//! Linear, finite-temperature Landau two-fluid propagation in a uniform channel.
//!
//! Pressure, temperature, normal velocity and superfluid velocity all evolve.
//! Thermal expansion couples first and second sound. The equations are
//! `A z_t = -D u_x`, `K u_t = -D z_x`, where `z=(p',T')`, `u=(v,v_n-v_s)`,
//! `A=[[kappa_T,-alpha],[-alpha,rho cp/T]]`, `K=diag(rho,rho f_s f_n)` and
//! `D=diag(1,rho f_s s)`. Thermodynamics comes from the original immutable
//! [`TwoFluidState`]; this is neither a gas nor a zero-viscosity LBM material.
//!
//! A staggered finite-volume incidence and implicit midpoint preserve quadratic
//! wave energy up to reported roundoff. A 2x2 energy normalization reduces the
//! coupled equations to two scalar tridiagonal systems, solved by the existing
//! `fs-la` owner in O(cells) work/storage. Buffers and factors are reused.
//! Each step checks the original equations, physical limits and energy balance
//! before publishing. Cancellation polls are bounded by 128 cells.
//!
//! Explicit mutual friction uses exact local material updates in a symmetric
//! split. Lost relative kinetic energy is exported as heat once; it is not
//! silently injected into this first-order thermal perturbation equation.
//! Reported storage is wave availability, not absolute internal energy.
//! The boundary thermal wave work is `T' q/T0`, distinct from transported heat.
//! Nonlinear heating/phase feedback, vortices, critical velocities, normal
//! viscous boundary layers, longitudinal viscosity, second viscosities,
//! Kapitza resistance and arbitrary geometry are outside this linear model.

use fs_la::tridiagonal::{TridiagonalError, TridiagonalWorkspace};
use fs_material::MaterialError;
use fs_material::superfluid::TwoFluidState;

mod boundary;
mod modes;
mod momentum;
mod step;
pub use boundary::{ChannelEnd, ImpermeableEndResponse, ImpermeableEndsResponse};
use modes::Modes;
pub use momentum::{BodyForce, ComponentForceDensity, ForcedChannelStep};

/// Source admission, caller domain, cancellation or numerical refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum ChannelError {
    Invalid(&'static str),
    Domain { field: &'static str, index: usize },
    Capacity,
    Cancelled,
    Material(MaterialError),
    Solve(TridiagonalError),
    Balance { defect_j: f64, tolerance_j: f64 },
    MomentumBalance { defect_n_s: f64, tolerance_n_s: f64 },
}
impl core::fmt::Display for ChannelError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "two-fluid channel: {self:?}")
    }
}
impl std::error::Error for ChannelError {}
impl From<MaterialError> for ChannelError {
    fn from(error: MaterialError) -> Self {
        Self::Material(error)
    }
}
impl From<TridiagonalError> for ChannelError {
    fn from(error: TridiagonalError) -> Self {
        if error == TridiagonalError::Cancelled {
            Self::Cancelled
        } else {
            Self::Solve(error)
        }
    }
}

/// A rejected channel proposal or the original thermal owner's refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum HeatStepError<E> {
    Channel(ChannelError),
    Heat(E),
}
impl<E: core::fmt::Display> core::fmt::Display for HeatStepError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Channel(error) => error.fmt(f),
            Self::Heat(error) => write!(f, "two-fluid heat acceptance: {error}"),
        }
    }
}
impl<E: std::error::Error + 'static> std::error::Error for HeatStepError<E> {}

/// Caller-selected bounds for the retained reference-state linearization.
/// These do not assert phase stability, a critical velocity or source validity.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LinearizationLimits {
    pub max_pressure_perturbation_pa: f64,
    pub max_temperature_perturbation_k: f64,
    pub max_component_speed_m_s: f64,
}

/// Read-only physical cell observables from the same accepted two-fluid state.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CellSample {
    pub pressure_perturbation_pa: f64,
    pub temperature_perturbation_k: f64,
    pub density_perturbation_kg_m3: f64,
    pub specific_entropy_perturbation_j_kg_k: f64,
}

/// Read-only velocity and transport observables. Positive is toward +x.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FaceSample {
    pub normal_velocity_m_s: f64,
    pub superfluid_velocity_m_s: f64,
    pub barycentric_velocity_m_s: f64,
    pub counterflow_velocity_m_s: f64,
    pub mass_flux_kg_m2_s: f64,
    pub entropy_flux_w_m2_k: f64,
    /// Heat relative to barycentric advection, not total enthalpy flux.
    pub counterflow_heat_flux_w_m2: f64,
}

/// Velocities prescribed for the full interval at one channel end.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct EndVelocities {
    pub normal_m_s: f64,
    pub superfluid_m_s: f64,
}

/// Both ends use +x orientation; inward flow is negative at the right end.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct BoundaryVelocities {
    pub left: EndVelocities,
    pub right: EndVelocities,
}
impl BoundaryVelocities {
    pub const CLOSED: Self = Self {
        left: EndVelocities {
            normal_m_s: 0.0,
            superfluid_m_s: 0.0,
        },
        right: EndVelocities {
            normal_m_s: 0.0,
            superfluid_m_s: 0.0,
        },
    };
}

/// Explicit local counterflow drag. No vortex density or coefficient is inferred.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum MutualFriction {
    #[default]
    None,
    Linear {
        gamma_kg_m3_s: f64,
    },
    GorterMellink {
        a_m_s_kg: f64,
    },
}

/// One accepted interval. Heat and small-signal boundary work have distinct roles.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ChannelStep {
    pub time_s: f64,
    pub wave_energy_before_j: f64,
    pub wave_energy_after_j: f64,
    pub boundary_wave_work_j: f64,
    /// Signed work supplied by authored body forces; friction heat is separate.
    pub body_force_work_j: f64,
    /// Midpoint pressure-perturbation impulses from the adjacent end CELLS,
    /// signed into the resolved interior: left +A integral(p), right -A integral(p).
    /// These are conjugate numerical port efforts. Physical end-half-cell
    /// inertia and body-force corrections belong to the boundary composition.
    pub boundary_pressure_impulse_n_s: [f64; 2],
    /// Actual interior-face momentum before and after the interval. Endpoint
    /// half-cell momentum is supplied by the mechanical boundary owner.
    pub resolved_momentum_before_kg_m_s: f64,
    pub resolved_momentum_after_kg_m_s: f64,
    /// Authored force impulse on the N-1 interior dual volumes, in channel +x.
    pub interior_body_impulse_n_s: f64,
    /// Actual resolved momentum change minus pressure and body impulses.
    pub resolved_momentum_defect_n_s: f64,
    pub mutual_friction_heat_j: f64,
    pub wave_balance_defect_j: f64,
    /// First-order mass imported through the two explicit ports.
    pub boundary_mass_kg: f64,
    /// First-order entropy imported, excluding second-order friction production.
    pub boundary_entropy_j_k: f64,
    /// Signed counterflow heat imported at the reference temperature.
    pub boundary_counterflow_heat_j: f64,
    /// Infinity-norm residual divided by the scale of the original linear solve.
    pub maximum_relative_linear_residual: f64,
}

#[derive(Debug)]
struct Field {
    cells: Vec<f64>,
    faces: Vec<f64>,
}
impl Field {
    fn new(n: usize) -> Result<Self, ChannelError> {
        Ok(Self {
            cells: filled(n)?,
            faces: filled(n + 1)?,
        })
    }
}

#[derive(Debug)]
struct SolveWorkspace {
    diagonal: Vec<f64>,
    off_diagonal: Vec<f64>,
    rhs: Vec<f64>,
    factor: TridiagonalWorkspace,
    factor_step: Option<f64>,
    boundary_inverse: [Option<[f64; 2]>; 2],
}
impl SolveWorkspace {
    fn new(n: usize) -> Result<Self, ChannelError> {
        Ok(Self {
            diagonal: filled(n)?,
            off_diagonal: filled(n - 1)?,
            rhs: filled(n)?,
            factor: TridiagonalWorkspace::new(n)?,
            factor_step: None,
            boundary_inverse: [None; 2],
        })
    }
}

/// Constant-area channel retaining one original thermodynamic reference state.
/// No allocation occurs in `step`; all proposed fields remain private on refusal.
#[derive(Debug)]
pub struct TwoFluidChannel {
    material: TwoFluidState,
    length_m: f64,
    area_m2: f64,
    dx: f64,
    cell_volume: f64,
    limits: LinearizationLimits,
    modes: Modes,
    accepted: [Field; 2],
    pending: [Field; 2],
    workspace: [SolveWorkspace; 2],
    last_heat_j: Vec<f64>,
    pending_heat_j: Vec<f64>,
    time_s: f64,
    total_heat_j: f64,
}

impl TwoFluidChannel {
    /// Create a quiescent uniform channel. Geometry is SI; cells must be >= 2.
    /// Allocation is fallible. Source-domain and linearization bounds are explicit.
    pub fn try_new(
        material: TwoFluidState,
        length_m: f64,
        area_m2: f64,
        cells: usize,
        limits: LinearizationLimits,
    ) -> Result<Self, ChannelError> {
        if cells < 2
            || cells
                .checked_mul(256)
                .is_none_or(|b| b > isize::MAX as usize)
            || ![
                length_m,
                area_m2,
                limits.max_pressure_perturbation_pa,
                limits.max_temperature_perturbation_k,
                limits.max_component_speed_m_s,
            ]
            .into_iter()
            .all(positive)
        {
            return Err(ChannelError::Invalid(
                "positive geometry and perturbation bounds, with at least two cells, are required",
            ));
        }
        let dx = length_m / cells as f64;
        let cell_volume = dx * area_m2;
        let modes = Modes::new(&material)?;
        if ![dx, cell_volume, length_m * area_m2, dx / modes.speed[0]]
            .into_iter()
            .all(positive)
        {
            return Err(ChannelError::Invalid(
                "channel geometry or wave time is not representable",
            ));
        }
        Ok(Self {
            material,
            length_m,
            area_m2,
            dx,
            cell_volume,
            limits,
            modes,
            accepted: [Field::new(cells)?, Field::new(cells)?],
            pending: [Field::new(cells)?, Field::new(cells)?],
            workspace: [SolveWorkspace::new(cells)?, SolveWorkspace::new(cells)?],
            last_heat_j: filled(cells)?,
            pending_heat_j: filled(cells)?,
            time_s: 0.0,
            total_heat_j: 0.0,
        })
    }

    /// Initialize physical perturbations and both component velocities.
    /// There are N cell values and N+1 face values. Allowed only before stepping.
    /// Refusal/cancellation leaves the original fields intact.
    pub fn initialize(
        &mut self,
        pressure_pa: &[f64],
        temperature_k: &[f64],
        normal_m_s: &[f64],
        superfluid_m_s: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), ChannelError> {
        poll(&mut cancelled)?;
        let n = self.cells();
        if self.time_s != 0.0
            || pressure_pa.len() != n
            || temperature_k.len() != n
            || normal_m_s.len() != n + 1
            || superfluid_m_s.len() != n + 1
        {
            return Err(ChannelError::Invalid(
                "initialization needs N cells and N+1 faces before the first step",
            ));
        }
        for i in 0..=n {
            if i % 128 == 0 {
                poll(&mut cancelled)?;
            }
            if i < n {
                let value = self.modes.cell_to_modes(pressure_pa[i], temperature_k[i]);
                self.pending[0].cells[i] = value[0];
                self.pending[1].cells[i] = value[1];
            }
            let v = self.material.normal_fraction() * normal_m_s[i]
                + self.material.superfluid_fraction() * superfluid_m_s[i];
            let w = normal_m_s[i] - superfluid_m_s[i];
            let value = self.modes.face_to_modes(v, w);
            self.pending[0].faces[i] = value[0];
            self.pending[1].faces[i] = value[1];
        }
        self.validate_pending(&mut cancelled)?;
        energy(&self.pending, self.cell_volume, &mut cancelled)?;
        poll(&mut cancelled)?;
        core::mem::swap(&mut self.accepted, &mut self.pending);
        Ok(())
    }

    pub const fn material(&self) -> &TwoFluidState {
        &self.material
    }
    pub fn cells(&self) -> usize {
        self.accepted[0].cells.len()
    }
    pub const fn length_m(&self) -> f64 {
        self.length_m
    }
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub const fn cumulative_mutual_friction_heat_j(&self) -> f64 {
        self.total_heat_j
    }
    /// Per-cell heat exported by the last accepted step, ready for a thermal owner.
    /// An interior face's heat is split equally between its two adjacent cells.
    pub fn last_mutual_friction_heat_j(&self) -> &[f64] {
        &self.last_heat_j
    }
    /// Courant <= 1 accuracy policy; implicit midpoint is not a time-error estimator.
    pub fn maximum_step_s(&self) -> f64 {
        self.dx / self.modes.speed[0]
    }

    pub fn cell(&self, index: usize) -> Option<CellSample> {
        Some(self.modes.cell(
            &self.material,
            [
                *self.accepted[0].cells.get(index)?,
                *self.accepted[1].cells.get(index)?,
            ],
        ))
    }
    pub fn face(&self, index: usize) -> Option<FaceSample> {
        Some(self.modes.face(
            &self.material,
            [
                *self.accepted[0].faces.get(index)?,
                *self.accepted[1].faces.get(index)?,
            ],
        ))
    }
    pub fn wave_energy_j(&self) -> Result<f64, ChannelError> {
        energy(&self.accepted, self.cell_volume, &mut || false)
    }

    /// Advance both sound modes, honoring imposed component flows at both ends.
    /// Callers can bridge the cancellation callback to their existing execution Cx.
    pub fn step(
        &mut self,
        dt_s: f64,
        boundary: BoundaryVelocities,
        friction: MutualFriction,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ChannelStep, ChannelError> {
        self.step_forced(dt_s, boundary, BodyForce::None, friction, cancelled)
    }

    /// Prepare the original thermal owner's candidate before publishing a wave step.
    ///
    /// The callback receives the proposed report and per-cell frictional joules.
    /// It must return a private candidate without changing external accepted state.
    /// A thermal refusal or the final cancellation drops that candidate and leaves
    /// channel fields, clock and heat unchanged. On success the returned candidate
    /// can be installed by an infallible assignment, completing an atomic coupled
    /// transaction. This does not change the channel's frozen reference properties.
    pub fn step_with_heat<T, E>(
        &mut self,
        dt_s: f64,
        boundary: BoundaryVelocities,
        friction: MutualFriction,
        cancelled: impl FnMut() -> bool,
        prepare_heat: impl FnOnce(&ChannelStep, &[f64]) -> Result<T, E>,
    ) -> Result<(ChannelStep, T), HeatStepError<E>> {
        self.step_forced_with_heat(
            dt_s,
            boundary,
            BodyForce::None,
            friction,
            cancelled,
            prepare_heat,
        )
    }

    fn publish_step(&mut self, report: ChannelStep) {
        // All arithmetic, capacity and cancellation checks precede this point.
        // stage_step already admitted the cumulative sum; publication cannot fail.
        core::mem::swap(&mut self.accepted, &mut self.pending);
        core::mem::swap(&mut self.last_heat_j, &mut self.pending_heat_j);
        self.time_s = report.time_s;
        self.total_heat_j += report.mutual_friction_heat_j;
    }

    fn validate_pending(&self, cancelled: &mut impl FnMut() -> bool) -> Result<(), ChannelError> {
        self.validate_pending_with_limits(true, cancelled)
    }

    fn validate_pending_with_limits(
        &self,
        physical_limits: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), ChannelError> {
        let n = self.cells();
        for i in 0..=n {
            if i % 128 == 0 {
                poll(cancelled)?;
            }
            if i < n {
                let c = self.modes.cell(
                    &self.material,
                    [self.pending[0].cells[i], self.pending[1].cells[i]],
                );
                if !c.pressure_perturbation_pa.is_finite()
                    || (physical_limits
                        && c.pressure_perturbation_pa.abs()
                            > self.limits.max_pressure_perturbation_pa)
                {
                    return Err(ChannelError::Domain {
                        field: "pressure",
                        index: i,
                    });
                }
                if !c.temperature_perturbation_k.is_finite()
                    || (physical_limits
                        && (c.temperature_perturbation_k.abs()
                            > self.limits.max_temperature_perturbation_k
                            || c.temperature_perturbation_k.abs()
                                > 0.1 * self.material.temperature_k()))
                {
                    return Err(ChannelError::Domain {
                        field: "temperature",
                        index: i,
                    });
                }
                if !c.density_perturbation_kg_m3.is_finite()
                    || (physical_limits
                        && c.density_perturbation_kg_m3.abs() > 0.1 * self.material.density_kg_m3())
                    || !c.specific_entropy_perturbation_j_kg_k.is_finite()
                {
                    return Err(ChannelError::Domain {
                        field: "linearized density/entropy",
                        index: i,
                    });
                }
            }
            let face = self.modes.face(
                &self.material,
                [self.pending[0].faces[i], self.pending[1].faces[i]],
            );
            for speed in [face.normal_velocity_m_s, face.superfluid_velocity_m_s] {
                if !speed.is_finite()
                    || (physical_limits
                        && (speed.abs() > self.limits.max_component_speed_m_s
                            || speed.abs() > 0.1 * self.modes.speed[1]))
                {
                    return Err(ChannelError::Domain {
                        field: "component velocity",
                        index: i,
                    });
                }
            }
            if [
                face.mass_flux_kg_m2_s,
                face.entropy_flux_w_m2_k,
                face.counterflow_heat_flux_w_m2,
            ]
            .into_iter()
            .any(|x| !x.is_finite())
            {
                return Err(ChannelError::Domain {
                    field: "transport flux",
                    index: i,
                });
            }
        }
        Ok(())
    }
}

fn filled(n: usize) -> Result<Vec<f64>, ChannelError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(n)
        .map_err(|_| ChannelError::Capacity)?;
    values.resize(n, 0.0);
    Ok(values)
}
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), ChannelError> {
    if cancelled() {
        Err(ChannelError::Cancelled)
    } else {
        Ok(())
    }
}
fn accumulate(total: &mut f64, correction: &mut f64, value: f64) {
    let next = value - *correction;
    let sum = *total + next;
    *correction = (sum - *total) - next;
    *total = sum;
}
fn energy(
    fields: &[Field; 2],
    volume: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, ChannelError> {
    let mut total = 0.0;
    let mut correction = 0.0;
    for field in fields {
        for (i, &value) in field
            .cells
            .iter()
            .chain(&field.faces[1..field.faces.len() - 1])
            .enumerate()
        {
            if i % 128 == 0 {
                poll(cancelled)?;
            }
            accumulate(&mut total, &mut correction, (0.5 * value) * value * volume);
        }
    }
    if total.is_finite() {
        Ok(total)
    } else {
        Err(ChannelError::Invalid("wave energy is not representable"))
    }
}
