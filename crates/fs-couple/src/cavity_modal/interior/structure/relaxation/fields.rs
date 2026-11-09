//! Original physical recovery around the shared complete-state harmonic owner.

use super::super::{poles::relaxation::stress_square_over, response};
use fs_math::{c64::C64, dd::Dd, det};
use fs_modal::relaxation_modes::MaxwellDampingBranch;
use fs_phs::harmonic::{HarmonicBudget, HarmonicError, HarmonicFactor, HarmonicModel};
use fs_phs::mechanical::{LinearMechanicalSystem, MAX_MECHANICAL_STATES};
use fs_phs::{PhsError, RelaxationBranch};

/// Complete dense source admission; no budget silently discards coordinates.
#[derive(Clone, Copy, Debug)]
pub struct StructuralInteriorRelaxationHarmonicLimits {
    pub maximum_coordinates: usize,
    pub maximum_branches: usize,
    /// Adapter preparation/recovery payload, separate from the native harmonic
    /// budget. Original borrowed arrays, earlier results and allocator metadata
    /// are excluded. For d=2*n+b, bound is
    /// 8*(32*d*d+16*n*n+8*b*n+64*d+64*n+32*b) bytes.
    pub maximum_adapter_bytes: u128,
    /// Conservative adapter work score 1024*d^3+128*n^2+128*b*n,
    /// covering mechanical admission and original-row/physical recovery.
    pub maximum_adapter_work: u128,
    pub numerical: HarmonicBudget,
    /// Additional componentwise check of the ORIGINAL mechanical, kinematic
    /// and branch equations, including their differentiated equations.
    pub physical_relative_residual_tolerance: f64,
}

impl Default for StructuralInteriorRelaxationHarmonicLimits {
    fn default() -> Self {
        Self {
            maximum_coordinates: 64,
            maximum_branches: 32,
            maximum_adapter_bytes: 64 * 1024 * 1024,
            maximum_adapter_work: 32_000_000_000,
            numerical: HarmonicBudget {
                max_dimension: MAX_MECHANICAL_STATES,
                max_dense_bytes: 64 * 1024 * 1024,
                max_dense_work: 1_000_000_000,
                relative_residual_tolerance: 1e-10,
            },
            physical_relative_residual_tolerance: 1e-9,
        }
    }
}

/// Original-source, native-solver or physical-row refusal; no partial response.
#[derive(Debug)]
pub enum StructuralInteriorRelaxationHarmonicError {
    Invalid(&'static str),
    Limit {
        quantity: &'static str,
        requested: u128,
        maximum: u128,
    },
    Dynamics(PhsError),
    Harmonic(HarmonicError),
    PhysicalResidual {
        relative: f64,
        tolerance: f64,
    },
    Cancelled,
}
impl core::fmt::Display for StructuralInteriorRelaxationHarmonicError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(what) => write!(f, "structural relaxation response: {what}"),
            Self::Limit {
                quantity,
                requested,
                maximum,
            } => write!(
                f,
                "structural relaxation {quantity}: {requested} exceeds {maximum}"
            ),
            Self::Dynamics(error) => error.fmt(f),
            Self::Harmonic(error) => error.fmt(f),
            Self::PhysicalResidual {
                relative,
                tolerance,
            } => write!(
                f,
                "original relaxation residual {relative:e} exceeds {tolerance:e}"
            ),
            Self::Cancelled => f.write_str("structural relaxation response cancelled"),
        }
    }
}
impl std::error::Error for StructuralInteriorRelaxationHarmonicError {}
impl From<PhsError> for StructuralInteriorRelaxationHarmonicError {
    fn from(error: PhsError) -> Self {
        Self::Dynamics(error)
    }
}
impl From<HarmonicError> for StructuralInteriorRelaxationHarmonicError {
    fn from(error: HarmonicError) -> Self {
        Self::Harmonic(error)
    }
}
impl From<response::ResponseError> for StructuralInteriorRelaxationHarmonicError {
    fn from(error: response::ResponseError) -> Self {
        match error {
            response::ResponseError::Invalid(what) => Self::Invalid(what),
            response::ResponseError::Cancelled => Self::Cancelled,
        }
    }
}
type Error = StructuralInteriorRelaxationHarmonicError;
type Limits = StructuralInteriorRelaxationHarmonicLimits;

mod curvature;
mod parameter;
pub(super) mod transfer;
pub use transfer::{
    StructuralInteriorRelaxationFixedReceiver, StructuralInteriorRelaxationReceiverJet,
    StructuralInteriorRelaxationReceiverSpectrum,
};
pub(super) mod power;
pub use curvature::{
    StructuralInteriorRelaxationCurvatureResponse, StructuralInteriorRelaxationMemoryCurvature,
    StructuralInteriorRelaxationPowerCurvature,
};
pub(super) use parameter::PreparedDirection;
pub use parameter::{
    StructuralInteriorRelaxationFrequencyParameterResponse,
    StructuralInteriorRelaxationMemoryParameterJet, StructuralInteriorRelaxationParameterDirection,
    StructuralInteriorRelaxationPowerParameterJet,
};
pub use power::{
    StructuralInteriorRelaxationMemoryFrequency, StructuralInteriorRelaxationPowerFrequency,
};

/// One original Maxwell arm in caller order, using peak phasors.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StructuralInteriorRelaxationMemoryResponse {
    pub stress: C64,
    /// Cycle-mean |sigma|^2/(4*k) [J], separate from relaxed spatial storage.
    pub mean_stored_energy_j: f64,
    /// Independent positive heat |sigma|^2/(2*k*tau) [W].
    pub average_dissipation_w: f64,
    /// Signed Re(conj(sigma)*b_structure^T*v)/2 [W].
    pub average_structural_work_w: f64,
    /// Signed Re(conj(sigma)*b_acoustic^T*p)/2 [W].
    pub average_acoustic_work_w: f64,
}
impl StructuralInteriorRelaxationMemoryResponse {
    pub fn power_balance_residual_w(self) -> f64 {
        self.average_structural_work_w + self.average_acoustic_work_w - self.average_dissipation_w
    }
}

/// Complete original fields and independent cycle-averaged physical accounts.
#[derive(Clone, Debug)]
pub struct StructuralInteriorRelaxationResponse {
    pub angular_frequency_rad_s: f64,
    pub structural_displacement: Vec<C64>,
    pub structural_velocity: Vec<C64>,
    pub pressure_potential_pa_s: Vec<C64>,
    /// Physical pressure, equal to -i*omega*psi within the original row tolerance.
    pub pressure_pa: Vec<C64>,
    /// B^T*p; branch reactions are accounted separately below.
    pub fluid_force: Vec<C64>,
    /// B*v; this is the original interface flow, without explicit branch flow.
    pub displaced_flow_m3_s: Vec<C64>,
    pub spatial_energy: response::StructuralInteriorEnergy,
    pub memory: Vec<StructuralInteriorRelaxationMemoryResponse>,
    /// Sum of the listed arm energies; not an additional energy debit.
    pub relaxation_energy_j: f64,
    pub average_structural_input_power_w: f64,
    pub average_acoustic_input_power_w: f64,
    /// Positive from fluid into structure through the original B/B^T interface.
    pub average_fluid_to_structure_power_w: f64,
    pub average_structural_dissipation_w: f64,
    pub average_wall_dissipation_w: f64,
    pub average_relaxation_dissipation_w: f64,
    pub average_structural_to_relaxation_power_w: f64,
    pub average_acoustic_to_relaxation_power_w: f64,
    /// Independently evaluated source minus all physical heat [W].
    pub power_balance_residual_w: f64,
    /// Maximum native canonical and original physical backward residual.
    pub relative_residual: f64,
    pub original_relative_residual: f64,
}
impl StructuralInteriorRelaxationResponse {
    pub fn mean_stored_energy_j(&self) -> f64 {
        self.spatial_energy.total_j() + self.relaxation_energy_j
    }
    pub fn structural_power_balance_residual_w(&self) -> f64 {
        self.average_structural_input_power_w + self.average_fluid_to_structure_power_w
            - self.average_structural_dissipation_w
            - self.average_structural_to_relaxation_power_w
    }
    pub fn acoustic_power_balance_residual_w(&self) -> f64 {
        self.average_acoustic_input_power_w
            - self.average_fluid_to_structure_power_w
            - self.average_wall_dissipation_w
            - self.average_acoustic_to_relaxation_power_w
    }
    pub fn relaxation_power_balance_residual_w(&self) -> f64 {
        self.average_structural_to_relaxation_power_w + self.average_acoustic_to_relaxation_power_w
            - self.average_relaxation_dissipation_w
    }
}

/// Same-factor continuous frequency tangent, in original field units per rad/s.
#[derive(Clone, Debug)]
pub struct StructuralInteriorRelaxationFrequencyResponse {
    pub response: StructuralInteriorRelaxationResponse,
    pub structural_displacement_derivative_per_rad_s: Vec<C64>,
    pub structural_velocity_derivative_per_rad_s: Vec<C64>,
    pub pressure_potential_derivative_pa_s_per_rad_s: Vec<C64>,
    pub pressure_derivative_pa_per_rad_s: Vec<C64>,
    pub fluid_force_derivative_per_rad_s: Vec<C64>,
    pub displaced_flow_derivative_m3_s_per_rad_s: Vec<C64>,
    pub relaxation_stress_derivative_per_rad_s: Vec<C64>,
    /// Independent physical heat/work, evaluated with the same two state solves.
    pub power: StructuralInteriorRelaxationPowerFrequency,
    pub memory: Vec<StructuralInteriorRelaxationMemoryFrequency>,
    /// Maximum native and original differentiated-equation backward residual.
    pub relative_derivative_residual: f64,
    pub original_relative_derivative_residual: f64,
}

#[derive(Clone, Copy)]
pub(super) struct Source<'s> {
    pub structural: usize,
    pub coordinates: usize,
    pub mass: &'s [f64],
    pub damping: &'s [f64],
    pub stiffness: &'s [f64],
    pub coupling: &'s [f64],
}
pub(super) struct Prepared<'s> {
    source: Source<'s>,
    branches: &'s [MaxwellDampingBranch],
    roots: Vec<f64>,
    harmonic: HarmonicModel,
    mechanical: LinearMechanicalSystem,
    limits: Limits,
}
pub(super) struct Factor<'p, 's> {
    prepared: &'p Prepared<'s>,
    harmonic: HarmonicFactor<'p>,
    omega: f64,
}
struct State {
    q: Vec<C64>,
    velocity: Vec<C64>,
    stress: Vec<C64>,
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}

fn admit(n: usize, b: usize, limits: Limits) -> Result<(), Error> {
    let (n, b) = (n as u128, b as u128);
    let d = 2 * n + b;
    for (quantity, requested, maximum) in [
        ("coordinates", n, limits.maximum_coordinates as u128),
        ("branches", b, limits.maximum_branches as u128),
        ("native states", d, MAX_MECHANICAL_STATES as u128),
    ] {
        if requested > maximum {
            return Err(Error::Limit {
                quantity,
                requested,
                maximum,
            });
        }
    }
    if n == 0
        || !limits.physical_relative_residual_tolerance.is_finite()
        || limits.physical_relative_residual_tolerance <= 0.0
        || limits.physical_relative_residual_tolerance >= 1.0
    {
        return Err(Error::Invalid(
            "nonempty coordinates and physical residual tolerance in (0,1) required",
        ));
    }
    let bytes = 8 * (32 * d * d + 16 * n * n + 8 * b * n + 64 * d + 64 * n + 32 * b);
    let work = 1024 * d * d * d + 128 * n * n + 128 * b * n;
    for (quantity, requested, maximum) in [
        ("adapter bytes", bytes, limits.maximum_adapter_bytes),
        ("adapter work", work, limits.maximum_adapter_work),
    ] {
        if requested > maximum {
            return Err(Error::Limit {
                quantity,
                requested,
                maximum,
            });
        }
    }
    Ok(())
}

impl<'s> Prepared<'s> {
    pub(super) fn new(
        source: Source<'s>,
        branches: &'s [MaxwellDampingBranch],
        limits: Limits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(cancelled)?;
        let (n, ns) = (source.coordinates, source.structural);
        admit(n, branches.len(), limits)?;
        if ns == 0
            || ns >= n
            || source.coupling.len() != ns * (n - ns)
            || source.coupling.iter().any(|v| !v.is_finite())
        {
            return Err(Error::Invalid(
                "complete finite original interface required",
            ));
        }
        if source.mass.len() != n * n
            || source.stiffness.len() != n * n
            || source.damping.len() != n * n
        {
            return Err(Error::Invalid("complete original M/C/K arrays required"));
        }
        for i in 0..n {
            poll(cancelled)?;
            if (0..i).any(|j| {
                source.mass[i * n + j] != source.mass[j * n + i]
                    || source.stiffness[i * n + j] != source.stiffness[j * n + i]
            }) {
                return Err(Error::Invalid(
                    "original mass and relaxed stiffness must be symmetric",
                ));
            }
        }
        let mut identity = vec![0.0; n * n];
        for i in 0..n {
            identity[i * n + i] = 1.0;
        }
        let mut native = Vec::with_capacity(branches.len());
        let mut roots = Vec::with_capacity(branches.len());
        for branch in branches {
            poll(cancelled)?;
            if branch.strain.len() != n
                || branch.strain.iter().any(|v| !v.is_finite())
                || !branch.stiffness.is_finite()
                || branch.stiffness <= 0.0
                || !branch.relaxation_time_s.is_finite()
                || branch.relaxation_time_s <= 0.0
            {
                return Err(Error::Invalid(
                    "complete original strains and positive finite branch coefficients required",
                ));
            }
            native.push(RelaxationBranch {
                projection: branch.strain.clone(),
                stiffness: branch.stiffness,
                relaxation_time_s: branch.relaxation_time_s,
            });
            roots.push(det::sqrt(branch.stiffness));
        }
        let mechanical = LinearMechanicalSystem::new(
            source.mass,
            source.damping,
            source.stiffness,
            &identity,
            n,
            n,
            native,
        )?;
        poll(cancelled)?;
        let harmonic =
            HarmonicModel::from_ode(mechanical.system(), limits.numerical, &mut *cancelled)?;
        poll(cancelled)?;
        Ok(Self {
            source,
            branches,
            roots,
            harmonic,
            mechanical,
            limits,
        })
    }

    pub(super) fn factor(
        &self,
        omega: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Factor<'_, 's>, Error> {
        let harmonic = self.harmonic.factor(omega, &mut *cancelled)?;
        Ok(Factor {
            prepared: self,
            harmonic,
            omega,
        })
    }

    fn load(&self, structural: &[C64], acoustic: &[C64]) -> Result<Vec<C64>, Error> {
        if structural.len() != self.source.structural
            || acoustic.len() != self.source.coordinates - self.source.structural
            || structural.iter().chain(acoustic).any(|&v| !finite(v))
        {
            return Err(Error::Invalid(
                "one finite force or inward-flow phasor per original coordinate required",
            ));
        }
        let mut result = Vec::with_capacity(self.source.coordinates);
        result.extend_from_slice(structural);
        result.extend_from_slice(acoustic);
        Ok(result)
    }

    fn state(
        &self,
        state: &[C64],
        effort: &[C64],
        output: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<State, Error> {
        let n = self.source.coordinates;
        let mut stress = Vec::with_capacity(self.branches.len());
        for (j, &root) in self.roots.iter().enumerate() {
            poll(cancelled)?;
            let value = effort[2 * n + j].scale(-root);
            if !finite(value) {
                return Err(Error::Invalid("original branch stress is unrepresentable"));
            }
            stress.push(value);
        }
        Ok(State {
            q: state[..n].to_vec(),
            velocity: output.to_vec(),
            stress,
        })
    }

    fn residual(
        &self,
        omega: f64,
        state: &State,
        input: &[C64],
        baseline: Option<(&State, f64)>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, Error> {
        let n = self.source.coordinates;
        let mut residual = 0.0_f64;
        for i in 0..n {
            poll(cancelled)?;
            let mut row = Row::new();
            row.add(state.q[i], C64::new(0.0, -omega));
            row.add(state.velocity[i], C64::from_re(-1.0));
            if let Some((base, order)) = baseline {
                row.add(base.q[i], C64::new(0.0, -order));
            }
            residual = residual.max(row.relative()?);
            let mut row = Row::new();
            row.add(input[i], C64::from_re(-1.0));
            for j in 0..n {
                row.add(state.q[j], C64::from_re(self.source.stiffness[i * n + j]));
                row.add(
                    state.velocity[j],
                    C64::new(
                        self.source.damping[i * n + j],
                        -omega * self.source.mass[i * n + j],
                    ),
                );
                if let Some((base, order)) = baseline {
                    row.add(
                        base.velocity[j],
                        C64::new(0.0, -order * self.source.mass[i * n + j]),
                    );
                }
            }
            for (branch, &stress) in self.branches.iter().zip(&state.stress) {
                row.add(stress, C64::from_re(branch.strain[i]));
            }
            residual = residual.max(row.relative()?);
        }
        for (j, branch) in self.branches.iter().enumerate() {
            poll(cancelled)?;
            let mut row = Row::new();
            row.add(
                state.stress[j],
                C64::new(branch.relaxation_time_s.recip(), -omega),
            );
            for i in 0..n {
                row.add(
                    state.velocity[i],
                    C64::from_re(-branch.stiffness * branch.strain[i]),
                );
            }
            if let Some((base, order)) = baseline {
                row.add(base.stress[j], C64::new(0.0, -order));
            }
            residual = residual.max(row.relative()?);
        }
        let tolerance = self.limits.physical_relative_residual_tolerance;
        if residual > tolerance {
            return Err(Error::PhysicalResidual {
                relative: residual,
                tolerance,
            });
        }
        Ok(residual)
    }

    fn observe(
        &self,
        omega: f64,
        mut state: State,
        input: &[C64],
        native_residual: f64,
        original_residual: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<StructuralInteriorRelaxationResponse, Error> {
        let (n, ns) = (self.source.coordinates, self.source.structural);
        let potential = state.q.split_off(ns);
        let pressure = state.velocity.split_off(ns);
        let (fluid_force, flow) =
            response::interface(self.source.coupling, &state.velocity, &pressure, cancelled)?;
        let energy = response::energy(
            self.source.mass,
            self.source.stiffness,
            &state.q,
            &state.velocity,
            &potential,
            &pressure,
            cancelled,
        )?;
        let structural_loss =
            response::dissipation(self.source.damping, n, 0, &state.velocity, cancelled)?;
        let wall_loss = response::dissipation(self.source.damping, n, ns, &pressure, cancelled)?;
        let mut memory = Vec::with_capacity(self.branches.len());
        let mut totals = [Dd::from_f64(0.0); 4];
        for (branch, &stress) in self.branches.iter().zip(&state.stress) {
            poll(cancelled)?;
            let mean_stored_energy_j =
                stress_square_over(stress, branch.stiffness, 4.0).map_err(Error::Invalid)?;
            let average_dissipation_w = 0.5
                * stress_square_over(stress, branch.stiffness, branch.relaxation_time_s)
                    .map_err(Error::Invalid)?;
            let structural_work = branch_work(stress, &branch.strain[..ns], &state.velocity)?;
            let acoustic_work = branch_work(stress, &branch.strain[ns..], &pressure)?;
            for (sum, value) in totals.iter_mut().zip([
                mean_stored_energy_j,
                average_dissipation_w,
                structural_work,
                acoustic_work,
            ]) {
                *sum = *sum + Dd::from_f64(value);
            }
            memory.push(StructuralInteriorRelaxationMemoryResponse {
                stress,
                mean_stored_energy_j,
                average_dissipation_w,
                average_structural_work_w: structural_work,
                average_acoustic_work_w: acoustic_work,
            });
        }
        let totals = totals.map(Dd::to_f64);
        let structural_input = response::work(&input[..ns], &state.velocity);
        let acoustic_input = response::work(&input[ns..], &pressure);
        let exchange = response::work(&fluid_force, &state.velocity);
        let residual = (Dd::from_f64(structural_input) + Dd::from_f64(acoustic_input)
            - Dd::from_f64(structural_loss)
            - Dd::from_f64(wall_loss)
            - Dd::from_f64(totals[1]))
        .to_f64();
        if totals
            .iter()
            .chain(&[
                structural_input,
                acoustic_input,
                exchange,
                structural_loss,
                wall_loss,
                residual,
            ])
            .any(|v| !v.is_finite())
        {
            return Err(Error::Invalid(
                "original storage or physical power is unrepresentable",
            ));
        }
        poll(cancelled)?;
        Ok(StructuralInteriorRelaxationResponse {
            angular_frequency_rad_s: omega,
            structural_displacement: state.q,
            structural_velocity: state.velocity,
            pressure_potential_pa_s: potential,
            pressure_pa: pressure,
            fluid_force,
            displaced_flow_m3_s: flow,
            spatial_energy: energy,
            memory,
            relaxation_energy_j: totals[0],
            average_structural_input_power_w: structural_input,
            average_acoustic_input_power_w: acoustic_input,
            average_fluid_to_structure_power_w: exchange,
            average_structural_dissipation_w: structural_loss,
            average_wall_dissipation_w: wall_loss,
            average_relaxation_dissipation_w: totals[1],
            average_structural_to_relaxation_power_w: totals[2],
            average_acoustic_to_relaxation_power_w: totals[3],
            power_balance_residual_w: residual,
            relative_residual: native_residual.max(original_residual),
            original_relative_residual: original_residual,
        })
    }
}

impl Factor<'_, '_> {
    pub(super) fn solve(
        &self,
        structural: &[C64],
        acoustic: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<StructuralInteriorRelaxationResponse, Error> {
        poll(cancelled)?;
        let input = self.prepared.load(structural, acoustic)?;
        let native = self.harmonic.response(&input, &mut *cancelled)?;
        let state =
            self.prepared
                .state(&native.state, &native.effort, &native.output, cancelled)?;
        let residual = self
            .prepared
            .residual(self.omega, &state, &input, None, cancelled)?;
        self.prepared.observe(
            self.omega,
            state,
            &input,
            native.relative_residual,
            residual,
            cancelled,
        )
    }

    pub(super) fn solve_with_frequency_derivative(
        &self,
        structural: &[C64],
        acoustic: &[C64],
        source_derivative: Option<(&[C64], &[C64])>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<StructuralInteriorRelaxationFrequencyResponse, Error> {
        poll(cancelled)?;
        let p = self.prepared;
        let input = p.load(structural, acoustic)?;
        let input_derivative = source_derivative.map(|(s, a)| p.load(s, a)).transpose()?;
        let native = self.harmonic.response_with_frequency_derivative(
            &input,
            input_derivative.as_deref(),
            &mut *cancelled,
        )?;
        self.recover_frequency(native, &input, input_derivative.as_deref(), cancelled)
    }

    // Reuse the same original physical recovery after a native curvature solve.
    fn recover_frequency(
        &self,
        native: fs_phs::harmonic::HarmonicFrequencyResponse,
        input: &[C64],
        input_derivative: Option<&[C64]>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<StructuralInteriorRelaxationFrequencyResponse, Error> {
        let p = self.prepared;
        let state = p.state(
            &native.response.state,
            &native.response.effort,
            &native.response.output,
            cancelled,
        )?;
        let mut derivative = p.state(
            &native.state_derivative,
            &native.effort_derivative,
            &native.output_derivative,
            cancelled,
        )?;
        let residual = p.residual(self.omega, &state, &input, None, cancelled)?;
        let zero = vec![C64::ZERO; p.source.coordinates];
        let slope_input = input_derivative.unwrap_or(&zero);
        let derivative_residual = p.residual(
            self.omega,
            &derivative,
            &slope_input,
            Some((&state, 1.0)),
            cancelled,
        )?;
        let response = p.observe(
            self.omega,
            state,
            &input,
            native.response.relative_residual,
            residual,
            cancelled,
        )?;
        let (force, flow) = response::interface(
            p.source.coupling,
            &derivative.velocity[..p.source.structural],
            &derivative.velocity[p.source.structural..],
            cancelled,
        )?;
        let (power, memory) = power::observe(
            p,
            &response,
            &derivative,
            &force,
            input,
            slope_input,
            &native.response.effort,
            &native.effort_derivative,
            cancelled,
        )?;
        let potential = derivative.q.split_off(p.source.structural);
        let pressure = derivative.velocity.split_off(p.source.structural);
        poll(cancelled)?;
        Ok(StructuralInteriorRelaxationFrequencyResponse {
            response,
            structural_displacement_derivative_per_rad_s: derivative.q,
            structural_velocity_derivative_per_rad_s: derivative.velocity,
            pressure_potential_derivative_pa_s_per_rad_s: potential,
            pressure_derivative_pa_per_rad_s: pressure,
            fluid_force_derivative_per_rad_s: force,
            displaced_flow_derivative_m3_s_per_rad_s: flow,
            relaxation_stress_derivative_per_rad_s: derivative.stress,
            power,
            memory,
            relative_derivative_residual: native
                .relative_derivative_residual
                .max(derivative_residual),
            original_relative_derivative_residual: derivative_residual,
        })
    }
}

fn branch_work(stress: C64, weights: &[f64], field: &[C64]) -> Result<f64, Error> {
    let mut work = Dd::from_f64(0.0);
    for (&w, &z) in weights.iter().zip(field) {
        work = work
            + Dd::from_f64(w)
                * (Dd::from_f64(stress.re) * Dd::from_f64(z.re)
                    + Dd::from_f64(stress.im) * Dd::from_f64(z.im));
    }
    let value = (Dd::from_f64(0.5) * work).to_f64();
    if value.is_finite() {
        Ok(value)
    } else {
        Err(Error::Invalid("original branch work is unrepresentable"))
    }
}

struct Row {
    re: Dd,
    im: Dd,
    scale: f64,
}
impl Row {
    fn new() -> Self {
        Self {
            re: Dd::from_f64(0.0),
            im: Dd::from_f64(0.0),
            scale: 0.0,
        }
    }
    fn add(&mut self, z: C64, a: C64) {
        self.re = self.re + Dd::from_f64(a.re) * Dd::from_f64(z.re)
            - Dd::from_f64(a.im) * Dd::from_f64(z.im);
        self.im = self.im
            + Dd::from_f64(a.re) * Dd::from_f64(z.im)
            + Dd::from_f64(a.im) * Dd::from_f64(z.re);
        self.scale += a.abs() * z.abs();
    }
    fn relative(self) -> Result<f64, Error> {
        let error = C64::new(self.re.to_f64(), self.im.to_f64()).abs();
        if !error.is_finite() || !self.scale.is_finite() {
            return Err(Error::Invalid(
                "original equation residual is unrepresentable",
            ));
        }
        Ok(if self.scale == 0.0 {
            if error == 0.0 { 0.0 } else { 1.0 }
        } else {
            error / self.scale
        })
    }
}

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
