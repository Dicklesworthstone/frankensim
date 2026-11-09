//! Resolved unsteady Newtonian films on a tangentially moving wall.
//!
//! The locally parallel equation is `rho*u_t = rho*g + (mu*u_y)_y`, with
//! `u(0)-U_wall=b*u_y(0)` and `mu*u_y(h)=tau_surface`. The default slip length
//! `b=0` is no-slip. Uniform finite volumes retain the
//! velocity through the film; implicit midpoint advances the original viscous
//! diffusion, rather than instantaneously replacing it with a Nusselt profile.
//! The flat-film reduction follows the tangential momentum/stress conditions
//! in <https://www.damp.tottori-u.ac.jp/~ooshida/article/RIMS0007/en_node1.html>.
//!
//! Each interval holds its authored acceleration, wall speed and surface
//! traction. For a smooth drive, supply their midpoint values and resolve its
//! bandwidth. The discrete work identity is body + wall + surface work =
//! kinetic change + bulk viscous heat + wall slip loss, up to the reported
//! linear/roundoff residual. Slip length is prescribed, not inferred from a
//! liquid material or a contact angle. Interface loss is not bulk fluid heating.
//! Heat is reported, not deposited into another material or thermal inventory.
//!
//! Geometry, density and viscosity are fixed. This is a laterally wide,
//! streamwise uniform laminar reduction: no evolving thickness, entrance flow,
//! capillary waves, turbulent-transition criterion or thermal evolution.
//! Midpoint is second order but can oscillate on unresolved diffusion modes;
//! unconditional linear stability is not a timestep-independent accuracy claim.

use fs_la::tridiagonal::{TridiagonalError, TridiagonalWorkspace};

use super::{FilmError, FilmMaterial, ShearedFilm};

mod carreau;
/// Original thermal storage coupled to resolved film momentum and local heat.
#[cfg(feature = "implicit-enthalpy-coupling")]
pub use carreau::thermal;
pub use carreau::{
    TransientCarreauControls, TransientCarreauError, TransientCarreauFilm,
    TransientCarreauFilmSpec, TransientCarreauStep,
};

/// Maximum retained wall-normal cells; both work and storage are linear in it.
pub const MAX_FILM_CELLS: usize = 65_536;

/// Signed tangential forcing, held over one interval (SI units).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct FilmDrive {
    /// Body acceleration, including any explicitly prescribed pressure gradient.
    pub acceleration_m_s2: f64,
    /// Translating wall velocity; normal wall motion is not represented.
    pub wall_velocity_m_s: f64,
    /// Traction applied to the liquid at the free surface, in Pa.
    pub surface_shear_pa: f64,
}

/// Explicit physical-speed and algebraic-solve admission limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FilmStepLimits {
    /// Bounds initial/final cell speeds, wall speed and the surface traces
    /// reconstructed with this interval's traction at both ends and midpoint.
    /// The caller qualifies the laminar regime; this is not a transition model.
    pub maximum_speed_m_s: f64,
    /// Maximum cell momentum defect divided by that cell's mass [m/s].
    pub maximum_velocity_residual_m_s: f64,
}

/// A refused step leaves accepted velocity and time unchanged.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum TransientFilmError {
    /// Native film material/geometry admission failure.
    Film(FilmError),
    /// A history, drive or numerical limit is invalid.
    Invalid(&'static str),
    /// Storage cannot be represented or allocated.
    Capacity,
    /// An intermediate or a required positive coefficient is unrepresentable.
    NonRepresentable,
    /// The prescribed physical-speed limit was exceeded.
    SpeedLimit,
    /// The original momentum equation exceeds its declared residual allowance.
    ResidualLimit,
    /// Cancellation observed before publication.
    Cancelled,
    /// Failure in the original native tridiagonal solve.
    Linear(TridiagonalError),
}
impl core::fmt::Display for TransientFilmError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "transient laminar film: {self:?}")
    }
}
impl std::error::Error for TransientFilmError {}
impl From<FilmError> for TransientFilmError {
    fn from(value: FilmError) -> Self {
        Self::Film(value)
    }
}
impl From<TridiagonalError> for TransientFilmError {
    fn from(value: TridiagonalError) -> Self {
        if value == TridiagonalError::Cancelled {
            Self::Cancelled
        } else {
            Self::Linear(value)
        }
    }
}

/// One accepted interval. Extensive quantities are per unit wetted wall area.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FilmStep {
    /// Duration of the accepted interval [s].
    pub duration_s: f64,
    /// Endpoint volume flux per width [m²/s].
    pub final_flow_per_width_m2_s: f64,
    /// Midpoint volume flux per width [m²/s].
    pub midpoint_flow_per_width_m2_s: f64,
    /// Velocity extrapolated across the surface half-cell using its traction.
    pub midpoint_surface_velocity_m_s: f64,
    /// Fluid-on-wall tangential traction; wall-on-fluid has the opposite sign.
    pub midpoint_wall_shear_pa: f64,
    /// Actual fluid-side wall trace, including the admitted Navier slip jump.
    pub midpoint_wall_fluid_velocity_m_s: f64,
    /// Change in integrated fluid momentum [kg/(m s)].
    pub momentum_change_kg_m_s: f64,
    /// Momentum change minus body and boundary impulses [kg/(m s)].
    pub momentum_residual_kg_m_s: f64,
    /// Evaluated as sum(m*u_mid*delta_u), avoiding subtraction of energies.
    pub kinetic_energy_change_j_m2: f64,
    /// Signed body-force work [J/m²].
    pub body_work_j_m2: f64,
    /// Signed work of the translating wall on the liquid [J/m²].
    pub wall_work_j_m2: f64,
    /// Signed work of the surface traction on the liquid [J/m²].
    pub surface_work_j_m2: f64,
    /// Nonnegative bulk viscous loss, including both fluid boundary half-cells.
    /// Excludes the separate wall/fluid interfacial loss.
    pub viscous_heat_j_m2: f64,
    /// Nonnegative interface loss tau_wall*(u_fluid_wall-U_wall)*dt [J/m²].
    /// No fraction is automatically deposited into either adjacent material.
    pub wall_slip_heat_j_m2: f64,
    /// Input work minus kinetic change, bulk heat and interface loss [J/m²].
    pub work_residual_j_m2: f64,
    /// Maximum original cell momentum residual divided by cell mass [m/s].
    pub maximum_velocity_residual_m_s: f64,
}

/// One original constant-property film, its resolved history and reusable O(n)
/// solve storage. Warm steps allocate nothing; equal durations reuse the native
/// factorization. Scratch may change on refusal, but accepted history does not.
#[derive(Debug)]
pub struct TransientFilm {
    material: FilmMaterial,
    thickness: f64,
    spacing: f64,
    wall_slip_length: f64,
    mass_per_cell_area: f64,
    shear_conductance: f64,
    diffusion_time: f64,
    time: f64,
    velocity: Vec<f64>,
    candidate: Vec<f64>,
    midpoint: Vec<f64>,
    delta: Vec<f64>,
    rhs: Vec<f64>,
    diagonal: Vec<f64>,
    adjacent: Vec<f64>,
    heat: Vec<f64>,
    linear: TridiagonalWorkspace,
    factored_duration: Option<f64>,
}

fn filled(n: usize) -> Result<Vec<f64>, TransientFilmError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(n)
        .map_err(|_| TransientFilmError::Capacity)?;
    values.resize(n, 0.0);
    Ok(values)
}
fn finite(value: f64) -> Result<f64, TransientFilmError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(TransientFilmError::NonRepresentable)
    }
}
fn positive(value: f64) -> Result<f64, TransientFilmError> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(TransientFilmError::NonRepresentable)
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), TransientFilmError> {
    if cancelled() {
        Err(TransientFilmError::Cancelled)
    } else {
        Ok(())
    }
}

// The dual shear segments tile the film. Each interior segment supplies equal
// halves to its neighbours; a boundary half-segment heats only its adjacent
// cell. Use the exact stresses from the accepted momentum work account.
fn stage_cell_heat(
    midpoint: &[f64],
    duration_s: f64,
    drive: FilmDrive,
    wall_bulk_difference: f64,
    surface_offset: f64,
    stress: impl Fn(usize) -> f64,
    heat: &mut [f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), TransientFilmError> {
    let n = midpoint.len();
    let segment = |traction: f64, difference: f64| -> Result<f64, TransientFilmError> {
        let value = finite(duration_s * traction * difference)?;
        if value < 0.0 {
            Err(TransientFilmError::NonRepresentable)
        } else {
            Ok(value)
        }
    };
    let mut left = segment(stress(0), wall_bulk_difference)?;
    for i in 0..n {
        if i % 128 == 0 {
            poll(cancelled)?;
        }
        let right = if i + 1 == n {
            segment(drive.surface_shear_pa, surface_offset)?
        } else {
            0.5 * segment(stress(i + 1), midpoint[i + 1] - midpoint[i])?
        };
        heat[i] = finite(left + right)?;
        left = right;
    }
    Ok(())
}

impl TransientFilm {
    /// Admit a fixed-thickness film initially at rest. The existing steady
    /// film owner supplies the same material/geometry admission rules.
    pub fn new(
        material: FilmMaterial,
        thickness_m: f64,
        cells: usize,
    ) -> Result<Self, TransientFilmError> {
        if !(2..=MAX_FILM_CELLS).contains(&cells) {
            return Err(TransientFilmError::Invalid("film cell count"));
        }
        Self::from_profile(material, thickness_m, filled(cells)?, 0.0)
    }

    /// Import actual cell-center velocities at a nonnegative physical time.
    /// This is an initial-condition operation, not a steady-profile reset.
    pub fn from_profile(
        material: FilmMaterial,
        thickness_m: f64,
        velocity_m_s: Vec<f64>,
        time_s: f64,
    ) -> Result<Self, TransientFilmError> {
        ShearedFilm::new(material, thickness_m, 0.0, 0.0)?;
        let n = velocity_m_s.len();
        if !(2..=MAX_FILM_CELLS).contains(&n)
            || !time_s.is_finite()
            || time_s < 0.0
            || velocity_m_s.iter().any(|u| !u.is_finite())
        {
            return Err(TransientFilmError::Invalid("film initial history"));
        }
        let spacing = positive(thickness_m / n as f64)?;
        let mass_per_cell_area = positive(material.density_kg_m3 * spacing)?;
        let shear_conductance =
            positive(material.density_kg_m3 * (material.kinematic_viscosity_m2_s / spacing))?;
        positive(2.0 * shear_conductance)?;
        let diffusion_time = positive(mass_per_cell_area / shear_conductance)?;
        Ok(Self {
            material,
            thickness: thickness_m,
            spacing,
            wall_slip_length: 0.0,
            mass_per_cell_area,
            shear_conductance,
            diffusion_time,
            time: time_s,
            velocity: velocity_m_s,
            candidate: filled(n)?,
            midpoint: filled(n)?,
            delta: filled(n)?,
            rhs: filled(n)?,
            diagonal: filled(n)?,
            adjacent: filled(n - 1)?,
            heat: filled(n)?,
            linear: TridiagonalWorkspace::new(n)?,
            factored_duration: None,
        })
    }

    /// Original constant SI material coefficients.
    pub const fn material(&self) -> FilmMaterial {
        self.material
    }
    /// Fixed wall-to-surface thickness.
    pub const fn thickness_m(&self) -> f64 {
        self.thickness
    }
    /// Uniform wall-normal control-volume width.
    pub const fn spacing_m(&self) -> f64 {
        self.spacing
    }
    /// Prescribed Navier length b in u_fluid_wall-U_wall = b*du/dy at the wall.
    /// Zero retains the original no-slip arithmetic. Infinite/perfect slip is
    /// not represented, and this value carries no material-card authority.
    pub const fn wall_slip_length_m(&self) -> f64 {
        self.wall_slip_length
    }
    /// Admit a new finite nonnegative slip length without replacing accepted
    /// velocity or time. The next interval uses it with its original wall
    /// speed and liquid viscosity. Refusal leaves the prior boundary intact.
    pub fn set_wall_slip_length_m(&mut self, length_m: f64) -> Result<(), TransientFilmError> {
        if !length_m.is_finite() || length_m < 0.0 {
            return Err(TransientFilmError::Invalid(
                "finite nonnegative wall slip length required",
            ));
        }
        let distance = positive(0.5 * self.spacing + length_m)?;
        positive(self.spacing / distance)?;
        positive(self.shear_conductance * (self.spacing / distance))?;
        self.wall_slip_length = length_m;
        self.factored_duration = None;
        Ok(())
    }

    fn wall_distance(&self) -> f64 {
        0.5 * self.spacing + self.wall_slip_length
    }

    fn wall_factor(&self) -> f64 {
        if self.wall_slip_length == 0.0 {
            2.0
        } else {
            self.spacing / self.wall_distance()
        }
    }

    // Preserve the no-slip subtraction exactly. For slip, each physical
    // distance multiplies the SAME original wall gradient; subtracting nearly
    // equal wall/fluid speeds must not erase a small bulk or interface loss.
    fn wall_differences(&self, difference: f64) -> Result<[f64; 2], TransientFilmError> {
        if self.wall_slip_length == 0.0 {
            Ok([difference, 0.0])
        } else {
            let rate = finite(difference / self.wall_distance())?;
            Ok([
                finite((0.5 * self.spacing) * rate)?,
                finite(self.wall_slip_length * rate)?,
            ])
        }
    }
    /// Physical time of the last fully admitted profile.
    pub const fn accepted_time_s(&self) -> f64 {
        self.time
    }
    /// Retained cell-center velocities, ordered from the wall to the surface.
    pub fn velocities_m_s(&self) -> &[f64] {
        &self.velocity
    }
    /// Cell-scale viscous diffusion time; resolve it for rapid forcing changes.
    pub fn diffusion_time_s(&self) -> f64 {
        self.diffusion_time
    }

    /// Advance the retained profile with one implicit-midpoint step. All
    /// numerical traversals poll at least every 128 cells. No accepted state
    /// changes until finite-value, speed, residual and cancellation checks pass.
    pub fn step(
        &mut self,
        duration_s: f64,
        drive: FilmDrive,
        limits: FilmStepLimits,
        cancelled: impl FnMut() -> bool,
    ) -> Result<FilmStep, TransientFilmError> {
        self.step_inner(duration_s, drive, limits, None, cancelled)
    }

    /// Advance and export this SAME interval's viscous heat into N original
    /// wall-normal cells [J/m² of wetted wall area]. Interior shear-segment work
    /// is split equally between adjacent cells; boundary half-segments heat
    /// only their adjacent cells. Interface slip loss is reported separately
    /// and never included in this output. Its sum is the native bulk loss,
    /// within summation roundoff. Supply the actual wall area when depositing
    /// into a caloric owner's joule inventory. Do not add mechanical work again.
    /// Failed/cancelled steps leave both accepted flow and output unchanged.
    pub fn step_with_heat(
        &mut self,
        duration_s: f64,
        drive: FilmDrive,
        limits: FilmStepLimits,
        cell_heat_j_m2: &mut [f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<FilmStep, TransientFilmError> {
        self.step_inner(duration_s, drive, limits, Some(cell_heat_j_m2), cancelled)
    }

    fn step_inner(
        &mut self,
        duration_s: f64,
        drive: FilmDrive,
        limits: FilmStepLimits,
        cell_heat_j_m2: Option<&mut [f64]>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FilmStep, TransientFilmError> {
        poll(&mut cancelled)?;
        if cell_heat_j_m2
            .as_ref()
            .is_some_and(|heat| heat.len() != self.velocity.len())
        {
            return Err(TransientFilmError::Invalid(
                "film needs N heat output cells",
            ));
        }
        if !duration_s.is_finite()
            || duration_s <= 0.0
            || !limits.maximum_speed_m_s.is_finite()
            || limits.maximum_speed_m_s <= 0.0
            || !limits.maximum_velocity_residual_m_s.is_finite()
            || limits.maximum_velocity_residual_m_s <= 0.0
            || ![
                drive.acceleration_m_s2,
                drive.wall_velocity_m_s,
                drive.surface_shear_pa,
            ]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(TransientFilmError::Invalid("film step inputs"));
        }
        let end_time = finite(self.time + duration_s)?;
        if end_time <= self.time {
            return Err(TransientFilmError::NonRepresentable);
        }
        if drive.wall_velocity_m_s.abs() > limits.maximum_speed_m_s {
            return Err(TransientFilmError::SpeedLimit);
        }
        let n = self.velocity.len();
        let weight = positive(duration_s / self.mass_per_cell_area)?;
        let alpha = positive(0.5 * weight * self.shear_conductance)?;
        let wall_factor = self.wall_factor();
        let wall_conductance = positive(wall_factor * self.shear_conductance)?;
        let body_increment = finite(duration_s * drive.acceleration_m_s2)?;
        if self.factored_duration != Some(duration_s) {
            self.factored_duration = None;
            for i in 0..n {
                if i % 128 == 0 {
                    poll(&mut cancelled)?;
                }
                // The wall half-cell and Navier length are in series. The top
                // face has specified stress, contributing no matrix unknown.
                let left = if i == 0 { wall_factor } else { 1.0 };
                let right = if i + 1 == n { 0.0 } else { 1.0 };
                self.diagonal[i] = finite(1.0 + alpha * (left + right))?;
                if i + 1 < n {
                    self.adjacent[i] = -alpha;
                }
            }
            self.linear.factor(
                &self.adjacent,
                &self.diagonal,
                &self.adjacent,
                &mut cancelled,
            )?;
            self.factored_duration = Some(duration_s);
        }
        for i in 0..n {
            if i % 128 == 0 {
                poll(&mut cancelled)?;
            }
            if self.velocity[i].abs() > limits.maximum_speed_m_s {
                return Err(TransientFilmError::SpeedLimit);
            }
            let left = if i == 0 {
                wall_conductance * (self.velocity[0] - drive.wall_velocity_m_s)
            } else {
                self.shear_conductance * (self.velocity[i] - self.velocity[i - 1])
            };
            let right = if i + 1 == n {
                drive.surface_shear_pa
            } else {
                self.shear_conductance * (self.velocity[i + 1] - self.velocity[i])
            };
            self.rhs[i] = finite(weight.mul_add(right - left, body_increment))?;
        }
        // Solve for the increment so small time steps do not subtract two
        // large equilibrium right-hand sides before resolving the dynamics.
        self.linear
            .solve_factored_into(&self.rhs, &mut self.delta, &mut cancelled)?;
        for i in 0..n {
            if i % 128 == 0 {
                poll(&mut cancelled)?;
            }
            self.candidate[i] = finite(self.velocity[i] + self.delta[i])?;
            // Audit the REPRESENTED update, including any rounding in addition.
            self.delta[i] = self.candidate[i] - self.velocity[i];
            self.midpoint[i] = finite(0.5f64.mul_add(self.delta[i], self.velocity[i]))?;
            if self.candidate[i].abs() > limits.maximum_speed_m_s {
                return Err(TransientFilmError::SpeedLimit);
            }
        }
        let surface_offset = finite(0.5 * drive.surface_shear_pa / self.shear_conductance)?;
        let receipt = self.report_candidate(
            duration_s,
            drive,
            limits,
            surface_offset,
            |face| {
                if face == 0 {
                    wall_conductance * (self.midpoint[0] - drive.wall_velocity_m_s)
                } else if face == n {
                    drive.surface_shear_pa
                } else {
                    self.shear_conductance * (self.midpoint[face] - self.midpoint[face - 1])
                }
            },
            &mut cancelled,
        )?;
        if cell_heat_j_m2.is_some() {
            let conductance = self.shear_conductance;
            let midpoint = &self.midpoint;
            let wall_bulk_difference =
                self.wall_differences(midpoint[0] - drive.wall_velocity_m_s)?[0];
            stage_cell_heat(
                midpoint,
                duration_s,
                drive,
                wall_bulk_difference,
                surface_offset,
                |face| {
                    if face == 0 {
                        wall_conductance * (midpoint[0] - drive.wall_velocity_m_s)
                    } else {
                        conductance * (midpoint[face] - midpoint[face - 1])
                    }
                },
                &mut self.heat,
                &mut cancelled,
            )?;
        }
        poll(&mut cancelled)?;
        if let Some(heat) = cell_heat_j_m2 {
            heat.copy_from_slice(&self.heat);
        }
        core::mem::swap(&mut self.velocity, &mut self.candidate);
        self.time = end_time;
        Ok(receipt)
    }

    // One finite-volume summation-by-parts account for Newtonian and nonlinear
    // films. The constitutive owner supplies its actual midpoint face stresses;
    // numerical defects are reported separately from the physical viscous heat.
    fn report_candidate(
        &self,
        duration_s: f64,
        drive: FilmDrive,
        limits: FilmStepLimits,
        surface_offset: f64,
        stress: impl Fn(usize) -> f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<FilmStep, TransientFilmError> {
        let n = self.velocity.len();
        let weight = positive(duration_s / self.mass_per_cell_area)?;
        let body_increment = finite(duration_s * drive.acceleration_m_s2)?;
        let wall_shear = finite(stress(0))?;
        let [wall_bulk_difference, wall_slip_velocity] =
            self.wall_differences(self.midpoint[0] - drive.wall_velocity_m_s)?;
        let wall_fluid_velocity = finite(drive.wall_velocity_m_s + wall_slip_velocity)?;
        let surface_velocity = finite(self.midpoint[n - 1] + surface_offset)?;
        for trace in [
            finite(self.velocity[n - 1] + surface_offset)?,
            surface_velocity,
            finite(self.candidate[n - 1] + surface_offset)?,
        ] {
            if trace.abs() > limits.maximum_speed_m_s {
                return Err(TransientFilmError::SpeedLimit);
            }
        }
        // With b>=0 the trace lies between the wall and first cell, but check
        // its represented value at all three profiles as part of speed admission.
        for velocity in [self.velocity[0], self.midpoint[0], self.candidate[0]] {
            let slip = self.wall_differences(velocity - drive.wall_velocity_m_s)?[1];
            if finite(drive.wall_velocity_m_s + slip)?.abs() > limits.maximum_speed_m_s {
                return Err(TransientFilmError::SpeedLimit);
            }
        }
        let mut q = 0.0;
        let mut q_final = 0.0;
        let mut momentum = 0.0;
        let mut kinetic = 0.0;
        let mut viscous_power =
            finite(wall_shear * wall_bulk_difference + drive.surface_shear_pa * surface_offset)?;
        let mut maximum_residual = 0.0_f64;
        for i in 0..n {
            if i % 128 == 0 {
                poll(cancelled)?;
            }
            let left = finite(stress(i))?;
            let right = finite(stress(i + 1))?;
            if i + 1 < n {
                let difference = self.midpoint[i + 1] - self.midpoint[i];
                viscous_power = finite(right.mul_add(difference, viscous_power))?;
            }
            let residual = finite(self.delta[i] - weight.mul_add(right - left, body_increment))?;
            maximum_residual = maximum_residual.max(residual.abs());
            q = finite(self.spacing.mul_add(self.midpoint[i], q))?;
            q_final = finite(self.spacing.mul_add(self.candidate[i], q_final))?;
            momentum = finite(self.mass_per_cell_area.mul_add(self.delta[i], momentum))?;
            kinetic = finite(
                (self.mass_per_cell_area * self.midpoint[i]).mul_add(self.delta[i], kinetic),
            )?;
        }
        if maximum_residual > limits.maximum_velocity_residual_m_s {
            return Err(TransientFilmError::ResidualLimit);
        }
        let body_work =
            finite(duration_s * self.material.density_kg_m3 * drive.acceleration_m_s2 * q)?;
        let wall_work = finite(-duration_s * wall_shear * drive.wall_velocity_m_s)?;
        let surface_work = finite(duration_s * drive.surface_shear_pa * surface_velocity)?;
        let heat = finite(duration_s * viscous_power)?;
        let slip_heat = if self.wall_slip_length == 0.0 {
            0.0
        } else {
            finite(duration_s * wall_shear * wall_slip_velocity)?
        };
        if slip_heat < 0.0 {
            return Err(TransientFilmError::NonRepresentable);
        }
        let receipt = FilmStep {
            duration_s,
            final_flow_per_width_m2_s: q_final,
            midpoint_flow_per_width_m2_s: q,
            midpoint_surface_velocity_m_s: surface_velocity,
            midpoint_wall_shear_pa: wall_shear,
            midpoint_wall_fluid_velocity_m_s: wall_fluid_velocity,
            momentum_change_kg_m_s: momentum,
            momentum_residual_kg_m_s: finite(
                momentum
                    - duration_s
                        * (self.material.density_kg_m3 * self.thickness * drive.acceleration_m_s2
                            + drive.surface_shear_pa
                            - wall_shear),
            )?,
            kinetic_energy_change_j_m2: kinetic,
            body_work_j_m2: body_work,
            wall_work_j_m2: wall_work,
            surface_work_j_m2: surface_work,
            viscous_heat_j_m2: heat,
            wall_slip_heat_j_m2: slip_heat,
            work_residual_j_m2: if self.wall_slip_length == 0.0 {
                finite(((body_work + wall_work) + surface_work) - kinetic - heat)?
            } else {
                finite(((body_work + wall_work) + surface_work) - kinetic - heat - slip_heat)?
            },
            maximum_velocity_residual_m_s: maximum_residual,
        };
        Ok(receipt)
    }
}

#[cfg(test)]
mod tests;
