//! Radially resolved startup and reversal of the original thermal Carreau law.
//!
//! The axially developed incompressible equation is
//! `rho u_t = delta_p/L + (r tau)_r/r`, `tau = mu(|u_r|) u_r`.
//! Annular finite volumes retain the entire radial velocity field. Implicit
//! midpoint uses the ORIGINAL local stress and positive stress tangent;
//! no mean-flow inertance or frozen apparent viscosity replaces radial diffusion.
//!
//! Summation by parts gives pressure work = kinetic change + positive viscous
//! heat, apart from the separately reported nonlinear/roundoff residual. Time
//! and radial discretization errors remain: midpoint is second order and can
//! oscillate on under-resolved stiff modes. Its energy balance is not a claim
//! of timestep-independent accuracy or positivity. Refine time and radius for
//! the actual forcing bandwidth; a radial diffusion timescale is exposed.
//!
//! Density, geometry and the original material calibration are fixed.
//! Temperature is uniform and held during each interval;
//! [`RadialCarreauPipe::at_temperature`] continues the original velocity
//! history under a new prescribed temperature.
//! Walls are rigid, stationary and no-slip. There is no axial development,
//! convective transport, compressibility, cavitation, transition prediction,
//! thermal evolution or junction/network model. Heat is reported, not deposited.

use core::f64::consts::PI;
use std::sync::Arc;

use fs_la::tridiagonal::{TridiagonalError, TridiagonalWorkspace};

use super::{LaminarCarreauPipe, MaterialError, positive};

/// Explicit spatial storage bound, independent of time and nonlinear budgets.
pub const MAX_RADIAL_CELLS: usize = 65_536;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RadialCarreauControls {
    /// Maximum cell momentum residual divided by that cell's mass [m/s].
    pub absolute_velocity_tolerance_m_s: f64,
    /// Relative to the largest initial or tentative axial velocity.
    pub relative_velocity_tolerance: f64,
    pub maximum_newton_iterations: u32,
    pub maximum_backtracks: u32,
    /// Total canonical local viscosity/tangent queries during one step.
    pub maximum_material_evaluations: u64,
}

impl Default for RadialCarreauControls {
    fn default() -> Self {
        Self {
            absolute_velocity_tolerance_m_s: 1e-12,
            relative_velocity_tolerance: 1e-10,
            maximum_newton_iterations: 64,
            maximum_backtracks: 24,
            maximum_material_evaluations: 1_000_000,
        }
    }
}

impl RadialCarreauControls {
    fn admit(self) -> Result<(), RadialCarreauError> {
        if !positive(self.absolute_velocity_tolerance_m_s)
            || !positive(self.relative_velocity_tolerance)
            || self.relative_velocity_tolerance >= 1.0
            || !(1..=128).contains(&self.maximum_newton_iterations)
            || self.maximum_backtracks > 64
            || self.maximum_material_evaluations == 0
        {
            return Err(RadialCarreauError::Input("invalid radial solve controls"));
        }
        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum RadialCarreauError {
    Input(&'static str),
    Material(MaterialError),
    Linear(TridiagonalError),
    Capacity,
    StateMismatch,
    Interrupted,
    MaterialBudget,
    NewtonBudget,
    BacktrackBudget,
    Unrepresentable,
}
impl core::fmt::Display for RadialCarreauError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "radial Carreau pipe: {self:?}")
    }
}
impl std::error::Error for RadialCarreauError {}
impl From<MaterialError> for RadialCarreauError {
    fn from(value: MaterialError) -> Self {
        Self::Material(value)
    }
}
impl From<TridiagonalError> for RadialCarreauError {
    fn from(value: TridiagonalError) -> Self {
        if value == TridiagonalError::Cancelled {
            Self::Interrupted
        } else {
            Self::Linear(value)
        }
    }
}

/// Uniform radial cells with exact annular areas and a prepared original law.
#[derive(Clone, Debug)]
pub struct RadialCarreauPipe {
    pipe: LaminarCarreauPipe,
    density_kg_m3: f64,
    cells: usize,
    spacing_m: f64,
    area_unit_m2: f64,
    mass_unit_kg: f64,
    diffusion_time_s: f64,
    identity: Arc<()>,
}

/// Accepted physical radial history; only a fully admitted step may replace it.
#[derive(Clone, Debug)]
pub struct RadialCarreauState {
    velocity_m_s: Vec<f64>,
    accepted_time_s: f64,
    identity: Arc<()>,
}
impl RadialCarreauState {
    pub fn velocities_m_s(&self) -> &[f64] {
        &self.velocity_m_s
    }
    pub const fn accepted_time_s(&self) -> f64 {
        self.accepted_time_s
    }
}

/// Reusable Newton, trial and native tridiagonal storage. Warm steps allocate
/// nothing; rejected trials may alter this private scratch, never the state.
#[derive(Debug)]
pub struct RadialCarreauWorkspace {
    candidate: Vec<f64>,
    trial: Vec<f64>,
    residual: Vec<f64>,
    diagonal: Vec<f64>,
    off_diagonal: Vec<f64>,
    correction: Vec<f64>,
    linear: TridiagonalWorkspace,
    identity: Arc<()>,
}

/// Physical quantities use the same accepted midpoint field and face stresses.
/// Work defects and finite discretization errors are not additional heat.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RadialCarreauStep {
    pub duration_s: f64,
    pub pressure_drop_pa: f64,
    pub initial_volume_flow_m3_s: f64,
    pub final_volume_flow_m3_s: f64,
    /// Midpoint quadrature of the time-average flow, not an exact trajectory integral.
    pub midpoint_volume_flow_m3_s: f64,
    pub displaced_volume_m3: f64,
    pub initial_kinetic_energy_j: f64,
    pub final_kinetic_energy_j: f64,
    /// Evaluated as sum(m*u_mid*delta_u), avoiding subtraction of large energies.
    pub kinetic_energy_change_j: f64,
    pub pressure_work_j: f64,
    /// Positive sum of original tau*gamma at the physical radial faces.
    pub viscous_heat_j: f64,
    pub work_residual_j: f64,
    /// Negative sum(u_mid * original momentum residual), in joules.
    pub momentum_residual_work_j: f64,
    /// Original sigma_rz at the wall, signed in the increasing-r direction.
    pub midpoint_wall_stress_pa: f64,
    /// Maximum face rate over the numerical linear-in-time trajectory.
    pub maximum_shear_rate_per_s: f64,
    pub max_velocity_residual_m_s: f64,
    /// d(midpoint Q)/d(delta_p), at fixed initial history and duration.
    pub midpoint_pressure_conductance_m3_s_pa: f64,
    pub newton_iterations: u32,
    pub backtracks: u32,
    pub material_evaluations: u64,
}

#[derive(Clone, Copy)]
struct Evaluation {
    norm: f64,
    scale: f64,
    dissipation_w: f64,
    wall_stress_pa: f64,
    maximum_rate: f64,
}

fn filled(n: usize) -> Result<Vec<f64>, RadialCarreauError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(n)
        .map_err(|_| RadialCarreauError::Capacity)?;
    values.resize(n, 0.0);
    Ok(values)
}
fn finite(x: f64) -> Result<f64, RadialCarreauError> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(RadialCarreauError::Unrepresentable)
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), RadialCarreauError> {
    if cancelled() {
        Err(RadialCarreauError::Interrupted)
    } else {
        Ok(())
    }
}
fn tile(index: usize, cancelled: &mut impl FnMut() -> bool) -> Result<(), RadialCarreauError> {
    if index % 64 == 0 {
        poll(cancelled)?;
    }
    Ok(())
}

// Positive and signed physical sums share a fixed compensated order.
#[derive(Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, x: f64) -> Result<(), RadialCarreauError> {
        let next = finite(self.value + x)?;
        self.correction = finite(
            self.correction
                + if self.value.abs() >= x.abs() {
                    (self.value - next) + x
                } else {
                    (x - next) + self.value
                },
        )?;
        self.value = next;
        Ok(())
    }
    fn finish(self) -> Result<f64, RadialCarreauError> {
        finite(self.value + self.correction)
    }
}

impl RadialCarreauPipe {
    pub fn new(
        pipe: LaminarCarreauPipe,
        density_kg_m3: f64,
        radial_cells: usize,
    ) -> Result<Self, RadialCarreauError> {
        if !positive(density_kg_m3) || !(2..=MAX_RADIAL_CELLS).contains(&radial_cells) {
            return Err(RadialCarreauError::Input(
                "invalid density or radial cell count",
            ));
        }
        let spacing_m = pipe.radius_m / radial_cells as f64;
        let area_unit_m2 = (PI * spacing_m) * spacing_m;
        let mass_unit_kg = (density_kg_m3 * area_unit_m2) * pipe.length_m;
        let diffusion_time_s = super::product_ratio(
            &[density_kg_m3, spacing_m, spacing_m],
            &[pipe.prepared.rest.viscosity],
        )?;
        if ![spacing_m, area_unit_m2, mass_unit_kg, 1.0 / spacing_m]
            .into_iter()
            .all(positive)
            || !positive(mass_unit_kg * (2 * radial_cells - 1) as f64)
        {
            return Err(RadialCarreauError::Unrepresentable);
        }
        Ok(Self {
            pipe,
            density_kg_m3,
            cells: radial_cells,
            spacing_m,
            area_unit_m2,
            mass_unit_kg,
            diffusion_time_s,
            identity: Arc::new(()),
        })
    }

    pub const fn pipe(&self) -> &LaminarCarreauPipe {
        &self.pipe
    }
    /// Continue the original radial history at a new prescribed temperature.
    ///
    /// The original thermal law admits the complete physical rate envelope at
    /// this temperature. Density, geometry, annular masses and every accepted
    /// velocity remain unchanged, so changing viscosity creates no kinetic
    /// energy or instantaneous momentum jump. Existing states, output buffers
    /// and workspaces remain compatible; their identity denotes this source
    /// calibration and spatial discretization, not a temperature fingerprint.
    /// Only the prepared constitutive response and diffusion timescale change.
    /// This operation allocates no radial storage and never mutates history.
    ///
    /// The caller supplies a spatially uniform temperature held over the next
    /// interval. The original fixed-density assumption must remain valid;
    /// thermal expansion, a density change, thermal storage and the energy used
    /// to change temperature require separate physical owners. This method
    /// neither evolves temperature nor deposits the reported viscous heat.
    pub fn at_temperature(&self, temperature_k: f64) -> Result<Self, RadialCarreauError> {
        let pipe = self.pipe.at_temperature(temperature_k)?;
        let diffusion_time_s = super::product_ratio(
            &[self.density_kg_m3, self.spacing_m, self.spacing_m],
            &[pipe.prepared.rest.viscosity],
        )?;
        Ok(Self {
            pipe,
            diffusion_time_s,
            ..self.clone()
        })
    }
    pub const fn density_kg_m3(&self) -> f64 {
        self.density_kg_m3
    }
    pub const fn radial_cells(&self) -> usize {
        self.cells
    }
    pub const fn spacing_m(&self) -> f64 {
        self.spacing_m
    }
    pub fn cell_radius_m(&self, index: usize) -> Option<f64> {
        (index < self.cells).then(|| (index as f64 + 0.5) * self.spacing_m)
    }
    /// rho*dr²/mu0 for resolving the fastest radial diffusion. This is a
    /// discretization timescale, not a time-error or transition certificate.
    pub fn radial_diffusion_time_s(&self) -> f64 {
        self.diffusion_time_s
    }

    pub fn workspace(&self) -> Result<RadialCarreauWorkspace, RadialCarreauError> {
        Ok(RadialCarreauWorkspace {
            candidate: filled(self.cells)?,
            trial: filled(self.cells)?,
            residual: filled(self.cells)?,
            diagonal: filled(self.cells)?,
            off_diagonal: filled(self.cells - 1)?,
            correction: filled(self.cells)?,
            linear: TridiagonalWorkspace::new(self.cells)?,
            identity: self.identity.clone(),
        })
    }

    pub fn rest_state(&self) -> Result<RadialCarreauState, RadialCarreauError> {
        Ok(RadialCarreauState {
            velocity_m_s: filled(self.cells)?,
            accepted_time_s: 0.0,
            identity: self.identity.clone(),
        })
    }

    /// Import an actual resolved radial history, with cell-center velocities.
    /// Source support is checked at every interior face and the no-slip wall.
    pub fn state_from_velocities(
        &self,
        velocities_m_s: Vec<f64>,
        accepted_time_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RadialCarreauState, RadialCarreauError> {
        poll(&mut cancelled)?;
        if velocities_m_s.len() != self.cells
            || !accepted_time_s.is_finite()
            || accepted_time_s < 0.0
        {
            return Err(RadialCarreauError::Input("invalid radial history"));
        }
        self.admit_velocities(&velocities_m_s, &mut cancelled)?;
        for face in 0..self.cells {
            tile(face, &mut cancelled)?;
            self.pipe
                .prepared
                .sample_rate(self.gradient(&velocities_m_s, face)?.abs())?;
        }
        poll(&mut cancelled)?;
        Ok(RadialCarreauState {
            velocity_m_s: velocities_m_s,
            accepted_time_s,
            identity: self.identity.clone(),
        })
    }

    fn gradient(&self, values: &[f64], face: usize) -> Result<f64, RadialCarreauError> {
        finite(if face + 1 < self.cells {
            (values[face + 1] - values[face]) / self.spacing_m
        } else {
            -values[face] / (0.5 * self.spacing_m)
        })
    }

    fn admit_velocities(
        &self,
        values: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, RadialCarreauError> {
        let mut maximum = 0.0_f64;
        for face in 0..self.cells {
            tile(face, cancelled)?;
            finite(values[face])?;
            let rate = self.gradient(values, face)?.abs();
            if rate > self.pipe.prepared.ceiling.rate {
                return Err(RadialCarreauError::Material(super::state_error(
                    "radial Carreau shear leaves the original source/rate envelope",
                )));
            }
            maximum = maximum.max(rate);
        }
        Ok(maximum)
    }

    #[allow(clippy::too_many_arguments)]
    fn evaluate(
        &self,
        initial: &[f64],
        candidate: &[f64],
        pressure: f64,
        dt: f64,
        controls: RadialCarreauControls,
        evaluations: &mut u64,
        residual: &mut [f64],
        diagonal: &mut [f64],
        off: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Evaluation, RadialCarreauError> {
        // Both endpoints are admitted; their linear interpolation cannot have
        // a larger face-rate magnitude. No unsupported endpoint is published.
        let maximum_rate = self.admit_velocities(candidate, cancelled)?;
        let mut scale = 0.0_f64;
        for i in 0..self.cells {
            tile(i, cancelled)?;
            let ring = (2 * i + 1) as f64;
            let mass = self.mass_unit_kg * ring;
            residual[i] = finite(
                mass * (candidate[i] - initial[i]) - ((self.area_unit_m2 * ring) * pressure) * dt,
            )?;
            diagonal[i] = mass;
            scale = scale.max(initial[i].abs()).max(candidate[i].abs());
        }
        let mut heat = Sum::default();
        let mut wall_stress = 0.0;
        for face in 0..self.cells {
            tile(face, cancelled)?;
            if *evaluations >= controls.maximum_material_evaluations {
                return Err(RadialCarreauError::MaterialBudget);
            }
            *evaluations += 1;
            let gradient =
                0.5 * self.gradient(initial, face)? + 0.5 * self.gradient(candidate, face)?;
            let material = self.pipe.prepared.sample_rate(gradient.abs())?;
            let stress = finite(material.viscosity * gradient)?;
            let radius = (face + 1) as f64 * self.spacing_m;
            let surface = (2.0 * PI * radius) * self.pipe.length_m;
            let distance = if face + 1 == self.cells {
                0.5 * self.spacing_m
            } else {
                self.spacing_m
            };
            let impulse = finite((surface * stress) * dt)?;
            let tangent = finite((0.5 * dt * surface) * (material.tangent / distance))?;
            residual[face] = finite(residual[face] - impulse)?;
            diagonal[face] = finite(diagonal[face] + tangent)?;
            if face + 1 < self.cells {
                residual[face + 1] = finite(residual[face + 1] + impulse)?;
                diagonal[face + 1] = finite(diagonal[face + 1] + tangent)?;
                off[face] = -tangent;
            } else {
                wall_stress = stress;
            }
            heat.add(((surface * distance) * stress) * gradient)?;
        }
        let mut norm = 0.0_f64;
        for (i, &r) in residual.iter().enumerate() {
            tile(i, cancelled)?;
            norm = norm.max(finite(r / (self.mass_unit_kg * (2 * i + 1) as f64))?.abs());
        }
        let dissipation_w = heat.finish()?;
        if dissipation_w < 0.0 {
            return Err(RadialCarreauError::Unrepresentable);
        }
        Ok(Evaluation {
            norm,
            scale,
            dissipation_w,
            wall_stress_pa: wall_stress,
            maximum_rate,
        })
    }

    /// Advance one held-pressure interval with a convergent implicit-midpoint
    /// radial discretization. Construction allocates; warm stepping does not.
    /// Cancellation, source refusal, budgets or unresolved arithmetic preserve
    /// accepted history/time exactly, and the same workspace remains reusable.
    #[allow(clippy::too_many_arguments)]
    pub fn step(
        &self,
        state: &mut RadialCarreauState,
        workspace: &mut RadialCarreauWorkspace,
        pressure_drop_pa: f64,
        duration_s: f64,
        controls: RadialCarreauControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RadialCarreauStep, RadialCarreauError> {
        let report = self.prepare_step(
            state,
            workspace,
            pressure_drop_pa,
            duration_s,
            controls,
            &mut cancelled,
        )?;
        core::mem::swap(&mut state.velocity_m_s, &mut workspace.candidate);
        state.accepted_time_s += duration_s;
        Ok(report)
    }

    /// Resolve a candidate from immutable initial history into a reusable
    /// output state. Network pressure trials can share the original radial
    /// solve without advancing accepted branches or allocating new profiles.
    /// Both states must belong to this model. Refusal or cancellation preserves
    /// the output exactly; success replaces its whole history and time.
    #[allow(clippy::too_many_arguments)]
    pub fn step_from(
        &self,
        initial: &RadialCarreauState,
        output: &mut RadialCarreauState,
        workspace: &mut RadialCarreauWorkspace,
        pressure_drop_pa: f64,
        duration_s: f64,
        controls: RadialCarreauControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RadialCarreauStep, RadialCarreauError> {
        if !Arc::ptr_eq(&self.identity, &output.identity) {
            return Err(RadialCarreauError::StateMismatch);
        }
        let report = self.prepare_step(
            initial,
            workspace,
            pressure_drop_pa,
            duration_s,
            controls,
            &mut cancelled,
        )?;
        core::mem::swap(&mut output.velocity_m_s, &mut workspace.candidate);
        output.accepted_time_s = initial.accepted_time_s + duration_s;
        Ok(report)
    }

    #[allow(clippy::too_many_arguments)]
    fn prepare_step(
        &self,
        state: &RadialCarreauState,
        workspace: &mut RadialCarreauWorkspace,
        pressure_drop_pa: f64,
        duration_s: f64,
        controls: RadialCarreauControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RadialCarreauStep, RadialCarreauError> {
        poll(&mut cancelled)?;
        controls.admit()?;
        if !Arc::ptr_eq(&self.identity, &state.identity)
            || !Arc::ptr_eq(&self.identity, &workspace.identity)
        {
            return Err(RadialCarreauError::StateMismatch);
        }
        if !pressure_drop_pa.is_finite() || !positive(duration_s) {
            return Err(RadialCarreauError::Input("invalid pressure or duration"));
        }
        let next_time = finite(state.accepted_time_s + duration_s)?;
        if next_time <= state.accepted_time_s {
            return Err(RadialCarreauError::Unrepresentable);
        }
        let initial_rate = self.admit_velocities(&state.velocity_m_s, &mut cancelled)?;
        for (i, &value) in state.velocity_m_s.iter().enumerate() {
            tile(i, &mut cancelled)?;
            workspace.candidate[i] = value;
        }
        let (mut evaluations, mut iterations, mut backtracks) = (0, 0, 0);
        let mut accepted = self.evaluate(
            &state.velocity_m_s,
            &workspace.candidate,
            pressure_drop_pa,
            duration_s,
            controls,
            &mut evaluations,
            &mut workspace.residual,
            &mut workspace.diagonal,
            &mut workspace.off_diagonal,
            &mut cancelled,
        )?;
        loop {
            poll(&mut cancelled)?;
            let tolerance = finite(
                controls.absolute_velocity_tolerance_m_s
                    + controls.relative_velocity_tolerance * accepted.scale,
            )?;
            if accepted.norm <= tolerance {
                break;
            }
            if iterations >= controls.maximum_newton_iterations {
                return Err(RadialCarreauError::NewtonBudget);
            }
            for (i, r) in workspace.residual.iter_mut().enumerate() {
                tile(i, &mut cancelled)?;
                *r = -*r;
            }
            workspace.linear.solve_into(
                &workspace.off_diagonal,
                &workspace.diagonal,
                &workspace.off_diagonal,
                &workspace.residual,
                &mut workspace.correction,
                &mut cancelled,
            )?;
            let mut step = 1.0_f64;
            let mut found = None;
            for reduction in 0..=controls.maximum_backtracks {
                poll(&mut cancelled)?;
                for i in 0..self.cells {
                    tile(i, &mut cancelled)?;
                    workspace.trial[i] =
                        step.mul_add(workspace.correction[i], workspace.candidate[i]);
                }
                match self.evaluate(
                    &state.velocity_m_s,
                    &workspace.trial,
                    pressure_drop_pa,
                    duration_s,
                    controls,
                    &mut evaluations,
                    &mut workspace.residual,
                    &mut workspace.diagonal,
                    &mut workspace.off_diagonal,
                    &mut cancelled,
                ) {
                    Ok(value) => {
                        let bound = controls.absolute_velocity_tolerance_m_s
                            + controls.relative_velocity_tolerance * value.scale;
                        if value.norm < accepted.norm || value.norm <= bound {
                            found = Some(value);
                            break;
                        }
                    }
                    Err(
                        error @ (RadialCarreauError::Interrupted
                        | RadialCarreauError::MaterialBudget),
                    ) => return Err(error),
                    Err(_) => {}
                }
                if reduction < controls.maximum_backtracks {
                    step *= 0.5;
                    backtracks += 1;
                }
            }
            accepted = found.ok_or(RadialCarreauError::BacktrackBudget)?;
            core::mem::swap(&mut workspace.candidate, &mut workspace.trial);
            iterations += 1;
        }
        // Differentiate the accepted original midpoint equations. The same
        // positive tangent matrix gives a pressure-conjugate flow response.
        for (i, rhs) in workspace.correction.iter_mut().enumerate() {
            tile(i, &mut cancelled)?;
            *rhs = finite(duration_s * self.area_unit_m2 * (2 * i + 1) as f64)?;
        }
        workspace.linear.factor(
            &workspace.off_diagonal,
            &workspace.diagonal,
            &workspace.off_diagonal,
            &mut cancelled,
        )?;
        workspace
            .linear
            .solve_factored_in_place(&mut workspace.correction, &mut cancelled)?;
        let (mut initial_flow, mut final_flow, mut mid_flow) =
            (Sum::default(), Sum::default(), Sum::default());
        let (mut initial_energy, mut final_energy, mut change) =
            (Sum::default(), Sum::default(), Sum::default());
        let (mut residual_work, mut conductance) = (Sum::default(), Sum::default());
        for i in 0..self.cells {
            tile(i, &mut cancelled)?;
            let old = state.velocity_m_s[i];
            let new = workspace.candidate[i];
            let mid = 0.5 * old + 0.5 * new;
            let area = self.area_unit_m2 * (2 * i + 1) as f64;
            let mass = self.mass_unit_kg * (2 * i + 1) as f64;
            initial_flow.add(area * old)?;
            final_flow.add(area * new)?;
            mid_flow.add(area * mid)?;
            initial_energy.add((0.5 * mass * old) * old)?;
            final_energy.add((0.5 * mass * new) * new)?;
            change.add((mass * mid) * (new - old))?;
            residual_work.add(-mid * workspace.residual[i])?;
            conductance.add(0.5 * area * workspace.correction[i])?;
        }
        let midpoint_volume_flow_m3_s = mid_flow.finish()?;
        let displaced_volume_m3 = finite(duration_s * midpoint_volume_flow_m3_s)?;
        let pressure_work_j = finite(pressure_drop_pa * displaced_volume_m3)?;
        let kinetic_energy_change_j = change.finish()?;
        let viscous_heat_j = finite(duration_s * accepted.dissipation_w)?;
        let work_residual_j = finite(pressure_work_j - kinetic_energy_change_j - viscous_heat_j)?;
        let midpoint_pressure_conductance_m3_s_pa = conductance.finish()?;
        if !positive(midpoint_pressure_conductance_m3_s_pa) {
            return Err(RadialCarreauError::Unrepresentable);
        }
        let report = RadialCarreauStep {
            duration_s,
            pressure_drop_pa,
            initial_volume_flow_m3_s: initial_flow.finish()?,
            final_volume_flow_m3_s: final_flow.finish()?,
            midpoint_volume_flow_m3_s,
            displaced_volume_m3,
            initial_kinetic_energy_j: initial_energy.finish()?,
            final_kinetic_energy_j: final_energy.finish()?,
            kinetic_energy_change_j,
            pressure_work_j,
            viscous_heat_j,
            work_residual_j,
            momentum_residual_work_j: residual_work.finish()?,
            midpoint_wall_stress_pa: accepted.wall_stress_pa,
            maximum_shear_rate_per_s: initial_rate.max(accepted.maximum_rate),
            max_velocity_residual_m_s: accepted.norm,
            midpoint_pressure_conductance_m3_s_pa,
            newton_iterations: iterations,
            backtracks,
            material_evaluations: evaluations,
        };
        poll(&mut cancelled)?;
        Ok(report)
    }
}

#[cfg(test)]
mod tests;
