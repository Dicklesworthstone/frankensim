//! Unsteady incompressible flow, and transient conjugate heat transfer on
//! it, with the finite-volume SIMPLEC discretization of [`super::simple`].
//!
//! # Time discretization
//!
//! The momentum equations gain `rho V du/dt` over each staggered volume,
//! implicit by backward Euler or BDF2 ([`TimeScheme`]; BDF2's first step is
//! backward Euler):
//!
//! ```text
//! BDF2:  rho V (3 u^(n+1) - 4 u^n + u^(n-1)) / (2 dt) + (steady terms at n+1) = 0
//! ```
//!
//! Every step iterates SIMPLEC sweeps to the step's own steady residuals
//! (`inner_tolerance` on the mass and momentum residuals of the
//! time-discrete equations; [`super::simple_flow`] defines them), then
//! removes the remaining divergence by one tight pressure correction, so
//! each step's fluxes are divergence-free to 1e-8 of the converged
//! imbalance. A step that does not converge within `inner_iterations`
//! refuses (it is never silently accepted). The implicit scheme has no
//! CFL stability limit; the step size is the caller's accuracy decision,
//! and halving it is the convergence evidence.
//!
//! [`march_conjugate`] advances flow and energy together: each step solves
//! the flow at `t^(n+1)` (with the Boussinesq force of `T^n` when declared,
//! a first-order lag), then the backward-Euler energy step of
//! [`super::march_energy`] on that step's fluxes, with its energy closure.
//!
//! # No-claim boundaries
//!
//! The SIMPLEC no-claims apply (staircase walls, power-law convection). An
//! LVEL transient remains an equilibrium eddy-viscosity solution. The
//! explicit `*_with_sgs` entries instead compute Smagorinsky/WALE from each
//! iterate's full velocity gradient, for exploration only; they do not
//! establish resolved turbulent spectra or thermal wall resolution. The
//! Energy coupling is first order in time. The native implicit-buoyancy and
//! material endpoint entries close force and transport on the returned
//! temperature; ordinary marches retain
//! the explicitly documented one-step source lag.
//! Face fans are steady models and refuse here; use
//! internal fans. The opt-in thermal_wall::march_conjugate_with_wall_heat
//! entry uses this same march with freshly prepared equilibrium wall transport.

use core::cell::Cell;
use fs_exec::CancelGate;

mod endpoint;
mod native_buoyancy;
mod pressure_program;
pub(crate) use endpoint::MaterialEndpointCoupling;
pub use endpoint::{MaterialEndpointConfig, MaterialEndpointMarch, MaterialEndpointStepReport};
pub use native_buoyancy::{
    march_conjugate_implicit_buoyancy, march_conjugate_with_sgs_implicit_buoyancy,
};
pub use pressure_program::{PressureOffsetProgram, PressureProgram};

#[cfg(feature = "enthalpy-coupling")]
mod material;
pub(in crate::conjugate) mod wall_heat;
#[cfg(feature = "enthalpy-coupling")]
pub(super) use material::MaterialMomentum;

use super::domain::{FluidProperties, SolidMaterial, VoxelDomain};
pub(crate) use super::energy::compact_owners;
use super::energy::{AdvectiveEnthalpyField, EnergyConfig, EnergySolution, ThermalSetup};
use super::flow::FlowField;
use super::radiation::{RadiationConfig, Radiators};
use super::simple::sgs::FvSgsConfig;
use super::simple::{FvBoundary, FvFlow, SimpleConfig, Solver, TimeScheme, Turbulence, admit};
use super::transient::{
    TransientRecord, energy_step_with_conductivity, heat_capacity, junction_stored_energy_change_j,
    record_energy_step, scale_power, solve_step_with_fields, storage_record,
};
use super::{ChtError, finite, finite_positive, poll};
#[cfg(feature = "enthalpy-coupling")]
use super::{UnsteadyHeatingAudit, UnsteadyMechanicalStep, ViscousHeatControls};

/// Optional physical-work account on the original source-caloric clock.
/// It remains private until the complete thermal/source march succeeds.
#[cfg(feature = "enthalpy-coupling")]
pub(crate) struct MaterialViscousHeat {
    pub controls: ViscousHeatControls,
    pub audit: UnsteadyHeatingAudit,
    pub records: Vec<UnsteadyMechanicalStep>,
}

#[cfg(feature = "enthalpy-coupling")]
impl MaterialViscousHeat {
    pub(crate) fn new(controls: ViscousHeatControls) -> Self {
        Self {
            controls,
            audit: UnsteadyHeatingAudit::default(),
            records: Vec::new(),
        }
    }

    /// Replace only the physical deposit while retaining the current authored
    /// power schedule. Fixed-cell heat has a reciprocal thermostat rejection.
    fn deposit(
        thermal: &mut ThermalSetup,
        authored_power_w: &[f64],
        receipt: &UnsteadyMechanicalStep,
        gate: &CancelGate,
    ) -> Result<f64, ChtError> {
        let cells = receipt.power.cell_heat_w.len();
        if !authored_power_w.is_empty() && authored_power_w.len() != cells {
            return Err(ChtError::InvalidInput {
                field: "thermal.power_w",
                reason: "expected one source per original cell".into(),
            });
        }
        thermal.power_w.clear();
        thermal.power_w.extend_from_slice(authored_power_w);
        if thermal.power_w.is_empty() {
            thermal.power_w.resize(cells, 0.0);
        }
        for (cell, (source, &heat)) in thermal.power_w.iter_mut()
            .zip(&receipt.power.cell_heat_w).enumerate()
        {
            if cell.is_multiple_of(256) { poll(gate)?; }
            *source += heat;
            finite("thermal.viscous_heat.power_w", *source)?;
        }
        let mut fixed: Vec<_> = thermal.fixed_temperature.iter()
            .map(|&(cell, _)| cell).collect();
        fixed.sort_unstable();
        fixed.dedup();
        let mut fixed_heat_w = 0.0;
        for cell in fixed {
            fixed_heat_w += *receipt.power.cell_heat_w.get(cell).ok_or_else(|| {
                ChtError::InvalidInput {
                    field: "thermal.fixed_temperature",
                    reason: "fixed cell outside original domain".into(),
                }
            })?;
        }
        let fixed_heat_j = fixed_heat_w * receipt.duration_s;
        finite("thermal.viscous_heat.fixed_energy_j", fixed_heat_j)?;
        Ok(fixed_heat_j)
    }
}

/// Time-march controls.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UnsteadyConfig {
    /// Time step, s.
    pub time_step_s: f64,
    /// Steps.
    pub steps: usize,
    /// Momentum time discretization.
    pub scheme: TimeScheme,
    /// SIMPLEC sweeps allowed per step.
    pub inner_iterations: usize,
    /// Mass and momentum residual each step must reach.
    pub inner_tolerance: f64,
    /// First step (1-based) of the time-averaging window.
    pub average_from_step: usize,
    /// Cell whose velocity every step records.
    pub probe: Option<usize>,
    /// Differential reservoir-pressure factor sampled at each step endpoint.
    /// This retained program reaches ordinary, SGS, wall-heat and material
    /// conjugate marches, including implicit endpoint feedback. A flow-only
    /// pressure callback, when supplied, multiplies this factor.
    pub pressure_program: PressureProgram,
    /// Optional signed-Pa offset programs on the six boundary planes in
    /// `Face3::ALL` order. Each selected plane must be Pressure/Outlet.
    /// Offsets add after the scalar factor and preserve authored profiles.
    pub opening_pressure_programs: [Option<PressureOffsetProgram>; 6],
}

impl UnsteadyConfig {
    /// BDF2 defaults: 100 sweeps per step to 1e-6, averaging the second
    /// half of the march.
    #[must_use]
    pub const fn new(time_step_s: f64, steps: usize) -> Self {
        Self {
            time_step_s,
            steps,
            scheme: TimeScheme::Bdf2,
            inner_iterations: 100,
            inner_tolerance: 1e-6,
            average_from_step: steps / 2 + 1,
            probe: None,
            pressure_program: PressureProgram::IDENTITY,
            opening_pressure_programs: [None; 6],
        }
    }

    fn opening_pressure_offsets_pa(&self, time_s: f64) -> Result<[f64; 6], ChtError> {
        let mut offsets = [0.0; 6];
        for (side, program) in self.opening_pressure_programs.iter().enumerate() {
            if let Some(program) = program {
                offsets[side] = program.sample_pa(time_s)?;
            }
        }
        Ok(offsets)
    }
}

/// Evidence of one flow step.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FlowStepRecord {
    /// Time at the end of the step, s.
    pub time_s: f64,
    /// SIMPLEC sweeps the step took.
    pub inner_iterations: usize,
    /// Mass residual at convergence of the step.
    pub mass_residual: f64,
    /// Momentum residual at convergence of the step.
    pub momentum_residual: f64,
    /// Velocity of the probe cell, m/s (zeros without a probe).
    pub probe_velocity_m_s: [f64; 3],
    /// Staggered kinetic energy `sum m_face * u_face^2 / 2`, J. Each face
    /// owns the adjacent fluid half cells; exposed faces own only the inner
    /// half cell. Cell-centred velocity averaging is not applied to this norm.
    pub kinetic_energy_j: f64,
}

/// An unsteady flow march.
#[derive(Debug, Clone, PartialEq)]
pub struct UnsteadyFlow {
    /// The flow at the final time (its report covers the whole march).
    pub flow: FvFlow,
    /// One record per step.
    pub records: Vec<FlowStepRecord>,
    /// Cell velocity averaged over the steps from `average_from_step`, m/s.
    pub mean_velocity_m_s: Vec<[f64; 3]>,
    /// Face fluxes averaged over the same steps: divergence-free like every
    /// step's (an average of projected fields), so the steady energy
    /// equation may run on it (neglecting the unsteady correlation
    /// `<u' T'>`).
    pub mean_field: FlowField,
    /// Steps in the average.
    pub averaged_steps: usize,
}

/// A coupled flow and energy march.
#[derive(Debug, Clone, PartialEq)]
pub struct ConjugateMarch {
    /// The flow march (final flow, flow records, mean velocity).
    pub flow: UnsteadyFlow,
    /// Temperature per cell at the final time, K.
    pub temperature: Vec<f64>,
    /// Cell temperature averaged over the same window, K.
    pub mean_temperature: Vec<f64>,
    /// One energy record per step (closure, peak temperatures).
    pub energy_records: Vec<TransientRecord>,
}

/// Boussinesq buoyancy for [`march_conjugate`] and
/// [`march_conjugate_implicit_buoyancy`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Boussinesq {
    /// Gravitational acceleration, m/s^2.
    pub gravity_m_s2: [f64; 3],
    /// Volumetric expansion coefficient, 1/K.
    pub expansion_per_k: f64,
    /// Reference (ambient) temperature, K.
    pub reference_temperature_k: f64,
}

fn admit_unsteady(
    domain: &VoxelDomain,
    config: &SimpleConfig,
    unsteady: &UnsteadyConfig,
) -> Result<(), ChtError> {
    finite_positive("unsteady.time_step_s", unsteady.time_step_s)?;
    finite_positive("unsteady.inner_tolerance", unsteady.inner_tolerance)?;
    let refuse = |field: &'static str, reason: &str| {
        Err(ChtError::InvalidInput {
            field,
            reason: reason.into(),
        })
    };
    if unsteady.steps == 0 || unsteady.inner_iterations == 0 {
        return refuse(
            "unsteady.steps",
            "steps and inner_iterations must be at least one",
        );
    }
    finite_positive(
        "unsteady.end_time_s",
        unsteady.steps as f64 * unsteady.time_step_s,
    )?;
    if unsteady.average_from_step == 0 || unsteady.average_from_step > unsteady.steps {
        return refuse(
            "unsteady.average_from_step",
            "must lie in 1..=steps (1-based)",
        );
    }
    if unsteady.probe.is_some_and(|c| c >= domain.cell_count()) {
        return refuse("unsteady.probe", "cell outside the domain");
    }
    if config.fan.is_some() {
        return refuse(
            "simple.fan",
            "a face fan's operating-point iteration is steady; use an internal fan",
        );
    }
    for (side, program) in unsteady.opening_pressure_programs.iter().enumerate() {
        if program.is_some()
            && !matches!(
                config.faces[side],
                FvBoundary::Pressure { .. } | FvBoundary::Outlet
            )
        {
            return refuse(
                "unsteady.opening_pressure_programs",
                "each selected plane requires an original pressure opening",
            );
        }
    }
    Ok(())
}

/// Iterate one step's sweeps to its tolerance.
fn converge_step(
    solver: &mut Solver<'_>,
    unsteady: &UnsteadyConfig,
    step: usize,
    gate: &CancelGate,
) -> Result<(usize, (f64, f64)), ChtError> {
    let mut residuals = (f64::INFINITY, f64::INFINITY);
    let mut tolerance = unsteady.inner_tolerance;
    for inner in 1..=unsteady.inner_iterations {
        residuals = solver.sweep(gate)?;
        if !(residuals.0.is_finite() && residuals.1.is_finite()) {
            return Err(ChtError::FlowDiverged { step });
        }
        if residuals.0 <= unsteady.inner_tolerance && residuals.1 <= unsteady.inner_tolerance {
            if solver.has_coriolis() {
                let (projected, required) =
                    solver.project_with_momentum_tolerance(unsteady.inner_tolerance, gate)?;
                residuals.1 = projected;
                tolerance = required;
                if projected > required {
                    continue;
                }
            }
            return Ok((inner, residuals));
        }
    }
    Err(ChtError::SolverNotConverged {
        system: "unsteady flow step",
        iterations: unsteady.inner_iterations,
        relative_residual: residuals.0.max(residuals.1),
        tolerance,
    })
}

/// Shared march state: per-step flow solve, record, and running mean.
struct Marcher<'a> {
    solver: Solver<'a>,
    domain: &'a VoxelDomain,
    unsteady: &'a UnsteadyConfig,
    records: Vec<FlowStepRecord>,
    mean: Vec<[f64; 3]>,
    mean_flux: Option<[Vec<f64>; 3]>,
    averaged: usize,
    sweeps: usize,
    residuals: (f64, f64),
}

impl<'a> Marcher<'a> {
    fn new(
        domain: &'a VoxelDomain,
        fluid: &FluidProperties,
        config: &'a SimpleConfig,
        unsteady: &'a UnsteadyConfig,
    ) -> Self {
        let mut solver = Solver::new(domain, fluid, config);
        solver.start_time(unsteady.time_step_s, unsteady.scheme);
        Self {
            solver,
            domain,
            unsteady,
            records: Vec::with_capacity(unsteady.steps),
            mean: vec![[0.0; 3]; domain.cell_count()],
            mean_flux: None,
            averaged: 0,
            sweeps: 0,
            residuals: (f64::INFINITY, f64::INFINITY),
        }
    }

    /// Solve step `step` (1-based) with the inlets scaled by `inflow`,
    /// project, record, and accumulate.
    fn step(&mut self, step: usize, inflow: f64, gate: &CancelGate) -> Result<(), ChtError> {
        let inner = self.solve_step(step, inflow, gate)?;
        self.record_step(step, inner, gate)
    }

    /// Solve the current native time-discrete equations without advancing
    /// their accepted history or recording a physical sample. Coupled trials
    /// may reuse the current iterate while previous/older time levels stay fixed.
    fn solve_step(
        &mut self,
        step: usize,
        inflow: f64,
        gate: &CancelGate,
    ) -> Result<usize, ChtError> {
        finite("unsteady.inlet_schedule", inflow)?;
        self.solver.scale_inlets(inflow);
        let (inner, residuals) = converge_step(&mut self.solver, self.unsteady, step, gate)?;
        if !self.solver.has_coriolis() {
            self.solver.project(gate)?;
            self.solver.update_sgs(gate)?;
        }
        self.sweeps += inner;
        self.residuals = residuals;
        Ok(inner)
    }

    fn record_step(
        &mut self,
        step: usize,
        inner: usize,
        gate: &CancelGate,
    ) -> Result<(), ChtError> {
        let kinetic = self.solver.kinetic_energy_j(gate)?;
        let averaging = step >= self.unsteady.average_from_step;
        if averaging {
            for c in 0..self.domain.cell_count() {
                if c.is_multiple_of(256) {
                    poll(gate)?;
                }
                if !self.domain.is_fluid(c) {
                    continue;
                }
                let v = self.solver.cell_velocity(c);
                for a in 0..3 {
                    self.mean[c][a] += v[a];
                }
            }
            self.averaged += 1;
            let field = self.solver.field();
            let faces = field.face_arrays();
            match &mut self.mean_flux {
                Some(sum) => {
                    for (total, now) in sum.iter_mut().zip(faces) {
                        for (t, q) in total.iter_mut().zip(now) {
                            *t += q;
                        }
                    }
                }
                None => self.mean_flux = Some(faces.map(<[f64]>::to_vec)),
            }
        }
        self.records.push(FlowStepRecord {
            time_s: step as f64 * self.unsteady.time_step_s,
            inner_iterations: inner,
            mass_residual: self.residuals.0,
            momentum_residual: self.residuals.1,
            probe_velocity_m_s: self
                .unsteady
                .probe
                .map_or([0.0; 3], |c| self.solver.cell_velocity(c)),
            kinetic_energy_j: kinetic,
        });
        Ok(())
    }

    fn finish(self, fluid: &FluidProperties, gate: &CancelGate) -> Result<UnsteadyFlow, ChtError> {
        let scale = 1.0 / self.averaged.max(1) as f64;
        let mean_velocity_m_s = self.mean.iter().map(|v| v.map(|x| x * scale)).collect();
        let averaged_steps = self.averaged;
        let records = self.records;
        let mean_field = match self.mean_flux {
            Some([fx, fy, fz]) => FlowField::from_face_arrays(
                self.domain,
                fx.into_iter().map(|q| q * scale).collect(),
                fy.into_iter().map(|q| q * scale).collect(),
                fz.into_iter().map(|q| q * scale).collect(),
            ),
            None => self.solver.field(),
        };
        let flow = if self.solver.has_coriolis() {
            // The accepted endpoint already passed tight continuity and its
            // frozen momentum gate. History has since rolled forward for the
            // next interval; reassembling now would audit an unadvanced step.
            poll(gate)?;
            self.solver.guard_finite()?;
            self.solver
                .finish_projected(fluid, self.sweeps, self.residuals)
        } else {
            self.solver
                .finish(fluid, self.sweeps, self.residuals, gate)?
        };
        Ok(UnsteadyFlow {
            flow,
            records,
            mean_velocity_m_s,
            mean_field,
            averaged_steps,
        })
    }
}

/// March unsteady incompressible flow from rest; the declared inlet
/// velocities are scaled by `inlet_schedule(t)` at each step's end time
/// (`|_| 1.0`: an impulsive start), wall velocities apply throughout.
///
/// # Errors
/// Input refusals (as [`super::simple_flow`], a face fan, an empty march),
/// [`ChtError::SolverNotConverged`] (system `"unsteady flow step"`) when a
/// step misses its tolerance, [`ChtError::FlowDiverged`], or
/// [`ChtError::Cancelled`].
pub fn simple_unsteady(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    config: &SimpleConfig,
    unsteady: &UnsteadyConfig,
    inlet_schedule: impl Fn(f64) -> f64,
    gate: &CancelGate,
) -> Result<UnsteadyFlow, ChtError> {
    simple_unsteady_inner(
        domain,
        fluid,
        config,
        unsteady,
        None,
        inlet_schedule,
        |_| 1.0,
        gate,
    )
}

/// March from rest with a time-dependent differential reservoir pressure.
///
/// At each step's end time, the finite signed `pressure_schedule(t)` multiplies
/// [`UnsteadyConfig::pressure_program`]. Their product modulates every authored
/// Pressure/Outlet value, including spatial profiles, before adding that plane's
/// optional [`UnsteadyConfig::opening_pressure_programs`] offset in Pa:
/// `p_face(t) = p_ref0 + factor(t)*(p_authored - p_ref0) + offset_side(t)`.
/// Each disconnected component keeps its original anchor opening; numerical
/// pressure is centered on `p_ref0 + offset_anchor(t)`. Neither the factor
/// nor independently scheduled offsets multiply a common atmospheric pressure.
///
/// A factor of one exactly reproduces [`simple_unsteady`]. Without opening
/// offset programs, zero equalizes each component's reservoir pressures and a
/// negative factor reverses their head. Offsets remain independent of the factor.
/// Existing inertia, viscosity, body loads and internal
/// fans still determine the actual flow, including backflow. The same native
/// solver, BE/BDF2 velocity history and pressure projection persist throughout
/// the march. The callback is sampled once per endpoint, never inside a
/// nonlinear sweep; resolve authored jumps with the time grid.
///
/// The callback remains scalar modulation of the original differential profile.
/// Independent per-plane mechanical pressure histories use the offset programs;
/// neither control evolves thermodynamic pressure. Reservoir tangential
/// velocities and other loads retain their own controls.
/// Energy and coupled thermal marches are not advanced by this flow-only entry.
///
/// # Errors
/// The refusals of [`simple_unsteady`], a non-finite schedule sample, or an
/// unrepresentable scaled pressure difference or reservoir pressure. No partial
/// march is returned on refusal or cancellation; caller-owned inputs are unchanged.
pub fn simple_unsteady_with_pressure_schedule(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    config: &SimpleConfig,
    unsteady: &UnsteadyConfig,
    inlet_schedule: impl Fn(f64) -> f64,
    pressure_schedule: impl Fn(f64) -> f64,
    gate: &CancelGate,
) -> Result<UnsteadyFlow, ChtError> {
    simple_unsteady_inner(
        domain,
        fluid,
        config,
        unsteady,
        None,
        inlet_schedule,
        pressure_schedule,
        gate,
    )
}

/// March transient flow with an explicitly selected Smagorinsky/WALE closure.
///
/// The filter width is `filter_width_factor * volume^(1/3)`. The closure is
/// refreshed from the latest resolved velocity at each nonlinear sweep and
/// after each step's pressure projection. It never enables itself through
/// [`SimpleConfig`], and a zero coefficient exactly recovers [`simple_unsteady`].
/// This is an exploration entry; turbulent DNS and wall resolution require
/// independent spatial and temporal convergence evidence.
///
/// # Errors
/// The refusals of [`simple_unsteady`], invalid SGS controls, acceptance use,
/// an eddy/molecular ratio above the declared cap, or simultaneous LVEL.
pub fn simple_unsteady_with_sgs(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    config: &SimpleConfig,
    unsteady: &UnsteadyConfig,
    sgs: FvSgsConfig,
    inlet_schedule: impl Fn(f64) -> f64,
    gate: &CancelGate,
) -> Result<UnsteadyFlow, ChtError> {
    simple_unsteady_inner(
        domain,
        fluid,
        config,
        unsteady,
        Some(sgs),
        inlet_schedule,
        |_| 1.0,
        gate,
    )
}

/// Combine the original transient SGS closure with differential-pressure forcing.
///
/// Pressure, sampling and rollback semantics are those of
/// [`simple_unsteady_with_pressure_schedule`]. The original
/// [`simple_unsteady_with_sgs`] closure, resolved-gradient update and regime
/// restrictions continue to apply; pressure modulation adds no new SGS model.
///
/// # Errors
/// The refusals of both original entries, including invalid pressure samples.
#[allow(clippy::too_many_arguments)]
pub fn simple_unsteady_with_sgs_and_pressure_schedule(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    config: &SimpleConfig,
    unsteady: &UnsteadyConfig,
    sgs: FvSgsConfig,
    inlet_schedule: impl Fn(f64) -> f64,
    pressure_schedule: impl Fn(f64) -> f64,
    gate: &CancelGate,
) -> Result<UnsteadyFlow, ChtError> {
    simple_unsteady_inner(
        domain,
        fluid,
        config,
        unsteady,
        Some(sgs),
        inlet_schedule,
        pressure_schedule,
        gate,
    )
}

#[allow(clippy::too_many_arguments)]
fn simple_unsteady_inner(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    config: &SimpleConfig,
    unsteady: &UnsteadyConfig,
    sgs: Option<FvSgsConfig>,
    inlet_schedule: impl Fn(f64) -> f64,
    pressure_schedule: impl Fn(f64) -> f64,
    gate: &CancelGate,
) -> Result<UnsteadyFlow, ChtError> {
    admit(domain, fluid, config)?;
    admit_unsteady(domain, config, unsteady)?;
    let mut marcher = Marcher::new(domain, fluid, config, unsteady);
    if let Some(sgs) = sgs {
        marcher.solver.enable_sgs(sgs, gate)?;
    }
    for step in 1..=unsteady.steps {
        poll(gate)?;
        let time_s = step as f64 * unsteady.time_step_s;
        let inflow = inlet_schedule(time_s);
        marcher.solver.set_pressure_loads(
            unsteady.pressure_program.sample(time_s)? * pressure_schedule(time_s),
            unsteady.opening_pressure_offsets_pa(time_s)?,
            gate,
        )?;
        marcher.step(step, inflow, gate)?;
        marcher.solver.advance_time();
    }
    marcher.finish(fluid, gate)
}

/// March flow and energy together from rest and `initial_temperature`;
/// inlets follow `inlet_schedule(t)` and sources `power_schedule(t)`;
/// buoyancy (when declared) follows the previous step's temperature, and
/// radiation (when declared) converges emission and irradiation at the new endpoint.
/// [`UnsteadyConfig::pressure_program`] scales the original differential
/// reservoir profile and [`UnsteadyConfig::opening_pressure_programs`] add
/// independent plane offsets once per endpoint, before momentum/thermal feedback.
/// The projected signed flux drives the original thermal backflow boundaries.
/// Physical-heating applicability is unchanged: this does not admit open
/// reservoirs into a transient heat owner that requires a closed box.
///
/// # Errors
/// The refusals of [`simple_unsteady`] and [`super::march_energy`]
/// (including ordinary solids or compact junctions without heat capacity), or
/// [`ChtError::Cancelled`].
#[allow(clippy::too_many_arguments)]// physics inputs + both controls
pub fn march_conjugate(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    solids: &[SolidMaterial],
    setup: &ThermalSetup,
    config: &SimpleConfig,
    unsteady: &UnsteadyConfig,
    buoyancy: Option<&Boussinesq>,
    initial_temperature: &[f64],
    inlet_schedule: impl Fn(f64) -> f64,
    power_schedule: impl Fn(f64) -> f64,
    energy: &EnergyConfig,
    radiation: Option<&RadiationConfig>,
    gate: &CancelGate,
) -> Result<ConjugateMarch, ChtError> {
    march_conjugate_inner(
        domain,
        fluid,
        solids,
        setup,
        config,
        unsteady,
        None,
        None,
        buoyancy,
        initial_temperature,
        inlet_schedule,
        power_schedule,
        energy,
        radiation,
        gate,
    )
}

/// March flow and heat with explicit transient Smagorinsky/WALE transport.
///
/// The same SGS viscosity enters complete symmetric momentum stress and
/// fluid thermal conductivity `rho * cp * nu_t / Pr_t`, with the positive
/// `Pr_t` declared in [`SimpleConfig::turbulent_prandtl`] (default 0.9).
/// Conductivity is refreshed from the projected flow before every energy
/// solve. This gradient-diffusion model adds neither SGS dissipation heating
/// nor an equilibrium wall heat-transfer law. Caller-authored additional
/// eddy conductivity refuses so that two thermal closures cannot be combined
/// or silently overwritten. Otherwise coupling follows [`march_conjugate`].
///
/// # Errors
/// The refusals of [`march_conjugate`] and [`simple_unsteady_with_sgs`], or a
/// supplied nonzero eddy-conductivity field.
#[allow(clippy::too_many_arguments)]// same physics inputs + explicit closure
pub fn march_conjugate_with_sgs(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    solids: &[SolidMaterial],
    setup: &ThermalSetup,
    config: &SimpleConfig,
    unsteady: &UnsteadyConfig,
    sgs: FvSgsConfig,
    buoyancy: Option<&Boussinesq>,
    initial_temperature: &[f64],
    inlet_schedule: impl Fn(f64) -> f64,
    power_schedule: impl Fn(f64) -> f64,
    energy: &EnergyConfig,
    radiation: Option<&RadiationConfig>,
    gate: &CancelGate,
) -> Result<ConjugateMarch, ChtError> {
    march_conjugate_inner(
        domain,
        fluid,
        solids,
        setup,
        config,
        unsteady,
        Some(sgs),
        None,
        buoyancy,
        initial_temperature,
        inlet_schedule,
        power_schedule,
        energy,
        radiation,
        gate,
    )
}

/// Reused per-voxel coefficient buffers; material sources remain in their
/// adapter, while the native solver owns momentum, temperature and time.
#[derive(Clone)]
pub(crate) struct MaterialCoefficients {
    pub viscosity_pa_s: Vec<f64>,
    pub conductivity_w_m_k: Vec<[f64; 3]>,
    pub force_n_m3: Vec<[f64; 3]>,
}

/// One native thermal interval. A constitutive iteration may replace its
/// storage linearization and conductivity, never its clock, old state or fluxes.
pub(crate) struct MaterialHeatStep<'a> {
    pub domain: &'a VoxelDomain,
    pub setup: &'a ThermalSetup,
    pub previous_temperature_k: &'a [f64],
    pub reference_capacity_j_k: &'a [f64],
    reference_storage_w_k: &'a [f64],
    pub time_s: f64,
    pub time_step_s: f64,
    /// Generated physical heat immediately rejected by fixed-cell thermostats.
    /// It is both a real source and reservoir outflow, never stored twice.
    physical_fixed_heat_j: f64,
    fluid: &'a FluidProperties,
    solids: &'a [SolidMaterial],
    pub flow: &'a FlowField,
    energy: &'a EnergyConfig,
    solids_present: bool,
    radiation: Option<(&'a Radiators, &'a RadiationConfig)>,
    /// Optional outer endpoint budget, shared by every caloric/radiative trial.
    endpoint_energy_budget: Option<&'a Cell<usize>>,
    gate: &'a CancelGate,
}

impl MaterialHeatStep<'_> {
    /// The exact native face scheme used by this thermal interval. Material
    /// acceptance must re-evaluate that same source-enthalpy reconstruction.
    pub(crate) fn convection_scheme(&self) -> super::ConvectionScheme {
        self.energy.scheme
    }

    /// Only validated compact aliases have zero native cell capacity.
    pub(crate) fn is_compact_cell(&self, cell: usize) -> bool {
        self.reference_capacity_j_k[cell] == 0.0
    }

    /// Append the independent junction rows to a material owner's cell rates.
    /// Constant capacities are final; C(T) entries are only reference rates.
    /// The shared endpoint owner refreshes their exact enthalpy secants.
    pub(crate) fn append_junction_storage(
        &self,
        storage_w_k: &mut Vec<f64>,
    ) -> Result<(), ChtError> {
        let cells = self.domain.cell_count();
        if storage_w_k.len() != cells
            || self.reference_capacity_j_k.len() != cells + self.setup.compact_components.len()
        {
            return Err(ChtError::InvalidInput {
                field: "thermal.storage",
                reason: "material storage requires one rate per cell before junction rates".into(),
            });
        }
        for (j, part) in self.setup.compact_components.iter().enumerate() {
            if j.is_multiple_of(256) {
                poll(self.gate)?;
            }
            let rate = self.reference_capacity_j_k[cells + j] / self.time_step_s;
            let old =
                self.previous_temperature_k[self.domain.index(part.lo[0], part.lo[1], part.lo[2])];
            finite_positive("thermal.storage.junction_coefficient_w_k", rate)?;
            finite("thermal.storage.junction_rhs_w", rate * old)?;
            storage_w_k.push(rate);
        }
        Ok(())
    }

    pub(crate) fn junction_stored_energy_change_j(
        &self,
        temperature: &[f64],
    ) -> Result<f64, ChtError> {
        junction_stored_energy_change_j(
            self.domain,
            self.setup,
            self.reference_capacity_j_k,
            self.previous_temperature_k,
            temperature,
            self.gate,
        )
    }

    /// Exact ordinary native storage for source-fluid marches. Source solid
    /// owners replace their own cells separately; compact aliases remain zero.
    pub(crate) fn native_stored_energy_change_j(
        &self, cell: usize, temperature_k: f64,
    ) -> Result<f64, ChtError> {
        super::transient::native_cell_stored_energy_change_j(
            self.domain, self.solids, self.reference_capacity_j_k,
            self.previous_temperature_k, cell, temperature_k,
        )
    }

    pub(crate) fn solve(
        &self,
        capacity_j_k: &[f64],
        storage_w_k: &[f64],
        conductivity_w_m_k: &[[f64; 3]],
    ) -> Result<(Vec<f64>, TransientRecord), ChtError> {
        self.solve_with_optional_conductivity(
            capacity_j_k, storage_w_k, Some(conductivity_w_m_k),
        )
    }

    /// Native k(T), wall/contact and compact laws remain inside the heat owner.
    fn solve_native(&self) -> Result<(Vec<f64>, TransientRecord), ChtError> {
        self.solve_with_optional_conductivity(
            self.reference_capacity_j_k, self.reference_storage_w_k, None,
        )
    }

    fn solve_with_optional_conductivity(
        &self,
        capacity_j_k: &[f64],
        storage_w_k: &[f64],
        conductivity_w_m_k: Option<&[[f64; 3]]>,
    ) -> Result<(Vec<f64>, TransientRecord), ChtError> {
        let (temperature, record) = if let Some(budget) = self.endpoint_energy_budget {
            let mut remaining = budget.get();
            let solved = solve_step_with_fields(
                self.domain,
                self.fluid,
                self.solids,
                self.flow,
                self.setup,
                self.energy,
                storage_w_k,
                self.previous_temperature_k,
                self.time_step_s,
                self.radiation,
                conductivity_w_m_k,
                None,
                Some(&mut remaining),
                self.gate,
            );
            budget.set(remaining);
            let (next, radiated_w, correction_w) = solved?;
            let record = storage_record(
                self.domain,
                self.solids,
                self.setup,
                capacity_j_k,
                self.previous_temperature_k,
                &next,
                (self.time_s, self.time_step_s, self.solids_present),
                (radiated_w, correction_w),
                self.gate,
            )?;
            (next.temperature, record)
        } else {
            energy_step_with_conductivity(
                self.domain,
                self.fluid,
                self.solids,
                self.flow,
                self.setup,
                self.energy,
                (capacity_j_k, storage_w_k),
                self.previous_temperature_k,
                (self.time_s, self.time_step_s, self.solids_present),
                self.radiation,
                conductivity_w_m_k,
                self.gate,
            )?
        };
        Ok((temperature, self.record_physical_heat(record)?))
    }

    /// A fixed-temperature row bypasses ordinary source assembly. Preserve
    /// its actual generated molecular heat as equal source/reservoir entries
    /// before the source-caloric owner's nonlinear closure test.
    fn record_physical_heat(
        &self,
        mut record: TransientRecord,
    ) -> Result<TransientRecord, ChtError> {
        if self.physical_fixed_heat_j != 0.0 {
            record.source_j += self.physical_fixed_heat_j;
            record.boundary_outflow_j += self.physical_fixed_heat_j;
            // Equal physical source/reservoir terms cancel exactly; retain the
            // original closure so cancellation cannot hide a nonlinear defect.
            for value in [record.source_j, record.boundary_outflow_j, record.closure_j] {
                finite("thermal.viscous_heat.energy_j", value)?;
            }
        }
        Ok(record)
    }

    /// Source-caloric face transport on this same accepted native interval.
    /// The original endpoint radiation owner also consumes the solve budget.
    pub(crate) fn solve_with_enthalpy(
        &self,
        storage_w_k: &[f64],
        conductivity_w_m_k: &[[f64; 3]],
        enthalpy: AdvectiveEnthalpyField<'_>,
        thermal: &ThermalSetup,
        remaining_solves: &mut usize,
    ) -> Result<(EnergySolution, f64, f64), ChtError> {
        let mut allowed = self
            .endpoint_energy_budget
            .map_or(*remaining_solves, |budget| {
                (*remaining_solves).min(budget.get())
            });
        let opening = allowed;
        let result = solve_step_with_fields(
            self.domain,
            self.fluid,
            self.solids,
            self.flow,
            thermal,
            self.energy,
            storage_w_k,
            self.previous_temperature_k,
            self.time_step_s,
            self.radiation,
            Some(conductivity_w_m_k),
            Some(enthalpy),
            Some(&mut allowed),
            self.gate,
        );
        let spent = opening - allowed;
        *remaining_solves -= spent;
        if let Some(budget) = self.endpoint_energy_budget {
            budget.set(budget.get() - spent);
        }
        result
    }

    /// Endpoint source-Cp SGS conductivity on the SAME native face/link laws.
    /// This observes a candidate only; it neither solves nor advances time.
    pub(crate) fn conductivity_flux_defect(
        &self,
        solution: &EnergySolution,
        molecular_w_m_k: &[[f64; 3]],
        thermal: &ThermalSetup,
        net_defects_w: &mut Vec<f64>,
    ) -> Result<(f64, f64), ChtError> {
        solution.conductivity_flux_defect(
            self.domain,
            self.solids,
            self.flow,
            thermal,
            molecular_w_m_k,
            net_defects_w,
            self.gate,
        )
    }

    /// Actual constitutive storage with the original boundary/reservoir ledger.
    pub(crate) fn record_with_enthalpy(
        &self,
        solution: &EnergySolution,
        stored_by_cell_j: &[f64],
        radiated_w: f64,
        outflow_correction_w: f64,
    ) -> Result<TransientRecord, ChtError> {
        if stored_by_cell_j.len() != self.domain.cell_count() {
            return Err(ChtError::InvalidInput {
                field: "thermal.storage",
                reason: "one actual source storage increment per cell is required".into(),
            });
        }
        for (cell, &stored) in stored_by_cell_j.iter().enumerate() {
            if cell.is_multiple_of(4096) {
                poll(self.gate)?;
            }
            if self.is_compact_cell(cell) && stored != 0.0 {
                return Err(ChtError::InvalidInput {
                    field: "thermal.storage",
                    reason: "compact aliases cannot also carry source-material storage".into(),
                });
            }
        }
        let mut record = record_energy_step(
            self.domain,
            self.setup,
            solution,
            |cell| stored_by_cell_j[cell],
            (self.time_s, self.time_step_s, self.solids_present),
            radiated_w,
            outflow_correction_w,
        )?;
        record.stored_energy_change_j +=
            self.junction_stored_energy_change_j(&solution.temperature)?;
        record.closure_j =
            record.stored_energy_change_j + record.boundary_outflow_j - record.source_j;
        finite("thermal.storage.energy_j", record.stored_energy_change_j)?;
        finite("thermal.storage.energy_j", record.closure_j)?;
        self.record_physical_heat(record)
    }
}

struct MaterialCoupling<'a, E> {
    fields: &'a mut MaterialCoefficients,
    refresh: &'a mut dyn FnMut(&[f64], &mut MaterialCoefficients) -> Result<(), E>,
}

/// Internal source-material seam on the same march. Typed source errors return
/// unchanged; neither a second time owner nor a replacement energy ledger exists.
#[allow(clippy::too_many_arguments)]
pub(crate) fn march_conjugate_with_material<E: From<ChtError>>(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    solids: &[SolidMaterial],
    setup: &ThermalSetup,
    config: &SimpleConfig,
    unsteady: &UnsteadyConfig,
    sgs: Option<FvSgsConfig>,
    #[cfg(feature = "enthalpy-coupling")] viscous_heat: Option<&mut MaterialViscousHeat>,
    initial_temperature: &[f64],
    inlet_schedule: impl Fn(f64) -> f64,
    power_schedule: impl Fn(f64) -> f64,
    energy: &EnergyConfig,
    radiation: Option<&RadiationConfig>,
    fields: &mut MaterialCoefficients,
    refresh: &mut dyn FnMut(&[f64], &mut MaterialCoefficients) -> Result<(), E>,
    material_heat: Option<
        &mut dyn FnMut(
            &MaterialHeatStep<'_>,
            &mut MaterialCoefficients,
        ) -> Result<(Vec<f64>, TransientRecord), E>,
    >,
    endpoint: Option<&mut MaterialEndpointCoupling>,
    gate: &CancelGate,
) -> Result<ConjugateMarch, E> {
    march_conjugate_inner_with_material(
        domain,
        fluid,
        solids,
        setup,
        config,
        unsteady,
        sgs,
        #[cfg(feature = "enthalpy-coupling")]
        viscous_heat,
        None,
        None,
        Some(MaterialCoupling { fields, refresh }),
        material_heat,
        initial_temperature,
        inlet_schedule,
        power_schedule,
        energy,
        radiation,
        endpoint,
        gate,
    )
}

#[allow(clippy::too_many_arguments)]
fn march_conjugate_inner(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    solids: &[SolidMaterial],
    setup: &ThermalSetup,
    config: &SimpleConfig,
    unsteady: &UnsteadyConfig,
    sgs: Option<FvSgsConfig>,
    wall_heat: Option<&mut wall_heat::Coupling<'_>>,
    buoyancy: Option<&Boussinesq>,
    initial_temperature: &[f64],
    inlet_schedule: impl Fn(f64) -> f64,
    power_schedule: impl Fn(f64) -> f64,
    energy: &EnergyConfig,
    radiation: Option<&RadiationConfig>,
    gate: &CancelGate,
) -> Result<ConjugateMarch, ChtError> {
    march_conjugate_inner_with_material(
        domain,
        fluid,
        solids,
        setup,
        config,
        unsteady,
        sgs,
        #[cfg(feature = "enthalpy-coupling")]
        None,
        wall_heat,
        buoyancy,
        None,
        None,
        initial_temperature,
        inlet_schedule,
        power_schedule,
        energy,
        radiation,
        None,
        gate,
    )
}

#[allow(clippy::too_many_arguments)] // one owner for all transient entry paths
fn march_conjugate_inner_with_material<E: From<ChtError>>(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    solids: &[SolidMaterial],
    setup: &ThermalSetup,
    config: &SimpleConfig,
    unsteady: &UnsteadyConfig,
    sgs: Option<FvSgsConfig>,
    #[cfg(feature = "enthalpy-coupling")] mut viscous_heat: Option<&mut MaterialViscousHeat>,
    mut wall_heat: Option<&mut wall_heat::Coupling<'_>>,
    buoyancy: Option<&Boussinesq>,
    mut material: Option<MaterialCoupling<'_, E>>,
    mut material_heat: Option<
        &mut dyn FnMut(
            &MaterialHeatStep<'_>,
            &mut MaterialCoefficients,
        ) -> Result<(Vec<f64>, TransientRecord), E>,
    >,
    initial_temperature: &[f64],
    inlet_schedule: impl Fn(f64) -> f64,
    power_schedule: impl Fn(f64) -> f64,
    energy: &EnergyConfig,
    radiation: Option<&RadiationConfig>,
    mut endpoint: Option<&mut MaterialEndpointCoupling>,
    gate: &CancelGate,
) -> Result<ConjugateMarch, E> {
    admit(domain, fluid, config)?;
    admit_unsteady(domain, config, unsteady)?;
    if let Some(coupling) = endpoint.as_deref() {
        coupling.admit(unsteady)?;
        if material.is_none() || wall_heat.is_some() {
            return Err(ChtError::InvalidInput {
                field: "unsteady.material_endpoint",
                reason: "implicit material feedback requires its original molecular material owner without a wall closure".into(),
            }.into());
        }
    }
    #[cfg(feature = "enthalpy-coupling")]
    if viscous_heat.is_some()
        && (unsteady.scheme != TimeScheme::BackwardEuler
            || sgs.is_some()
            || material.is_none()
            || material_heat.is_none()
            || wall_heat.is_some())
    {
        return Err(ChtError::InvalidInput {
            field: "unsteady.material_viscous_heat",
            reason: "source-caloric physical heating requires native backward-Euler Stokes without SGS or wall closures".into(),
        }.into());
    }
    if material.is_some()
        && (buoyancy.is_some()
            || wall_heat.is_some()
            || config.turbulence != Turbulence::Laminar
            || !setup.eddy_conductivity_w_m_k.is_empty())
    {
        return Err(ChtError::InvalidInput {
            field: "unsteady.material",
            reason: "material transport owns buoyancy and molecular conductivity; scalar beta, equilibrium LVEL and wall closures cannot be combined".into(),
        }.into());
    }
    if let Some(sgs) = sgs {
        sgs.validate()?;
        if setup.eddy_conductivity_w_m_k.iter().any(|&k| k != 0.0) {
            return Err(ChtError::InvalidInput {
                field: "thermal.eddy_conductivity_w_m_k",
                reason: "explicit SGS supplies the fluid eddy conductivity; a second field cannot be combined".into(),
            }.into());
        }
    }
    let cells = domain.cell_count();
    if initial_temperature.len() != cells {
        return Err(ChtError::InvalidInput {
            field: "unsteady.initial_temperature",
            reason: format!(
                "expected {cells} entries, got {}",
                initial_temperature.len()
            ),
        }
        .into());
    }
    for &t in initial_temperature {
        finite("unsteady.initial_temperature", t)?;
    }
    if let Some(b) = buoyancy {
        for g in b.gravity_m_s2 {
            finite("buoyancy.gravity_m_s2", g)?;
        }
        finite("buoyancy.expansion_per_k", b.expansion_per_k)?;
        finite(
            "buoyancy.reference_temperature_k",
            b.reference_temperature_k,
        )?;
    }
    if material.is_none() || radiation.is_some() {
        // The source-field heat path can omit unused compact material slots.
        // Radiation still requires its separately declared surface materials.
        domain.check_materials(solids.len())?;
    }
    let capacity = heat_capacity(domain, fluid, solids, setup, initial_temperature, gate)?;
    let storage: Vec<f64> = capacity.iter().map(|c| c / unsteady.time_step_s).collect();
    let solids_present = (0..cells).any(|c| !domain.is_fluid(c));
    let radiators = match radiation {
        Some(r) => Some(Radiators::build(domain, solids, r, gate)?),
        None => None,
    };
    let mut marcher = Marcher::new(domain, fluid, config, unsteady);
    #[cfg(feature = "enthalpy-coupling")]
    if let Some(heating) = &viscous_heat {
        marcher
            .solver
            .enable_unsteady_stokes_heat(heating.controls)?;
    }
    if let Some(sgs) = sgs {
        marcher.solver.enable_sgs(sgs, gate)?;
    }
    let mut temperature = initial_temperature.to_vec();
    let mut mean_temperature = vec![0.0; cells];
    let mut energy_records = Vec::with_capacity(unsteady.steps);
    let mut stepped = setup.clone();
    if let Some(coupling) = material.as_mut() {
        if let Some(feedback) = endpoint.as_deref_mut() {
            feedback.refresh(coupling, &temperature)?;
        } else {
            (coupling.refresh)(&temperature, coupling.fields)?;
        }
    }
    for step in 1..=unsteady.steps {
        poll(gate)?;
        let time_s = step as f64 * unsteady.time_step_s;
        // Hold the same authored reservoir load through every coupled trial.
        // The original solver, pressure datum and time history remain intact.
        marcher.solver.set_pressure_loads(
            unsteady.pressure_program.sample(time_s)?,
            unsteady.opening_pressure_offsets_pa(time_s)?,
            gate,
        )?;
        let (next, record) = if let Some(feedback) = endpoint.as_deref_mut() {
            let inflow = inlet_schedule(time_s);
            let scale = power_schedule(time_s);
            finite("unsteady.power_schedule", scale)?;
            scale_power(&mut stepped, setup, scale);
            #[cfg(feature = "enthalpy-coupling")]
            if viscous_heat.is_some() && setup.power_w.is_empty() {
                // An empty authored field must not retain last step's deposit.
                stepped.power_w.fill(0.0);
            }
            feedback.step(
                &mut marcher,
                material.as_mut().expect("admitted original material owner"),
                &mut material_heat,
                #[cfg(feature = "enthalpy-coupling")]
                viscous_heat.as_deref_mut(),
                fluid,
                solids,
                &mut stepped,
                energy,
                &capacity,
                &storage,
                &temperature,
                step,
                inflow,
                solids_present,
                radiators.as_ref().zip(radiation),
                gate,
            )?
        } else {
            if let Some(coupling) = material.as_ref() {
                marcher.solver.set_material_transport(
                    &coupling.fields.viscosity_pa_s,
                    Some(&coupling.fields.force_n_m3),
                    gate,
                )?;
            }
            if let Some(b) = buoyancy {
                marcher.solver.set_force(
                    (0..cells)
                        .map(|c| {
                            if !domain.is_fluid(c) {
                                return [0.0; 3];
                            }
                            let scale = -fluid.density_kg_m3
                                * b.expansion_per_k
                                * (temperature[c] - b.reference_temperature_k);
                            b.gravity_m_s2.map(|g| scale * g)
                        })
                        .collect(),
                );
            }
            marcher.step(step, inlet_schedule(time_s), gate)?;
            let scale = power_schedule(time_s);
            finite("unsteady.power_schedule", scale)?;
            scale_power(&mut stepped, setup, scale);
            #[cfg(feature = "enthalpy-coupling")]
            let mechanical_step = viscous_heat
                .as_ref()
                .map(|heating| {
                    marcher.solver.unsteady_mechanical_step(
                        (step - 1) as f64 * unsteady.time_step_s,
                        time_s,
                        material
                            .as_ref()
                            .map(|coupling| coupling.fields.force_n_m3.as_slice()),
                        heating.controls,
                        gate,
                    )
                })
                .transpose()?;
            let physical_fixed_heat_j = 0.0;
            #[cfg(feature = "enthalpy-coupling")]
            let physical_fixed_heat_j = if let Some(receipt) = &mechanical_step {
                let authored = if setup.power_w.is_empty() {
                    Vec::new()
                } else {
                    stepped.power_w.clone()
                };
                MaterialViscousHeat::deposit(&mut stepped, &authored, receipt, gate)?
            } else {
                physical_fixed_heat_j
            };

            if wall_heat.is_none()
                && (config.turbulence != Turbulence::Laminar || marcher.solver.has_sgs())
            {
                if marcher.solver.has_sgs() {
                    stepped.eddy_wall_treatment = super::EddyWallTreatment::Molecular;
                }
                stepped.eddy_conductivity_w_m_k = marcher
                    .solver
                    .heat_eddy_viscosity()
                    .iter()
                    .map(|nu_t| fluid.volumetric_heat_capacity() * nu_t / config.turbulent_prandtl)
                    .collect();
            }
            // One projected flux image is shared with wall preparation and energy.
            // `stepped` never receives generated wall corrections: it remains the
            // authored setup plus this step's source schedule.
            let field = marcher.solver.field();
            let prepared_walls = wall_heat
                .as_deref_mut()
                .map(|wall| wall.prepare(&marcher.solver, &field, &stepped, gate))
                .transpose()?;
            let thermal = prepared_walls
                .as_ref()
                .map_or(&stepped, |walls| &walls.thermal);
            let (next, record) = if let (Some(advance), Some(coupling)) =
                (material_heat.as_deref_mut(), material.as_mut())
            {
                let interval = MaterialHeatStep {
                    domain,
                    setup: thermal,
                    previous_temperature_k: &temperature,
                    reference_capacity_j_k: &capacity,
                    reference_storage_w_k: &storage,
                    time_s,
                    time_step_s: unsteady.time_step_s,
                    physical_fixed_heat_j,
                    fluid,
                    solids,
                    flow: &field,
                    energy,
                    solids_present,
                    radiation: radiators.as_ref().zip(radiation),
                    endpoint_energy_budget: None,
                    gate,
                };
                advance(&interval, coupling.fields)?
            } else {
                energy_step_with_conductivity(
                    domain,
                    fluid,
                    solids,
                    &field,
                    thermal,
                    energy,
                    (&capacity, &storage),
                    &temperature,
                    (time_s, unsteady.time_step_s, solids_present),
                    radiators.as_ref().zip(radiation),
                    material
                        .as_ref()
                        .map(|coupling| coupling.fields.conductivity_w_m_k.as_slice()),
                    gate,
                )?
            };
            // Validate every source at the accepted endpoint, including the final
            // step, before exposing results. Reuse these fields for the next step.
            if let Some(coupling) = material.as_mut() {
                (coupling.refresh)(&next, coupling.fields)?;
            }
            #[cfg(feature = "enthalpy-coupling")]
            if let (Some(heating), Some(receipt)) = (viscous_heat.as_deref_mut(), mechanical_step) {
                heating.audit.accept(&receipt, heating.controls)?;
                heating.records.push(receipt);
            }
            if let (Some(wall), Some(prepared)) = (wall_heat.as_deref_mut(), prepared_walls) {
                wall.record(prepared, &next, time_s, unsteady.time_step_s, gate)?;
            }
            (next, record)
        };
        energy_records.push(record);
        temperature = next;
        if step >= unsteady.average_from_step {
            for (m, t) in mean_temperature.iter_mut().zip(&temperature) {
                *m += t;
            }
        }
        #[cfg(feature = "enthalpy-coupling")]
        if viscous_heat.is_some() {
            // The physical work owner uses one BE interval. Start its next
            // interval from the exact accepted faces without a stale BDF level.
            marcher
                .solver
                .start_time(unsteady.time_step_s, TimeScheme::BackwardEuler);
        } else {
            marcher.solver.advance_time();
        }
        #[cfg(not(feature = "enthalpy-coupling"))]
        marcher.solver.advance_time();
    }
    let window = (unsteady.steps + 1 - unsteady.average_from_step) as f64;
    for m in &mut mean_temperature {
        *m /= window;
    }
    Ok(ConjugateMarch {
        flow: marcher.finish(fluid, gate)?,
        temperature,
        mean_temperature,
        energy_records,
    })
}
