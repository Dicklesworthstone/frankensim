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
//! The SIMPLEC no-claims apply (staircase walls, power-law convection, LVEL
//! as the only turbulence closure: an LVEL "transient" is an unsteady
//! eddy-viscosity solution, not a turbulence-resolving simulation). The
//! energy coupling is first order in time and lags buoyancy by one step.
//! Face fans and compact components are steady models and refuse here; use
//! internal fans.

use fs_exec::CancelGate;

use super::domain::{FluidProperties, SolidMaterial, VoxelDomain};
use super::energy::{EnergyConfig, ThermalSetup};
use super::simple::{FvFlow, SimpleConfig, Solver, TimeScheme, Turbulence, admit};
use super::transient::{TransientRecord, energy_step, heat_capacity};
use super::turbulence::TURBULENT_PRANDTL;
use super::{ChtError, finite, finite_positive, poll};

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
        }
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
    /// Kinetic energy `sum rho |u|^2 V / 2` of the fluid, J.
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

/// Boussinesq buoyancy for [`march_conjugate`].
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
    for inner in 1..=unsteady.inner_iterations {
        residuals = solver.sweep(gate)?;
        if !(residuals.0.is_finite() && residuals.1.is_finite()) {
            return Err(ChtError::FlowDiverged { step });
        }
        if residuals.0 <= unsteady.inner_tolerance && residuals.1 <= unsteady.inner_tolerance {
            return Ok((inner, residuals));
        }
    }
    Err(ChtError::SolverNotConverged {
        system: "unsteady flow step",
        iterations: unsteady.inner_iterations,
        relative_residual: residuals.0.max(residuals.1),
        tolerance: unsteady.inner_tolerance,
    })
}

/// Shared march state: per-step flow solve, record, and running mean.
struct Marcher<'a> {
    solver: Solver<'a>,
    domain: &'a VoxelDomain,
    rho: f64,
    unsteady: &'a UnsteadyConfig,
    records: Vec<FlowStepRecord>,
    mean: Vec<[f64; 3]>,
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
            rho: fluid.density_kg_m3,
            unsteady,
            records: Vec::with_capacity(unsteady.steps),
            mean: vec![[0.0; 3]; domain.cell_count()],
            averaged: 0,
            sweeps: 0,
            residuals: (f64::INFINITY, f64::INFINITY),
        }
    }

    /// Solve step `step` (1-based) with the inlets scaled by `inflow`,
    /// project, record, and accumulate.
    fn step(&mut self, step: usize, inflow: f64, gate: &CancelGate) -> Result<(), ChtError> {
        finite("unsteady.inlet_schedule", inflow)?;
        self.solver.scale_inlets(inflow);
        let (inner, residuals) = converge_step(&mut self.solver, self.unsteady, step, gate)?;
        self.solver.project(gate)?;
        self.sweeps += inner;
        self.residuals = residuals;
        let volume = self.domain.dx().powi(3);
        let mut kinetic = 0.0;
        let averaging = step >= self.unsteady.average_from_step;
        for c in 0..self.domain.cell_count() {
            if !self.domain.is_fluid(c) {
                continue;
            }
            let v = self.solver.cell_velocity(c);
            kinetic += 0.5 * self.rho * volume * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
            if averaging {
                for a in 0..3 {
                    self.mean[c][a] += v[a];
                }
            }
        }
        if averaging {
            self.averaged += 1;
        }
        self.records.push(FlowStepRecord {
            time_s: step as f64 * self.unsteady.time_step_s,
            inner_iterations: inner,
            mass_residual: residuals.0,
            momentum_residual: residuals.1,
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
        let flow = self
            .solver
            .finish(fluid, self.sweeps, self.residuals, gate)?;
        Ok(UnsteadyFlow {
            flow,
            records,
            mean_velocity_m_s,
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
    admit(domain, fluid, config)?;
    admit_unsteady(domain, config, unsteady)?;
    let mut marcher = Marcher::new(domain, fluid, config, unsteady);
    for step in 1..=unsteady.steps {
        poll(gate)?;
        marcher.step(
            step,
            inlet_schedule(step as f64 * unsteady.time_step_s),
            gate,
        )?;
        marcher.solver.advance_time();
    }
    marcher.finish(fluid, gate)
}

/// March flow and energy together from rest and `initial_temperature`;
/// inlets follow `inlet_schedule(t)` and sources `power_schedule(t)`;
/// buoyancy (when declared) follows the previous step's temperature.
///
/// # Errors
/// The refusals of [`simple_unsteady`] and [`super::march_energy`]
/// (including solids without heat capacity and compact components), or
/// [`ChtError::Cancelled`].
#[allow(clippy::too_many_arguments)] // physics inputs + both controls
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
    gate: &CancelGate,
) -> Result<ConjugateMarch, ChtError> {
    admit(domain, fluid, config)?;
    admit_unsteady(domain, config, unsteady)?;
    let cells = domain.cell_count();
    if initial_temperature.len() != cells {
        return Err(ChtError::InvalidInput {
            field: "unsteady.initial_temperature",
            reason: format!(
                "expected {cells} entries, got {}",
                initial_temperature.len()
            ),
        });
    }
    for &t in initial_temperature {
        finite("unsteady.initial_temperature", t)?;
    }
    if !setup.compact_components.is_empty() {
        return Err(ChtError::InvalidInput {
            field: "thermal.compact_components",
            reason: "two-resistor compact models are steady; the transient march refuses them"
                .into(),
        });
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
    domain.check_materials(solids.len())?;
    let capacity = heat_capacity(domain, fluid, solids)?;
    let storage: Vec<f64> = capacity.iter().map(|c| c / unsteady.time_step_s).collect();
    let solids_present = (0..cells).any(|c| !domain.is_fluid(c));
    let mut marcher = Marcher::new(domain, fluid, config, unsteady);
    let mut temperature = initial_temperature.to_vec();
    let mut mean_temperature = vec![0.0; cells];
    let mut energy_records = Vec::with_capacity(unsteady.steps);
    let mut stepped = setup.clone();
    for step in 1..=unsteady.steps {
        poll(gate)?;
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
        let time_s = step as f64 * unsteady.time_step_s;
        marcher.step(step, inlet_schedule(time_s), gate)?;
        let scale = power_schedule(time_s);
        finite("unsteady.power_schedule", scale)?;
        if !setup.power_w.is_empty() {
            for (out, base) in stepped.power_w.iter_mut().zip(&setup.power_w) {
                *out = base * scale;
            }
        }
        if config.turbulence != Turbulence::Laminar {
            stepped.eddy_conductivity_w_m_k = marcher
                .solver
                .heat_eddy_viscosity()
                .iter()
                .map(|nu_t| fluid.volumetric_heat_capacity() * nu_t / TURBULENT_PRANDTL)
                .collect();
        }
        let (next, record) = energy_step(
            domain,
            fluid,
            solids,
            &marcher.solver.field(),
            &stepped,
            energy,
            (&capacity, &storage),
            &temperature,
            (time_s, unsteady.time_step_s, solids_present),
            gate,
        )?;
        energy_records.push(record);
        temperature = next;
        if step >= unsteady.average_from_step {
            for (m, t) in mean_temperature.iter_mut().zip(&temperature) {
                *m += t;
            }
        }
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
