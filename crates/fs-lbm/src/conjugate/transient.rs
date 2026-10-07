//! Transient conjugate heat transfer over a frozen steady flow field.
//!
//! The energy equation `C dT/dt + A T = b(t)` (the steady conservative
//! operator of [`super::solve_energy`] plus storage) is advanced by backward
//! Euler, which is unconditionally stable for the M-matrix `A` and so admits
//! the ~10^4 conductivity contrast and the large solid heat capacities of
//! electronics without a stability limit on the step:
//!
//! ```text
//! (C / dt + A) T^(n+1) = (C / dt) T^n + b(t_(n+1))
//! ```
//!
//! `C_c = (rho c)_c V_c` uses the fluid's volumetric heat capacity in fluid
//! cells and each solid's declared `volumetric_heat_capacity_j_m3_k` (a
//! solid without one refuses). Sources are `power_w` scaled by
//! `power_schedule(t)`. Each step retains an energy closure computed two
//! independent ways: the stored-energy change `sum C (T^(n+1) - T^n)` and
//! `dt (source - boundary outflow)` from boundary fluxes alone; their
//! difference is the solver residual, never a modelling term.
//!
//! No-claims: the flow is frozen (forced convection whose velocity does not
//! depend on temperature; for buoyant flows use the steady
//! `natural_convection`); first-order in time (the step size is the
//! caller's accuracy decision; halving it is the convergence evidence); no
//! temperature-dependent properties.

use fs_exec::CancelGate;

use super::domain::{FluidProperties, SolidMaterial, Voxel, VoxelDomain};
use super::energy::{EnergyConfig, PseudoStep, ThermalSetup, solve_energy_inner};
use super::flow::FlowField;
use super::{ChtError, finite, finite_positive, poll};

/// Time-march configuration.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TransientConfig {
    /// Backward-Euler step, s.
    pub time_step_s: f64,
    /// Number of steps.
    pub steps: usize,
    /// Energy solver settings.
    pub energy: EnergyConfig,
}

/// Evidence of one time step.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TransientRecord {
    /// Time at the end of the step, s.
    pub time_s: f64,
    /// Largest cell temperature, K.
    pub max_temperature_k: f64,
    /// Largest solid-cell temperature, K (NaN without solids).
    pub max_solid_temperature_k: f64,
    /// `sum C (T^(n+1) - T^n)`, J.
    pub stored_energy_change_j: f64,
    /// `dt * sources`, J.
    pub source_j: f64,
    /// `dt * net boundary outflow`, J.
    pub boundary_outflow_j: f64,
    /// `stored - (source - outflow)`, J: the step's solver residual.
    pub closure_j: f64,
    /// Krylov iterations of the step.
    pub iterations: usize,
}

/// Final temperature and the step records.
#[derive(Debug, Clone, PartialEq)]
pub struct TransientSolution {
    /// Temperature per cell at the final time, K.
    pub temperature: Vec<f64>,
    /// One record per step.
    pub records: Vec<TransientRecord>,
}

/// March the conjugate energy equation from `initial_temperature`.
///
/// # Errors
/// Input refusals (including a solid without heat capacity), solver
/// refusals, or [`ChtError::Cancelled`].
#[allow(clippy::too_many_arguments)] // physics inputs + schedule + controls
pub fn march_energy(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    solids: &[SolidMaterial],
    flow: &FlowField,
    setup: &ThermalSetup,
    initial_temperature: &[f64],
    power_schedule: impl Fn(f64) -> f64,
    config: &TransientConfig,
    gate: &CancelGate,
) -> Result<TransientSolution, ChtError> {
    fluid.validate()?;
    finite_positive("transient.time_step_s", config.time_step_s)?;
    if !setup.compact_components.is_empty() {
        // A two-resistor model has no heat capacity: it is a steady model.
        return Err(ChtError::InvalidInput {
            field: "thermal.compact_components",
            reason: "two-resistor compact models are steady; the transient march refuses them"
                .into(),
        });
    }
    let cells = domain.cell_count();
    if initial_temperature.len() != cells {
        return Err(ChtError::InvalidInput {
            field: "transient.initial_temperature",
            reason: format!(
                "expected {cells} entries, got {}",
                initial_temperature.len()
            ),
        });
    }
    for &t in initial_temperature {
        finite("transient.initial_temperature", t)?;
    }
    domain.check_materials(solids.len())?;
    let capacity = heat_capacity(domain, fluid, solids)?;
    let storage: Vec<f64> = capacity.iter().map(|c| c / config.time_step_s).collect();
    let solids_present = (0..cells).any(|c| !domain.is_fluid(c));
    let mut temperature = initial_temperature.to_vec();
    let mut records = Vec::with_capacity(config.steps);
    let mut stepped = setup.clone();
    for step in 1..=config.steps {
        poll(gate)?;
        let time_s = step as f64 * config.time_step_s;
        let scale = power_schedule(time_s);
        finite("transient.power_schedule", scale)?;
        if !setup.power_w.is_empty() {
            for (out, base) in stepped.power_w.iter_mut().zip(&setup.power_w) {
                *out = base * scale;
            }
        }
        let (next, record) = energy_step(
            domain,
            fluid,
            solids,
            flow,
            &stepped,
            &config.energy,
            (&capacity, &storage),
            &temperature,
            (time_s, config.time_step_s, solids_present),
            gate,
        )?;
        records.push(record);
        temperature = next;
    }
    Ok(TransientSolution {
        temperature,
        records,
    })
}

/// Heat capacity `(rho c)_c V_c` of every cell, J/K (a solid without a
/// declared heat capacity refuses).
pub(crate) fn heat_capacity(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    solids: &[SolidMaterial],
) -> Result<Vec<f64>, ChtError> {
    let volume = domain.dx() * domain.dx() * domain.dx();
    (0..domain.cell_count())
        .map(|c| {
            let rho_c = match domain.voxel_at(c) {
                Voxel::Fluid => fluid.volumetric_heat_capacity(),
                Voxel::Solid(m) => {
                    let solid = &solids[usize::from(m)];
                    let value = solid.volumetric_heat_capacity_j_m3_k.ok_or_else(|| {
                        ChtError::InvalidInput {
                            field: "solid.volumetric_heat_capacity_j_m3_k",
                            reason: format!("solid `{}` declares no heat capacity", solid.label),
                        }
                    })?;
                    finite_positive("solid.volumetric_heat_capacity_j_m3_k", value)?;
                    value
                }
            };
            Ok(rho_c * volume)
        })
        .collect()
}

/// One backward-Euler energy step on `flow` from `previous`, with its
/// closure record. `(capacity, storage)` are `C` and `C / dt` per cell;
/// `(time_s, dt, solids_present)` label the record.
#[allow(clippy::too_many_arguments)] // physics inputs + step state
pub(crate) fn energy_step(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    solids: &[SolidMaterial],
    flow: &FlowField,
    setup: &ThermalSetup,
    energy: &EnergyConfig,
    (capacity, storage): (&[f64], &[f64]),
    previous: &[f64],
    (time_s, dt, solids_present): (f64, f64, bool),
    gate: &CancelGate,
) -> Result<(Vec<f64>, TransientRecord), ChtError> {
    let next = solve_energy_inner(
        domain,
        fluid,
        solids,
        flow,
        setup,
        energy,
        Some(&PseudoStep {
            coefficient_w_k: storage,
            previous,
        }),
        gate,
    )?;
    let stored: f64 = capacity
        .iter()
        .zip(next.temperature.iter().zip(previous))
        .map(|(c, (a, b))| c * (a - b))
        .sum();
    let balance = next.report.balance;
    let source_j = balance.source_w * dt;
    let boundary_outflow_j =
        (balance.boundary_outflow_w + balance.sink_outflow_w - balance.fixed_cell_injection_w) * dt;
    let max_temperature_k = next
        .temperature
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let max_solid_temperature_k = if solids_present {
        next.max_where(|c| !domain.is_fluid(c))
            .map_or(f64::NAN, |(_, t)| t)
    } else {
        f64::NAN
    };
    let record = TransientRecord {
        time_s,
        max_temperature_k,
        max_solid_temperature_k,
        stored_energy_change_j: stored,
        source_j,
        boundary_outflow_j,
        closure_j: stored - (source_j - boundary_outflow_j),
        iterations: next.report.iterations,
    };
    Ok((next.temperature, record))
}
