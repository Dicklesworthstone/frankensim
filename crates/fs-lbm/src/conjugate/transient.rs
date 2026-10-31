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
//! Radiation, when declared, is converged at the backward-Euler endpoint:
//! both emission and surface irradiation use the new temperature. Every
//! nonlinear trial retains the original step's storage state. The closure
//! uses the physical face law, and an unconverged endpoint refuses.
//!
//! No-claims: the flow is frozen (forced convection whose velocity does not
//! depend on temperature; buoyant transients march through
//! `march_conjugate`); first-order in time (the step size is the caller's
//! accuracy decision; halving it is the convergence evidence).

use fs_exec::CancelGate;

use super::domain::{FluidProperties, SolidMaterial, Voxel, VoxelDomain};
use super::energy::{EnergyConfig, EnergySolution, PseudoStep, ThermalSetup, solve_energy_inner};
use super::flow::FlowField;
use super::radiation::{RadiationConfig, Radiators};
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
    /// Stored-energy increment, J. Ordinary constant-storage steps use
    /// `sum C (T^(n+1) - T^n)`; sourced solid coupling uses its exact Cp integral.
    pub stored_energy_change_j: f64,
    /// `dt * sources`, J.
    pub source_j: f64,
    /// `dt * net boundary outflow`, J.
    pub boundary_outflow_j: f64,
    /// `stored - (source - outflow)`, J: the step's solver residual.
    pub closure_j: f64,
    /// Net heat radiated to the surroundings at the step's end temperature,
    /// W (zero without radiation; radiation's share of the outflow).
    pub radiated_w: f64,
    /// Total Krylov iterations over all nonlinear solves of the step.
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
    radiation: Option<&RadiationConfig>,
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
    let radiators = match radiation {
        Some(r) => Some(Radiators::build(domain, solids, r, gate)?),
        None => None,
    };
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
            radiators.as_ref().zip(radiation),
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
            Ok(rho_c * domain.volume(c))
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
    radiation: Option<(&Radiators, &RadiationConfig)>,
    gate: &CancelGate,
) -> Result<(Vec<f64>, TransientRecord), ChtError> {
    energy_step_with_conductivity(
        domain,
        fluid,
        solids,
        flow,
        setup,
        energy,
        (capacity, storage),
        previous,
        (time_s, dt, solids_present),
        radiation,
        None,
        gate,
    )
}

/// The same storage, radiation and closure owner with an optional admitted
/// source-material conductivity field frozen at the preceding temperature.
#[allow(clippy::too_many_arguments)]
pub(crate) fn energy_step_with_conductivity(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    solids: &[SolidMaterial],
    flow: &FlowField,
    setup: &ThermalSetup,
    energy: &EnergyConfig,
    (capacity, storage): (&[f64], &[f64]),
    previous: &[f64],
    (time_s, dt, solids_present): (f64, f64, bool),
    radiation: Option<(&Radiators, &RadiationConfig)>,
    molecular_w_m_k: Option<&[[f64; 3]]>,
    gate: &CancelGate,
) -> Result<(Vec<f64>, TransientRecord), ChtError> {
    let (next, radiated_w, radiation_outflow_correction_w) = solve_step(
        domain, fluid, solids, flow, setup, energy, storage, previous, radiation,
        molecular_w_m_k, gate,
    )?;
    let stored: f64 = capacity
        .iter()
        .zip(next.temperature.iter().zip(previous))
        .map(|(c, (a, b))| c * (a - b))
        .sum();
    let balance = next.report.balance;
    let source_j = balance.source_w * dt;
    // Prescribed cells may jump from the initial state. Their reservoir
    // supplies that storage change as well as the steady face/sink reaction.
    // Repeated declarations use the solver's last value and count only once.
    let mut fixed_cells: Vec<usize> = setup.fixed_temperature.iter().map(|&(cell, _)| cell).collect();
    fixed_cells.sort_unstable();
    fixed_cells.dedup();
    let fixed_storage_j: f64 = fixed_cells
        .iter()
        .map(|&cell| capacity[cell] * (next.temperature[cell] - previous[cell]))
        .sum();
    let boundary_outflow_j = (balance.boundary_outflow_w + balance.sink_outflow_w
        - balance.fixed_cell_injection_w
        + radiation_outflow_correction_w)
        * dt
        - fixed_storage_j;
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
        radiated_w,
        iterations: next.report.iterations,
    };
    Ok((next.temperature, record))
}

/// Solve one endpoint, returning its external radiation and the correction
/// from tangent sinks to the physical radiation law in the net boundary
/// balance. The non-radiating path remains one unchanged energy solve.
#[allow(clippy::too_many_arguments)]
fn solve_step(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    solids: &[SolidMaterial],
    flow: &FlowField,
    setup: &ThermalSetup,
    energy: &EnergyConfig,
    storage: &[f64],
    previous: &[f64],
    radiation: Option<(&Radiators, &RadiationConfig)>,
    molecular_w_m_k: Option<&[[f64; 3]]>,
    gate: &CancelGate,
) -> Result<(EnergySolution, f64, f64), ChtError> {
    let solve = |setup: &ThermalSetup| {
        // This is the time level, not the previous radiation iterate.
        let step = PseudoStep {
            coefficient_w_k: storage,
            previous,
        };
        match molecular_w_m_k {
            Some(field) => energy.solve_step_with_conductivity(
                domain, fluid, solids, flow, setup, field, &step, gate,
            ),
            None => solve_energy_inner(
                domain,
                fluid,
                solids,
                flow,
                setup,
                energy,
                Some(&step),
                gate,
            ),
        }
    };
    let Some((radiators, controls)) = radiation else {
        return Ok((solve(setup)?, 0.0, 0.0));
    };
    finite_positive("radiation.tolerance", controls.tolerance)?;
    if controls.max_iterations == 0 {
        return Err(ChtError::InvalidInput {
            field: "radiation.max_iterations",
            reason: "must be at least one".into(),
        });
    }
    let cells = domain.cell_count();
    let mut fixed = vec![false; cells];
    for &(cell, _) in &setup.fixed_temperature {
        // The energy solver supplies the input refusal for invalid indices.
        if let Some(value) = fixed.get_mut(cell) {
            *value = true;
        }
    }
    let mut iterate = previous.to_vec();
    let mut endpoint = radiators.endpoint(domain, &iterate, gate)?;
    let mut step = setup.clone();
    let base_sinks = setup.cell_sinks.len();
    let mut tangent_power = vec![0.0; cells];
    let mut tangent_scale = vec![0.0; cells];
    let mut total_iterations = 0usize;
    let mut relative_residual = f64::INFINITY;
    for _ in 0..controls.max_iterations {
        poll(gate)?;
        step.cell_sinks.truncate(base_sinks);
        step.cell_sinks.extend_from_slice(&endpoint.sinks);
        let mut next = solve(&step)?;
        total_iterations = total_iterations.saturating_add(next.report.iterations);
        let physical = radiators.endpoint(domain, &next.temperature, gate)?;
        tangent_power.fill(0.0);
        tangent_scale.fill(0.0);
        for sink in &endpoint.sinks {
            tangent_power[sink.cell] +=
                sink.conductance_w_k * (next.temperature[sink.cell] - sink.temperature_k);
            tangent_scale[sink.cell] += sink.conductance_w_k
                * (next.temperature[sink.cell].abs() + sink.temperature_k.abs());
        }
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        let mut change = 0.0f64;
        let mut power_residual = 0.0f64;
        let mut correction_w = 0.0;
        for cell in 0..cells {
            if cell % 1024 == 0 {
                poll(gate)?;
            }
            let t = next.temperature[cell];
            lo = lo.min(t);
            hi = hi.max(t);
            change = change.max((t - iterate[cell]).abs());
            let defect = physical.cell_power_w[cell] - tangent_power[cell];
            finite("radiation.power_defect_w", defect)?;
            let scale = physical.cell_power_w[cell]
                .abs()
                .max(tangent_power[cell].abs())
                .max(storage[cell] * (t - previous[cell]).abs())
                .max(setup.power_w.get(cell).copied().unwrap_or(0.0).abs());
            // Near radiative equilibrium, emission and absorption cancel;
            // distinguish that roundoff from a real local power imbalance.
            finite("radiation.tangent_scale_w", tangent_scale[cell])?;
            finite("radiation.power_scale_w", scale)?;
            let roundoff = 64.0 * f64::EPSILON * physical.cell_scale_w[cell]
                + 64.0 * f64::EPSILON * tangent_scale[cell];
            finite("radiation.power_roundoff_w", roundoff)?;
            power_residual = power_residual.max(
                (defect.abs() - roundoff).max(0.0) / scale.max(f64::MIN_POSITIVE),
            );
            if !fixed[cell] {
                correction_w += defect;
            }
        }
        finite("radiation.outflow_correction_w", correction_w)?;
        let temperature_residual =
            (change - 64.0 * f64::EPSILON * hi.abs()).max(0.0) / (hi - lo).max(1.0);
        relative_residual = temperature_residual.max(power_residual);
        if relative_residual <= controls.tolerance {
            // A fixed cell's physical sink correction also corrects the
            // reservoir injection by the same amount. Only free-cell
            // corrections therefore enter the net boundary balance.
            next.report.iterations = total_iterations;
            return Ok((next, physical.radiated_w, correction_w));
        }
        iterate = next.temperature;
        endpoint = physical;
    }
    Err(ChtError::SolverNotConverged {
        system: "transient radiation",
        iterations: controls.max_iterations,
        relative_residual,
        tolerance: controls.tolerance,
    })
}
