//! Cell-centered gray radiation on the original tensor/caloric timeline.
//!
//! The existing ray geometry and radiosity owner supplies both the tangent
//! sinks and the physical radiative powers. Every retry solves the same
//! original storage equation; a radiation iterate is never a time advance.

use fs_exec::CancelGate;

use super::super::coordinate::CellCoordinate;
use super::super::sdirk::remaining_controls;
use super::{
    CaloricStage, PreparedCaloric, TensorCaloricEndpoint, TensorCaloricStep, TensorCaloricStorage,
    TensorConductivityLaw, TensorNonlinearConfig, end_time, solve_stage,
};
use crate::conjugate::radiation::{RadiationConfig, Radiators};
use crate::conjugate::tensor_conduction::invalid;
use crate::conjugate::transient::record_conduction_step;
use crate::conjugate::{
    ChtError, SolidMaterial, ThermalSetup, VoxelDomain, finite, finite_positive, poll,
};

/// Original radiative evidence at one accepted caloric stage.
///
/// The native tensor endpoint retains its last tangent-sink equations and
/// their linear/nonlinear diagnostics. The accompanying transient record
/// instead uses the physical radiation, including this explicitly reported
/// difference. Ray factors remain estimates, not geometric certificates.
#[derive(Clone, Debug, PartialEq)]
pub struct TensorCaloricRadiationReport {
    /// Emitting faces retained by the original escape/enclosure geometry.
    pub exposed_faces: usize,
    /// Original rays traced once for the complete physical interval.
    pub rays: usize,
    /// Original exchange patches; zero in escape-only mode.
    pub patches: usize,
    /// Native heat solves used by this stage, including zero-correction solves.
    pub iterations: usize,
    /// Physical power escaping to the declared surroundings at this stage [W].
    /// Surface-to-surface transfers cancel within the enclosing system.
    pub radiated_w: f64,
    /// Last iterate's temperature change divided by the temperature span,
    /// with the original radiation owner's one-Kelvin minimum scale.
    pub temperature_change: f64,
    /// Largest absolute physical-minus-tangent cell radiation power [W].
    pub maximum_cell_power_defect_w: f64,
    /// Largest original roundoff-aware, heat-scaled power discrepancy.
    pub relative_power_defect: f64,
    /// Physical-minus-tangent outflow summed over free cells [W].
    /// A fixed cell's radiative correction changes its thermostat reaction
    /// by the same amount, so it cancels from the net heat account.
    pub net_outflow_correction_w: f64,
}

/// Advance one original-caloric backward-Euler interval with gray radiation.
///
/// The existing radiation owner evaluates emission at each emitting solid
/// voxel's CELL temperature, not at the tensor face traces. Its transparent
/// gas, gray/diffuse/opaque surface, patch averaging, ray geometry and source
/// limitations are unchanged. Conduction and radiation are distinct physical
/// routes on those solid/fluid interfaces. Authored cell sinks must not
/// already contain the same radiation; no duplicate route is inferred.
///
/// All radiation retries retain the original mass and specific energy at
/// the physical start. The original tensor owner reevaluates calorics,
/// conductivity integrals and contacts at each trial. A latent plateau
/// therefore remains an energy interval with zero dT/de.
///
/// The nonlinear, Krylov and backtracking budgets are shared across EVERY
/// radiative retry. RadiationConfig::max_iterations bounds the total native
/// heat solves. Both the temperature update and the original nonlinear
/// radiative powers must pass RadiationConfig::tolerance. The returned heat
/// record uses physical radiation, not only the last tangent sink. Source
/// caloric/conductivity budgets remain with those same original law objects.
///
/// # Errors
/// Original caloric/tensor/radiation admissions apply. A failed source query,
/// exhausted work budget, nonconvergence or cancellation returns no step.
#[allow(clippy::too_many_arguments)]
pub fn advance_caloric_tensor_conduction_radiating(
    domain: &VoxelDomain,
    solids: &[SolidMaterial],
    setup: &ThermalSetup,
    conductivity_laws: &[&dyn TensorConductivityLaw],
    storage: TensorCaloricStorage<'_>,
    prescribed_specific_energy_j_kg: &[(usize, f64)],
    start_time_s: f64,
    time_step_s: f64,
    config: &TensorNonlinearConfig,
    radiation: &RadiationConfig,
    gate: &CancelGate,
) -> Result<(TensorCaloricStep, TensorCaloricRadiationReport), ChtError> {
    poll(gate)?;
    domain.check_materials(solids.len())?;
    let time_s = end_time(start_time_s, time_step_s)?;
    let mut work = RadiationWork::new(config, radiation)?;
    let prepared =
        PreparedCaloric::new(domain, conductivity_laws.len(), storage, time_step_s, gate)?;
    let stage = prepared.stage(setup, prescribed_specific_energy_j_kg, gate)?;
    let radiators = Radiators::build(domain, solids, radiation, gate)?;
    let accepted = solve_radiating_stage(
        domain,
        solids,
        setup,
        conductivity_laws,
        &stage,
        storage.previous_specific_energy_j_kg,
        &vec![0.0; domain.cell_count()],
        None,
        &radiators,
        radiation,
        &mut work,
        gate,
    )?;
    let record = physical_stage_record(domain, setup, &prepared, &accepted, time_s, time_step_s)?;
    poll(gate)?;
    Ok((
        TensorCaloricStep {
            endpoint: accepted.endpoint,
            record,
        },
        accepted.report,
    ))
}

/// Work belongs to the whole interval, never to an individual radiation retry.
struct RadiationWork {
    remaining: TensorNonlinearConfig,
    heat_solves: usize,
    maximum_heat_solves: usize,
}

impl RadiationWork {
    fn new(config: &TensorNonlinearConfig, radiation: &RadiationConfig) -> Result<Self, ChtError> {
        finite_positive("tensor.caloric.radiation.tolerance", radiation.tolerance)?;
        if radiation.max_iterations == 0 {
            return Err(invalid(
                "caloric radiation needs a positive heat-solve budget",
            ));
        }
        Ok(Self {
            remaining: *config,
            heat_solves: 0,
            maximum_heat_solves: radiation.max_iterations,
        })
    }

    fn begin_solve(&mut self, residual: f64, tolerance: f64) -> Result<(), ChtError> {
        if self.heat_solves == self.maximum_heat_solves {
            return Err(ChtError::SolverNotConverged {
                system: "caloric tensor radiation",
                iterations: self.heat_solves,
                relative_residual: residual,
                tolerance,
            });
        }
        self.heat_solves += 1;
        Ok(())
    }
}

struct RadiatingStage {
    endpoint: TensorCaloricEndpoint,
    report: TensorCaloricRadiationReport,
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn solve_radiating_stage(
    domain: &VoxelDomain,
    solids: &[SolidMaterial],
    setup: &ThermalSetup,
    conductivity_laws: &[&dyn TensorConductivityLaw],
    stage: &CaloricStage<'_, '_>,
    previous_specific_energy_j_kg: &[f64],
    initial: &[f64],
    carried_power_w: Option<&[f64]>,
    radiators: &Radiators,
    radiation: &RadiationConfig,
    work: &mut RadiationWork,
    gate: &CancelGate,
) -> Result<RadiatingStage, ChtError> {
    let cells = domain.cell_count();
    if initial.len() != cells || previous_specific_energy_j_kg.len() != cells {
        return Err(invalid(
            "caloric radiation needs one original energy and seed per cell",
        ));
    }
    let mut seed = initial.to_vec();
    let mut iterate = Vec::with_capacity(cells);
    for (cell, &value) in seed.iter().enumerate() {
        if cell.is_multiple_of(256) {
            poll(gate)?;
        }
        iterate.push(stage.temperature_and_slope(cell, value)?[0]);
    }
    let mut tangent = radiators.endpoint(domain, &iterate, gate)?;
    let mut step = setup.clone();
    let original_sinks = step.cell_sinks.len();
    let mut fixed = vec![false; cells];
    for &(cell, _) in &setup.fixed_temperature {
        let fixed = fixed.get_mut(cell).ok_or_else(|| {
            invalid("caloric radiation fixed cell is outside the original domain")
        })?;
        *fixed = true;
    }
    let mut tangent_power = vec![0.0; cells];
    let mut tangent_scale = vec![0.0; cells];
    let mut stage_iterations = 0_usize;
    let mut stage_newton = 0_usize;
    let mut stage_backtracks = 0_usize;
    let mut solves = 0_usize;
    let mut residual = f64::INFINITY;
    loop {
        poll(gate)?;
        work.begin_solve(residual, radiation.tolerance)?;
        step.cell_sinks.truncate(original_sinks);
        step.cell_sinks.extend_from_slice(&tangent.sinks);
        let mut next = solve_stage(
            domain,
            solids,
            &step,
            conductivity_laws,
            stage,
            &seed,
            carried_power_w,
            &work.remaining,
            gate,
        )?;
        work.remaining = remaining_controls(&work.remaining, &next.solution)?;
        stage_iterations = stage_iterations
            .checked_add(next.solution.solution.report.iterations)
            .ok_or_else(|| invalid("caloric radiation Krylov work count overflowed"))?;
        stage_newton = stage_newton
            .checked_add(next.solution.newton_iterations)
            .ok_or_else(|| invalid("caloric radiation Newton work count overflowed"))?;
        stage_backtracks = stage_backtracks
            .checked_add(next.solution.backtracks)
            .ok_or_else(|| invalid("caloric radiation backtrack work count overflowed"))?;
        solves += 1;

        let temperature = &next.solution.solution.temperature;
        let physical = radiators.endpoint(domain, temperature, gate)?;
        tangent.write_tangent_power(temperature, &mut tangent_power, &mut tangent_scale, gate)?;
        let mut low = f64::INFINITY;
        let mut high = f64::NEG_INFINITY;
        let mut change = 0.0_f64;
        let mut power_residual = 0.0_f64;
        let mut maximum_defect = 0.0_f64;
        let mut correction = 0.0;
        for cell in 0..cells {
            if cell.is_multiple_of(512) {
                poll(gate)?;
            }
            let t = temperature[cell];
            low = low.min(t);
            high = high.max(t);
            change = change.max((t - iterate[cell]).abs());
            let storage_rate = stage
                .storage
                .storage_rate_w(cell, next.specific_energy_j_kg[cell]);
            let scale = storage_rate
                .abs()
                .max(setup.power_w.get(cell).copied().unwrap_or(0.0).abs())
                .max(carried_power_w.map_or(0.0, |carry| carry[cell].abs()));
            let defect =
                physical.power_defect(cell, tangent_power[cell], tangent_scale[cell], scale)?;
            power_residual = power_residual.max(defect.relative);
            maximum_defect = maximum_defect.max(defect.difference_w.abs());
            if !fixed[cell] {
                correction += defect.difference_w;
            }
        }
        finite("tensor.caloric.radiation.outflow_correction_w", correction)?;
        let temperature_change =
            (change - 64.0 * f64::EPSILON * high.abs()).max(0.0) / (high - low).max(1.0);
        finite(
            "tensor.caloric.radiation.temperature_change",
            temperature_change,
        )?;
        residual = temperature_change.max(power_residual);
        if residual <= radiation.tolerance {
            // Retain total stage work, not only its final successful tangent.
            next.solution.solution.report.iterations = stage_iterations;
            next.solution.newton_iterations = stage_newton;
            next.solution.backtracks = stage_backtracks;
            let (exposed_faces, rays, patches) = radiators.counts(radiation.rays_per_face);
            let report = TensorCaloricRadiationReport {
                exposed_faces,
                rays,
                patches,
                iterations: solves,
                radiated_w: physical.radiated_w,
                temperature_change,
                maximum_cell_power_defect_w: maximum_defect,
                relative_power_defect: power_residual,
                net_outflow_correction_w: correction,
            };
            poll(gate)?;
            return Ok(RadiatingStage {
                endpoint: next,
                report,
            });
        }
        // Only the numerical seed changes. The stage still borrows the SAME
        // physical start energy and mass, including on an isothermal plateau.
        for cell in 0..cells {
            if cell.is_multiple_of(512) {
                poll(gate)?;
            }
            seed[cell] = next.specific_energy_j_kg[cell] - previous_specific_energy_j_kg[cell];
            finite("tensor.caloric.radiation.energy_seed_j_kg", seed[cell])?;
        }
        iterate.copy_from_slice(temperature);
        tangent = physical;
    }
}

fn physical_stage_record(
    domain: &VoxelDomain,
    setup: &ThermalSetup,
    storage: &PreparedCaloric<'_>,
    accepted: &RadiatingStage,
    time_s: f64,
    time_step_s: f64,
) -> Result<crate::conjugate::transient::TransientRecord, ChtError> {
    let endpoint = &accepted.endpoint;
    record_conduction_step(
        domain,
        setup,
        &endpoint.solution.solution.temperature,
        endpoint.solution.solution.report.balance,
        endpoint.solution.solution.report.iterations,
        |cell| storage.stored_j(cell, endpoint.specific_energy_j_kg[cell]),
        (
            time_s,
            time_step_s,
            (0..domain.cell_count()).any(|cell| !domain.is_fluid(cell)),
        ),
        accepted.report.radiated_w,
        accepted.report.net_outflow_correction_w,
    )
}
