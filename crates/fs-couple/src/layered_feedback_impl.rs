//! Nonlinear solve for the public layered thermo-mechanical coupling.
use std::cell::{Cell,RefCell};
use std::collections::{BTreeMap,BTreeSet};
use fs_material::thermovisco::thermal::CaloricPoint;
use fs_solid::deposition::{self,AxialControl,BendingControl,Control,LayeredStrip};
use fs_solver::{NewtonKrylovState,SolverRunProgress};
use fs_solver::finite_difference::{FallibleFiniteDifference,FiniteDifferenceError};
use fs_thermal::enthalpy::{DepositionError,DepositionSimulation};
use fs_thermal::enthalpy::backward_euler::{BackwardEulerStep,StepCoefficients};
use super::types::*;

fn finite(value:f64)->Result<f64,DepositionError>{
    if value.is_finite(){Ok(value)}else{Err(DepositionError::Arithmetic)}
}
fn same_positive(a:f64,b:f64)->bool{
    a.is_finite()&&b.is_finite()&&a>0.0&&b>0.0
        &&(a-b).abs()<=64.0*f64::EPSILON*a.abs().max(b.abs())
}

#[allow(clippy::too_many_arguments)]
pub fn step_candidate<'t, 'm, S, E>(
    thermal: &DepositionSimulation<'t>, strip: &LayeredStrip<'m>,
    bindings: &[LayerThermalBinding<'_>], target_s: f64, control: Control,
    config: StripEnthalpyConfig, initial_energy_increments_j: &[f64],
    resolve: impl Fn(&DepositionSimulation<'t>) -> Result<ThermalCoefficients<S>, E>,
    cancelled: impl FnMut() -> bool,
) -> Result<StripEnthalpyCandidate<'t, 'm, S>, StripEnthalpyError<E>> {
    let numerical = |e:FiniteDifferenceError<StripEnthalpyError<E>>| match e {
        FiniteDifferenceError::Source(e)=>e,
        FiniteDifferenceError::Invalid(s)=>StripEnthalpyError::Invalid(s),
    };
    let cancel = RefCell::new(cancelled);
    let checkpoint = || if (*cancel.borrow_mut())() { Err(DepositionError::Cancelled) } else { Ok(()) };
    checkpoint()?;
    config.newton.validate()?;
    if strip.time_s() != thermal.time_s() || strip.layers().is_empty()
        || bindings.len() != strip.layers().len() || config.maximum_cells == 0
        || thermal.cells().len() > config.maximum_cells
        || initial_energy_increments_j.len() != thermal.cells().len()
        || config.maximum_structural_trials == 0 || config.maximum_newton_iterations == 0
        || [config.cell_energy_tolerance_j, config.thermal_energy_tolerance_j, config.total_energy_tolerance_j]
            .iter().any(|v| !v.is_finite() || *v < 0.0) {
        return Err(StripEnthalpyError::Invalid("complete active bindings, common clocks and valid work/energy controls required"));
    }
    let mut layers = BTreeMap::new();
    let mut cells = BTreeSet::new();
    for binding in bindings {
        checkpoint()?;
        let layer = strip.layers().get(&binding.layer)
            .ok_or(StripEnthalpyError::Invalid("unknown mechanical layer"))?;
        let cell = thermal.cells().get(&binding.cell)
            .ok_or(StripEnthalpyError::Invalid("inactive thermal destination"))?;
        if layers.insert(binding.layer, binding).is_some() || !cells.insert(binding.cell)
            || !thermal.uses_material_law(cell.material_id, binding.caloric.curve()) {
            return Err(StripEnthalpyError::Invalid("duplicate layer/cell or different caloric object"));
        }
        let density = binding.material.reference_density_kg_m3(binding.caloric)
            .map_err(StripEnthalpyError::Binding)?;
        let volume = finite(layer.geometry.area_m2 * strip.length_m())?;
        if !same_positive(volume, cell.geometry.volume_m3)
            || !same_positive(finite(volume*density)?, cell.mass_kg) {
            return Err(StripEnthalpyError::Invalid("thermal mass/volume differ from the source-bound reference layer"));
        }
        for point in layer.material_points() {
            if !std::ptr::eq(point.law(), binding.material.law()) {
                return Err(StripEnthalpyError::Invalid("layer does not retain the declared mechanical source law"));
            }
            CaloricPoint::try_new(point.clone(), binding.caloric.curve(), cell.specific_h_j_kg,
                cell.mass_kg, volume, cell.geometry.surface_area_m2)
                .map_err(StripEnthalpyError::CaloricAdmission)?;
        }
    }
    let step = BackwardEulerStep::new(thermal, target_s, checkpoint)?;
    let dt = step.duration_s();
    let trials = Cell::new(0usize);
    let evaluate = |x: &[f64]| step.evaluate(x, |endpoint| {
        checkpoint()?;
        if trials.get() >= config.maximum_structural_trials { return Err(StripEnthalpyError::StructuralBudget); }
        trials.set(trials.get()+1);
        let mut temperatures = BTreeMap::new();
        for (&id, binding) in &layers {
            checkpoint()?;
            temperatures.insert(id, endpoint.temperature_k(binding.cell)?);
        }
        let (next, mechanical, layer_energy) = strip.advance_with_layer_energy(target_s,
            &temperatures, control, config.mechanics,
            || checkpoint().map_err(|_| deposition::Error::Cancelled))?;
        let coefficients = resolve(endpoint).map_err(StripEnthalpyError::Source)?;
        checkpoint()?;
        let other_heat = finite(dt * endpoint.external_step_power_w(&coefficients.sources_w,
            &coefficients.conductances_w_k, checkpoint)?)?;
        let mut sources = coefficients.sources_w;
        let mut deposited = 0.0;
        for (&id, energy) in &layer_energy {
            checkpoint()?;
            if !energy.dissipated.is_finite() || energy.dissipated < 0.0 {
                return Err(StripEnthalpyError::Invalid("nonphysical layer dissipation"));
            }
            let power = finite(energy.dissipated/dt)?;
            if energy.dissipated != 0.0 && power == 0.0 {
                return Err(StripEnthalpyError::Invalid("layer heat underflowed in power conversion"));
            }
            deposited = finite(deposited + finite(dt*power)?)?;
            if power != 0.0 {
                let total = sources.entry(layers[&id].cell).or_insert(0.0);
                *total = finite(*total + power)?;
            }
        }
        if (deposited-mechanical.energy_j.dissipated).abs() > config.total_energy_tolerance_j {
            return Err(StripEnthalpyError::Invalid("layer heat differs from the complete material account"));
        }
        Ok::<_, StripEnthalpyError<E>>(StepCoefficients { sources_w: sources,
            conductances_w_k: coefficients.conductances_w_k,
            payload: (next, mechanical, layer_energy, coefficients.source, deposited, other_heat) })
    }, checkpoint);
    let problem = FallibleFiniteDifference::new(step.cell_ids().len(), config.difference_step_j,
        |x: &[f64]| evaluate(x).map(|trial| trial.residual_j().to_vec())).map_err(numerical)?;
    let construction = NewtonKrylovState::new(&problem, initial_energy_increments_j.to_vec(), config.newton);
    if let Some(e) = problem.take_failure() { return Err(numerical(e)); }
    let mut solver = construction?;
    let result = solver.run_controlled(&problem, config.maximum_newton_iterations,
        &mut || (*cancel.borrow_mut())() || problem.has_failed());
    if let Some(e) = problem.take_failure() { return Err(numerical(e)); }
    checkpoint()?;
    if result.progress != SolverRunProgress::Complete { return Err(StripEnthalpyError::Interrupted(result.progress)); }
    if !result.report.converged { return Err(StripEnthalpyError::NotConverged(result.report)); }
    let trial = evaluate(&solver.x)?;
    let maximum_cell_residual_j = trial.maximum_cell_residual_j();
    let (thermal_next, (strip_next, mechanical, layer_energy_j, source, deposited_heat_j,
        non_mechanical_heat_j), audit) = trial.accept(config.cell_energy_tolerance_j,
            config.thermal_energy_tolerance_j, checkpoint)?;
    let force = match control.axial { AxialControl::Strain(_) => mechanical.axial_force_n, AxialControl::ForceN(n) => n };
    let moment = match control.bending { BendingControl::CurvaturePerM(_) => mechanical.bending_moment_nm, BendingControl::MomentNm(m) => m };
    let work = finite(strip.length_m() * finite(force * (strip_next.axial_strain()-strip.axial_strain())
        + moment * (strip_next.curvature_per_m()-strip.curvature_per_m()))?)?;
    let total = finite(finite(strip_next.stored_energy_j()-strip.stored_energy_j())?
        + audit.stored_energy_change_j + mechanical.energy_j.numerical_dissipation
        - work - non_mechanical_heat_j)?;
    if total.abs() > config.total_energy_tolerance_j
        || finite(audit.external_heat_j-deposited_heat_j-non_mechanical_heat_j)?.abs() > config.total_energy_tolerance_j {
        return Err(StripEnthalpyError::Invalid("thermal, mechanical and independently imposed work accounts do not close"));
    }
    checkpoint()?;
    Ok(StripEnthalpyCandidate { thermal: thermal_next, strip: strip_next,
        report: StripEnthalpyReport { mechanical, layer_energy_j, thermal: audit, source,
            nonlinear: result.report, structural_trials: trials.get(), maximum_cell_residual_j,
            deposited_heat_j, non_mechanical_heat_j, external_mechanical_work_j: work,
            total_energy_residual_j: total } })
}