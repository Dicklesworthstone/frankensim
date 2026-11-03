//! Native adhesive damage and viscosity with finite adhesive heat storage.
//!
//! One original sourced cohesive patch represents the through-thickness reduced
//! mechanics of one existing homogeneous adhesive cell. Only its surviving-bond
//! viscosity heats that cell. Fracture keeps its nonthermal surface account;
//! damaged-surface friction refuses because its location is not the layer bulk.
//! Actual layer mass, caloric law, geometry and two source contact resistances
//! remain with `ThirdBodyHeatStep`. No heat partition or extra capacity is added.
//!
//! This is a small-deformation, fixed material-reference, closed-normal-path
//! split with frozen mechanical coefficients. The declared source must describe
//! this layer thickness; matching IDs/areas cannot establish that calibration.
//! The original card's named third body is explicitly mapped to the thermal
//! layer material. Evolving thickness, open cracks and thermal weakening need
//! their own geometry/source transfer. The optional entropy account covers
//! caloric storage, conduction and viscous generation, not fracture/surface
//! entropy or a temperature-dependent cohesive free energy.

use super::{
    ThirdBodyHeatInput, ThirdBodyHeatReport, ThirdBodyHeatRoute, ThirdBodyHeatStep,
    ThirdBodyHeatTolerances,
};
use crate::friction_heat::cohesive::CohesiveHeatError;
use fs_exec::CancelGate;
use fs_matdb::{InterfaceSystemCard, MaterialStateId};
use fs_solid::cohesive::dynamics::{CohesiveMotion, DynamicsStep};
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress};
use fs_thermal::enthalpy::backward_euler::{StepCoefficients, resolvable_increment_j};
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation, MaterialId};
use std::cell::Cell;

/// Explicit ordered assignment of one original surface to its thermal body.
#[derive(Clone, Copy)]
pub struct CohesiveLayerHeatBody<'a> {
    pub body_id: u64,
    pub thermal_material: MaterialId,
    pub mechanical_material: &'a MaterialStateId,
}

/// One original joint and its actual finite adhesive material.
/// Body order is the original card's surface A/B order. The named third body
/// maps to `layer_material`; a bare material pair cannot invent an adhesive.
#[derive(Clone, Copy)]
pub struct CohesiveLayerHeatBinding<'a> {
    pub patch_id: u64,
    pub route: &'a ThirdBodyHeatRoute,
    pub card: &'a InterfaceSystemCard,
    pub bodies: [CohesiveLayerHeatBody<'a>; 2],
    pub layer_material: MaterialId,
    pub layer_name: &'a str,
}

/// Original endpoint transport. The caller cannot author the mechanical loss.
pub struct CohesiveLayerHeatSource<P> {
    /// Body A, body B, adhesive, in W/(m K).
    pub conductivity_w_m_k: [f64; 3],
    /// Body A/adhesive, body B/adhesive, in m² K/W.
    pub area_resistance_m2_k_w: [f64; 2],
    /// OTHER heat and links; both original layer contacts are already owned.
    pub transport: StepCoefficients<P>,
}

#[derive(Clone, Copy, Debug)]
pub struct CohesiveLayerHeatControls {
    pub maximum_cells: usize,
    pub maximum_trials: usize,
    pub maximum_newton_iterations: usize,
    pub difference_step_j: f64,
    pub heat: ThirdBodyHeatTolerances,
    pub combined_energy_tolerance_j: f64,
    pub newton: NewtonKrylovConfig,
    pub audit_entropy: bool,
}

/// Immutable native joint/thermal proposal; neither accepted owner is mutated.
pub struct HeatedCohesiveLayerCandidate<'s, 'a, 't, P> {
    pub mechanical: DynamicsStep<'s, 'a>,
    pub thermal: DepositionSimulation<'t>,
    pub heat: ThirdBodyHeatReport,
    pub source: P,
    pub nonlinear: NewtonReport,
    pub thermal_trials: usize,
    pub fracture_work_j: f64,
    /// Delta(K + bulk/bond storage + H) + fracture - applied work - other heat.
    pub combined_energy_residual_j: f64,
    pub summed_energy_residual_j: f64,
}

/// Advance one original finite adhesive joint and its surrounding heat network.
///
/// Every native patch must be covered: this rung admits exactly one patch and
/// one layer, with matching full layer cross-section and both contact areas.
/// The mechanical owner advances once; thermal trials reuse its actual local
/// viscosity loss and accepted damage. The transport callback may use that
/// damage and the trial layer temperature to query original material sources.
/// It must not publish its payload before this function returns success.
///
/// A nonzero damaged-surface friction loss refuses, even when small. It needs
/// separately located surface heat routes, not an invented fraction into this
/// adhesive volume. Intrinsic fracture, reversible stored energy and source
/// parameter work never enter the layer heat source. Mechanical temperature
/// transfer is not performed here. Ordinary outside heat remains separately
/// accounted for; numerical residuals never become material heat.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn propose_deformable_cohesive_layer_heat<'s, 'a, 't, P, E>(
    motion: &CohesiveMotion<'s, 'a>,
    mean_force_n: &[f64],
    endpoint_prescribed_jumps_m: &[[f64; 3]],
    dt_s: f64,
    binding: CohesiveLayerHeatBinding<'_>,
    thermal: &DepositionSimulation<'t>,
    controls: CohesiveLayerHeatControls,
    energy_guess_j: &[f64],
    resolve: impl Fn(
        &DepositionSimulation<'t>,
        f64,
        &DynamicsStep<'s, 'a>,
    ) -> Result<CohesiveLayerHeatSource<P>, E>,
    gate: &CancelGate,
) -> Result<HeatedCohesiveLayerCandidate<'s, 'a, 't, P>, CohesiveHeatError<E>> {
    let check = || {
        if gate.is_requested() {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    };
    check()?;
    controls.newton.validate()?;
    let cells = thermal.cells().len();
    let target = motion.time_s() + dt_s;
    let patches = motion.state().system().patches();
    if motion.time_s() != thermal.time_s()
        || !dt_s.is_finite()
        || dt_s <= 0.0
        || !target.is_finite()
        || target <= motion.time_s()
        || cells == 0
        || cells > controls.maximum_cells
        || controls.maximum_trials == 0
        || controls.maximum_newton_iterations == 0
        || energy_guess_j.len() != cells
        || energy_guess_j.iter().any(|v| !v.is_finite())
        || !controls.difference_step_j.is_finite()
        || controls.difference_step_j <= 0.0
        || [
            controls.heat.cell_j,
            controls.heat.thermal_j,
            controls.heat.layer_j,
            controls.heat.total_j,
            controls.combined_energy_tolerance_j,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v < 0.0)
        || patches.len() != 1
        || patches[0].id() != binding.patch_id
    {
        return Err(CohesiveHeatError::Invalid(
            "one complete sourced joint and bounded synchronized finite thermal state required",
        ));
    }
    // This also checks unchanged original masses, materials, geometry and links.
    let step = ThirdBodyHeatStep::new(thermal, target, binding.route, check)?;
    let source = patches[0].source().ok_or(CohesiveHeatError::Invalid(
        "finite adhesive heat requires the original sourced cohesive patch",
    ))?;
    let area = patches[0].area_m2();
    let layer = &thermal.cells()[&binding.route.layer_cell()];
    let layer_area = layer.geometry.volume_m3 / (2.0 * binding.route.distances_m()[0][1]);
    if source.resolved().card_identity() != binding.card.content_hash()
        || binding.card.third_body() != Some(binding.layer_name)
        || binding.layer_name.trim().is_empty()
        || layer.material_id != binding.layer_material
        || !same_area(area, layer_area)
        || binding
            .route
            .patches
            .iter()
            .any(|p| !same_area(area, p.link.area_m2))
        || motion.state().points()[0].jump_m()[0] > 0.0
    {
        return Err(CohesiveHeatError::Invalid(
            "original named adhesive, finite layer material, full patch areas and closed path required",
        ));
    }
    for (side, cell_id) in binding.route.body_cells().into_iter().enumerate() {
        let cell = &thermal.cells()[&cell_id];
        let body = binding.bodies[side];
        if cell.body_id != body.body_id
            || cell.material_id != body.thermal_material
            || &source.surfaces()[side].material != body.mechanical_material
        {
            return Err(CohesiveHeatError::Invalid(
                "original ordered cohesive surfaces and finite thermal body materials differ",
            ));
        }
    }
    let mechanical = motion.advance_controlled(
        mean_force_n,
        endpoint_prescribed_jumps_m,
        step.duration_s(),
        &mut || gate.is_requested(),
    )?;
    if mechanical.motion.time_s() != target
        || mechanical.patch_dissipation.len() != 1
        || mechanical.motion.state().points()[0].jump_m()[0] > 0.0
        || mechanical.friction_dissipated_j != 0.0
    {
        return Err(CohesiveHeatError::Invalid(
            "finite adhesive heat needs a closed native path without unlocated surface friction",
        ));
    }
    let loss = mechanical.patch_dissipation[0];
    if loss.patch_id != binding.patch_id
        || loss.friction_j != 0.0
        || [loss.viscous_j, loss.fracture_j]
            .iter()
            .any(|v| !v.is_finite() || *v < 0.0)
        || (loss.viscous_j - mechanical.viscous_dissipated_j).abs()
            > controls.combined_energy_tolerance_j
        || (loss.fracture_j - mechanical.cohesive_energy_j.dissipated).abs()
            > controls.combined_energy_tolerance_j
    {
        return Err(CohesiveHeatError::Invalid(
            "original adhesive loss accounts differ",
        ));
    }
    let trials = Cell::new(0usize);
    let evaluate = |x: &[f64], entropy: bool| {
        check()?;
        if trials.get() >= controls.maximum_trials {
            return Err(CohesiveHeatError::TrialBudget);
        }
        trials.set(trials.get() + 1);
        if x.len() != cells || x.iter().any(|v| !v.is_finite()) {
            return Err(CohesiveHeatError::Invalid(
                "complete finite layer heat iterate required",
            ));
        }
        let increments: Vec<_> = thermal
            .cells()
            .values()
            .zip(x)
            .map(|(cell, &q)| resolvable_increment_j(cell.specific_h_j_kg, q, cell.mass_kg))
            .collect();
        let endpoint = |state: &DepositionSimulation<'t>, temperature| {
            let source =
                resolve(state, temperature, &mechanical).map_err(CohesiveHeatError::Source)?;
            Ok(ThirdBodyHeatInput {
                dissipated_j: loss.viscous_j,
                conductivity_w_m_k: source.conductivity_w_m_k,
                area_resistance_m2_k_w: source.area_resistance_m2_k_w,
                transport: source.transport,
            })
        };
        if entropy {
            step.evaluate_with_entropy(&increments, endpoint, check)
        } else {
            step.evaluate(&increments, endpoint, check)
        }
    };
    let problem =
        FallibleFiniteDifference::new(cells, controls.difference_step_j, |x: &[f64]| {
            evaluate(x, false).map(|trial| trial.residual_j().to_vec())
        })
        .map_err(numerical)?;
    let construction = NewtonKrylovState::new(&problem, energy_guess_j.to_vec(), controls.newton);
    if let Some(error) = problem.take_failure() {
        return Err(numerical(error));
    }
    let mut solver = construction?;
    let result = solver.run_controlled(&problem, controls.maximum_newton_iterations, &mut || {
        gate.is_requested() || problem.has_failed()
    });
    if let Some(error) = problem.take_failure() {
        return Err(numerical(error));
    }
    check()?;
    if result.progress != SolverRunProgress::Complete {
        return Err(CohesiveHeatError::Interrupted(result.progress));
    }
    if !result.report.converged {
        return Err(CohesiveHeatError::NotConverged(result.report));
    }
    let (thermal, source, heat) =
        evaluate(&solver.x, controls.audit_entropy)?.accept(controls.heat, check)?;
    let combined = mechanical.kinetic_after_j - mechanical.kinetic_before_j
        + mechanical.bulk_stored_after_j
        - mechanical.bulk_stored_before_j
        + mechanical.cohesive_stored_after_j
        - mechanical.cohesive_stored_before_j
        + heat.thermal.stored_energy_change_j
        + mechanical.cohesive_energy_j.dissipated
        - mechanical.external_work_j
        - mechanical.prescribed_jump_work_j
        - heat.external_heat_j;
    let summed = heat.total_energy_residual_j - mechanical.energy_residual_j;
    if [combined, summed]
        .iter()
        .any(|v| !v.is_finite() || v.abs() > controls.combined_energy_tolerance_j)
    {
        return Err(CohesiveHeatError::Invalid(
            "native adhesive and finite layer energy do not balance",
        ));
    }
    check()?;
    Ok(HeatedCohesiveLayerCandidate {
        fracture_work_j: loss.fracture_j,
        mechanical,
        thermal,
        heat,
        source,
        nonlinear: result.report,
        thermal_trials: trials.get(),
        combined_energy_residual_j: combined,
        summed_energy_residual_j: summed,
    })
}

fn same_area(a: f64, b: f64) -> bool {
    a.is_finite()
        && b.is_finite()
        && a > 0.0
        && b > 0.0
        && (a - b).abs() <= 64.0 * f64::EPSILON * a.max(b)
}

fn numerical<E>(error: FiniteDifferenceError<CohesiveHeatError<E>>) -> CohesiveHeatError<E> {
    match error {
        FiniteDifferenceError::Source(e) => e,
        FiniteDifferenceError::Invalid(field) => CohesiveHeatError::Invalid(field),
    }
}
