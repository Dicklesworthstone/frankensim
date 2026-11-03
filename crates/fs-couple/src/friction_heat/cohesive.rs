//! Deformable adhesive reclosure and original finite-body caloric transport.
//!
//! One native mechanical proposal supplies the actual per-patch friction and
//! surviving-bond viscous losses. The existing massless-face heat network solves
//! their temperatures and signed exchange with finite bodies. Intrinsic fracture
//! work remains outside the caloric account: creating damaged surfaces is not
//! automatically heat. No prescribed partition or added heat capacity is used.
//!
//! This is a fixed-reference, closed-contact, frozen-mechanical-coefficient split.
//! Mechanical source-temperature transfer remains explicit between increments.
//! The caller declares each mechanical patch's coincidence and ordered A/B
//! material assignment with its thermal face; a reduced jump map alone cannot
//! reconstruct world placement. Opening, moving contact coverage, finite bond
//! thermal resistance/storage and strong thermal feedback need their own owners.

use super::network::{
    FrictionHeatNetworkInput, FrictionHeatNetworkReport, FrictionHeatNetworkStep,
    InterfaceHeatInput,
};
use super::{FrictionHeatRoute, FrictionHeatTolerances};
use fs_exec::CancelGate;
use fs_solid::cohesive::dynamics::{CohesiveMotion, DynamicsError, DynamicsStep};
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress,
};
use fs_thermal::enthalpy::backward_euler::resolvable_increment_j;
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};
use std::cell::Cell;
use std::collections::{BTreeMap, BTreeSet};

/// Same canonical link-keyed conductivity/other-source inputs as rigid heat.
pub use super::hertz::multiple::{
    RigidHertzIslandHeatConductivity as CohesiveHeatConductivity,
    RigidHertzIslandHeatSource as CohesiveHeatSource,
};

/// Explicit correspondence between one original cohesive patch and heat face.
/// The route's A/B cell order must be the original source's surface A/B order.
#[derive(Clone, Copy)]
pub struct CohesiveHeatContact<'a> {
    pub patch_id: u64,
    pub route: &'a FrictionHeatRoute,
    /// Positive numerical seed; the common face temperature is solved.
    pub temperature_guess_k: f64,
}

/// Thermal solve bounds; mechanics retains its original independent controls.
#[derive(Clone, Copy, Debug)]
pub struct CohesiveHeatControls {
    pub maximum_cells: usize,
    pub maximum_interfaces: usize,
    pub maximum_unknowns: usize,
    pub maximum_trials: usize,
    pub maximum_newton_iterations: usize,
    pub difference_step_j: f64,
    /// Numerical J/K coordinate scaling, never a physical interface capacity.
    pub temperature_scale_j_per_k: f64,
    pub heat: FrictionHeatTolerances,
    pub combined_energy_tolerance_j: f64,
    pub newton: NewtonKrylovConfig,
    /// Audit the original caloric entropy on the accepted root only.
    pub audit_entropy: bool,
}

#[derive(Debug)]
pub enum CohesiveHeatError<E> {
    Invalid(&'static str),
    Mechanical(DynamicsError),
    Thermal(DepositionError),
    Source(E),
    Newton(NewtonError),
    NotConverged(NewtonReport),
    Interrupted(SolverRunProgress),
    TrialBudget,
}
impl<E> From<DynamicsError> for CohesiveHeatError<E> {
    fn from(e: DynamicsError) -> Self {
        Self::Mechanical(e)
    }
}
impl<E> From<DepositionError> for CohesiveHeatError<E> {
    fn from(e: DepositionError) -> Self {
        Self::Thermal(e)
    }
}
impl<E> From<NewtonError> for CohesiveHeatError<E> {
    fn from(e: NewtonError) -> Self {
        Self::Newton(e)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for CohesiveHeatError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "deformable cohesive heat refused: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for CohesiveHeatError<E> {}

/// The native mechanical and caloric endpoints publish together, exactly once.
pub struct HeatedCohesiveCandidate<'s, 'a, 't, P> {
    pub mechanical: DynamicsStep<'s, 'a>,
    pub thermal: DepositionSimulation<'t>,
    pub heat: FrictionHeatNetworkReport,
    pub source: P,
    pub nonlinear: NewtonReport,
    pub thermal_trials: usize,
    /// In the same ascending link order as heat.interfaces.
    pub patch_ids: Vec<u64>,
    /// Physical fracture work remains in its original nonthermal account [J].
    pub fracture_work_j: f64,
    /// Delta(K + bulk/bond storage + H) + fracture - applied work - other heat.
    pub combined_energy_residual_j: f64,
    /// Independent thermal residual minus the native mechanical input residual.
    pub summed_energy_residual_j: f64,
}

/// Advance the actual deformable joint and all its finite thermal receivers.
///
/// Every original patch requires one distinct, matching-area heat route and
/// source binding. Shared receiving cells are supported. The complete straight
/// normal-jump path must remain closed; otherwise the entire proposal refuses.
/// Per-patch losses come from the native accepted mechanical replay, never a
/// caller-authored receipt or a second material update. Fracture work is retained
/// separately while friction and physical bond viscosity enter heat once.
///
/// Mechanics is solved once at its original source coefficients; only thermal
/// transport is iterated. The callback supplies endpoint conductivities and
/// OTHER heat, and must stage its payload until acceptance. Update source
/// temperatures explicitly before a subsequent mechanical phase and refine the
/// split timestep for feedback accuracy. All native caloric domains, events,
/// conduction links and external reservoirs remain with their original owner.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn propose_deformable_cohesive_heat<'s, 'a, 't, P, E>(
    motion: &CohesiveMotion<'s, 'a>,
    mean_force_n: &[f64],
    endpoint_prescribed_jumps_m: &[[f64; 3]],
    dt_s: f64,
    contacts: &[CohesiveHeatContact<'_>],
    thermal: &DepositionSimulation<'t>,
    controls: CohesiveHeatControls,
    energy_guess_j: &[f64],
    resolve: impl Fn(&DepositionSimulation<'t>, &[u64], &[f64]) -> Result<CohesiveHeatSource<P>, E>,
    gate: &CancelGate,
) -> Result<HeatedCohesiveCandidate<'s, 'a, 't, P>, CohesiveHeatError<E>> {
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
    let count = contacts.len();
    let dimension = cells
        .checked_add(count)
        .ok_or(CohesiveHeatError::Invalid("thermal dimension overflow"))?;
    let target = motion.time_s() + dt_s;
    let patches = motion.state().system().patches();
    if motion.time_s() != thermal.time_s()
        || !dt_s.is_finite()
        || dt_s <= 0.0
        || !target.is_finite()
        || target <= motion.time_s()
        || cells == 0
        || cells > controls.maximum_cells
        || count == 0
        || count != patches.len()
        || count > controls.maximum_interfaces
        || dimension > controls.maximum_unknowns
        || controls.maximum_trials == 0
        || controls.maximum_newton_iterations == 0
        || energy_guess_j.len() != cells
        || energy_guess_j.iter().any(|v| !v.is_finite())
        || !controls.difference_step_j.is_finite()
        || controls.difference_step_j <= 0.0
        || !controls.temperature_scale_j_per_k.is_finite()
        || controls.temperature_scale_j_per_k <= 0.0
        || [
            controls.heat.cell_j,
            controls.heat.thermal_j,
            controls.heat.interface_j,
            controls.heat.total_j,
            controls.combined_energy_tolerance_j,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v < 0.0)
    {
        return Err(CohesiveHeatError::Invalid(
            "bounded complete synchronized mechanical/thermal state required",
        ));
    }
    let indices: BTreeMap<_, _> = patches
        .iter()
        .enumerate()
        .map(|(i, p)| (p.id(), i))
        .collect();
    let mut ordered: Vec<_> = contacts.iter().collect();
    ordered.sort_by_key(|c| c.route.link_id());
    let mut seen = BTreeSet::new();
    let mut mechanical_indices = Vec::with_capacity(count);
    for contact in &ordered {
        check()?;
        let &index = indices
            .get(&contact.patch_id)
            .ok_or(CohesiveHeatError::Invalid("unknown cohesive heat patch"))?;
        let patch = &patches[index];
        let area = patch.area_m2();
        let other = contact.route.contact_area_m2();
        if !seen.insert(contact.patch_id)
            || patch.source().is_none()
            || (area - other).abs() > 64.0 * f64::EPSILON * area.max(other)
            || !contact.temperature_guess_k.is_finite()
            || contact.temperature_guess_k <= 0.0
            || motion.state().points()[index].jump_m()[0] > 0.0
        {
            return Err(CohesiveHeatError::Invalid(
                "heat needs unique sourced closed patches with matching original face area",
            ));
        }
        mechanical_indices.push(index);
    }
    let routes: Vec<_> = ordered.iter().map(|c| c.route).collect();
    let step =
        FrictionHeatNetworkStep::new(thermal, target, &routes, controls.maximum_interfaces, check)?;
    let mechanical = motion.advance_controlled(
        mean_force_n,
        endpoint_prescribed_jumps_m,
        step.duration_s(),
        &mut || gate.is_requested(),
    )?;
    if mechanical.motion.time_s() != target || mechanical.patch_dissipation.len() != count {
        return Err(CohesiveHeatError::Invalid(
            "native cohesive endpoint changed clock or patch coverage",
        ));
    }
    let mut losses = Vec::with_capacity(count);
    let mut physical_heat = 0.0;
    for (contact, &index) in ordered.iter().zip(&mechanical_indices) {
        check()?;
        let loss = mechanical.patch_dissipation[index];
        if loss.patch_id != contact.patch_id
            || mechanical.motion.state().points()[index].jump_m()[0] > 0.0
            || [loss.fracture_j, loss.viscous_j, loss.friction_j]
                .iter()
                .any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(CohesiveHeatError::Invalid(
                "native closed-contact loss or original patch identity refused",
            ));
        }
        let heat = loss.viscous_j + loss.friction_j;
        physical_heat += heat;
        losses.push(heat);
    }
    let native_heat = mechanical.viscous_dissipated_j + mechanical.friction_dissipated_j;
    if !physical_heat.is_finite()
        || !native_heat.is_finite()
        || (physical_heat - native_heat).abs() > controls.combined_energy_tolerance_j
    {
        return Err(CohesiveHeatError::Invalid(
            "per-patch heat does not match original mechanical loss",
        ));
    }
    let trials = Cell::new(0usize);
    let evaluate = |x: &[f64], entropy: bool| {
        check()?;
        if trials.get() >= controls.maximum_trials {
            return Err(CohesiveHeatError::TrialBudget);
        }
        trials.set(trials.get() + 1);
        if x.len() != dimension || x.iter().any(|v| !v.is_finite()) {
            return Err(CohesiveHeatError::Invalid(
                "complete finite thermal iterate required",
            ));
        }
        let temperatures: Vec<_> = ordered
            .iter()
            .enumerate()
            .map(|(i, c)| c.temperature_guess_k + x[cells + i] / controls.temperature_scale_j_per_k)
            .collect();
        // Project only sub-ULP iterate tails. The original residual still
        // charges every source and flux, and its acceptance gates are unchanged.
        let increments: Vec<_> = thermal
            .cells()
            .values()
            .zip(&x[..cells])
            .map(|(cell, &q)| resolvable_increment_j(cell.specific_h_j_kg, q, cell.mass_kg))
            .collect();
        let endpoint = |state: &DepositionSimulation<'t>, temperatures: &[f64]| {
            let source =
                resolve(state, step.link_ids(), temperatures).map_err(CohesiveHeatError::Source)?;
            if source.interfaces.len() != count
                || source
                    .interfaces
                    .iter()
                    .zip(step.link_ids())
                    .any(|(input, id)| input.link_id != *id)
            {
                return Err(CohesiveHeatError::Invalid(
                    "conductivity must retain every original link in canonical order",
                ));
            }
            let mut interfaces = Vec::with_capacity(count);
            for (input, &loss) in source.interfaces.iter().zip(&losses) {
                check()?;
                interfaces.push(InterfaceHeatInput {
                    dissipated_j: loss,
                    conductivity_w_m_k: input.conductivity_w_m_k,
                });
            }
            Ok(FrictionHeatNetworkInput {
                interfaces,
                transport: source.transport,
            })
        };
        if entropy {
            step.evaluate_with_entropy(&increments, &temperatures, endpoint, check)
        } else {
            step.evaluate(&increments, &temperatures, endpoint, check)
        }
    };
    let problem =
        FallibleFiniteDifference::new(dimension, controls.difference_step_j, |x: &[f64]| {
            evaluate(x, false).map(|trial| trial.residual_j().to_vec())
        })
        .map_err(numerical)?;
    let mut initial = energy_guess_j.to_vec();
    initial.resize(dimension, 0.0);
    let construction = NewtonKrylovState::new(&problem, initial, controls.newton);
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
    let trial = evaluate(&solver.x, controls.audit_entropy)?;
    let heat = trial.report();
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
            "original cohesive mechanics and finite heat do not balance",
        ));
    }
    let (thermal, source, heat) = trial.accept(controls.heat, check)?;
    check()?;
    Ok(HeatedCohesiveCandidate {
        fracture_work_j: mechanical.cohesive_energy_j.dissipated,
        mechanical,
        thermal,
        heat,
        source,
        nonlinear: result.report,
        thermal_trials: trials.get(),
        patch_ids: ordered.iter().map(|c| c.patch_id).collect(),
        combined_energy_residual_j: combined,
        summed_energy_residual_j: summed,
    })
}

fn numerical<E>(error: FiniteDifferenceError<CohesiveHeatError<E>>) -> CohesiveHeatError<E> {
    match error {
        FiniteDifferenceError::Source(e) => e,
        FiniteDifferenceError::Invalid(field) => CohesiveHeatError::Invalid(field),
    }
}
