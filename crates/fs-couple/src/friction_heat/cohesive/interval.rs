//! Complete deformable cohesive stopping paths and original finite-body heat.
//!
//! The original immutable interval owner retains inertia, bulk elasticity,
//! damage, viscosity, friction and prescribed-jump work. Only its accepted
//! per-patch losses feed the shared heat network; no trial is advanced twice.

use super::*;
use fs_solid::cohesive::dynamics::interval::{DynamicsInterval, IntervalConfig};

/// Both complete endpoints, returned only after the independent joint audit.
pub struct HeatedCohesiveIntervalCandidate<'s, 'a, 't, P> {
    pub mechanical: DynamicsInterval<'s, 'a>,
    pub thermal: DepositionSimulation<'t>,
    pub heat: FrictionHeatNetworkReport,
    pub source: P,
    pub nonlinear: NewtonReport,
    pub thermal_trials: usize,
    /// Original patch IDs in canonical ascending thermal-link order.
    pub patch_ids: Vec<u64>,
    pub fracture_work_j: f64,
    /// Delta(K + bulk + bond + H) + fracture - load/prescribed work - other heat.
    pub combined_energy_residual_j: f64,
    pub summed_energy_residual_j: f64,
}

/// Complete a native interval under constant reduced force and affine
/// prescribed jump offsets, including stops, static recapture and renewed slip.
/// Every original sourced patch has one matching reference thermal face; shared
/// finite receiving cells retain their original network and caloric laws.
///
/// Actual accepted per-patch friction and surviving-bond viscosity supply heat
/// once. Fracture remains a separate surface-creation account. Initial and all
/// accepted local normal paths must stay closed; intermediate opening cannot
/// be hidden by final reclosure. Source, domain, budget, nonlinear, energy and
/// cancellation refusal leave both input owners unchanged.
///
/// This retains the parent's fixed material-reference geometry and explicit
/// ordered surface/thermal-face correspondence. Mechanical coefficients stay
/// frozen. Transport uses interval-average physical source power and endpoint
/// backward Euler; refine complete intervals to resolve heat timing/feedback.
/// The caloric entropy audit does not infer cohesive free-energy derivatives,
/// surface entropy, or a closed thermomechanical model.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn propose_deformable_cohesive_interval_heat<'s, 'a, 't, P, E>(
    motion: &CohesiveMotion<'s, 'a>,
    constant_force_n: &[f64],
    endpoint_prescribed_jumps_m: &[[f64; 3]],
    duration_s: f64,
    interval_controls: IntervalConfig,
    contacts: &[CohesiveHeatContact<'_>],
    thermal: &DepositionSimulation<'t>,
    controls: CohesiveHeatControls,
    energy_guess_j: &[f64],
    resolve: impl Fn(&DepositionSimulation<'t>, &[u64], &[f64]) -> Result<CohesiveHeatSource<P>, E>,
    gate: &CancelGate,
) -> Result<HeatedCohesiveIntervalCandidate<'s, 'a, 't, P>, CohesiveHeatError<E>> {
    let check = || {
        if gate.is_requested() {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    };
    let (ordered, indices) = ordered_heat_contacts(
        motion,
        thermal,
        contacts,
        duration_s,
        controls,
        energy_guess_j,
        gate,
    )?;
    let target = motion.time_s() + duration_s;
    let routes: Vec<_> = ordered.iter().map(|c| c.route).collect();
    let network =
        FrictionHeatNetworkStep::new(thermal, target, &routes, controls.maximum_interfaces, check)?;
    let mechanical = motion.advance_interval_controlled(
        constant_force_n,
        endpoint_prescribed_jumps_m,
        duration_s,
        interval_controls,
        &mut || gate.is_requested(),
    )?;
    if mechanical.motion.time_s() != target
        || mechanical.patch_dissipation.len() != contacts.len()
        || mechanical.maximum_normal_jump_m.len() != contacts.len()
    {
        return Err(CohesiveHeatError::Invalid(
            "complete native heat interval required",
        ));
    }
    let mut losses = Vec::with_capacity(contacts.len());
    let mut physical_heat = 0.0;
    for (contact, &index) in ordered.iter().zip(&indices) {
        check()?;
        let loss = mechanical.patch_dissipation[index];
        let maximum_jump = mechanical.maximum_normal_jump_m[index];
        if loss.patch_id != contact.patch_id
            || !maximum_jump.is_finite()
            || maximum_jump > 0.0
            || [loss.fracture_j, loss.viscous_j, loss.friction_j]
                .iter()
                .any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(CohesiveHeatError::Invalid(
                "every accepted native patch path must stay closed with nonnegative physical losses",
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
            "native interval loss differs from per-patch heat",
        ));
    }
    let guesses: Vec<_> = ordered.iter().map(|c| c.temperature_guess_k).collect();
    let solved = solve_frozen_cohesive_heat(
        thermal,
        &network,
        &guesses,
        &losses,
        controls,
        energy_guess_j,
        resolve,
        gate,
    )?;
    let combined = mechanical.kinetic_after_j - mechanical.kinetic_before_j
        + mechanical.bulk_stored_after_j
        - mechanical.bulk_stored_before_j
        + mechanical.cohesive_stored_after_j
        - mechanical.cohesive_stored_before_j
        + solved.heat.thermal.stored_energy_change_j
        + mechanical.cohesive_energy_j.dissipated
        - mechanical.external_work_j
        - mechanical.prescribed_jump_work_j
        - solved.heat.external_heat_j;
    let summed = solved.heat.total_energy_residual_j - mechanical.energy_residual_j;
    if [combined, summed]
        .iter()
        .any(|v| !v.is_finite() || v.abs() > controls.combined_energy_tolerance_j)
    {
        return Err(CohesiveHeatError::Invalid(
            "native cohesive interval and finite heat do not balance",
        ));
    }
    check()?;
    Ok(HeatedCohesiveIntervalCandidate {
        fracture_work_j: mechanical.cohesive_energy_j.dissipated,
        mechanical,
        thermal: solved.thermal,
        heat: solved.heat,
        source: solved.source,
        nonlinear: solved.nonlinear,
        thermal_trials: solved.trials,
        patch_ids: ordered.iter().map(|c| c.patch_id).collect(),
        combined_energy_residual_j: combined,
        summed_energy_residual_j: summed,
    })
}
