//! Original rigid stop/restick paths and finite-body heat accept together.
//!
//! Accepted native substeps own all mechanical history and physical losses.
//! Their per-site sums supply one finite enthalpy solve over the requested
//! interval. Rejected stopping trials supply no heat. The final sticking
//! substep cannot erase heat already generated while the bodies were sliding.
//! Mechanical coefficients remain frozen; the first-order thermal split needs
//! interval refinement and is not a closed temperature-feedback calculation.

use super::*;
use crate::cohesive_contact::rigid::{RigidCohesiveInterval, RigidCohesiveIntervalConfig};

/// Full interval mechanics and the corresponding original caloric endpoint.
pub struct HeatedRigidCohesiveIntervalReport<P> {
    pub mechanical: RigidCohesiveInterval,
    pub heat: FrictionHeatNetworkReport,
    pub source: P,
    pub nonlinear: NewtonReport,
    pub thermal_trials: usize,
    /// Native site indices in ascending original thermal-link order.
    pub site_indices: Vec<usize>,
    pub fracture_work_j: f64,
    /// Delta(native kinetic/gravity/bond + H) + fracture - work - other heat.
    pub combined_energy_residual_j: f64,
    pub summed_energy_residual_j: f64,
}

/// Advance original rigid bodies through stops, resticking and finite heating.
///
/// Loads stay constant and every source/history remains owned by the native
/// interval stepper. All original sites require their same unique, matching-area
/// thermal routes and explicit ordered body/material assignments. Every accepted
/// straight local constitutive path must stay closed, including intermediate
/// substeps; a closed final endpoint cannot hide an earlier open contact.
/// Geometry remains the parent's fixed material-reference thermal attachment.
///
/// Exact accepted per-site friction plus surviving-bond viscosity is deposited
/// once. Fracture retains its separate surface account, with no numerical or
/// parameter-work heating. Transport uses the interval-average source power and
/// endpoint backward Euler; shorter complete intervals refine heat timing.
/// The optional entropy account covers original caloric transport and physical
/// heat, without inferring cohesive free-energy derivatives or surface entropy.
///
/// Thermal/source/energy refusal or cancellation restores the entire incoming
/// native state, observations and clock, even after provisional stop/restick
/// substeps. The caloric endpoint publishes only after the native final gate.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn apply_rigid_cohesive_interval_heat<'a, 't, P, E>(
    mechanical: &mut RigidCohesiveSystem<'a>,
    loads: [RigidCohesiveLoad; 2],
    duration_s: f64,
    interval_controls: RigidCohesiveIntervalConfig,
    thermal: &mut DepositionSimulation<'t>,
    bodies: [RigidCohesiveHeatBody<'_>; 2],
    contacts: &[RigidCohesiveHeatContact<'_>],
    controls: CohesiveHeatControls,
    energy_guess_j: &[f64],
    resolve: impl Fn(&DepositionSimulation<'t>, &[u64], &[f64]) -> Result<CohesiveHeatSource<P>, E>,
    gate: &CancelGate,
) -> Result<HeatedRigidCohesiveIntervalReport<P>, RigidCohesiveHeatError<E>> {
    let check = || {
        if gate.is_requested() {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    };
    check()?;
    let ordered = ordered_heat_contacts(
        mechanical,
        thermal,
        bodies,
        contacts,
        controls.maximum_interfaces,
        gate,
    )?;
    let start = mechanical.time_s();
    let target = start + duration_s;
    if !duration_s.is_finite() || duration_s <= 0.0 || !target.is_finite() || target <= start {
        return Err(RigidCohesiveHeatError::Invalid(
            "finite advancing heat interval required",
        ));
    }
    let routes: Vec<_> = ordered.iter().map(|c| c.route).collect();
    let guesses: Vec<_> = ordered.iter().map(|c| c.temperature_guess_k).collect();
    let network =
        FrictionHeatNetworkStep::new(thermal, target, &routes, controls.maximum_interfaces, check)?;
    let (interval, (solved, combined, summed)) = mechanical.advance_interval_with_acceptance(
        loads, duration_s, interval_controls, Some(gate), |interval| {
            if interval.initial_time_s != start || interval.time_s != target
                || interval.sites.len() != ordered.len()
            {
                return Err(RigidCohesiveHeatError::Invalid("complete native heat interval required"));
            }
            let mut losses = Vec::with_capacity(ordered.len());
            let mut total_heat = 0.0;
            for contact in &ordered {
                check()?;
                let site = &interval.sites[contact.site_index];
                if !site.maximum_normal_jump_m.is_finite() || site.maximum_normal_jump_m > 0.0
                    || [site.friction_dissipation_j, site.viscous_dissipation_j,
                        site.fracture_dissipation_j].iter().any(|v| !v.is_finite() || *v < 0.0)
                {
                    return Err(RigidCohesiveHeatError::Invalid(
                        "every accepted native material path must remain closed with physical nonnegative losses",
                    ));
                }
                let loss = site.friction_dissipation_j + site.viscous_dissipation_j;
                total_heat += loss;
                losses.push(loss);
            }
            if !total_heat.is_finite()
                || (total_heat - interval.friction_dissipation_j - interval.viscous_dissipation_j).abs()
                    > controls.combined_energy_tolerance_j
            {
                return Err(RigidCohesiveHeatError::Invalid(
                    "accepted site heat differs from original complete interval loss",
                ));
            }
            let solved = solve_frozen_cohesive_heat(
                thermal, &network, &guesses, &losses, controls, energy_guess_j, resolve, gate,
            )?;
            let combined = interval.energy_after_j - interval.energy_before_j
                + solved.heat.thermal.stored_energy_change_j + interval.fracture_dissipation_j
                - interval.external_work_j - solved.heat.external_heat_j;
            let summed = interval.energy_residual_j + solved.heat.total_energy_residual_j;
            if [combined, summed].iter().any(|r| !r.is_finite() || r.abs() > controls.combined_energy_tolerance_j) {
                return Err(RigidCohesiveHeatError::Invalid(
                    "original complete native interval and finite heat do not balance",
                ));
            }
            Ok((solved, combined, summed))
        },
    )?;
    // The native rollback owner has accepted heat and passed its final gate.
    // No fallible operation separates the corresponding caloric publication.
    *thermal = solved.thermal;
    Ok(HeatedRigidCohesiveIntervalReport {
        fracture_work_j: interval.fracture_dissipation_j,
        mechanical: interval,
        heat: solved.heat,
        source: solved.source,
        nonlinear: solved.nonlinear,
        thermal_trials: solved.trials,
        site_indices: ordered.iter().map(|c| c.site_index).collect(),
        combined_energy_residual_j: combined,
        summed_energy_residual_j: summed,
    })
}
