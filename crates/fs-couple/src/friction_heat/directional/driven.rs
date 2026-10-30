//! Explicit mechanical drive, sourced directional friction and original heat.
//!
//! This is the driven consumer of the SAME material/geometry/thermal transaction
//! as the unforced trajectory. The original rigid owner accounts for applied
//! forces and torques, including common-body acceleration during static hold.
//! Only physical friction loss enters the original finite enthalpy cells;
//! external mechanical work remains a separate supply in the joint balance.

use super::trajectory::march::{
    DirectionalHeatMarchError, DirectionalHeatMarchReport, DirectionalMarchControls,
    advance_directional_heat_phases,
};
use super::trajectory::{
    DirectionalHeatMechanics, DirectionalTrajectoryHeatError, HeatedDirectionalTrajectoryReport,
    solve_directional_heat_phase,
};
use super::*;
use fs_contact::interface_binding::thermal_friction::elliptic::driven::SourceDrivenError;
use fs_mbd::Wrench;
use fs_mbd::tangential_friction::directional::driven::{
    RigidDirectionalDrivenError, RigidDirectionalDrivenReport, apply_rigid_directional_driven_with,
};
use fs_tribo::tangential::elliptic::driven::{EllipticDrivenControls, EllipticDrivenMode};

pub type HeatedDirectionalDrivenReport<P> =
    HeatedDirectionalTrajectoryReport<P, RigidDirectionalDrivenReport>;
pub type DirectionalDrivenHeatError<E> =
    DirectionalTrajectoryHeatError<E, RigidDirectionalDrivenError<SourceDrivenError>>;

pub type DrivenDirectionalHeatMarchReport<P> =
    DirectionalHeatMarchReport<P, RigidDirectionalDrivenReport>;
pub type DrivenDirectionalHeatMarchError<E> =
    DirectionalHeatMarchError<E, RigidDirectionalDrivenError<SourceDrivenError>>;

/// Refine driven material/thermal phases using the original shared controller.
///
/// The explicit external wrenches are held throughout this fixed-geometry span.
/// Accepted mode persists through private steps and publishes with original
/// bodies, source history and finite heat only after the whole span succeeds.
/// Work supplied by the drive is included in each phase and whole-span balance.
/// The controller retains original transport payloads and bounds accumulated
/// source-path error. A later failure rolls back earlier private accepted phases.
///
/// Outer refinement resolves changing source/thermal samples. The native curved
/// driven discretization has its own maximum_step_s and must also be refined;
/// an outer comparison can share a common inner discretization error. Local
/// estimates do not certify a global bound or a unique discontinuous-law branch.
#[allow(clippy::too_many_arguments)]
pub fn advance_source_directional_driven_heat<'m, P, E>(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    history: &mut SourceThermalEllipticState,
    mode: &mut EllipticDrivenMode,
    law: &SourceThermalEllipticCoulomb,
    contact: RigidDirectionalContact,
    external: [Wrench; 2],
    thermal: &mut DepositionSimulation<'m>,
    route: &FrictionHeatRoute,
    material: DirectionalHeatMaterial<'_>,
    duration_s: f64,
    phase_controls: DirectionalHeatSolve,
    driven: EllipticDrivenControls,
    controls: DirectionalMarchControls,
    interface_temperature_guess_k: f64,
    resolve: impl Fn(
        &DepositionSimulation<'m>,
        [EquilibriumPhaseState; 2],
    ) -> Result<DirectionalHeatSource<P>, E>,
    gate: &CancelGate,
) -> Result<DrivenDirectionalHeatMarchReport<P>, DrivenDirectionalHeatMarchError<E>> {
    let mut joint = (history.clone(), *mode);
    let result = advance_directional_heat_phases(
        bodies,
        properties,
        &mut joint,
        law,
        contact,
        thermal,
        duration_s,
        phase_controls,
        controls,
        interface_temperature_guess_k,
        |state| &state.0,
        |pair, properties, contact, state, thermal, dt, controls, guess, ti| {
            apply_source_directional_driven_heat(
                pair,
                properties,
                &mut state.0,
                &mut state.1,
                law,
                contact,
                external,
                thermal,
                route,
                material,
                dt,
                controls,
                driven,
                guess,
                ti,
                &resolve,
                gate,
            )
        },
        gate,
    )?;
    *history = joint.0;
    *mode = joint.1;
    Ok(result)
}

/// Advance one driven source/rigid/heat phase with accepted mode continuation.
///
/// `external` uses contact A/B order. Forces are world vectors and torques are
/// body vectors about each original centre of mass. Poses, arms, contact frame,
/// inertia and normal load are held over the phase; actual normal recoil is
/// reported by mechanics. Normal load and external wrenches remain distinct
/// physical inputs. Temperature and pressure fields on `contact` must be None.
///
/// Resting/Sliding is accepted mechanical state, not a speed cutoff. Its native
/// residual admission prevents roundoff after arrest from causing breakaway
/// between static and kinetic thresholds. The same original source claims can
/// change those thresholds as the solved finite surface temperatures evolve.
///
/// Original body momenta, source history, mode and thermal state accept together.
/// The energy gate includes explicit drive work; it never credits work as heat
/// merely because friction holds a relative contact stationary. Native curved
/// drive steps and outer source/thermal sampling have separate refinement axes.
#[allow(clippy::too_many_arguments)]
pub fn apply_source_directional_driven_heat<'m, P, E>(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    history: &mut SourceThermalEllipticState,
    mode: &mut EllipticDrivenMode,
    law: &SourceThermalEllipticCoulomb,
    contact: RigidDirectionalContact,
    external: [Wrench; 2],
    thermal: &mut DepositionSimulation<'m>,
    route: &FrictionHeatRoute,
    material: DirectionalHeatMaterial<'_>,
    dt_s: f64,
    controls: DirectionalHeatSolve,
    driven: EllipticDrivenControls,
    energy_guess_j: &[f64],
    interface_temperature_guess_k: f64,
    resolve: impl Fn(
        &DepositionSimulation<'m>,
        [EquilibriumPhaseState; 2],
    ) -> Result<DirectionalHeatSource<P>, E>,
    gate: &CancelGate,
) -> Result<HeatedDirectionalDrivenReport<P>, DirectionalDrivenHeatError<E>> {
    let mut joint = (history.clone(), *mode);
    let result = solve_directional_heat_phase(
        bodies,
        properties,
        &mut joint,
        law,
        contact,
        thermal,
        route,
        material,
        dt_s,
        controls,
        energy_guess_j,
        interface_temperature_guess_k,
        |pair, properties, contact, accepted, dt, query| {
            let mut motion = None;
            let report = apply_rigid_directional_driven_with(
                pair,
                properties,
                contact,
                external,
                &mut accepted.0,
                &mut accepted.1,
                dt,
                controls.mechanical,
                |state, input| {
                    let source = law.advance_driven_motion_consistent(
                        state,
                        input,
                        query,
                        controls.motion,
                        driven,
                        &mut || gate.is_requested(),
                    )?;
                    motion = Some((
                        [
                            source.queried_distance_m,
                            source.actual_distance_m,
                            source.distance_residual_m,
                        ],
                        source.evaluations,
                    ));
                    Ok::<_, SourceDrivenError>((source.driven, source.state))
                },
                gate,
            )?;
            let (motion, source_evaluations) = motion.ok_or_else(|| {
                RigidDirectionalDrivenError::Driven(SourceDrivenError::Query(
                    MotionQueryError::Invalid("missing driven source-motion receipt"),
                ))
            })?;
            Ok(DirectionalHeatMechanics {
                kinetic_energy_change_j: report.kinetic_energy_change_j,
                external_work_j: report.external_work_j,
                dissipated_j: report.driven.dissipated_j,
                net_slip_m: report.driven.slip_increment_m,
                report,
                motion,
                source_evaluations,
            })
        },
        resolve,
        gate,
    )?;
    // No fallible work or callback follows the common publication boundary.
    *history = joint.0;
    *mode = joint.1;
    Ok(result)
}
