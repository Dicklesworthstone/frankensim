//! One native sphere trajectory with separate finite-patch event and force laws.
//!
//! The original event engine detects impact, applies the original Hertz impulse,
//! and continues with the original Hertz sustained-support driver. Their radii,
//! coefficients, dimensional controls and source histories remain independent.
//! Neither an impact duration nor a pressure impulse is inferred from support.
//! Only actual native motion crosses the transition; no native momentum is reset
//! or external split clock inserted to join the two physical owners.

use super::hertz_impacts::{
    self, SpherePlaneHertzImpactError, SpherePlaneHertzImpactEvent, SpherePlaneHertzImpactPhase,
};
use super::hertz_support::{
    self, SpherePlaneHertzState, SpherePlaneHertzSupportAdvance, SpherePlaneHertzSupportControls,
    SpherePlaneHertzSupportError, SpherePlaneHertzSupportPhase,
};
use super::*;
use core::convert::Infallible;
use fs_tribo::sliding_patch::contact::HertzWrenchContact;
use fs_tribo::sliding_patch::contact::impulse::{HertzImpulseControls, HertzImpulsePatch};

/// Original numerical and physical budgets, without any inferred conversion
/// between instantaneous pressure impulse and a sustained normal force.
#[derive(Clone, Copy, Debug)]
pub struct SpherePlaneHertzContactControls {
    pub support: SpherePlaneHertzSupportControls,
    pub impulse: HertzImpulseControls,
}

/// One accepted original trajectory, plus its actual native impact responses.
/// `support.motion.dissipated_energy` already includes both impact and support
/// loss. The event receipts and support subaccounts are not additional heat.
#[derive(Clone, Debug)]
pub struct SpherePlaneHertzContactAdvance {
    pub support: SpherePlaneHertzSupportAdvance,
    pub impacts: Vec<SpherePlaneHertzImpactEvent>,
}

/// Either original owner may refuse. Both accepted source histories and the
/// incoming motion remain unchanged, including after an earlier private event.
#[derive(Debug)]
pub enum SpherePlaneHertzContactError<SE = Infallible, IE = Infallible> {
    Support(SpherePlaneHertzSupportError<SE>),
    Impact(SpherePlaneHertzImpactError<IE>),
}
impl<SE: core::fmt::Debug, IE: core::fmt::Debug> core::fmt::Display
    for SpherePlaneHertzContactError<SE, IE>
{
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite-patch sphere impact and support: {self:?}")
    }
}
impl<SE: core::fmt::Debug, IE: core::fmt::Debug> std::error::Error
    for SpherePlaneHertzContactError<SE, IE>
{
}

/// Advance with separately declared, held impact and sustained-support laws.
///
/// The same normal restitution, event times, flight, force/motor/plane work,
/// angular balance and attitude evolution remain with the original engine.
/// The actual post-impact shear/spin reaches support directly. Support's exact
/// pure-motion events and first-order mixed-motion refinement are unchanged.
#[allow(clippy::too_many_arguments)]
pub fn advance_sphere_plane_hertz_contact(
    scene: &SpherePlane,
    state: SpherePlaneHertzState,
    torque_world_nm: Vec3,
    support_law: &HertzWrenchContact,
    impact_patch: &HertzImpulsePatch,
    duration_s: f64,
    controls: SpherePlaneHertzContactControls,
    cancelled: impl FnMut() -> bool,
) -> Result<SpherePlaneHertzContactAdvance, SpherePlaneHertzContactError> {
    advance_sphere_plane_hertz_contact_with(
        scene,
        state,
        torque_world_nm,
        &mut (),
        &mut (),
        duration_s,
        controls,
        |_, _, _| Ok((*support_law, ())),
        |_, _, _| Ok((*impact_patch, ())),
        cancelled,
    )
}

/// Resolve support at its original physical phase origins and impacts at the
/// actual detected events. Each callback sees its own accepted candidate history
/// and returns its original constitutive owner, not a prescribed response.
///
/// Independent source types, error types and state progress are preserved. The
/// histories are cloned once for the whole proposal. Both publish only after
/// the complete native trajectory and a final cancellation check pass; swaps
/// finish before either displaced caller history can run a destructor. The
/// callbacks' own external effects retain their original ownership boundary.
/// This composition adds no law, impact duration, pressure model, heat split,
/// or additional momentum application.
#[allow(clippy::too_many_arguments)]
pub fn advance_sphere_plane_hertz_contact_with<S: Clone, I: Clone, SE, IE>(
    scene: &SpherePlane,
    state: SpherePlaneHertzState,
    torque_world_nm: Vec3,
    support_history: &mut S,
    impact_history: &mut I,
    duration_s: f64,
    controls: SpherePlaneHertzContactControls,
    resolve_support: impl FnMut(
        &S,
        &SpherePlaneHertzSupportPhase,
        &mut dyn FnMut() -> bool,
    ) -> Result<(HertzWrenchContact, S), SE>,
    mut resolve_impact: impl FnMut(
        &I,
        &SpherePlaneHertzImpactPhase,
        &mut dyn FnMut() -> bool,
    ) -> Result<(HertzImpulsePatch, I), IE>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<SpherePlaneHertzContactAdvance, SpherePlaneHertzContactError<SE, IE>> {
    let mechanical =
        |e| SpherePlaneHertzContactError::Support(SpherePlaneHertzSupportError::Mechanics(e));
    poll(&mut cancelled).map_err(mechanical)?;
    controls
        .impulse
        .validate()
        .map_err(|e| SpherePlaneHertzContactError::Impact(SpherePlaneHertzImpactError::Law(e)))?;
    hertz_support::validate_controls(controls.support)
        .map_err(SpherePlaneHertzContactError::Support)?;
    let mut impact_candidate = impact_history.clone();
    let mut impact_error = None;
    let mut impacts = Vec::new();
    let result = hertz_support::candidate_with_impact(
        scene,
        state,
        torque_world_nm,
        support_history.clone(),
        duration_s,
        controls.support,
        resolve_support,
        |engine, body, time, normal_before, normal, slip, cancel| {
            let phase = match hertz_impacts::phase(engine, body, time, normal_before, normal, slip)
            {
                Ok(phase) => phase,
                Err(error) => {
                    impact_error = Some(SpherePlaneHertzImpactError::Mechanics(error));
                    return Err(SpherePlaneError::Unresolved("Hertz impact phase refused"));
                }
            };
            let (patch, next) = match resolve_impact(&impact_candidate, &phase, &mut || cancel()) {
                Ok(value) => value,
                Err(error) => {
                    impact_error = Some(SpherePlaneHertzImpactError::Source(error));
                    return Err(SpherePlaneError::Unresolved("Hertz impact source refused"));
                }
            };
            let (accepted, event) =
                match hertz_impacts::solve_event(engine, phase, patch, controls.impulse, cancel) {
                    Ok(value) => value,
                    Err(SpherePlaneHertzImpactError::Mechanics(error)) => {
                        impact_error = Some(SpherePlaneHertzImpactError::Mechanics(error));
                        return Err(SpherePlaneError::Unresolved(
                            "Hertz impact mechanics refused",
                        ));
                    }
                    Err(SpherePlaneHertzImpactError::Law(error)) => {
                        impact_error = Some(SpherePlaneHertzImpactError::Law(error));
                        return Err(SpherePlaneError::Unresolved("Hertz impact law refused"));
                    }
                    Err(SpherePlaneHertzImpactError::Source(never)) => match never {},
                };
            impacts.push(event);
            impact_candidate = next;
            Ok(Some(accepted))
        },
        &mut cancelled,
    );
    if let Some(error) = impact_error {
        return Err(SpherePlaneHertzContactError::Impact(error));
    }
    let (support, mut support_candidate) = result.map_err(SpherePlaneHertzContactError::Support)?;
    poll(&mut cancelled).map_err(mechanical)?;
    core::mem::swap(support_history, &mut support_candidate);
    core::mem::swap(impact_history, &mut impact_candidate);
    Ok(SpherePlaneHertzContactAdvance { support, impacts })
}
