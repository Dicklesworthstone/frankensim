//! Sourced directional friction through moving-sphere impacts and support.
//!
//! Original material claims feed the original native impulse and driven laws.
//! The sphere engine owns event times, pose, spin, normal response and work.
//! Conditions come from explicit physical owners at each real event/phase;
//! this bridge does not infer contact pressure from a point or allocate impact
//! loss to heat. Its impact receipts survive subsequent support-source updates.

use fs_contact::interface_binding::thermal_friction::elliptic::driven::SourceDrivenSegmentError;
use fs_contact::interface_binding::thermal_friction::elliptic::impulse::{
    SourceImpulseCoefficientRole, SourceImpulseError, SourceImpulseReport,
};
use fs_contact::interface_binding::thermal_friction::elliptic::motion::MotionQuerySolve;
use fs_contact::interface_binding::thermal_friction::elliptic::{
    SourceFrictionConditions, SourceThermalEllipticCoulomb, SourceThermalEllipticState,
};
use fs_mbd::sphere_plane::directional::impact::{
    SpherePlaneDirectionalImpactPhase, advance_sphere_plane_directional_with_impacts,
};
use fs_mbd::sphere_plane::directional::{
    SpherePlaneDirectionalAdvance, SpherePlaneDirectionalControls, SpherePlaneDirectionalError,
    SpherePlaneDirectionalPhase,
};
use fs_mbd::sphere_plane::{SpherePlane, SpherePlaneError, SpherePlaneState};
use fs_tribo::partial_slip::TangentFrame;
use fs_tribo::tangential::elliptic::impulse::EllipticImpulseControls;

/// Conditions and source-distance bracket for this actual support phase.
/// A speed-limited source can use [0, maximum_speed * phase.input.dt_s];
/// its accepted path and all bracket endpoints must remain in the source domain.
#[derive(Clone, Debug)]
pub struct SourceSpherePlaneSupportConditions {
    pub friction: SourceFrictionConditions,
    pub motion: MotionQuerySolve,
}

/// Independent motion, instantaneous-law and numerical acceptance controls.
#[derive(Clone, Copy, Debug)]
pub struct SourceSpherePlaneControls {
    pub motion: SpherePlaneDirectionalControls,
    pub impulse: EllipticImpulseControls,
    /// One explicitly chosen source coefficient ellipse for each whole impact.
    /// This is never an automatic static/kinetic switch during collision.
    pub impact_coefficient_role: SourceImpulseCoefficientRole,
}

/// Original material claims and mechanics at one accepted collision.
/// The response's dissipated work is tangential impact loss, not friction heat.
/// The original trajectory separately reports the complete normal+tangent loss.
#[derive(Clone, Debug)]
pub struct SourceSpherePlaneImpactReceipt {
    pub phase: SpherePlaneDirectionalImpactPhase,
    pub response: SourceImpulseReport,
}

/// All receipts belong to the same accepted physical trajectory.
#[derive(Clone, Debug)]
pub struct SourceSpherePlaneAdvance {
    pub trajectory: SpherePlaneDirectionalAdvance,
    /// Chronological, bounded by the original motion.max_events control.
    pub impacts: Vec<SourceSpherePlaneImpactReceipt>,
    /// Original source/native support trials, excluding accepted replay.
    /// Each such trial resolves both principal source bundles.
    pub support_source_evaluations: usize,
}

#[derive(Debug)]
pub enum SourceSpherePlaneError<E> {
    SupportConditions(E),
    ImpactConditions(E),
    Support(SourceDrivenSegmentError<SpherePlaneError>),
    Impact(SourceImpulseError),
    Invalid(&'static str),
}
impl<E: core::fmt::Debug> core::fmt::Display for SourceSpherePlaneError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source sphere/plane friction: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for SourceSpherePlaneError<E> {}

/// Advance the original sphere with one ordered source history across impacts
/// and sustained friction. The scalar scene friction is superseded at every
/// impact and support phase; original normal restitution remains authoritative.
///
/// Query callbacks receive actual pre-phase mechanics and the latest private
/// accepted history. They supply physical temperature/pressure and a complete
/// query without replacing velocity, mobility, normal reaction or time.
/// Support closes speed/path against actual native travel; impacts sample
/// preimpact speed and preserve accepted path exactly. Full static/kinetic
/// source coefficients survive impact for subsequent supported motion.
///
/// Callbacks must not publish external state: a later query, mechanics or
/// cancellation refusal discards the complete candidate and all receipts.
/// Body state is returned by value, and history is replaced only on success.
/// Ballistic flight alone neither queries friction nor advances its history.
///
/// The frame is fixed to the prescribed plane. Source phase refinement and
/// native sliding/attitude refinement remain separate controls. No finite-area
/// thermal receiver or impact-to-heat partition is supplied by this point model.
#[allow(clippy::too_many_arguments)]
pub fn advance_source_sphere_plane<E>(
    scene: &SpherePlane,
    state: SpherePlaneState,
    law: &SourceThermalEllipticCoulomb,
    history: &mut SourceThermalEllipticState,
    frame: TangentFrame,
    duration_s: f64,
    controls: SourceSpherePlaneControls,
    mut support_conditions: impl FnMut(
        &SourceThermalEllipticState,
        &SpherePlaneDirectionalPhase,
    ) -> Result<SourceSpherePlaneSupportConditions, E>,
    mut impact_conditions: impl FnMut(
        &SourceThermalEllipticState,
        &SpherePlaneDirectionalImpactPhase,
    ) -> Result<SourceFrictionConditions, E>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<SourceSpherePlaneAdvance, SpherePlaneDirectionalError<SourceSpherePlaneError<E>>> {
    let mut candidate = history.clone();
    let mut impacts = Vec::new();
    let mut support_source_evaluations = 0usize;
    let trajectory = advance_sphere_plane_directional_with_impacts(
        scene,
        state,
        &mut candidate,
        frame,
        duration_s,
        controls.motion,
        |accepted, phase, visit, cancel| {
            let conditions = support_conditions(accepted, phase)
                .map_err(SourceSpherePlaneError::SupportConditions)?;
            let mut input = phase.input;
            input.contact_temperature_k = conditions.friction.contact_temperature_k;
            input.nominal_pressure_pa = conditions.friction.nominal_pressure_pa;
            let report = law
                .advance_driven_motion_consistent_with_segments(
                    accepted,
                    &input,
                    &conditions.friction.query,
                    conditions.motion,
                    controls.motion.driven,
                    &mut || cancel(),
                    &mut |segment| visit(segment),
                )
                .map_err(SourceSpherePlaneError::Support)?;
            support_source_evaluations = support_source_evaluations
                .checked_add(report.evaluations)
                .ok_or(SourceSpherePlaneError::Invalid(
                    "support source evaluation count overflow",
                ))?;
            Ok((report.driven, report.state))
        },
        |accepted, phase, cancel| {
            if impacts.len() >= controls.motion.motion.max_events {
                return Err(SourceSpherePlaneError::Invalid(
                    "impact source receipt event budget",
                ));
            }
            let conditions = impact_conditions(accepted, phase)
                .map_err(SourceSpherePlaneError::ImpactConditions)?;
            let report = law
                .solve_preimpact_impulse(
                    accepted,
                    &phase.input,
                    &conditions,
                    controls.impact_coefficient_role,
                    controls.impulse,
                    &mut || cancel(),
                )
                .map_err(SourceSpherePlaneError::Impact)?;
            let next = report.state.clone();
            let impulse = report.impulse;
            impacts.push(SourceSpherePlaneImpactReceipt {
                phase: *phase,
                response: report,
            });
            Ok((impulse, next))
        },
        &mut cancelled,
    )?;
    if cancelled() {
        return Err(SpherePlaneDirectionalError::Mechanics(
            SpherePlaneError::Cancelled,
        ));
    }
    *history = candidate;
    Ok(SourceSpherePlaneAdvance {
        trajectory,
        impacts,
        support_source_evaluations,
    })
}
