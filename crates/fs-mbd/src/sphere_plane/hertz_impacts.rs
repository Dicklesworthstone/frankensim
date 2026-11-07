//! Finite-patch shear and spin impulses through the original sphere trajectory.
//!
//! The instantaneous law uses one supplied Hertz pressure-impulse distribution
//! and one coupled friction budget. The sphere's actual isotropic inertia,
//! radius and contact arm determine every recoil coordinate. Normal response,
//! event timing, pose transport, motor/force work and moving-plane work remain
//! with the original engine. The normal and tangent/spin mobilities decouple
//! for this geometry; no general coupled-body impact claim is made.
//!
//! The effective impact radius and pressure-impulse shape are explicit source
//! assumptions, not a footprint inferred from rigid sphere geometry. No impact
//! duration, pressure, travel, normal compliance or heat partition is invented.
//! Persistent support uses the separately supplied point-Coulomb law. Its
//! zero torsional resistance does not extrapolate the transient impact patch.

use super::*;
use core::convert::Infallible;
use fs_tribo::sliding_patch::contact::impulse::{
    HertzImpulseControls, HertzImpulseError, HertzImpulseInput, HertzImpulsePatch,
    HertzImpulseReport,
};

/// Independent original motion and native impulse-solve limits.
#[derive(Clone, Copy, Debug)]
pub struct SpherePlaneHertzImpactControls {
    pub motion: SpherePlaneConfig,
    pub impulse: HertzImpulseControls,
}

/// Actual event state supplied to an impact source before any impulse is applied.
#[derive(Clone, Copy, Debug)]
pub struct SpherePlaneHertzImpactPhase {
    pub time_s: f64,
    pub body: RigidBodyState,
    pub contact_point_world_m: Vec3,
    pub normal_world: Vec3,
    /// Orthonormal coordinate axes for the reported shear components.
    pub tangent_axes_world: [Vec3; 2],
    pub normal_velocity_before_m_s: f64,
    pub normal_velocity_after_m_s: f64,
    pub normal_impulse_ns: f64,
    pub slip_velocity_before_m_s: [f64; 2],
    /// The prescribed translating plane has no angular velocity.
    pub spin_before_rad_s: f64,
    pub tangent_mobility_per_kg: f64,
    pub spin_mobility_per_kg_m2: f64,
}

/// Original constitutive event and its real mechanical inputs, in time order.
#[derive(Clone, Copy, Debug)]
pub struct SpherePlaneHertzImpactEvent {
    pub phase: SpherePlaneHertzImpactPhase,
    pub patch: HertzImpulsePatch,
    pub input: HertzImpulseInput,
    pub response: HertzImpulseReport,
}

/// Full native trajectory plus the bounded list of finite-patch impact receipts.
#[derive(Clone, Debug)]
pub struct SpherePlaneHertzImpactAdvance {
    pub motion: SpherePlaneAdvance,
    pub impacts: Vec<SpherePlaneHertzImpactEvent>,
    /// Directly accumulated normal plus shear/spin impact loss [J]. This is
    /// already included in motion.dissipated_energy; it is not additional heat.
    pub impact_loss_j: f64,
}

/// Every failure preserves the accepted source history and publishes no motion.
#[derive(Debug)]
pub enum SpherePlaneHertzImpactError<E = Infallible> {
    Mechanics(SpherePlaneError),
    Law(HertzImpulseError),
    Source(E),
}
impl<E> From<SpherePlaneError> for SpherePlaneHertzImpactError<E> {
    fn from(error: SpherePlaneError) -> Self {
        Self::Mechanics(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for SpherePlaneHertzImpactError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite-patch sphere/plane impacts: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for SpherePlaneHertzImpactError<E> {}

/// Advance the real sphere with one explicitly declared Hertz impulse law.
///
/// The held world couple acts about the COM during flight and support. Every
/// detected impact receives the joint sliding/spin impulse, so its changed
/// spin also drives the subsequent original attitude trajectory. Support
/// friction remains a separate static/kinetic point law. The scene's default
/// impact coefficient is replaced only at the actual impact events.
#[allow(clippy::too_many_arguments)]
pub fn advance_sphere_plane_hertz_impacts(
    scene: &SpherePlane,
    state: SpherePlaneState,
    torque_world_nm: Vec3,
    support_friction: &FrictionLaw,
    impact_patch: &HertzImpulsePatch,
    duration_s: f64,
    controls: SpherePlaneHertzImpactControls,
    cancelled: impl FnMut() -> bool,
) -> Result<SpherePlaneHertzImpactAdvance, SpherePlaneHertzImpactError> {
    advance_sphere_plane_hertz_impacts_with(
        scene,
        state,
        torque_world_nm,
        support_friction,
        &mut (),
        duration_s,
        controls,
        |_, _, _| Ok((*impact_patch, ())),
        cancelled,
    )
}

/// Query an explicit radius/impulse-coefficient source at each detected event.
///
/// The callback receives actual incoming slip/spin, the independently admitted
/// normal impulse, event time and the original body. It returns the original
/// Hertz impulse law and a candidate source history; the native solver owns
/// the impulse response. No finite-time material query or normal pressure is
/// manufactured. All source updates remain private until the complete motion
/// passes geometry, momentum, work and cancellation checks. A later event's
/// refusal cannot publish the history from an earlier private impact.
#[allow(clippy::too_many_arguments)]
pub fn advance_sphere_plane_hertz_impacts_with<S: Clone, E>(
    scene: &SpherePlane,
    state: SpherePlaneState,
    torque_world_nm: Vec3,
    support_friction: &FrictionLaw,
    history: &mut S,
    duration_s: f64,
    controls: SpherePlaneHertzImpactControls,
    mut resolve: impl FnMut(
        &S,
        &SpherePlaneHertzImpactPhase,
        &mut dyn FnMut() -> bool,
    ) -> Result<(HertzImpulsePatch, S), E>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<SpherePlaneHertzImpactAdvance, SpherePlaneHertzImpactError<E>> {
    poll(&mut cancelled)?;
    controls
        .impulse
        .validate()
        .map_err(SpherePlaneHertzImpactError::Law)?;
    let mut candidate = history.clone();
    let mut source_error = None;
    let mut impacts = Vec::new();
    let motion = advance_with_contact(
        scene,
        state,
        torque_world_nm,
        support_friction,
        duration_s,
        controls.motion,
        None,
        &mut cancelled,
        |_, _, _, _, _, _, _| Ok(None),
        |engine, body, time, normal_before, normal, slip, cancel| {
            let phase = SpherePlaneHertzImpactPhase {
                time_s: time,
                body,
                contact_point_world_m: body.pose().position_world().add(engine.arm),
                normal_world: engine.n,
                tangent_axes_world: engine.tangent,
                normal_velocity_before_m_s: normal_before,
                normal_velocity_after_m_s: normal.relative_normal_speed_after,
                normal_impulse_ns: normal.normal_impulse,
                slip_velocity_before_m_s: slip,
                spin_before_rad_s: finite(world_spin(body).dot(engine.n) / engine.inertia)?,
                tangent_mobility_per_kg: engine.inverse_tangent_mass,
                spin_mobility_per_kg_m2: engine.inertia.recip(),
            };
            let (patch, next) = match resolve(&candidate, &phase, &mut || cancel()) {
                Ok(result) => result,
                Err(error) => {
                    source_error = Some(SpherePlaneHertzImpactError::Source(error));
                    return Err(SpherePlaneError::Unresolved("Hertz impact source refused"));
                }
            };
            poll(cancel)?;
            let radius = patch.radius_m();
            let spin_response = finite((radius / engine.inertia) * radius)?;
            if spin_response <= 0.0 {
                return Err(SpherePlaneError::Invalid(
                    "unrepresentable Hertz spin response",
                ));
            }
            let w = engine.inverse_tangent_mass;
            let input = HertzImpulseInput {
                normal_impulse_ns: normal.normal_impulse,
                velocity_before_m_s: slip,
                spin_before_rad_s: phase.spin_before_rad_s,
                scaled_mobility_per_kg: [[w, 0.0, 0.0], [0.0, w, 0.0], [0.0, 0.0, spin_response]],
            };
            let response = match patch.solve(&input, controls.impulse, &mut || cancel()) {
                Ok(result) => result,
                Err(error) => {
                    source_error = Some(SpherePlaneHertzImpactError::Law(error));
                    return Err(SpherePlaneError::Unresolved("Hertz impact law refused"));
                }
            };
            let accepted = admit_impact(&phase, radius, &response, controls.impulse)?;
            poll(cancel)?;
            impacts.push(SpherePlaneHertzImpactEvent {
                phase,
                patch,
                input,
                response,
            });
            candidate = next;
            Ok(Some(accepted))
        },
    );
    if let Some(error) = source_error {
        return Err(error);
    }
    let (motion, impact_loss_j) = motion?;
    poll(&mut cancelled)?;
    *history = candidate;
    Ok(SpherePlaneHertzImpactAdvance {
        motion,
        impacts,
        impact_loss_j,
    })
}

/// Reconstruct the unscaled conjugate ports independently of radius conversion
/// and the law's inverse-mobility energy metric. The original engine then checks
/// actual body endpoint rates and complete world-origin momentum and energy.
fn admit_impact(
    phase: &SpherePlaneHertzImpactPhase,
    radius: f64,
    response: &HertzImpulseReport,
    controls: HertzImpulseControls,
) -> Result<TangentialImpactAdvance, SpherePlaneError> {
    let before = phase.slip_velocity_before_m_s;
    let w = phase.tangent_mobility_per_kg;
    let after = [
        finite(before[0] + w * response.impulse_ns[0])?,
        finite(before[1] + w * response.impulse_ns[1])?,
    ];
    let spin_after = finite(
        phase.spin_before_rad_s + phase.spin_mobility_per_kg_m2 * response.spin_impulse_nms,
    )?;
    let middle = [0, 1].map(|i| 0.5 * before[i] + 0.5 * after[i]);
    let middle_spin = 0.5 * phase.spin_before_rad_s + 0.5 * spin_after;
    let physical_loss = finite(
        -response.impulse_ns[0] * middle[0]
            - response.impulse_ns[1] * middle[1]
            - response.spin_impulse_nms * middle_spin,
    )?;
    let kinetic_change = finite(
        ((after[0] - before[0]) * middle[0] + (after[1] - before[1]) * middle[1]) / w
            + (spin_after - phase.spin_before_rad_s) * middle_spin / phase.spin_mobility_per_kg_m2,
    )?;
    if [0, 1].into_iter().any(|i| {
        (after[i] - response.velocity_after_m_s[i]).abs() > controls.velocity_tolerance_m_s
    }) || radius * (spin_after - response.spin_after_rad_s).abs()
        > controls.velocity_tolerance_m_s
        || physical_loss < -controls.energy_tolerance_j
        || (physical_loss - response.dissipated_j).abs() > controls.energy_tolerance_j
        || (kinetic_change - response.kinetic_energy_change_j).abs() > controls.energy_tolerance_j
        || (kinetic_change + physical_loss).abs() > controls.energy_tolerance_j
    {
        return Err(SpherePlaneError::Unresolved(
            "native Hertz impact work or endpoint",
        ));
    }
    Ok(TangentialImpactAdvance {
        impulse: response.impulse_ns,
        spin_impulse: response.spin_impulse_nms,
        velocity_after: after,
        dissipated: response.dissipated_j,
    })
}
