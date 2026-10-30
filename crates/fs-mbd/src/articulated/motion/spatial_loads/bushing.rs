//! Finite, anisotropic translational and rotational attachments on the native tree.
//!
//! For material frames A/B, `q = R_A^T (x_B-x_A) - q_rest` and
//! `r = Log(R_A^T R_B)` store `U = (q.Kt.q + r.Kr.r)/2`. Damping is
//! nonnegative in these exact strain rates. The physical couple uses the
//! inverse-transpose SO(3) left Jacobian; treating `r_dot` as relative angular
//! velocity would be wrong for anisotropic finite rotations. The log differential
//! follows Sola, Deray and Atchuthan, *A micro Lie theory*, equations 74 and 146:
//! <https://arxiv.org/abs/1812.01537>.
//!
//! The A-side couple also contains `(x_B-x_A) cross F_A`, required by the
//! rotating A-frame translation law. Thus both total force and total spatial
//! moment vanish, and power is exactly `-dU/dt - loss` in real arithmetic.
//! A prescribed rest orientation is represented by the two attachment frames.
//! This is finite compliance, not a rigid six-axis joint or unilateral contact.

use super::super::{
    ArticulatedModel, EvaluatedStage, FreeFloatingMotionError, FreeFloatingMotionState,
    MAX_FREE_FLOATING_STEP_LINKS, MomentumMap, Se3, StageLoad, Vec3, Wrench, poll, step_evaluated,
    validate_external_wrenches, validate_slice, validate_twist, validate_vec3,
};
use super::implicit::{ImplicitSpatialConfig, ImplicitSpatialError, ImplicitSpatialStep};
use super::{
    MaterialPoint, SpatialMotionLoads, SpatialMotionStep, add_wrench, finish_with_bushings,
    position_base, shape_point_velocity,
};
use fs_ga::{So3Tangent, so3_left_jacobian_inverse};

/// Principal-log chart boundary. Active rotational elements refuse angles at
/// or above pi minus one microradian; no quaternion sign can select a smooth
/// continuation through pi for this single-well potential.
pub const MAX_BUSHING_ROTATION_RADIANS: f64 = core::f64::consts::PI - 1.0e-6;

/// A fixed pose of a material attachment relative to its original link.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MaterialFrame {
    /// Original model link index.
    pub link: usize,
    /// Maps attachment coordinates into the link frame, in metres/radians.
    pub link_from_attachment: Se3,
}

/// Finite Kelvin--Voigt bushing with diagonal, anisotropic strain coefficients.
///
/// All stiffness/damping components must be finite and nonnegative; zero axes
/// are allowed. Translation coefficients are expressed in attachment A axes;
/// rotation coefficients act on the principal rotation-vector coordinates.
/// This permits a freely moving direction or torsion-only/translational-only
/// element without a fake small stiffness. Cross-axis material coupling is not
/// inferred. Attachment frames may be coincident, but links must be distinct.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpatialBushing {
    /// First material frame.
    pub a: MaterialFrame,
    /// Second material frame.
    pub b: MaterialFrame,
    /// Unstrained B-origin position in A, in metres.
    pub rest_translation_a: Vec3,
    /// Diagonal translational stiffness, in N/m.
    pub translation_stiffness: Vec3,
    /// Diagonal translational damping, in N s/m.
    pub translation_damping: Vec3,
    /// Diagonal rotation-vector stiffness, in N m/rad.
    pub rotation_stiffness: Vec3,
    /// Diagonal rotation-vector damping, in N m s/rad.
    pub rotation_damping: Vec3,
}

/// Instantaneous bushing loads and constitutive energy, before time integration.
#[derive(Debug, Clone, PartialEq)]
pub struct BushingForces {
    /// Per-link internal wrenches about original link origins in link frames.
    pub body_wrench: Vec<Wrench>,
    /// Sum of the actual quadratic strain storage, in joules.
    pub stored_energy: f64,
    /// Work rate of these internal wrenches, in watts.
    pub internal_power: f64,
    /// Nonnegative quadratic damping loss, in watts.
    pub dissipation_power: f64,
}

/// Evaluate finite bushings on the actual original articulated state.
///
/// Relative strains/rates are formed in the shape frame before adding any base
/// motion. Common rigid rotation, translation and boosts therefore cancel
/// algebraically instead of contaminating strain through large world numbers.
/// Work is O(links+bushings), with cancellation between native passes and every
/// 64 elements. Inputs are immutable on every success or refusal.
///
/// # Errors
/// Refuses the original model/state admission, excessive work, invalid material
/// data, an active rotational chart boundary, nonfinite arithmetic or cancellation.
pub fn evaluate_bushing_forces(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    bushings: &[SpatialBushing],
    mut cancelled: impl FnMut() -> bool,
) -> Result<BushingForces, FreeFloatingMotionError> {
    poll(&mut cancelled)?;
    if model.link_count() > MAX_FREE_FLOATING_STEP_LINKS {
        return Err(FreeFloatingMotionError::TooManyLinks {
            got: model.link_count(),
            maximum: MAX_FREE_FLOATING_STEP_LINKS,
        });
    }
    model.free_floating_complexity()?;
    validate_twist(state.base.twist_body, "bushing.base_twist")?;
    admit_bushings(model, bushings, &mut cancelled)?;
    let map = MomentumMap::new(model, state, &mut cancelled)?;
    let mut body_wrench = vec![Wrench::new(Vec3::default(), Vec3::default()); model.link_count()];
    let mut load = StageLoad {
        force_world: Vec3::default(),
        torque_com_world: Vec3::default(),
        body_wrench_power: 0.0,
        internal_body_power: 0.0,
        dissipation_power: 0.0,
        stored_energy: 0.0,
    };
    accumulate(&map, bushings, &mut body_wrench, &mut load, &mut cancelled)?;
    poll(&mut cancelled)?;
    Ok(BushingForces {
        body_wrench,
        stored_energy: load.stored_energy,
        internal_power: load.internal_body_power,
        dissipation_power: load.dissipation_power,
    })
}

/// Advance axial connections, spatial loads and bushings in one explicit Lie
/// midpoint step. Every internal law is reevaluated at the actual stage state.
/// All storage/loss fields in the result include both axial links and bushings.
///
/// Resolve the fastest compliant mode by timestep refinement. Rotational chart
/// admission is checked at the initial, midpoint and final state; these checks
/// do not enclose the continuous path or permit traversal through the pi cut.
///
/// # Errors
/// Returns the original motion/load refusals and bushing material/chart refusals.
pub fn step_free_floating_bushings(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    loads: SpatialMotionLoads<'_>,
    bushings: &[SpatialBushing],
    duration: f64,
    mut cancelled: impl FnMut() -> bool,
) -> Result<SpatialMotionStep, FreeFloatingMotionError> {
    poll(&mut cancelled)?;
    admit(model, loads, bushings, &mut cancelled)?;
    let evaluated = step_evaluated(
        model,
        state,
        loads.held,
        duration,
        &mut cancelled,
        |map, base, cancelled| evaluate(map, base, loads, bushings, cancelled),
    )?;
    finish_with_bushings(model, evaluated, loads, bushings, &mut cancelled)
}

/// Advance the same bushing law with the native simultaneous implicit midpoint
/// solve. This admits stiff linear compliance without the explicit stability
/// bound, while retaining all original inertias, base recoil and momenta.
/// Nonlinear finite-step energy defects are reported, not hidden as damping.
/// The implicit rule is not L-stable; refine for phase accuracy and geometry.
///
/// # Errors
/// Includes bounded-solve, cancellation, original-state and bushing refusals.
pub fn step_free_floating_bushings_implicit(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    loads: SpatialMotionLoads<'_>,
    bushings: &[SpatialBushing],
    duration: f64,
    config: ImplicitSpatialConfig,
    cancelled: impl FnMut() -> bool,
) -> Result<ImplicitSpatialStep, ImplicitSpatialError> {
    super::implicit::step_with_bushings(model, state, loads, bushings, duration, config, cancelled)
}

pub(super) fn admit(
    model: &ArticulatedModel,
    loads: SpatialMotionLoads<'_>,
    bushings: &[SpatialBushing],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), FreeFloatingMotionError> {
    super::admit(model, loads, cancelled)?;
    if loads
        .world_points
        .len()
        .saturating_add(loads.connections.len())
        .saturating_add(bushings.len())
        > super::MAX_SPATIAL_MOTION_LOADS
    {
        return Err(FreeFloatingMotionError::InvalidLoad(
            "combined spatial load count",
        ));
    }
    admit_bushings(model, bushings, cancelled)
}

fn admit_bushings(
    model: &ArticulatedModel,
    bushings: &[SpatialBushing],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), FreeFloatingMotionError> {
    if bushings.len() > super::MAX_SPATIAL_MOTION_LOADS {
        return Err(FreeFloatingMotionError::InvalidLoad("bushing count"));
    }
    for (i, b) in bushings.iter().enumerate() {
        if i % 64 == 0 {
            poll(cancelled)?;
        }
        if b.a.link >= model.link_count() || b.b.link >= model.link_count() || b.a.link == b.b.link
        {
            return Err(FreeFloatingMotionError::InvalidLoad(
                "bushing attachment links",
            ));
        }
        validate_vec3(b.rest_translation_a, "bushing.rest_translation")?;
        for v in [
            b.translation_stiffness,
            b.translation_damping,
            b.rotation_stiffness,
            b.rotation_damping,
        ] {
            if [v.x, v.y, v.z].iter().any(|x| !x.is_finite() || *x < 0.0) {
                return Err(FreeFloatingMotionError::InvalidLoad(
                    "nonnegative finite bushing coefficients",
                ));
            }
        }
    }
    Ok(())
}

pub(super) fn evaluate<'a>(
    map: &MomentumMap,
    base: super::super::FreeFloatingBaseState,
    loads: SpatialMotionLoads<'a>,
    bushings: &[SpatialBushing],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<EvaluatedStage<'a>, FreeFloatingMotionError> {
    let mut stage = super::evaluate(map, base, loads, cancelled)?;
    if !bushings.is_empty() {
        accumulate(
            map,
            bushings,
            stage.body_wrench.to_mut(),
            &mut stage.load,
            cancelled,
        )?;
    }
    Ok(stage)
}

struct Strain {
    translation: Vec3,
    rotation: Vec3,
    frame_a_rotation: fs_ga::So3,
    delta_base: Vec3,
}

fn strain(map: &MomentumMap, b: SpatialBushing) -> Result<Strain, FreeFloatingMotionError> {
    let a = map.shape.world_from_link[b.a.link]
        .rotation()
        .compose(b.a.link_from_attachment.rotation())?;
    let delta_base = position_base(map, point(b.b))? - position_base(map, point(b.a))?;
    let translation = a.inverse().rotate(delta_base)? - b.rest_translation_a;
    let mut rotation = Vec3::default();
    if b.rotation_stiffness != Vec3::default() || b.rotation_damping != Vec3::default() {
        let r_b = map.shape.world_from_link[b.b.link]
            .rotation()
            .compose(b.b.link_from_attachment.rotation())?;
        rotation = a.inverse().compose(r_b)?.log().angular;
        if !rotation.norm().is_finite() || rotation.norm() >= MAX_BUSHING_ROTATION_RADIANS {
            return Err(FreeFloatingMotionError::InvalidLoad(
                "bushing rotation reaches principal-log chart boundary",
            ));
        }
    }
    Ok(Strain {
        translation,
        rotation,
        frame_a_rotation: a,
        delta_base,
    })
}

fn accumulate(
    map: &MomentumMap,
    bushings: &[SpatialBushing],
    wrenches: &mut [Wrench],
    load: &mut StageLoad,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), FreeFloatingMotionError> {
    for (i, b) in bushings.iter().copied().enumerate() {
        if i % 64 == 0 {
            poll(cancelled)?;
        }
        let s = strain(map, b)?;
        let a_inverse = s.frame_a_rotation.inverse();
        let angular_a = map.shape.world_from_link[b.a.link]
            .rotation()
            .rotate(map.shape.body_twist[b.a.link].angular)?;
        let angular_b = map.shape.world_from_link[b.b.link]
            .rotation()
            .rotate(map.shape.body_twist[b.b.link].angular)?;
        // Cancel both common base translation AND common base rotation exactly.
        let translation_rate = a_inverse.rotate(
            shape_point_velocity(map, point(b.b))?
                - shape_point_velocity(map, point(b.a))?
                - angular_a.cross(s.delta_base),
        )?;
        let inverse_jacobian = so3_left_jacobian_inverse(So3Tangent::new(s.rotation))?;
        let rotation_rate = inverse_jacobian.apply(a_inverse.rotate(angular_b - angular_a)?);
        let force_a = diagonal(b.translation_stiffness, s.translation)
            + diagonal(b.translation_damping, translation_rate);
        let log_force = diagonal(b.rotation_stiffness, s.rotation)
            + diagonal(b.rotation_damping, rotation_rate);
        let force = s.frame_a_rotation.rotate(force_a)?;
        let couple = s
            .frame_a_rotation
            .rotate(inverse_jacobian.transpose().apply(log_force))?;
        add_wrench(
            map,
            wrenches,
            point(b.a),
            force,
            s.delta_base.cross(force) + couple,
        )?;
        add_wrench(
            map,
            wrenches,
            point(b.b),
            force.scale(-1.0),
            couple.scale(-1.0),
        )?;
        load.internal_body_power -= force_a.dot(translation_rate) + log_force.dot(rotation_rate);
        load.dissipation_power += translation_rate
            .dot(diagonal(b.translation_damping, translation_rate))
            + rotation_rate.dot(diagonal(b.rotation_damping, rotation_rate));
        load.stored_energy += storage(b, &s);
    }
    validate_external_wrenches(wrenches)?;
    validate_slice(
        &[
            load.internal_body_power,
            load.dissipation_power,
            load.stored_energy,
        ],
        "bushing.energy",
    )?;
    Ok(())
}

pub(super) fn stored_energy(
    map: &MomentumMap,
    bushings: &[SpatialBushing],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, FreeFloatingMotionError> {
    let mut energy = 0.0;
    for (i, b) in bushings.iter().copied().enumerate() {
        if i % 64 == 0 {
            poll(cancelled)?;
        }
        energy += storage(b, &strain(map, b)?);
    }
    validate_slice(&[energy], "bushing.final_storage")?;
    Ok(energy)
}

fn storage(b: SpatialBushing, s: &Strain) -> f64 {
    0.5 * (s
        .translation
        .dot(diagonal(b.translation_stiffness, s.translation))
        + s.rotation.dot(diagonal(b.rotation_stiffness, s.rotation)))
}

fn diagonal(d: Vec3, x: Vec3) -> Vec3 {
    Vec3::new(d.x * x.x, d.y * x.y, d.z * x.z)
}

fn point(f: MaterialFrame) -> MaterialPoint {
    MaterialPoint {
        link: f.link,
        point_body: f.link_from_attachment.translation(),
    }
}
