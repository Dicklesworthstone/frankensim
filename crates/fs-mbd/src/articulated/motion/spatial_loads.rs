//! Inertial point loads and finite axial links on a freely moving body tree.
//!
//! Forces are rebuilt at both original Lie-midpoint stages. Axial connections
//! use the true material-point separation and relative velocity, including base
//! rotation and every intervening joint. Their equal/opposite central forces
//! enter native ABA but contribute no external spatial impulse. They are finite
//! spring/damper elements, not holonomic constraints or a contact regularization.

use super::{
    ArticulatedModel, Cow, EvaluatedStage, FreeFloatingBaseState, FreeFloatingMotionError,
    FreeFloatingMotionLoads, FreeFloatingMotionState, FreeFloatingMotionStep, MomentumMap, Vec3,
    Wrench, poll, step_evaluated, validate_slice, validate_vec3,
};

/// Combined point-load/connection bound; evaluations poll every 64 entries.
pub const MAX_SPATIAL_MOTION_LOADS: usize = 16384;

/// A fixed material point in one original link frame.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MaterialPoint {
    /// Original model link index.
    pub link: usize,
    /// Position relative to that link's declared origin, in metres.
    pub point_body: Vec3,
}

/// Held inertial force and free couple acting at a moving material point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WorldPointLoad {
    /// Material application point; it moves with the original link.
    pub point: MaterialPoint,
    /// Force held in world coordinates over the step, in newtons.
    pub force_world: Vec3,
    /// Additional free couple held in world coordinates, in newton metres.
    pub couple_world: Vec3,
}

/// Bilateral central Kelvin--Voigt connection between two original links.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AxialConnection {
    /// First material endpoint.
    pub a: MaterialPoint,
    /// Second material endpoint on a distinct original link.
    pub b: MaterialPoint,
    /// Unstretched length, in metres; finite and nonnegative.
    pub rest_length_m: f64,
    /// Axial stiffness in N/m; finite and nonnegative.
    pub stiffness_n_per_m: f64,
    /// Axial damping in N s/m; finite and nonnegative.
    pub damping_ns_per_m: f64,
}

/// Loads evaluated on the actual original tree at each numerical stage.
#[derive(Debug, Clone, Copy)]
pub struct SpatialMotionLoads<'a> {
    /// Existing held joint efforts, link-frame wrenches, and uniform gravity.
    pub held: FreeFloatingMotionLoads<'a>,
    /// Additional inertial forces/couples at material points.
    pub world_points: &'a [WorldPointLoad],
    /// Internal finite connections; endpoints remain on their original links.
    pub connections: &'a [AxialConnection],
}

/// Accepted mechanics and explicit finite-connection energy accounting.
#[derive(Debug, Clone, PartialEq)]
pub struct SpatialMotionStep {
    /// Original native state and impulse/work accounting. External work and
    /// impulse exclude the internal connections; their work has its own field.
    pub motion: FreeFloatingMotionStep,
    /// Sum of k (length-rest_length)^2 / 2 at the actual supplied state, in J.
    pub stored_energy_before: f64,
    /// The same storage at the actual returned state, in J.
    pub stored_energy_after: f64,
    /// Nonnegative midpoint quadrature of c length_rate^2, in J. This is not an
    /// exact continuous loss integral or an assigned thermal destination.
    pub dissipated_energy: f64,
    /// Connection work + actual stored-energy change + quadrature loss, in J.
    /// This finite-step storage/work defect is not physical dissipation.
    pub connection_work_residual: f64,
    /// Actual change of kinetic+connection energy plus loss, minus external
    /// and actuator work, in J. No exact finite-step passivity is implied.
    pub mechanical_work_residual: f64,
}

/// Advance world-loaded articulated motion and finite internal attachments.
///
/// The original O(links) ABA and momentum reconstruction remain authoritative;
/// connections add O(points+connections) work and O(links) wrench storage. No
/// dense generalized mass matrix or separate connection solver is introduced.
/// The explicit second-order lane requires a step resolving connection modes.
/// It is not suitable for treating arbitrarily stiff links as rigid constraints.
///
/// All points, coefficients, link indices, stage loads, final geometry and
/// native joint limits must be admitted before returning an updated state.
/// Coincident connection endpoints refuse because their axial direction is
/// undefined. All caller states and load descriptions are borrowed unchanged.
pub fn step_free_floating_spatial(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    loads: SpatialMotionLoads<'_>,
    duration: f64,
    mut cancelled: impl FnMut() -> bool,
) -> Result<SpatialMotionStep, FreeFloatingMotionError> {
    poll(&mut cancelled)?;
    admit(model, loads, &mut cancelled)?;
    let evaluated = step_evaluated(
        model,
        state,
        loads.held,
        duration,
        &mut cancelled,
        |map, base, cancelled| evaluate(map, base, loads, cancelled),
    )?;
    let mut stored_energy_after = 0.0;
    if !loads.connections.is_empty() {
        let final_map = MomentumMap::new(model, &evaluated.step.state_after, &mut cancelled)?;
        for (i, connection) in loads.connections.iter().enumerate() {
            if i % 64 == 0 {
                poll(&mut cancelled)?;
            }
            let (_, length) = separation(&final_map, *connection)?;
            stored_energy_after += storage(*connection, length);
        }
    }
    let connection_work_residual = evaluated.step.internal_body_work + stored_energy_after
        - evaluated.stored_energy_before
        + evaluated.dissipated_energy;
    let mechanical_work_residual = evaluated.step.kinetic_work_residual + connection_work_residual;
    validate_slice(
        &[
            stored_energy_after,
            connection_work_residual,
            mechanical_work_residual,
        ],
        "spatial_motion.final_energy",
    )?;
    poll(&mut cancelled)?;
    Ok(SpatialMotionStep {
        motion: evaluated.step,
        stored_energy_before: evaluated.stored_energy_before,
        stored_energy_after,
        dissipated_energy: evaluated.dissipated_energy,
        connection_work_residual,
        mechanical_work_residual,
    })
}

fn admit(
    model: &ArticulatedModel,
    loads: SpatialMotionLoads<'_>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), FreeFloatingMotionError> {
    if loads
        .world_points
        .len()
        .saturating_add(loads.connections.len())
        > MAX_SPATIAL_MOTION_LOADS
    {
        return Err(FreeFloatingMotionError::InvalidLoad("spatial load count"));
    }
    let point = |p: MaterialPoint| {
        if p.link >= model.link_count() {
            return Err(FreeFloatingMotionError::InvalidLoad("material point link"));
        }
        validate_vec3(p.point_body, "spatial_motion.material_point")?;
        Ok(())
    };
    for (i, load) in loads.world_points.iter().enumerate() {
        if i % 64 == 0 {
            poll(cancelled)?;
        }
        point(load.point)?;
        validate_vec3(load.force_world, "spatial_motion.world_force")?;
        validate_vec3(load.couple_world, "spatial_motion.world_couple")?;
    }
    for (i, c) in loads.connections.iter().enumerate() {
        if i % 64 == 0 {
            poll(cancelled)?;
        }
        point(c.a)?;
        point(c.b)?;
        if c.a.link == c.b.link
            || [c.rest_length_m, c.stiffness_n_per_m, c.damping_ns_per_m]
                .iter()
                .any(|x| !x.is_finite() || *x < 0.0)
        {
            return Err(FreeFloatingMotionError::InvalidLoad("axial connection"));
        }
    }
    Ok(())
}

fn evaluate<'a>(
    map: &MomentumMap,
    base: FreeFloatingBaseState,
    loads: SpatialMotionLoads<'a>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<EvaluatedStage<'a>, FreeFloatingMotionError> {
    // The original external load account is evaluated before adding internal
    // forces. Thus central connection rounding cannot become an external source.
    let mut load = map.loads(base, loads.held, cancelled)?;
    let mut wrenches = loads.held.external_body_wrench.to_vec();
    let rotation = base.world_from_base.rotation();
    for (i, applied) in loads.world_points.iter().enumerate() {
        if i % 64 == 0 {
            poll(cancelled)?;
        }
        let point = applied.point;
        let position = position_base(map, point)?;
        let shape_velocity = shape_point_velocity(map, point)?;
        let velocity = rotation.rotate(
            base.twist_body.linear + base.twist_body.angular.cross(position) + shape_velocity,
        )?;
        let angular = rotation.rotate(
            base.twist_body.angular
                + map.shape.world_from_link[point.link]
                    .rotation()
                    .rotate(map.shape.body_twist[point.link].angular)?,
        )?;
        let force_base = rotation.inverse().rotate(applied.force_world)?;
        let couple_base = rotation.inverse().rotate(applied.couple_world)?;
        add_wrench(map, &mut wrenches, point, force_base, couple_base)?;
        load.force_world = load.force_world + applied.force_world;
        load.torque_com_world = load.torque_com_world
            + rotation
                .rotate(position - map.center_base)?
                .cross(applied.force_world)
            + applied.couple_world;
        load.body_wrench_power +=
            applied.force_world.dot(velocity) + applied.couple_world.dot(angular);
    }
    for (i, connection) in loads.connections.iter().enumerate() {
        if i % 64 == 0 {
            poll(cancelled)?;
        }
        let c = *connection;
        let (delta, length) = separation(map, c)?;
        let direction = delta.scale(1.0 / length);
        // Subtract shape velocities and use omega cross delta directly: a
        // large common translation never contaminates the relative strain rate.
        let relative = base.twist_body.angular.cross(delta) + shape_point_velocity(map, c.b)?
            - shape_point_velocity(map, c.a)?;
        let rate = direction.dot(relative);
        let magnitude =
            c.stiffness_n_per_m * (length - c.rest_length_m) + c.damping_ns_per_m * rate;
        let force = direction.scale(magnitude);
        add_wrench(map, &mut wrenches, c.a, force, Vec3::default())?;
        add_wrench(map, &mut wrenches, c.b, force.scale(-1.0), Vec3::default())?;
        load.internal_body_power -= magnitude * rate;
        load.dissipation_power += c.damping_ns_per_m * rate * rate;
        load.stored_energy += storage(c, length);
    }
    validate_vec3(load.force_world, "spatial_motion.force")?;
    validate_vec3(load.torque_com_world, "spatial_motion.torque")?;
    validate_slice(
        &[
            load.body_wrench_power,
            load.internal_body_power,
            load.dissipation_power,
            load.stored_energy,
        ],
        "spatial_motion.stage_energy",
    )?;
    Ok(EvaluatedStage {
        body_wrench: Cow::Owned(wrenches),
        load,
    })
}

fn add_wrench(
    map: &MomentumMap,
    wrenches: &mut [Wrench],
    point: MaterialPoint,
    force_base: Vec3,
    couple_base: Vec3,
) -> Result<(), FreeFloatingMotionError> {
    let inverse = map.shape.world_from_link[point.link].rotation().inverse();
    let force = inverse.rotate(force_base)?;
    let torque = point.point_body.cross(force) + inverse.rotate(couple_base)?;
    let old = wrenches[point.link];
    wrenches[point.link] = Wrench::new(old.torque + torque, old.force + force);
    Ok(())
}

fn position_base(map: &MomentumMap, point: MaterialPoint) -> Result<Vec3, FreeFloatingMotionError> {
    Ok(map.shape.world_from_link[point.link].transform_point(point.point_body)?)
}

fn shape_point_velocity(
    map: &MomentumMap,
    point: MaterialPoint,
) -> Result<Vec3, FreeFloatingMotionError> {
    let twist = map.shape.body_twist[point.link];
    Ok(map.shape.world_from_link[point.link]
        .rotation()
        .rotate(twist.linear + twist.angular.cross(point.point_body))?)
}

fn separation(
    map: &MomentumMap,
    connection: AxialConnection,
) -> Result<(Vec3, f64), FreeFloatingMotionError> {
    let delta = position_base(map, connection.b)? - position_base(map, connection.a)?;
    let length = delta.norm();
    if !length.is_finite() || length <= 0.0 {
        return Err(FreeFloatingMotionError::InvalidLoad(
            "coincident or nonfinite connection",
        ));
    }
    Ok((delta, length))
}

fn storage(connection: AxialConnection, length: f64) -> f64 {
    let extension = length - connection.rest_length_m;
    0.5 * connection.stiffness_n_per_m * extension * extension
}
