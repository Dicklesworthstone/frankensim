//! Applied momentum impulses through the original finite-body joint graph.

use super::{
    ArticulatedJoint, Dd, DistanceGraphControls, DistanceGraphError, DistanceResidual,
    DistanceStepError, Geometry, JointRows, MAX_DISTANCE_ROWS, MassProperties, Mobility, Problem,
    RigidBodyIntegrator, RigidBodyState, Row, Rows, Vec3, Wrench, limit, norm, poll, positive,
    reserved, zeros,
};
use crate::{Gravity, MomentumImpulseReceipt, MomentumWork, momentum_difference, rounded_vector};

/// One external resultant at the native body's fixed pose.
/// Linear impulse is world-frame [N s]; angular impulse is about its COM in
/// body coordinates [N m s]. A point impulse contributes
/// `arm_body cross orientation.rotate_world_to_body(linear_world_ns)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GraphBodyImpulse {
    /// World-frame linear impulse [N s].
    pub linear_world_ns: Vec3,
    /// COM angular impulse in the body's original coordinates [N m s].
    pub angular_body_nms: Vec3,
}
impl GraphBodyImpulse {
    /// No external impulse.
    pub const ZERO: Self = Self {
        linear_world_ns: Vec3::ZERO,
        angular_body_nms: Vec3::ZERO,
    };

    /// Lift a world impulse at an original COM-relative material point.
    /// The pose is only used for the frame transform; no body is advanced.
    pub fn at_body_point(
        state: RigidBodyState,
        arm_body_m: Vec3,
        impulse_world_ns: Vec3,
    ) -> Result<Self, crate::DynamicsError> {
        if !arm_body_m.is_finite() || !impulse_world_ns.is_finite() {
            return Err(crate::DynamicsError::NonFinite("graph point impulse"));
        }
        let local = state
            .pose()
            .orientation()
            .rotate_world_to_body(impulse_world_ns);
        let angular = arm_body_m.cross(local);
        if !angular.is_finite() {
            return Err(crate::DynamicsError::NonFinite(
                "graph point angular impulse",
            ));
        }
        Ok(Self {
            linear_world_ns: impulse_world_ns,
            angular_body_nms: angular,
        })
    }
}

/// Physical admission for one fixed-configuration impulse transaction.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GraphImpulseControls {
    /// Original geometry, tangent-speed, row/body, pivot and storage bounds.
    /// Only `step.maximum_impulse_ns` is used from the time-step controls:
    /// the event has no duration, angular travel or nonlinear trials.
    pub graph: DistanceGraphControls,
    /// Maximum source/reaction work balance and ideal-reaction work defect [J].
    pub work_tolerance_j: f64,
    /// Per-body and total linear impulse representability tolerance [N s].
    pub linear_impulse_tolerance_ns: f64,
    /// Per-body and total angular impulse representability tolerance [N m s].
    pub angular_impulse_tolerance_nms: f64,
}

/// One simultaneous source-and-reaction event; all poses are unchanged.
#[derive(Debug, Clone, PartialEq)]
pub struct GraphImpulseResponse {
    /// Accepted states in original body order.
    pub states_after: Vec<RigidBodyState>,
    /// Exactly one original resultant update per body, including free bodies.
    pub native_impulses: Vec<MomentumImpulseReceipt>,
    /// Reactions in original scalar-row order. Angular rows use their existing
    /// length-scaled coordinates and physical couple bounds.
    pub multipliers_ns: Vec<f64>,
    /// Authored joint to original scalar-row map.
    pub joint_rows: Vec<JointRows>,
    /// Original geometric and tangent residuals before the event.
    pub residuals_before: Vec<DistanceResidual>,
    /// Residuals at the same poses with accepted endpoint momenta.
    pub residuals_after: Vec<DistanceResidual>,
    /// Original external sources paired with whole-event mean velocities [J].
    pub applied_work_j: f64,
    /// Ideal reactions paired with the same mean velocities [J]. Finite
    /// admitted tangency/rounding defects remain visible; this is not heat.
    pub constraint_work_j: f64,
    /// Actual kinetic change from the stored momentum increments [J].
    pub kinetic_energy_change_j: f64,
    /// Actual kinetic change minus both independently paired works [J].
    pub work_energy_residual_j: f64,
    /// Original source linear impulse, summed in world coordinates [N s].
    pub applied_linear_impulse_world_ns: Vec3,
    /// Original source spin plus orbital impulse about the world origin [N m s].
    pub applied_angular_impulse_world_nms: Vec3,
    /// Requested internal angular impulse about the world origin [N m s].
    /// A finite accepted material-point closure can leave a small real couple.
    pub constraint_angular_impulse_world_nms: Vec3,
    /// Actual stored total linear-momentum change [N s].
    pub linear_momentum_change_world_ns: Vec3,
    /// Actual stored spin-plus-orbital momentum change [N m s].
    pub angular_momentum_change_world_nms: Vec3,
    /// Actual change minus original source and reaction impulses [N s].
    pub linear_momentum_residual_world_ns: Vec3,
    /// Actual change minus original source and reaction impulses [N m s].
    pub angular_momentum_residual_world_nms: Vec3,
    /// Conservative logical owned-payload bound admitted before allocation.
    pub planned_live_bytes: usize,
}

/// Apply external impulses through rods, balls, hinges, finite planes and sliders.
///
/// At the original fixed poses the existing graph owner forms the complete
/// `G = J M^-1 J^T`. Its original normalized Cholesky solves
/// `G lambda = -J M^-1 source`. Source response is evaluated on impulse-only
/// native states: no subtraction of boosted velocities or free-source endpoint
/// tangency admission is used. An already admitted initial tangent defect is
/// not projected away; zero source leaves the motion unchanged after geometry
/// and rank admission.
///
/// Each body's external source and all its force/couple reactions are summed
/// before exactly one native momentum update. Thus shared-body cross terms
/// belong to one physical simultaneous event. Applied and constraint work use
/// the actual complete event's mean velocities, never intermediate free-source
/// or sequential-reaction velocities. No kinetic loss is inferred as heat.
///
/// Uses the original finite full inertia, row geometry, scaled mobility, grouped
/// impulse caps and tangent admission. Cost is O(B+R²+R³), with R at most 32,
/// and no generalized mass matrix or free-motion integration. Every refusal,
/// including late cancellation or undeliverable impulses, publishes no state.
///
/// This is a bilateral impulse response at already satisfied joint geometry.
/// It is not collision detection, restitution, constraint engagement, friction,
/// a prescribed joint-rate jump or a continuous-time load step.
#[allow(clippy::too_many_lines)]
pub fn apply_articulated_graph_impulses(
    states: &[RigidBodyState],
    properties: &[MassProperties],
    joints: &[ArticulatedJoint],
    sources: &[GraphBodyImpulse],
    controls: GraphImpulseControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<GraphImpulseResponse, DistanceGraphError> {
    poll(&mut cancelled)?;
    let graph = controls.graph;
    if states.len() != properties.len()
        || states.len() != sources.len()
        || [
            graph.position_tolerance_m,
            graph.velocity_tolerance_mps,
            graph.step.maximum_impulse_ns,
            controls.work_tolerance_j,
            controls.linear_impulse_tolerance_ns,
            controls.angular_impulse_tolerance_nms,
        ]
        .iter()
        .any(|&x| !positive(x))
        || !graph.minimum_scaled_pivot.is_finite()
        || !(4096.0 * f64::EPSILON..1.0).contains(&graph.minimum_scaled_pivot)
    {
        return Err(DistanceStepError::Invalid(
            "complete graph impulse inputs and positive controls",
        )
        .into());
    }
    limit("bodies", states.len() as u128, graph.maximum_bodies as u128)?;
    let maximum = graph.maximum_constraints.min(MAX_DISTANCE_ROWS);
    limit("joints", joints.len() as u128, maximum as u128)?;
    let mut bank = [Row::EMPTY; MAX_DISTANCE_ROWS];
    let mut count = 0;
    for &joint in joints {
        poll(&mut cancelled)?;
        let width = joint.row_count();
        limit("constraints", (count + width) as u128, maximum as u128)?;
        joint.write_rows(&mut bank[count..count + width]);
        count += width;
    }
    let rows = Rows::Mixed(&bank[..count]);
    use core::mem::size_of as size;
    let bytes = states.len() as u128
        * (2 * size::<RigidBodyState>()
            + size::<MomentumImpulseReceipt>()
            + size::<AccumulatedImpulse>()) as u128
        + count as u128
            * (2 * size::<Geometry>()
                + 2 * size::<DistanceResidual>()
                + size::<super::AngularRow>()
                + 2 * size::<f64>()) as u128
        + joints.len() as u128 * size::<JointRows>() as u128
        + 4 * count as u128 * count as u128 * size::<f64>() as u128;
    limit("live bytes", bytes, graph.maximum_live_bytes as u128)?;

    let mut source_states = reserved(states.len())?;
    for ((&state, &property), &source) in states.iter().zip(properties).zip(sources) {
        poll(&mut cancelled)?;
        state.validate()?;
        property.validate()?;
        source_states.push(RigidBodyState::new(
            state.pose(),
            source.linear_world_ns,
            source.angular_body_nms,
        )?);
    }
    // Only the original geometry and inverse-mass response are used.
    let problem = Problem {
        states,
        properties,
        links: rows,
        integrator: RigidBodyIntegrator::new(Gravity::ZERO),
        flow: &[] as &[Wrench],
        duration_s: 0.0,
        controls: graph,
        targets: None,
    };
    let mut before = reserved(count)?;
    let mut multipliers = zeros(count)?;
    let mut residuals_before = reserved(count)?;
    for (index, row) in rows.iter().enumerate() {
        poll(&mut cancelled)?;
        if row.bodies[0] == row.bodies[1] || row.bodies.iter().any(|&b| b >= states.len()) {
            return Err(DistanceStepError::Invalid("distinct existing joint bodies").into());
        }
        let geometry = problem.geometry(row, row.bodies.map(|b| states[b]), None)?;
        if !geometry.gap.is_finite()
            || !geometry.rate.is_finite()
            || geometry.gap.abs() > graph.position_tolerance_m
            || geometry.rate.abs() > graph.velocity_tolerance_mps
        {
            return Err(DistanceGraphError::InitialConstraint { constraint: index });
        }
        let response = problem.geometry(row, row.bodies.map(|b| source_states[b]), None)?;
        multipliers[index] = -response.rate;
        residuals_before.push(residual(geometry));
        before.push(geometry);
    }
    if let Some(constraint) = rows.grouped_violation(
        &before,
        graph.position_tolerance_m,
        graph.velocity_tolerance_mps,
    ) {
        return Err(DistanceGraphError::InitialConstraint { constraint });
    }
    if count != 0 {
        Mobility::new(&problem, states, &before, &mut cancelled)?.solve(&mut multipliers)?;
    }
    if !rows.admits_impulses(&multipliers, &before, graph.step.maximum_impulse_ns) {
        return Err(DistanceStepError::LocalBranch.into());
    }
    let mut reactions = reserved(states.len())?;
    reactions.resize(states.len(), AccumulatedImpulse::ZERO);
    for ((row, geometry), &multiplier) in rows.iter().zip(&before).zip(&multipliers) {
        poll(&mut cancelled)?;
        let world = geometry.normal.scale(multiplier);
        if !world.is_finite() {
            return Err(DistanceStepError::Invalid("graph reaction impulse").into());
        }
        for side in 0..2 {
            let body = row.bodies[side];
            let sign = if side == 0 { -1.0 } else { 1.0 };
            let local = states[body]
                .pose()
                .orientation()
                .rotate_world_to_body(world.scale(sign));
            let reaction = if row.is_angular() {
                GraphBodyImpulse {
                    linear_world_ns: Vec3::ZERO,
                    angular_body_nms: local,
                }
            } else {
                GraphBodyImpulse {
                    linear_world_ns: world.scale(sign),
                    angular_body_nms: geometry.arms_body_m[side].cross(local),
                }
            };
            reactions[body].add(AccumulatedImpulse::from_source(reaction));
        }
    }
    let mut states_after = reserved(states.len())?;
    let mut native_impulses = reserved(states.len())?;
    let mut kinetic = MomentumWork::default();
    let mut applied_work = MomentumWork::default();
    let mut constraint_work = MomentumWork::default();
    let mut applied_total = AccumulatedImpulse::ZERO;
    let mut reaction_total = AccumulatedImpulse::ZERO;
    let mut actual_total = AccumulatedImpulse::ZERO;
    let mut defect_total = AccumulatedImpulse::ZERO;
    for body in 0..states.len() {
        poll(&mut cancelled)?;
        let old = states[body];
        let property = properties[body];
        let source = AccumulatedImpulse::from_source(sources[body]);
        let reaction = reactions[body];
        let mut total = source;
        total.add(reaction);
        let event = old.apply_momentum_impulse(
            property,
            rounded_vector(total.linear),
            rounded_vector(total.angular),
        )?;
        let new = event.state_after;
        let actual = AccumulatedImpulse {
            linear: momentum_difference(old.linear_momentum_world(), new.linear_momentum_world()),
            angular: momentum_difference(old.angular_momentum_body(), new.angular_momentum_body()),
        };
        let defect = actual.minus(total);
        admit(
            norm(rounded_vector(defect.linear)),
            controls.linear_impulse_tolerance_ns,
            "unresolved body linear impulse",
        )?;
        admit(
            norm(rounded_vector(defect.angular)),
            controls.angular_impulse_tolerance_nms,
            "unresolved body angular impulse",
        )?;
        let velocity = [
            old.center_of_mass_velocity_world(property)?,
            new.center_of_mass_velocity_world(property)?,
        ];
        let omega = [
            property.angular_velocity_body_checked(old.angular_momentum_body())?,
            property.angular_velocity_body_checked(new.angular_momentum_body())?,
        ];
        for (work, impulse) in [
            (&mut kinetic, actual),
            (&mut applied_work, source),
            (&mut constraint_work, reaction),
        ] {
            work.add(velocity[0], velocity[1], impulse.linear);
            work.add(omega[0], omega[1], impulse.angular);
        }
        let mut local_defect = MomentumWork::default();
        local_defect.add(velocity[0], velocity[1], defect.linear);
        local_defect.add(omega[0], omega[1], defect.angular);
        admit(
            local_defect.total().to_f64().abs(),
            controls.work_tolerance_j,
            "unresolved body impulse work",
        )?;
        applied_total.add(source.world_about_origin(old));
        reaction_total.add(reaction.world_about_origin(old));
        actual_total.add(actual.world_about_origin(old));
        defect_total.add(defect.world_about_origin(old));
        states_after.push(new);
        native_impulses.push(event);
    }
    let mut after = reserved(count)?;
    let mut residuals_after = reserved(count)?;
    for row in rows.iter() {
        poll(&mut cancelled)?;
        let geometry = problem.geometry(row, row.bodies.map(|b| states_after[b]), None)?;
        if !geometry.gap.is_finite()
            || !geometry.rate.is_finite()
            || geometry.gap.abs() > graph.position_tolerance_m
            || geometry.rate.abs() > graph.velocity_tolerance_mps
        {
            return Err(DistanceStepError::NotConverged.into());
        }
        residuals_after.push(residual(geometry));
        after.push(geometry);
    }
    if rows
        .grouped_violation(
            &after,
            graph.position_tolerance_m,
            graph.velocity_tolerance_mps,
        )
        .is_some()
    {
        return Err(DistanceStepError::NotConverged.into());
    }
    let work_defect = kinetic.total() - applied_work.total() - constraint_work.total();
    admit(
        work_defect.to_f64().abs(),
        controls.work_tolerance_j,
        "graph impulse work balance",
    )?;
    admit(
        constraint_work.total().to_f64().abs(),
        controls.work_tolerance_j,
        "ideal graph reaction work",
    )?;
    admit(
        norm(rounded_vector(defect_total.linear)),
        controls.linear_impulse_tolerance_ns,
        "graph linear impulse balance",
    )?;
    admit(
        norm(rounded_vector(defect_total.angular)),
        controls.angular_impulse_tolerance_nms,
        "graph angular impulse balance",
    )?;
    // Refuse unrepresentable accounts even when a smaller residual is finite.
    for value in [
        kinetic.total(),
        applied_work.total(),
        constraint_work.total(),
    ] {
        if !value.to_f64().is_finite() {
            return Err(DistanceStepError::Invalid("graph impulse work inventory").into());
        }
    }
    for value in [applied_total, reaction_total, actual_total] {
        if !rounded_vector(value.linear).is_finite() || !rounded_vector(value.angular).is_finite() {
            return Err(DistanceStepError::Invalid("graph impulse momentum inventory").into());
        }
    }
    let mut joint_rows = reserved(joints.len())?;
    let mut first_row = 0;
    for &joint in joints {
        poll(&mut cancelled)?;
        let row_count = joint.row_count();
        joint_rows.push(JointRows {
            first_row,
            row_count,
        });
        first_row += row_count;
    }
    poll(&mut cancelled)?;
    Ok(GraphImpulseResponse {
        states_after,
        native_impulses,
        multipliers_ns: multipliers,
        joint_rows,
        residuals_before,
        residuals_after,
        applied_work_j: applied_work.total().to_f64(),
        constraint_work_j: constraint_work.total().to_f64(),
        kinetic_energy_change_j: kinetic.total().to_f64(),
        work_energy_residual_j: work_defect.to_f64(),
        applied_linear_impulse_world_ns: rounded_vector(applied_total.linear),
        applied_angular_impulse_world_nms: rounded_vector(applied_total.angular),
        constraint_angular_impulse_world_nms: rounded_vector(reaction_total.angular),
        linear_momentum_change_world_ns: rounded_vector(actual_total.linear),
        angular_momentum_change_world_nms: rounded_vector(actual_total.angular),
        linear_momentum_residual_world_ns: rounded_vector(defect_total.linear),
        angular_momentum_residual_world_nms: rounded_vector(defect_total.angular),
        planned_live_bytes: bytes as usize,
    })
}

fn residual(geometry: Geometry) -> DistanceResidual {
    DistanceResidual {
        gap_m: geometry.gap,
        range_rate_mps: geometry.rate,
    }
}
fn admit(value: f64, tolerance: f64, reason: &'static str) -> Result<(), DistanceGraphError> {
    if !value.is_finite() || value > tolerance {
        Err(DistanceStepError::Invalid(reason).into())
    } else {
        Ok(())
    }
}
#[derive(Clone, Copy)]
struct AccumulatedImpulse {
    linear: [Dd; 3],
    angular: [Dd; 3],
}
impl AccumulatedImpulse {
    const ZERO: Self = Self {
        linear: [Dd::ZERO; 3],
        angular: [Dd::ZERO; 3],
    };
    fn from_source(source: GraphBodyImpulse) -> Self {
        Self {
            linear: components(source.linear_world_ns).map(Dd::from_f64),
            angular: components(source.angular_body_nms).map(Dd::from_f64),
        }
    }
    fn add(&mut self, other: Self) {
        for axis in 0..3 {
            self.linear[axis] = self.linear[axis] + other.linear[axis];
            self.angular[axis] = self.angular[axis] + other.angular[axis];
        }
    }
    fn minus(self, other: Self) -> Self {
        Self {
            linear: core::array::from_fn(|d| self.linear[d] - other.linear[d]),
            angular: core::array::from_fn(|d| self.angular[d] - other.angular[d]),
        }
    }
    fn world_about_origin(self, state: RigidBodyState) -> Self {
        let mut angular = [Dd::ZERO; 3];
        for part in [
            Vec3::new(self.angular[0].hi, self.angular[1].hi, self.angular[2].hi),
            Vec3::new(self.angular[0].lo, self.angular[1].lo, self.angular[2].lo),
        ] {
            let rotated = state.pose().orientation().rotate_body_to_world(part);
            for (sum, value) in angular.iter_mut().zip(components(rotated)) {
                *sum = *sum + Dd::from_f64(value);
            }
        }
        let position = components(state.pose().position_world()).map(Dd::from_f64);
        for d in 0..3 {
            let b = (d + 1) % 3;
            let c = (d + 2) % 3;
            angular[d] = angular[d] + position[b] * self.linear[c] - position[c] * self.linear[b];
        }
        Self {
            linear: self.linear,
            angular,
        }
    }
}
fn components(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
