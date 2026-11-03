//! Simultaneous material-point closures between independently moving native bodies.
//!
//! The original pair constraint supplies positive-distance geometry. Spherical
//! joints add three Cartesian closure rows through [`advance_joint_graph`]. All
//! scalar reactions share the full native body-point impulse owner and mobility
//! `J M^-1 J^T`; a chain is never approximated by independent pair solves.
//! Central rod reactions conserve internal angular momentum; spherical rows
//! retain and report the small couple from accepted finite closure residuals.
//! Native free steps, actual impulses and finite constraint work remain visible,
//! and poses are not repaired.

use super::{
    DistanceStepControls, DistanceStepError, Geometry, RigidDistanceConstraint, norm, poll,
    positive,
};
use crate::{
    ActionReactionImpulseReceipt, DynamicsError, MassProperties, RigidBodyIntegrator,
    RigidBodyState, StepReceipt, Vec3, Wrench,
};
use fs_la::factor::{Cholesky, cholesky};
use fs_math::dd::Dd;

#[path = "graph/driven.rs"]
mod driven;
pub use driven::{
    DistanceDrive, DrivenDistanceGraphStep, DrivenDistanceLinkStep, advance_driven_distance_graph,
};

#[path = "graph/joints.rs"]
mod joints;
pub use joints::{GraphJoint, JointGraphStep, JointRows, SphericalJoint, advance_joint_graph};

/// Maximum simultaneous independent scalar rows in the bounded dense solve.
pub const MAX_DISTANCE_ROWS: usize = 32;

/// One positive bilateral rod, with endpoints in the original COM body frames.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DistanceLink {
    /// Distinct indices into the supplied native body/property/wrench slices.
    pub bodies: [usize; 2],
    /// Body-fixed material points, in metres from the respective COMs.
    pub points_body_m: [Vec3; 2],
    /// Fixed positive separation [m]; tension and compression are admitted.
    pub length_m: f64,
}

// Typed scalar rows feed one numerical owner. Existing distance slices are
// borrowed without flattening or allocation; only mixed joints need a bounded
// stack bank of rows. A Cartesian row never masquerades as a zero-length rod.
#[derive(Clone, Copy)]
enum RowKind {
    Distance(f64),
    Cartesian(usize),
}
#[derive(Clone, Copy)]
struct Row {
    bodies: [usize; 2],
    points_body_m: [Vec3; 2],
    kind: RowKind,
}
impl Row {
    const EMPTY: Self = Self {
        bodies: [0; 2],
        points_body_m: [Vec3::ZERO; 2],
        kind: RowKind::Distance(1.0),
    };
    fn distance(link: DistanceLink) -> Self {
        Self {
            bodies: link.bodies,
            points_body_m: link.points_body_m,
            kind: RowKind::Distance(link.length_m),
        }
    }
}
#[derive(Clone, Copy)]
enum Rows<'a> {
    Distances(&'a [DistanceLink]),
    Mixed(&'a [Row]),
}
impl<'a> Rows<'a> {
    fn len(self) -> usize {
        match self {
            Self::Distances(links) => links.len(),
            Self::Mixed(rows) => rows.len(),
        }
    }
    fn get(self, index: usize) -> Row {
        match self {
            Self::Distances(links) => Row::distance(links[index]),
            Self::Mixed(rows) => rows[index],
        }
    }
    fn iter(self) -> impl ExactSizeIterator<Item = Row> + 'a {
        (0..self.len()).map(move |i| self.get(i))
    }
    fn admits_impulses(self, impulses: &[f64], maximum: f64) -> bool {
        let mut i = 0;
        while i < self.len() {
            let (magnitude, width) = if matches!(self.get(i).kind, RowKind::Cartesian(0)) {
                (impulses[i].hypot(impulses[i + 1]).hypot(impulses[i + 2]), 3)
            } else {
                (impulses[i].abs(), 1)
            };
            if !magnitude.is_finite() || magnitude > maximum {
                return false;
            }
            i += width;
        }
        true
    }
    fn error(self, geometry: &[Geometry]) -> f64 {
        let mut i = 0;
        let mut error: f64 = 0.0;
        while i < self.len() {
            let (magnitude, width) = if matches!(self.get(i).kind, RowKind::Cartesian(0)) {
                (
                    geometry[i]
                        .gap
                        .hypot(geometry[i + 1].gap)
                        .hypot(geometry[i + 2].gap),
                    3,
                )
            } else {
                (geometry[i].gap.abs(), 1)
            };
            error = error.max(magnitude);
            i += width;
        }
        error
    }
    fn spherical_violation(
        self,
        geometry: &[Geometry],
        position: f64,
        velocity: f64,
    ) -> Option<usize> {
        for (i, row) in self.iter().enumerate() {
            if matches!(row.kind, RowKind::Cartesian(0))
                && (geometry[i]
                    .gap
                    .hypot(geometry[i + 1].gap)
                    .hypot(geometry[i + 2].gap)
                    > position
                    || geometry[i]
                        .rate
                        .hypot(geometry[i + 1].rate)
                        .hypot(geometry[i + 2].rate)
                        > velocity)
            {
                return Some(i);
            }
        }
        None
    }
}

/// Local physical and resource bounds, admitted before transaction allocation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DistanceGraphControls {
    /// Original impulse, angular-travel, nonlinear-iteration and backtrack caps.
    pub step: DistanceStepControls,
    /// Absolute initial/final closure tolerance [m], below every rod length.
    /// Spherical joints use their complete Euclidean separation norm.
    pub position_tolerance_m: f64,
    /// Absolute initial/final range-rate tolerance [m/s], or the Euclidean
    /// relative point-speed tolerance for a spherical joint.
    pub velocity_tolerance_mps: f64,
    /// Squared Cholesky-pivot floor after unit-diagonal scaling, in
    /// `[4096*EPSILON, 1)`. This admits regular rows, not a certified condition bound.
    pub minimum_scaled_pivot: f64,
    /// Maximum borrowed bodies, including unconstrained bodies.
    pub maximum_bodies: usize,
    /// Maximum simultaneous scalar rows, at most [`MAX_DISTANCE_ROWS`]. A rod
    /// consumes one row; a spherical joint in [`advance_joint_graph`] consumes three.
    pub maximum_constraints: usize,
    /// Maximum complete nonlinear trial evaluations, including the initial one.
    pub maximum_trial_evaluations: u32,
    /// Bound on logical owned payload, excluding borrowed inputs, stack storage,
    /// allocator metadata and capacity rounding. Native Cholesky allocation
    /// semantics remain unchanged; transaction buffers reserve fallibly.
    pub maximum_live_bytes: usize,
}

/// Actual accepted geometric and hidden-constraint residual of one scalar row.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DistanceResidual {
    /// Endpoint length minus the authored distance [m], or the selected world
    /// Cartesian separation component of a spherical-joint row.
    pub gap_m: f64,
    /// Endpoint range rate minus its prescribed value [m/s], or the selected
    /// relative point-velocity component of a spherical-joint row.
    pub range_rate_mps: f64,
}

/// One atomic native-body graph step, with every actual reaction event.
#[derive(Debug, Clone, PartialEq)]
pub struct DistanceGraphStep {
    /// Accepted bodies in original body-slice order.
    pub states_after: Vec<RigidBodyState>,
    /// Initial native pair impulses in original scalar-row order.
    pub position_impulses: Vec<ActionReactionImpulseReceipt>,
    /// Original unconstrained native steps between the two impulse stages.
    pub free_steps: Vec<StepReceipt>,
    /// Endpoint native pair impulses in original scalar-row order.
    pub velocity_impulses: Vec<ActionReactionImpulseReceipt>,
    /// Sum of actual initial and endpoint pair work [J], including all cross
    /// terms on shared bodies. Individual pair works use the reported order.
    pub constraint_work_j: f64,
    /// Actual endpoint residuals in original scalar-row order.
    pub residuals: Vec<DistanceResidual>,
    /// Accepted coupled position corrections.
    pub position_iterations: u32,
    /// Complete native graph trials, including rejected backtracks.
    pub trial_evaluations: u32,
    /// Admitted conservative bound on owned logical payload [bytes].
    pub planned_live_bytes: usize,
}

/// A failed graph transaction publishes no body or partial reaction inventory.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum DistanceGraphError {
    /// Original pair geometry, native dynamics, local branch or cancellation refusal.
    Step(DistanceStepError),
    /// Original initial gap or range-rate tolerance is not satisfied.
    InitialConstraint { constraint: usize },
    /// A prescribed length path or its represented continuation is inadmissible.
    PrescribedMotion {
        constraint: usize,
        reason: &'static str,
    },
    /// The normalized mobility is dependent or below its declared pivot floor.
    DependentConstraints { constraint: usize },
    /// Caller resource cap exceeded before the corresponding work/allocation.
    Limit {
        resource: &'static str,
        required: u128,
        maximum: u128,
    },
    /// A transaction buffer could not reserve its admitted capacity.
    Allocation,
}
impl From<DistanceStepError> for DistanceGraphError {
    fn from(error: DistanceStepError) -> Self {
        Self::Step(error)
    }
}
impl From<DynamicsError> for DistanceGraphError {
    fn from(error: DynamicsError) -> Self {
        Self::Step(error.into())
    }
}
impl core::fmt::Display for DistanceGraphError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "rigid distance graph: {self:?}")
    }
}
impl std::error::Error for DistanceGraphError {}

/// Evolve a regular distance graph using simultaneous native constraint impulses.
///
/// Every position correction solves the complete initial shared-body mobility
/// and backtracks against the maximum actual length residual after native free
/// evolution. The endpoint solve uses the complete actual endpoint mobility.
/// Bodies may participate in several rods; disconnected islands and unlinked
/// bodies are admitted. Redundant or insufficiently regular rows refuse.
///
/// Trial buffers are reused. For B bodies and R rods, a native trial is O(B+R),
/// mobility assembly O(R²), and the two bounded factorizations O(R³). Cancellation
/// is checked for each body/row and before publication; no intermediate states
/// escape. Reactions are solved simultaneously, then applied at unchanged poses
/// in caller order by the original pair owner. That ordered application retains
/// all shared-body work cross terms and only affects floating-point roundoff.
///
/// This extends the original locally connected SHAKE/RATTLE construction. It
/// does not claim exact energy conservation, symplecticity, a global branch
/// certificate, unilateral cables, impacts, or moving length commands.
#[allow(clippy::too_many_arguments)]
pub fn advance_distance_graph(
    states: &[RigidBodyState],
    properties: &[MassProperties],
    links: &[DistanceLink],
    integrator: RigidBodyIntegrator,
    wrenches: &[Wrench],
    duration_s: f64,
    controls: DistanceGraphControls,
    cancelled: impl FnMut() -> bool,
) -> Result<DistanceGraphStep, DistanceGraphError> {
    advance_between(
        states,
        properties,
        Rows::Distances(links),
        integrator,
        wrenches,
        duration_s,
        controls,
        None,
        cancelled,
    )
    .map(|(step, _)| step)
}

#[allow(clippy::too_many_arguments)]
fn advance_between(
    states: &[RigidBodyState],
    properties: &[MassProperties],
    links: Rows<'_>,
    integrator: RigidBodyIntegrator,
    wrenches: &[Wrench],
    duration_s: f64,
    controls: DistanceGraphControls,
    targets: Option<&[driven::Target]>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<(DistanceGraphStep, driven::Work), DistanceGraphError> {
    poll(&mut cancelled)?;
    let local = controls.step;
    if states.len() != properties.len()
        || states.len() != wrenches.len()
        || !positive(duration_s)
        || !positive(local.maximum_impulse_ns)
        || !positive(local.maximum_rotation_increment_rad)
        || local.maximum_rotation_increment_rad > 0.25
        || local.maximum_position_iterations == 0
        || local.maximum_backtracks == 0
        || !positive(controls.position_tolerance_m)
        || !positive(controls.velocity_tolerance_mps)
        || !controls.minimum_scaled_pivot.is_finite()
        || controls.minimum_scaled_pivot < 4096.0 * f64::EPSILON
        || controls.minimum_scaled_pivot >= 1.0
    {
        return Err(DistanceStepError::Invalid("graph shapes or finite local controls").into());
    }
    limit(
        "bodies",
        states.len() as u128,
        controls.maximum_bodies as u128,
    )?;
    limit(
        "constraints",
        links.len() as u128,
        controls.maximum_constraints.min(MAX_DISTANCE_ROWS) as u128,
    )?;
    limit(
        "trial evaluations",
        1,
        u128::from(controls.maximum_trial_evaluations),
    )?;
    let bytes = planned_bytes(states.len(), links.len());
    limit("live bytes", bytes, controls.maximum_live_bytes as u128)?;
    let problem = Problem {
        states,
        properties,
        links,
        integrator,
        wrenches,
        duration_s,
        controls,
        targets,
    };
    let mut initial = reserved(links.len())?;
    for (i, link) in links.iter().enumerate() {
        poll(&mut cancelled)?;
        if link.bodies[0] == link.bodies[1] || link.bodies.iter().any(|&b| b >= states.len()) {
            return Err(
                DistanceStepError::Invalid("rod endpoints need distinct existing bodies").into(),
            );
        }
        let mut g = problem.geometry(link, link.bodies.map(|b| states[b]))?;
        if let Some(targets) = targets {
            g.rate = (driven::range_rate(
                link.bodies.map(|b| states[b]),
                link.bodies.map(|b| properties[b]),
                link.points_body_m,
                g.normal,
            )? - Dd::from_f64(targets[i].initial_rate))
            .to_f64();
        }
        if !g.rate.is_finite()
            || g.gap.abs() > controls.position_tolerance_m
            || g.rate.abs() > controls.velocity_tolerance_mps
        {
            return Err(DistanceGraphError::InitialConstraint { constraint: i });
        }
        initial.push(g);
    }
    if let Some(constraint) = links.spherical_violation(
        &initial,
        controls.position_tolerance_m,
        controls.velocity_tolerance_mps,
    ) {
        return Err(DistanceGraphError::InitialConstraint { constraint });
    }
    let mut mobility = Mobility::new(&problem, states, &initial, &mut cancelled)?;
    let mut impulse = zeros(links.len())?;
    let mut correction = zeros(links.len())?;
    let mut next = zeros(links.len())?;
    let mut accepted = Trial::new(states.len(), links.len())?;
    let mut candidate = Trial::new(states.len(), links.len())?;
    let mut evaluations = 0;
    problem.evaluate(
        &impulse,
        &initial,
        &mut accepted,
        &mut evaluations,
        &mut cancelled,
    )?;
    let mut iterations = 0;
    while links.error(&accepted.geometry) > controls.position_tolerance_m {
        if iterations == local.maximum_position_iterations {
            return Err(DistanceStepError::NotConverged.into());
        }
        for (value, g) in correction.iter_mut().zip(&accepted.geometry) {
            *value = -g.gap / duration_s;
        }
        mobility.solve(&mut correction)?;
        let mut scale = 1.0;
        let mut improved = false;
        for _ in 0..local.maximum_backtracks {
            poll(&mut cancelled)?;
            for i in 0..links.len() {
                next[i] = impulse[i] + scale * correction[i];
            }
            if next.iter().zip(&impulse).all(|(a, b)| a == b) {
                break;
            }
            if links.admits_impulses(&next, local.maximum_impulse_ns) {
                match problem.evaluate(
                    &next,
                    &initial,
                    &mut candidate,
                    &mut evaluations,
                    &mut cancelled,
                ) {
                    Ok(())
                        if links.error(&candidate.geometry) < links.error(&accepted.geometry) =>
                    {
                        core::mem::swap(&mut accepted, &mut candidate);
                        impulse.copy_from_slice(&next);
                        improved = true;
                        break;
                    }
                    Ok(()) | Err(DistanceGraphError::Step(DistanceStepError::LocalBranch)) => {}
                    Err(error) => return Err(error),
                }
            }
            scale *= 0.5;
        }
        if !improved {
            return Err(DistanceStepError::LocalBranch.into());
        }
        iterations += 1;
    }
    let mut drive_work = driven::Work::default();
    if let Some(targets) = targets {
        drive_work.add(&accepted.impulses, &impulse, targets, false, &mut cancelled)?;
    }
    mobility = Mobility::new(
        &problem,
        &accepted.states,
        &accepted.geometry,
        &mut cancelled,
    )?;
    for (value, g) in impulse.iter_mut().zip(&accepted.geometry) {
        *value = -g.rate;
    }
    mobility.solve(&mut impulse)?;
    if !links.admits_impulses(&impulse, local.maximum_impulse_ns) {
        return Err(DistanceStepError::LocalBranch.into());
    }
    problem.apply(
        &mut accepted.states,
        &impulse,
        &accepted.geometry,
        &mut candidate.impulses,
        &mut cancelled,
    )?;
    if let Some(targets) = targets {
        drive_work.add(&candidate.impulses, &impulse, targets, true, &mut cancelled)?;
    }
    let mut residuals = reserved(links.len())?;
    for (i, link) in links.iter().enumerate() {
        poll(&mut cancelled)?;
        let g = problem.endpoint_geometry(i, link.bodies.map(|b| accepted.states[b]))?;
        if g.gap.abs() > controls.position_tolerance_m
            || g.rate.abs() > controls.velocity_tolerance_mps
        {
            return Err(DistanceStepError::NotConverged.into());
        }
        if let Some(targets) = targets {
            let target = targets[i];
            let represented_gap = Dd::from_f64(g.gap) - target.length_residual;
            let represented_rate = Dd::from_f64(g.rate)
                + (target.final_rate - Dd::from_f64(target.final_rate.to_f64()));
            if represented_gap.abs().to_f64() > controls.position_tolerance_m
                || represented_rate.abs().to_f64() > controls.velocity_tolerance_mps
            {
                return Err(DistanceStepError::NotConverged.into());
            }
        }
        residuals.push(DistanceResidual {
            gap_m: g.gap,
            range_rate_mps: g.rate,
        });
        accepted.geometry[i] = g;
    }
    if links
        .spherical_violation(
            &accepted.geometry,
            controls.position_tolerance_m,
            controls.velocity_tolerance_mps,
        )
        .is_some()
    {
        return Err(DistanceStepError::NotConverged.into());
    }
    let mut work = fs_math::dd::Dd::ZERO;
    for event in accepted.impulses.iter().chain(&candidate.impulses) {
        work = work + fs_math::dd::Dd::from_f64(event.impulse_work_j);
    }
    let constraint_work_j = work.to_f64();
    if !constraint_work_j.is_finite() {
        return Err(DistanceStepError::Invalid("graph impulse work").into());
    }
    poll(&mut cancelled)?;
    Ok((
        DistanceGraphStep {
            states_after: accepted.states,
            position_impulses: accepted.impulses,
            free_steps: accepted.steps,
            velocity_impulses: candidate.impulses,
            constraint_work_j,
            residuals,
            position_iterations: iterations,
            trial_evaluations: evaluations,
            planned_live_bytes: bytes as usize,
        },
        drive_work,
    ))
}

struct Problem<'a> {
    states: &'a [RigidBodyState],
    properties: &'a [MassProperties],
    links: Rows<'a>,
    integrator: RigidBodyIntegrator,
    wrenches: &'a [Wrench],
    duration_s: f64,
    controls: DistanceGraphControls,
    targets: Option<&'a [driven::Target]>,
}
impl Problem<'_> {
    fn endpoint_geometry(
        &self,
        index: usize,
        states: [RigidBodyState; 2],
    ) -> Result<Geometry, DistanceStepError> {
        let mut link = self.links.get(index);
        if let Some(targets) = self.targets {
            link.kind = RowKind::Distance(targets[index].length);
        }
        let mut g = self.geometry(link, states)?;
        if let Some(targets) = self.targets {
            let target = targets[index];
            // Solve against the original polynomial, retaining the part lost
            // when the endpoint length/rate are rounded for continuation.
            g.gap = (Dd::from_f64(g.gap) + target.length_residual).to_f64();
            let link = self.links.get(index);
            g.rate = (driven::range_rate(
                states,
                link.bodies.map(|b| self.properties[b]),
                link.points_body_m,
                g.normal,
            )? - target.final_rate)
                .to_f64();
            if !g.gap.is_finite() || !g.rate.is_finite() {
                return Err(DistanceStepError::Invalid("prescribed graph residual"));
            }
        }
        Ok(g)
    }
    fn constraint(&self, link: Row) -> Result<RigidDistanceConstraint, DistanceStepError> {
        let RowKind::Distance(length_m) = link.kind else {
            return Err(DistanceStepError::Invalid(
                "a Cartesian closure is not a radial rod",
            ));
        };
        RigidDistanceConstraint::new(
            link.bodies.map(|b| self.properties[b]),
            link.points_body_m,
            length_m,
            self.controls.position_tolerance_m,
            self.controls.velocity_tolerance_mps,
        )
    }
    fn geometry(
        &self,
        link: Row,
        states: [RigidBodyState; 2],
    ) -> Result<Geometry, DistanceStepError> {
        match link.kind {
            RowKind::Distance(_) => self.constraint(link)?.geometry(states),
            RowKind::Cartesian(axis) => joints::geometry(
                states,
                link.bodies.map(|b| self.properties[b]),
                link.points_body_m,
                axis,
            ),
        }
    }
    fn apply(
        &self,
        states: &mut [RigidBodyState],
        impulses: &[f64],
        geometry: &[Geometry],
        events: &mut Vec<ActionReactionImpulseReceipt>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), DistanceGraphError> {
        events.clear();
        for ((link, &j), g) in self.links.iter().zip(impulses).zip(geometry) {
            poll(cancelled)?;
            let event = match link.kind {
                RowKind::Distance(_) => self
                    .constraint(link)?
                    .pair_impulse(link.bodies.map(|b| states[b]), g.normal.scale(j))?,
                RowKind::Cartesian(_) => crate::apply_equal_and_opposite_impulse_at_body_points(
                    states[link.bodies[0]],
                    self.properties[link.bodies[0]],
                    link.points_body_m[0],
                    states[link.bodies[1]],
                    self.properties[link.bodies[1]],
                    link.points_body_m[1],
                    g.normal.scale(j),
                )?,
            };
            states[link.bodies[0]] = event.first.state_after;
            states[link.bodies[1]] = event.second.state_after;
            events.push(event);
        }
        Ok(())
    }
    fn evaluate(
        &self,
        impulses: &[f64],
        initial: &[Geometry],
        out: &mut Trial,
        evaluations: &mut u32,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), DistanceGraphError> {
        poll(cancelled)?;
        limit(
            "trial evaluations",
            u128::from(*evaluations) + 1,
            u128::from(self.controls.maximum_trial_evaluations),
        )?;
        *evaluations += 1;
        out.states.clear();
        out.states.extend_from_slice(self.states);
        self.apply(
            &mut out.states,
            impulses,
            initial,
            &mut out.impulses,
            cancelled,
        )?;
        out.steps.clear();
        out.geometry.clear();
        for (i, state) in out.states.iter_mut().enumerate() {
            poll(cancelled)?;
            let inertia = self.properties[i].principal_inertia_body();
            let travel = (norm(state.angular_momentum_body())
                + self.duration_s * norm(self.wrenches[i].torque_body))
                * (self.duration_s / inertia.x.min(inertia.y).min(inertia.z));
            if !travel.is_finite() || travel > self.controls.step.maximum_rotation_increment_rad {
                return Err(DistanceStepError::LocalBranch.into());
            }
            let step = self.integrator.step(
                *state,
                self.properties[i],
                self.wrenches[i],
                self.duration_s,
            )?;
            *state = step.state_after;
            out.steps.push(step);
        }
        let chord =
            2.0 * fs_math::det::sin(0.5 * self.controls.step.maximum_rotation_increment_rad);
        for (i, (link, before)) in self.links.iter().zip(initial).enumerate() {
            poll(cancelled)?;
            let g = self.endpoint_geometry(i, link.bodies.map(|b| out.states[b]))?;
            if norm(g.normal.sub(before.normal)) > chord {
                return Err(DistanceStepError::LocalBranch.into());
            }
            out.geometry.push(g);
        }
        Ok(())
    }
}
struct Trial {
    states: Vec<RigidBodyState>,
    steps: Vec<StepReceipt>,
    impulses: Vec<ActionReactionImpulseReceipt>,
    geometry: Vec<Geometry>,
}
impl Trial {
    fn new(bodies: usize, rows: usize) -> Result<Self, DistanceGraphError> {
        Ok(Self {
            states: reserved(bodies)?,
            steps: reserved(bodies)?,
            impulses: reserved(rows)?,
            geometry: reserved(rows)?,
        })
    }
}

struct AngularRow {
    torque: [Vec3; 2],
    response: [Vec3; 2],
}
struct Mobility {
    diagonal: Vec<f64>,
    factor: Cholesky,
}
impl Mobility {
    fn new(
        problem: &Problem<'_>,
        states: &[RigidBodyState],
        geometry: &[Geometry],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, DistanceGraphError> {
        let count = problem.links.len();
        let mut angular: Vec<AngularRow> = reserved(count)?;
        for (link, g) in problem.links.iter().zip(geometry) {
            poll(cancelled)?;
            let torque = core::array::from_fn(|k| {
                link.points_body_m[k].cross(
                    states[link.bodies[k]]
                        .pose()
                        .orientation()
                        .rotate_world_to_body(g.normal),
                )
            });
            angular.push(AngularRow {
                torque,
                response: [
                    problem.properties[link.bodies[0]].angular_velocity_body_checked(torque[0])?,
                    problem.properties[link.bodies[1]].angular_velocity_body_checked(torque[1])?,
                ],
            });
        }
        let mut matrix = zeros(count * count)?;
        let mut diagonal = reserved(count)?;
        for i in 0..count {
            poll(cancelled)?;
            for j in 0..=i {
                let mut value = 0.0;
                for a in 0..2 {
                    for b in 0..2 {
                        let body = problem.links.get(i).bodies[a];
                        if body == problem.links.get(j).bodies[b] {
                            let entry = geometry[i].normal.dot(geometry[j].normal)
                                / problem.properties[body].mass()
                                + angular[i].torque[a].dot(angular[j].response[b]);
                            value += if a == b { entry } else { -entry };
                        }
                    }
                }
                if !value.is_finite() {
                    return Err(DistanceStepError::Invalid("graph mobility").into());
                }
                matrix[i * count + j] = value;
                matrix[j * count + i] = value;
            }
            if !positive(matrix[i * count + i]) {
                return Err(DistanceGraphError::DependentConstraints { constraint: i });
            }
            diagonal.push(fs_math::det::sqrt(matrix[i * count + i]));
        }
        for i in 0..count {
            for j in 0..count {
                matrix[i * count + j] = (matrix[i * count + j] / diagonal[i]) / diagonal[j];
            }
        }
        let factor =
            cholesky(&matrix, count).map_err(|error| DistanceGraphError::DependentConstraints {
                constraint: match error {
                    fs_la::factor::FactorError::NotSpd { index }
                    | fs_la::factor::FactorError::Singular { index } => index,
                },
            })?;
        for i in 0..count {
            if factor.l(i, i).powi(2) < problem.controls.minimum_scaled_pivot {
                return Err(DistanceGraphError::DependentConstraints { constraint: i });
            }
        }
        poll(cancelled)?;
        Ok(Self { diagonal, factor })
    }
    fn solve(&self, rhs: &mut [f64]) -> Result<(), DistanceGraphError> {
        for (v, d) in rhs.iter_mut().zip(&self.diagonal) {
            *v /= d;
        }
        self.factor.solve(rhs);
        for (v, d) in rhs.iter_mut().zip(&self.diagonal) {
            *v /= d;
        }
        if rhs.iter().any(|v| !v.is_finite()) {
            return Err(DistanceStepError::Invalid("graph impulse solution").into());
        }
        Ok(())
    }
}

fn planned_bytes(bodies: usize, rows: usize) -> u128 {
    use core::mem::size_of as s;
    2 * bodies as u128 * (s::<RigidBodyState>() + s::<StepReceipt>()) as u128
        + rows as u128 * (2 * s::<ActionReactionImpulseReceipt>() + 3 * s::<Geometry>()
            + s::<AngularRow>() + 5 * s::<f64>() + s::<DistanceResidual>()) as u128
        // Matrix, old/new factors and conservative dense-factor scratch. At
        // this 32-row bound the original factor uses its unblocked panel.
        + 6 * rows as u128 * rows as u128 * s::<f64>() as u128
}
fn reserved<T>(count: usize) -> Result<Vec<T>, DistanceGraphError> {
    let mut out = Vec::new();
    out.try_reserve_exact(count)
        .map_err(|_| DistanceGraphError::Allocation)?;
    Ok(out)
}
fn zeros(count: usize) -> Result<Vec<f64>, DistanceGraphError> {
    let mut out = reserved(count)?;
    out.resize(count, 0.0);
    Ok(out)
}
fn limit(resource: &'static str, required: u128, maximum: u128) -> Result<(), DistanceGraphError> {
    if required > maximum {
        Err(DistanceGraphError::Limit {
            resource,
            required,
            maximum,
        })
    } else {
        Ok(())
    }
}
