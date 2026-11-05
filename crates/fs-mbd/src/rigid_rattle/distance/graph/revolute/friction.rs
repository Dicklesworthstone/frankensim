//! Dry axial friction on a finite hinge with the complete graph recoil.
//!
//! At the supplied poses, let `J` contain the original bilateral rows and `a`
//! the hinge's relative-spin row. With `A = J M^-1 J^T` and `b = J M^-1 a^T`,
//! an axial impulse `p` requires the simultaneous reactions `-A^-1 b p`.
//! Its physical inverse inertia is `a M^-1 a^T - b^T A^-1 b`. The original
//! graph owns the geometry, factorization and native impulse applications;
//! independent pair inertias are not substituted for the constrained response.
//!
//! This is a fixed-pose impulse operation, not a time integrator. Its midpoint
//! impulse work is not a resolved finite-slip heat source. One hinge has dry
//! friction; all other supplied joints retain their original bilateral laws.

use super::super::{
    ArticulatedJoint, DistanceGraphControls, DistanceGraphError, DistanceStepError, Geometry,
    MAX_DISTANCE_ROWS, Mobility, Problem, Row, Rows, limit, planned_bytes_for, poll, reserved,
    zeros,
};
use super::{CoupleImpulseReceipt, JointReactionImpulse, NativeReaction, components, difference};
use crate::{
    DynamicsError, Gravity, MassProperties, MomentumWork, RigidBodyIntegrator, RigidBodyState,
    Vec3, Wrench, momentum_difference, rounded_vector,
};
use fs_contact::tangential::impact::{
    CoulombImpulseError, CoulombImpulseInput, FrictionLaw, FrictionRegime,
    solve_isotropic_coulomb_impulse,
};
use fs_math::dd::Dd;

/// Physical phase and caller-supplied torque-impulse capacities [N m s].
/// No coefficient, normal load, contact radius or elapsed time is inferred.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RevoluteFrictionInput {
    /// An already moving hinge uses only kinetic capacity, including arrest.
    /// The caller owns any preceding drive or smooth mechanical step.
    Sliding {
        /// Nonnegative admitted kinetic torque-impulse capacity [N m s].
        kinetic_capacity_nms: f64,
    },
    /// An exactly resting axial DOF receives an equal-and-opposite drive couple.
    /// The drive and friction act together. Within the static capacity the
    /// unchanged state has exactly zero work; beyond it only kinetic capacity
    /// opposes breakaway. A nonzero initial axial rate refuses this operation.
    RestingDrive {
        /// Signed axial drive on body two; body one receives its negative.
        drive_impulse_nms: f64,
        /// Nonnegative static bound, at least the kinetic bound [N m s].
        static_capacity_nms: f64,
        /// Nonnegative kinetic capacity after breakaway [N m s].
        kinetic_capacity_nms: f64,
    },
}

/// Original graph admission plus dimensional bounds for the angular port.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RevoluteFrictionControls {
    /// Uses the original body/row/byte, position/rate, reaction-vector and
    /// scaled-pivot bounds. Time-stepping iteration/travel fields are unused.
    pub graph: DistanceGraphControls,
    /// Bounds each actual drive, friction and resultant axial impulse [N m s].
    pub maximum_axial_impulse_nms: f64,
    /// Actual native axial rate versus the original scalar solution [rad/s].
    pub angular_rate_tolerance_rad_s: f64,
    /// Bounds both native work closure and bilateral residual work [J].
    pub work_tolerance_j: f64,
    /// Bounds the actual total linear-momentum increment [N s].
    pub linear_impulse_tolerance_ns: f64,
    /// Bounds the actual total angular-momentum increment [N m s]. This also
    /// admits or rejects any couple caused by finite joint closure errors.
    pub angular_impulse_tolerance_nms: f64,
}

/// Outcome of the selected physical phase; no speed-based capture threshold.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RevoluteFrictionRegime {
    /// An initially resting hinge held its authored drive without motion.
    Held,
    /// Kinetic capacity left nonzero axial motion, or drive caused breakaway.
    Sliding,
    /// Kinetic capacity removed the initial relative spin. The actual residual
    /// is still reported; no pose, spin or momentum is snapped to zero.
    Arrested,
}

/// Original states and the independently evaluated physical impulse account.
#[derive(Debug, Clone, PartialEq)]
#[must_use]
pub struct RevoluteFrictionReport {
    /// All original graph bodies, with poses unchanged and native momenta updated.
    pub states_after: Vec<RigidBodyState>,
    /// Static holding, continuing kinetic motion, or kinetic arrest.
    pub regime: RevoluteFrictionRegime,
    /// Normalized mean of the two admitted, positively aligned native axes.
    pub axis_world: Vec3,
    /// Second-minus-first angular velocity projected on the actual axis [rad/s].
    pub relative_rate_before_rad_s: f64,
    /// The same relative rate reconstructed from accepted native momenta [rad/s].
    pub relative_rate_after_rad_s: f64,
    /// Actual native axial rate minus the predicted constitutive endpoint [rad/s].
    pub axial_rate_residual_rad_s: f64,
    /// The full graph's Schur-complement response [1/(kg m²)].
    pub inverse_effective_inertia: f64,
    /// Authored drive, or zero for the already sliding phase [N m s].
    pub drive_impulse_nms: f64,
    /// Signed impulse on the second body about `axis_world`.
    pub friction_impulse_nms: f64,
    /// Coupled constraint coordinates [N s], in the original scalar-row order.
    pub constraint_impulses_ns: Vec<f64>,
    /// Native ordered applications retain all shared-body work cross terms.
    /// Their individual work allocation differs from simultaneous port work.
    pub native_axial_impulse: CoupleImpulseReceipt,
    /// Original point and couple reactions, in the authored scalar-row order.
    pub native_constraint_impulses: Vec<JointReactionImpulse>,
    /// Each physical port uses the original and final relative rates, including
    /// the full simultaneous constraint response. These are signed work [J].
    pub drive_work_j: f64,
    /// Nonpositive kinetic-friction impulse work, or zero during static holding [J].
    pub friction_work_j: f64,
    /// Bilateral constraint residual work, separately bounded and never heat [J].
    pub constraint_work_j: f64,
    /// Reconstructed from actual native momentum differences, without
    /// subtracting large total energies. No dissipation is relabeled as heat.
    pub kinetic_energy_change_j: f64,
    /// Native kinetic change minus drive, friction and constraint impulse work [J].
    pub work_energy_residual_j: f64,
    /// Actual total linear-momentum increment, whose internal commanded sum is zero.
    pub linear_momentum_change_ns: Vec3,
    /// Actual spin plus orbital angular-momentum increment about the world origin.
    pub angular_momentum_change_nms: Vec3,
    /// Conservative original graph bound on the transaction's owned payload [bytes].
    pub planned_live_bytes: usize,
}

/// Every refusal leaves all borrowed input states unchanged.
#[derive(Debug, Clone, PartialEq)]
pub enum RevoluteFrictionError {
    /// Original graph geometry, inertia, reaction, resource or cancellation refusal.
    Native(DistanceGraphError),
    /// Original work-conjugate Coulomb projection refused its angular inputs.
    Coulomb(CoulombImpulseError),
    /// Malformed physical input or dimensional control.
    Invalid(&'static str),
    /// Static capacity cannot be used to stop an already moving hinge.
    NotResting {
        /// Actual supplied initial relative spin, without a capture threshold.
        relative_rate_rad_s: f64,
    },
    /// The original mechanical response or native endpoint did not resolve.
    Unresolved(&'static str),
}
impl From<DistanceGraphError> for RevoluteFrictionError {
    fn from(error: DistanceGraphError) -> Self {
        Self::Native(error)
    }
}
impl From<DistanceStepError> for RevoluteFrictionError {
    fn from(error: DistanceStepError) -> Self {
        Self::Native(error.into())
    }
}
impl From<DynamicsError> for RevoluteFrictionError {
    fn from(error: DynamicsError) -> Self {
        Self::Native(error.into())
    }
}
impl core::fmt::Display for RevoluteFrictionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite revolute friction: {self:?}")
    }
}
impl std::error::Error for RevoluteFrictionError {}

/// Apply dry friction to one existing revolute joint, with every graph reaction.
///
/// `hinge_index` selects an `ArticulatedJoint::Revolute` in the original graph.
/// Positions and all bilateral velocity constraints must already be admissible.
/// The original full-inertia graph factorization admits only independent rows;
/// a locked or unresolved axial DOF refuses without numerical compliance.
/// Sliding uses the shared Coulomb projection in angular work-conjugate units.
/// The resting-drive branch separately enforces its static capacity before
/// selecting kinetic breakaway, never spending static capacity on moving spin.
///
/// Cost is the original bounded O(rows³) dense factorization, O(rows²) assembly
/// and O(bodies + rows) native application, with at most 32 bilateral rows.
/// No new graph dependency, constitutive stiffness, pose projection, hidden
/// normal load, drive trajectory, arrest time or thermal allocation is added.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn apply_revolute_friction_impulse(
    states: &[RigidBodyState],
    properties: &[MassProperties],
    joints: &[ArticulatedJoint],
    hinge_index: usize,
    input: RevoluteFrictionInput,
    controls: RevoluteFrictionControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<RevoluteFrictionReport, RevoluteFrictionError> {
    poll(&mut cancelled)?;
    validate(input, controls)?;
    if states.len() != properties.len() {
        return Err(RevoluteFrictionError::Invalid("body/property lengths"));
    }
    let Some(ArticulatedJoint::Revolute(hinge)) = joints.get(hinge_index).copied() else {
        return Err(RevoluteFrictionError::Invalid("selected revolute joint"));
    };
    let graph = controls.graph;
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
    // The original graph's larger two-trial payload is a conservative bound
    // for this one-state, two-geometry, one-factor impulse transaction.
    let bytes = planned_bytes_for::<JointReactionImpulse>(states.len(), count);
    limit("live bytes", bytes, graph.maximum_live_bytes as u128)?;
    for (&state, &property) in states.iter().zip(properties) {
        poll(&mut cancelled)?;
        state.validate()?;
        property.validate()?;
    }
    // Only the original geometry and fixed-pose impulse methods are used.
    // No free-flow step, time interval, gravity or wrench is evaluated here.
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
    let mut before: Vec<Geometry> = reserved(count)?;
    for (i, row) in rows.iter().enumerate() {
        poll(&mut cancelled)?;
        if row.bodies[0] == row.bodies[1] || row.bodies.iter().any(|&b| b >= states.len()) {
            return Err(RevoluteFrictionError::Invalid(
                "distinct existing joint bodies",
            ));
        }
        let g = problem.geometry(row, row.bodies.map(|b| states[b]), None)?;
        if !g.rate.is_finite()
            || g.gap.abs() > graph.position_tolerance_m
            || g.rate.abs() > graph.velocity_tolerance_mps
        {
            return Err(DistanceGraphError::InitialConstraint { constraint: i }.into());
        }
        before.push(g);
    }
    if let Some(constraint) = rows.grouped_violation(
        &before,
        graph.position_tolerance_m,
        graph.velocity_tolerance_mps,
    ) {
        return Err(DistanceGraphError::InitialConstraint { constraint }.into());
    }
    let bodies = hinge.bodies();
    let axes = core::array::from_fn::<_, 2, _>(|i| {
        states[bodies[i]]
            .pose()
            .orientation()
            .rotate_body_to_world(hinge.axes_body()[i])
    });
    let axis = axes[0]
        .add(axes[1])
        .normalized("finite hinge friction axis")?;
    let rate_before = relative_rate(states, properties, bodies, axis)?;
    if matches!(input, RevoluteFrictionInput::RestingDrive { .. }) && rate_before != 0.0 {
        return Err(RevoluteFrictionError::NotResting {
            relative_rate_rad_s: rate_before,
        });
    }
    let factor = Mobility::new(&problem, states, &before, &mut cancelled)?;
    // Analytic unit-couple response: do not subtract two large trial momenta
    // to obtain mobility, or probe unrelated unconstrained pair masses.
    let mut spin = [Vec3::ZERO; 2];
    for i in 0..2 {
        let q = states[bodies[i]].pose().orientation();
        spin[i] = q.rotate_body_to_world(properties[bodies[i]].angular_velocity_body_checked(
            q.rotate_world_to_body(axis.scale(if i == 0 { -1.0 } else { 1.0 })),
        )?);
    }
    let raw = project(difference(spin[1], spin[0]), axis);
    let mut coupling = [0.0; MAX_DISTANCE_ROWS];
    let mut reactions = zeros(count)?;
    for (i, (row, g)) in rows.iter().zip(&before).enumerate() {
        poll(&mut cancelled)?;
        let mut value = Dd::ZERO;
        for side in 0..2 {
            for h in 0..2 {
                if row.bodies[side] == bodies[h] {
                    let response = if row.is_angular() {
                        spin[h]
                    } else {
                        let arm = states[bodies[h]]
                            .pose()
                            .orientation()
                            .rotate_body_to_world(g.arms_body_m[side]);
                        spin[h].cross(arm)
                    };
                    let sign = if side == 0 { -1.0 } else { 1.0 };
                    value = value + dot(response, g.normal) * Dd::from_f64(sign);
                }
            }
        }
        coupling[i] = value.to_f64();
        reactions[i] = -coupling[i];
    }
    factor.solve(&mut reactions)?;
    let inverse = coupling[..count]
        .iter()
        .zip(&reactions)
        .fold(raw, |sum, (&b, &r)| sum + Dd::from_f64(b) * Dd::from_f64(r))
        .to_f64();
    if !inverse.is_finite()
        || inverse <= 0.0
        || !raw.to_f64().is_finite()
        || raw.to_f64() <= 0.0
        || inverse / raw.to_f64() < graph.minimum_scaled_pivot
    {
        return Err(RevoluteFrictionError::Unresolved(
            "locked or ill-conditioned hinge axis",
        ));
    }
    let (drive, friction, target_rate, regime) = response(input, rate_before, inverse)?;
    let total = (Dd::from_f64(drive) + Dd::from_f64(friction)).to_f64();
    if [drive, friction, total]
        .iter()
        .any(|v| !v.is_finite() || v.abs() > controls.maximum_axial_impulse_nms)
    {
        return Err(RevoluteFrictionError::Invalid("axial impulse bound"));
    }
    for r in &mut reactions {
        *r *= total;
    }
    if !rows.admits_impulses(&reactions, &before, graph.step.maximum_impulse_ns) {
        return Err(RevoluteFrictionError::Unresolved(
            "bilateral reaction bound",
        ));
    }
    let mut after_states = reserved(states.len())?;
    after_states.extend_from_slice(states);
    let native = JointReactionImpulse::couple(
        bodies.map(|b| states[b]),
        bodies.map(|b| properties[b]),
        axis.scale(total),
    )?;
    let pair = native.states_after();
    for i in 0..2 {
        after_states[bodies[i]] = pair[i];
    }
    let JointReactionImpulse::Couple(native_axial_impulse) = native else {
        return Err(RevoluteFrictionError::Invalid("native axial couple"));
    };
    let mut native_constraints = reserved(count)?;
    problem.apply(
        &mut after_states,
        &reactions,
        &before,
        &mut native_constraints,
        &mut cancelled,
    )?;
    let rate_after = relative_rate(&after_states, properties, bodies, axis)?;
    let rate_residual = (Dd::from_f64(rate_after) - Dd::from_f64(target_rate)).to_f64();
    if !rate_residual.is_finite() || rate_residual.abs() > controls.angular_rate_tolerance_rad_s {
        return Err(RevoluteFrictionError::Unresolved("native axial endpoint"));
    }
    let mean_rate = (Dd::from_f64(rate_before) + Dd::from_f64(rate_after)) * Dd::from_f64(0.5);
    let drive_work = Dd::from_f64(drive) * mean_rate;
    let friction_work = Dd::from_f64(friction) * mean_rate;
    let mut constraint_work = Dd::ZERO;
    let mut after = reserved(count)?;
    for (i, row) in rows.iter().enumerate() {
        poll(&mut cancelled)?;
        let g = problem.geometry(row, row.bodies.map(|b| after_states[b]), None)?;
        constraint_work = constraint_work
            + Dd::from_f64(reactions[i])
                * (Dd::from_f64(before[i].rate) + Dd::from_f64(g.rate))
                * Dd::from_f64(0.5);
        after.push(g);
    }
    if rows
        .grouped_violation(
            &after,
            graph.position_tolerance_m,
            graph.velocity_tolerance_mps,
        )
        .is_some()
        || after
            .iter()
            .any(|g| !g.rate.is_finite() || g.rate.abs() > graph.velocity_tolerance_mps)
    {
        return Err(RevoluteFrictionError::Unresolved(
            "native bilateral endpoint",
        ));
    }
    let (energy, linear, angular) = account(states, &after_states, properties, &mut cancelled)?;
    let work_residual = (energy - drive_work - friction_work - constraint_work).to_f64();
    if [
        energy.to_f64(),
        drive_work.to_f64(),
        friction_work.to_f64(),
        constraint_work.to_f64(),
        work_residual,
    ]
    .iter()
    .any(|x| !x.is_finite())
        || !linear.is_finite()
        || !angular.is_finite()
    {
        return Err(RevoluteFrictionError::Unresolved(
            "nonfinite impulse account",
        ));
    }
    if friction_work.to_f64() > 0.0
        || work_residual.abs() > controls.work_tolerance_j
        || constraint_work.to_f64().abs() > controls.work_tolerance_j
    {
        return Err(RevoluteFrictionError::Unresolved(
            "native work and passivity",
        ));
    }
    if norm(linear) > controls.linear_impulse_tolerance_ns
        || norm(angular) > controls.angular_impulse_tolerance_nms
    {
        return Err(RevoluteFrictionError::Unresolved("native momentum balance"));
    }
    poll(&mut cancelled)?;
    Ok(RevoluteFrictionReport {
        states_after: after_states,
        regime,
        axis_world: axis,
        relative_rate_before_rad_s: rate_before,
        relative_rate_after_rad_s: rate_after,
        axial_rate_residual_rad_s: rate_residual,
        inverse_effective_inertia: inverse,
        drive_impulse_nms: drive,
        friction_impulse_nms: friction,
        constraint_impulses_ns: reactions,
        native_axial_impulse,
        native_constraint_impulses: native_constraints,
        drive_work_j: drive_work.to_f64(),
        friction_work_j: friction_work.to_f64(),
        constraint_work_j: constraint_work.to_f64(),
        kinetic_energy_change_j: energy.to_f64(),
        work_energy_residual_j: work_residual,
        linear_momentum_change_ns: linear,
        angular_momentum_change_nms: angular,
        planned_live_bytes: bytes as usize,
    })
}

fn validate(
    input: RevoluteFrictionInput,
    controls: RevoluteFrictionControls,
) -> Result<(), RevoluteFrictionError> {
    let graph = controls.graph;
    if [
        graph.position_tolerance_m,
        graph.velocity_tolerance_mps,
        graph.step.maximum_impulse_ns,
        controls.maximum_axial_impulse_nms,
        controls.angular_rate_tolerance_rad_s,
        controls.work_tolerance_j,
        controls.linear_impulse_tolerance_ns,
        controls.angular_impulse_tolerance_nms,
    ]
    .iter()
    .any(|x| !x.is_finite() || *x <= 0.0)
        || !graph.minimum_scaled_pivot.is_finite()
        || !(4096.0 * f64::EPSILON..1.0).contains(&graph.minimum_scaled_pivot)
    {
        return Err(RevoluteFrictionError::Invalid(
            "finite positive impulse controls",
        ));
    }
    let kinetic = match input {
        RevoluteFrictionInput::Sliding {
            kinetic_capacity_nms,
        } => kinetic_capacity_nms,
        RevoluteFrictionInput::RestingDrive {
            drive_impulse_nms,
            static_capacity_nms,
            kinetic_capacity_nms,
        } => {
            if !drive_impulse_nms.is_finite()
                || !static_capacity_nms.is_finite()
                || static_capacity_nms < kinetic_capacity_nms
            {
                return Err(RevoluteFrictionError::Invalid(
                    "resting drive and static capacity",
                ));
            }
            kinetic_capacity_nms
        }
    };
    if !kinetic.is_finite() || kinetic < 0.0 {
        return Err(RevoluteFrictionError::Invalid("kinetic impulse capacity"));
    }
    Ok(())
}

fn response(
    input: RevoluteFrictionInput,
    rate: f64,
    inverse: f64,
) -> Result<(f64, f64, f64, RevoluteFrictionRegime), RevoluteFrictionError> {
    match input {
        RevoluteFrictionInput::Sliding {
            kinetic_capacity_nms,
        } => {
            // The shared projection accepts consistent work-conjugate units:
            // angular impulse, angular rate and inverse inertia. Unit radius
            // coefficient merely passes the authored capacity, not a material.
            let law = FrictionLaw::Coulomb {
                static_mu: 1.0,
                kinetic_mu: 1.0,
            };
            let out = solve_isotropic_coulomb_impulse(
                CoulombImpulseInput {
                    normal_impulse: kinetic_capacity_nms,
                    relative_tangent_velocity: [rate, 0.0],
                    inverse_tangent_mass: inverse,
                },
                &law,
            )
            .map_err(RevoluteFrictionError::Coulomb)?;
            Ok((
                0.0,
                out.tangent_impulse[0],
                out.relative_tangent_velocity_after[0],
                if out.regime == FrictionRegime::Sticking {
                    RevoluteFrictionRegime::Arrested
                } else {
                    RevoluteFrictionRegime::Sliding
                },
            ))
        }
        RevoluteFrictionInput::RestingDrive {
            drive_impulse_nms: drive,
            static_capacity_nms,
            kinetic_capacity_nms,
        } => {
            if drive.abs() <= static_capacity_nms {
                Ok((drive, -drive, 0.0, RevoluteFrictionRegime::Held))
            } else {
                let friction = -drive.signum() * kinetic_capacity_nms;
                let rate = inverse * (drive + friction);
                if !rate.is_finite() {
                    return Err(RevoluteFrictionError::Unresolved("breakaway rate"));
                }
                Ok((drive, friction, rate, RevoluteFrictionRegime::Sliding))
            }
        }
    }
}

fn dot(a: Vec3, b: Vec3) -> Dd {
    project(components(a).map(Dd::from_f64), b)
}
fn project(a: [Dd; 3], b: Vec3) -> Dd {
    a.into_iter()
        .zip(components(b))
        .fold(Dd::ZERO, |sum, (a, b)| sum + a * Dd::from_f64(b))
}
fn norm(v: Vec3) -> f64 {
    v.x.hypot(v.y).hypot(v.z)
}
fn relative_rate(
    states: &[RigidBodyState],
    properties: &[MassProperties],
    bodies: [usize; 2],
    axis: Vec3,
) -> Result<f64, RevoluteFrictionError> {
    let mut spin = [Vec3::ZERO; 2];
    for i in 0..2 {
        spin[i] = states[bodies[i]]
            .point_kinematics(properties[bodies[i]], Vec3::ZERO)?
            .angular_velocity_world;
    }
    let value = project(difference(spin[1], spin[0]), axis).to_f64();
    if !value.is_finite() {
        return Err(RevoluteFrictionError::Invalid("relative axial rate"));
    }
    Ok(value)
}

// Shared fixed-pose native momentum account for the single hinge and the
// coupled revolute/prismatic friction owner. Keep the measured native momentum
// differences and polarized work in one implementation.
pub(in crate::rigid_rattle::distance::graph) fn account(
    before: &[RigidBodyState],
    after: &[RigidBodyState],
    properties: &[MassProperties],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Dd, Vec3, Vec3), DistanceGraphError> {
    let mut work = MomentumWork::default();
    let mut linear = [Dd::ZERO; 3];
    let mut angular = [Dd::ZERO; 3];
    for ((&old, &new), &property) in before.iter().zip(after).zip(properties) {
        poll(cancelled)?;
        let dp = momentum_difference(old.linear_momentum_world(), new.linear_momentum_world());
        let dl = momentum_difference(old.angular_momentum_body(), new.angular_momentum_body());
        work.add(
            old.center_of_mass_velocity_world(property)?,
            new.center_of_mass_velocity_world(property)?,
            dp,
        );
        work.add(
            property.angular_velocity_body_checked(old.angular_momentum_body())?,
            property.angular_velocity_body_checked(new.angular_momentum_body())?,
            dl,
        );
        let r = components(old.pose().position_world()).map(Dd::from_f64);
        for d in 0..3 {
            linear[d] = linear[d] + dp[d];
            let b = (d + 1) % 3;
            let c = (d + 2) % 3;
            angular[d] = angular[d] + r[b] * dp[c] - r[c] * dp[b];
        }
        for part in [
            Vec3::new(dl[0].hi, dl[1].hi, dl[2].hi),
            Vec3::new(dl[0].lo, dl[1].lo, dl[2].lo),
        ] {
            let rotated = old.pose().orientation().rotate_body_to_world(part);
            for (sum, value) in angular.iter_mut().zip(components(rotated)) {
                *sum = *sum + Dd::from_f64(value);
            }
        }
    }
    Ok((
        work.total(),
        rounded_vector(linear),
        rounded_vector(angular),
    ))
}
