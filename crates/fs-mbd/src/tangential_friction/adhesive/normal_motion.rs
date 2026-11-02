//! Native normal adhesion and two-dimensional junction sliding on a DMT pair.
//!
//! The original source-owned DMT load line supplies exact normal path work and
//! the analytic mean REAL area, including gap/contact crossings. All nine
//! entries of the actual point-force mobility enter the same midpoint solve.
//! Attraction across a finite gap changes native momentum and normal storage;
//! it never creates a solid shear area or friction heat.
//!
//! This is a fixed-Jacobian momentum phase. The binding's collinear geometry
//! markers and immutable approach datum authenticate the initial normal state.
//! Paired tractions act at their common midpoint reaction plane. That plane is
//! the declared contact discretization, not a solved partition of the two
//! materials' surface deformation. Actual marker and COM midpoint travel are
//! returned for the surrounding pose owner. Poses are not advanced here, and an
//! advanced normal state cannot be reused against unchanged marker geometry.
//!
//! The discrete source path is affine; it does not resolve curved trajectories,
//! finite rotations, impacts or contact detection. Torsion during a path with
//! real contact, and unresolved sliding reversals, explicitly refuse. Gap-only
//! spin is free. No torque, artificial stopping impulse or numerical heat is
//! supplied for an unrepresented physical coordinate.

use crate::{DynamicsError, MassProperties, PointKinematics, RigidBodyState, Vec3};
use fs_contact::interface_binding::adhesive_shear::SourceAdhesiveShear;
use fs_contact::normal_patch::dmt::{
    SourceDmtSphereContact,
    coupled::{
        DmtCoupledConfig, DmtCoupledError, DmtCoupledLoadLine, DmtCoupledRegime, DmtCoupledState,
        DmtCoupledStep,
    },
};
use fs_exec::CancelGate;
use fs_tribo::partial_slip::TangentFrame;

use super::super::hertz::{
    RigidHertzError, apply_pair_wrench, points, relative_point, relative_spin,
};
use super::super::{RigidContactBinding, RigidFrictionError, components, magnitude, poll, vector};

/// Original source ordering, geometry markers and held contact discretization.
#[derive(Clone, Copy, Debug)]
pub struct RigidDmtContact<'a> {
    /// Body indices and COM-relative geometry-marker arms, not reaction arms.
    pub binding: RigidContactBinding,
    /// Normal points from A toward B. Positive A-minus-B motion closes the gap.
    pub frame: TangentFrame,
    /// Approach when the two geometry markers coincide [m]. This fixed datum
    /// belongs to the declared geometry; never update it to hide a stale pose.
    pub reference_approach_m: f64,
    /// Original source-owned finite-gap normal contact.
    pub normal: &'a SourceDmtSphereContact,
    /// Original ordered junction-strength source, matched by the load line.
    pub shear: &'a SourceAdhesiveShear,
}

/// Original constitutive bounds and independent dimensional native-body gates.
#[derive(Clone, Copy, Debug)]
pub struct RigidDmtControls {
    /// Existing joint normal/shear solve, source domains and work-path budget.
    pub contact: DmtCoupledConfig,
    /// Maximum size of the original body array.
    pub maximum_bodies: usize,
    /// Marker collinearity, reaction coincidence and initial approach gate [m].
    pub geometry_tolerance_m: f64,
    /// Actual midpoint motion minus the original solved coordinate path [m].
    pub displacement_tolerance_m: f64,
    /// Native recoil/normal-Jacobian gate and unresolved static-motion gate.
    pub velocity_tolerance_m_s: f64,
    /// Unsupported relative normal-axis spin gate on real-contact paths.
    pub angular_velocity_tolerance_rad_s: f64,
    /// Individual and complete pair linear-impulse residual gate [N s].
    pub linear_momentum_tolerance_ns: f64,
    /// Individual and complete pair angular-impulse residual gate [N m s].
    pub angular_momentum_tolerance_nms: f64,
    /// Native work, total storage/heat and assisting endpoint work gate [J].
    pub energy_tolerance_j: f64,
}

/// Independent physical checks on actual body endpoints and the original law.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidDmtAudit {
    /// Actual displacement minus the source trial, in normal/shear order [m].
    pub displacement_residual_m: [f64; 3],
    /// Actual endpoint rates minus the full original mobility response [m/s].
    pub velocity_residual_m_s: [f64; 3],
    /// Marker approach travel minus work-conjugate reaction-point travel [m].
    pub marker_normal_displacement_residual_m: f64,
    /// Largest initial/final marker-versus-reaction normal-rate difference.
    pub marker_normal_velocity_residual_m_s: f64,
    /// Actual native kinetic change minus actual contact-impulse work [J].
    pub work_energy_residual_j: f64,
    /// Delta K + original normal potential change + physical shear heat [J].
    pub total_energy_residual_j: f64,
    /// Largest positive shear-force work at either physical endpoint [J].
    pub maximum_assisting_work_j: f64,
    /// Actual minus commanded linear impulse, in original A/B order.
    pub body_linear_momentum_residual_ns: [Vec3; 2],
    /// Actual minus commanded body-frame angular impulse, in A/B order.
    pub body_angular_momentum_residual_nms: [Vec3; 2],
    /// Complete pair linear-momentum change [N s].
    pub linear_momentum_change_ns: Vec3,
    /// Complete orbital plus spin angular-momentum change [N m s].
    pub angular_momentum_change_nms: Vec3,
}

/// Native normal/shear motion and independently reevaluated source path.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidDmtReport {
    /// Original law reevaluated on actual native midpoint coordinates.
    pub contact: DmtCoupledStep,
    /// Actual `J M^-1 J^T` for [approach, shear 1, shear 2], in kg^-1.
    pub mobility_per_kg: [[f64; 3]; 3],
    /// Constant efforts used for the actual impulse, positive resisting motion.
    pub applied_mean_force_n: [f64; 3],
    /// Paired point impulse applied to body A; body B receives its opposite.
    pub force_impulse_on_a_world_ns: Vec3,
    /// Actual initial/final relative reaction-point rates in normal/shear order.
    pub relative_velocity_m_s: [[f64; 3]; 2],
    /// Actual initial/final relative spin about the original normal.
    pub relative_spin_rad_s: [f64; 2],
    /// Actual midpoint coordinate increments, in normal/shear order [m].
    pub relative_displacement_m: [f64; 3],
    /// COM travel for the surrounding pose owner, in original A/B order [m].
    pub body_displacement_world_m: [Vec3; 2],
    /// Fixed-Jacobian geometry-marker travel for that owner, in A/B order [m].
    pub marker_displacement_world_m: [Vec3; 2],
    /// Fixed-Jacobian reaction-point travel, in original A/B order [m].
    pub reaction_displacement_world_m: [Vec3; 2],
    /// Differential native kinetic change, without subtracting large totals.
    pub kinetic_energy_change_j: f64,
    /// Signed contact impulse work on actual mean endpoint velocities [J].
    pub midpoint_impulse_work_j: f64,
    /// Original normal-potential change on the actual approach path [J].
    pub normal_storage_change_j: f64,
    /// Physical source-owned junction loss, unpartitioned [J].
    pub dissipated_energy_j: f64,
    /// Constitutive nonlinear iterations, excluding the independent path replay.
    pub iterations: usize,
    /// Trial constitutive probes, excluding the one independent path replay.
    pub evaluations: usize,
    /// All independent native acceptance residuals.
    pub audit: RigidDmtAudit,
}

/// Immutable proposal; source history and both body momenta publish together.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidDmtCandidate {
    /// Original A/B endpoints. Poses remain unchanged during this phase.
    pub bodies_after: [RigidBodyState; 2],
    /// Authenticated original source state with actual endpoint occupancy.
    pub state_after: DmtCoupledState,
    /// Physical source path, native motion and balance checks.
    pub report: RigidDmtReport,
}

/// Source refusal, unrepresented motion and native balance remain distinct.
#[derive(Debug)]
pub enum RigidDmtError {
    /// Invalid physical controls, geometry or nonfinite derived quantity.
    Invalid(&'static str),
    /// Original normal/shear constitutive owner refused the path.
    Contact(DmtCoupledError),
    /// Actual native kinematics or impulse owner refused.
    Mechanics(RigidHertzError),
    /// Original normal state disagrees with the binding's actual marker gap.
    GeometryMismatch {
        /// Original source approach minus the geometry-authenticated approach.
        residual_m: f64,
    },
    /// The source law has no torsional traction on a real solid contact.
    UnsupportedTorsion {
        /// Actual initial relative spin about the original normal.
        before_rad_s: f64,
        /// Actual final relative spin about the original normal.
        after_rad_s: f64,
    },
    /// A held path conceals continued slip or would assist an endpoint motion.
    UnresolvedSlipEvent {
        /// Largest actual tangential speed at the two endpoints.
        maximum_slip_speed_m_s: f64,
        /// Largest positive endpoint shear-force work over the requested time.
        maximum_assisting_work_j: f64,
    },
    /// Actual endpoints did not satisfy the original dimensional balance gates.
    Audit(Box<RigidDmtAudit>),
}
impl From<DmtCoupledError> for RigidDmtError {
    fn from(value: DmtCoupledError) -> Self {
        Self::Contact(value)
    }
}
impl From<RigidHertzError> for RigidDmtError {
    fn from(value: RigidHertzError) -> Self {
        Self::Mechanics(value)
    }
}
impl From<DynamicsError> for RigidDmtError {
    fn from(value: DynamicsError) -> Self {
        Self::Mechanics(value.into())
    }
}
impl From<RigidFrictionError> for RigidDmtError {
    fn from(value: RigidFrictionError) -> Self {
        Self::Mechanics(value.into())
    }
}
impl core::fmt::Display for RigidDmtError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native DMT normal/shear phase: {self:?}")
    }
}
impl std::error::Error for RigidDmtError {}

/// Captured physical pair and reusable original constitutive workspace.
/// Every trial starts from the same native endpoints; retry does not evolve
/// the captured bodies, the source properties or an accepted history token.
pub struct PreparedRigidDmtPair<'a> {
    contact: RigidDmtContact<'a>,
    before: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    reaction_arms: [Vec3; 2],
    axes: [Vec3; 3],
    markers_before: [PointKinematics; 2],
    reaction_before: [PointKinematics; 2],
    approach_m: f64,
    mobility: [[f64; 3]; 3],
    load_line: DmtCoupledLoadLine<'a>,
    dt_s: f64,
    controls: RigidDmtControls,
}

impl<'a> PreparedRigidDmtPair<'a> {
    /// Capture actual geometry and construct the original full load-line solve.
    pub fn new(
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        contact: RigidDmtContact<'a>,
        dt_s: f64,
        controls: RigidDmtControls,
        gate: &CancelGate,
    ) -> Result<Self, RigidDmtError> {
        poll(gate)?;
        let binding = contact.binding;
        if bodies.len() != properties.len()
            || bodies.len() > controls.maximum_bodies
            || binding.body_a >= bodies.len()
            || binding.body_b >= bodies.len()
            || binding.body_a == binding.body_b
            || !dt_s.is_finite()
            || dt_s <= 0.0
            || dt_s > controls.contact.maximum_step_s
            || !contact.reference_approach_m.is_finite()
            || [
                controls.geometry_tolerance_m,
                controls.displacement_tolerance_m,
                controls.velocity_tolerance_m_s,
                controls.angular_velocity_tolerance_rad_s,
                controls.linear_momentum_tolerance_ns,
                controls.angular_momentum_tolerance_nms,
                controls.energy_tolerance_j,
            ]
            .into_iter()
            .any(|v| !v.is_finite() || v < 0.0)
        {
            return Err(RigidDmtError::Invalid(
                "bounded distinct bodies, finite controls and positive time required",
            ));
        }
        let before = [bodies[binding.body_a], bodies[binding.body_b]];
        let properties = [properties[binding.body_a], properties[binding.body_b]];
        let marker_arms = [binding.arm_a_body, binding.arm_b_body];
        let axes = [
            contact.frame.normal(),
            contact.frame.longitudinal(),
            contact.frame.lateral(),
        ]
        .map(vector);
        let markers_before = points(before, properties, marker_arms)?;
        let separation = markers_before[1]
            .point_world
            .sub(markers_before[0].point_world);
        let normal_separation = finite(separation.dot(axes[0]))?;
        if magnitude(separation.sub(axes[0].scale(normal_separation)))
            > controls.geometry_tolerance_m
        {
            return Err(RigidDmtError::Invalid(
                "geometry markers must lie on the declared normal line",
            ));
        }
        let approach_m = finite(contact.reference_approach_m - normal_separation)?;
        // The original marker normal row is unchanged by shifting along the
        // normal. The common reaction plane makes every traction work-dual and
        // avoids a spurious internal couple when a gap becomes a solid contact.
        let reaction_arms = core::array::from_fn(|i| {
            let sign = if i == 0 { 0.5 } else { -0.5 };
            marker_arms[i].add(
                before[i]
                    .pose()
                    .orientation()
                    .rotate_world_to_body(separation.scale(sign)),
            )
        });
        let reaction_before = points(before, properties, reaction_arms)?;
        if magnitude(
            reaction_before[0]
                .point_world
                .sub(reaction_before[1].point_world),
        ) > controls.geometry_tolerance_m
        {
            return Err(RigidDmtError::Invalid(
                "unrepresentable common reaction plane",
            ));
        }
        let mobility = point_response(before, properties, reaction_arms, axes)?;
        let half_dt_squared = finite(0.5 * dt_s * dt_s)?;
        if half_dt_squared == 0.0 {
            return Err(RigidDmtError::Invalid("unrepresentable midpoint duration"));
        }
        let compliance = mobility.map(|row| row.map(|v| v * half_dt_squared));
        let load_line =
            DmtCoupledLoadLine::new(contact.normal, contact.shear, compliance, controls.contact)?;
        poll(gate)?;
        Ok(Self {
            contact,
            before,
            properties,
            reaction_arms,
            axes,
            markers_before,
            reaction_before,
            approach_m,
            mobility,
            load_line,
            dt_s,
            controls,
        })
    }

    /// Original sources, marker binding and fixed approach datum.
    pub const fn contact(&self) -> RigidDmtContact<'a> {
        self.contact
    }

    /// Actual signed approach authenticated by the captured geometry [m].
    pub const fn approach_m(&self) -> f64 {
        self.approach_m
    }

    /// Complete native normal/shear response, without a caller effective mass.
    pub const fn mobility_per_kg(&self) -> [[f64; 3]; 3] {
        self.mobility
    }

    /// Original midpoint load line, `dt^2 J M^-1 J^T / 2` [m/N].
    pub const fn compliance_m_per_n(&self) -> [[f64; 3]; 3] {
        self.load_line.compliance_m_per_n()
    }

    /// Actual COM-relative reaction arms, distinct from geometry-marker arms.
    pub const fn reaction_arms_body(&self) -> [Vec3; 2] {
        self.reaction_arms
    }

    /// Initialize source history from the ACTUAL normal geometry and native
    /// slip. Only the arbitrary shear-coordinate origins are caller supplied.
    pub fn initialize(&self, shear_position_m: [f64; 2]) -> Result<DmtCoupledState, RigidDmtError> {
        let velocity = relative_point(self.reaction_before);
        Ok(self.load_line.initialize_with_velocity(
            [self.approach_m, shear_position_m[0], shear_position_m[1]],
            [velocity.dot(self.axes[1]), velocity.dot(self.axes[2])],
        )?)
    }

    /// Resolve one source-owned path and independently replay it on actual body
    /// motion. The enclosing owner must apply the returned geometry travel and
    /// refresh the Jacobian before continuing from the advanced normal state.
    #[allow(clippy::too_many_lines)]
    pub fn advance(
        &mut self,
        previous: DmtCoupledState,
        gate: &CancelGate,
    ) -> Result<RigidDmtCandidate, RigidDmtError> {
        poll(gate)?;
        let axes = self.axes;
        let dt = self.dt_s;
        let controls = self.controls;
        let initial = relative_point(self.reaction_before);
        let rates_before = axes.map(|axis| initial.dot(axis));
        let old = previous.position_m();
        let geometry_residual = finite(old[0] - self.approach_m)?;
        if geometry_residual.abs() > controls.geometry_tolerance_m {
            return Err(RigidDmtError::GeometryMismatch {
                residual_m: geometry_residual,
            });
        }
        // A dynamic endpoint, not a zero previous increment, selects occupancy.
        let previous = self
            .load_line
            .with_endpoint_velocity(previous, [rates_before[1], rates_before[2]])?;
        let free = core::array::from_fn(|i| old[i] + dt * rates_before[i]);
        let trial = self.load_line.trial(previous, free, dt, gate)?;
        let applied = trial.mean_force_n;
        let impulse = axes[0]
            .scale(-dt * applied[0])
            .add(axes[1].scale(-dt * applied[1]))
            .add(axes[2].scale(-dt * applied[2]));
        let after = apply_pair_wrench(
            self.before,
            self.properties,
            self.reaction_arms,
            impulse,
            Vec3::ZERO,
            gate,
        )?;
        let reaction_after = points(after, self.properties, self.reaction_arms)?;
        let marker_arms = [
            self.contact.binding.arm_a_body,
            self.contact.binding.arm_b_body,
        ];
        let markers_after = points(after, self.properties, marker_arms)?;
        let endpoint = relative_point(reaction_after);
        let rates_after = axes.map(|axis| endpoint.dot(axis));
        let delta = core::array::from_fn(|i| dt * (0.5 * rates_before[i] + 0.5 * rates_after[i]));
        let position = core::array::from_fn(|i| old[i] + delta[i]);
        let actual =
            self.load_line
                .advance_actual(previous, position, applied, trial.regime, dt, gate)?;
        let spin = [
            relative_spin(self.reaction_before).dot(axes[0]),
            relative_spin(reaction_after).dot(axes[0]),
        ];
        if spin.into_iter().any(|v| !v.is_finite()) {
            return Err(RigidDmtError::Invalid("unrepresentable normal-axis spin"));
        }
        if actual.path.mean_real_area_m2 > 0.0
            && spin
                .into_iter()
                .any(|v| v.abs() > controls.angular_velocity_tolerance_rad_s)
        {
            return Err(RigidDmtError::UnsupportedTorsion {
                before_rad_s: spin[0],
                after_rad_s: spin[1],
            });
        }
        let maximum_slip = rates_before[1]
            .hypot(rates_before[2])
            .max(rates_after[1].hypot(rates_after[2]));
        let assisting = finite(
            dt * [rates_before, rates_after]
                .into_iter()
                .map(|v| -applied[1] * v[1] - applied[2] * v[2])
                .fold(0.0, f64::max),
        )?;
        if (actual.path.mean_real_area_m2 > 0.0
            && trial.regime == DmtCoupledRegime::Sticking
            && maximum_slip > controls.velocity_tolerance_m_s)
            || assisting > controls.energy_tolerance_j
        {
            return Err(RigidDmtError::UnresolvedSlipEvent {
                maximum_slip_speed_m_s: maximum_slip,
                maximum_assisting_work_j: assisting,
            });
        }
        let state_after = self
            .load_line
            .with_endpoint_velocity(actual.state(), [rates_after[1], rates_after[2]])?;
        let mut kinetic = 0.0;
        let mut linear = Vec3::ZERO;
        let mut angular = Vec3::ZERO;
        let mut body_linear = [Vec3::ZERO; 2];
        let mut body_angular = [Vec3::ZERO; 2];
        let mut body_travel = [Vec3::ZERO; 2];
        let mut marker_travel = [Vec3::ZERO; 2];
        let mut reaction_travel = [Vec3::ZERO; 2];
        for (i, sign) in [(0, 1.0), (1, -1.0)] {
            poll(gate)?;
            let dp = after[i]
                .linear_momentum_world()
                .sub(self.before[i].linear_momentum_world());
            let dl = after[i]
                .angular_momentum_body()
                .sub(self.before[i].angular_momentum_body());
            let mean_v = mean(
                self.reaction_before[i].center_of_mass_velocity_world,
                reaction_after[i].center_of_mass_velocity_world,
            );
            let mean_w = mean(
                self.reaction_before[i].angular_velocity_body,
                reaction_after[i].angular_velocity_body,
            );
            kinetic += dp.dot(mean_v) + dl.dot(mean_w);
            body_travel[i] = mean_v.scale(dt);
            marker_travel[i] = mean(
                self.markers_before[i].point_velocity_world,
                markers_after[i].point_velocity_world,
            )
            .scale(dt);
            reaction_travel[i] = mean(
                self.reaction_before[i].point_velocity_world,
                reaction_after[i].point_velocity_world,
            )
            .scale(dt);
            let orientation = self.before[i].pose().orientation();
            let expected = impulse.scale(sign);
            body_linear[i] = dp.sub(expected);
            body_angular[i] =
                dl.sub(self.reaction_arms[i].cross(orientation.rotate_world_to_body(expected)));
            linear = linear.add(dp);
            angular = angular
                .add(self.before[i].pose().position_world().cross(dp))
                .add(orientation.rotate_body_to_world(dl));
        }
        let work = finite(impulse.dot(mean(initial, endpoint)))?;
        let marker_normal = marker_travel[0].sub(marker_travel[1]).dot(axes[0]);
        let marker_rate_residual = (relative_point(self.markers_before).dot(axes[0])
            - rates_before[0])
            .abs()
            .max((relative_point(markers_after).dot(axes[0]) - rates_after[0]).abs());
        let audit = RigidDmtAudit {
            displacement_residual_m: core::array::from_fn(|i| {
                delta[i] - (trial.state().position_m()[i] - old[i])
            }),
            velocity_residual_m_s: core::array::from_fn(|i| {
                rates_after[i] - rates_before[i]
                    + dt * (0..3)
                        .map(|j| self.mobility[i][j] * applied[j])
                        .sum::<f64>()
            }),
            marker_normal_displacement_residual_m: marker_normal - delta[0],
            marker_normal_velocity_residual_m_s: marker_rate_residual,
            work_energy_residual_j: kinetic - work,
            total_energy_residual_j: kinetic
                + actual.path.normal_work_j
                + actual.friction_dissipation_j,
            maximum_assisting_work_j: assisting,
            body_linear_momentum_residual_ns: body_linear,
            body_angular_momentum_residual_nms: body_angular,
            linear_momentum_change_ns: linear,
            angular_momentum_change_nms: angular,
        };
        if audit
            .displacement_residual_m
            .into_iter()
            .chain([audit.marker_normal_displacement_residual_m])
            .any(|v| !v.is_finite() || v.abs() > controls.displacement_tolerance_m)
            || audit
                .velocity_residual_m_s
                .into_iter()
                .chain([marker_rate_residual])
                .any(|v| !v.is_finite() || v.abs() > controls.velocity_tolerance_m_s)
            || [audit.work_energy_residual_j, audit.total_energy_residual_j]
                .into_iter()
                .any(|v| !v.is_finite() || v.abs() > controls.energy_tolerance_j)
            || body_linear
                .into_iter()
                .chain([linear])
                .any(|v| !v.is_finite() || magnitude(v) > controls.linear_momentum_tolerance_ns)
            || body_angular
                .into_iter()
                .chain([angular])
                .any(|v| !v.is_finite() || magnitude(v) > controls.angular_momentum_tolerance_nms)
            || body_travel
                .into_iter()
                .chain(marker_travel)
                .chain(reaction_travel)
                .any(|v| !v.is_finite())
        {
            return Err(RigidDmtError::Audit(Box::new(audit)));
        }
        poll(gate)?;
        Ok(RigidDmtCandidate {
            bodies_after: after,
            state_after,
            report: RigidDmtReport {
                contact: actual,
                mobility_per_kg: self.mobility,
                applied_mean_force_n: applied,
                force_impulse_on_a_world_ns: impulse,
                relative_velocity_m_s: [rates_before, rates_after],
                relative_spin_rad_s: spin,
                relative_displacement_m: delta,
                body_displacement_world_m: body_travel,
                marker_displacement_world_m: marker_travel,
                reaction_displacement_world_m: reaction_travel,
                kinetic_energy_change_j: kinetic,
                midpoint_impulse_work_j: work,
                normal_storage_change_j: actual.path.normal_work_j,
                dissipated_energy_j: actual.friction_dissipation_j,
                iterations: trial.iterations,
                evaluations: trial.evaluations,
                audit,
            },
        })
    }
}

/// Stage a source-preserving phase without changing accepted bodies or history.
#[allow(clippy::too_many_arguments)]
pub fn prepare_rigid_dmt_contact(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contact: RigidDmtContact<'_>,
    previous: DmtCoupledState,
    dt_s: f64,
    controls: RigidDmtControls,
    gate: &CancelGate,
) -> Result<RigidDmtCandidate, RigidDmtError> {
    PreparedRigidDmtPair::new(bodies, properties, contact, dt_s, controls, gate)?
        .advance(previous, gate)
}

/// Publish both original body momenta and source history after final admission.
/// Poses remain unchanged: the surrounding pose owner must consume the returned
/// marker/COM travel before a continuation from the advanced normal state.
#[allow(clippy::too_many_arguments)]
pub fn apply_rigid_dmt_contact(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    contact: RigidDmtContact<'_>,
    history: &mut DmtCoupledState,
    dt_s: f64,
    controls: RigidDmtControls,
    gate: &CancelGate,
) -> Result<RigidDmtReport, RigidDmtError> {
    let candidate =
        prepare_rigid_dmt_contact(bodies, properties, contact, *history, dt_s, controls, gate)?;
    poll(gate)?;
    bodies[contact.binding.body_a] = candidate.bodies_after[0];
    bodies[contact.binding.body_b] = candidate.bodies_after[1];
    *history = candidate.state_after;
    Ok(candidate.report)
}

/// All three rows are POINT forces; the sliding/spin helper's third row is a
/// free couple and is deliberately not used for normal motion.
fn point_response(
    bodies: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    arms: [Vec3; 2],
    axes: [Vec3; 3],
) -> Result<[[f64; 3]; 3], RigidDmtError> {
    let mut factor = [[0.0; 12]; 3];
    for (body, sign) in [(0, 1.0), (1, -1.0)] {
        let inverse_mass_root = properties[body].mass().sqrt().recip();
        let orientation = bodies[body].pose().orientation();
        for i in 0..3 {
            let translation = components(axes[i].scale(sign * inverse_mass_root));
            let torque = arms[body]
                .cross(orientation.rotate_world_to_body(axes[i]))
                .scale(sign);
            let rotation = components(properties[body].angular_mobility_factor(torque));
            for j in 0..3 {
                factor[i][6 * body + j] = translation[j];
                factor[i][6 * body + 3 + j] = rotation[j];
            }
        }
    }
    let mut response = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..=i {
            response[i][j] = finite(factor[i].iter().zip(&factor[j]).map(|(a, b)| a * b).sum())?;
            response[j][i] = response[i][j];
        }
    }
    Ok(response)
}
fn mean(a: Vec3, b: Vec3) -> Vec3 {
    a.scale(0.5).add(b.scale(0.5))
}
fn finite(value: f64) -> Result<f64, RigidDmtError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(RigidDmtError::Invalid(
            "unrepresentable native DMT mechanics",
        ))
    }
}
