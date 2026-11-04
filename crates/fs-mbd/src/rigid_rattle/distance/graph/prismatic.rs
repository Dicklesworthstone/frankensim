//! Finite sliding rails through original moving-plane and pure-couple owners.

use super::{
    Dd, DistanceStepError, Geometry, MassProperties, PointPlaneJoint, RigidBodyState, Row, RowKind,
    Vec3, positive,
};
use crate::UnitQuaternion;

/// A freely sliding joint between a finite carriage and finite rail carrier.
///
/// The first body's material point lies on the line through the second body's
/// rail anchor along its joint-frame Z axis. The two joint frames have the same
/// world orientation. Thus axial translation is free, while both transverse
/// translations and all three relative rotations are constrained.
///
/// Frames map joint coordinates into the respective native COM body frames.
/// They need not coincide with either body's principal axes. Both bodies retain
/// their actual full inertia and recoil; there is no prescribed support.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PrismaticJoint {
    bodies: [usize; 2],
    points_body_m: [Vec3; 2],
    frames_body: [UnitQuaternion; 2],
    transverse: [PointPlaneJoint; 2],
    angular_length_m: f64,
}

impl PrismaticJoint {
    /// Prepare a carriage/rail joint with a free local-Z displacement.
    ///
    /// `angular_length_m = L` converts the relative-frame sine vector to metres.
    /// On the admitted relative-angle branch `|theta| < pi/2`, graph tolerances
    /// `(p,v)` require `L*|sin(theta)| <= p` and the norm of the derivative of
    /// that same vector at most `v`. At alignment this rate norm is
    /// `L*|omega_second-omega_first|`. The complete three-row couple is bounded
    /// by `L*maximum_impulse_ns`; the complete two-row transverse force impulse
    /// has its separate original `maximum_impulse_ns` bound.
    ///
    /// The plane normals and offsets are prepared once. Their reaction arms
    /// follow the current sliding point, including the moving rail's torque.
    /// The joint does not infer travel limits, friction, an axial motor law or
    /// a zero-length radial constraint.
    ///
    /// # Errors
    /// Refuses equal body indices, nonfinite points and a nonpositive or
    /// nonfinite angular scale. The graph also requires actual distinct bodies,
    /// `L > position_tolerance_m`, the local orientation branch and independent
    /// rows. Supplied frame values already satisfy native quaternion admission.
    pub fn new(
        bodies: [usize; 2],
        points_body_m: [Vec3; 2],
        frames_body: [UnitQuaternion; 2],
        angular_length_m: f64,
    ) -> Result<Self, DistanceStepError> {
        if !positive(angular_length_m) {
            return Err(DistanceStepError::Invalid("prismatic angular length scale"));
        }
        let normals = [Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0)]
            .map(|axis| frames_body[1].rotate_body_to_world(axis));
        let transverse = [
            PointPlaneJoint::new(bodies, points_body_m[0], points_body_m[1], normals[0])?,
            PointPlaneJoint::new(bodies, points_body_m[0], points_body_m[1], normals[1])?,
        ];
        Ok(Self {
            bodies,
            points_body_m,
            frames_body,
            transverse,
            angular_length_m,
        })
    }

    /// Carriage followed by finite carrier, in the original graph body slices.
    #[must_use]
    pub const fn bodies(self) -> [usize; 2] {
        self.bodies
    }

    /// Carriage material point and authored rail anchor in their COM frames [m].
    #[must_use]
    pub const fn points_body_m(self) -> [Vec3; 2] {
        self.points_body_m
    }

    /// Joint-to-body frame maps. Their local Z axes define the free rail axis.
    #[must_use]
    pub const fn frames_body(self) -> [UnitQuaternion; 2] {
        self.frames_body
    }

    /// Authored scale for the three angular equations [m].
    #[must_use]
    pub const fn angular_length_m(self) -> f64 {
        self.angular_length_m
    }

    pub(super) fn write_rows(self, rows: &mut [Row]) {
        for (component, plane) in self.transverse.into_iter().enumerate() {
            rows[component] = plane.row(component, 2);
        }
        for component in 0..3 {
            rows[2 + component] = Row {
                bodies: self.bodies,
                points_body_m: [Vec3::ZERO; 2],
                kind: RowKind::FrameLock(AngularConstraint {
                    frames_body: self.frames_body,
                    length_m: self.angular_length_m,
                    component,
                }),
            };
        }
    }
}

#[derive(Clone, Copy)]
pub(super) struct AngularConstraint {
    pub frames_body: [UnitQuaternion; 2],
    pub length_m: f64,
    pub component: usize,
}

/// The scaled skew of R0^T R1, expressed in joint frame zero. With relative
/// quaternion q=(w,v), g=2 L w v. For a relative world angular velocity dw,
/// dg_i/dt = n_i.dot(dw), where
/// n_i = L R0 [w² e_i - v_i v + w(v cross e_i)].
/// These exact moving-frame Jacobians receive equal/opposite native couples.
/// The positive-cosine branch excludes the skew chart's pi-angle false zero
/// and its pi/2 rank loss. No angular finite difference is used in production.
pub(super) fn geometry(
    states: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    row: AngularConstraint,
) -> Result<Geometry, DistanceStepError> {
    geometry_with_spin(states, properties, row, None)
}

// A revolute actuator advances the first joint frame about its local Z axis.
// Subtract that commanded motion before the final DD dot product; a large
// shared or commanded spin must not erase its small hidden-constraint error.
pub(super) fn geometry_with_spin(
    states: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    row: AngularConstraint,
    spin: Option<(Vec3, Dd)>,
) -> Result<Geometry, DistanceStepError> {
    let frames = core::array::from_fn::<_, 2, _>(|i| {
        states[i].pose().orientation().multiply(row.frames_body[i])
    });
    let relative = frames[0].conjugate().multiply(frames[1]);
    let [w, x, y, z] = relative.components();
    if w * w <= 0.5 {
        return Err(DistanceStepError::LocalBranch);
    }
    let v = Vec3::new(x, y, z);
    let component = [x, y, z][row.component];
    let axis = [
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    ][row.component];
    let local = axis
        .scale(w * w)
        .sub(v.scale(component))
        .add(v.cross(axis).scale(w));
    let normal = frames[0].rotate_body_to_world(local).scale(row.length_m);
    let gap =
        (Dd::from_f64(2.0 * w) * Dd::from_f64(component) * Dd::from_f64(row.length_m)).to_f64();
    let angular = [
        states[0]
            .point_kinematics(properties[0], Vec3::ZERO)?
            .angular_velocity_world,
        states[1]
            .point_kinematics(properties[1], Vec3::ZERO)?
            .angular_velocity_world,
    ];
    let command = spin.map(|(axis, rate)| {
        (
            states[0].pose().orientation().rotate_body_to_world(axis),
            rate,
        )
    });
    let axis = command.map_or([0.0; 3], |(axis, _)| [axis.x, axis.y, axis.z]);
    let rate = [angular[0].x, angular[0].y, angular[0].z]
        .into_iter()
        .zip([angular[1].x, angular[1].y, angular[1].z])
        .zip([normal.x, normal.y, normal.z])
        .zip(axis)
        .fold(Dd::ZERO, |sum, (((before, after), n), axis)| {
            let mut relative = Dd::from_f64(after) - Dd::from_f64(before);
            if let Some((_, rate)) = command {
                relative = relative - Dd::from_f64(axis) * rate;
            }
            sum + relative * Dd::from_f64(n)
        })
        .to_f64();
    if !normal.is_finite() || !gap.is_finite() || !rate.is_finite() {
        return Err(DistanceStepError::Invalid("prismatic angular geometry"));
    }
    Ok(Geometry {
        normal,
        arms_body_m: [Vec3::ZERO; 2],
        gap,
        rate,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Pose;

    fn near(a: f64, b: f64, tolerance: f64) {
        assert!(
            (a - b).abs() <= tolerance,
            "{a:e} != {b:e}, tol={tolerance:e}"
        );
    }

    #[test]
    fn g1_frame_lock_jacobian_matches_independent_axis_rotation_differences() {
        let frames = [
            UnitQuaternion::from_axis_angle(Vec3::new(0.2, 0.3, 1.0), 0.21).unwrap(),
            UnitQuaternion::from_axis_angle(Vec3::new(-0.5, 1.0, 0.1), -0.13).unwrap(),
        ];
        let orientation = [
            UnitQuaternion::from_axis_angle(Vec3::new(1.0, -2.0, 0.4), 0.35).unwrap(),
            UnitQuaternion::from_axis_angle(Vec3::new(-0.4, 1.0, 0.7), -0.23).unwrap(),
        ];
        let properties = core::array::from_fn::<_, 2, _>(|i| {
            MassProperties::from_principal_frame(
                1.0 + i as f64,
                Vec3::ZERO,
                Vec3::new(0.7, 1.0, 1.3),
                frames[i],
            )
            .unwrap()
        });
        let omega = [Vec3::new(0.2, -0.4, 0.3), Vec3::new(-0.5, 0.1, 0.4)];
        let states = core::array::from_fn(|i| {
            RigidBodyState::new(
                Pose::new(Vec3::ZERO, orientation[i]).unwrap(),
                Vec3::ZERO,
                properties[i]
                    .angular_momentum_body_checked(orientation[i].rotate_world_to_body(omega[i]))
                    .unwrap(),
            )
            .unwrap()
        });
        let basis = [
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ];
        // This oracle uses four transformed axes, not relative-quaternion
        // components or the production Jacobian expression.
        let axis_gap = |states: [RigidBodyState; 2], component: usize| {
            let i = (component + 1) % 3;
            let j = (component + 2) % 3;
            let a = |axis| {
                states[0]
                    .pose()
                    .orientation()
                    .rotate_body_to_world(frames[0].rotate_body_to_world(axis))
            };
            let b = |axis| {
                states[1]
                    .pose()
                    .orientation()
                    .rotate_body_to_world(frames[1].rotate_body_to_world(axis))
            };
            0.65 * (a(basis[j]).dot(b(basis[i])) - a(basis[i]).dot(b(basis[j])))
        };
        let perturb = |mut states: [RigidBodyState; 2], body: usize, axis: Vec3, amount: f64| {
            let s = states[body];
            let q = s.pose().orientation();
            let q = q
                .right_exp(q.rotate_world_to_body(axis).scale(amount))
                .unwrap();
            states[body] = RigidBodyState::new(
                Pose::new(s.pose().position_world(), q).unwrap(),
                s.linear_momentum_world(),
                s.angular_momentum_body(),
            )
            .unwrap();
            states
        };
        for component in 0..3 {
            let g = geometry(
                states,
                properties,
                AngularConstraint {
                    frames_body: frames,
                    length_m: 1.3,
                    component,
                },
            )
            .unwrap();
            near(g.gap, axis_gap(states, component), 8e-16);
            let h = 3e-6;
            for (axis, normal) in basis.into_iter().zip([g.normal.x, g.normal.y, g.normal.z]) {
                for body in 0..2 {
                    let high = axis_gap(perturb(states, body, axis, h), component);
                    let low = axis_gap(perturb(states, body, axis, -h), component);
                    near(
                        (high - low) / (2.0 * h),
                        if body == 0 { -normal } else { normal },
                        2e-10,
                    );
                }
            }
            let high = perturb(perturb(states, 0, omega[0], h), 1, omega[1], h);
            let low = perturb(perturb(states, 0, omega[0], -h), 1, omega[1], -h);
            near(
                (axis_gap(high, component) - axis_gap(low, component)) / (2.0 * h),
                g.rate,
                2e-10,
            );
        }
    }
}
