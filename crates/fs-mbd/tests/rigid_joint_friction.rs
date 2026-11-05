use fs_contact::tangential::impact::multiple::{
    CoupledFrictionImpulseControls, CoupledFrictionImpulseError,
};
use fs_mbd::rigid_rattle::distance::graph::joint_friction::{
    JointFrictionControls, JointFrictionError, JointFrictionMotion, JointFrictionPhase,
    JointFrictionPort, JointFrictionReport, PrismaticFrictionInput, apply_joint_friction_impulses,
    apply_prismatic_friction_impulse,
};
use fs_mbd::rigid_rattle::distance::graph::revolute_friction::{
    RevoluteFrictionControls, RevoluteFrictionInput, apply_revolute_friction_impulse,
};
use fs_mbd::rigid_rattle::distance::graph::{
    ArticulatedJoint, DistanceGraphControls, DistanceGraphError, DistanceLink, PrismaticJoint,
    RevoluteJoint, advance_articulated_graph,
};
use fs_mbd::rigid_rattle::distance::{DistanceStepControls, DistanceStepError};
use fs_mbd::{
    Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState, UnitQuaternion, Vec3,
    Wrench,
};

const X: Vec3 = Vec3::new(1.0, 0.0, 0.0);
const Y: Vec3 = Vec3::new(0.0, 1.0, 0.0);
const Z: Vec3 = Vec3::new(0.0, 0.0, 1.0);

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tolerance {tolerance:e}"
    );
}
fn vector(actual: Vec3, expected: Vec3, tolerance: f64) {
    near(actual.x, expected.x, tolerance);
    near(actual.y, expected.y, tolerance);
    near(actual.z, expected.z, tolerance);
}
fn properties(mass: f64, inertia: f64) -> MassProperties {
    MassProperties::new(mass, Vec3::ZERO, Vec3::new(inertia, inertia, inertia)).unwrap()
}
fn state(property: MassProperties, position: Vec3, velocity: Vec3, spin: Vec3) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(position, UnitQuaternion::IDENTITY).unwrap(),
        velocity.scale(property.mass()),
        property.angular_momentum_body_checked(spin).unwrap(),
    )
    .unwrap()
}
fn controls() -> JointFrictionControls {
    JointFrictionControls {
        graph: DistanceGraphControls {
            step: DistanceStepControls {
                maximum_impulse_ns: 100.0,
                maximum_rotation_increment_rad: 0.25,
                maximum_position_iterations: 30,
                maximum_backtracks: 16,
            },
            position_tolerance_m: 1e-10,
            velocity_tolerance_mps: 1e-10,
            minimum_scaled_pivot: 1e-10,
            maximum_bodies: 32,
            maximum_constraints: 32,
            maximum_trial_evaluations: 160,
            maximum_live_bytes: 2_000_000,
        },
        solve: CoupledFrictionImpulseControls {
            maximum_iterations: 256,
            impulse_tolerance_ns: 1e-10,
            velocity_tolerance_m_s: 1e-10,
            energy_tolerance_j: 1e-9,
            relative_matrix_tolerance: 1e-12,
            maximum_impulse_ns: 100.0,
        },
        work_tolerance_j: 1e-9,
        linear_impulse_tolerance_ns: 1e-9,
        angular_impulse_tolerance_nms: 1e-9,
    }
}
fn hinge(bodies: [usize; 2], axis: Vec3, length: f64) -> ArticulatedJoint {
    ArticulatedJoint::Revolute(
        RevoluteJoint::new(bodies, [Vec3::ZERO; 2], [axis; 2], length).unwrap(),
    )
}
fn slider(bodies: [usize; 2], points: [Vec3; 2], length: f64) -> ArticulatedJoint {
    ArticulatedJoint::Prismatic(
        PrismaticJoint::new(bodies, points, [UnitQuaternion::IDENTITY; 2], length).unwrap(),
    )
}
fn revolute_port(joint_index: usize, capacity_nms: f64) -> JointFrictionPort {
    JointFrictionPort::Revolute {
        joint_index,
        capacity_nms,
    }
}
fn prismatic_port(joint_index: usize, capacity_ns: f64) -> JointFrictionPort {
    JointFrictionPort::Prismatic {
        joint_index,
        capacity_ns,
    }
}
fn run(
    states: &[RigidBodyState],
    properties: &[MassProperties],
    joints: &[ArticulatedJoint],
    ports: &[JointFrictionPort],
) -> JointFrictionReport {
    apply_joint_friction_impulses(states, properties, joints, ports, controls(), || false).unwrap()
}
fn spin(state: RigidBodyState, property: MassProperties) -> Vec3 {
    state
        .point_kinematics(property, Vec3::ZERO)
        .unwrap()
        .angular_velocity_world
}
fn velocity(state: RigidBodyState, property: MassProperties) -> Vec3 {
    state.center_of_mass_velocity_world(property).unwrap()
}
fn audit(before: &[RigidBodyState], properties: &[MassProperties], out: &JointFrictionReport) {
    let integrator = RigidBodyIntegrator::new(Gravity::ZERO);
    let mut energy = 0.0;
    let mut linear = Vec3::ZERO;
    let mut angular = Vec3::ZERO;
    for ((&before, &after), &p) in before.iter().zip(&out.states_after).zip(properties) {
        assert_eq!(before.pose(), after.pose());
        let a = integrator.diagnostics(before, p).unwrap();
        let b = integrator.diagnostics(after, p).unwrap();
        energy += b.mechanical_energy - a.mechanical_energy;
        linear = linear.add(b.linear_momentum_world.sub(a.linear_momentum_world));
        angular = angular.add(b.angular_momentum_world.sub(a.angular_momentum_world));
    }
    near(energy, out.kinetic_energy_change_j, 2e-10);
    near(
        energy,
        out.drive_work_j + out.friction_work_j + out.constraint_work_j,
        2e-10,
    );
    assert!(out.friction_work_j <= 0.0);
    near(out.work_energy_residual_j, 0.0, 2e-10);
    near(out.constraint_work_j, 0.0, 2e-10);
    near(out.native_duality_gap_j, 0.0, 2e-10);
    vector(linear, Vec3::ZERO, 2e-10);
    vector(angular, Vec3::ZERO, 2e-10);
    vector(out.linear_momentum_change_ns, linear, 2e-10);
    vector(out.angular_momentum_change_nms, angular, 2e-10);
    for port in &out.ports {
        near(port.scaled_velocity_residual_m_s, 0.0, 2e-10);
    }
}

#[test]
fn g1_three_rotors_solve_together_and_preserve_physical_angular_scaling() {
    let p = [
        properties(1.0, 2.0),
        properties(2.0, 3.0),
        properties(3.0, 5.0),
    ];
    let s = core::array::from_fn::<_, 3, _>(|i| {
        state(p[i], Vec3::ZERO, Vec3::ZERO, Z.scale([3.0, -1.0, 2.0][i]))
    });
    let mut reference = None;
    for lengths in [[1.0, 1.0], [0.25, 2.0], [3.0, 0.3]] {
        let joints = [hinge([0, 1], Z, lengths[0]), hinge([1, 2], Z, lengths[1])];
        let ports = [revolute_port(0, 1.0), revolute_port(1, 10.0)];
        let out = run(&s, &p, &joints, &ports);
        assert_eq!(out.phase, JointFrictionPhase::FrozenCapacities);
        assert_eq!(out, run(&s, &p, &joints, &[ports[1], ports[0]]));
        for (i, expected) in [2.5, 1.0, 1.0].into_iter().enumerate() {
            vector(spin(out.states_after[i], p[i]), Z.scale(expected), 5e-13);
            vector(velocity(out.states_after[i], p[i]), Vec3::ZERO, 5e-13);
            if let Some(previous) = &reference {
                let previous: &Vec<RigidBodyState> = previous;
                vector(
                    out.states_after[i].angular_momentum_body(),
                    previous[i].angular_momentum_body(),
                    5e-13,
                );
            }
        }
        for i in 0..2 {
            let JointFrictionMotion::Revolute {
                angular_length_m,
                rate_before_rad_s,
                rate_after_rad_s,
                friction_impulse_nms,
            } = out.ports[i].motion
            else {
                panic!("hinge port")
            };
            near(angular_length_m, lengths[i], 0.0);
            near(rate_before_rad_s, [-4.0, 3.0][i], 2e-14);
            near(rate_after_rad_s, [-1.5, 0.0][i], 5e-13);
            near(friction_impulse_nms, [1.0, -5.0][i], 5e-13);
            for j in 0..2 {
                let physical = [[5.0 / 6.0, -1.0 / 3.0], [-1.0 / 3.0, 8.0 / 15.0]][i][j];
                near(
                    out.scaled_mobility_per_kg[i][j],
                    lengths[i] * lengths[j] * physical,
                    5e-13,
                );
            }
        }
        near(out.friction_work_j, -10.25, 5e-13);
        audit(&s, &p, &out);
        reference = Some(out.states_after);
    }

    // A single pass of individually correct joint solves does not solve the
    // same frozen event: AB's later impulse reopens the previously stopped BC.
    let joints = [hinge([0, 1], Z, 1.0), hinge([1, 2], Z, 1.0)];
    let c = controls();
    let scalar = RevoluteFrictionControls {
        graph: c.graph,
        maximum_axial_impulse_nms: 100.0,
        angular_rate_tolerance_rad_s: c.solve.velocity_tolerance_m_s,
        work_tolerance_j: c.work_tolerance_j,
        linear_impulse_tolerance_ns: c.linear_impulse_tolerance_ns,
        angular_impulse_tolerance_nms: c.angular_impulse_tolerance_nms,
    };
    let bc = apply_revolute_friction_impulse(
        &s,
        &p,
        &joints,
        1,
        RevoluteFrictionInput::Sliding {
            kinetic_capacity_nms: 10.0,
        },
        scalar,
        || false,
    )
    .unwrap();
    let ab = apply_revolute_friction_impulse(
        &bc.states_after,
        &p,
        &joints,
        0,
        RevoluteFrictionInput::Sliding {
            kinetic_capacity_nms: 1.0,
        },
        scalar,
        || false,
    )
    .unwrap();
    let reopened = spin(ab.states_after[2], p[2]).z - spin(ab.states_after[1], p[1]).z;
    near(reopened, -1.0 / 3.0, 5e-13);
}

#[test]
fn g1_free_off_center_slider_uses_translational_and_rotational_recoil() {
    let p = [properties(2.0, 2.0), properties(3.0, 3.0)];
    let reduced_mass = 6.0 / 5.0;
    for axial in [0.0, 3.0] {
        let d = Vec3::new(2.0, 0.0, axial);
        let s = [
            state(p[0], d, Z.scale(9.0 / 5.0), Vec3::ZERO),
            state(p[1], Vec3::ZERO, Z.scale(-6.0 / 5.0), Vec3::ZERO),
        ];
        let joints = [slider([0, 1], [X.scale(-2.0), Vec3::ZERO], 0.4)];
        // Eliminate the common translation and conserve the complete angular
        // momentum of two co-rotating bodies. The reduced kinetic metric gives
        // W = 1/mu + d_perp^2/(I0+I1+mu*d_axial^2), independent of graph rows.
        let denominator = 5.0 + reduced_mass * axial * axial;
        let mobility = 1.0 / reduced_mass + 4.0 / denominator;
        for capacity in [0.0, 1.0, 10.0] {
            let out = apply_prismatic_friction_impulse(
                &s,
                &p,
                &joints,
                0,
                PrismaticFrictionInput::Sliding {
                    kinetic_capacity_ns: capacity,
                },
                controls(),
                || false,
            )
            .unwrap();
            let impulse = -capacity.min(3.0 / mobility);
            let rate = 3.0 + mobility * impulse;
            let omega = Y.scale(2.0 * impulse / denominator);
            let relative = omega.cross(d).add(Z.scale(rate));
            near(out.scaled_mobility_per_kg[0][0], mobility, 5e-13);
            vector(spin(out.states_after[0], p[0]), omega, 5e-13);
            vector(spin(out.states_after[1], p[1]), omega, 5e-13);
            vector(
                velocity(out.states_after[0], p[0]),
                relative.scale(3.0 / 5.0),
                5e-13,
            );
            vector(
                velocity(out.states_after[1], p[1]),
                relative.scale(-2.0 / 5.0),
                5e-13,
            );
            let JointFrictionMotion::Prismatic {
                velocity_before_m_s,
                velocity_after_m_s,
                friction_impulse_ns,
                drive_impulse_ns,
            } = out.ports[0].motion
            else {
                panic!("slider port")
            };
            near(velocity_before_m_s, 3.0, 5e-14);
            near(velocity_after_m_s, rate, 5e-13);
            near(friction_impulse_ns, impulse, 5e-13);
            assert_eq!(drive_impulse_ns, 0.0);
            assert_eq!(
                out.phase,
                if capacity < 3.0 / mobility {
                    JointFrictionPhase::Sliding
                } else {
                    JointFrictionPhase::Arrested
                }
            );
            near(out.friction_work_j, impulse * (3.0 + rate) * 0.5, 5e-13);
            audit(&s, &p, &out);
            let coupled = run(&s, &p, &joints, &[prismatic_port(0, capacity)]);
            for i in 0..2 {
                vector(
                    coupled.states_after[i].linear_momentum_world(),
                    out.states_after[i].linear_momentum_world(),
                    5e-13,
                );
                vector(
                    coupled.states_after[i].angular_momentum_body(),
                    out.states_after[i].angular_momentum_body(),
                    5e-13,
                );
            }
            // These are the existing body's actual states, immediately usable
            // by the native finite-motion owner and its constraint admission.
            let next = advance_articulated_graph(
                &out.states_after,
                &p,
                &joints,
                RigidBodyIntegrator::new(Gravity::ZERO),
                &[Wrench::ZERO; 2],
                0.0001,
                controls().graph,
                || false,
            )
            .unwrap();
            run(
                &next.step.states_after,
                &p,
                &joints,
                &[prismatic_port(0, 0.0)],
            );
        }
    }
}

#[test]
fn g1_slider_retains_products_of_inertia_in_its_native_rotational_response() {
    // The summed full inertia is [[5,1,0],[1,4,0],[0,0,6]]. Its off-diagonal
    // entry creates real X spin from the Y torque needed by an eccentric rail.
    let p = [0.4, 0.6].map(|share| {
        MassProperties::from_inertia_tensor(
            5.0 * share,
            Vec3::ZERO,
            fs_ga::Mat3 {
                m: [5.0, 1.0, 0.0, 1.0, 4.0, 0.0, 0.0, 0.0, 6.0].map(|x| x * share),
            },
        )
        .unwrap()
    });
    let s = [
        state(p[0], X.scale(2.0), Z.scale(9.0 / 5.0), Vec3::ZERO),
        state(p[1], Vec3::ZERO, Z.scale(-6.0 / 5.0), Vec3::ZERO),
    ];
    let joints = [slider([0, 1], [X.scale(-2.0), Vec3::ZERO], 1.0)];
    let out = run(&s, &p, &joints, &[prismatic_port(0, 10.0)]);
    // W=1/mu + 4*(I_sum^-1)_yy = 5/6 + 20/19.
    let mobility = 215.0 / 114.0;
    let impulse = -3.0 / mobility;
    let omega = Vec3::new(-2.0 / 19.0, 10.0 / 19.0, 0.0).scale(impulse);
    near(out.scaled_mobility_per_kg[0][0], mobility, 5e-13);
    for i in 0..2 {
        vector(spin(out.states_after[i], p[i]), omega, 5e-13);
    }
    assert!(omega.x.abs() > 0.1);
    near(out.friction_work_j, 1.5 * impulse, 5e-13);
    audit(&s, &p, &out);
}

fn mixed(
    lengths: [f64; 2],
) -> (
    [RigidBodyState; 3],
    [MassProperties; 3],
    [ArticulatedJoint; 2],
) {
    let p = [
        properties(1.0, 2.0),
        properties(1.0, 3.0),
        properties(2.0, 2.0),
    ];
    let s = [
        state(p[0], Vec3::ZERO, Vec3::ZERO, Y.scale(3.0)),
        state(p[1], Vec3::ZERO, Vec3::ZERO, Y.scale(-1.0)),
        state(
            p[2],
            Vec3::new(2.0, 0.0, 1.0),
            Vec3::new(-1.0, 0.0, 5.0),
            Y.scale(-1.0),
        ),
    ];
    (
        s,
        p,
        [
            hinge([0, 1], Y, lengths[0]),
            slider([2, 1], [X.scale(-2.0), Vec3::ZERO], lengths[1]),
        ],
    )
}

#[test]
fn g1_mixed_hinge_and_slider_share_the_complete_constrained_mobility() {
    for lengths in [[1.0, 1.0], [0.2, 3.0]] {
        let (s, p, joints) = mixed(lengths);
        let ports = [revolute_port(0, 1.0), prismatic_port(1, 10.0)];
        let out = run(&s, &p, &joints, &ports);
        // Independently reduce kinetic energy to (a,b,u):
        // T_internal = a^2 + 5*b^2 - 2*b*u + u^2/2.
        // Port velocities are (b-a,u), hence physical W=[[2/3,1/3],[1/3,5/3]].
        let expected = [
            [2.0 / 3.0 * lengths[0].powi(2), lengths[0] / 3.0],
            [lengths[0] / 3.0, 5.0 / 3.0],
        ];
        for i in 0..2 {
            for j in 0..2 {
                near(out.scaled_mobility_per_kg[i][j], expected[i][j], 8e-13);
            }
        }
        let JointFrictionMotion::Revolute {
            rate_after_rad_s,
            friction_impulse_nms,
            ..
        } = out.ports[0].motion
        else {
            panic!("hinge port")
        };
        near(friction_impulse_nms, 1.0, 5e-13);
        near(rate_after_rad_s, -4.0, 5e-13);
        let JointFrictionMotion::Prismatic {
            velocity_after_m_s,
            friction_impulse_ns,
            ..
        } = out.ports[1].motion
        else {
            panic!("slider port")
        };
        near(friction_impulse_ns, -2.0, 8e-13);
        near(velocity_after_m_s, 0.0, 8e-13);
        for i in 0..3 {
            vector(
                spin(out.states_after[i], p[i]),
                Y.scale([2.5, -1.5, -1.5][i]),
                8e-13,
            );
            let v = if i == 2 {
                Vec3::new(-1.25, 0.0, 4.0)
            } else {
                Vec3::new(0.25, 0.0, 1.0)
            };
            vector(velocity(out.states_after[i], p[i]), v, 8e-13);
        }
        near(out.friction_work_j, -7.0, 8e-13);
        audit(&s, &p, &out);
        assert_eq!(out, run(&s, &p, &joints, &[ports[1], ports[0]]));
        let reversed = run(
            &s,
            &p,
            &[joints[1], joints[0]],
            &[prismatic_port(0, 10.0), revolute_port(1, 1.0)],
        );
        for i in 0..3 {
            vector(
                reversed.states_after[i].linear_momentum_world(),
                out.states_after[i].linear_momentum_world(),
                8e-13,
            );
            vector(
                reversed.states_after[i].angular_momentum_body(),
                out.states_after[i].angular_momentum_body(),
                8e-13,
            );
        }
    }
}

#[test]
fn g1_coupled_port_work_can_transfer_energy_while_total_friction_is_passive() {
    let (mut s, p, joints) = mixed([1.0, 1.0]);
    s[0] = state(p[0], Vec3::ZERO, Vec3::ZERO, Y.scale(-1.25));
    let out = run(
        &s,
        &p,
        &joints,
        &[revolute_port(0, 1.0), prismatic_port(1, 10.0)],
    );
    // Initial port rates (1/4,3) and the independent mixed metric give
    // p=(7/12,-23/12). Both stop together: the first impulse temporarily
    // transfers energy into its relative mode, with larger loss in the second.
    let JointFrictionMotion::Revolute {
        friction_impulse_nms,
        rate_after_rad_s,
        ..
    } = out.ports[0].motion
    else {
        panic!("hinge port")
    };
    let JointFrictionMotion::Prismatic {
        friction_impulse_ns,
        velocity_after_m_s,
        ..
    } = out.ports[1].motion
    else {
        panic!("slider port")
    };
    near(friction_impulse_nms, 7.0 / 12.0, 5e-13);
    near(friction_impulse_ns, -23.0 / 12.0, 5e-13);
    near(rate_after_rad_s, 0.0, 5e-13);
    near(velocity_after_m_s, 0.0, 5e-13);
    near(out.ports[0].friction_work_j, 7.0 / 96.0, 5e-13);
    assert!(out.ports[0].friction_work_j > 0.0);
    near(out.friction_work_j, -269.0 / 96.0, 5e-13);
    audit(&s, &p, &out);
}

#[test]
fn g1_mixed_native_solution_is_covariant_under_world_rotation_origin_and_boost() {
    let (s, p, joints) = mixed([0.7, 1.3]);
    let ports = [revolute_port(0, 1.0), prismatic_port(1, 10.0)];
    let base = run(&s, &p, &joints, &ports);
    let rotation = UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, 3.0), 0.61).unwrap();
    let offset = Vec3::new(5.0, -7.0, 11.0);
    let boost = Vec3::new(3.0, -4.0, 2.0);
    let transform = |s: RigidBodyState, p: MassProperties| {
        RigidBodyState::new(
            Pose::new(
                rotation
                    .rotate_body_to_world(s.pose().position_world())
                    .add(offset),
                rotation,
            )
            .unwrap(),
            rotation
                .rotate_body_to_world(s.linear_momentum_world())
                .add(boost.scale(p.mass())),
            s.angular_momentum_body(),
        )
        .unwrap()
    };
    let transformed = core::array::from_fn::<_, 3, _>(|i| transform(s[i], p[i]));
    let out = run(&transformed, &p, &joints, &ports);
    for i in 0..3 {
        let expected = transform(base.states_after[i], p[i]);
        assert_eq!(out.states_after[i].pose(), expected.pose());
        vector(
            out.states_after[i].linear_momentum_world(),
            expected.linear_momentum_world(),
            5e-12,
        );
        vector(
            out.states_after[i].angular_momentum_body(),
            expected.angular_momentum_body(),
            5e-12,
        );
    }
    for i in 0..2 {
        vector(
            out.ports[i].axis_world,
            rotation.rotate_body_to_world(base.ports[i].axis_world),
            5e-14,
        );
        for j in 0..2 {
            near(
                out.scaled_mobility_per_kg[i][j],
                base.scaled_mobility_per_kg[i][j],
                2e-13,
            );
        }
    }
    near(out.friction_work_j, base.friction_work_j, 5e-12);
    audit(&transformed, &p, &out);
}

#[test]
fn g1_slider_resting_drive_holds_breaks_away_and_never_captures_existing_slip() {
    let p = [properties(2.0, 2.0), properties(3.0, 3.0)];
    let s = p.map(|p| state(p, Vec3::ZERO, Vec3::ZERO, Vec3::ZERO));
    let joints = [slider([0, 1], [Vec3::ZERO; 2], 1.0)];
    for drive in [-1.2, -1.0, -0.3, 0.0, 0.3, 1.0, 1.2] {
        let out = apply_prismatic_friction_impulse(
            &s,
            &p,
            &joints,
            0,
            PrismaticFrictionInput::RestingDrive {
                drive_impulse_ns: drive,
                static_capacity_ns: 1.0,
                kinetic_capacity_ns: 0.4,
            },
            controls(),
            || false,
        )
        .unwrap();
        let held = drive.abs() <= 1.0;
        let friction = if held { -drive } else { -drive.signum() * 0.4 };
        let resultant = drive + friction;
        assert_eq!(
            out.phase,
            if held {
                JointFrictionPhase::Held
            } else {
                JointFrictionPhase::Sliding
            }
        );
        if held {
            assert_eq!(out.states_after, s);
            assert_eq!(out.drive_work_j, 0.0);
            assert_eq!(out.friction_work_j, 0.0);
        }
        let JointFrictionMotion::Prismatic {
            velocity_after_m_s,
            friction_impulse_ns,
            ..
        } = out.ports[0].motion
        else {
            panic!("slider port")
        };
        near(velocity_after_m_s, 5.0 / 6.0 * resultant, 2e-14);
        near(friction_impulse_ns, friction, 2e-14);
        near(out.drive_work_j, drive * velocity_after_m_s * 0.5, 2e-14);
        near(
            out.friction_work_j,
            friction * velocity_after_m_s * 0.5,
            2e-14,
        );
        audit(&s, &p, &out);
    }
    let moving = [state(p[0], Vec3::ZERO, Z.scale(1e-12), Vec3::ZERO), s[1]];
    assert!(matches!(
        apply_prismatic_friction_impulse(
            &moving,
            &p,
            &joints,
            0,
            PrismaticFrictionInput::RestingDrive {
                drive_impulse_ns: 0.0,
                static_capacity_ns: 1.0,
                kinetic_capacity_ns: 0.4,
            },
            controls(),
            || false,
        ),
        Err(JointFrictionError::NotResting { .. })
    ));
    let stopped = apply_prismatic_friction_impulse(
        &moving,
        &p,
        &joints,
        0,
        PrismaticFrictionInput::Sliding {
            kinetic_capacity_ns: 0.4,
        },
        controls(),
        || false,
    )
    .unwrap();
    let JointFrictionMotion::Prismatic {
        friction_impulse_ns,
        velocity_after_m_s,
        ..
    } = stopped.ports[0].motion
    else {
        panic!("slider port")
    };
    near(friction_impulse_ns, -1.2e-12, 1e-26);
    near(velocity_after_m_s, 0.0, 1e-26);
    assert_eq!(stopped.phase, JointFrictionPhase::Arrested);
}

#[test]
fn g1_original_graph_bounds_locked_ports_and_cancellation_refuse_atomically() {
    let (s, p, joints) = mixed([1.0, 1.0]);
    let ports = [revolute_port(0, 1.0), prismatic_port(1, 10.0)];
    let good = run(&s, &p, &joints, &ports);
    for invalid in [
        vec![ports[0], ports[0]],
        vec![prismatic_port(0, 1.0)],
        vec![revolute_port(2, 1.0)],
        vec![revolute_port(0, -1.0)],
        vec![prismatic_port(1, f64::NAN)],
    ] {
        assert!(matches!(
            apply_joint_friction_impulses(&s, &p, &joints, &invalid, controls(), || false),
            Err(JointFrictionError::Invalid(_))
        ));
    }
    for changed in 0..4 {
        let mut c = controls();
        match changed {
            0 => c.graph.maximum_bodies = 2,
            1 => c.graph.maximum_constraints = 9,
            2 => c.graph.maximum_live_bytes = good.planned_live_bytes - 1,
            _ => c.graph.step.maximum_impulse_ns = 1e-4,
        }
        assert!(apply_joint_friction_impulses(&s, &p, &joints, &ports, c, || false).is_err());
    }
    let mut c = controls();
    c.solve.maximum_impulse_ns = 1.5;
    assert!(apply_joint_friction_impulses(&s, &p, &joints, &ports, c, || false).is_err());
    c = controls();
    c.graph.maximum_live_bytes = good.planned_live_bytes;
    assert_eq!(
        good,
        apply_joint_friction_impulses(&s, &p, &joints, &ports, c, || false).unwrap()
    );

    let locked_p = [properties(1.0, 1.0); 2];
    let locked_s = [
        state(locked_p[0], Z, Vec3::ZERO, Vec3::ZERO),
        state(locked_p[1], Vec3::ZERO, Vec3::ZERO, Vec3::ZERO),
    ];
    let locked_joints = [
        slider([0, 1], [Vec3::ZERO; 2], 1.0),
        ArticulatedJoint::Distance(DistanceLink {
            bodies: [0, 1],
            points_body_m: [Vec3::ZERO; 2],
            length_m: 1.0,
        }),
    ];
    assert!(matches!(
        run_error(
            &locked_s,
            &locked_p,
            &locked_joints,
            &[prismatic_port(0, 1.0)]
        ),
        JointFrictionError::Unresolved(_)
    ));

    let mut polls = 0;
    apply_joint_friction_impulses(&s, &p, &joints, &ports, controls(), || {
        polls += 1;
        false
    })
    .unwrap();
    assert!(polls > 100);
    for stop in 1..=polls {
        let mut seen = 0;
        let error = apply_joint_friction_impulses(&s, &p, &joints, &ports, controls(), || {
            seen += 1;
            seen == stop
        })
        .unwrap_err();
        assert!(
            matches!(
                error,
                JointFrictionError::Native(DistanceGraphError::Step(DistanceStepError::Cancelled))
                    | JointFrictionError::Coupled(CoupledFrictionImpulseError::Cancelled)
            ),
            "checkpoint {stop}/{polls}: {error:?}"
        );
    }
    assert_eq!(good, run(&s, &p, &joints, &ports));
}

fn run_error(
    states: &[RigidBodyState],
    properties: &[MassProperties],
    joints: &[ArticulatedJoint],
    ports: &[JointFrictionPort],
) -> JointFrictionError {
    apply_joint_friction_impulses(states, properties, joints, ports, controls(), || false)
        .unwrap_err()
}

#[test]
fn g1_empty_and_zero_capacity_events_are_native_identities() {
    let (s, p, joints) = mixed([1.0, 1.0]);
    for ports in [vec![], vec![revolute_port(0, 0.0), prismatic_port(1, 0.0)]] {
        let out = run(&s, &p, &joints, &ports);
        assert_eq!(out.states_after, s);
        assert_eq!(out.kinetic_energy_change_j, 0.0);
        assert_eq!(out.friction_work_j, 0.0);
        assert_eq!(out.constraint_work_j, 0.0);
        assert_eq!(out.ports.len(), ports.len());
        if ports.is_empty() {
            assert_eq!(out.projection.as_ref().unwrap().factorizations, 0);
        }
    }
}
