use fs_mbd::rigid_rattle::distance::graph::{
    ArticulatedJoint, DistanceGraphControls, DistanceGraphError, DistanceLink, GraphBodyImpulse,
    GraphImpulseControls, GraphImpulseResponse, PointPlaneJoint, PrismaticJoint, RevoluteJoint,
    SphericalJoint, apply_articulated_graph_impulses,
};
use fs_mbd::rigid_rattle::distance::{DistanceStepControls, DistanceStepError};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};

const Z: Vec3 = Vec3::new(0.0, 0.0, 1.0);

fn controls() -> GraphImpulseControls {
    GraphImpulseControls {
        graph: DistanceGraphControls {
            step: DistanceStepControls {
                maximum_impulse_ns: 100.0,
                maximum_rotation_increment_rad: 0.1,
                maximum_position_iterations: 20,
                maximum_backtracks: 8,
            },
            position_tolerance_m: 1e-10,
            velocity_tolerance_mps: 1e-10,
            minimum_scaled_pivot: 1e-10,
            maximum_bodies: 64,
            maximum_constraints: 32,
            maximum_trial_evaluations: 100,
            maximum_live_bytes: 2_000_000,
        },
        work_tolerance_j: 2e-8,
        linear_impulse_tolerance_ns: 2e-9,
        angular_impulse_tolerance_nms: 2e-8,
    }
}
fn mass(m: f64, inertia: f64) -> MassProperties {
    MassProperties::new(m, Vec3::ZERO, Vec3::new(inertia, inertia, inertia)).unwrap()
}
fn state(
    property: MassProperties,
    position: Vec3,
    q: UnitQuaternion,
    velocity: Vec3,
) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(position, q).unwrap(),
        velocity.scale(property.mass()),
        Vec3::ZERO,
    )
    .unwrap()
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() <= 3e-9, "{a:e} != {b:e}");
}
fn vector(a: Vec3, b: Vec3) {
    near(a.x, b.x);
    near(a.y, b.y);
    near(a.z, b.z);
}
fn run(
    states: &[RigidBodyState],
    properties: &[MassProperties],
    joints: &[ArticulatedJoint],
    sources: &[GraphBodyImpulse],
) -> GraphImpulseResponse {
    apply_articulated_graph_impulses(states, properties, joints, sources, controls(), || false)
        .unwrap()
}
fn rod(a: usize, b: usize) -> ArticulatedJoint {
    ArticulatedJoint::Distance(DistanceLink {
        bodies: [a, b],
        points_body_m: [Vec3::ZERO; 2],
        length_m: 1.0,
    })
}

// G0/G3: one shock transfers through the complete finite chain in one event;
// an independent total-mass solution also checks frame, origin and boost changes.
#[test]
fn chain_impulse_shares_momentum_and_whole_event_work() {
    let properties = [mass(1.0, 1.0), mass(2.0, 2.0), mass(3.0, 3.0)];
    for (q, origin, boost) in [
        (UnitQuaternion::IDENTITY, Vec3::ZERO, Vec3::ZERO),
        (
            UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, 3.0), 0.7).unwrap(),
            Vec3::new(13.0, -29.0, 41.0),
            Vec3::new(2.0, -3.0, 0.5),
        ),
    ] {
        let direction = q.rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0));
        let states: Vec<_> = properties
            .iter()
            .enumerate()
            .map(|(i, &p)| state(p, origin.add(direction.scale(i as f64)), q, boost))
            .collect();
        let sources = [
            GraphBodyImpulse {
                linear_world_ns: direction.scale(3.0),
                angular_body_nms: Vec3::ZERO,
            },
            GraphBodyImpulse::ZERO,
            GraphBodyImpulse::ZERO,
        ];
        let result = run(&states, &properties, &[rod(0, 1), rod(1, 2)], &sources);
        for i in 0..3 {
            assert_eq!(result.states_after[i].pose(), states[i].pose());
            vector(
                result.states_after[i]
                    .center_of_mass_velocity_world(properties[i])
                    .unwrap(),
                boost.add(direction.scale(0.5)),
            );
            vector(result.states_after[i].angular_momentum_body(), Vec3::ZERO);
            assert_eq!(result.native_impulses[i].state_before, states[i]);
            assert_eq!(
                result.native_impulses[i].state_after,
                result.states_after[i]
            );
        }
        near(result.multipliers_ns[0], 2.5);
        near(result.multipliers_ns[1], 1.5);
        near(
            result.applied_work_j,
            0.75 + boost.dot(direction.scale(3.0)),
        );
        near(result.kinetic_energy_change_j, result.applied_work_j);
        near(result.constraint_work_j, 0.0);
        near(result.work_energy_residual_j, 0.0);
        vector(result.linear_momentum_change_world_ns, direction.scale(3.0));
        vector(
            result.angular_momentum_change_world_nms,
            origin.cross(direction.scale(3.0)),
        );
        vector(result.constraint_angular_impulse_world_nms, Vec3::ZERO);
        assert!(
            result
                .residuals_after
                .iter()
                .all(|r| r.range_rate_mps.abs() < 1e-10)
        );
    }
}

// G0: a ball transfers all translational impulses but does not weld the spin.
// The source is an actual off-center point impulse, not an invented COM force.
#[test]
fn spherical_joint_transmits_point_load_recoil_without_welding_spin() {
    let properties = [mass(1.0, 2.0), mass(2.0, 4.0)];
    let states = properties.map(|p| state(p, Vec3::ZERO, UnitQuaternion::IDENTITY, Vec3::ZERO));
    let source = GraphBodyImpulse::at_body_point(
        states[0],
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(3.0, 6.0, 0.0),
    )
    .unwrap();
    let joints = [ArticulatedJoint::Spherical(SphericalJoint {
        bodies: [0, 1],
        points_body_m: [Vec3::ZERO; 2],
    })];
    let result = run(
        &states,
        &properties,
        &joints,
        &[source, GraphBodyImpulse::ZERO],
    );
    for i in 0..2 {
        vector(
            result.states_after[i]
                .center_of_mass_velocity_world(properties[i])
                .unwrap(),
            Vec3::new(1.0, 2.0, 0.0),
        );
    }
    vector(
        result.states_after[0].angular_momentum_body(),
        Vec3::new(0.0, 0.0, -3.0),
    );
    vector(result.states_after[1].angular_momentum_body(), Vec3::ZERO);
    // Translation: J²/(2*(m0+m1)) = 7.5; free spin: 3²/(2*I0) = 2.25.
    near(result.kinetic_energy_change_j, 9.75);
    near(result.applied_work_j, 9.75);
    near(result.constraint_work_j, 0.0);
    vector(
        result.angular_momentum_change_world_nms,
        Vec3::new(0.0, 0.0, -3.0),
    );
}

// G0/G3: shared hinge tilt uses both full inertia tensors. The independently
// summed 2x2 tensor gives omega_x=omega_y=1; axial spin remains free.
#[test]
fn hinge_uses_full_tensor_reaction_and_preserves_free_axial_spin() {
    let properties = [
        MassProperties::from_inertia_tensor(
            1.0,
            Vec3::ZERO,
            fs_ga::Mat3 {
                m: [2.0, 0.5, 0.0, 0.5, 3.0, 0.0, 0.0, 0.0, 4.0],
            },
        )
        .unwrap(),
        MassProperties::from_inertia_tensor(
            2.0,
            Vec3::ZERO,
            fs_ga::Mat3 {
                m: [5.0, -0.25, 0.0, -0.25, 4.0, 0.0, 0.0, 0.0, 6.0],
            },
        )
        .unwrap(),
    ];
    for q in [
        UnitQuaternion::IDENTITY,
        UnitQuaternion::from_axis_angle(Vec3::new(2.0, -1.0, 3.0), 0.9).unwrap(),
    ] {
        let states = properties.map(|p| state(p, Vec3::ZERO, q, Vec3::ZERO));
        let joint = ArticulatedJoint::Revolute(
            RevoluteJoint::new([0, 1], [Vec3::ZERO; 2], [Z; 2], 1.0).unwrap(),
        );
        let source = GraphBodyImpulse {
            linear_world_ns: Vec3::ZERO,
            angular_body_nms: Vec3::new(7.25, 7.25, 8.0),
        };
        let result = run(
            &states,
            &properties,
            &[joint],
            &[source, GraphBodyImpulse::ZERO],
        );
        vector(
            result.states_after[0].angular_momentum_body(),
            Vec3::new(2.5, 3.5, 8.0),
        );
        vector(
            result.states_after[1].angular_momentum_body(),
            Vec3::new(4.75, 3.75, 0.0),
        );
        near(result.kinetic_energy_change_j, 15.25);
        near(result.applied_work_j, 15.25);
        near(result.constraint_work_j, 0.0);
        vector(
            result.angular_momentum_change_world_nms,
            q.rotate_body_to_world(source.angular_body_nms),
        );
    }
}

// G0: carried planes and rails react only in their authored constrained
// directions; finite carriers recoil and the free rail coordinate survives.
#[test]
fn plane_and_prismatic_joints_retain_their_physical_freedoms() {
    let properties = [mass(1.0, 2.0), mass(2.0, 4.0)];
    let states = properties.map(|p| state(p, Vec3::ZERO, UnitQuaternion::IDENTITY, Vec3::ZERO));
    let source = GraphBodyImpulse {
        linear_world_ns: Vec3::new(0.0, 6.0, 3.0),
        angular_body_nms: Vec3::new(0.0, 0.0, 12.0),
    };
    let plane = ArticulatedJoint::PointPlane(
        PointPlaneJoint::new([0, 1], Vec3::ZERO, Vec3::ZERO, Vec3::new(0.0, 1.0, 0.0)).unwrap(),
    );
    let rail = ArticulatedJoint::Prismatic(
        PrismaticJoint::new([0, 1], [Vec3::ZERO; 2], [UnitQuaternion::IDENTITY; 2], 1.0).unwrap(),
    );
    for (joint, expected_spins, energy) in [(plane, [12.0, 0.0], 46.5), (rail, [4.0, 8.0], 22.5)] {
        let result = run(
            &states,
            &properties,
            &[joint],
            &[source, GraphBodyImpulse::ZERO],
        );
        vector(
            result.states_after[0]
                .center_of_mass_velocity_world(properties[0])
                .unwrap(),
            Vec3::new(0.0, 2.0, 3.0),
        );
        vector(
            result.states_after[1]
                .center_of_mass_velocity_world(properties[1])
                .unwrap(),
            Vec3::new(0.0, 2.0, 0.0),
        );
        for i in 0..2 {
            vector(
                result.states_after[i].angular_momentum_body(),
                Z.scale(expected_spins[i]),
            );
        }
        near(result.kinetic_energy_change_j, energy);
        near(result.applied_work_j, energy);
        near(result.constraint_work_j, 0.0);
        vector(
            result.linear_momentum_change_world_ns,
            source.linear_world_ns,
        );
        vector(
            result.angular_momentum_change_world_nms,
            source.angular_body_nms,
        );
    }
}

// G0: an impulse response is not a projection of existing admitted defects,
// and the source probe must not be rejected as if it were the final endpoint.
#[test]
fn zero_source_is_identity_and_unconstrained_sources_use_native_resultants() {
    let properties = [mass(1.0, 1.0), mass(2.0, 2.0)];
    let states = [
        state(
            properties[0],
            Vec3::ZERO,
            UnitQuaternion::IDENTITY,
            Vec3::ZERO,
        ),
        state(
            properties[1],
            Vec3::new(1.0, 0.0, 0.0),
            UnitQuaternion::IDENTITY,
            Vec3::new(5e-11, 0.0, 0.0),
        ),
    ];
    let identity = run(
        &states,
        &properties,
        &[rod(0, 1)],
        &[GraphBodyImpulse::ZERO; 2],
    );
    assert_eq!(identity.states_after, states);
    assert_eq!(identity.multipliers_ns, [0.0]);
    assert_eq!(identity.residuals_before, identity.residuals_after);
    let source = GraphBodyImpulse {
        linear_world_ns: Vec3::new(3.0, -2.0, 1.0),
        angular_body_nms: Z,
    };
    let free = run(&states, &properties, &[], &[source, GraphBodyImpulse::ZERO]);
    assert_eq!(
        free.native_impulses[0],
        states[0]
            .apply_momentum_impulse(
                properties[0],
                source.linear_world_ns,
                source.angular_body_nms
            )
            .unwrap()
    );
    assert_eq!(free.states_after[1], states[1]);
    assert!(run(&[], &[], &[], &[]).states_after.is_empty());
}

// G0: dependent rows, undeliverable local momentum, resource limits and every
// cancellation checkpoint refuse atomically; a retry gives the same event.
#[test]
fn rank_resource_momentum_and_cancellation_refusals_are_atomic() {
    let properties = [mass(1.0, 1.0), mass(2.0, 2.0)];
    let states = [
        state(
            properties[0],
            Vec3::ZERO,
            UnitQuaternion::IDENTITY,
            Vec3::ZERO,
        ),
        state(
            properties[1],
            Vec3::new(1.0, 0.0, 0.0),
            UnitQuaternion::IDENTITY,
            Vec3::ZERO,
        ),
    ];
    let sources = [
        GraphBodyImpulse {
            linear_world_ns: Vec3::new(3.0, 0.0, 0.0),
            angular_body_nms: Vec3::ZERO,
        },
        GraphBodyImpulse::ZERO,
    ];
    assert!(matches!(
        apply_articulated_graph_impulses(
            &states,
            &properties,
            &[rod(0, 1), rod(0, 1)],
            &sources,
            controls(),
            || false
        ),
        Err(DistanceGraphError::DependentConstraints { .. })
    ));
    let mut small = controls();
    small.graph.maximum_live_bytes = 0;
    assert!(matches!(
        apply_articulated_graph_impulses(
            &states,
            &properties,
            &[rod(0, 1)],
            &sources,
            small,
            || false
        ),
        Err(DistanceGraphError::Limit {
            resource: "live bytes",
            ..
        })
    ));
    let huge = state(
        properties[0],
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
        Vec3::new(2.0_f64.powi(54), 0.0, 0.0),
    );
    assert!(matches!(
        apply_articulated_graph_impulses(
            &[huge],
            &properties[..1],
            &[],
            &[GraphBodyImpulse {
                linear_world_ns: Vec3::new(1.0, 0.0, 0.0),
                angular_body_nms: Vec3::ZERO
            }],
            controls(),
            || false
        ),
        Err(DistanceGraphError::Step(DistanceStepError::Invalid(_)))
    ));
    let mut checkpoints = 0;
    let expected = apply_articulated_graph_impulses(
        &states,
        &properties,
        &[rod(0, 1)],
        &sources,
        controls(),
        || {
            checkpoints += 1;
            false
        },
    )
    .unwrap();
    for stop in 1..=checkpoints {
        let mut seen = 0;
        assert!(matches!(
            apply_articulated_graph_impulses(
                &states,
                &properties,
                &[rod(0, 1)],
                &sources,
                controls(),
                || {
                    seen += 1;
                    seen == stop
                }
            ),
            Err(DistanceGraphError::Step(DistanceStepError::Cancelled))
        ));
    }
    assert_eq!(expected, run(&states, &properties, &[rod(0, 1)], &sources));
}

// G0: finite off-center joint mobility shares a transverse source between
// linear recoil and spin. Independent scalar mobility is 1+1+1+1=4.
#[test]
fn offset_ball_reactions_exchange_orbital_and_spin_momentum() {
    let properties = [mass(1.0, 1.0); 2];
    let states = [
        state(
            properties[0],
            Vec3::new(-1.0, 0.0, 0.0),
            UnitQuaternion::IDENTITY,
            Vec3::ZERO,
        ),
        state(
            properties[1],
            Vec3::new(1.0, 0.0, 0.0),
            UnitQuaternion::IDENTITY,
            Vec3::ZERO,
        ),
    ];
    let joint = ArticulatedJoint::Spherical(SphericalJoint {
        bodies: [0, 1],
        points_body_m: [Vec3::new(1.0, 0.0, 0.0), Vec3::new(-1.0, 0.0, 0.0)],
    });
    let source = GraphBodyImpulse {
        linear_world_ns: Vec3::new(0.0, 4.0, 0.0),
        angular_body_nms: Vec3::ZERO,
    };
    let result = run(
        &states,
        &properties,
        &[joint],
        &[source, GraphBodyImpulse::ZERO],
    );
    vector(
        result.states_after[0].linear_momentum_world(),
        Vec3::new(0.0, 3.0, 0.0),
    );
    vector(
        result.states_after[1].linear_momentum_world(),
        Vec3::new(0.0, 1.0, 0.0),
    );
    for after in &result.states_after {
        vector(after.angular_momentum_body(), Vec3::new(0.0, 0.0, -1.0));
    }
    vector(
        result.angular_momentum_change_world_nms,
        Vec3::new(0.0, 0.0, -4.0),
    );
    vector(result.constraint_angular_impulse_world_nms, Vec3::ZERO);
    near(result.applied_work_j, 6.0);
    near(result.kinetic_energy_change_j, 6.0);
    near(result.constraint_work_j, 0.0);
}
