use fs_mbd::assembly::release::PointMassRelease;
use fs_mbd::assembly::transfer::network::{
    PointMassTransferEdge, PointMassTransferNetworkResult, transfer_point_mass_network,
};
use fs_mbd::assembly::transfer::{
    PointMassTransfer, PointMassTransferControls, transfer_point_masses,
};
use fs_mbd::assembly::{AssemblyControls, AssemblyError};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};

fn controls(bodies: usize, edges: usize) -> PointMassTransferControls {
    PointMassTransferControls {
        assembly: AssemblyControls {
            maximum_parts: bodies + edges,
            maximum_impulses: 2 * edges,
            linear_momentum_tolerance_ns: 1e-9,
            angular_momentum_tolerance_nms: 1e-8,
            energy_tolerance_j: 1e-8,
        },
        attachment_tolerance_m: 1e-11,
        mass_tolerance_kg: 1e-12,
    }
}
fn close(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:.17e} != {b:.17e}");
}
fn vector(a: Vec3, b: Vec3, tol: f64) {
    close(a.x, b.x, tol);
    close(a.y, b.y, tol);
    close(a.z, b.z, tol);
}
fn properties(mass: f64) -> MassProperties {
    MassProperties::new(mass, Vec3::ZERO, Vec3::new(5.0, 6.0, 7.0)).unwrap()
}
fn body(p: MassProperties, x: Vec3, q: UnitQuaternion, v: Vec3, w: Vec3) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(x, q).unwrap(),
        v.scale(p.mass()),
        p.angular_momentum_body_checked(w).unwrap(),
    )
    .unwrap()
}
fn central(p: MassProperties, speed: f64) -> RigidBodyState {
    body(
        p,
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
        Vec3::new(speed, 0.0, 0.0),
        Vec3::ZERO,
    )
}
fn edge(
    donor: usize,
    receiver: usize,
    mass: f64,
    from: Vec3,
    to: Vec3,
    impulse: Vec3,
) -> PointMassTransferEdge {
    PointMassTransferEdge {
        donor,
        receiver,
        transfer: PointMassTransfer {
            release: PointMassRelease {
                mass_kg: mass,
                attachment_arm_body_m: from,
                impulse_on_parcel_body_ns: impulse,
            },
            receiver_attachment_arm_body_m: to,
        },
    }
}
fn totals(b: &[RigidBodyState], p: &[MassProperties]) -> (f64, Vec3, Vec3, f64) {
    let observer = Vec3::new(-0.7, 1.3, -0.9);
    let mut mass = 0.0;
    let mut linear = Vec3::ZERO;
    let mut angular = Vec3::ZERO;
    let mut kinetic = 0.0;
    for (b, p) in b.iter().zip(p) {
        mass += p.mass();
        linear = linear.add(b.linear_momentum_world());
        angular = angular
            .add(
                b.pose()
                    .orientation()
                    .rotate_body_to_world(b.angular_momentum_body()),
            )
            .add(
                b.pose()
                    .position_world()
                    .sub(observer)
                    .cross(b.linear_momentum_world()),
            );
        kinetic += b.kinetic_energy(*p).unwrap();
    }
    (mass, linear, angular, kinetic)
}
fn conservation(
    b: &[RigidBodyState],
    p: &[MassProperties],
    result: &PointMassTransferNetworkResult,
) {
    let before = totals(b, p);
    let after = totals(
        &result.bodies.iter().map(|e| e.state).collect::<Vec<_>>(),
        &result
            .bodies
            .iter()
            .map(|e| e.properties)
            .collect::<Vec<_>>(),
    );
    close(after.0, before.0, 2e-12);
    vector(after.1, before.1, 2e-10);
    vector(after.2, before.2, 2e-9);
    close(
        after.3 - before.3,
        result.supplied_work_j - result.removed_relative_energy_j,
        3e-9,
    );
}

#[test]
fn simultaneous_cycle_keeps_original_donor_snapshots() {
    let p = [properties(4.0); 3];
    let b = [central(p[0], 4.0), central(p[1], -2.0), central(p[2], 1.0)];
    let edges = core::array::from_fn::<_, 3, _>(|i| {
        edge(i, (i + 1) % 3, 1.0, Vec3::ZERO, Vec3::ZERO, Vec3::ZERO)
    });
    let out = transfer_point_mass_network(&b, &p, &edges, controls(3, 3), || false).unwrap();
    // Each node retains 3/4 of its original momentum and receives 1/4
    // of its predecessor's ORIGINAL momentum, even around the cycle.
    for (endpoint, expected) in out.bodies.iter().zip([13.0, -2.0, 1.0]) {
        close(endpoint.properties.mass(), 4.0, 0.0);
        vector(
            endpoint.state.linear_momentum_world(),
            Vec3::new(expected, 0.0, 0.0),
            0.0,
        );
    }
    close(out.supplied_work_j, 0.0, 0.0);
    close(out.removed_relative_energy_j, 81.0 / 4.0, 2e-14);
    for (endpoint, loss) in out.bodies.iter().zip([27.0 / 8.0, 27.0 / 2.0, 27.0 / 8.0]) {
        close(endpoint.removed_relative_energy_j, loss, 1e-14);
    }
    conservation(&b, &p, &out);

    // Sequential pair handoffs are different events: newly received material
    // changes the velocity used by a later outgoing edge.
    let mut sequential = b;
    let mut sequential_properties = p;
    for e in edges {
        let result = transfer_point_masses(
            &sequential[e.donor],
            sequential_properties[e.donor],
            &sequential[e.receiver],
            sequential_properties[e.receiver],
            &[e.transfer],
            controls(2, 1),
            || false,
        )
        .unwrap();
        sequential[e.donor] = result.donor_state;
        sequential_properties[e.donor] = result.donor_properties;
        sequential[e.receiver] = result.receiver.state;
        sequential_properties[e.receiver] = result.receiver.properties;
    }
    close(sequential[0].linear_momentum_world().x, 12.64, 5e-14);
    assert!(
        (sequential[0].linear_momentum_world().x - out.bodies[0].state.linear_momentum_world().x)
            .abs()
            > 0.3
    );
}

#[test]
fn joint_incoming_streams_share_receiver_recoil() {
    let p = [properties(2.0); 3];
    let edges = [
        edge(0, 2, 1.0, Vec3::ZERO, Vec3::ZERO, Vec3::ZERO),
        edge(1, 2, 1.0, Vec3::ZERO, Vec3::ZERO, Vec3::ZERO),
    ];
    for (second, expected_momentum, expected_loss) in [(-2.0, 0.0, 4.0), (2.0, 4.0, 2.0)] {
        let b = [
            central(p[0], 2.0),
            central(p[1], second),
            central(p[2], 0.0),
        ];
        let out = transfer_point_mass_network(&b, &p, &edges, controls(3, 2), || false).unwrap();
        assert_eq!(out.bodies[2].properties.mass(), 4.0);
        vector(
            out.bodies[2].state.linear_momentum_world(),
            Vec3::new(expected_momentum, 0.0, 0.0),
            0.0,
        );
        close(
            out.bodies[2].removed_relative_energy_j,
            expected_loss,
            2e-14,
        );
        close(out.removed_relative_energy_j, expected_loss, 2e-14);
        conservation(&b, &p, &out);
    }
}

fn square() -> (
    [RigidBodyState; 4],
    [MassProperties; 4],
    [PointMassTransferEdge; 4],
) {
    let p = [properties(4.0); 4];
    let positions = [
        Vec3::ZERO,
        Vec3::new(2.0, 0.0, 0.0),
        Vec3::new(2.0, 2.0, 0.0),
        Vec3::new(0.0, 2.0, 0.0),
    ];
    let b = core::array::from_fn(|i| {
        body(
            p[i],
            positions[i],
            UnitQuaternion::from_axis_angle(Vec3::new(0.2, -0.4, 0.7), 0.3 * i as f64).unwrap(),
            Vec3::new(0.3 * i as f64 - 0.4, 0.2 - 0.1 * i as f64, 0.05 * i as f64),
            Vec3::new(0.03 * i as f64, -0.04, 0.1 + 0.02 * i as f64),
        )
    });
    let edges = core::array::from_fn(|i| {
        let j = (i + 1) % 4;
        let point = positions[i].add(positions[j]).scale(0.5);
        edge(
            i,
            j,
            1.0,
            b[i].pose()
                .orientation()
                .rotate_world_to_body(point.sub(positions[i])),
            b[j].pose()
                .orientation()
                .rotate_world_to_body(point.sub(positions[j])),
            Vec3::new(0.1, -0.04, 0.03).scale(1.0 + i as f64),
        )
    });
    (b, p, edges)
}

#[test]
fn off_center_cycle_rebases_each_donor_and_conserves_full_spatial_momentum() {
    let (b, p, edges) = square();
    let out = transfer_point_mass_network(&b, &p, &edges, controls(4, 4), || false).unwrap();
    conservation(&b, &p, &out);
    vector(
        out.bodies[0].center_shift_body_m,
        Vec3::new(-0.25, 0.25, 0.0),
        3e-15,
    );
    let expected = [
        23.0 / 4.0,
        -0.25,
        0.0,
        -0.25,
        19.0 / 4.0,
        0.0,
        0.0,
        0.0,
        13.0 / 2.0,
    ];
    for (a, e) in out.bodies[0]
        .properties
        .inertia_tensor_body()
        .m
        .into_iter()
        .zip(expected)
    {
        close(a, e, 3e-12);
    }
    for (original, endpoint) in b.iter().zip(&out.bodies) {
        assert_eq!(endpoint.properties.mass(), 4.0);
        assert_eq!(
            endpoint.state.pose().orientation(),
            original.pose().orientation()
        );
        vector(
            endpoint.state.pose().position_world(),
            original
                .pose()
                .point_world_from_body(endpoint.center_shift_body_m)
                .unwrap(),
            3e-14,
        );
    }
    assert!(out.supplied_work_j > 0.0 && out.removed_relative_energy_j > 0.0);
    close(out.mass_residual_kg, 0.0, 0.0);
}

#[test]
fn body_edge_permutations_and_split_parcels_preserve_the_network_event() {
    let (b, p, edges) = square();
    let run = |b: &[RigidBodyState], p: &[MassProperties], e: &[PointMassTransferEdge]| {
        transfer_point_mass_network(b, p, e, controls(b.len(), e.len()), || false).unwrap()
    };
    let out = run(&b, &p, &edges);
    let order = [2, 0, 3, 1];
    let inverse = [1, 3, 0, 2];
    let reordered_b = order.map(|i| b[i]);
    let reordered_p = order.map(|i| p[i]);
    let reordered_e = [edges[3], edges[1], edges[0], edges[2]].map(|e| PointMassTransferEdge {
        donor: inverse[e.donor],
        receiver: inverse[e.receiver],
        ..e
    });
    let changed = run(&reordered_b, &reordered_p, &reordered_e);
    for (original, endpoint) in order.into_iter().zip(&changed.bodies) {
        vector(
            endpoint.state.linear_momentum_world(),
            out.bodies[original].state.linear_momentum_world(),
            2e-13,
        );
        vector(
            endpoint.state.angular_momentum_body(),
            out.bodies[original].state.angular_momentum_body(),
            4e-12,
        );
    }
    close(changed.supplied_work_j, out.supplied_work_j, 4e-13);
    close(
        changed.removed_relative_energy_j,
        out.removed_relative_energy_j,
        4e-12,
    );
    let split = |fraction: f64| PointMassTransferEdge {
        transfer: PointMassTransfer {
            release: PointMassRelease {
                mass_kg: fraction,
                impulse_on_parcel_body_ns: edges[0]
                    .transfer
                    .release
                    .impulse_on_parcel_body_ns
                    .scale(fraction),
                ..edges[0].transfer.release
            },
            ..edges[0].transfer
        },
        ..edges[0]
    };
    let divided = run(
        &b,
        &p,
        &[split(0.25), edges[1], edges[2], edges[3], split(0.75)],
    );
    for (a, c) in out.bodies.iter().zip(&divided.bodies) {
        vector(
            a.state.linear_momentum_world(),
            c.state.linear_momentum_world(),
            2e-13,
        );
        vector(
            a.state.angular_momentum_body(),
            c.state.angular_momentum_body(),
            5e-12,
        );
        vector(
            a.state.pose().position_world(),
            c.state.pose().position_world(),
            3e-14,
        );
    }
    close(divided.supplied_work_j, out.supplied_work_j, 4e-13);
    close(
        divided.removed_relative_energy_j,
        out.removed_relative_energy_j,
        5e-12,
    );
    conservation(&b, &p, &divided);
}

#[test]
fn canonical_pair_parity_isolated_owners_and_common_boost_keep_small_exchanges() {
    let p = [properties(1.0), properties(98.0), properties(49.0)];
    let isolated = body(
        p[0],
        Vec3::new(1e300, -1e300, 0.0),
        UnitQuaternion::IDENTITY,
        Vec3::new(1e300, 0.0, 0.0),
        Vec3::ZERO,
    );
    let b = [isolated, central(p[1], 0.0), central(p[2], 0.0)];
    let e = edge(1, 2, 49.0, Vec3::ZERO, Vec3::ZERO, Vec3::new(1.0, 0.0, 0.0));
    let out = transfer_point_mass_network(&b, &p, &[e], controls(3, 1), || false).unwrap();
    let pair = transfer_point_masses(
        &b[1],
        p[1],
        &b[2],
        p[2],
        &[e.transfer],
        controls(2, 1),
        || false,
    )
    .unwrap();
    assert_eq!(out.bodies[0].state, isolated);
    assert_eq!(out.bodies[0].properties, p[0]);
    assert_eq!(out.bodies[1].state, pair.donor_state);
    assert_eq!(out.bodies[2].state, pair.receiver.state);
    assert_eq!(
        out.bodies[2].state.linear_momentum_world(),
        Vec3::new(1.0, 0.0, 0.0)
    );
    assert_eq!(out.supplied_work_j, pair.supplied_work_j);
    assert_eq!(
        out.removed_relative_energy_j,
        pair.receiver.removed_relative_energy_j
    );

    let p = [properties(4.0); 3];
    let b = [central(p[0], 4.0), central(p[1], -2.0), central(p[2], 1.0)];
    let edges = core::array::from_fn::<_, 3, _>(|i| {
        edge(i, (i + 1) % 3, 1.0, Vec3::ZERO, Vec3::ZERO, Vec3::ZERO)
    });
    let base = transfer_point_mass_network(&b, &p, &edges, controls(3, 3), || false).unwrap();
    let boost = Vec3::new(1_099_511_627_776.0, 0.0, 0.0);
    let moved = core::array::from_fn::<_, 3, _>(|i| {
        RigidBodyState::new(
            b[i].pose(),
            b[i].linear_momentum_world().add(boost.scale(p[i].mass())),
            b[i].angular_momentum_body(),
        )
        .unwrap()
    });
    let result = transfer_point_mass_network(&moved, &p, &edges, controls(3, 3), || false).unwrap();
    assert_eq!(
        result.removed_relative_energy_j,
        base.removed_relative_energy_j
    );
    assert_eq!(result.supplied_work_j, base.supplied_work_j);
    assert_eq!(result.balance.linear_momentum_residual_world_ns, Vec3::ZERO);
    assert_eq!(
        result.balance.angular_momentum_residual_world_nms,
        Vec3::ZERO
    );
    close(result.balance.energy_residual_j, 0.0, 2e-12);
}

#[test]
fn network_limits_invalid_intermediate_inventory_and_late_cancel_are_atomic() {
    let (b, p, edges) = square();
    let empty = transfer_point_mass_network(&b, &p, &[], controls(4, 0), || false).unwrap();
    for (i, e) in empty.bodies.iter().enumerate() {
        assert_eq!(e.state, b[i]);
        assert_eq!(e.properties, p[i]);
        assert_eq!(e.center_shift_body_m, Vec3::ZERO);
    }
    let zero = edge(0, 1, 0.0, Vec3::ZERO, Vec3::ZERO, Vec3::ZERO);
    let no_op = transfer_point_mass_network(&b, &p, &[zero], controls(4, 1), || false).unwrap();
    assert_eq!(no_op, empty);
    assert!(matches!(
        transfer_point_mass_network(&b, &p, &[zero], controls(4, 0), || false),
        Err(AssemblyError::Budget {
            requested: 5,
            maximum: 4
        })
    ));
    let mut c = controls(4, 4);
    c.assembly.maximum_impulses = 7;
    assert!(matches!(
        transfer_point_mass_network(&b, &p, &edges, c, || false),
        Err(AssemblyError::Budget {
            requested: 8,
            maximum: 7
        })
    ));
    for invalid in [
        PointMassTransferEdge {
            donor: 4,
            ..edges[0]
        },
        PointMassTransferEdge {
            receiver: 0,
            ..edges[0]
        },
        edge(0, 1, -1.0, Vec3::ZERO, Vec3::ZERO, Vec3::ZERO),
        edge(0, 1, 0.0, Vec3::ZERO, Vec3::ZERO, Vec3::new(1.0, 0.0, 0.0)),
        edge(
            0,
            1,
            0.0,
            Vec3::ZERO,
            Vec3::new(f64::NAN, 0.0, 0.0),
            Vec3::ZERO,
        ),
    ] {
        assert!(transfer_point_mass_network(&b, &p, &[invalid], controls(4, 1), || false).is_err());
    }
    assert!(transfer_point_mass_network(&[], &[], &[], controls(1, 0), || false).is_err());
    assert!(transfer_point_mass_network(&b, &p[..3], &edges, controls(4, 4), || false).is_err());

    // Incoming replacement does not authorize a zero-inertia intermediate
    // body: all original outgoing releases must be physically admitted first.
    let p2 = [properties(2.0); 2];
    let b2 = [central(p2[0], 1.0), central(p2[1], -1.0)];
    let draining = [
        edge(0, 1, 2.0, Vec3::ZERO, Vec3::ZERO, Vec3::ZERO),
        edge(1, 0, 1.0, Vec3::ZERO, Vec3::ZERO, Vec3::ZERO),
    ];
    assert!(transfer_point_mass_network(&b2, &p2, &draining, controls(2, 2), || false).is_err());
    let separated = PointMassTransferEdge {
        transfer: PointMassTransfer {
            receiver_attachment_arm_body_m: edges[0]
                .transfer
                .receiver_attachment_arm_body_m
                .add(Vec3::new(0.01, 0.0, 0.0)),
            ..edges[0].transfer
        },
        ..edges[0]
    };
    assert!(transfer_point_mass_network(&b, &p, &[separated], controls(4, 1), || false).is_err());
    let mut count = 0;
    let expected = transfer_point_mass_network(&b, &p, &edges, controls(4, 4), || {
        count += 1;
        false
    })
    .unwrap();
    assert!(count > 12 * edges.len());
    let mut polls = 0;
    assert_eq!(
        transfer_point_mass_network(&b, &p, &edges, controls(4, 4), || {
            polls += 1;
            polls == count
        }),
        Err(AssemblyError::Cancelled)
    );
    assert_eq!(polls, count);
    assert_eq!(
        transfer_point_mass_network(&b, &p, &edges, controls(4, 4), || false).unwrap(),
        expected
    );
    assert_eq!(b, square().0);
    assert_eq!(p, square().1);
}
