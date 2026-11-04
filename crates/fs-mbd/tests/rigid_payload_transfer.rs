use fs_mbd::assembly::transfer::PointMassTransferControls;
use fs_mbd::assembly::transfer::rigid::{PreparedRigidPayloadTransfer, RigidPayloadTransferResult};
use fs_mbd::assembly::{
    AssemblyControls, AssemblyError, AssemblyPart, RigidAssembly, SeparationImpulse,
};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};

fn controls() -> PointMassTransferControls {
    PointMassTransferControls {
        assembly: AssemblyControls {
            maximum_parts: 8,
            maximum_impulses: 8,
            linear_momentum_tolerance_ns: 1e-10,
            angular_momentum_tolerance_nms: 1e-9,
            energy_tolerance_j: 1e-9,
        },
        attachment_tolerance_m: 2e-12,
        mass_tolerance_kg: 1e-12,
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}
fn vector(a: Vec3, b: Vec3, tolerance: f64) {
    close(a.x, b.x, tolerance);
    close(a.y, b.y, tolerance);
    close(a.z, b.z, tolerance);
}
fn rotation(axis: Vec3, angle: f64) -> UnitQuaternion {
    UnitQuaternion::from_axis_angle(axis, angle).unwrap()
}
fn compose(a: UnitQuaternion, b: UnitQuaternion) -> UnitQuaternion {
    let [aw, ax, ay, az] = a.components();
    let [bw, bx, by, bz] = b.components();
    UnitQuaternion::new(
        aw * bw - ax * bx - ay * by - az * bz,
        aw * bx + ax * bw + ay * bz - az * by,
        aw * by + ay * bw + az * bx - ax * bz,
        aw * bz + az * bw + ax * by - ay * bx,
    )
    .unwrap()
}
fn state(p: MassProperties, x: Vec3, q: UnitQuaternion, v: Vec3, w: Vec3) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(x, q).unwrap(),
        v.scale(p.mass()),
        p.angular_momentum_body_checked(w).unwrap(),
    )
    .unwrap()
}
fn diagonal(m: f64, inertia: Vec3) -> MassProperties {
    MassProperties::new(m, Vec3::ZERO, inertia).unwrap()
}
fn components(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn diagonal_product(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(a.x * b.x, a.y * b.y, a.z * b.z)
}
fn diagonal_divide(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(a.x / b.x, a.y / b.y, a.z / b.z)
}
fn world_spin(body: RigidBodyState) -> Vec3 {
    body.pose()
        .orientation()
        .rotate_body_to_world(body.angular_momentum_body())
}

// Independent 3x3 tensor arithmetic and cofactor inverse; no native inertia
// construction, diagonalization, angular-velocity solve or balance receipt.
type Matrix = [[f64; 3]; 3];
fn add(a: Matrix, b: Matrix) -> Matrix {
    std::array::from_fn(|i| std::array::from_fn(|j| a[i][j] + b[i][j]))
}
fn multiply(a: Matrix, v: Vec3) -> Vec3 {
    let v = components(v);
    let w: [f64; 3] = std::array::from_fn(|i| (0..3).map(|j| a[i][j] * v[j]).sum());
    Vec3::new(w[0], w[1], w[2])
}
fn inverse(a: Matrix) -> Matrix {
    let cofactor: Matrix = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            a[(i + 1) % 3][(j + 1) % 3] * a[(i + 2) % 3][(j + 2) % 3]
                - a[(i + 1) % 3][(j + 2) % 3] * a[(i + 2) % 3][(j + 1) % 3]
        })
    });
    let determinant: f64 = (0..3).map(|j| a[0][j] * cofactor[0][j]).sum();
    assert!(determinant > 0.0);
    std::array::from_fn(|i| std::array::from_fn(|j| cofactor[j][i] / determinant))
}
fn shifted(m: f64, r: Vec3) -> Matrix {
    let norm = r.dot(r);
    let r = components(r);
    std::array::from_fn(|i| {
        std::array::from_fn(|j| m * (if i == j { norm } else { 0.0 } - r[i] * r[j]))
    })
}
fn tensor(properties: MassProperties, orientation: UnitQuaternion) -> Matrix {
    let principal = components(properties.principal_inertia_body());
    let axes = [
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    ]
    .map(|v| {
        components(
            orientation
                .rotate_body_to_world(properties.principal_axes_to_body().rotate_body_to_world(v)),
        )
    });
    std::array::from_fn(|i| {
        std::array::from_fn(|j| (0..3).map(|k| principal[k] * axes[k][i] * axes[k][j]).sum())
    })
}
fn matrix(a: Matrix, b: Matrix, tolerance: f64) {
    for i in 0..3 {
        for j in 0..3 {
            close(a[i][j], b[i][j], tolerance);
        }
    }
}
fn totals(bodies: &[(RigidBodyState, MassProperties)]) -> (f64, Vec3, Vec3, f64) {
    let origin = Vec3::new(-0.75, 0.5, 1.25);
    let mut mass = 0.0;
    let mut linear = Vec3::ZERO;
    let mut angular = Vec3::ZERO;
    let mut energy = 0.0;
    for &(body, properties) in bodies {
        let p = body.linear_momentum_world();
        let spin = world_spin(body);
        mass += properties.mass();
        linear = linear.add(p);
        angular = angular
            .add(spin)
            .add(body.pose().position_world().sub(origin).cross(p));
        energy += 0.5
            * (p.dot(p) / properties.mass()
                + spin.dot(multiply(
                    inverse(tensor(properties, body.pose().orientation())),
                    spin,
                )));
    }
    (mass, linear, angular, energy)
}
fn conserve(before: &[(RigidBodyState, MassProperties)], out: &RigidPayloadTransferResult) {
    let a = totals(before);
    let b = totals(&[
        (out.donor_state, out.donor_properties),
        (out.receiver.state, out.receiver.assembly.properties()),
    ]);
    close(a.0, b.0, 2e-13);
    vector(a.1, b.1, 2e-11);
    vector(a.2, b.2, 2e-10);
    close(
        b.3 - a.3,
        out.supplied_work_j - out.receiver.removed_relative_energy_j,
        2e-10,
    );
}

#[test]
fn finite_spin_recoil_work_and_sticking_loss_match_closed_form() {
    let c = controls();
    let id = Vec3::new(3.0, 4.0, 5.0);
    let ip = Vec3::new(1.0, 1.5, 2.0);
    let ir = Vec3::new(4.0, 5.0, 6.0);
    let retained = diagonal(4.0, id);
    let payload = diagonal(2.0, ip);
    let receiver_properties = diagonal(6.0, ir);
    let parts = [retained, payload].map(|properties| AssemblyPart {
        properties,
        pose_in_reference: Pose::new(Vec3::ZERO, UnitQuaternion::IDENTITY).unwrap(),
    });
    let assembly = RigidAssembly::new(&parts, c.assembly, || false).unwrap();
    let plan = PreparedRigidPayloadTransfer::new(
        &assembly,
        &[1],
        Vec3::ZERO,
        Vec3::ZERO,
        c.assembly,
        || false,
    )
    .unwrap();
    let vd = Vec3::new(1.0, -0.5, 0.25);
    let wd = Vec3::new(0.3, -0.4, 0.6);
    let vr = Vec3::new(-0.2, 0.6, 0.4);
    let wr = Vec3::new(-0.1, 0.2, -0.3);
    let donor = state(
        assembly.properties(),
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
        vd,
        wd,
    );
    let receiver = state(
        receiver_properties,
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
        vr,
        wr,
    );
    let j = Vec3::new(0.4, -0.6, 0.2);
    let k = Vec3::new(0.3, -0.1, 0.2);
    let impulse = SeparationImpulse {
        part_a: 1,
        part_b: 0,
        point_assembly_m: Vec3::ZERO,
        impulse_on_a_assembly_ns: j,
        couple_on_a_assembly_nms: k,
    };
    let out = plan
        .transfer(
            &donor,
            &receiver,
            receiver_properties,
            &[impulse],
            c,
            || false,
        )
        .unwrap();
    vector(
        out.donor_state.linear_momentum_world(),
        vd.scale(4.0).sub(j),
        2e-14,
    );
    vector(
        out.donor_state.angular_momentum_body(),
        diagonal_product(id, wd).sub(k),
        2e-14,
    );
    vector(
        out.receiver.state.linear_momentum_world(),
        vr.scale(6.0).add(vd.scale(2.0)).add(j),
        2e-14,
    );
    let expected_spin = diagonal_product(ir, wr)
        .add(diagonal_product(ip, wd))
        .add(k);
    vector(world_spin(out.receiver.state), expected_spin, 2e-14);
    vector(
        out.receiver
            .assembly
            .properties()
            .angular_velocity_body_checked(expected_spin)
            .unwrap(),
        diagonal_divide(expected_spin, ir.add(ip)),
        2e-14,
    );
    let work = 0.5
        * (j.dot(j) * (1.0 / 4.0 + 1.0 / 2.0)
            + k.dot(diagonal_divide(k, id).add(diagonal_divide(k, ip))));
    let dv = vd.add(j.scale(0.5)).sub(vr);
    let dw = wd.add(diagonal_divide(k, ip)).sub(wr);
    let reduced_inertia = diagonal_divide(diagonal_product(ip, ir), ip.add(ir));
    let loss = 0.5 * (1.5 * dv.dot(dv) + dw.dot(diagonal_product(reduced_inertia, dw)));
    close(out.supplied_work_j, work, 2e-14);
    close(out.receiver.removed_relative_energy_j, loss, 2e-14);
    close(out.donor_properties.mass(), 4.0, 0.0);
    close(out.receiver.assembly.properties().mass(), 8.0, 0.0);
    conserve(
        &[
            (donor, assembly.properties()),
            (receiver, receiver_properties),
        ],
        &out,
    );
}

struct Scene {
    assembly: RigidAssembly,
    plan: PreparedRigidPayloadTransfer,
    donor: RigidBodyState,
    receiver: RigidBodyState,
    receiver_properties: MassProperties,
    impulses: [SeparationImpulse; 2],
}
fn scene() -> Scene {
    let c = controls();
    let parts: Vec<_> = [2.0, 3.0, 1.5, 2.5]
        .into_iter()
        .enumerate()
        .map(|(i, mass)| {
            let t = i as f64;
            AssemblyPart {
                properties: MassProperties::from_principal_frame(
                    mass,
                    Vec3::ZERO,
                    Vec3::new(0.4, 0.6, 0.8).scale(mass),
                    rotation(Vec3::new(1.0, -0.2, 0.5), 0.25 + t * 0.13),
                )
                .unwrap(),
                pose_in_reference: Pose::new(
                    Vec3::new(-0.6 + 0.5 * t, 0.15 * (t - 1.0).powi(2), -0.1 * t),
                    rotation(Vec3::new(0.2, 1.0, -0.4), -0.35 + 0.2 * t),
                )
                .unwrap(),
            }
        })
        .collect();
    let assembly = RigidAssembly::new(&parts, c.assembly, || false).unwrap();
    let donor_arm = Vec3::new(0.25, 0.5, -0.3);
    let receiver_arm = Vec3::new(-0.35, 0.15, 0.4);
    let plan = PreparedRigidPayloadTransfer::new(
        &assembly,
        &[3, 1],
        donor_arm,
        receiver_arm,
        c.assembly,
        || false,
    )
    .unwrap();
    let donor = state(
        assembly.properties(),
        Vec3::new(1.0, -2.0, 0.5),
        rotation(Vec3::new(0.4, -0.7, 1.0), 0.6),
        Vec3::new(0.7, -1.2, 0.4),
        Vec3::new(0.3, 0.6, -0.8),
    );
    let receiver_properties = MassProperties::from_principal_frame(
        5.0,
        Vec3::ZERO,
        Vec3::new(1.1, 1.4, 1.9),
        rotation(Vec3::new(-0.6, 1.0, 0.2), -0.4),
    )
    .unwrap();
    let qr = rotation(Vec3::new(0.5, 0.1, 1.0), -0.5);
    let xr = donor
        .pose()
        .point_world_from_body(donor_arm)
        .unwrap()
        .sub(qr.rotate_body_to_world(receiver_arm));
    let receiver = state(
        receiver_properties,
        xr,
        qr,
        Vec3::new(-0.2, 0.4, 0.7),
        Vec3::new(-0.4, 0.2, 0.1),
    );
    let impulses = [
        SeparationImpulse {
            part_a: 1,
            part_b: 0,
            point_assembly_m: Vec3::new(0.1, 0.3, -0.4),
            impulse_on_a_assembly_ns: Vec3::new(0.6, -0.25, 0.3),
            couple_on_a_assembly_nms: Vec3::new(0.1, 0.05, -0.08),
        },
        SeparationImpulse {
            part_a: 3,
            part_b: 2,
            point_assembly_m: Vec3::new(-0.1, 0.6, 0.2),
            impulse_on_a_assembly_ns: Vec3::new(-0.2, 0.35, 0.1),
            couple_on_a_assembly_nms: Vec3::new(-0.03, 0.02, 0.07),
        },
    ];
    Scene {
        assembly,
        plan,
        donor,
        receiver,
        receiver_properties,
        impulses,
    }
}

#[test]
fn off_center_full_tensors_close_spin_orbit_and_material_inventory() {
    let s = scene();
    let c = controls();
    let out = s
        .plan
        .transfer(
            &s.donor,
            &s.receiver,
            s.receiver_properties,
            &s.impulses,
            c,
            || false,
        )
        .unwrap();
    let q = s.donor.pose().orientation();
    let vd = s
        .donor
        .center_of_mass_velocity_world(s.assembly.properties())
        .unwrap();
    let wd = q.rotate_body_to_world(Vec3::new(0.3, 0.6, -0.8));
    let mut group_mass = [0.0; 2];
    let mut group_center = [Vec3::ZERO; 2];
    let mut group_tensor = [[[0.0; 3]; 3]; 2];
    // Reconstruct both groups directly from the original material components.
    for (i, part) in s.assembly.parts().iter().enumerate() {
        let g = usize::from(i == 1 || i == 3);
        group_mass[g] += part.properties.mass();
        group_center[g] = group_center[g].add(
            part.pose_in_reference
                .position_world()
                .scale(part.properties.mass()),
        );
    }
    for g in 0..2 {
        group_center[g] = group_center[g].scale(1.0 / group_mass[g]);
    }
    for (i, part) in s.assembly.parts().iter().enumerate() {
        let g = usize::from(i == 1 || i == 3);
        let r =
            q.rotate_body_to_world(part.pose_in_reference.position_world().sub(group_center[g]));
        group_tensor[g] = add(
            group_tensor[g],
            add(
                tensor(
                    part.properties,
                    compose(q, part.pose_in_reference.orientation()),
                ),
                shifted(part.properties.mass(), r),
            ),
        );
    }
    let mut impulses = [Vec3::ZERO; 2];
    let mut couples = [Vec3::ZERO; 2];
    for impulse in s.impulses {
        for (g, sign) in [(1, 1.0), (0, -1.0)] {
            let j = impulse.impulse_on_a_assembly_ns.scale(sign);
            impulses[g] = impulses[g].add(q.rotate_body_to_world(j));
            couples[g] = couples[g].add(
                q.rotate_body_to_world(
                    impulse
                        .point_assembly_m
                        .sub(group_center[g])
                        .cross(j)
                        .add(impulse.couple_on_a_assembly_nms.scale(sign)),
                ),
            );
        }
    }
    let group_positions = group_center.map(|r| {
        s.donor
            .pose()
            .position_world()
            .add(q.rotate_body_to_world(r))
    });
    let group_p: [Vec3; 2] = std::array::from_fn(|g| {
        vd.add(wd.cross(q.rotate_body_to_world(group_center[g])))
            .scale(group_mass[g])
            .add(impulses[g])
    });
    let group_h: [Vec3; 2] = std::array::from_fn(|g| multiply(group_tensor[g], wd).add(couples[g]));
    let expected_work: f64 = (0..2)
        .map(|g| {
            0.5 * (impulses[g].dot(impulses[g]) / group_mass[g]
                + couples[g].dot(multiply(inverse(group_tensor[g]), couples[g])))
        })
        .sum();
    vector(out.donor_state.linear_momentum_world(), group_p[0], 2e-13);
    vector(world_spin(out.donor_state), group_h[0], 2e-13);
    matrix(
        tensor(out.donor_properties, out.donor_state.pose().orientation()),
        group_tensor[0],
        2e-12,
    );
    let mr = s.receiver_properties.mass();
    let xr = s.receiver.pose().position_world();
    let combined_mass = mr + group_mass[1];
    let center = xr.add(
        group_positions[1]
            .sub(xr)
            .scale(group_mass[1] / combined_mass),
    );
    let receiver_tensor = tensor(s.receiver_properties, s.receiver.pose().orientation());
    let combined_tensor = add(
        add(receiver_tensor, group_tensor[1]),
        add(
            shifted(mr, xr.sub(center)),
            shifted(group_mass[1], group_positions[1].sub(center)),
        ),
    );
    let p = s.receiver.linear_momentum_world().add(group_p[1]);
    let h = world_spin(s.receiver)
        .add(group_h[1])
        .add(xr.sub(center).cross(s.receiver.linear_momentum_world()))
        .add(group_positions[1].sub(center).cross(group_p[1]));
    let omega = multiply(inverse(combined_tensor), h);
    let velocity = p.scale(1.0 / combined_mass);
    let expected_loss: f64 = [
        (
            mr,
            xr,
            s.receiver.linear_momentum_world(),
            world_spin(s.receiver),
            receiver_tensor,
        ),
        (
            group_mass[1],
            group_positions[1],
            group_p[1],
            group_h[1],
            group_tensor[1],
        ),
    ]
    .into_iter()
    .map(|(m, x, p, spin, inertia)| {
        let dv = p
            .scale(1.0 / m)
            .sub(velocity.add(omega.cross(x.sub(center))));
        let dw = multiply(inverse(inertia), spin).sub(omega);
        0.5 * (m * dv.dot(dv) + dw.dot(multiply(inertia, dw)))
    })
    .sum();
    vector(out.receiver.state.pose().position_world(), center, 2e-13);
    vector(out.receiver.state.linear_momentum_world(), p, 2e-13);
    vector(world_spin(out.receiver.state), h, 2e-12);
    matrix(
        tensor(
            out.receiver.assembly.properties(),
            out.receiver.state.pose().orientation(),
        ),
        combined_tensor,
        2e-12,
    );
    close(out.supplied_work_j, expected_work, 2e-13);
    close(out.receiver.removed_relative_energy_j, expected_loss, 2e-12);
    conserve(
        &[
            (s.donor, s.assembly.properties()),
            (s.receiver, s.receiver_properties),
        ],
        &out,
    );

    // Exact native-chain parity and subsequent material expansion ensure the
    // payload appears only in the receiver, with every original leaf preserved.
    let released = s
        .plan
        .partition()
        .release(s.donor, &s.impulses, c.assembly, || false)
        .unwrap();
    let captured = RigidAssembly::latch(
        &[s.receiver, released.bodies[1]],
        &[
            s.receiver_properties,
            s.plan.partition().groups()[1].properties(),
        ],
        s.receiver.pose().orientation(),
        c.assembly,
        || false,
    )
    .unwrap();
    assert_eq!(out.donor_state, released.bodies[0]);
    assert_eq!(out.receiver, captured);
    assert_eq!(out.supplied_work_j, released.supplied_work_j);
    assert_eq!(
        s.plan.partition().original_part_indices(0),
        Some([0, 2].as_slice())
    );
    assert_eq!(
        s.plan.partition().original_part_indices(1),
        Some([1, 3].as_slice())
    );
    let receiving_parts = out
        .receiver
        .assembly
        .release(out.receiver.state, &[], c.assembly, || false)
        .unwrap();
    let original_parts = s
        .plan
        .partition()
        .expand(
            &[out.donor_state, receiving_parts.bodies[1]],
            c.assembly,
            || false,
        )
        .unwrap();
    let mut all: Vec<_> = original_parts
        .into_iter()
        .zip(s.assembly.parts().iter().map(|p| p.properties))
        .collect();
    all.push((receiving_parts.bodies[0], s.receiver_properties));
    assert_eq!(all.len(), s.assembly.parts().len() + 1);
    let expanded = totals(&all);
    let accepted = totals(&[
        (out.donor_state, out.donor_properties),
        (out.receiver.state, out.receiver.assembly.properties()),
    ]);
    close(expanded.0, accepted.0, 2e-13);
    vector(expanded.1, accepted.1, 2e-12);
    vector(expanded.2, accepted.2, 2e-11);
    close(expanded.3, accepted.3, 2e-11);
}

#[test]
fn moving_observer_rotation_translation_and_common_boost_preserve_handoff() {
    let s = scene();
    let c = controls();
    let original = s
        .plan
        .transfer(
            &s.donor,
            &s.receiver,
            s.receiver_properties,
            &s.impulses,
            c,
            || false,
        )
        .unwrap();
    let rotation = rotation(Vec3::new(-0.5, 0.3, 1.0), 0.9);
    let translation = Vec3::new(14.0, -9.0, 3.0);
    let boost = Vec3::new(1024.0, -512.0, 256.0);
    let transform = |body: RigidBodyState, mass: f64| {
        RigidBodyState::new(
            Pose::new(
                rotation
                    .rotate_body_to_world(body.pose().position_world())
                    .add(translation),
                compose(rotation, body.pose().orientation()),
            )
            .unwrap(),
            rotation
                .rotate_body_to_world(body.linear_momentum_world())
                .add(boost.scale(mass)),
            body.angular_momentum_body(),
        )
        .unwrap()
    };
    let donor = transform(s.donor, s.assembly.properties().mass());
    let receiver = transform(s.receiver, s.receiver_properties.mass());
    let transformed = s
        .plan
        .transfer(
            &donor,
            &receiver,
            s.receiver_properties,
            &s.impulses,
            c,
            || false,
        )
        .unwrap();
    for (a, b, mass) in [
        (
            original.donor_state,
            transformed.donor_state,
            original.donor_properties.mass(),
        ),
        (
            original.receiver.state,
            transformed.receiver.state,
            original.receiver.assembly.properties().mass(),
        ),
    ] {
        let expected = transform(a, mass);
        vector(
            b.pose().position_world(),
            expected.pose().position_world(),
            2e-13,
        );
        vector(
            b.linear_momentum_world(),
            expected.linear_momentum_world(),
            4e-12,
        );
        vector(
            world_spin(b),
            rotation.rotate_body_to_world(world_spin(a)),
            4e-11,
        );
    }
    close(original.supplied_work_j, transformed.supplied_work_j, 2e-13);
    close(
        original.receiver.removed_relative_energy_j,
        transformed.receiver.removed_relative_energy_j,
        4e-11,
    );
    close(transformed.mass_residual_kg, 0.0, 2e-13);
}

#[test]
fn preparation_and_event_caps_keep_positive_finite_inventory() {
    let s = scene();
    let c = controls();
    for indices in [&[][..], &[0, 1, 2, 3][..], &[1, 1][..], &[7][..]] {
        assert!(matches!(
            PreparedRigidPayloadTransfer::new(
                &s.assembly,
                indices,
                Vec3::ZERO,
                Vec3::ZERO,
                c.assembly,
                || false
            ),
            Err(AssemblyError::Invalid(_))
        ));
    }
    let mut tight = c;
    tight.assembly.maximum_parts = s.assembly.parts().len();
    assert!(matches!(
        PreparedRigidPayloadTransfer::new(
            &s.assembly,
            &[1],
            Vec3::ZERO,
            Vec3::ZERO,
            tight.assembly,
            || false
        ),
        Err(AssemblyError::Budget {
            requested: 5,
            maximum: 4
        })
    ));
    assert!(matches!(
        s.plan.transfer(
            &s.donor,
            &s.receiver,
            s.receiver_properties,
            &s.impulses,
            tight,
            || false
        ),
        Err(AssemblyError::Budget {
            requested: 5,
            maximum: 4
        })
    ));
    tight = c;
    tight.assembly.maximum_parts = 5;
    tight.assembly.maximum_impulses = 3;
    s.plan
        .transfer(
            &s.donor,
            &s.receiver,
            s.receiver_properties,
            &s.impulses,
            tight,
            || false,
        )
        .unwrap();
    tight.assembly.maximum_impulses = 2;
    assert!(matches!(
        s.plan.transfer(
            &s.donor,
            &s.receiver,
            s.receiver_properties,
            &s.impulses,
            tight,
            || false
        ),
        Err(AssemblyError::Budget {
            requested: 3,
            maximum: 2
        })
    ));
    let bad = SeparationImpulse {
        part_a: 1,
        part_b: 3,
        ..s.impulses[0]
    };
    assert!(matches!(
        s.plan.transfer(
            &s.donor,
            &s.receiver,
            s.receiver_properties,
            &[bad],
            c,
            || false
        ),
        Err(AssemblyError::Invalid(
            "separation impulse must cross retained groups"
        ))
    ));
    for invalid in [f64::NAN, f64::INFINITY, -1.0] {
        let mut bad = c;
        bad.mass_tolerance_kg = invalid;
        assert!(
            s.plan
                .transfer(
                    &s.donor,
                    &s.receiver,
                    s.receiver_properties,
                    &[],
                    bad,
                    || false
                )
                .is_err()
        );
        bad = c;
        bad.attachment_tolerance_m = invalid;
        assert!(
            s.plan
                .transfer(
                    &s.donor,
                    &s.receiver,
                    s.receiver_properties,
                    &[],
                    bad,
                    || false
                )
                .is_err()
        );
    }
    assert!(
        s.plan
            .partition()
            .groups()
            .iter()
            .all(|g| g.properties().mass() > 0.0)
    );
}

#[test]
fn local_geometry_cannot_be_hidden_by_large_world_position_rounding() {
    let c = controls();
    let p = diagonal(1.0, Vec3::new(1.0, 1.0, 1.0));
    let parts = [-1.0, 1.0].map(|x| AssemblyPart {
        properties: p,
        pose_in_reference: Pose::new(Vec3::new(x, 0.0, 0.0), UnitQuaternion::IDENTITY).unwrap(),
    });
    let assembly = RigidAssembly::new(&parts, c.assembly, || false).unwrap();
    let x = Vec3::new(1e16, 0.0, 0.0);
    let donor = state(
        assembly.properties(),
        x,
        UnitQuaternion::IDENTITY,
        Vec3::ZERO,
        Vec3::ZERO,
    );
    let receiver = state(p, x, UnitQuaternion::IDENTITY, Vec3::ZERO, Vec3::ZERO);
    let before = (donor, receiver);
    // The original endpoints coincide. Release rounds the payload COM, making
    // its retained local attachment one metre away even though the two rounded
    // world attachment coordinates still compare equal. Refuse after release.
    let plan = PreparedRigidPayloadTransfer::new(
        &assembly,
        &[1],
        Vec3::ZERO,
        Vec3::ZERO,
        c.assembly,
        || false,
    )
    .unwrap();
    assert_eq!(
        plan.transfer(&donor, &receiver, p, &[], c, || false),
        Err(AssemblyError::Invalid(
            "unrepresented rigid transfer attachment coincidence"
        ))
    );
    let wrong = PreparedRigidPayloadTransfer::new(
        &assembly,
        &[1],
        Vec3::new(0.25, 0.0, 0.0),
        Vec3::ZERO,
        c.assembly,
        || false,
    )
    .unwrap();
    assert_eq!(
        wrong.transfer(&donor, &receiver, p, &[], c, || false),
        Err(AssemblyError::Invalid(
            "noncoincident rigid transfer attachment points"
        ))
    );
    assert_eq!((donor, receiver), before);
}

#[test]
fn cancellation_through_final_publication_and_late_balance_refusal_are_atomic() {
    let s = scene();
    let c = controls();
    let snapshot = (s.donor, s.receiver, s.plan.clone(), s.assembly.clone());
    let mut polls = 0;
    let expected = s
        .plan
        .transfer(
            &s.donor,
            &s.receiver,
            s.receiver_properties,
            &s.impulses,
            c,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert!(polls > 1);
    for at in 1..=polls {
        let mut seen = 0;
        assert_eq!(
            s.plan.transfer(
                &s.donor,
                &s.receiver,
                s.receiver_properties,
                &s.impulses,
                c,
                || {
                    seen += 1;
                    seen == at
                }
            ),
            Err(AssemblyError::Cancelled)
        );
        assert_eq!(
            (s.donor, s.receiver, s.plan.clone(), s.assembly.clone()),
            snapshot
        );
    }
    let mut strict = c;
    strict.assembly.linear_momentum_tolerance_ns = 0.0;
    strict.assembly.angular_momentum_tolerance_nms = 0.0;
    strict.assembly.energy_tolerance_j = 0.0;
    assert!(matches!(
        s.plan.transfer(
            &s.donor,
            &s.receiver,
            s.receiver_properties,
            &s.impulses,
            strict,
            || false
        ),
        Err(AssemblyError::Balance(_))
    ));
    assert_eq!(
        s.plan
            .transfer(
                &s.donor,
                &s.receiver,
                s.receiver_properties,
                &s.impulses,
                c,
                || false
            )
            .unwrap(),
        expected
    );
    let mut preparation_polls = 0;
    PreparedRigidPayloadTransfer::new(
        &s.assembly,
        &[1, 3],
        Vec3::ZERO,
        Vec3::ZERO,
        c.assembly,
        || {
            preparation_polls += 1;
            false
        },
    )
    .unwrap();
    for at in 1..=preparation_polls {
        let mut seen = 0;
        assert_eq!(
            PreparedRigidPayloadTransfer::new(
                &s.assembly,
                &[1, 3],
                Vec3::ZERO,
                Vec3::ZERO,
                c.assembly,
                || {
                    seen += 1;
                    seen == at
                }
            ),
            Err(AssemblyError::Cancelled)
        );
    }
    assert_eq!(
        (s.donor, s.receiver, s.plan.clone(), s.assembly.clone()),
        snapshot
    );
}
