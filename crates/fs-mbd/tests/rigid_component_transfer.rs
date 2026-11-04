//! Full-inertia transfers: independent mechanical and inventory oracles.
use fs_mbd::assembly::transfer::rigid::{
    RigidTransferControls, RigidTransferInput, RigidTransferResult, transfer_rigid_components,
};
use fs_mbd::assembly::{
    AssemblyControls, AssemblyError, AssemblyPart, RigidAssembly, SeparationImpulse,
};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};

fn controls() -> RigidTransferControls {
    RigidTransferControls {
        assembly: AssemblyControls {
            maximum_parts: 64,
            maximum_impulses: 16,
            linear_momentum_tolerance_ns: 2e-8,
            angular_momentum_tolerance_nms: 2e-8,
            energy_tolerance_j: 2e-8,
        },
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
fn compose(q: UnitQuaternion, r: UnitQuaternion) -> UnitQuaternion {
    let [a, b, c, d] = q.components();
    let [e, f, g, h] = r.components();
    UnitQuaternion::new(
        a * e - b * f - c * g - d * h,
        a * f + b * e + c * h - d * g,
        a * g - b * h + c * e + d * f,
        a * h + b * g - c * f + d * e,
    )
    .unwrap()
}
fn state(
    assembly: &RigidAssembly,
    position: Vec3,
    axes: UnitQuaternion,
    velocity: Vec3,
    angular_velocity_body: Vec3,
) -> RigidBodyState {
    let p = assembly.properties();
    RigidBodyState::new(
        Pose::new(position, axes).unwrap(),
        velocity.scale(p.mass()),
        p.angular_momentum_body_checked(angular_velocity_body)
            .unwrap(),
    )
    .unwrap()
}
fn spin(body: RigidBodyState) -> Vec3 {
    body.pose()
        .orientation()
        .rotate_body_to_world(body.angular_momentum_body())
}
fn basis() -> [Vec3; 3] {
    [
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    ]
}
fn same_pose(a: Pose, b: Pose, tolerance: f64) {
    vector(a.position_world(), b.position_world(), tolerance);
    for axis in basis() {
        vector(
            a.orientation().rotate_body_to_world(axis),
            b.orientation().rotate_body_to_world(axis),
            tolerance,
        );
    }
}
#[derive(Clone)]
struct Fixture {
    assemblies: Vec<RigidAssembly>,
    states: Vec<RigidBodyState>,
    destinations: Vec<Vec<usize>>,
    impulses: Vec<Vec<SeparationImpulse>>,
}
fn run(
    f: &Fixture,
    controls: RigidTransferControls,
    cancelled: impl FnMut() -> bool,
) -> Result<RigidTransferResult, AssemblyError> {
    let inputs: Vec<_> = (0..f.assemblies.len())
        .map(|i| RigidTransferInput {
            assembly: &f.assemblies[i],
            state: f.states[i],
            destinations: &f.destinations[i],
            impulses: &f.impulses[i],
        })
        .collect();
    transfer_rigid_components(&inputs, controls, cancelled)
}
fn original_pose(f: &Fixture, body: usize, part: usize) -> Pose {
    let reference = f.states[body].pose();
    let local = f.assemblies[body].parts()[part].pose_in_reference;
    Pose::new(
        reference
            .point_world_from_body(local.position_world())
            .unwrap(),
        compose(reference.orientation(), local.orientation()),
    )
    .unwrap()
}
fn inventory(f: &Fixture, result: &RigidTransferResult, tolerance: f64) {
    let mut seen: Vec<_> = f
        .assemblies
        .iter()
        .map(|a| vec![false; a.parts().len()])
        .collect();
    for (destination, endpoint) in result.bodies.iter().enumerate() {
        assert_eq!(endpoint.assembly.parts().len(), endpoint.origins.len());
        for (part, origin) in endpoint.assembly.parts().iter().zip(&endpoint.origins) {
            assert!(!seen[origin.body][origin.part], "duplicated material");
            seen[origin.body][origin.part] = true;
            assert_eq!(f.destinations[origin.body][origin.part], destination);
            assert_eq!(
                part.properties,
                f.assemblies[origin.body].parts()[origin.part].properties
            );
            let pose = Pose::new(
                endpoint
                    .state
                    .pose()
                    .point_world_from_body(part.pose_in_reference.position_world())
                    .unwrap(),
                compose(
                    endpoint.state.pose().orientation(),
                    part.pose_in_reference.orientation(),
                ),
            )
            .unwrap();
            same_pose(pose, original_pose(f, origin.body, origin.part), tolerance);
        }
    }
    assert!(seen.into_iter().flatten().all(|present| present));
}
fn simple(boost: Vec3) -> Fixture {
    let unit = MassProperties::new(1.0, Vec3::ZERO, Vec3::new(0.25, 0.25, 0.25)).unwrap();
    let parts: Vec<_> = [(-1.0, -1.0), (-1.0, 1.0), (1.0, -1.0), (1.0, 1.0)]
        .into_iter()
        .map(|(x, y)| AssemblyPart {
            properties: unit,
            pose_in_reference: Pose::new(Vec3::new(x, y, 0.0), UnitQuaternion::IDENTITY).unwrap(),
        })
        .collect();
    let donor = RigidAssembly::new(&parts, controls().assembly, || false).unwrap();
    let receiver = RigidAssembly::new(
        &[AssemblyPart {
            properties: MassProperties::new(2.0, Vec3::ZERO, Vec3::new(0.5, 0.5, 0.5)).unwrap(),
            pose_in_reference: Pose::new(Vec3::ZERO, UnitQuaternion::IDENTITY).unwrap(),
        }],
        controls().assembly,
        || false,
    )
    .unwrap();
    let states = vec![
        state(
            &donor,
            Vec3::ZERO,
            UnitQuaternion::IDENTITY,
            boost,
            Vec3::ZERO,
        ),
        state(
            &receiver,
            Vec3::new(3.0, 0.0, 0.0),
            UnitQuaternion::IDENTITY,
            boost.add(Vec3::new(0.0, -0.5, 0.0)),
            Vec3::new(0.0, 0.0, 2.0),
        ),
    ];
    Fixture {
        assemblies: vec![donor, receiver],
        states,
        destinations: vec![vec![0, 0, 1, 1], vec![1]],
        impulses: vec![
            vec![
                SeparationImpulse {
                    part_a: 2,
                    part_b: 0,
                    point_assembly_m: Vec3::ZERO,
                    impulse_on_a_assembly_ns: Vec3::new(0.0, 1.0, 0.0),
                    couple_on_a_assembly_nms: Vec3::new(0.0, 0.0, 0.2),
                },
                SeparationImpulse {
                    part_a: 3,
                    part_b: 1,
                    point_assembly_m: Vec3::ZERO,
                    impulse_on_a_assembly_ns: Vec3::new(0.0, 2.0, 0.0),
                    couple_on_a_assembly_nms: Vec3::new(0.0, 0.0, 0.3),
                },
            ],
            vec![],
        ],
    }
}

#[test]
fn g1_shared_group_recoil_and_intrinsic_spin_have_closed_form_work_and_loss() {
    for boost in [Vec3::ZERO, Vec3::new(1e9, -2e9, 3e9)] {
        let f = simple(boost);
        let result = run(&f, controls(), || false).unwrap();
        let donor = &result.bodies[0];
        let receiver = &result.bodies[1];
        // Both separated groups have m=2 and Iz=5/2. The two impulses
        // combine to J_y=3, K_out=-5/2 and K_retained=-7/2.
        close(result.supplied_work_j, 41.0 / 5.0, 2e-12);
        assert!(
            result.supplied_work_j > 4.552 + 3.0,
            "shared-work cross terms"
        );
        close(donor.supplied_work_j, 41.0 / 5.0, 2e-12);
        close(receiver.supplied_work_j, 0.0, 0.0);
        close(donor.assembly.properties().mass(), 2.0, 0.0);
        vector(
            donor.state.pose().position_world(),
            Vec3::new(-1.0, 0.0, 0.0),
            1e-13,
        );
        vector(
            donor.state.linear_momentum_world(),
            boost.scale(2.0).add(Vec3::new(0.0, -3.0, 0.0)),
            1e-12,
        );
        vector(spin(donor.state), Vec3::new(0.0, 0.0, -3.5), 2e-12);
        // The incoming group's spin (-5/2) and both orbital contributions
        // combine with the receiver's original spin (+1) to Lz=-11/2.
        close(receiver.assembly.properties().mass(), 4.0, 0.0);
        close(
            receiver.assembly.properties().inertia_tensor_body().m[8],
            7.0,
            2e-12,
        );
        vector(
            receiver.state.pose().position_world(),
            Vec3::new(2.0, 0.0, 0.0),
            1e-13,
        );
        vector(
            receiver.state.linear_momentum_world(),
            boost.scale(4.0).add(Vec3::new(0.0, 2.0, 0.0)),
            1e-12,
        );
        vector(spin(receiver.state), Vec3::new(0.0, 0.0, -5.5), 2e-12);
        close(receiver.removed_relative_energy_j, 117.0 / 56.0, 2e-11);
        close(result.removed_relative_energy_j, 117.0 / 56.0, 2e-11);
        close(result.mass_residual_kg, 0.0, 0.0);
        vector(
            result.balance.linear_momentum_residual_world_ns,
            Vec3::ZERO,
            1e-12,
        );
        vector(
            result.balance.angular_momentum_residual_world_nms,
            Vec3::ZERO,
            2e-10,
        );
        close(result.balance.energy_residual_j, 0.0, 2e-10);
        inventory(&f, &result, 2e-12);
    }
}

#[test]
fn g1_retained_packets_do_not_acquire_artificial_leaf_latching_loss() {
    let mut f = simple(Vec3::ZERO);
    // The original receiver follows the retained packet's FINAL rigid field.
    // Its material moves to body 0; body 1 receives only the outgoing packet.
    // Neither final connection removes any relative motion.
    f.destinations[1][0] = 0;
    f.states[1] = state(
        &f.assemblies[1],
        Vec3::new(3.0, 0.0, 0.0),
        UnitQuaternion::IDENTITY,
        Vec3::new(0.0, -7.1, 0.0),
        Vec3::new(0.0, 0.0, -1.4),
    );
    let result = run(&f, controls(), || false).unwrap();
    close(result.supplied_work_j, 8.2, 2e-12);
    close(result.removed_relative_energy_j, 0.0, 2e-25);
    let retained = &result.bodies[0];
    let outgoing = &result.bodies[1];
    vector(
        retained
            .state
            .center_of_mass_velocity_world(retained.assembly.properties())
            .unwrap(),
        Vec3::new(0.0, -4.3, 0.0),
        2e-12,
    );
    vector(
        retained
            .assembly
            .properties()
            .angular_velocity_body_checked(retained.state.angular_momentum_body())
            .unwrap(),
        Vec3::new(0.0, 0.0, -1.4),
        2e-12,
    );
    vector(
        outgoing
            .state
            .center_of_mass_velocity_world(outgoing.assembly.properties())
            .unwrap(),
        Vec3::new(0.0, 1.5, 0.0),
        2e-12,
    );
    vector(
        outgoing
            .assembly
            .properties()
            .angular_velocity_body_checked(outgoing.state.angular_momentum_body())
            .unwrap(),
        Vec3::new(0.0, 0.0, -1.0),
        2e-12,
    );
    inventory(&f, &result, 2e-12);
}

// Independent dense 3x3 mechanics on original material leaves. These helpers
// never call native assembly partition, release, latch or transfer routines.
fn mul(a: [f64; 9], v: Vec3) -> Vec3 {
    Vec3::new(
        a[0] * v.x + a[1] * v.y + a[2] * v.z,
        a[3] * v.x + a[4] * v.y + a[5] * v.z,
        a[6] * v.x + a[7] * v.y + a[8] * v.z,
    )
}
fn determinant(a: [f64; 9]) -> f64 {
    a[0] * (a[4] * a[8] - a[5] * a[7]) - a[1] * (a[3] * a[8] - a[5] * a[6])
        + a[2] * (a[3] * a[7] - a[4] * a[6])
}
fn solve(a: [f64; 9], v: Vec3) -> Vec3 {
    let rhs = [v.x, v.y, v.z];
    let d = determinant(a);
    assert!(d > 0.0);
    let mut answer = [0.0; 3];
    for column in 0..3 {
        let mut replaced = a;
        for row in 0..3 {
            replaced[3 * row + column] = rhs[row];
        }
        answer[column] = determinant(replaced) / d;
    }
    Vec3::new(answer[0], answer[1], answer[2])
}
fn rotated_tensor(a: [f64; 9], q: UnitQuaternion) -> [f64; 9] {
    let mut result = [0.0; 9];
    for (column, axis) in basis().into_iter().enumerate() {
        let image = q.rotate_body_to_world(mul(a, q.rotate_world_to_body(axis)));
        for (row, value) in [image.x, image.y, image.z].into_iter().enumerate() {
            result[3 * row + column] = value;
        }
    }
    result
}
struct Leaf {
    mass: f64,
    position: Vec3,
    velocity: Vec3,
    omega: Vec3,
    inertia_world: [f64; 9],
}
fn leaves(f: &Fixture, destination: usize) -> Vec<Leaf> {
    let mut result = Vec::new();
    for body in 0..f.assemblies.len() {
        let native = f.states[body];
        let assembly = &f.assemblies[body];
        let velocity = native
            .linear_momentum_world()
            .scale(1.0 / assembly.properties().mass());
        let omega = native.pose().orientation().rotate_body_to_world(solve(
            assembly.properties().inertia_tensor_body().m,
            native.angular_momentum_body(),
        ));
        for (part, original) in assembly.parts().iter().enumerate() {
            if f.destinations[body][part] != destination {
                continue;
            }
            let pose = original_pose(f, body, part);
            let arm = pose.position_world().sub(native.pose().position_world());
            result.push(Leaf {
                mass: original.properties.mass(),
                position: pose.position_world(),
                velocity: velocity.add(omega.cross(arm)),
                omega,
                inertia_world: rotated_tensor(
                    original.properties.inertia_tensor_body().m,
                    pose.orientation(),
                ),
            });
        }
    }
    result
}
fn independent_passive_oracle(f: &Fixture, result: &RigidTransferResult) {
    let mut loss_sum = 0.0;
    for (target, endpoint) in result.bodies.iter().enumerate() {
        let leaves = leaves(f, target);
        let mass: f64 = leaves.iter().map(|p| p.mass).sum();
        let center = leaves.iter().fold(Vec3::ZERO, |sum, p| {
            sum.add(p.position.scale(p.mass / mass))
        });
        let linear = leaves
            .iter()
            .fold(Vec3::ZERO, |sum, p| sum.add(p.velocity.scale(p.mass)));
        let mut angular = Vec3::ZERO;
        let mut inertia = [0.0; 9];
        for leaf in &leaves {
            let arm = leaf.position.sub(center);
            let r = [arm.x, arm.y, arm.z];
            angular = angular
                .add(mul(leaf.inertia_world, leaf.omega))
                .add(arm.cross(leaf.velocity.scale(leaf.mass)));
            for row in 0..3 {
                for col in 0..3 {
                    inertia[3 * row + col] += leaf.inertia_world[3 * row + col]
                        + leaf.mass
                            * (if row == col { arm.dot(arm) } else { 0.0 } - r[row] * r[col]);
                }
            }
        }
        let velocity = linear.scale(1.0 / mass);
        let omega = solve(inertia, angular);
        let loss: f64 = leaves
            .iter()
            .map(|leaf| {
                let dv = leaf
                    .velocity
                    .sub(velocity.add(omega.cross(leaf.position.sub(center))));
                let dw = leaf.omega.sub(omega);
                0.5 * (leaf.mass * dv.dot(dv) + dw.dot(mul(leaf.inertia_world, dw)))
            })
            .sum();
        loss_sum += loss;
        close(endpoint.assembly.properties().mass(), mass, 2e-12);
        vector(endpoint.state.pose().position_world(), center, 2e-12);
        vector(endpoint.state.linear_momentum_world(), linear, 2e-11);
        vector(spin(endpoint.state), angular, 2e-10);
        let actual_inertia = rotated_tensor(
            endpoint.assembly.properties().inertia_tensor_body().m,
            endpoint.state.pose().orientation(),
        );
        for (actual, expected) in actual_inertia.into_iter().zip(inertia) {
            close(actual, expected, 2e-10);
        }
        vector(
            endpoint.state.pose().orientation().rotate_body_to_world(
                endpoint
                    .assembly
                    .properties()
                    .angular_velocity_body_checked(endpoint.state.angular_momentum_body())
                    .unwrap(),
            ),
            omega,
            2e-11,
        );
        close(endpoint.removed_relative_energy_j, loss, 2e-10);
    }
    close(result.supplied_work_j, 0.0, 0.0);
    close(result.removed_relative_energy_j, loss_sum, 4e-10);
}
fn tensor_cycle() -> Fixture {
    let mut assemblies = Vec::new();
    let mut states = Vec::new();
    for body in 0..3 {
        let b = body as f64;
        let parts: Vec<_> = (0..3)
            .map(|part| {
                let p = part as f64;
                let mass = 1.0 + b + 0.5 * p;
                AssemblyPart {
                    properties: MassProperties::from_principal_frame(
                        mass,
                        Vec3::ZERO,
                        Vec3::new(0.4, 0.7, 0.9).scale(mass),
                        rotation(Vec3::new(0.3, 0.8, -0.4), 0.13 + 0.2 * b + 0.17 * p),
                    )
                    .unwrap(),
                    pose_in_reference: Pose::new(
                        Vec3::new(0.4 * p - 0.2, 0.3 * p * p - 0.4, 0.2 - 0.13 * p),
                        rotation(Vec3::new(-0.5, 0.4, 0.7), -0.1 + 0.3 * p - 0.2 * b),
                    )
                    .unwrap(),
                }
            })
            .collect();
        let assembly = RigidAssembly::new(&parts, controls().assembly, || false).unwrap();
        states.push(state(
            &assembly,
            Vec3::new(1.5 * b, 0.4 - 0.3 * b, -0.2 * b),
            rotation(Vec3::new(0.4, -0.7, 0.2), 0.27 + 0.31 * b),
            Vec3::new(0.8 - 0.3 * b, -0.2 + 0.4 * b, 0.1 - 0.15 * b),
            Vec3::new(0.3 + 0.2 * b, -0.4 + 0.1 * b, 0.8 - 0.3 * b),
        ));
        assemblies.push(assembly);
    }
    Fixture {
        assemblies,
        states,
        destinations: (0..3).map(|i| vec![i, (i + 1) % 3, (i + 1) % 3]).collect(),
        impulses: vec![vec![]; 3],
    }
}
fn world_totals(
    states: impl IntoIterator<Item = (RigidBodyState, MassProperties)>,
) -> (Vec3, Vec3, f64) {
    let observer = Vec3::new(-0.6, 1.2, -0.4);
    let mut linear = Vec3::ZERO;
    let mut angular = Vec3::ZERO;
    let mut energy = 0.0;
    for (state, properties) in states {
        let p = state.linear_momentum_world();
        let l = state.angular_momentum_body();
        linear = linear.add(p);
        angular = angular
            .add(spin(state))
            .add(state.pose().position_world().sub(observer).cross(p));
        energy += 0.5
            * (p.dot(p) / properties.mass() + l.dot(solve(properties.inertia_tensor_body().m, l)));
    }
    (linear, angular, energy)
}

#[test]
fn g1_g3_full_tensor_cycle_matches_leaf_mechanics_and_scene_covariance() {
    let f = tensor_cycle();
    let result = run(&f, controls(), || false).unwrap();
    independent_passive_oracle(&f, &result);
    inventory(&f, &result, 2e-12);
    let before = world_totals(
        f.states
            .iter()
            .copied()
            .zip(f.assemblies.iter().map(|a| a.properties())),
    );
    let after = world_totals(
        result
            .bodies
            .iter()
            .map(|b| (b.state, b.assembly.properties())),
    );
    vector(before.0, after.0, 2e-11);
    vector(before.1, after.1, 2e-10);
    close(
        after.2 - before.2 + result.removed_relative_energy_j,
        0.0,
        4e-10,
    );
    assert!(result.removed_relative_energy_j > 0.1);
    assert!(
        result
            .bodies
            .iter()
            .any(|b| b.assembly.properties().inertia_tensor_body().m[1].abs() > 0.1)
    );

    let q = rotation(Vec3::new(0.3, -0.4, 0.8), 0.61);
    let shift = Vec3::new(13.0, -7.0, 4.0);
    let boost = Vec3::new(2e4, -3e4, 4e4);
    let mut moved = f.clone();
    for (i, old) in f.states.iter().enumerate() {
        moved.states[i] = RigidBodyState::new(
            Pose::new(
                q.rotate_body_to_world(old.pose().position_world())
                    .add(shift),
                compose(q, old.pose().orientation()),
            )
            .unwrap(),
            q.rotate_body_to_world(old.linear_momentum_world())
                .add(boost.scale(f.assemblies[i].properties().mass())),
            old.angular_momentum_body(),
        )
        .unwrap();
    }
    let mut loose = controls();
    loose.assembly.linear_momentum_tolerance_ns = 2e-6;
    loose.assembly.angular_momentum_tolerance_nms = 2e-6;
    loose.assembly.energy_tolerance_j = 2e-6;
    let transformed = run(&moved, loose, || false).unwrap();
    for (a, b) in result.bodies.iter().zip(&transformed.bodies) {
        vector(
            b.state.pose().position_world(),
            q.rotate_body_to_world(a.state.pose().position_world())
                .add(shift),
            2e-11,
        );
        vector(
            b.state.linear_momentum_world(),
            q.rotate_body_to_world(a.state.linear_momentum_world())
                .add(boost.scale(a.assembly.properties().mass())),
            2e-9,
        );
        vector(spin(b.state), q.rotate_body_to_world(spin(a.state)), 2e-8);
        close(
            b.removed_relative_energy_j,
            a.removed_relative_energy_j,
            2e-8,
        );
    }
    inventory(&moved, &transformed, 2e-11);

    // Permute physical body identities and remap every destination, then
    // compare the corresponding actual endpoints, not only their receipts.
    let order = [2, 0, 1];
    let inverse = [1, 2, 0];
    let permuted = Fixture {
        assemblies: order.iter().map(|&i| f.assemblies[i].clone()).collect(),
        states: order.iter().map(|&i| f.states[i]).collect(),
        destinations: order
            .iter()
            .map(|&i| f.destinations[i].iter().map(|&j| inverse[j]).collect())
            .collect(),
        impulses: vec![vec![]; 3],
    };
    let p = run(&permuted, controls(), || false).unwrap();
    for (index, &original) in order.iter().enumerate() {
        let a = &result.bodies[original];
        let b = &p.bodies[index];
        same_pose(a.state.pose(), b.state.pose(), 2e-12);
        vector(
            a.state.linear_momentum_world(),
            b.state.linear_momentum_world(),
            2e-11,
        );
        vector(spin(a.state), spin(b.state), 2e-10);
        close(
            a.removed_relative_energy_j,
            b.removed_relative_energy_j,
            2e-10,
        );
    }
    inventory(&permuted, &p, 2e-12);
}

#[test]
fn g3_exact_noop_whole_inventory_cycles_and_next_event_preserve_material_identity() {
    let mut identity = tensor_cycle();
    for (i, row) in identity.destinations.iter_mut().enumerate() {
        row.fill(i);
    }
    let unchanged = run(&identity, controls(), || false).unwrap();
    for (i, endpoint) in unchanged.bodies.iter().enumerate() {
        assert_eq!(endpoint.assembly, identity.assemblies[i]);
        assert_eq!(endpoint.state, identity.states[i]);
        assert_eq!(endpoint.supplied_work_j, 0.0);
        assert_eq!(endpoint.removed_relative_energy_j, 0.0);
    }
    inventory(&identity, &unchanged, 0.0);
    assert_eq!(unchanged.mass_residual_kg, 0.0);
    assert_eq!(unchanged.balance.energy_residual_j, 0.0);

    // Every donor gives away its entire inventory; original replacement
    // material arrives from another donor in the SAME closed cycle.
    let mut f = tensor_cycle();
    for (i, row) in f.destinations.iter_mut().enumerate() {
        row.fill((i + 1) % 3);
    }
    let transferred = run(&f, controls(), || false).unwrap();
    close(transferred.supplied_work_j, 0.0, 0.0);
    close(transferred.removed_relative_energy_j, 0.0, 2e-25);
    inventory(&f, &transferred, 2e-12);
    for (destination, &source) in [2, 0, 1].iter().enumerate() {
        let endpoint = &transferred.bodies[destination];
        vector(
            endpoint.state.pose().position_world(),
            f.states[source].pose().position_world(),
            2e-12,
        );
        vector(
            endpoint.state.linear_momentum_world(),
            f.states[source].linear_momentum_world(),
            2e-11,
        );
        vector(spin(endpoint.state), spin(f.states[source]), 2e-11);
        assert_eq!(
            endpoint.assembly.parts().len(),
            f.assemblies[source].parts().len()
        );
        close(
            endpoint.assembly.properties().mass(),
            f.assemblies[source].properties().mass(),
            0.0,
        );
    }

    // Reuse the published flat leaf inventory. No expansion, external parcel
    // cache or hidden original assembly is needed to transfer it back.
    let second = Fixture {
        assemblies: transferred
            .bodies
            .iter()
            .map(|b| b.assembly.clone())
            .collect(),
        states: transferred.bodies.iter().map(|b| b.state).collect(),
        destinations: transferred
            .bodies
            .iter()
            .map(|b| b.origins.iter().map(|o| o.body).collect())
            .collect(),
        impulses: vec![vec![]; 3],
    };
    let restored = run(&second, controls(), || false).unwrap();
    inventory(&second, &restored, 2e-12);
    close(restored.supplied_work_j, 0.0, 0.0);
    close(restored.removed_relative_energy_j, 0.0, 2e-25);
    for (i, endpoint) in restored.bodies.iter().enumerate() {
        same_pose(endpoint.state.pose(), f.states[i].pose(), 2e-12);
        vector(
            endpoint.state.linear_momentum_world(),
            f.states[i].linear_momentum_world(),
            2e-11,
        );
        vector(spin(endpoint.state), spin(f.states[i]), 4e-11);
        for (actual, expected) in endpoint
            .assembly
            .properties()
            .inertia_tensor_body()
            .m
            .into_iter()
            .zip(f.assemblies[i].properties().inertia_tensor_body().m)
        {
            close(actual, expected, 4e-11);
        }
    }
}

#[test]
fn g4_late_inventory_refusal_caps_and_final_cancellation_leave_no_partial_event() {
    let f = simple(Vec3::ZERO);
    let states = f.states.clone();
    let inventories = f.assemblies.clone();
    let mut invalid = f.clone();
    invalid.destinations[1][0] = 2;
    assert!(matches!(
        run(&invalid, controls(), || false),
        Err(AssemblyError::Invalid(_))
    ));
    invalid = f.clone();
    invalid.destinations[1].clear();
    assert!(matches!(
        run(&invalid, controls(), || false),
        Err(AssemblyError::Invalid(_))
    ));
    invalid = f.clone();
    for row in &mut invalid.destinations {
        row.fill(0);
    }
    invalid.impulses[0].clear();
    assert!(matches!(
        run(&invalid, controls(), || false),
        Err(AssemblyError::Invalid(_))
    ));
    invalid = f.clone();
    invalid.impulses[0][0].part_a = 0;
    invalid.impulses[0][0].part_b = 1;
    assert!(matches!(
        run(&invalid, controls(), || false),
        Err(AssemblyError::Invalid(_))
    ));

    let mut limited = controls();
    limited.assembly.maximum_parts = 6; // Two bodies plus five original leaves.
    assert!(matches!(
        run(&f, limited, || false),
        Err(AssemblyError::Budget {
            requested: 7,
            maximum: 6
        })
    ));
    limited = controls();
    limited.assembly.maximum_impulses = 1;
    assert!(matches!(
        run(&f, limited, || false),
        Err(AssemblyError::Budget {
            requested: 2,
            maximum: 1
        })
    ));

    let mut polls = 0;
    let accepted = run(&f, controls(), || {
        polls += 1;
        false
    })
    .unwrap();
    assert!(polls > 1);
    let mut cancelled_polls = 0;
    assert!(matches!(
        run(&f, controls(), || {
            cancelled_polls += 1;
            cancelled_polls == polls
        }),
        Err(AssemblyError::Cancelled)
    ));
    assert_eq!(cancelled_polls, polls);
    assert_eq!(f.states, states);
    assert_eq!(f.assemblies, inventories);
    assert_eq!(run(&f, controls(), || false).unwrap(), accepted);
}
