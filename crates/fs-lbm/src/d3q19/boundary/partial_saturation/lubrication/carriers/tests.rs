//! G0/G3 public shared-carrier packing and unique physical inventories.
use super::*;

fn body(center: [f64; 3], mass: f64, velocity: [f64; 3]) -> PartialSaturationBody3 {
    PartialSaturationBody3 {
        center,
        mass,
        velocity,
        inertia_world: [
            [0.12, 0.009, -0.004],
            [0.009, 0.17, 0.006],
            [-0.004, 0.006, 0.21],
        ]
        .map(|row| row.map(|v| mass * v)),
        angular_velocity: [0.002, -0.003, 0.001],
    }
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 3e-15, "{a:.17e} != {b:.17e}");
}
fn difference(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    core::array::from_fn(|i| a[i] - b[i])
}
fn energy(b: PartialSaturationBody3) -> f64 {
    0.5 * b.mass * dot(b.velocity, b.velocity)
        + 0.5
            * dot(
                b.angular_velocity,
                b.inertia_world.map(|r| dot(r, b.angular_velocity)),
            )
}
fn momentum(b: PartialSaturationBody3, origin: [f64; 3]) -> [f64; 6] {
    let p = b.velocity.map(|u| b.mass * u);
    let orbital = cross(difference(b.center, origin), p);
    let spin = b.inertia_world.map(|r| dot(r, b.angular_velocity));
    core::array::from_fn(|i| {
        if i < 3 {
            p[i]
        } else {
            orbital[i - 3] + spin[i - 3]
        }
    })
}
fn world_load(b: PartialSaturationBody3, load: [f64; 6], origin: [f64; 3]) -> [f64; 6] {
    let orbital = cross(
        difference(b.center, origin),
        core::array::from_fn(|i| load[i]),
    );
    core::array::from_fn(|i| {
        if i < 3 {
            load[i]
        } else {
            load[i] + orbital[i - 3]
        }
    })
}
fn work(a: PartialSaturationBody3, b: PartialSaturationBody3, load: [f64; 6]) -> f64 {
    (0..3)
        .map(|i| {
            0.5 * load[i] * (a.velocity[i] + b.velocity[i])
                + 0.5 * load[i + 3] * (a.angular_velocity[i] + b.angular_velocity[i])
        })
        .sum()
}
fn face(carrier: usize, point: [f64; 3], normal: [f64; 3]) -> SharedLubricationSurface3 {
    SharedLubricationSurface3::Carried(
        CarrierLubricationFace3::new(carrier, point, normal, 0.012, 0.001)
            .unwrap()
            .with_tangential(0.01, 0.001)
            .unwrap(),
    )
}
fn model() -> SphereLubrication3 {
    SphereLubrication3::new(0.11, 0.2, 0.001, 3)
        .unwrap()
        .with_resistance_mode_limit(15)
        .unwrap()
}
fn prescribed() -> SharedLubricationSurface3 {
    SharedLubricationSurface3::Prescribed(
        LubricationPlane3::new(
            [0.0, 0.0, 4.505],
            [0.0, 0.0, -1.0],
            [0.002, -0.001, 0.0],
            0.012,
            0.001,
        )
        .unwrap()
        .with_tangential(0.01, 0.001)
        .unwrap()
        .with_axial_spin(0.003)
        .unwrap(),
    )
}

#[test]
fn two_faces_and_one_actuator_return_one_carrier_inventory_with_complete_wrench_balance() {
    let sphere = body([3.0, 4.0, 4.0], 2.0, [-0.012, 0.007, -0.003]);
    let carrier = body([2.2, 3.3, 4.7], 4.0, [-0.001, 0.003, 0.002]);
    let surfaces = [
        face(0, [2.495, 0.0, 0.0], [1.0, 0.0, 0.0]),
        face(0, [0.0, 3.495, 0.0], [0.0, 1.0, 0.0]),
        prescribed(),
    ];
    let external = [0.0004, -0.0002, 0.0001, 0.00003, -0.00001, 0.00002];
    let carrier_external = [0.0002, -0.0001, 0.0003, 0.00001, 0.00002, -0.00003];
    let origin = [-2.0, 0.7, 1.2];
    let gate = CancelGate::new();
    let native = prepare(&surfaces, &[carrier], &[carrier_external], &gate).unwrap();
    assert_eq!(native[0].recoil.unwrap().carrier_anchor, Some(0));
    assert_eq!(native[1].recoil.unwrap().carrier_anchor, Some(0));
    let mut proposed = [sphere];
    let (film, reports) = solve(
        model(),
        [0.5],
        &native,
        &[sphere],
        &[Response::new(external)],
        &mut proposed,
        origin,
        &gate,
    )
    .unwrap();
    let (faces, carriers) = collect(&surfaces, &[carrier], &reports, &gate).unwrap();
    assert_eq!(faces.len(), 3);
    assert_eq!(carriers.len(), 1);
    assert_eq!(
        faces.iter().map(|f| f.carrier).collect::<Vec<_>>(),
        [Some(0), Some(0), None]
    );
    assert_eq!(carriers[0].active_films, 2);
    assert_eq!(carriers[0].tangential_films, 2);
    assert_eq!(film.active_resistance_modes, 15);
    let after = carriers[0].motion.body_after;
    let before_s = momentum(sphere, origin);
    let after_s = momentum(proposed[0], origin);
    let before_c = momentum(carrier, origin);
    let after_c = momentum(after, origin);
    let js = world_load(sphere, external, origin);
    let jc = world_load(carrier, carrier_external, origin);
    for i in 0..6 {
        let actuator = if i < 3 {
            faces[2].impulse[i]
        } else {
            faces[2].angular_impulse[i - 3]
        };
        close(
            after_s[i] - before_s[i] + after_c[i] - before_c[i] - js[i] - jc[i] + actuator,
            0.0,
        );
        close(carriers[0].motion.momentum_residual[i], 0.0);
    }
    for i in 0..3 {
        close(
            carriers[0].impulse[i],
            faces[0].impulse[i] + faces[1].impulse[i],
        );
        close(
            carriers[0].motion.spin_impulse[i],
            reports[0].finite_body.unwrap().spin_impulse[i]
                + reports[1].finite_body.unwrap().spin_impulse[i],
        );
    }
    close(carriers[0].work, faces[0].work + faces[1].work);
    close(
        carriers[0].motion.kinetic_energy_change,
        energy(after) - energy(carrier),
    );
    close(
        carriers[0].motion.external_work,
        work(carrier, after, carrier_external),
    );
    close(carriers[0].motion.work_residual, 0.0);
    close(
        energy(proposed[0]) - energy(sphere) + energy(after)
            - energy(carrier)
            - work(sphere, proposed[0], external)
            - work(carrier, after, carrier_external)
            + faces[2].work
            + film.dissipated_energy,
        0.0,
    );
}

#[test]
fn one_face_per_carrier_and_prescribed_only_keep_legacy_arithmetic_and_caller_order() {
    let sphere = body([3.0, 4.0, 4.0], 2.0, [-0.012, 0.007, -0.003]);
    let carriers = [
        body([2.2, 3.3, 4.7], 4.0, [-0.001, 0.003, 0.002]),
        body([3.3, 3.2, 3.7], 3.0, [0.003, -0.002, 0.001]),
    ];
    let surfaces = [
        face(1, [2.495, 0.0, 0.0], [1.0, 0.0, 0.0]),
        prescribed(),
        face(0, [0.0, 3.495, 0.0], [0.0, 1.0, 0.0]),
    ];
    let gate = CancelGate::new();
    let native = prepare(&surfaces, &carriers, &[[0.0; 6]; 2], &gate).unwrap();
    let mut legacy = native.clone();
    for plane in &mut legacy {
        if let Some(input) = &mut plane.recoil {
            input.carrier_anchor = None;
        }
    }
    let mut proposed = [sphere];
    let mut old = [sphere];
    let response = [Response::default()];
    let actual = solve(
        model(),
        [0.5],
        &native,
        &[sphere],
        &response,
        &mut proposed,
        [0.0; 3],
        &gate,
    )
    .unwrap();
    let expected = solve(
        model(),
        [0.5],
        &legacy,
        &[sphere],
        &response,
        &mut old,
        [0.0; 3],
        &gate,
    )
    .unwrap();
    assert_eq!(proposed, old);
    assert_eq!(actual, expected);
    let (faces, groups) = collect(&surfaces, &carriers, &actual.1, &gate).unwrap();
    assert_eq!(faces[0].carrier, Some(1));
    assert_eq!(faces[2].carrier, Some(0));
    assert_eq!(
        groups[1].motion.body_after,
        actual.1[0].finite_body.unwrap().body_after
    );
    assert_eq!(
        groups[0].motion.body_after,
        actual.1[2].finite_body.unwrap().body_after
    );
    let SharedLubricationSurface3::Prescribed(plane) = prescribed() else {
        unreachable!()
    };
    assert_eq!(prepare(&[prescribed()], &[], &[], &gate).unwrap(), [plane]);
    assert!(prepare(&[], &[], &[], &gate).unwrap().is_empty());
    let (faces, groups) = collect(&[], &[], &[], &gate).unwrap();
    assert!(faces.is_empty() && groups.is_empty());
}

#[test]
fn shared_surface_admission_and_receipt_cancellation_return_no_partial_inventory() {
    let carrier = body([2.2, 3.3, 4.7], 4.0, [-0.001, 0.003, 0.002]);
    let good = face(0, [2.495, 0.0, 0.0], [1.0, 0.0, 0.0]);
    let gate = CancelGate::new();
    assert!(prepare(&[good], &[], &[], &gate).is_err());
    assert!(prepare(&[], &[carrier], &[[0.0; 6]], &gate).is_err());
    assert!(prepare(&[good], &[carrier], &[], &gate).is_err());
    assert!(prepare(&[good], &[carrier], &[[f64::NAN; 6]], &gate).is_err());
    assert!(prepare(&[good; 65], &[carrier], &[[0.0; 6]], &gate).is_err());
    assert!(prepare(&[good], &[carrier; 65], &[[0.0; 6]; 65], &gate).is_err());
    let missing = face(1, [2.495, 0.0, 0.0], [1.0, 0.0, 0.0]);
    assert!(prepare(&[missing], &[carrier], &[[0.0; 6]], &gate).is_err());
    let original = prepare(&[good], &[carrier], &[[0.0; 6]], &gate).unwrap();
    gate.request();
    assert!(prepare(&[good], &[carrier], &[[0.0; 6]], &gate).is_err());
    assert!(collect(&[], &[], &[], &gate).is_err());
    assert_eq!(
        prepare(&[good], &[carrier], &[[0.0; 6]], &CancelGate::new()).unwrap(),
        original
    );
}
