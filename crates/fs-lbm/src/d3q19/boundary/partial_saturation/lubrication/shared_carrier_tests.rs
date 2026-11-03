//! G1/G3 common spatial inertia across distinct film faces; G4 refusal.
use super::*;

const RADIUS: f64 = 0.5;
const CUTOFF: f64 = 0.04;
const ROUGHNESS: f64 = 0.0001;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}
fn vector_close(actual: [f64; 3], expected: [f64; 3], tolerance: f64) {
    for i in 0..3 {
        close(actual[i], expected[i], tolerance);
    }
}
fn body(center: [f64; 3]) -> PartialSaturationBody3 {
    PartialSaturationBody3 {
        center,
        mass: 2.0,
        inertia_world: [[0.2, 0.01, 0.02], [0.01, 0.25, -0.01], [0.02, -0.01, 0.3]],
        velocity: [0.0; 3],
        angular_velocity: [0.0; 3],
    }
}
fn face(
    point: [f64; 3],
    normal: [f64; 3],
    body: PartialSaturationBody3,
    external: [f64; 6],
    anchor: usize,
    tangent: bool,
) -> LubricationPlane3 {
    let plane = FiniteLubricationPlane3::new(point, normal, body, CUTOFF, ROUGHNESS).unwrap();
    let plane = if tangent {
        plane.with_tangential(0.01, ROUGHNESS).unwrap()
    } else {
        plane
    };
    let mut native = plane.driven(external).unwrap();
    native.recoil.as_mut().unwrap().carrier_anchor = Some(anchor);
    native
}
fn model() -> SphereLubrication3 {
    SphereLubrication3::new(0.08, CUTOFF, ROUGHNESS, 64)
        .unwrap()
        .with_resistance_mode_limit(64)
        .unwrap()
}
fn run<const N: usize>(
    bodies: [PartialSaturationBody3; N],
    external: [[f64; 6]; N],
    planes: &[LubricationPlane3],
    origin: [f64; 3],
) -> (
    [PartialSaturationBody3; N],
    SphereLubricationExchange3<N>,
    Vec<PlaneLubricationExchange3>,
) {
    run_response(bodies, external.map(Response::new), planes, origin)
}
fn run_response<const N: usize>(
    bodies: [PartialSaturationBody3; N],
    response: [Response; N],
    planes: &[LubricationPlane3],
    origin: [f64; 3],
) -> (
    [PartialSaturationBody3; N],
    SphereLubricationExchange3<N>,
    Vec<PlaneLubricationExchange3>,
) {
    let model = model();
    model.admit([RADIUS; N]).unwrap();
    planes::admit_planes(planes).unwrap();
    let mut after = bodies;
    let (report, planes) = solve(
        model,
        [RADIUS; N],
        planes,
        &bodies,
        &response,
        &mut after,
        origin,
        &CancelGate::new(),
    )
    .unwrap();
    (after, report, planes)
}
fn total<const P: usize>(items: [[f64; 3]; P]) -> [f64; 3] {
    core::array::from_fn(|i| items.into_iter().map(|v| v[i]).sum())
}
fn rotate(v: [f64; 3]) -> [f64; 3] {
    [0.8 * v[0] + 0.6 * v[2], v[1], -0.6 * v[0] + 0.8 * v[2]]
}
fn rotated_body(b: PartialSaturationBody3) -> PartialSaturationBody3 {
    let r = [[0.8, 0.0, 0.6], [0.0, 1.0, 0.0], [-0.6, 0.0, 0.8]];
    let mut inertia_world = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..=i {
            let entry = (0..3)
                .flat_map(|a| {
                    (0..3).map(move |b_axis| r[i][a] * b.inertia_world[a][b_axis] * r[j][b_axis])
                })
                .sum();
            // The native contract admits exactly symmetric world tensors.
            inertia_world[i][j] = entry;
            inertia_world[j][i] = entry;
        }
    }
    PartialSaturationBody3 {
        center: rotate(b.center),
        velocity: rotate(b.velocity),
        angular_velocity: rotate(b.angular_velocity),
        inertia_world,
        ..b
    }
}
fn rotated_load(load: [f64; 6]) -> [f64; 6] {
    let force = rotate(core::array::from_fn(|i| load[i]));
    let torque = rotate(core::array::from_fn(|i| load[i + 3]));
    core::array::from_fn(|i| if i < 3 { force[i] } else { torque[i - 3] })
}

#[test]
fn two_nonparallel_faces_match_one_mass_and_fluid_independent_midpoint_oracle() {
    let mut sphere = body([0.8, 0.7, 0.3]);
    sphere.velocity = [0.01, -0.004, 0.0];
    let mut carrier = body(sphere.center);
    carrier.mass = 5.0;
    carrier.velocity = [-0.003, 0.002, 0.0];
    let external = [0.001, -0.002, 0.0, 0.0, 0.0, 0.0];
    let load = [-0.002, 0.001, 0.0, 0.0, 0.0, 0.0];
    let fluid_weight = 0.8;
    let fluid_slip = [0.004, -0.002, 0.0];
    let mut response = Response::new(external);
    response.add(fluid_weight, [0.0; 3], fluid_slip);
    let sphere_response = sphere.mass + 0.5 * fluid_weight;
    let drive: [f64; 3] = core::array::from_fn(|i| external[i] + fluid_weight * fluid_slip[i]);
    let normals = [
        [1.0, 0.0, 0.0],
        [
            core::f64::consts::FRAC_1_SQRT_2,
            core::f64::consts::FRAC_1_SQRT_2,
            0.0,
        ],
    ];
    let planes = core::array::from_fn::<_, 2, _>(|i| {
        face(
            core::array::from_fn(|j| {
                sphere.center[j] - (RADIUS + [0.01, 0.015][i]) * normals[i][j]
            }),
            normals[i],
            carrier,
            load,
            0,
            false,
        )
    });
    let directions = planes.map(|p| p.normal);
    let zeta = planes.map(|p| {
        6.0 * core::f64::consts::PI
            * 0.08
            * RADIUS
            * RADIUS
            * (1.0 / (p.distance(sphere.center) - RADIUS) - 1.0 / CUTOFF)
    });
    // Independent two-mode equation including native midpoint fluid mobility.
    // Both diagonal AND cross-face mobility contain the same one carrier.
    let c = 0.5 * (1.0 / sphere_response + 1.0 / carrier.mass);
    let diagonal = [1.0 / zeta[0] + c, 1.0 / zeta[1] + c];
    let off = c * dot(directions[0], directions[1]);
    let free_slip = core::array::from_fn(|i| {
        sphere.velocity[i] - carrier.velocity[i]
            + 0.5 * (drive[i] / sphere_response - load[i] / carrier.mass)
    });
    let rhs = directions.map(|n| -dot(n, free_slip));
    let det = diagonal[0] * diagonal[1] - off * off;
    let lambda = [
        (diagonal[1] * rhs[0] - off * rhs[1]) / det,
        (diagonal[0] * rhs[1] - off * rhs[0]) / det,
    ];
    let film = total(core::array::from_fn::<_, 2, _>(|i| {
        directions[i].map(|n| lambda[i] * n)
    }));
    let expected_sphere =
        core::array::from_fn(|i| sphere.velocity[i] + (drive[i] + film[i]) / sphere_response);
    let expected_carrier =
        core::array::from_fn(|i| carrier.velocity[i] + (load[i] - film[i]) / carrier.mass);
    let (after, report, faces) = run_response([sphere], [response], &planes, [1.0, -2.0, 0.7]);
    vector_close(after[0].velocity, expected_sphere, 2e-17);
    for i in 0..3 {
        let delta = after[0].velocity[i] - sphere.velocity[i];
        let fluid_impulse = fluid_weight * (fluid_slip[i] - 0.5 * delta);
        close(
            sphere.mass * delta,
            external[i] + fluid_impulse + film[i],
            4e-17,
        );
    }
    let shared = faces[0].finite_body.unwrap();
    vector_close(shared.body_after.velocity, expected_carrier, 2e-17);
    assert_eq!(shared.body_after, faces[1].finite_body.unwrap().body_after);
    for i in 0..2 {
        vector_close(
            faces[i].impulse,
            directions[i].map(|n| -lambda[i] * n),
            5e-17,
        );
        vector_close(faces[i].finite_body.unwrap().spin_impulse, [0.0; 3], 1e-18);
        close(
            faces[i].work,
            dot(faces[i].impulse, shared.velocity_midpoint),
            2e-19,
        );
    }
    vector_close(report.body_impulses[0], film, 4e-17);
    close(shared.work_residual, 0.0, 3e-19);
    close(report.work_residual, 0.0, 4e-19);
    close(
        shared.kinetic_energy_change,
        shared.external_work + faces.iter().map(|f| f.work).sum::<f64>(),
        3e-19,
    );
    assert_eq!(report.active_resistance_modes, 2);
    // Independently-owned copies give a different cross-face response.
    let independent = planes.map(|mut p| {
        p.recoil.as_mut().unwrap().carrier_anchor = None;
        p
    });
    let (different, _, _) = run_response([sphere], [response], &independent, [1.0, -2.0, 0.7]);
    assert!((different[0].velocity[0] - after[0].velocity[0]).abs() > 1e-5);
}

fn spatial_case() -> (
    PartialSaturationBody3,
    PartialSaturationBody3,
    [[f64; 3]; 2],
    [[f64; 3]; 2],
    [f64; 6],
    [f64; 6],
) {
    let mut sphere = body([0.505, 0.35, 0.505]);
    sphere.velocity = [-0.005, 0.006, -0.004];
    sphere.angular_velocity = [0.005, -0.004, 0.002];
    let mut carrier = body([-0.2, 0.1, -0.3]);
    carrier.mass = 3.5;
    carrier.inertia_world = [[0.9, 0.11, -0.07], [0.11, 1.2, 0.05], [-0.07, 0.05, 0.7]];
    carrier.velocity = [0.003, -0.002, 0.004];
    carrier.angular_velocity = [-0.004, 0.002, 0.003];
    (
        sphere,
        carrier,
        [[0.0, 0.0, 0.0]; 2],
        [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0]],
        [0.0004, -0.0003, 0.0002, 0.0001, 0.0002, -0.0001],
        [-0.0003, 0.0002, 0.0005, -0.0001, 0.0001, 0.0002],
    )
}

#[test]
fn shared_full_tensor_tangent_faces_close_unique_wrench_work_and_transform() {
    let (sphere, carrier, points, normals, external, load) = spatial_case();
    let planes =
        core::array::from_fn::<_, 2, _>(|i| face(points[i], normals[i], carrier, load, 0, true));
    let origin = [-1.0, 2.0, 0.7];
    let (after, report, faces) = run([sphere], [external], &planes, origin);
    let full = faces[0].finite_body.unwrap();
    let force = total([faces[0].impulse, faces[1].impulse]);
    let torque = total([
        faces[0].finite_body.unwrap().spin_impulse,
        faces[1].finite_body.unwrap().spin_impulse,
    ]);
    let dp = core::array::from_fn(|i| {
        carrier.mass * (full.body_after.velocity[i] - carrier.velocity[i])
    });
    let dw =
        core::array::from_fn(|i| full.body_after.angular_velocity[i] - carrier.angular_velocity[i]);
    let dl = carrier.inertia_world.map(|row| dot(row, dw));
    vector_close(dp, core::array::from_fn(|i| load[i] + force[i]), 5e-17);
    vector_close(dl, core::array::from_fn(|i| load[i + 3] + torque[i]), 5e-17);
    vector_close(
        total([
            report.body_angular_impulses[0],
            faces[0].angular_impulse,
            faces[1].angular_impulse,
        ]),
        [0.0; 3],
        5e-17,
    );
    let kinetic = dot(dp, full.velocity_midpoint) + dot(dl, full.angular_velocity_midpoint);
    close(full.kinetic_energy_change, kinetic, 1e-18);
    close(
        kinetic - full.external_work - faces.iter().map(|f| f.work).sum::<f64>(),
        0.0,
        2e-18,
    );
    close(report.work_residual, 0.0, 3e-18);
    assert!(report.dissipated_energy > 0.0);
    assert!(
        faces
            .iter()
            .all(|f| f.finite_body.unwrap().body_after == full.body_after)
    );
    let rotated = core::array::from_fn::<_, 2, _>(|i| {
        face(
            rotate(points[i]),
            rotate(normals[i]),
            rotated_body(carrier),
            rotated_load(load),
            0,
            true,
        )
    });
    let (rotated_after, rotated_report, rotated_faces) = run(
        [rotated_body(sphere)],
        [rotated_load(external)],
        &rotated,
        rotate(origin),
    );
    vector_close(rotated_after[0].velocity, rotate(after[0].velocity), 3e-16);
    vector_close(
        rotated_after[0].angular_velocity,
        rotate(after[0].angular_velocity),
        3e-16,
    );
    vector_close(
        rotated_faces[0].finite_body.unwrap().body_after.velocity,
        rotate(full.body_after.velocity),
        3e-16,
    );
    vector_close(
        rotated_faces[0]
            .finite_body
            .unwrap()
            .body_after
            .angular_velocity,
        rotate(full.body_after.angular_velocity),
        3e-16,
    );
    close(
        rotated_report.dissipated_energy,
        report.dissipated_energy,
        2e-18,
    );
    close(rotated_report.plane_work, report.plane_work, 2e-18);
    let shifted_origin = [1e6, -2e6, 3e6];
    let (shifted_after, shifted_report, shifted_faces) =
        run([sphere], [external], &planes, shifted_origin);
    assert_eq!(shifted_after, after);
    assert_eq!(shifted_report.body_spin_impulses, report.body_spin_impulses);
    assert_eq!(shifted_report.plane_work, report.plane_work);
    assert_eq!(shifted_faces[0].finite_body, faces[0].finite_body);
}

#[test]
fn common_rigid_motion_is_null_and_face_permutation_keeps_one_carrier() {
    let (mut sphere, mut carrier, points, normals, _, _) = spatial_case();
    let translation = [0.003, -0.002, 0.001];
    let omega = [0.002, -0.001, 0.003];
    sphere.velocity = core::array::from_fn(|i| translation[i] + cross(omega, sphere.center)[i]);
    carrier.velocity = core::array::from_fn(|i| translation[i] + cross(omega, carrier.center)[i]);
    sphere.angular_velocity = omega;
    carrier.angular_velocity = omega;
    let planes = core::array::from_fn::<_, 2, _>(|i| {
        face(points[i], normals[i], carrier, [0.0; 6], 0, true)
    });
    let (after, report, faces) = run([sphere], [[0.0; 6]], &planes, [0.0; 3]);
    vector_close(after[0].velocity, sphere.velocity, 2e-18);
    vector_close(after[0].angular_velocity, sphere.angular_velocity, 2e-18);
    vector_close(
        faces[0].finite_body.unwrap().body_after.velocity,
        carrier.velocity,
        2e-18,
    );
    close(report.dissipated_energy, 0.0, 1e-30);
    vector_close(report.plane_impulse, [0.0; 3], 2e-18);
    let (sphere, carrier, points, normals, external, load) = spatial_case();
    let planes =
        core::array::from_fn::<_, 2, _>(|i| face(points[i], normals[i], carrier, load, 0, true));
    let (a, ra, pa) = run([sphere], [external], &planes, [0.2, -0.3, 0.1]);
    let reversed = [planes[1], planes[0]];
    let (b, rb, pb) = run([sphere], [external], &reversed, [0.2, -0.3, 0.1]);
    vector_close(a[0].velocity, b[0].velocity, 2e-17);
    vector_close(
        pa[0].finite_body.unwrap().body_after.angular_velocity,
        pb[0].finite_body.unwrap().body_after.angular_velocity,
        2e-17,
    );
    close(ra.plane_work, rb.plane_work, 1e-18);
    for i in 0..2 {
        vector_close(pa[i].impulse, pb[1 - i].impulse, 3e-17);
        vector_close(
            pa[i].finite_body.unwrap().spin_impulse,
            pb[1 - i].finite_body.unwrap().spin_impulse,
            3e-17,
        );
    }
}

#[test]
fn inactive_faces_receive_external_impulse_once_and_single_face_keeps_native_bits() {
    let sphere = body([0.8, 0.7, 0.6]);
    let mut carrier = body([0.0; 3]);
    carrier.mass = 3.0;
    carrier.velocity = [0.004, -0.002, 0.001];
    let load = [0.003, -0.002, 0.001, 0.0002, -0.0001, 0.0003];
    let planes = [
        face([0.0, 0.0, -10.0], [0.0, 0.0, 1.0], carrier, load, 0, true),
        face([-10.0, 0.0, 0.0], [1.0, 0.0, 0.0], carrier, load, 0, true),
    ];
    let (after, report, faces) = run([sphere], [[0.0; 6]], &planes, [0.0; 3]);
    assert_eq!(after, [sphere]);
    assert_eq!(report.active_resistance_modes, 0);
    let full = faces[0].finite_body.unwrap();
    let expected = Response::new(load).solve(carrier).unwrap();
    assert_eq!(full.body_after, expected);
    assert_eq!(faces[1].finite_body.unwrap().body_after, expected);
    assert_eq!(report.plane_work, 0.0);
    close(full.kinetic_energy_change - full.external_work, 0.0, 1e-19);
    let independent = planes.map(|mut p| {
        p.recoil.as_mut().unwrap().carrier_anchor = None;
        p
    });
    let single = run([sphere], [[0.0; 6]], &planes[..1], [0.0; 3]);
    let legacy = run([sphere], [[0.0; 6]], &independent[..1], [0.0; 3]);
    assert_eq!(single, legacy);
}

#[test]
fn invalid_carrier_identity_and_cancellation_refuse_before_any_proposal() {
    let sphere = body([0.505, 0.2, 0.505]);
    let carrier = body([0.0; 3]);
    let planes = [
        face([0.0; 3], [0.0, 0.0, 1.0], carrier, [0.0; 6], 0, true),
        face([0.0; 3], [1.0, 0.0, 0.0], carrier, [0.0; 6], 0, true),
    ];
    for kind in 0..4 {
        let mut bad = planes;
        match kind {
            0 => bad[0].recoil.as_mut().unwrap().carrier_anchor = Some(1),
            1 => bad[1].recoil.as_mut().unwrap().body.mass += 1.0,
            2 => bad[1].recoil.as_mut().unwrap().external[0] = 0.001,
            _ => bad[0].recoil = None,
        }
        let mut after = [sphere];
        assert!(
            solve(
                model(),
                [RADIUS],
                &bad,
                &[sphere],
                &[Response::default()],
                &mut after,
                [0.0; 3],
                &CancelGate::new()
            )
            .is_err()
        );
        assert_eq!(after, [sphere]);
    }
    let gate = CancelGate::new();
    gate.request();
    let mut after = [sphere];
    assert!(
        solve(
            model(),
            [RADIUS],
            &planes,
            &[sphere],
            &[Response::default()],
            &mut after,
            [0.0; 3],
            &gate
        )
        .is_err()
    );
    assert_eq!(after, [sphere]);
    let a = run([sphere], [[0.0; 6]], &planes, [0.0; 3]);
    let b = run([sphere], [[0.0; 6]], &planes, [0.0; 3]);
    assert_eq!(a, b);
}
