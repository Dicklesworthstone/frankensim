//! Original resistance solve with an explicitly rotating prescribed surface.
use super::*;

const RADIUS: f64 = 0.5;
const AXIS: [f64; 3] = [0.0, 0.0, 2.0];

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}
fn vector_close(a: [f64; 3], b: [f64; 3], tolerance: f64) {
    for i in 0..3 {
        close(a[i], b[i], tolerance);
    }
}
fn rotate(v: [f64; 3]) -> [f64; 3] {
    [0.8 * v[0] + 0.6 * v[2], v[1], -0.6 * v[0] + 0.8 * v[2]]
}
fn body(center: [f64; 3]) -> PartialSaturationBody3 {
    PartialSaturationBody3 {
        center,
        mass: 2.0,
        inertia_world: [[0.2, 0.0, 0.0], [0.0, 0.2, 0.0], [0.0, 0.0, 0.2]],
        velocity: [0.0; 3],
        angular_velocity: [0.0; 3],
    }
}
fn plane(velocity: [f64; 3], spin: f64) -> LubricationPlane3 {
    LubricationPlane3::new(AXIS, [0.0, 0.0, 1.0], velocity, 0.04, 0.0001)
        .unwrap()
        .with_tangential(0.01, 0.0001)
        .unwrap()
        .with_axial_spin(spin)
        .unwrap()
}
fn run<const N: usize>(
    bodies: [PartialSaturationBody3; N],
    responses: [Response; N],
    plane: LubricationPlane3,
    origin: [f64; 3],
) -> (
    [PartialSaturationBody3; N],
    SphereLubricationExchange3<N>,
    PlaneLubricationExchange3,
) {
    let model = SphereLubrication3::new(0.11, 0.04, 0.0001, 64).unwrap();
    model.admit([RADIUS; N]).unwrap();
    planes::admit_planes(&[plane]).unwrap();
    let mut proposed = bodies;
    let (report, planes) = solve(
        model,
        [RADIUS; N],
        &[plane],
        &bodies,
        &responses,
        &mut proposed,
        origin,
        &CancelGate::new(),
    )
    .unwrap();
    (proposed, report, planes[0])
}

#[test]
fn axial_surface_velocity_follows_its_translating_axis_and_rotates_covariantly() {
    let translation = [0.003, -0.002, 0.001];
    let original = plane(translation, 0.006);
    let point = [2.0, 0.5, 2.0];
    let expected = [0.0, 0.010, 0.001];
    vector_close(original.velocity_at(point).unwrap(), expected, 3e-18);
    assert_eq!(original.axial_spin(), 0.006);
    let shifted_point = core::array::from_fn(|i| point[i] + 0.5 * translation[i]);
    vector_close(
        original
            .translated(0.5)
            .unwrap()
            .velocity_at(shifted_point)
            .unwrap(),
        expected,
        4e-18,
    );
    let rotated = LubricationPlane3::new(
        rotate(AXIS),
        rotate([0.0, 0.0, 1.0]),
        rotate(translation),
        0.04,
        0.0001,
    )
    .unwrap()
    .with_axial_spin(0.006)
    .unwrap();
    vector_close(
        rotated.velocity_at(rotate(point)).unwrap(),
        rotate(expected),
        8e-18,
    );
    assert!(original.velocity_at([f64::NAN, 0.0, 0.0]).is_err());
    assert!(original.with_axial_spin(f64::INFINITY).is_err());
    assert!(original.with_axial_spin(f64::NAN).is_err());
}

#[test]
fn local_rotating_wall_matches_translation_without_erasing_rolling_resistance() {
    for translation in [
        [0.002, 0.0, 0.0],
        [0.0, -0.012, 0.0],
        [0.0, -0.012 + 1e-14, 0.0],
    ] {
        let spinning = plane(translation, 0.006);
        let foot = [2.0, 0.0, 2.0];
        let local = spinning.velocity_at(foot).unwrap();
        let mut initial = body([2.0, 0.0, 2.505]);
        initial.velocity = [local[0] + 0.007, local[1], local[2]];
        initial.angular_velocity = [0.0, 0.014, 0.005];
        // Relative rolling has zero point slip, but the native pressure-driven
        // rolling mode still dissipates. Axial torsion is not fabricated.
        close(
            initial.velocity[0] - local[0] - RADIUS * initial.angular_velocity[1],
            0.0,
            2e-18,
        );
        let translating = plane(local, 0.0);
        let (a, film, support) = run([initial], [Response::default()], spinning, [-2.0, 0.7, 1.2]);
        let (b, reference, equivalent) = run(
            [initial],
            [Response::default()],
            translating,
            [-2.0, 0.7, 1.2],
        );
        assert_eq!(a, b);
        assert_eq!(film.body_impulses, reference.body_impulses);
        assert_eq!(film.body_spin_impulses, reference.body_spin_impulses);
        assert_eq!(film.dissipated_energy, reference.dissipated_energy);
        assert!(film.dissipated_energy > 0.0);
        assert_eq!(a[0].angular_velocity[2], initial.angular_velocity[2]);
        close(support.work, equivalent.work, 2e-19);
        close(support.work, dot(support.impulse, local), 2e-19);
        close(film.work_residual, 0.0, 2e-18);
        vector_close(
            core::array::from_fn(|i| film.body_angular_impulses[0][i] + support.angular_impulse[i]),
            [0.0; 3],
            2e-16,
        );
    }
}

#[test]
fn opposite_footpoints_transfer_spin_work_with_equal_opposite_world_impulse() {
    let spinning = plane([0.0; 3], 0.008);
    let initial = [body([-1.2, 0.0, 2.505]), body([1.2, 0.0, 2.505])];
    let (after, film, support) = run(initial, [Response::default(); 2], spinning, [0.0; 3]);
    assert_eq!(film.active_pairs, 0);
    assert_eq!(film.active_plane_pairs, 2);
    assert_eq!(film.active_resistance_modes, 10);
    assert!(after[0].velocity[1] < 0.0 && after[1].velocity[1] > 0.0);
    assert!(support.work < 0.0 && support.angular_impulse[2] < 0.0);
    assert!(film.dissipated_energy > 0.0);
    close(
        support.work,
        spinning.axial_spin() * support.angular_impulse[2],
        2e-19,
    );
    for i in 0..3 {
        close(
            film.body_impulses[0][i] + film.body_impulses[1][i] + support.impulse[i],
            0.0,
            2e-17,
        );
        close(
            film.body_angular_impulses[0][i]
                + film.body_angular_impulses[1][i]
                + support.angular_impulse[i],
            0.0,
            2e-17,
        );
    }
    close(film.work_residual, 0.0, 2e-18);
    // Reporting far away must not subtract large torques to recover work.
    let origin = [1024.0, -2048.0, 4096.0];
    let (shifted_after, shifted, shifted_support) =
        run(initial, [Response::default(); 2], spinning, origin);
    assert_eq!(shifted_after, after);
    assert_eq!(shifted.plane_work, film.plane_work);
    assert_eq!(shifted_support.work, support.work);
    let origin_moment = cross(origin, support.impulse);
    vector_close(
        shifted_support.angular_impulse,
        core::array::from_fn(|i| support.angular_impulse[i] - origin_moment[i]),
        3e-14,
    );
}

#[test]
fn spinning_plane_uses_the_original_full_fluid_body_response_under_rotation() {
    let spinning = plane([0.002, -0.012, 0.0], 0.006);
    let mut initial = body([2.0, 0.2, 2.505]);
    initial.velocity = [0.020, -0.010, 0.002];
    initial.angular_velocity = [0.010, 0.020, 0.005];
    let external = [0.001, -0.002, 0.0003, 0.0001, -0.0002, 0.0003];
    let samples = [
        (0.3, [0.2, -0.1, 0.1], [-0.01, 0.003, -0.004]),
        (0.15, [-0.1, 0.3, 0.2], [0.002, -0.004, 0.001]),
    ];
    let mut response = Response::new(external);
    for (weight, arm, slip) in samples {
        response.add(weight, arm, slip);
    }
    let origin = [-2.0, 0.7, 1.2];
    let (after, film, support) = run([initial], [response], spinning, origin);
    let rotated_body = PartialSaturationBody3 {
        center: rotate(initial.center),
        velocity: rotate(initial.velocity),
        angular_velocity: rotate(initial.angular_velocity),
        ..initial
    };
    let rotated_plane = LubricationPlane3::new(
        rotate(AXIS),
        rotate([0.0, 0.0, 1.0]),
        rotate(spinning.velocity()),
        0.04,
        0.0001,
    )
    .unwrap()
    .with_tangential(0.01, 0.0001)
    .unwrap()
    .with_axial_spin(spinning.axial_spin())
    .unwrap();
    let force = rotate(external[..3].try_into().unwrap());
    let torque = rotate(external[3..].try_into().unwrap());
    let mut rotated_response = Response::new([
        force[0], force[1], force[2], torque[0], torque[1], torque[2],
    ]);
    for (weight, arm, slip) in samples {
        rotated_response.add(weight, rotate(arm), rotate(slip));
    }
    let (rotated_after, rotated_film, rotated_support) = run(
        [rotated_body],
        [rotated_response],
        rotated_plane,
        rotate(origin),
    );
    vector_close(rotated_after[0].velocity, rotate(after[0].velocity), 3e-15);
    vector_close(
        rotated_after[0].angular_velocity,
        rotate(after[0].angular_velocity),
        3e-15,
    );
    vector_close(
        rotated_film.body_impulses[0],
        rotate(film.body_impulses[0]),
        3e-15,
    );
    vector_close(
        rotated_film.body_spin_impulses[0],
        rotate(film.body_spin_impulses[0]),
        3e-15,
    );
    vector_close(
        rotated_support.angular_impulse,
        rotate(support.angular_impulse),
        2e-14,
    );
    close(rotated_support.work, support.work, 2e-17);
    close(
        rotated_film.dissipated_energy,
        film.dissipated_energy,
        2e-17,
    );
    close(rotated_film.work_residual, 0.0, 3e-18);
    let arm = core::array::from_fn(|i| AXIS[i] - origin[i]);
    let shift = cross(arm, support.impulse);
    let work_from_wrench = dot(support.impulse, spinning.velocity())
        + spinning.axial_spin() * (support.angular_impulse[2] - shift[2]);
    close(support.work, work_from_wrench, 3e-19);
}

#[test]
fn common_axial_motion_is_null_and_no_finite_torsional_drag_is_invented() {
    let spinning = plane([0.0; 3], 0.01);
    for center in [[1.2, 0.0, 2.505], [0.0, 0.0, 2.505]] {
        let mut initial = body(center);
        initial.velocity = spinning.velocity_at([center[0], center[1], 2.0]).unwrap();
        initial.angular_velocity = [0.0, 0.0, 0.01];
        let (after, film, support) = run([initial], [Response::default()], spinning, [0.0; 3]);
        assert_eq!(after, [initial]);
        assert_eq!(film.body_impulses, [[0.0; 3]]);
        assert_eq!(film.body_spin_impulses, [[0.0; 3]]);
        assert_eq!(film.dissipated_energy, 0.0);
        assert_eq!(support.work, 0.0);
    }
    let initial = body([0.0, 0.0, 2.505]);
    let (after, film, _) = run([initial], [Response::default()], spinning, [0.0; 3]);
    assert_eq!(after, [initial]);
    assert_eq!(film.body_spin_impulses, [[0.0; 3]]);
}

#[test]
fn zero_spin_retains_native_bits_and_local_speed_refuses_before_proposal() {
    let original = LubricationPlane3::new(AXIS, [0.0, 0.0, 1.0], [0.003, 0.0, 0.0], 0.04, 0.0001)
        .unwrap()
        .with_tangential(0.01, 0.0001)
        .unwrap();
    assert_eq!(original, original.with_axial_spin(0.0).unwrap());
    let mut finite = original;
    finite.recoil = Some(finite_planes::FinitePlaneInput {
        body: body(AXIS),
        external: [0.0; 6],
    });
    assert!(finite.with_axial_spin(0.01).is_err());
    assert_eq!(finite.with_axial_spin(0.0).unwrap(), finite);
    let initial = [body([1.2, 0.0, 2.505])];
    let response = [Response::new([0.002, 0.0, 0.0, 0.0, 0.0001, 0.0])];
    assert_eq!(
        run(initial, response, original, [0.0; 3]),
        run(
            initial,
            response,
            original.with_axial_spin(0.0).unwrap(),
            [0.0; 3]
        )
    );
    let fast = original.with_axial_spin(1.0).unwrap();
    assert!(fast.velocity_at([1.2, 0.0, 2.0]).is_err());
    let mut proposed = initial;
    let model = SphereLubrication3::new(0.11, 0.04, 0.0001, 1).unwrap();
    assert!(
        solve(
            model,
            [RADIUS],
            &[fast],
            &initial,
            &response,
            &mut proposed,
            [0.0; 3],
            &CancelGate::new()
        )
        .is_err()
    );
    assert_eq!(proposed, initial);
    let gate = CancelGate::new();
    gate.request();
    assert!(
        solve(
            model,
            [RADIUS],
            &[original],
            &initial,
            &response,
            &mut proposed,
            [0.0; 3],
            &gate
        )
        .is_err()
    );
    assert_eq!(proposed, initial);
}
