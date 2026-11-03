use super::*;

fn near(a: [f64; 3], b: [f64; 3]) {
    for (a, b) in a.into_iter().zip(b) {
        assert!((a - b).abs() < 2e-11, "{a} != {b}");
    }
}
fn elbow() -> PipeLoadGeometry {
    PipeLoadGeometry::try_new(
        vec![[0.0; 3], [2.0, 0.0, 0.0], [2.0, 1.0, 0.0]],
        Area::new(0.5),
        2,
        || false,
    )
    .unwrap()
}
fn flow(q: [f64; 2]) -> FlowInterval {
    FlowInterval {
        duration_s: 0.25,
        volume_flows_m3_s: q,
        density_kg_m3: 3.0,
        pressures_pa: [10.0, 4.0],
        exterior_pressure_pa: 2.0,
        gravity_m_s2: [0.0; 3],
    }
}

#[test]
fn g1_bend_pressure_traction_and_both_directions_of_reversal() {
    let pipe = elbow();
    for q in [2.0, -2.0] {
        let load = pipe
            .evaluate(flow([q, q]), [0.0; 3], |_| 3.0, || false)
            .unwrap();
        // Pressure contributes (4,-1) N; momentum turning contributes (24,-24) N.
        near(load.force_n, [28.0, -25.0, 0.0]);
        near(load.torque_nm, [0.0, 0.0, -50.0]);
    }
    let reversing = pipe
        .evaluate(flow([2.0, -2.0]), [0.0; 3], |_| 3.0, || false)
        .unwrap();
    // Integral Q^2 dt = 1/3 m^6/s, despite zero net transported volume.
    // Stored fluid momentum change is (-24,-12) N s and spin change -24 N m s.
    near(reversing.force_n, [108.0, 39.0, 0.0]);
    near(reversing.torque_nm, [0.0, 0.0, 78.0]);
    // Splitting the SAME affine trajectory preserves total delivered impulse.
    let mut a = flow([2.0, 0.0]);
    a.duration_s = 0.125;
    let mut b = flow([0.0, -2.0]);
    b.duration_s = 0.125;
    let a = pipe.evaluate(a, [0.0; 3], |_| 3.0, || false).unwrap();
    let b = pipe.evaluate(b, [0.0; 3], |_| 3.0, || false).unwrap();
    near(add(a.impulse_ns, b.impulse_ns), reversing.impulse_ns);
    near(
        add(a.angular_impulse_nms, b.angular_impulse_nms),
        reversing.angular_impulse_nms,
    );
}

#[test]
fn g1_straight_acceleration_retains_fluid_impulse_before_loading_the_pipe() {
    let pipe =
        PipeLoadGeometry::try_new(vec![[0.0; 3], [4.0, 0.0, 0.0]], Area::new(0.25), 1, || {
            false
        })
        .unwrap();
    let mut sample = flow([0.0, 3.0]);
    sample.density_kg_m3 = 2.0;
    sample.duration_s = 2.0;
    sample.pressures_pa = [61.0, 7.0];
    let load = pipe.evaluate(sample, [0.0; 3], |_| 2.0, || false).unwrap();
    // 13.5 N pressure force - 12 N accelerating retained fluid = 1.5 N wall load.
    near(load.force_n, [1.5, 0.0, 0.0]);
    near(load.impulse_ns, [3.0, 0.0, 0.0]);
    near(load.torque_nm, [0.0; 3]);
}

#[test]
fn g1_original_local_weight_and_arbitrary_moment_origin() {
    let pipe = elbow();
    let mut sample = flow([0.0; 2]);
    sample.pressures_pa = [2.0; 2];
    sample.gravity_m_s2 = [0.0, -10.0, 0.0];
    let load = pipe
        .evaluate(
            sample,
            [5.0, -2.0, 1.0],
            |section| [4.0, 8.0][section],
            || false,
        )
        .unwrap();
    // Each original section weighs 40 N; their centroids are (1,0,0),(2,0.5,0).
    near(load.force_n, [0.0, -80.0, 0.0]);
    near(load.torque_nm, [-80.0, 0.0, 280.0]);
}

#[test]
fn g3_pressure_datum_translation_and_rotation_preserve_the_physical_wrench() {
    let pipe = elbow();
    let original = pipe
        .evaluate(flow([2.0, -2.0]), [1.0, 2.0, 3.0], |_| 3.0, || false)
        .unwrap();
    let rotate = |x: [f64; 3]| [-x[1], x[0], x[2]];
    let shift = [10.0, -20.0, 30.0];
    let moved = PipeLoadGeometry::try_new(
        pipe.vertices_m
            .iter()
            .map(|&x| add(rotate(x), shift))
            .collect(),
        Area::new(0.5),
        2,
        || false,
    )
    .unwrap();
    let mut sample = flow([2.0, -2.0]);
    sample.pressures_pa = sample.pressures_pa.map(|x| x + 1e5);
    sample.exterior_pressure_pa += 1e5;
    let result = moved
        .evaluate(
            sample,
            add(rotate(original.about_m), shift),
            |_| 3.0,
            || false,
        )
        .unwrap();
    near(result.force_n, rotate(original.force_n));
    near(result.torque_nm, rotate(original.torque_nm));
}

#[test]
fn g0_invalid_geometry_interval_and_cancellation_refuse() {
    for vertices in [vec![], vec![[0.0; 3]], vec![[0.0; 3]; 2]] {
        assert!(PipeLoadGeometry::try_new(vertices, Area::new(1.0), 2, || false).is_err());
    }
    assert!(
        PipeLoadGeometry::try_new(vec![[0.0; 3], [1.0, 0.0, 0.0]], Area::new(1.0), 0, || false)
            .is_err()
    );
    let pipe = elbow();
    let mut invalid = flow([0.0; 2]);
    invalid.duration_s = 0.0;
    assert!(pipe.evaluate(invalid, [0.0; 3], |_| 3.0, || false).is_err());
    assert!(
        pipe.evaluate(flow([0.0; 2]), [0.0; 3], |_| 0.0, || false)
            .is_err()
    );
    assert_eq!(
        pipe.evaluate(flow([0.0; 2]), [0.0; 3], |_| 3.0, || true),
        Err(PipeLoadError::Cancelled)
    );
}

#[test]
fn g1_pipe_wrench_drives_the_native_rigid_body_impulse_and_work() {
    use fs_ga::{Mat3, Se3, Vec3, Wrench};
    use fs_mbd::articulated::{
        ArticulatedModel, FreeFloatingBaseState, JointModel, Link, SpatialInertia,
        apply_free_floating_impulses,
    };
    let pipe = elbow();
    let load = pipe
        .evaluate(flow([2.0; 2]), [0.0; 3], |_| 3.0, || false)
        .unwrap();
    // Pipe coordinates already coincide with the body's axes; the chosen load
    // reference is its COM. This is a one-way mechanical kick, not moving-pipe FSI.
    let inertia = SpatialInertia::new(
        2.0,
        Vec3::new(0.0, 0.0, 0.0),
        Mat3 {
            m: [3.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 5.0],
        },
    )
    .unwrap();
    let model = ArticulatedModel::new(vec![Link::new(
        "pipe support",
        None,
        Se3::identity(),
        JointModel::FIXED,
        inertia,
    )])
    .unwrap();
    let v = |a: [f64; 3]| Vec3::new(a[0], a[1], a[2]);
    let response = apply_free_floating_impulses(
        &model,
        FreeFloatingBaseState::stationary(Se3::identity()),
        &[],
        &[],
        &[],
        &[Wrench::new(v(load.angular_impulse_nms), v(load.impulse_ns))],
    )
    .unwrap();
    let array = |v: Vec3| [v.x, v.y, v.z];
    near(
        array(response.linear_momentum_change_world),
        [7.0, -6.25, 0.0],
    );
    near(
        array(response.angular_momentum_change_world),
        [0.0, 0.0, -12.5],
    );
    // J^2/(2 M) + K_z^2/(2 I_z), independently from the original native tensor.
    assert!((response.impulse_work - 37.640625).abs() < 2e-11);
    assert!(response.kinetic_work_residual.abs() < 2e-11);
}
