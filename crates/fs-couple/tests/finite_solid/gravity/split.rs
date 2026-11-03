//! G1/G3/G4: original gravity, consistent mass, contact impulses and split work.
use super::*;

fn near_vector(actual: [f64; 3], expected: [f64; 3], tolerance: f64) {
    for axis in 0..3 {
        near(actual[axis], expected[axis], tolerance);
    }
}

#[test]
fn uniform_gravity_survives_contact_impulses_and_repeated_native_split_motion() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let rho = [3.0, 15.0, 3.0, 3.0];
        let (mass, rows, _) = moments(rho);
        let reference_com = [rows[1] / mass, rows[2] / mass, rows[3] / mass];
        let g = [1.5, -2.0, 0.75];
        let source = uniform(g);
        let origin = [0.1, -0.4, 0.2];
        let mut u = [0.2, -0.1, 0.3];
        let mut v = [-0.3, 0.2, 0.1];
        let mut body = FiniteSolidBody::new_with_gravity_and_sources(
            problem(&material, &[], &[]),
            Some(graded_density(rho)),
            None,
            field(&source),
            &[u; 4],
            &[v; 4],
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        let point = FiniteSolidPoint::new(0, rows.map(|row| row / mass)).unwrap();
        let event_impulse = [0.003, -0.002, 0.001];
        let h = 0.01;
        for step in 0..6 {
            let old = body.native_state().clone();
            let before_com = std::array::from_fn(|axis| reference_com[axis] + u[axis]);
            let kicked_v: [f64; 3] = std::array::from_fn(|axis| v[axis] + h * g[axis]);
            let expected_work = mass
                * h
                * vector(g).dot(vector(std::array::from_fn(|axis| {
                    0.5 * v[axis] + 0.5 * kicked_v[axis]
                })));
            let (kicked, kick) = body
                .force_kick_candidate(h, &[], [0.0; 3], origin, Default::default(), cx)
                .unwrap();
            assert_eq!(body.native_state(), &old);
            assert_eq!(kicked.native_state().q, old.q);
            assert_eq!(kicked.native_state().history, old.history);
            near_vector(
                kick.source_impulse.linear_kg_m_s,
                g.map(|g| mass * h * g),
                2e-11,
            );
            near_vector(
                kick.source_impulse.angular_kg_m2_s,
                components(
                    vector(before_com)
                        .sub(vector(origin))
                        .cross(vector(g).scale(mass * h)),
                ),
                2e-11,
            );
            near(kick.gravity_impulse_work_j, expected_work, 2e-11);
            near(
                kick.gravity_potential_before_j.unwrap(),
                -mass * vector(g).dot(vector(before_com)),
                2e-11,
            );
            let (contacted, _) = kicked
                .point_impulse_candidate(point, event_impulse, origin, Default::default(), cx)
                .unwrap();
            v = std::array::from_fn(|axis| kicked_v[axis] + event_impulse[axis] / mass);
            let dx = v.map(|v| h * v);
            u = std::array::from_fn(|axis| u[axis] + dx[axis]);
            let work = mass * vector(g).dot(vector(dx));
            let defect = expected_work - work;
            let (next, drift) = contacted
                .contact_drift_candidate(&kick, controls(defect.abs() + 1e-9), cx)
                .unwrap();
            near(drift.gravity_work_j, work, 2e-11);
            near(drift.source_displacement_work_j, work, 2e-11);
            near(drift.gravity_splitting_work_j, defect, 2e-11);
            near(drift.integration_defect_j, defect, 2e-11);
            near(drift.absolute_integration_defect_j, defect.abs(), 2e-11);
            for node in 0..4 {
                near_vector(next.displacement_m()[node], u, 2e-11);
                near_vector(next.velocity_m_s()[node], v, 2e-11);
            }
            for (dof, &a) in next.native_state().a.iter().enumerate() {
                near(a, g[dof % 3], 2e-10);
            }
            near(next.native_state().t, (step + 1) as f64 * h, 2e-15);
            assert_eq!(next.native_state().history.len(), step + 1);
            body = next;
        }
        // A split history remains a valid starting point for the native
        // coupled midpoint owner, retaining the same field and mass.
        let (next, _) = body
            .step_candidate(h, &[], [0.0; 3], controls(1e-9), cx)
            .unwrap();
        for node in 0..4 {
            near_vector(
                next.velocity_m_s()[node],
                std::array::from_fn(|axis| v[axis] + h * g[axis]),
                2e-10,
            );
            near_vector(
                next.displacement_m()[node],
                std::array::from_fn(|axis| u[axis] + h * v[axis] + 0.5 * h * h * g[axis]),
                2e-10,
            );
        }
    });
}

#[test]
fn tidal_split_uses_original_p1_moments_pressure_supports_and_endpoint_acceleration() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let fixed: Vec<_> = (0..12)
            .filter(|&dof| dof != 3)
            .map(|dof| (dof, 0.0))
            .collect();
        let source = Tidal {
            g0: 100.0,
            k: 200.0,
        };
        let held = [(3, 0.2)];
        let pressure = 0.7;
        let surface = Arc::new(
            ClosedTetPressureSurface::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx)
                .unwrap(),
        );
        for graded in [false, true] {
            let rho = if graded {
                [3.0, 15.0, 3.0, 3.0]
            } else {
                [6.0; 4]
            };
            let (total_mass, rows, column) = moments(rho);
            let mass = column[1];
            let (stretch, v, h) = (1.2, 0.3, 0.003);
            let mut u0 = [[0.0; 3]; 4];
            let mut v0 = [[0.0; 3]; 4];
            u0[1][0] = stretch - 1.0;
            v0[1][0] = v;
            let body = FiniteSolidBody::new_with_gravity_and_sources(
                problem(&material, &fixed, &held),
                graded.then(|| graded_density(rho)),
                Some((Arc::clone(&surface), TetPressure::Uniform(pressure))),
                field(&source),
                &u0,
                &v0,
                0.0,
                Default::default(),
                cx,
            )
            .unwrap();
            let before = body.native_state().clone();
            let gravity_force = |lambda: f64| source.g0 * rows[1] - source.k * mass * lambda;
            let acceleration = |lambda: f64| {
                (held[0].1 + gravity_force(lambda)
                    - VOLUME * nominal_stress(lambda, 1.0)
                    - pressure * VOLUME)
                    / mass
            };
            let next_v = v + h * acceleration(stretch);
            let next_stretch = stretch + h * next_v;
            let average = 0.5 * (v + next_v);
            let potential = |lambda: f64| {
                0.5 * source.k * mass * lambda * lambda - source.g0 * rows[1] * lambda
            };
            let gravity_work = potential(stretch) - potential(next_stretch);
            let gravity_defect = h * gravity_force(stretch) * average - gravity_work;
            let dx = next_stretch - stretch;
            let dead_defect = held[0].1 * (h * average - dx);
            let pressure_defect = -pressure * VOLUME * (h * average - dx);
            let elastic_defect = strain_energy(next_stretch, 1.0)
                - strain_energy(stretch, 1.0)
                - h * VOLUME * nominal_stress(stretch, 1.0) * average;
            let absolute_defect = gravity_defect.abs()
                + dead_defect.abs()
                + pressure_defect.abs()
                + elastic_defect.abs();
            let (kicked, kick) = body
                .force_kick_candidate_with_pressure(
                    h,
                    &held,
                    [0.0; 3],
                    pressure,
                    [0.0; 3],
                    Default::default(),
                    cx,
                )
                .unwrap();
            near(kicked.native_state().v[0], next_v, 2e-11);
            let total_gravity = source.g0 * total_mass - source.k * rows[1] * stretch;
            near(
                kick.source_impulse.linear_kg_m_s[0],
                h * (held[0].1 + total_gravity),
                2e-11,
            );
            near(
                kick.impulse.support.linear_kg_m_s[0] + kick.source_impulse.linear_kg_m_s[0],
                rows[1] * (next_v - v),
                2e-11,
            );
            // With force along x, its z and y moments depend on the original
            // nodal density second moments, not the cell-mean density.
            near(
                kick.source_impulse.angular_kg_m2_s[1],
                h * (source.g0 * rows[3] - source.k * column[3] * stretch),
                2e-11,
            );
            near(
                kick.source_impulse.angular_kg_m2_s[2],
                -h * (source.g0 * rows[2] - source.k * column[2] * stretch),
                2e-11,
            );
            let (next, drift) = kicked
                .contact_drift_candidate(&kick, controls(absolute_defect + 1e-9), cx)
                .unwrap();
            assert_eq!(body.native_state(), &before);
            near(next.native_state().q[0] + 1.0, next_stretch, 2e-11);
            near(next.native_state().a[0], acceleration(next_stretch), 2e-9);
            near(
                kick.gravity_potential_before_j.unwrap(),
                potential(stretch),
                2e-11,
            );
            near(
                drift.gravity_potential_after_j.unwrap(),
                potential(next_stretch),
                2e-11,
            );
            near(drift.gravity_work_j, gravity_work, 2e-11);
            near(drift.pressure_work_j, -pressure * VOLUME * dx, 2e-11);
            near(
                drift.source_displacement_work_j,
                held[0].1 * dx - pressure * VOLUME * dx + gravity_work,
                2e-11,
            );
            near(drift.gravity_splitting_work_j, gravity_defect, 2e-11);
            near(
                drift.integration_defect_j,
                gravity_defect + dead_defect + pressure_defect + elastic_defect,
                2e-11,
            );
            near(drift.absolute_integration_defect_j, absolute_defect, 2e-11);
        }
    });
}

#[test]
fn opposing_source_defects_cannot_cancel_or_be_changed_through_public_kick_receipts() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let source = uniform([1.0, 0.0, 0.0]);
        let body = FiniteSolidBody::new_with_gravity(
            problem(&material, &[], &[]),
            field(&source),
            &[[0.0; 3]; 4],
            &[[0.2, 0.0, 0.0]; 4],
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        let (kicked, kick) = body
            .force_kick_candidate(0.1, &[], [-3.0, 0.0, 0.0], [0.0; 3], Default::default(), cx)
            .unwrap();
        let before = kicked.native_state().clone();
        // M=1kg, gravity=1m/s², dead acceleration=−0.5m/s², delta-v=.05m/s.
        // The gravity defect is −.0025J and the dead-load defect +.00125J.
        // A .002J gate admits the signed sum but must reject their magnitudes.
        assert!(matches!(
            kicked.contact_drift_candidate(&kick, controls(0.002), cx),
            Err(FiniteSolidError::Energy { .. })
        ));
        assert_eq!(kicked.native_state(), &before);
        let (next, report) = kicked
            .contact_drift_candidate(&kick, controls(0.004), cx)
            .unwrap();
        near(report.gravity_splitting_work_j, -0.0025, 2e-12);
        near(report.integration_defect_j, -0.00125, 2e-12);
        near(report.absolute_integration_defect_j, 0.00375, 2e-12);
        let mut edited = kick;
        edited.duration_s = 10.0;
        edited.gravity_impulse_work_j = 1e12;
        edited.gravity_potential_before_j = None;
        edited.source_impulse_work_j = -1e12;
        let (same, repeated) = kicked
            .contact_drift_candidate(&edited, controls(0.004), cx)
            .unwrap();
        assert_eq!(same.native_state(), next.native_state());
        near(repeated.gravity_work_j, report.gravity_work_j, 0.0);
        near(
            repeated.absolute_integration_defect_j,
            report.absolute_integration_defect_j,
            0.0,
        );
        assert_eq!(kicked.native_state(), &before);
    });
}

#[test]
fn gravity_split_cancellation_at_kick_or_endpoint_retains_accepted_state() {
    for at_endpoint in [false, true] {
        with_cx(|gate, cx| {
            let material = bind(&source_card(6.0, 1.0, false, false));
            let source = Cancelling {
                inner: uniform([0.0, -1.0, 0.0]),
                gate,
                armed: Cell::new(false),
            };
            let body = FiniteSolidBody::new_with_gravity(
                problem(&material, &[], &[]),
                field(&source),
                &[[0.0; 3]; 4],
                &[[0.0; 3]; 4],
                0.0,
                Default::default(),
                cx,
            )
            .unwrap();
            let (body, _) = body
                .step_candidate(0.02, &[], [0.0; 3], controls(1e-9), cx)
                .unwrap();
            let before = body.native_state().clone();
            if at_endpoint {
                let (kicked, kick) = body
                    .force_kick_candidate(0.01, &[], [0.0; 3], [0.0; 3], Default::default(), cx)
                    .unwrap();
                let before_drift = kicked.native_state().clone();
                source.armed.set(true);
                let error = kicked
                    .contact_drift_candidate(&kick, controls(1e-3), cx)
                    .unwrap_err();
                assert!(error.to_string().to_ascii_lowercase().contains("cancel"));
                assert_eq!(kicked.native_state(), &before_drift);
            } else {
                source.armed.set(true);
                let error = body
                    .force_kick_candidate(0.01, &[], [0.0; 3], [0.0; 3], Default::default(), cx)
                    .unwrap_err();
                assert!(error.to_string().to_ascii_lowercase().contains("cancel"));
            }
            assert_eq!(body.native_state(), &before);
            assert_eq!(body.native_state().history.len(), 1);
        });
    }
}
