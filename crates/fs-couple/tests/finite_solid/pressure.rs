//! Three-axis breathing under actual follower pressure. The independent
//! scalar equation uses the analytic volumetric Ogden energy on F=lambda I.
use super::*;
use fs_solid::hyper3::pressure::{ClosedTetPressureSurface, TetPressure};
use std::sync::Arc;

fn breathing_supports() -> Vec<(usize, f64)> {
    (0..12)
        .filter(|dof| ![3, 7, 11].contains(dof))
        .map(|dof| (dof, 0.0))
        .collect()
}

fn breathing_state(stretch: f64, velocity: f64) -> ([[f64; 3]; 4], [[f64; 3]; 4]) {
    let mut u = [[0.0; 3]; 4];
    let mut v = [[0.0; 3]; 4];
    for i in 0..3 {
        u[i + 1][i] = stretch - 1.0;
        v[i + 1][i] = velocity;
    }
    (u, v)
}

fn breathing_oracle(stretch: f64, velocity: f64, h: f64, pressure: f64) -> (f64, f64) {
    let mass_per_coordinate = 6.0 * VOLUME / 10.0;
    let residual = |next: f64| {
        let mid = 0.5 * (stretch + next);
        2.0 * mass_per_coordinate * (next - stretch - h * velocity) / (h * h)
            + VOLUME * (40.0 * (mid.powi(3) - 1.0) + pressure) * mid * mid
    };
    let (mut low, mut high) = (0.6, 1.5);
    assert!(residual(low) < 0.0 && residual(high) > 0.0);
    for _ in 0..100 {
        let middle = 0.5 * (low + high);
        if residual(middle) < 0.0 {
            low = middle;
        } else {
            high = middle;
        }
    }
    let next = 0.5 * (low + high);
    (next, 2.0 * (next - stretch) / h - velocity)
}

#[test]
fn g1_pressure_changes_acceleration_motion_and_exact_work_on_existing_mass() {
    with_cx(|_, cx| {
        let card = source_card(6.0, 1.0, false, false);
        let material = bind(&card);
        let fixed = breathing_supports();
        let p = problem(&material, &fixed, &[]);
        let surface = Arc::new(
            ClosedTetPressureSurface::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx)
                .unwrap(),
        );
        let (u, v) = breathing_state(1.1, 0.05);
        let initial = FiniteSolidBody::new_with_pressure(
            p,
            surface,
            0.8,
            &u,
            &v,
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        let mass = 6.0 * VOLUME / 10.0;
        let acceleration = -VOLUME * (40.0 * (1.1_f64.powi(3) - 1.0) + 0.8) * 1.1 * 1.1 / mass;
        for actual in &initial.native_state().a {
            near(*actual, acceleration, 2.0e-9);
        }
        let mut body = initial.clone();
        let (mut stretch, mut velocity) = (1.1, 0.05);
        for pressure in [0.8, 0.8, 1.6, 1.6] {
            let before = body.native_state().clone();
            let (expected_stretch, expected_velocity) =
                breathing_oracle(stretch, velocity, 0.01, pressure);
            let (next, report) = body
                .step_candidate_with_pressure(0.01, &[], [0.0; 3], pressure, controls(1.0e-3), cx)
                .unwrap();
            assert_eq!(body.native_state(), &before);
            for (q, v) in next.native_state().q.iter().zip(&next.native_state().v) {
                near(*q, expected_stretch - 1.0, 3.0e-9);
                near(*v, expected_velocity, 3.0e-8);
            }
            let work = -pressure * VOLUME * (expected_stretch.powi(3) - stretch.powi(3));
            near(report.pressure_work_j, work, 2.0e-10);
            near(report.source_work_j, work, 2.0e-10);
            let strain = |s: f64| VOLUME * 20.0 * (s.powi(3) - 1.0).powi(2);
            near(report.strain_after_j, strain(expected_stretch), 3.0e-9);
            near(
                report.kinetic_after_j,
                1.5 * mass * expected_velocity.powi(2),
                3.0e-9,
            );
            let defect = strain(expected_stretch) - strain(stretch)
                + 1.5 * mass * (expected_velocity.powi(2) - velocity.powi(2))
                - work;
            near(report.integration_defect_j, defect, 6.0e-9);
            // Fixed node zero receives both the pressure reaction and the
            // consistent mass off-diagonal contribution from the moving vertex.
            let actual_acceleration = (expected_velocity - velocity) / 0.01;
            let reaction = report
                .support_reactions_n
                .iter()
                .find(|(dof, _)| *dof == 0)
                .unwrap()
                .1;
            near(reaction, 1.5 * mass * actual_acceleration, 1.0e-7);
            near(next.reference_mass_kg(), initial.reference_mass_kg(), 0.0);
            let (replay, replay_report) = body
                .step_candidate_with_pressure(0.01, &[], [0.0; 3], pressure, controls(1.0e-3), cx)
                .unwrap();
            assert_eq!(replay.native_state(), next.native_state());
            assert_eq!(
                replay_report.pressure_work_j.to_bits(),
                report.pressure_work_j.to_bits()
            );
            (stretch, velocity) = (expected_stretch, expected_velocity);
            body = next;
        }
        assert_eq!(initial.native_state().steps, 0);
        assert_eq!(body.native_state().steps, 4);
    });
}

#[test]
fn g3_zero_pressure_recovers_original_dynamics_and_signed_pressure_reverses_initial_force() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let fixed = breathing_supports();
        let p = problem(&material, &fixed, &[]);
        let surface = Arc::new(
            ClosedTetPressureSurface::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx)
                .unwrap(),
        );
        let (u, v) = breathing_state(1.05, 0.01);
        let old = FiniteSolidBody::new(p, &u, &v, 0.0, Default::default(), cx).unwrap();
        let pressured = FiniteSolidBody::new_with_pressure(
            p,
            surface.clone(),
            0.0,
            &u,
            &v,
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        assert_eq!(old.native_state(), pressured.native_state());
        let (old_next, old_report) = old
            .step_candidate(0.01, &[], [0.0; 3], controls(1.0e-3), cx)
            .unwrap();
        let (next, report) = pressured
            .step_candidate_with_pressure(0.01, &[], [0.0; 3], 0.0, controls(1.0e-3), cx)
            .unwrap();
        assert_eq!(old_next.native_state(), next.native_state());
        assert_eq!(
            old_report.integration_defect_j.to_bits(),
            report.integration_defect_j.to_bits()
        );
        assert_eq!(report.pressure_work_j, 0.0);
        let (u, v) = breathing_state(1.0, 0.0);
        let inward = FiniteSolidBody::new_with_pressure(
            p,
            surface.clone(),
            0.7,
            &u,
            &v,
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        let outward = FiniteSolidBody::new_with_pressure(
            p,
            surface,
            -0.7,
            &u,
            &v,
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        for (a, b) in inward
            .native_state()
            .a
            .iter()
            .zip(&outward.native_state().a)
        {
            assert!(*a < 0.0 && *b > 0.0);
            near(*a, -*b, 1.0e-13);
        }
    });
}

#[test]
fn g4_pressure_binding_energy_and_cancel_refusals_leave_the_accepted_body_unchanged() {
    with_cx(|gate, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let fixed = breathing_supports();
        let p = problem(&material, &fixed, &[]);
        let surface = Arc::new(
            ClosedTetPressureSurface::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx)
                .unwrap(),
        );
        let changed_nodes = NODES.map(|node| node.map(|value| 2.0 * value));
        let unrelated = Arc::new(
            ClosedTetPressureSurface::new(
                &changed_nodes,
                &TETS,
                TetAssemblyBudget::standard(),
                4,
                cx,
            )
            .unwrap(),
        );
        let (u, v) = breathing_state(1.1, 0.05);
        assert!(
            FiniteSolidBody::new_with_pressure(
                p,
                unrelated,
                0.8,
                &u,
                &v,
                0.0,
                Default::default(),
                cx,
            )
            .is_err()
        );
        let body = FiniteSolidBody::new_with_pressure(
            p,
            surface,
            0.8,
            &u,
            &v,
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        let before = body.native_state().clone();
        assert!(
            body.step_candidate(0.01, &[], [0.0; 3], controls(1.0e-3), cx)
                .is_err(),
            "a pressure binding cannot silently disappear through the original API"
        );
        assert!(
            body.step_candidate_with_pressure(0.01, &[], [0.0; 3], f64::NAN, controls(1.0e-3), cx,)
                .is_err()
        );
        assert!(matches!(
            body.step_candidate_with_pressure(0.01, &[], [0.0; 3], 0.8, controls(1.0e-18), cx,),
            Err(FiniteSolidError::Energy { .. })
        ));
        assert_eq!(body.native_state(), &before);
        gate.request();
        assert!(
            body.step_candidate_with_pressure(0.01, &[], [0.0; 3], 0.8, controls(1.0e-3), cx,)
                .is_err()
        );
        assert_eq!(body.native_state(), &before);
    });
}

#[test]
fn g1_local_pressure_advances_the_original_material_body_with_independent_work() {
    with_cx(|gate, cx| {
        let fixed: Vec<_> = (0..12).filter(|dof| *dof != 3).map(|dof| (dof, 0.0)).collect();
        let surface = Arc::new(ClosedTetPressureSurface::new(
            &NODES, &TETS, TetAssemblyBudget::standard(), 4, cx,
        ).unwrap());
        let mut pressures = vec![0.0; surface.faces().len()];
        let sloping = surface.faces().iter().position(|face| !face.contains(&0)).unwrap();
        pressures[sloping] = 0.9;
        let pressure = TetPressure::PerFace(&pressures);
        let mut displacement = [[0.0; 3]; 4];
        let mut velocity = [[0.0; 3]; 4];
        displacement[1][0] = 0.1;
        velocity[1][0] = 0.04;
        let mut observed = Vec::new();
        for (rho, scale) in [(6.0, 1.0), (9.0, 2.0)] {
            let material = bind(&source_card(rho, scale, false, false));
            let p = problem(&material, &fixed, &[]);
            let body = FiniteSolidBody::new_with_pressure(
                p, surface.clone(), pressure, &displacement, &velocity, 0.0,
                Default::default(), cx,
            ).unwrap();
            // The projected x area of face (1,2,3) is exactly 1/2 regardless
            // of this axial stretch. Each of its three vertices receives -p/6.
            let force = -0.9 / 6.0;
            let mass = rho * VOLUME / 10.0;
            near(body.native_state().a[0],
                (force - VOLUME * nominal_stress(1.1, scale)) / mass, 2.0e-9);
            let before = body.native_state().clone();
            let (expected, expected_velocity) = midpoint_oracle(1.1, 0.04, rho, scale, 0.01, force);
            let (next, report) = body.step_candidate_with_pressure(
                0.01, &[], [0.0; 3], pressure, controls(1.0e-4), cx,
            ).unwrap();
            near(next.native_state().q[0], expected - 1.0, 2.0e-10);
            near(next.native_state().v[0], expected_velocity, 3.0e-8);
            near(report.pressure_work_j, force * (expected - 1.1), 1.0e-12);
            near(report.source_work_j, report.pressure_work_j, 0.0);
            near(report.integration_defect_j,
                strain_energy(expected, scale) - strain_energy(1.1, scale)
                    + 0.5 * mass * (expected_velocity.powi(2) - 0.04_f64.powi(2))
                    - force * (expected - 1.1), 1.0e-10);
            assert_eq!(body.native_state(), &before);
            let invalid = TetPressure::PerFace(&pressures[..3]);
            assert!(body.step_candidate_with_pressure(
                0.01, &[], [0.0; 3], invalid, controls(1.0e-4), cx,
            ).is_err());
            assert!(matches!(body.step_candidate_with_pressure(
                0.01, &[], [0.0; 3], pressure, controls(1.0e-18), cx,
            ), Err(FiniteSolidError::Energy { .. })));
            assert_eq!(body.native_state(), &before);
            let (retry, retry_report) = body.step_candidate_with_pressure(
                0.01, &[], [0.0; 3], pressure, controls(1.0e-4), cx,
            ).unwrap();
            assert_eq!(retry.native_state(), next.native_state());
            assert_eq!(retry_report.pressure_work_j.to_bits(), report.pressure_work_j.to_bits());
            observed.push(next.native_state().q[0]);
        }
        assert!((observed[0] - observed[1]).abs() > 1.0e-5);
        let material = bind(&source_card(6.0, 1.0, false, false));
        let body = FiniteSolidBody::new_with_pressure(
            problem(&material, &fixed, &[]), surface, pressure, &displacement,
            &velocity, 0.0, Default::default(), cx,
        ).unwrap();
        let before = body.native_state().clone();
        gate.request();
        assert!(body.step_candidate_with_pressure(
            0.01, &[], [0.0; 3], pressure, controls(1.0e-4), cx,
        ).is_err());
        assert_eq!(body.native_state(), &before);
    });
}

#[test]
fn g3_uniform_face_field_reproduces_uniform_motion_and_zero_field_reproduces_dead_loads() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let fixed = breathing_supports();
        let p = problem(&material, &fixed, &[]);
        let surface = Arc::new(ClosedTetPressureSurface::new(
            &NODES, &TETS, TetAssemblyBudget::standard(), 4, cx,
        ).unwrap());
        let (u, v) = breathing_state(1.05, 0.01);
        for value in [0.0, 0.8] {
            let field = vec![value; surface.faces().len()];
            let pressure = TetPressure::PerFace(&field);
            let body = FiniteSolidBody::new_with_pressure(
                p, surface.clone(), pressure, &u, &v, 0.0, Default::default(), cx,
            ).unwrap();
            let uniform = FiniteSolidBody::new_with_pressure(
                p, surface.clone(), value, &u, &v, 0.0, Default::default(), cx,
            ).unwrap();
            assert_eq!(body.native_state(), uniform.native_state());
            let (next, report) = body.step_candidate_with_pressure(
                0.01, &[], [0.0; 3], pressure, controls(1.0e-3), cx,
            ).unwrap();
            let (uniform_next, uniform_report) = uniform.step_candidate_with_pressure(
                0.01, &[], [0.0; 3], value, controls(1.0e-3), cx,
            ).unwrap();
            assert_eq!(next.native_state(), uniform_next.native_state());
            near(report.pressure_work_j, uniform_report.pressure_work_j, 2.0e-14);
            if value == 0.0 {
                let original = FiniteSolidBody::new(p, &u, &v, 0.0, Default::default(), cx).unwrap();
                let (old_next, _) = original.step_candidate(
                    0.01, &[], [0.0; 3], controls(1.0e-3), cx,
                ).unwrap();
                assert_eq!(next.native_state(), old_next.native_state());
            }
        }
    });
}
