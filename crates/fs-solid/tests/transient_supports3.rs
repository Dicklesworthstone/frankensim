//! G1/G3/G4: prescribed support inertia, coupled oscillation, work and refusal.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_solid::linear3::*;
use fs_solid::loads3::transient::adjoint::TransientAdjointConfig3;
use fs_solid::loads3::transient::{TransientConfig3, TransientError3};
use fs_solid::loads3::{TetLoadGeometry3, VolumeLoad3};

const NODES: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const ZERO: [[f64; 3]; 4] = [[0.0; 3]; 4];
const FRAME: ContentHash = ContentHash([79; 32]);

fn with_cx<R>(cancelled: bool, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    if cancelled {
        gate.request();
    }
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 79,
                kernel_id: 79,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}
fn material() -> TetElasticMaterial {
    TetElasticMaterial::try_new(6.0, 60.0, 0.0, [79; 32]).unwrap()
}
fn problem<'a>(m: &'a TetElasticMaterial, fixed: &'a [usize]) -> TetLinearElasticProblem<'a> {
    TetLinearElasticProblem {
        nodes_m: &NODES,
        tetrahedra: &TETS,
        materials: TetMaterialField::Uniform(m),
        fixed_dofs: fixed,
        budget: TetAssemblyBudget::default(),
    }
}
fn geometry(cx: &Cx<'_>) -> TetLoadGeometry3<'static> {
    TetLoadGeometry3::new(&NODES, &TETS, TetAssemblyBudget::default(), 4, cx).unwrap()
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-11, "{a:e} != {b:e}");
}

#[test]
fn g1_initial_acceleration_includes_mass_stiffness_and_both_damping_cross_terms() {
    with_cx(false, |cx| {
        let m = material();
        let fixed: Vec<_> = (0..12).filter(|d| *d != 3).collect();
        let p = problem(&m, &fixed);
        let g = geometry(cx);
        let load = g.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let mut u = ZERO;
        u[0][0] = 0.0005;
        u[1][0] = 0.0003;
        let mut v = ZERO;
        v[0][0] = 0.002;
        v[1][0] = -0.001;
        let mut a = ZERO;
        a[0][0] = 0.004;
        let config = TransientConfig3 {
            mass_damping_per_s: 0.4,
            stiffness_damping_s: 0.006,
            ..Default::default()
        };
        let mut run = load
            .begin_transient_with_supports(&p, None, 0.0, &u, &v, &a, config, cx)
            .unwrap();
        // V=1/6, E=60, rho=6: M_ff=.1, M_fc=.05, K_ff=10, K_fc=-10.
        // C_ff=.4*.1+.006*10=.1; C_fc=.4*.05+.006*(-10)=-.04.
        let expected = (-10.0 * u[1][0] + 10.0 * u[0][0] - 0.1 * v[1][0] + 0.04 * v[0][0]
            - 0.05 * a[0][0])
            / 0.1;
        close(run.state().acceleration[0], expected);
        let full = run.nodal_motion(cx).unwrap();
        assert_eq!(full[0], u);
        assert_eq!(full[1], v);
        close(full[2][0][0], a[0][0]);
        close(full[2][1][0], expected);
        let response = run.recover_response(FRAME, 0.01, cx).unwrap();
        assert!(response.free_force_residual_n < 1e-11);
        close(
            response.inertia_force_n[1][0],
            0.1 * expected + 0.05 * a[0][0],
        );
        close(
            response.damping_force_n[1][0],
            0.1 * v[1][0] - 0.04 * v[0][0],
        );
        close(
            response.kinetic_plus_elastic_energy_j,
            run.state().total_energy_j,
        );
        let energy = run.state().total_energy_j;
        let loaded = g
            .assemble(&[], VolumeLoad3::None, &[(0, 0.034), (3, 0.012)], cx)
            .unwrap();
        let jump = run.apply_load_jump(&loaded, cx).unwrap();
        close(jump.acceleration_change_l2_m_s2, 0.012 / 0.1);
        close(run.state().acceleration[0], expected + 0.012 / 0.1);
        assert_eq!(run.state().time_s, 0.0);
        assert_eq!(run.state().step, 0);
        assert_eq!(run.state().total_energy_j, energy);
        let jumped = run.nodal_motion(cx).unwrap();
        assert_eq!(jumped[0], u);
        assert_eq!(jumped[1], v);
        assert_eq!(jumped[2][0], a[0]);
        let right = run.recover_response(FRAME, 0.01, cx).unwrap();
        // Reaction change includes M_cf*delta(a_f) and the newly applied
        // fixed-node force, not merely minus the external force.
        close(
            right.support_reactions_n[0][0] - response.support_reactions_n[0][0],
            0.05 * (0.012 / 0.1) - 0.034,
        );
        let h = 0.01;
        let step = run
            .advance_with_supports(
                &loaded,
                &[(0, u[0][0] + h * v[0][0] + 0.5 * h * h * a[0][0])],
                h,
                cx,
            )
            .unwrap();
        close(step.balance_error_j, 0.0);
    });
}

#[test]
fn g1_accelerating_base_matches_an_independent_oscillator_with_second_order_error() {
    with_cx(false, |cx| {
        let m = material();
        let fixed: Vec<_> = (0..12).filter(|d| *d != 3).collect();
        let p = problem(&m, &fixed);
        let g = geometry(cx);
        let load = g.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let acceleration = 0.02;
        let mut a = ZERO;
        a[0][0] = acceleration;
        let mut errors = Vec::new();
        for steps in [25, 50, 100] {
            let mut run = load
                .begin_transient_with_supports(
                    &p,
                    None,
                    0.0,
                    &ZERO,
                    &ZERO,
                    &a,
                    TransientConfig3::default(),
                    cx,
                )
                .unwrap();
            close(run.state().acceleration[0], -0.5 * acceleration);
            let mut work = 0.0;
            for i in 1..=steps {
                let time = 0.5 * i as f64 / steps as f64;
                let report = run
                    .advance_with_supports(
                        &load,
                        &[(0, 0.5 * acceleration * time * time)],
                        0.5 / steps as f64,
                        cx,
                    )
                    .unwrap();
                assert!(report.balance_error_j.abs() < 1e-11);
                close(report.mechanical_work_j, 0.0);
                work += report.support_work_j;
            }
            let time = run.state().time_s;
            // .1 q''+10 q=10 s-.05 s'', s=.5*a*t², q(0)=q'(0)=0.
            let exact = 0.5 * acceleration * time * time
                - 0.015 * acceleration * (1.0 - (10.0 * time).cos());
            errors.push((run.state().displacement[0] - exact).abs());
            close(work, run.state().total_energy_j);
            let full = run.nodal_motion(cx).unwrap();
            close(full[1][0][0], acceleration * time);
            // Recovering acceleration from rounded displacement endpoints
            // amplifies their rounding by h^-2, without degrading q's order.
            assert!((full[2][0][0] - acceleration).abs() < 5e-11);
            assert!(
                run.recover_response(FRAME, 0.01, cx)
                    .unwrap()
                    .free_force_residual_n
                    < 1e-10
            );
        }
        for ratio in [errors[0] / errors[1], errors[1] / errors[2]] {
            assert!((3.7..4.3).contains(&ratio), "{errors:?}");
        }
        assert!(errors[2] < 5e-7);
    });
}

#[test]
fn g1_direct_acceleration_drives_the_complete_partitioned_equation() {
    with_cx(false, |cx| {
        let m = material();
        let fixed: Vec<_> = (0..12).filter(|d| *d != 3).collect();
        let p = problem(&m, &fixed);
        let g = geometry(cx);
        let load = g.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let mut u = ZERO;
        u[0][0] = 0.0005;
        u[1][0] = 0.0003;
        let mut v = ZERO;
        v[0][0] = 0.002;
        v[1][0] = -0.001;
        let mut a = ZERO;
        a[0][0] = 0.004;
        let config = TransientConfig3 {
            mass_damping_per_s: 0.4,
            stiffness_damping_s: 0.006,
            ..Default::default()
        };
        let mut run = load
            .begin_transient_with_supports(&p, None, 0.0, &u, &v, &a, config, cx)
            .unwrap();
        let mut support = [u[0][0], v[0][0], a[0][0]];
        let mut free = [u[1][0], v[1][0], run.state().acceleration[0]];
        let h = 0.01;
        let b = h * h / 4.0;
        let mut total_work = 0.0;
        let initial_energy = run.state().total_energy_j;
        for acceleration in [0.02, -0.01, 0.03, 0.0] {
            // Independent scalar partition: M_ff=.1, M_fc=.05,
            // K_ff=10, K_fc=-10, C_ff=.1 and C_fc=-.04.
            let support_velocity = support[1] + 0.5 * h * (support[2] + acceleration);
            let support_displacement = support[0] + 0.5 * h * (support[1] + support_velocity);
            let predictor = free[0] + h * free[1] + b * free[2];
            let velocity_predictor = free[1] + 0.5 * h * free[2];
            let force = 10.0 * support_displacement + 0.04 * support_velocity - 0.05 * acceleration;
            let free_acceleration =
                (force - 10.0 * predictor - 0.1 * velocity_predictor) / (0.1 + 0.05 * h + 10.0 * b);
            free = [
                predictor + b * free_acceleration,
                velocity_predictor + 0.5 * h * free_acceleration,
                free_acceleration,
            ];
            support = [support_displacement, support_velocity, acceleration];
            let step = run
                .advance_with_support_accelerations(&load, &[(0, acceleration)], h, cx)
                .unwrap();
            let full = run.nodal_motion(cx).unwrap();
            for field in 0..3 {
                close(full[field][0][0], support[field]);
                close(full[field][1][0], free[field]);
            }
            assert_eq!(full[2][0][0], acceleration);
            assert!(
                run.recover_response(FRAME, 0.01, cx)
                    .unwrap()
                    .free_force_residual_n
                    < 1e-11
            );
            close(step.balance_error_j, 0.0);
            total_work += step.support_work_j - step.damping_loss_j;
        }
        close(run.state().total_energy_j - initial_energy, total_work);
    });
}

#[test]
fn g1_direct_support_acceleration_and_work_survive_unresolved_displacement() {
    with_cx(false, |cx| {
        let m = material();
        let fixed: Vec<_> = (0..12).collect();
        let p = problem(&m, &fixed);
        let g = geometry(cx);
        let load = g.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let u = [[1.0 / 1024.0, 0.0, 0.0]; 4];
        let acceleration = [[1.0, 0.0, 0.0]; 4];
        let config = TransientConfig3 {
            absolute_energy_tolerance_j: 1e-30,
            ..Default::default()
        };
        let mut run = load
            .begin_transient_with_supports(&p, None, 0.0, &u, &ZERO, &acceleration, config, cx)
            .unwrap();
        let h = 2.0_f64.powi(-32);
        let step = run
            .advance_with_support_accelerations(
                &load,
                &[(0, 1.0), (3, 1.0), (6, 1.0), (9, 1.0)],
                h,
                cx,
            )
            .unwrap();
        let full = run.nodal_motion(cx).unwrap();
        assert_eq!(
            full[0], u,
            "the absolute displacement cannot represent this increment"
        );
        assert_eq!(full[1], [[h, 0.0, 0.0]; 4]);
        assert_eq!(full[2], acceleration);
        // rho*V = 1 kg: all support work is the exact rigid kinetic gain.
        let energy = 0.5 * h * h;
        assert!((step.support_work_j / energy - 1.0).abs() < 1e-12);
        assert!((run.state().kinetic_energy_j / energy - 1.0).abs() < 1e-12);
        assert!(step.balance_error_j.abs() < 1e-30);
        assert_eq!(step.iterations, 0);
    });
}

#[test]
fn g1_fixed_displacement_endpoint_retains_resolved_predictor_drift() {
    with_cx(false, |cx| {
        let m = material();
        let fixed: Vec<_> = (0..12).collect();
        let p = problem(&m, &fixed);
        let g = geometry(cx);
        let load = g.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let displacement = 2.0_f64.powi(-10);
        let h = 2.0_f64.powi(-32);
        let u = [[displacement, 0.0, 0.0]; 4];
        let velocity = [[h, 0.0, 0.0]; 4];
        let config = TransientConfig3 {
            absolute_energy_tolerance_j: 1e-30,
            ..Default::default()
        };
        let mut run = load
            .begin_transient_with_supports(&p, None, 0.0, &u, &velocity, &ZERO, config, cx)
            .unwrap();
        let initial_energy = run.state().total_energy_j;
        let step = run
            .advance_with_supports(
                &load,
                &[
                    (0, displacement),
                    (3, displacement),
                    (6, displacement),
                    (9, displacement),
                ],
                h,
                cx,
            )
            .unwrap();
        let full = run.nodal_motion(cx).unwrap();
        assert_eq!(full[0], u);
        // q1=q0 requires zero midpoint velocity. Although h*v0 is below
        // one displacement ULP, a1=-4*v0/h remains exactly representable.
        assert_eq!(full[1], [[-h, 0.0, 0.0]; 4]);
        assert_eq!(full[2], [[-4.0, 0.0, 0.0]; 4]);
        assert_eq!(step.support_work_j, 0.0);
        assert_eq!(run.state().total_energy_j, initial_energy);
        assert!(step.balance_error_j.abs() < 1e-30);
    });
}

#[test]
fn g4_direct_support_acceleration_refusals_preserve_all_motion() {
    with_cx(false, |cx| {
        let m = material();
        let fixed: Vec<_> = (0..12).filter(|d| *d != 3).collect();
        let p = problem(&m, &fixed);
        let g = geometry(cx);
        let load = g.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let mut run = load
            .begin_transient_with_supports(
                &p,
                None,
                0.0,
                &ZERO,
                &ZERO,
                &ZERO,
                Default::default(),
                cx,
            )
            .unwrap();
        let before = run.state().clone();
        let motion = run.nodal_motion(cx).unwrap();
        for values in [
            vec![(3, 1.0)],
            vec![(12, 1.0)],
            vec![(0, 1.0), (0, 2.0)],
            vec![(1, 1.0), (0, 1.0)],
            vec![(0, f64::NAN)],
            vec![(0, f64::INFINITY)],
            vec![(0, 1e20)],
        ] {
            assert!(
                run.advance_with_support_accelerations(&load, &values, 0.01, cx)
                    .is_err()
            );
            assert_eq!(*run.state(), before);
            assert_eq!(run.nodal_motion(cx).unwrap(), motion);
        }
        with_cx(true, |cancelled| {
            assert!(matches!(
                run.advance_with_support_accelerations(&load, &[(0, 0.01)], 0.01, cancelled),
                Err(TransientError3::Cancelled)
            ));
        });
        assert_eq!(*run.state(), before);
        assert_eq!(run.nodal_motion(cx).unwrap(), motion);
        run.advance_with_support_accelerations(&load, &[(0, 0.01)], 0.01, cx)
            .unwrap();
        let old = run.nodal_motion(cx).unwrap();
        run.advance_with_support_accelerations(&load, &[], 0.01, cx)
            .unwrap();
        let continued = run.nodal_motion(cx).unwrap();
        assert_eq!(continued[2][0][0], 0.0);
        assert!(continued[1][0][0] > old[1][0][0]);
        assert!(continued[0][0][0] > old[0][0][0]);
    });
}

#[test]
fn g1_fully_prescribed_deformation_closes_external_thermal_and_support_work() {
    with_cx(false, |cx| {
        let m = material();
        let fixed: Vec<_> = (0..12).collect();
        let p = problem(&m, &fixed);
        let g = geometry(cx);
        // Public spectral assembly still refuses a zero-dimensional eigenproblem.
        assert!(matches!(p.assemble(cx), Err(TetElasticError::NoFreeDofs)));
        let thermal = TetThermalStrainState::try_new(
            301.0,
            300.0,
            [0.001, 0.0, 0.0, 0.0, 0.0, 0.0],
            m.material_state_identity,
            ContentHash([80; 32]),
            0.1,
        )
        .unwrap();
        let load = g.assemble(&[], VolumeLoad3::None, &[(3, 0.2)], cx).unwrap();
        let speed = 1.0 / 1024.0;
        let mut velocity = ZERO;
        velocity[1][0] = speed;
        let mut run = load
            .begin_transient_with_supports(
                &p,
                Some(TetThermalStrainField::Uniform(&thermal)),
                0.0,
                &ZERO,
                &velocity,
                &ZERO,
                TransientConfig3::default(),
                cx,
            )
            .unwrap();
        assert!(run.free_dofs().is_empty());
        let initial = run.state().total_energy_j;
        for i in 1..=8 {
            let time = i as f64 / 32.0;
            let step = run
                .advance_with_supports(&load, &[(3, speed * time)], 1.0 / 32.0, cx)
                .unwrap();
            assert_eq!(step.iterations, 0);
            close(step.mechanical_work_j, 0.2 * speed / 32.0);
            close(step.thermal_equivalent_work_j, 0.01 * speed / 32.0);
            assert!(step.balance_error_j.abs() < 1e-12);
        }
        let q = speed * 0.25;
        close(run.state().kinetic_energy_j, 0.5 * 0.1 * speed * speed);
        close(run.state().total_energy_j - initial, 5.0 * q * q);
        let response = run.recover_response(FRAME, 0.01, cx).unwrap();
        close(response.support_reactions_n[1][0], 10.0 * (q - 0.001) - 0.2);
        close(
            response.kinetic_plus_elastic_energy_j,
            0.05 * speed * speed + 5.0 * (q - 0.001).powi(2),
        );
    });
}

#[test]
fn g4_invalid_supports_cancel_and_incompatible_operations_preserve_the_checkpoint() {
    with_cx(false, |cx| {
        let m = material();
        let fixed: Vec<_> = (0..12).filter(|d| *d != 3).collect();
        let p = problem(&m, &fixed);
        let g = geometry(cx);
        let load = g.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let mut bad_a = ZERO;
        bad_a[1][0] = 1.0;
        assert!(
            load.begin_transient_with_supports(
                &p,
                None,
                0.0,
                &ZERO,
                &ZERO,
                &bad_a,
                TransientConfig3::default(),
                cx
            )
            .is_err()
        );
        let mut run = load
            .begin_transient_with_supports(
                &p,
                None,
                0.0,
                &ZERO,
                &ZERO,
                &ZERO,
                TransientConfig3::default(),
                cx,
            )
            .unwrap();
        let before = run.state().clone();
        let motion = run.nodal_motion(cx).unwrap();
        for endpoint in [
            &[(3, 0.0)][..],
            &[(0, 0.0), (0, 0.0)],
            &[(12, 0.0)],
            &[(0, f64::NAN)],
            &[(0, 1.0)],
        ] {
            assert!(
                run.advance_with_supports(&load, endpoint, 0.01, cx)
                    .is_err()
            );
        }
        assert!(run.advance(&load, 0.01, cx).is_err());
        assert_eq!(
            run.apply_load_jump(&load, cx)
                .unwrap()
                .acceleration_change_l2_m_s2,
            0.0
        );
        assert!(run.adjoint(TransientAdjointConfig3::default(), cx).is_err());
        with_cx(true, |cancelled| {
            assert_eq!(
                run.advance_with_supports(&load, &[(0, 0.00001)], 0.01, cancelled),
                Err(TransientError3::Cancelled)
            )
        });
        assert_eq!(run.state(), &before);
        assert_eq!(run.nodal_motion(cx).unwrap(), motion);
        let mut fork = run.clone();
        for i in 1..=5 {
            let endpoint = [(0, i as f64 / 1_000_000.0)];
            assert_eq!(
                run.advance_with_supports(&load, &endpoint, 0.01, cx)
                    .unwrap(),
                fork.advance_with_supports(&load, &endpoint, 0.01, cx)
                    .unwrap()
            );
            assert_eq!(
                run.nodal_motion(cx).unwrap(),
                fork.nodal_motion(cx).unwrap()
            );
        }
    });
}
