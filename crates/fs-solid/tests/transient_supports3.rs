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
fn unresolved_support_velocity_cannot_hide_its_work_in_background_kinetic_energy() {
    with_cx(false, |cx| {
        // Numerical stress fixture: a large mass moving at an ordinary speed.
        // The requested acceleration's velocity increment is below one ULP,
        // while its independent reaction work is about half a joule.
        let m = TetElasticMaterial::try_new(6.0 * 2.0_f64.powi(54), 60.0, 0.0, [79; 32]).unwrap();
        let fixed: Vec<_> = (0..12).collect();
        let p = problem(&m, &fixed);
        let g = geometry(cx);
        let load = g.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let velocities = [[1.0, 0.0, 0.0]; 4];
        let config = TransientConfig3 {
            absolute_energy_tolerance_j: 1e-10, relative_energy_tolerance: 1e-10,
            ..Default::default()
        };
        let mut run = load.begin_transient_with_supports(&p, None, 0.0,
            &ZERO, &velocities, &ZERO, config, cx).unwrap();
        let before = run.nodal_motion(cx).unwrap();
        let result = run.advance_with_support_accelerations(&load,
            &[(0, 1.0), (3, 1.0), (6, 1.0), (9, 1.0)], 2.0_f64.powi(-54), cx);
        match result {
            Err(TransientError3::Energy { defect_j, allowance_j }) => {
                assert!((defect_j + 0.5).abs() < 1e-14, "{defect_j}");
                assert!(allowance_j < 1e-8);
            }
            other => panic!("rounded support velocity must refuse unresolved work: {other:?}"),
        }
        assert_eq!(run.nodal_motion(cx).unwrap(), before);
        assert_eq!(run.state().step, 0);
        assert_eq!(run.state().time_s, 0.0);
    });
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


fn eigenstrain(m: &TetElasticMaterial, strain: f64) -> TetThermalStrainState {
    TetThermalStrainState::try_new(300.0 + 1000.0 * strain, 300.0,
        [strain, 0.0, 0.0, 0.0, 0.0, 0.0], m.material_state_identity, FRAME, 0.1).unwrap()
}

#[test]
fn g1_g3_changing_eigenstrain_and_supports_share_full_damped_equations_and_work() {
    with_cx(false, |cx| {
        let m = material();
        let refreshed = TetElasticMaterial::try_new(6.0, 60.0, 0.0, [80; 32]).unwrap();
        let fixed: Vec<_> = (0..12).filter(|d| *d != 3).collect();
        let p = problem(&m, &fixed);
        let next_problem = problem(&refreshed, &fixed);
        let g = geometry(cx);
        let load = g.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let (mut qc, mut qf, mut vc, mut vf, mut ac, mut e) =
            (0.0005, 0.0003, 0.002, -0.001, 0.004, 0.0002);
        // Full consistent mass and Rayleigh cross terms, independently reduced:
        // M_ff=.1, M_fc=.05, K_ff=10, K_fc=-10, C_ff=.1, C_fc=-.04.
        let mut af = (10.0 * e - 10.0 * qf + 10.0 * qc - 0.1 * vf
            + 0.04 * vc - 0.05 * ac) / 0.1;
        let mut u = ZERO; u[0][0] = qc; u[1][0] = qf;
        let mut v = ZERO; v[0][0] = vc; v[1][0] = vf;
        let mut a = ZERO; a[0][0] = ac;
        let thermal = eigenstrain(&m, e);
        let config = TransientConfig3 { mass_damping_per_s: 0.4,
            stiffness_damping_s: 0.006, ..Default::default() };
        let mut run = load.begin_transient_with_supports(&p,
            Some(TetThermalStrainField::Uniform(&thermal)), 0.0, &u, &v, &a, config, cx).unwrap();
        close(run.state().acceleration[0], af);
        for (next_e, next_ac) in [(0.0007, 0.01), (-0.0001, -0.02), (0.0004, 0.0)] {
            let h = 0.01; let b = 0.25 * h * h;
            let next_vc = vc + 0.5 * h * (ac + next_ac);
            let next_qc = qc + 0.5 * h * (vc + next_vc);
            let qp = qf + h * vf + b * af;
            let vp = vf + 0.5 * h * af;
            let next_af = (10.0 * next_e + 10.0 * next_qc + 0.04 * next_vc
                - 0.05 * next_ac - 10.0 * qp - 0.1 * vp)
                / (0.1 + 0.05 * h + 10.0 * b);
            let next_qf = qp + b * next_af;
            let next_vf = vp + 0.5 * h * next_af;
            let reaction = |q_c: f64, q_f: f64, v_c: f64, v_f: f64,
                a_c: f64, a_f: f64, strain: f64| {
                0.1 * a_c + 0.05 * a_f + 0.16 * v_c - 0.04 * v_f
                    + 20.0 * q_c - 10.0 * q_f + 10.0 * strain
            };
            let old_reaction = reaction(qc, qf, vc, vf, ac, af, e);
            let next_reaction = reaction(next_qc, next_qf, next_vc, next_vf,
                next_ac, next_af, next_e);
            let target = eigenstrain(&refreshed, next_e);
            let field = TetThermalStrainField::Uniform(&target);
            let checkpoint = run.clone();
            let (next, report) = run.advance_thermal_with_support_accelerations_candidate(
                &next_problem, field, &load, &[(0, next_ac)], h, FRAME, 0.1, cx).unwrap();
            let (replayed, replay_report) = checkpoint.advance_thermal_with_support_accelerations_candidate(
                &next_problem, field, &load, &[(0, next_ac)], h, FRAME, 0.1, cx).unwrap();
            assert_eq!(report, replay_report);
            assert_eq!(next.nodal_motion(cx).unwrap(), replayed.nodal_motion(cx).unwrap());
            assert_eq!(run.state(), checkpoint.state());
            let motion = next.nodal_motion(cx).unwrap();
            for (actual, expected) in [
                (motion[0][0][0], next_qc), (motion[1][0][0], next_vc),
                (motion[2][0][0], next_ac), (motion[0][1][0], next_qf),
                (motion[1][1][0], next_vf), (motion[2][1][0], next_af),
            ] { close(actual, expected); }
            let response = next.recover_response(FRAME, 0.1, cx).unwrap();
            close(response.support_reactions_n[0][0], next_reaction);
            close(response.stress.elements()[0].stress_mandel_pa()[0],
                60.0 * (next_qf - next_qc - next_e));
            assert!(response.free_force_residual_n < 1e-11);
            close(report.native.thermal_equivalent_work_j,
                5.0 * (e + next_e) * ((next_qf - qf) - (next_qc - qc)));
            close(report.native.support_work_j,
                0.5 * (old_reaction + next_reaction) * (next_qc - qc));
            assert!(report.native.damping_loss_j > 0.0);
            assert!(report.native.balance_error_j.abs() < 1e-12);
            assert!(report.physical_balance_error_j.abs() < 1e-12);
            let (position_driven, position_report) = run.advance_thermal_with_supports_candidate(
                &next_problem, field, &load, &[(0, next_qc)], h, FRAME, 0.1, cx).unwrap();
            let positional = position_driven.nodal_motion(cx).unwrap();
            for component in 0..3 { for node in [0, 1] {
                close(positional[component][node][0], motion[component][node][0]);
            }}
            assert!(position_report.physical_balance_error_j.abs() < 1e-12);
            run = next;
            (qc, qf, vc, vf, ac, af, e) =
                (next_qc, next_qf, next_vc, next_vf, next_ac, next_af, next_e);
        }
    });
}

#[test]
fn g1_fully_prescribed_heating_averages_thermal_forces_and_retains_the_new_reactions() {
    with_cx(false, |cx| {
        let m = material();
        let fixed: Vec<_> = (0..12).collect();
        let p = problem(&m, &fixed);
        let g = geometry(cx);
        let load = g.assemble(&[], VolumeLoad3::None, &[(3, 0.2)], cx).unwrap();
        let speed = 1.0 / 128.0; let h = 1.0 / 32.0; let q = speed * h;
        let mut v = ZERO; v[1][0] = speed;
        let run = load.begin_transient_with_supports(&p, None, 0.0,
            &ZERO, &v, &ZERO, Default::default(), cx).unwrap();
        let target = eigenstrain(&m, 0.001);
        let (mut next, report) = run.advance_thermal_with_supports_candidate(&p,
            TetThermalStrainField::Uniform(&target), &load, &[(3, q)], h, FRAME, 0.1, cx).unwrap();
        assert!(next.free_dofs().is_empty());
        assert_eq!(report.native.iterations, 0);
        close(report.native.thermal_equivalent_work_j, 5.0 * 0.001 * q);
        close(report.native.mechanical_work_j, 0.2 * q);
        close(report.native.support_work_j, (5.0 * (q - 0.001) - 0.2) * q);
        close(report.thermal_material_work_j, -5.0 * (q - 0.001) * 0.001);
        close(report.elastic_energy_change_j, 5.0 * (q - 0.001).powi(2));
        assert!(report.native.balance_error_j.abs() < 1e-12);
        assert!(report.physical_balance_error_j.abs() < 1e-12);
        let response = next.recover_response(FRAME, 0.1, cx).unwrap();
        close(response.support_reactions_n[1][0], 10.0 * (q - 0.001) - 0.2);
        // Ordinary support advancement consumes the refreshed full thermal load.
        let frozen = next.advance_with_support_accelerations(&load, &[], h, cx).unwrap();
        close(frozen.thermal_equivalent_work_j, 10.0 * 0.001 * q);
        close(frozen.support_work_j, (15.0 * q - 10.0 * 0.001 - 0.2) * q);
        assert!(frozen.balance_error_j.abs() < 1e-12);
        close(next.nodal_motion(cx).unwrap()[0][1][0], 2.0 * q);
        close(next.recover_response(FRAME, 0.1, cx).unwrap().support_reactions_n[1][0],
            10.0 * (2.0 * q - 0.001) - 0.2);
        assert_eq!(run.state().step, 0);
        assert_eq!(run.nodal_motion(cx).unwrap()[0], ZERO);
    });
}

#[test]
fn g4_thermal_support_refusals_leave_the_accepted_checkpoint_usable() {
    with_cx(false, |cx| {
        let m = material();
        let fixed: Vec<_> = (0..12).filter(|d| *d != 3).collect();
        let p = problem(&m, &fixed);
        let g = geometry(cx);
        let load = g.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let run = load.begin_transient_with_supports(&p, None, 0.0,
            &ZERO, &ZERO, &ZERO, TransientConfig3 { maximum_steps: 1,
                ..Default::default() }, cx).unwrap();
        let before = run.state().clone();
        let motion = run.nodal_motion(cx).unwrap();
        let response = run.recover_response(FRAME, 0.1, cx).unwrap();
        let target = eigenstrain(&m, 0.001);
        let field = TetThermalStrainField::Uniform(&target);
        for supports in [&[(3, 0.0)][..], &[(0, 0.0), (0, 0.0)], &[(0, f64::NAN)]] {
            assert!(run.advance_thermal_with_support_accelerations_candidate(
                &p, field, &load, supports, 0.01, FRAME, 0.1, cx).is_err());
        }
        for changed in [
            TetElasticMaterial::try_new(7.0, 60.0, 0.0, [82; 32]).unwrap(),
        ] {
            let changed_problem = problem(&changed, &fixed);
            let changed_thermal = eigenstrain(&changed, 0.001);
            assert!(matches!(run.advance_thermal_with_supports_candidate(
                &changed_problem, TetThermalStrainField::Uniform(&changed_thermal),
                &load, &[(0, 0.00001)], 0.01, FRAME, 0.1, cx),
                Err(TransientError3::Invalid(_))));
        }
        // Late physical strain refusal must discard an otherwise solved candidate.
        assert!(run.advance_thermal_with_supports_candidate(&p, field, &load,
            &[(0, 0.00001)], 0.01, FRAME, 1e-8, cx).is_err());
        assert!(run.advance_thermal_candidate(&p, field, &load, 0.01, FRAME, 0.1, cx).is_err());
        with_cx(true, |cancelled| {
            assert!(matches!(run.advance_thermal_with_support_accelerations_candidate(
                &p, field, &load, &[(0, 0.01)], 0.01, FRAME, 0.1, cancelled),
                Err(TransientError3::Cancelled)));
        });
        assert_eq!(run.state(), &before);
        assert_eq!(run.nodal_motion(cx).unwrap(), motion);
        assert_eq!(run.recover_response(FRAME, 0.1, cx).unwrap(), response);
        assert!(run.thermal_strain_field().is_none());
        let (next, report) = run.advance_thermal_with_support_accelerations_candidate(
            &p, field, &load, &[(0, 0.01)], 0.01, FRAME, 0.1, cx).unwrap();
        assert!(report.physical_balance_error_j.abs() < 1e-12);
        assert!(matches!(next.advance_thermal_with_support_accelerations_candidate(
            &p, field, &load, &[(0, 0.01)], 0.01, FRAME, 0.1, cx),
            Err(TransientError3::Budget)));
        let stationary = load.begin_transient(&p, None, 0.0,
            &ZERO, &ZERO, Default::default(), cx).unwrap();
        assert!(matches!(stationary.advance_thermal_with_supports_candidate(
            &p, field, &load, &[(0, 0.0)], 0.01, FRAME, 0.1, cx),
            Err(TransientError3::Invalid(_))));
    });
}

#[test]
fn changing_stiffness_and_supports_follow_the_full_passive_midpoint_equation() {
    with_cx(false, |cx| {
      for zeta in [0.0, 0.006] {
        let materials = [
            material(),
            TetElasticMaterial::try_new(6.0, 30.0, 0.0, [91; 32]).unwrap(),
            TetElasticMaterial::try_new(6.0, 108.0, 0.0, [92; 32]).unwrap(),
            TetElasticMaterial::try_new(6.0, 42.0, 0.0, [93; 32]).unwrap(),
        ];
        let fixed: Vec<_> = (0..12).filter(|d| *d != 3).collect();
        let problems: Vec<_> = materials.iter().map(|m| problem(m, &fixed)).collect();
        let g = geometry(cx);
        let load = g.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let eta = 0.4; let mass = 0.1; let cross_mass = 0.05;
        let (mut qc, mut qf, mut vc, mut vf, mut ac, mut e, mut k) =
            (0.0005, 0.0003, 0.002, -0.001, 0.004, 0.0002, 10.0);
        let mut u = ZERO; u[0][0] = qc; u[1][0] = qf;
        let mut v = ZERO; v[0][0] = vc; v[1][0] = vf;
        let mut a = ZERO; a[0][0] = ac;
        let config = TransientConfig3 { mass_damping_per_s: eta,
            stiffness_damping_s: zeta, ..Default::default() };
        let initial = eigenstrain(&materials[0], e);
        let mut run = load.begin_transient_with_supports(&problems[0],
            Some(TetThermalStrainField::Uniform(&initial)), 0.0, &u, &v, &a, config, cx).unwrap();
        for (index, (k1, e1, ac1)) in [(5.0, 0.0007, -0.02),
            (18.0, -0.0001, 0.01), (7.0, 0.0004, 0.0)].into_iter().enumerate() {
            let h = 0.02;
            let vc1 = vc + 0.5*h*(ac+ac1);
            let qc1 = qc + 0.5*h*(vc+vc1);
            let vcbar = 0.5*(vc+vc1);
            let acbar = 0.5*(ac+ac1);
            let kbar = 0.5*(k+k1);
            let cff = eta*mass + zeta*kbar;
            let cfc = eta*cross_mass - zeta*kbar;
            let ccc = eta*mass + 2.0*zeta*kbar;
            // Independent two-coordinate reduction of the original tet:
            // M=[[.1,.05],[.05,.1]], K=k*[[2,-1],[-1,1]].
            let rate = (0.5*(k*e+k1*e1) - cff*vf - kbar*qf - 0.5*h*k1*vf
                - cross_mass*acbar - cfc*vcbar + 0.5*(k*qc+k1*qc1))
                / (mass + 0.5*h*cff + 0.25*h*h*k1);
            let vf1 = vf + h*rate;
            let qf1 = qf + 0.5*h*(vf+vf1);
            let vfbar = 0.5*(vf+vf1);
            let af1 = (k1*e1 - k1*(qf1-qc1) - (eta*mass+zeta*k1)*vf1
                - (eta*cross_mass-zeta*k1)*vc1 - cross_mass*ac1) / mass;
            let reaction = mass*ac1 + cross_mass*af1
                + (eta*mass+2.0*zeta*k1)*vc1 + (eta*cross_mass-zeta*k1)*vf1
                + 2.0*k1*qc1 - k1*qf1 + k1*e1;
            let reaction_midpoint = mass*acbar + cross_mass*rate + ccc*vcbar + cfc*vfbar
                + 0.5*(2.0*k*qc-k*qf+2.0*k1*qc1-k1*qf1+k*e+k1*e1);
            let expected_loss = h*(ccc*vcbar*vcbar + 2.0*cfc*vcbar*vfbar + cff*vfbar*vfbar);
            let expected_parameter = 0.5*(k1-k)
                *(2.0*qc*qc1-qc*qf1-qf*qc1+qf*qf1);
            let expected_material = 0.5*(k1-k)
                *((qf-qc-e)*(qf1-qc1-e1)+qc*qc1);
            let target = eigenstrain(&materials[index+1], e1);
            let field = TetThermalStrainField::Uniform(&target);
            let checkpoint = run.clone();
            let (next, report) = run.advance_thermal_with_support_accelerations_candidate(
                &problems[index+1], field, &load, &[(0, ac1)], h, FRAME, 0.1, cx).unwrap();
            let (replayed, replay_report) = checkpoint.advance_thermal_with_support_accelerations_candidate(
                &problems[index+1], field, &load, &[(0, ac1)], h, FRAME, 0.1, cx).unwrap();
            assert_eq!(report, replay_report);
            assert_eq!(next.nodal_motion(cx).unwrap(), replayed.nodal_motion(cx).unwrap());
            assert_eq!(run.state(), checkpoint.state());
            let motion = next.nodal_motion(cx).unwrap();
            for (actual, expected) in [(motion[0][0][0], qc1), (motion[1][0][0], vc1),
                (motion[2][0][0], ac1), (motion[0][1][0], qf1),
                (motion[1][1][0], vf1), (motion[2][1][0], af1)] { close(actual, expected); }
            let response = next.recover_response(FRAME, 0.1, cx).unwrap();
            close(response.support_reactions_n[0][0], reaction);
            assert!(response.free_force_residual_n < 1e-11);
            close(report.native.support_work_j, reaction_midpoint*(qc1-qc));
            close(report.native.damping_loss_j, expected_loss);
            close(report.native.stiffness_parameter_work_j, expected_parameter);
            close(report.stiffness_material_work_j, expected_material);
            close(report.thermal_material_work_j,
                -0.5*(k*(qf-qc-e)+k1*(qf1-qc1-e1))*(e1-e));
            assert!(report.native.damping_loss_j > 0.0);
            assert!(report.native.balance_error_j.abs() < 1e-12);
            assert!(report.physical_balance_error_j.abs() < 1e-12);
            let (position_driven, position_report) = run.advance_thermal_with_supports_candidate(
                &problems[index+1], field, &load, &[(0, qc1)], h, FRAME, 0.1, cx).unwrap();
            let by_position = position_driven.nodal_motion(cx).unwrap();
            for component in 0..3 { for node in [0,1] {
                close(by_position[component][node][0], motion[component][node][0]);
            }}
            assert!(position_report.physical_balance_error_j.abs() < 1e-12);
            // The next ordinary support step starts from the physical endpoint
            // acceleration at the NEW tensor, not the preceding midpoint rate.
            let mut continued = next.clone();
            let vc2 = vc1 + h*ac1;
            let qc2 = qc1 + 0.5*h*(vc1+vc2);
            let qp = qf1+h*vf1+0.25*h*h*af1;
            let vp = vf1+0.5*h*af1;
            let af2 = (k1*e1+k1*qc2-(eta*cross_mass-zeta*k1)*vc2-cross_mass*ac1
                - k1*qp-(eta*mass+zeta*k1)*vp)
                /(mass+0.5*h*(eta*mass+zeta*k1)+0.25*h*h*k1);
            continued.advance_with_support_accelerations(&load, &[(0,ac1)], h, cx).unwrap();
            close(continued.state().acceleration[0], af2);
            close(continued.state().displacement[0], qp+0.25*h*h*af2);
            run = next;
            (qc,qf,vc,vf,ac,e,k) = (qc1,qf1,vc1,vf1,ac1,e1,k1);
        }
      }
    });
}

#[test]
fn fully_prescribed_stiffness_changes_retain_material_work_at_zero_boundary_work() {
    with_cx(false, |cx| {
        let m0 = material();
        let m1 = TetElasticMaterial::try_new(6.0, 120.0, 0.0, [94;32]).unwrap();
        let fixed: Vec<_> = (0..12).collect();
        let p0 = problem(&m0, &fixed); let p1 = problem(&m1, &fixed);
        let g = geometry(cx);
        let load = g.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let q = 0.001; let e0 = 0.0002; let e1 = 0.0004;
        let mut u = ZERO; u[1][0] = q;
        let old = eigenstrain(&m0, e0); let new = eigenstrain(&m1, e1);
        for zeta in [0.0,0.006] {
            let run = load.begin_transient_with_supports(&p0,
                Some(TetThermalStrainField::Uniform(&old)), 0.0, &u, &ZERO, &ZERO,
                TransientConfig3 { stiffness_damping_s: zeta, ..Default::default() }, cx).unwrap();
            let (next, report) = run.advance_thermal_with_supports_candidate(&p1,
                TetThermalStrainField::Uniform(&new), &load, &[(3,q)], 0.01, FRAME, 0.1, cx).unwrap();
            assert!(next.free_dofs().is_empty());
            assert_eq!(report.native.iterations, 0);
            assert_eq!(report.native.support_work_j, 0.0);
            assert_eq!(report.native.mechanical_work_j, 0.0);
            assert_eq!(report.native.damping_loss_j, 0.0);
            close(report.native.stiffness_parameter_work_j, 5.0*q*q);
            close(report.stiffness_material_work_j, 5.0*(q-e0)*(q-e1));
            close(report.thermal_material_work_j, -0.5*(10.0*(q-e0)+20.0*(q-e1))*(e1-e0));
            close(report.elastic_energy_change_j, 10.0*(q-e1).powi(2)-5.0*(q-e0).powi(2));
            assert!(report.native.balance_error_j.abs() < 1e-12);
            assert!(report.physical_balance_error_j.abs() < 1e-12);
            close(next.recover_response(FRAME, 0.1, cx).unwrap().support_reactions_n[1][0], 20.0*(q-e1));
        }
    });
}
