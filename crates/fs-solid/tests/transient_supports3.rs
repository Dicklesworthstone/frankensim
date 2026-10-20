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
        let run = load
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
        assert!(run.apply_load_jump(&load, cx).is_err());
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
