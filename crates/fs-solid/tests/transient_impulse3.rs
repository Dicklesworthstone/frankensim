//! G1/G3/G4 physical impulses on the original consistent transient mass.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_solid::linear3::{TetElasticMaterial, TetLinearElasticProblem, TetMaterialField};
use fs_solid::loads3::transient::events::NodalImpulseControls3;
use fs_solid::loads3::transient::{TransientConfig3, TransientError3};
use fs_solid::loads3::{TetLoadGeometry3, VolumeLoad3};

const NODES: [[f64; 3]; 4] = [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const ZERO: [[f64; 3]; 4] = [[0.0; 3]; 4];

fn with_cx<R>(cancelled: bool, run: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    if cancelled {
        gate.request();
    }
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        run(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 41,
                kernel_id: 17,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn problem<'a>(
    material: &'a TetElasticMaterial,
    fixed: &'a [usize],
) -> TetLinearElasticProblem<'a> {
    TetLinearElasticProblem {
        nodes_m: &NODES,
        tetrahedra: &TETS,
        materials: TetMaterialField::Uniform(material),
        fixed_dofs: fixed,
        budget: Default::default(),
    }
}

fn config() -> TransientConfig3 {
    TransientConfig3 {
        maximum_displacement_m: 1.0,
        ..Default::default()
    }
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}

#[test]
fn g1_off_center_impulse_retains_consistent_inertia_linear_and_angular_momentum() {
    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(6.0, 60.0, 0.0, [41; 32]).unwrap();
        let p = problem(&material, &[]);
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, p.budget, 4, cx).unwrap();
        let empty = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let mut run = empty
            .begin_transient(&p, None, 2.0, &ZERO, &ZERO, config(), cx)
            .unwrap();
        let report = run
            .apply_nodal_impulse(&[(4, 0.2)], Default::default(), cx)
            .unwrap();
        // For total mass one, each Cartesian block is (I+11')/20.
        // Its inverse is 20*(I-11'/5), independent of any production matrix.
        let after = run.nodal_motion(cx).unwrap();
        for node in 0..4 {
            close(after[1][node][1], if node == 1 { 3.2 } else { -0.8 }, 2e-14);
            close(after[1][node][0], 0.0, 1e-15);
            close(after[1][node][2], 0.0, 1e-15);
        }
        let momentum = &report.nodal_momentum_change_kg_m_s;
        close(momentum.iter().map(|p| p[1]).sum::<f64>(), 0.2, 2e-15);
        let angular_z: f64 = NODES
            .iter()
            .zip(momentum)
            .map(|(x, p)| x[0] * p[1] - x[1] * p[0])
            .sum();
        close(angular_z, 0.2, 2e-15);
        assert_eq!(report.support_impulse_kg_m_s, ZERO);
        close(report.kinetic_change_j, 0.32, 2e-14);
        close(report.applied_impulse_work_j, 0.32, 2e-14);
        close(report.kinetic_energy_j[1], 0.32, 2e-14);
        assert_eq!(after[0], ZERO);
        assert_eq!(run.state().time_s, 2.0);
        assert_eq!(run.state().step, 0);
        assert!(std::ptr::eq(run.material_problem(), &p));
        // Native continuation starts with the actual impulsive velocity.
        run.advance(&empty, 0.001, cx).unwrap();
        assert_eq!(run.state().step, 1);
        assert_eq!(run.state().time_s, 2.001);
    });
}

#[test]
fn g1_supported_impulse_reports_cross_mass_reaction_and_refreshes_damped_acceleration() {
    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(6.0, 60.0, 0.0, [41; 32]).unwrap();
        let fixed: Vec<_> = (0..12).filter(|&d| d != 3).collect();
        let p = problem(&material, &fixed);
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, p.budget, 4, cx).unwrap();
        let load = geometry
            .assemble(&[], VolumeLoad3::None, &[(3, 0.2)], cx)
            .unwrap();
        let mut u = ZERO;
        u[1][0] = 0.01;
        let mut v = ZERO;
        v[1][0] = 0.2;
        let cfg = TransientConfig3 {
            mass_damping_per_s: 0.4,
            stiffness_damping_s: 0.02,
            ..config()
        };
        let mut run = load
            .begin_transient(&p, None, 0.0, &u, &v, cfg, cx)
            .unwrap();
        // Independent one-DOF values: m=.1, k=10, c=.4*.1+.02*10=.24.
        close(run.state().acceleration[0], 0.52, 1e-12);
        let report = run
            .apply_nodal_impulse(&[(0, 0.1), (3, 0.03)], Default::default(), cx)
            .unwrap();
        close(run.state().velocity[0], 0.5, 1e-14);
        close(run.state().acceleration[0], -0.2, 1e-12);
        assert_eq!(run.state().displacement, vec![0.01]);
        close(report.support_impulse_kg_m_s[0][0], -0.085, 1e-14);
        close(report.support_impulse_kg_m_s[2][0], 0.015, 1e-14);
        close(report.support_impulse_kg_m_s[3][0], 0.015, 1e-14);
        close(report.support_impulse_kg_m_s[1][0], 0.0, 1e-15);
        close(
            report
                .nodal_momentum_change_kg_m_s
                .iter()
                .map(|p| p[0])
                .sum::<f64>(),
            0.075,
            1e-14,
        );
        close(report.kinetic_change_j, 0.0105, 1e-14);
        assert_eq!(report.support_impulse_work_j, 0.0);
        close(report.work_residual_j, 0.0, 1e-14);
        let response = run
            .recover_response(fs_blake3::ContentHash([41; 32]), 0.1, cx)
            .unwrap();
        assert!(response.free_force_residual_n < 1e-11);
        close(response.external_force_n[1][0], 0.2, 0.0);
        run.advance(&load, 0.01, cx).unwrap();
        let qp = 0.01 + 0.01 * 0.5 - 0.000025 * 0.2;
        let vp = 0.5 - 0.005 * 0.2;
        let a = (0.2 - 10.0 * qp - 0.24 * vp) / (0.1 + 0.005 * 0.24 + 0.000025 * 10.0);
        close(run.state().acceleration[0], a, 1e-11);
    });
}

#[test]
fn g1_moving_support_impulse_work_uses_held_velocity_and_full_mass() {
    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(6.0, 60.0, 0.0, [41; 32]).unwrap();
        let fixed: Vec<_> = (0..12).filter(|&d| d != 3).collect();
        let p = problem(&material, &fixed);
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, p.budget, 4, cx).unwrap();
        let empty = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let velocity = [[1.0, 0.0, 0.0]; 4];
        let mut run = empty
            .begin_transient_with_supports(&p, None, 0.0, &ZERO, &velocity, &ZERO, config(), cx)
            .unwrap();
        let report = run
            .apply_nodal_impulse(&[(0, 0.04), (3, 0.02)], Default::default(), cx)
            .unwrap();
        let motion = run.nodal_motion(cx).unwrap();
        for node in [0, 2, 3] {
            assert_eq!(motion[1][node], velocity[node]);
        }
        close(motion[1][1][0], 1.2, 1e-14);
        assert_eq!(motion[0], ZERO);
        assert_eq!(motion[2], ZERO);
        close(report.kinetic_energy_j[0], 0.5, 1e-14);
        close(report.kinetic_energy_j[1], 0.552, 1e-14);
        close(report.kinetic_change_j, 0.052, 1e-14);
        close(report.applied_impulse_work_j, 0.062, 1e-14);
        close(report.support_impulse_work_j, -0.01, 1e-14);
        close(report.work_residual_j, 0.0, 1e-14);
        let supports: Vec<_> = fixed
            .iter()
            .map(|&d| (d, if d % 3 == 0 { 0.001 } else { 0.0 }))
            .collect();
        run.advance_with_supports(&empty, &supports, 0.001, cx)
            .unwrap();
        assert_eq!(run.state().step, 1);
    });
}

#[test]
fn g4_impulse_refusals_and_unrepresentable_velocity_leave_the_accepted_run_intact() {
    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(6.0, 60.0, 0.0, [41; 32]).unwrap();
        let p = problem(&material, &[]);
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, p.budget, 4, cx).unwrap();
        let empty = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let mut run = empty
            .begin_transient(&p, None, 0.0, &ZERO, &ZERO, config(), cx)
            .unwrap();
        let saved = run.state().clone();
        for bad in [
            vec![(1, 1.0), (0, 1.0)],
            vec![(0, 1.0), (0, 2.0)],
            vec![(12, 1.0)],
            vec![(0, f64::NAN)],
            vec![(0, f64::INFINITY)],
        ] {
            assert!(
                run.apply_nodal_impulse(&bad, Default::default(), cx)
                    .is_err()
            );
            assert_eq!(run.state(), &saved);
        }
        with_cx(true, |stopped| {
            assert!(matches!(
                run.apply_nodal_impulse(&[(0, 1.0)], Default::default(), stopped),
                Err(TransientError3::Cancelled)
            ));
        });
        assert_eq!(run.state(), &saved);
        let no_op = run
            .apply_nodal_impulse(&[], Default::default(), cx)
            .unwrap();
        assert_eq!(run.state(), &saved);
        assert_eq!(no_op.impulse_solve_iterations, 0);
        assert_eq!(no_op.kinetic_change_j, 0.0);
        // A solved increment that is swallowed when added to the accepted
        // velocity must fail the physical residual, not report a fictitious J.
        let large_velocity = [[1e16, 0.0, 0.0]; 4];
        let mut moving = empty
            .begin_transient(&p, None, 0.0, &ZERO, &large_velocity, config(), cx)
            .unwrap();
        let original = moving.state().clone();
        let strict = NodalImpulseControls3 {
            absolute_momentum_tolerance_kg_m_s: 1e-20,
            relative_tolerance: 1e-12,
            ..Default::default()
        };
        assert!(matches!(
            moving.apply_nodal_impulse(&[(0, 0.01)], strict, cx),
            Err(TransientError3::Momentum { .. })
        ));
        assert_eq!(moving.state(), &original);
    });
}
