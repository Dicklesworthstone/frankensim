//! G1/G3/G4: typed displacement sources drive actual 3-D inertial supports.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ga::Vec3;
use fs_qty::{Dims, QtyAny};
use fs_scenario::payload::{
    OrientationParity, OutsideDomainPolicy, PayloadId, PayloadMeta, QuantityContract,
    ReferenceSemantics, SampleSource, TableInterpolation, VectorPayload,
};
use fs_scenario::{Frame, FrameId, FrameMotion, FrameTree, WORLD};
use fs_solid::linear3::*;
use fs_solid::loads3::history::supports::{DISPLACEMENT_DIMS, SupportMotion3};
use fs_solid::loads3::history::{TRACTION_DIMS, TetLoadHistory3, VolumeHistory3};
use fs_solid::loads3::transient::TransientConfig3;
use fs_solid::loads3::{TetLoadGeometry3, VolumeLoad3};

const NODES: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const ZERO: [[f64; 3]; 4] = [[0.0; 3]; 4];
const ALL: [usize; 4] = [0, 1, 2, 3];
const FIXED: [usize; 12] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11];

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
                seed: 75,
                kernel_id: 75,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn source(
    dims: Dims,
    frame: FrameId,
    a: [f64; 3],
    b: [f64; 3],
    outside: OutsideDomainPolicy,
) -> VectorPayload {
    let basis = if frame == WORLD {
        "basis/world-cartesian"
    } else {
        "basis/frame-cartesian"
    };
    VectorPayload::new(
        PayloadMeta::new(
            QuantityContract::Dimensions(dims),
            PayloadId::new(basis).unwrap(),
            frame,
            OrientationParity::Even,
            ReferenceSemantics::Continuous,
        )
        .unwrap(),
        SampleSource::table(
            vec![
                QtyAny::new(0.0, Dims([0, 0, 1, 0, 0, 0])),
                QtyAny::new(1.0, Dims([0, 0, 1, 0, 0, 0])),
            ],
            vec![
                a.map(|v| QtyAny::new(v, dims)).to_vec(),
                b.map(|v| QtyAny::new(v, dims)).to_vec(),
            ],
            TableInterpolation::Linear,
            outside,
        )
        .unwrap(),
    )
    .unwrap()
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
        budget: TetAssemblyBudget::default(),
    }
}

fn history(frames: &FrameTree) -> TetLoadHistory3<'_> {
    TetLoadHistory3 {
        surfaces: &[],
        volume: VolumeHistory3::Fixed(VolumeLoad3::None),
        nodal_forces_n: &[],
        frames,
    }
}

#[test]
fn g1_uniform_translation_keeps_support_velocity_and_accounts_for_damping_work() {
    with_cx(false, |cx| {
        // Total mass = rho*V = 1 kg. A rigid translation has zero elastic
        // energy; prescribed supports supply the eta*M*v damping force.
        let material = TetElasticMaterial::try_new(6.0, 60.0, 0.0, [75; 32]).unwrap();
        let p = problem(&material, &FIXED);
        let g = TetLoadGeometry3::new(&NODES, &TETS, p.budget, 4, cx).unwrap();
        let frames = FrameTree::new();
        let h = history(&frames);
        let rate = 1.0 / 1024.0;
        let motion = source(
            DISPLACEMENT_DIMS,
            WORLD,
            [0.0; 3],
            [rate, 0.0, 0.0],
            OutsideDomainPolicy::Refuse,
        );
        let rows = [SupportMotion3 {
            nodes: &ALL,
            components: [true; 3],
            displacement: &motion,
        }];
        let config = TransientConfig3 {
            mass_damping_per_s: 0.25,
            ..Default::default()
        };
        let mut run = h
            .begin_transient_with_support_history(
                &g,
                &p,
                None,
                &rows,
                0.0,
                &ZERO,
                &[[rate, 0.0, 0.0]; 4],
                &ZERO,
                config,
                cx,
            )
            .unwrap();
        let mut work = 0.0;
        for _ in 0..8 {
            let step = run
                .advance_history_with_supports(&h, &rows, 1.0 / 32.0, cx)
                .unwrap();
            assert!(step.balance_error_j.abs() < 1e-12);
            assert!((step.support_work_j - step.damping_loss_j).abs() < 1e-12);
            work += step.support_work_j;
        }
        let state = run.nodal_motion(cx).unwrap();
        for n in 0..4 {
            assert!((state[0][n][0] - rate * 0.25).abs() < 1e-13);
            assert!((state[1][n][0] - rate).abs() < 1e-13);
            assert!(state[2][n][0].abs() < 1e-12);
        }
        assert!((run.state().kinetic_energy_j - 0.5 * rate * rate).abs() < 1e-12);
        assert!((work - 0.25 * rate * rate * 0.25).abs() < 1e-12);
    });
}

#[test]
fn g3_typed_force_and_support_samples_match_manual_coupled_endpoints() {
    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(6.0, 60.0, 0.0, [75; 32]).unwrap();
        let fixed: Vec<_> = FIXED.into_iter().filter(|d| *d != 3).collect();
        let p = problem(&material, &fixed);
        let g = TetLoadGeometry3::new(&NODES, &TETS, p.budget, 4, cx).unwrap();
        let frames = FrameTree::new();
        let h = TetLoadHistory3 {
            nodal_forces_n: &[(3, 0.01)],
            ..history(&frames)
        };
        let rate = 1.0 / 1024.0;
        let motion = source(
            DISPLACEMENT_DIMS,
            WORLD,
            [0.0; 3],
            [rate, 0.0, 0.0],
            OutsideDomainPolicy::Refuse,
        );
        let rows = [SupportMotion3 {
            nodes: &[0],
            components: [true, false, false],
            displacement: &motion,
        }];
        let mut velocity = ZERO;
        velocity[0][0] = rate;
        let load = g
            .assemble(&[], VolumeLoad3::None, &[(3, 0.01)], cx)
            .unwrap();
        let config = TransientConfig3::default();
        let mut scheduled = h
            .begin_transient_with_support_history(
                &g, &p, None, &rows, 0.0, &ZERO, &velocity, &ZERO, config, cx,
            )
            .unwrap();
        let mut direct = load
            .begin_transient_with_supports(&p, None, 0.0, &ZERO, &velocity, &ZERO, config, cx)
            .unwrap();
        for _ in 0..16 {
            let time = scheduled.state().time_s + 1.0 / 128.0;
            let a = scheduled
                .advance_history_with_supports(&h, &rows, 1.0 / 128.0, cx)
                .unwrap();
            let b = direct
                .advance_with_supports(&load, &[(0, rate * time)], 1.0 / 128.0, cx)
                .unwrap();
            assert_eq!(a, b);
            assert_eq!(scheduled.state(), direct.state());
            assert_eq!(
                scheduled.nodal_motion(cx).unwrap(),
                direct.nodal_motion(cx).unwrap()
            );
        }
        assert!(scheduled.state().velocity[0].abs() > 1e-6);
    });
}

#[test]
fn g3_rotating_support_keeps_physical_frame_time_beyond_clamped_value_table() {
    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(6.0, 60.0, 0.0, [75; 32]).unwrap();
        let p = problem(&material, &FIXED);
        let g = TetLoadGeometry3::new(&NODES, &TETS, p.budget, 4, cx).unwrap();
        let mut frames = FrameTree::new();
        frames.add(Frame {
            id: FrameId(1),
            name: "moving-support-vector".into(),
            parent: WORLD,
            motion: FrameMotion::Rotating {
                axis: [0.0, 0.0, 1.0],
                center: Vec3::new(1e6, -2e6, 0.0),
                rate: QtyAny::new(2.0, Dims([0, 0, -1, 0, 0, 0])),
            },
        });
        let h = history(&frames);
        let amplitude = 1.0 / 1024.0;
        let motion = source(
            DISPLACEMENT_DIMS,
            FrameId(1),
            [0.0; 3],
            [amplitude, 0.0, 0.0],
            OutsideDomainPolicy::Clamp,
        );
        let rows = [SupportMotion3 {
            nodes: &ALL,
            components: [true; 3],
            displacement: &motion,
        }];
        // Use the exact admitted endpoint for q0; the history adapter forbids
        // mixing an independently rounded initialization with the sampled BC.
        let start = 1.25;
        let sample = h.at_with_supports(&g, &p, &rows, start, 0.01, cx).unwrap();
        let mut displacement = ZERO;
        for &(d, u) in sample.prescribed_dofs_m() {
            displacement[d / 3][d % 3] = u;
        }
        let angle = 2.0_f64 * start;
        let velocity = [[
            -2.0 * amplitude * angle.sin(),
            2.0 * amplitude * angle.cos(),
            0.0,
        ]; 4];
        let acceleration = displacement.map(|row| row.map(|u| -4.0 * u));
        let mut run = h
            .begin_transient_with_support_history(
                &g,
                &p,
                None,
                &rows,
                start,
                &displacement,
                &velocity,
                &acceleration,
                Default::default(),
                cx,
            )
            .unwrap();
        run.advance_history_with_supports(&h, &rows, 1.0 / 128.0, cx)
            .unwrap();
        let endpoint = run.nodal_motion(cx).unwrap();
        let angle = 2.0 * run.state().time_s;
        for q in &endpoint[0] {
            assert!((q[0] - amplitude * angle.cos()).abs() < 1e-12);
            assert!((q[1] - amplitude * angle.sin()).abs() < 1e-12);
        }
        assert!((endpoint[0][0][0] - amplitude * 2.0_f64.cos()).abs() > 1e-4);
    });
}

#[test]
fn g4_invalid_histories_initial_mismatch_and_cancellation_publish_no_partial_motion() {
    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(6.0, 60.0, 0.0, [75; 32]).unwrap();
        let p = problem(&material, &FIXED);
        let g = TetLoadGeometry3::new(&NODES, &TETS, p.budget, 4, cx).unwrap();
        let frames = FrameTree::new();
        let h = history(&frames);
        let motion = source(
            DISPLACEMENT_DIMS,
            WORLD,
            [0.0; 3],
            [0.001, 0.0, 0.0],
            OutsideDomainPolicy::Refuse,
        );
        let wrong = source(
            TRACTION_DIMS,
            WORLD,
            [0.0; 3],
            [0.001, 0.0, 0.0],
            OutsideDomainPolicy::Refuse,
        );
        let valid = SupportMotion3 {
            nodes: &ALL,
            components: [true; 3],
            displacement: &motion,
        };
        assert!(
            h.begin_transient_with_support_history(
                &g,
                &p,
                None,
                &[valid],
                0.0,
                &[[0.001, 0.0, 0.0]; 4],
                &ZERO,
                &ZERO,
                Default::default(),
                cx,
            )
            .is_err()
        );
        let mut run = h
            .begin_transient_with_support_history(
                &g,
                &p,
                None,
                &[valid],
                0.0,
                &ZERO,
                &ZERO,
                &ZERO,
                Default::default(),
                cx,
            )
            .unwrap();
        let before = run.state().clone();
        let full_before = run.nodal_motion(cx).unwrap();
        for bad in [
            vec![SupportMotion3 {
                displacement: &wrong,
                ..valid
            }],
            vec![valid, valid],
            vec![SupportMotion3 {
                nodes: &[4],
                ..valid
            }],
        ] {
            assert!(
                run.advance_history_with_supports(&h, &bad, 0.01, cx)
                    .is_err()
            );
        }
        assert!(
            run.advance_history_with_supports(&h, &[valid], 2.0, cx)
                .is_err()
        );
        for step in [0.0, -0.01, f64::NAN, f64::INFINITY] {
            assert!(
                run.advance_history_with_supports(&h, &[valid], step, cx)
                    .is_err()
            );
        }
        with_cx(true, |cancelled| {
            assert!(
                run.advance_history_with_supports(&h, &[valid], 0.01, cancelled)
                    .is_err()
            );
        });
        assert_eq!(run.state(), &before);
        assert_eq!(run.nodal_motion(cx).unwrap(), full_before);
        let mut fork = run.clone();
        assert_eq!(
            run.advance_history_with_supports(&h, &[valid], 0.01, cx)
                .unwrap(),
            fork.advance_history_with_supports(&h, &[valid], 0.01, cx)
                .unwrap()
        );
        assert_eq!(
            run.nodal_motion(cx).unwrap(),
            fork.nodal_motion(cx).unwrap()
        );
    });
}
