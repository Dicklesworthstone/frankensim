//! Driven-support trajectory objectives, authored controls and atomic recording.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_qty::{Dims, QtyAny};
use fs_scenario::payload::{
    OrientationParity, OutsideDomainPolicy, PayloadId, PayloadMeta, QuantityContract,
    ReferenceSemantics, SampleSource, TableInterpolation, VectorPayload,
};
use fs_scenario::{FrameTree, WORLD};
use fs_solid::linear3::*;
use fs_solid::loads3::history::supports::{DISPLACEMENT_DIMS, SupportMotion3};
use fs_solid::loads3::history::{TetLoadHistory3, VolumeHistory3};
use fs_solid::loads3::transient::TransientConfig3;
use fs_solid::loads3::transient::adjoint::support::trajectory::{
    NodalMotionObservation3, PrescribedRecordingBudget3,
};
use fs_solid::loads3::transient::adjoint::support::{NodalMotionCotangent3, SupportAdjointBudget3};
use fs_solid::loads3::transient::adjoint::{TransientAdjointConfig3, TransientAdjointError3};
use fs_solid::loads3::{TetLoadGeometry3, VolumeLoad3};

const NODES: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const ZERO: [[f64; 3]; 4] = [[0.0; 3]; 4];
const TIMES: [f64; 3] = [0.025, 0.04, 0.02];
type FullMotion = [Vec<[f64; 3]>; 3];

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
                seed: 93,
                kernel_id: 93,
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
        budget: TetAssemblyBudget::default(),
    }
}
fn zero_seed() -> NodalMotionCotangent3 {
    NodalMotionCotangent3 {
        q: ZERO.to_vec(),
        v: ZERO.to_vec(),
        a: ZERO.to_vec(),
    }
}
fn config() -> TransientConfig3 {
    TransientConfig3 {
        mass_damping_per_s: 0.4,
        stiffness_damping_s: 0.006,
        ..Default::default()
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * (1.0 + a.abs().max(b.abs())),
        "{a:e} != {b:e}"
    );
}

// state, node, q/v/a channel, physical target, residual scale. These are actual
// least-squares observations, including direct support motion at t0 and later.
const TERMS: [(usize, usize, usize, f64, f64); 8] = [
    (0, 0, 0, 0.0002, 0.001),
    (0, 0, 2, 0.003, 0.3),
    (1, 1, 0, 0.0004, 0.001),
    (1, 1, 1, -0.002, 0.02),
    (2, 0, 1, 0.003, 0.02),
    (2, 1, 2, -0.004, 0.3),
    (3, 0, 2, 0.001, 0.3),
    (3, 1, 0, 0.0001, 0.001),
];

fn objective(states: &[FullMotion]) -> (f64, NodalMotionCotangent3, Vec<NodalMotionObservation3>) {
    let mut value = 0.0;
    let mut terminal = zero_seed();
    let mut observations = Vec::new();
    for (state, node, field, target, scale) in TERMS {
        let residual = (states[state][field][node][0] - target) / scale;
        value += 0.5 * residual * residual;
        let derivative = residual / scale;
        if state == 3 && node == 1 {
            terminal.q[node][0] += derivative;
        } else {
            let mut row = NodalMotionObservation3 {
                state,
                node,
                q: [0.0; 3],
                v: [0.0; 3],
                a: [0.0; 3],
            };
            match field {
                0 => row.q[0] = derivative,
                1 => row.v[0] = derivative,
                _ => row.a[0] = derivative,
            }
            observations.push(row);
        }
    }
    // Deliberately supply nonchronological observations. Equal-state rows must
    // retain their caller-defined summation order after endpoint indexing.
    observations.reverse();
    (value, terminal, observations)
}

#[test]
fn g3_complete_trajectory_gradient_matches_observed_physics_for_every_control_family() {
    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(6.0, 60.0, 0.0, [93; 32]).unwrap();
        let fixed: Vec<_> = (0..12).filter(|d| *d != 3).collect();
        let problem = problem(&material, &fixed);
        let geometry =
            TetLoadGeometry3::new(&NODES, &TETS, TetAssemblyBudget::default(), 4, cx).unwrap();
        // Four physical forces, four support displacements, initial free q/v,
        // initial support v/a. Initial free acceleration is solved, not varied.
        let controls = [
            0.0005, 0.002, -0.001, 0.0015, 0.0001, 0.00012, 0.00020, 0.00022, 0.00015, 0.001,
            -0.0005, 0.002,
        ];
        let evaluate = |x: &[f64; 12]| {
            let mut q = ZERO;
            q[0][0] = x[4];
            q[1][0] = x[8];
            let mut v = ZERO;
            v[0][0] = x[10];
            v[1][0] = x[9];
            let mut a = ZERO;
            a[0][0] = x[11];
            let initial = geometry
                .assemble(&[], VolumeLoad3::None, &[(3, x[0])], cx)
                .unwrap();
            let run = initial
                .begin_transient_with_supports(&problem, None, 0.0, &q, &v, &a, config(), cx)
                .unwrap();
            let mut recorded = run
                .record_prescribed_controls(
                    PrescribedRecordingBudget3::default(),
                    TransientAdjointConfig3::default(),
                    SupportAdjointBudget3::default(),
                    cx,
                )
                .unwrap();
            let mut states = vec![recorded.run().nodal_motion(cx).unwrap()];
            for (i, h) in TIMES.into_iter().enumerate() {
                let load = geometry
                    .assemble(&[], VolumeLoad3::None, &[(3, x[i + 1])], cx)
                    .unwrap();
                recorded
                    .advance_with_supports(&load, &[(0, x[5 + i])], h, cx)
                    .unwrap();
                states.push(recorded.run().nodal_motion(cx).unwrap());
            }
            let (value, terminal, observations) = objective(&states);
            (value, recorded, terminal, observations)
        };
        let (_, recorded, terminal, observations) = evaluate(&controls);
        let gradient = recorded
            .reverse_with_observations(&terminal, &observations, cx)
            .unwrap();
        assert_eq!(recorded.sample_count(), 4);
        assert_eq!(gradient.forces.len(), 4);
        assert_eq!(gradient.displacements.len(), 4);
        assert_eq!(gradient.initial.q[0], [0.0; 3]);
        assert_eq!(gradient.initial.a[1][0], 0.0);
        let support_index = recorded.adjoint().fixed_dofs().binary_search(&0).unwrap();
        let expected: Vec<_> = gradient
            .forces
            .iter()
            .map(|f| f.values[0])
            .chain(
                gradient
                    .displacements
                    .iter()
                    .map(|s| s.values[support_index]),
            )
            .chain([
                gradient.initial.q[1][0],
                gradient.initial.v[1][0],
                gradient.initial.v[0][0],
                gradient.initial.a[0][0],
            ])
            .collect();
        for control in 0..controls.len() {
            let epsilon = if (4..=8).contains(&control) {
                1e-7
            } else {
                1e-6
            };
            let mut plus = controls;
            let mut minus = controls;
            plus[control] += epsilon;
            minus[control] -= epsilon;
            let finite_difference = (evaluate(&plus).0 - evaluate(&minus).0) / (2.0 * epsilon);
            close(expected[control], finite_difference, 3e-6);
        }
        close(
            gradient.thermal_equivalent_force[0],
            gradient.forces.iter().map(|r| r.values[0]).sum(),
            1e-12,
        );
        assert_eq!(
            recorded.accepted_times().collect::<Vec<_>>(),
            gradient.forces.iter().map(|f| f.time_s).collect::<Vec<_>>()
        );
        assert_eq!(
            gradient,
            recorded
                .clone()
                .reverse_with_observations(&terminal, &observations, cx)
                .unwrap()
        );
        assert!(gradient.maximum_adjoint_relative_residual <= 1e-10);
    });
}

fn ramp(amplitude: f64) -> VectorPayload {
    VectorPayload::new(
        PayloadMeta::new(
            QuantityContract::Dimensions(DISPLACEMENT_DIMS),
            PayloadId::new("basis/world-cartesian").unwrap(),
            WORLD,
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
                vec![QtyAny::new(0.0, DISPLACEMENT_DIMS); 3],
                vec![
                    QtyAny::new(amplitude, DISPLACEMENT_DIMS),
                    QtyAny::new(0.0, DISPLACEMENT_DIMS),
                    QtyAny::new(0.0, DISPLACEMENT_DIMS),
                ],
            ],
            TableInterpolation::Linear,
            OutsideDomainPolicy::Reject,
        )
        .unwrap(),
    )
    .unwrap()
}

#[test]
fn g1_typed_history_gradient_chains_to_an_authored_support_ramp_amplitude() {
    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(6.0, 60.0, 0.0, [94; 32]).unwrap();
        let fixed: Vec<_> = (0..12).collect();
        let problem = problem(&material, &fixed);
        let geometry =
            TetLoadGeometry3::new(&NODES, &TETS, TetAssemblyBudget::default(), 4, cx).unwrap();
        let frames = FrameTree::new();
        let history = TetLoadHistory3 {
            surfaces: &[],
            volume: VolumeHistory3::Fixed(VolumeLoad3::None),
            nodal_forces_n: &[],
            frames: &frames,
        };
        let amplitude = 0.001;
        let payload = ramp(amplitude);
        let nodes = [0, 1, 2, 3];
        let supports = [SupportMotion3 {
            nodes: &nodes,
            components: [true, false, false],
            displacement: &payload,
        }];
        let velocity = [[amplitude, 0.0, 0.0]; 4];
        let run = history
            .begin_transient_with_support_history(
                &geometry,
                &problem,
                None,
                &supports,
                0.0,
                &ZERO,
                &velocity,
                &ZERO,
                config(),
                cx,
            )
            .unwrap();
        let mut recorded = run
            .record_prescribed_controls(
                PrescribedRecordingBudget3::default(),
                TransientAdjointConfig3::default(),
                SupportAdjointBudget3::default(),
                cx,
            )
            .unwrap();
        for h in [0.125, 0.25] {
            recorded
                .advance_history_with_supports(&history, &supports, h, cx)
                .unwrap();
        }
        let full = recorded.run().nodal_motion(cx).unwrap();
        let mut terminal = zero_seed();
        for node in 0..4 {
            terminal.q[node][0] = full[0][node][0];
            terminal.v[node][0] = 0.4 * full[1][node][0];
        }
        // J=.5*sum(q_x^2)+.2*sum(v_x^2); s=A*t and the explicitly supplied
        // initial support velocity is A. Both dependencies must be chained.
        let gradient = recorded.reverse(&terminal, cx).unwrap();
        assert!(gradient.forces.iter().all(|r| r.values.is_empty()));
        let mut derivative: f64 = gradient.initial.v.iter().map(|row| row[0]).sum();
        for row in &gradient.displacements {
            for (index, &dof) in recorded.adjoint().fixed_dofs().iter().enumerate() {
                if dof % 3 == 0 {
                    derivative += row.time_s * row.values[index];
                }
            }
        }
        close(
            derivative,
            4.0 * amplitude * (0.375_f64.powi(2) + 0.4),
            1e-11,
        );
        for row in gradient.initial.q {
            assert_eq!(row, [0.0; 3]);
        }
    });
}

#[test]
fn g4_recording_budgets_failures_and_cancellation_preserve_the_accepted_path() {
    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(6.0, 60.0, 0.0, [95; 32]).unwrap();
        let fixed: Vec<_> = (0..12).filter(|d| *d != 3).collect();
        let problem = problem(&material, &fixed);
        let geometry =
            TetLoadGeometry3::new(&NODES, &TETS, TetAssemblyBudget::default(), 4, cx).unwrap();
        let load = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let run = load
            .begin_transient_with_supports(&problem, None, 0.0, &ZERO, &ZERO, &ZERO, config(), cx)
            .unwrap();
        // 12 initial force/support scalars + 36 initial-motion + 1 thermal.
        let mut bounded = run
            .clone()
            .record_prescribed_controls(
                PrescribedRecordingBudget3 {
                    maximum_gradient_values: 49,
                    ..Default::default()
                },
                TransientAdjointConfig3::default(),
                SupportAdjointBudget3::default(),
                cx,
            )
            .unwrap();
        assert!(bounded.reverse(&zero_seed(), cx).is_ok());
        assert!(matches!(
            bounded.advance_with_supports(&load, &[], 0.01, cx),
            Err(TransientAdjointError3::Budget)
        ));
        assert_eq!(bounded.sample_count(), 1);
        assert_eq!(bounded.run().state(), run.state());
        let mut recorded = run
            .record_prescribed_controls(
                PrescribedRecordingBudget3 {
                    maximum_steps: 1,
                    maximum_observations: 1,
                    ..Default::default()
                },
                TransientAdjointConfig3::default(),
                SupportAdjointBudget3::default(),
                cx,
            )
            .unwrap();
        assert!(
            recorded
                .advance_with_supports(&load, &[(0, 0.001), (0, 0.001)], 0.01, cx)
                .is_err()
        );
        assert_eq!(recorded.sample_count(), 1);
        with_cx(true, |cancelled| {
            assert!(matches!(
                recorded.advance_with_supports(&load, &[], 0.01, cancelled),
                Err(TransientAdjointError3::Cancelled)
            ));
            assert!(matches!(
                recorded.reverse(&zero_seed(), cancelled),
                Err(TransientAdjointError3::Cancelled)
            ));
        });
        let mut fork = recorded.clone();
        recorded
            .advance_with_supports(&load, &[(0, 0.00001)], 0.01, cx)
            .unwrap();
        fork.advance_with_supports(&load, &[(0, 0.00001)], 0.01, cx)
            .unwrap();
        assert_eq!(
            recorded.run().nodal_motion(cx).unwrap(),
            fork.run().nodal_motion(cx).unwrap()
        );
        assert!(matches!(
            recorded.advance_with_supports(&load, &[], 0.01, cx),
            Err(TransientAdjointError3::Budget)
        ));
        let row = NodalMotionObservation3 {
            state: 2,
            node: 0,
            q: [0.0; 3],
            v: [0.0; 3],
            a: [0.0; 3],
        };
        assert!(
            recorded
                .reverse_with_observations(&zero_seed(), &[row], cx)
                .is_err()
        );
        assert!(matches!(
            recorded.reverse_with_observations(&zero_seed(), &[row, row], cx),
            Err(TransientAdjointError3::Budget)
        ));
        let mut bad = zero_seed();
        bad.a[0][0] = f64::NAN;
        assert!(recorded.reverse(&bad, cx).is_err());
        assert_eq!(recorded.sample_count(), 2);
        assert_eq!(
            recorded.reverse(&zero_seed(), cx).unwrap(),
            fork.reverse(&zero_seed(), cx).unwrap()
        );
    });
}
