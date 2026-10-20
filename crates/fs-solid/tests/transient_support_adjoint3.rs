//! G1/G3/G4: full-state/support-control transpose of the native Newmark owner.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_solid::linear3::*;
use fs_solid::loads3::transient::TransientConfig3;
use fs_solid::loads3::transient::adjoint::support::*;
use fs_solid::loads3::transient::adjoint::{TransientAdjointConfig3, TransientAdjointError3};
use fs_solid::loads3::{TetLoadGeometry3, VolumeLoad3};

const NODES: [[f64; 3]; 4] = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const ZERO: [[f64; 3]; 4] = [[0.; 3]; 4];

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
                seed: 83,
                kernel_id: 83,
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
fn geometry(cx: &Cx<'_>) -> TetLoadGeometry3<'static> {
    TetLoadGeometry3::new(&NODES, &TETS, TetAssemblyBudget::default(), 4, cx).unwrap()
}
fn config() -> TransientConfig3 {
    TransientConfig3 {
        maximum_displacement_m: 1.,
        mass_damping_per_s: 0.4,
        stiffness_damping_s: 0.006,
        solve: TetStaticSolveConfig {
            relative_residual_tolerance: 1e-12,
            ..Default::default()
        },
        ..Default::default()
    }
}
fn zero_seed() -> NodalMotionCotangent3 {
    NodalMotionCotangent3 {
        q: ZERO.to_vec(),
        v: ZERO.to_vec(),
        a: ZERO.to_vec(),
    }
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}
fn field_dot(a: &[[f64; 3]], b: &[[f64; 3]]) -> f64 {
    a.iter().zip(b).map(|(x, y)| dot(x, y)).sum()
}
fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() < 2e-10 * (1. + a.abs() + b.abs()),
        "{a:e} != {b:e}"
    );
}

#[test]
fn g1_support_transpose_retains_exact_mass_stiffness_and_both_rayleigh_cross_terms() {
    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(6., 60., 0., [83; 32]).unwrap();
        let fixed: Vec<_> = (0..12).filter(|d| *d != 3).rev().collect();
        let problem = problem(&material, &fixed);
        let geometry = geometry(cx);
        let load = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let run = load
            .begin_transient_with_supports(&problem, None, 0., &ZERO, &ZERO, &ZERO, config(), cx)
            .unwrap();
        let adjoint = run
            .prescribed_adjoint(Default::default(), Default::default(), cx)
            .unwrap();
        assert_eq!(adjoint.free_dofs(), &[3]);
        assert!(adjoint.fixed_dofs().windows(2).all(|w| w[0] < w[1]));
        let mut seed = zero_seed();
        seed.q[1][0] = 0.7;
        seed.v[1][0] = -0.2;
        seed.a[1][0] = 0.03;
        seed.q[0][0] = 0.11;
        seed.v[0][0] = -0.07;
        seed.a[0][0] = 0.02;
        // V=1/6, M_ff=.1, M_fc=.05, K_ff=10, K_fc=-10,
        // C_ff=.1, C_fc=-.04. Both consistent mass and stiffness affect support seeds.
        let initial = adjoint.pullback_initialization(&seed, cx).unwrap();
        let lambda = 0.03 / 0.1;
        close(initial.force[0], lambda);
        close(initial.initial.q[0][0], 0.11 + 10. * lambda);
        close(initial.initial.v[0][0], -0.07 + 0.04 * lambda);
        close(initial.initial.a[0][0], 0.02 - 0.05 * lambda);
        close(initial.initial.q[1][0], 0.7 - 10. * lambda);
        close(initial.initial.v[1][0], -0.2 - 0.1 * lambda);
        assert_eq!(initial.initial.a[1][0], 0.);
        let h = 0.1;
        let b = h * h / 4.;
        let d = h / 2.;
        let step = adjoint.pullback_step(h, &seed, cx).unwrap();
        let lambda = (b * 0.7 - d * 0.2 + 0.03) / (0.1 + d * 0.1 + b * 10.);
        let v = -0.07 + 0.04 * lambda;
        let a = 0.02 - 0.05 * lambda + d * v;
        close(step.force[0], lambda);
        close(step.prescribed_displacement[0], 0.11 + 10. * lambda + a / b);
        close(step.previous.q[0][0], -a / b);
        close(step.previous.v[0][0], v - h * a / b);
        close(step.previous.a[0][0], d * v - a);
        assert!(step.true_relative_residual < 1e-11);
        let constrained_force = geometry
            .assemble(&[], VolumeLoad3::None, &[(0, 1.)], cx)
            .unwrap();
        assert_eq!(
            adjoint
                .contract_load_direction(&step.force, &constrained_force, cx)
                .unwrap(),
            0.
        );
    });
}

#[test]
fn g3_full_nodal_objective_reverses_all_controls_through_multiple_native_steps() {
    with_cx(false, |cx| {
        let fixed = [0, 1, 2, 4, 5, 8];
        let free = [3, 6, 7, 9, 10, 11];
        let steps = [0.025, 0.031, 0.02];
        let geometry = geometry(cx);
        let direction =
            |group: usize, index: usize| 0.01 * ((11 * group + 7 * index + 3) as f64).sin();
        let field_direction = |group| {
            (0..4)
                .map(|node| core::array::from_fn(|c| direction(group, 3 * node + c)))
                .collect::<Vec<[f64; 3]>>()
        };
        let q_direction = field_direction(0);
        let v_direction = field_direction(1);
        let mut a_direction = field_direction(2);
        for &d in &free {
            a_direction[d / 3][d % 3] = 0.;
        }
        let force_directions: Vec<Vec<_>> = (0..4)
            .map(|step| (0..free.len()).map(|i| direction(3 + step, i)).collect())
            .collect();
        let support_directions: Vec<Vec<_>> = (0..3)
            .map(|step| (0..fixed.len()).map(|i| direction(7 + step, i)).collect())
            .collect();
        let mut terminal = zero_seed();
        for d in 0..12 {
            terminal.q[d / 3][d % 3] = ((d + 1) as f64).cos();
            terminal.v[d / 3][d % 3] = 0.03 * ((d + 2) as f64).sin();
            terminal.a[d / 3][d % 3] = 0.003 * ((d + 3) as f64).cos();
        }
        for nu in [0., 0.25] {
            let material = TetElasticMaterial::try_new(6., 60., nu, [84; 32]).unwrap();
            let problem = problem(&material, &fixed);
            let run_at = |epsilon: f64, group: usize| {
                let mut q = ZERO;
                let mut v = ZERO;
                let mut a = ZERO;
                for d in 0..12 {
                    q[d / 3][d % 3] = 0.00002 * (d as f64 - 5.)
                        + if group == 0 {
                            epsilon * q_direction[d / 3][d % 3]
                        } else {
                            0.
                        };
                    v[d / 3][d % 3] = 0.00003 * (3. - d as f64)
                        + if group == 1 {
                            epsilon * v_direction[d / 3][d % 3]
                        } else {
                            0.
                        };
                }
                for &d in &fixed {
                    a[d / 3][d % 3] = 0.0002 * (d as f64 - 4.)
                        + if group == 2 {
                            epsilon * a_direction[d / 3][d % 3]
                        } else {
                            0.
                        };
                }
                let load_at = |step| {
                    let force: Vec<_> = free
                        .iter()
                        .enumerate()
                        .map(|(i, &d)| {
                            (
                                d,
                                0.0001 * ((step + i + 1) as f64).sin()
                                    + if group == 3 + step {
                                        epsilon * force_directions[step][i]
                                    } else {
                                        0.
                                    },
                            )
                        })
                        .collect();
                    geometry
                        .assemble(&[], VolumeLoad3::None, &force, cx)
                        .unwrap()
                };
                let mut run = load_at(0)
                    .begin_transient_with_supports(&problem, None, 0., &q, &v, &a, config(), cx)
                    .unwrap();
                for (step, &h) in steps.iter().enumerate() {
                    let support: Vec<_> = fixed
                        .iter()
                        .enumerate()
                        .map(|(i, &d)| {
                            (
                                d,
                                0.00002 * (d as f64 - 5.)
                                    + 0.00003 * ((step + i + 1) as f64).sin()
                                    + if group == 7 + step {
                                        epsilon * support_directions[step][i]
                                    } else {
                                        0.
                                    },
                            )
                        })
                        .collect();
                    run.advance_with_supports(&load_at(step + 1), &support, h, cx)
                        .unwrap();
                }
                run
            };
            let baseline = run_at(0., usize::MAX);
            let adjoint = baseline
                .prescribed_adjoint(Default::default(), Default::default(), cx)
                .unwrap();
            let mut seed = terminal.clone();
            let mut force = vec![vec![]; 4];
            let mut support = vec![vec![]; 3];
            for step in (0..steps.len()).rev() {
                let pb = adjoint.pullback_step(steps[step], &seed, cx).unwrap();
                force[step + 1] = pb.force;
                support[step] = pb.prescribed_displacement;
                seed = pb.previous;
            }
            let initial = adjoint.pullback_initialization(&seed, cx).unwrap();
            force[0] = initial.force;
            assert!(free.iter().all(|&d| initial.initial.a[d / 3][d % 3] == 0.));
            for group in 0..10 {
                let predicted = match group {
                    0 => field_dot(&initial.initial.q, &q_direction),
                    1 => field_dot(&initial.initial.v, &v_direction),
                    2 => field_dot(&initial.initial.a, &a_direction),
                    3..=6 => dot(&force[group - 3], &force_directions[group - 3]),
                    _ => dot(&support[group - 7], &support_directions[group - 7]),
                };
                let evaluate = |epsilon| {
                    let run = run_at(epsilon, group);
                    let motion = run.nodal_motion(cx).unwrap();
                    field_dot(&terminal.q, &motion[0])
                        + field_dot(&terminal.v, &motion[1])
                        + field_dot(&terminal.a, &motion[2])
                };
                let eps = 1e-5;
                let finite_difference = (evaluate(eps) - evaluate(-eps)) / (2. * eps);
                assert!(
                    (predicted - finite_difference).abs() < 2e-7 * (1. + predicted.abs()),
                    "nu={nu}, group={group}: {predicted:e} != {finite_difference:e}"
                );
            }
        }
    });
}

#[test]
fn g1_fully_prescribed_body_has_the_kinematic_transpose_without_a_sparse_solve() {
    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(6., 60., 0.25, [85; 32]).unwrap();
        let fixed: Vec<_> = (0..12).collect();
        let problem = problem(&material, &fixed);
        let geometry = geometry(cx);
        let load = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let run = load
            .begin_transient_with_supports(&problem, None, 0., &ZERO, &ZERO, &ZERO, config(), cx)
            .unwrap();
        let adjoint = run
            .prescribed_adjoint(Default::default(), Default::default(), cx)
            .unwrap();
        let mut seed = zero_seed();
        seed.q[0][0] = 0.7;
        seed.v[0][0] = -0.2;
        seed.a[0][0] = 0.03;
        let h = 0.05;
        let b = h * h / 4.;
        let d = h / 2.;
        let step = adjoint.pullback_step(h, &seed, cx).unwrap();
        assert_eq!(step.iterations, 0);
        assert!(step.force.is_empty());
        let common = (0.03 - d * 0.2) / b;
        close(step.prescribed_displacement[0], 0.7 + common);
        close(step.previous.q[0][0], -common);
        close(step.previous.v[0][0], -0.2 - h * common);
        close(step.previous.a[0][0], -0.03);
        let initial = adjoint.pullback_initialization(&seed, cx).unwrap();
        assert_eq!(initial.initial, seed);
        assert_eq!(initial.iterations, 0);
        assert!(initial.force.is_empty());
    });
}

#[test]
fn g4_invalid_seed_limits_and_cancellation_publish_no_gradient_and_leave_run_unchanged() {
    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(6., 60., 0.25, [86; 32]).unwrap();
        let fixed = [0, 1, 2, 4, 5, 8];
        let problem = problem(&material, &fixed);
        let geometry = geometry(cx);
        let load = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let run = load
            .begin_transient_with_supports(&problem, None, 0., &ZERO, &ZERO, &ZERO, config(), cx)
            .unwrap();
        let before = run.nodal_motion(cx).unwrap();
        assert!(run.adjoint(Default::default(), cx).is_err());
        for budget in [
            SupportAdjointBudget3 {
                maximum_full_dofs: 11,
                maximum_elements: 1,
            },
            SupportAdjointBudget3 {
                maximum_full_dofs: 12,
                maximum_elements: 0,
            },
        ] {
            assert!(matches!(
                run.prescribed_adjoint(Default::default(), budget, cx),
                Err(TransientAdjointError3::Budget)
            ));
        }
        assert!(matches!(
            run.prescribed_adjoint(
                TransientAdjointConfig3 {
                    maximum_transpose_entries: 0,
                    ..Default::default()
                },
                Default::default(),
                cx
            ),
            Err(TransientAdjointError3::Budget)
        ));
        let adjoint = run
            .prescribed_adjoint(Default::default(), Default::default(), cx)
            .unwrap();
        let mut seed = zero_seed();
        seed.q[1][0] = 1.;
        seed.v[0][0] = 0.1;
        let good = adjoint.pullback_step(0.03, &seed, cx).unwrap();
        let mut short = seed.clone();
        short.q.pop();
        assert!(adjoint.pullback_step(0.03, &short, cx).is_err());
        let mut nonfinite = seed.clone();
        nonfinite.a[0][0] = f64::NAN;
        assert!(adjoint.pullback_initialization(&nonfinite, cx).is_err());
        for h in [0., -1., f64::NAN, f64::MIN_POSITIVE, f64::MAX] {
            assert!(adjoint.pullback_step(h, &seed, cx).is_err());
        }
        with_cx(true, |cancelled| {
            assert!(matches!(
                run.prescribed_adjoint(Default::default(), Default::default(), cancelled),
                Err(TransientAdjointError3::Cancelled)
            ));
            assert!(matches!(
                adjoint.pullback_step(0.03, &seed, cancelled),
                Err(TransientAdjointError3::Cancelled)
            ));
            assert!(matches!(
                adjoint.pullback_initialization(&seed, cancelled),
                Err(TransientAdjointError3::Cancelled)
            ));
        });
        assert_eq!(
            adjoint.clone().pullback_step(0.03, &seed, cx).unwrap(),
            good
        );
        assert_eq!(run.nodal_motion(cx).unwrap(), before);
        let homogeneous = load
            .begin_transient(&problem, None, 0., &ZERO, &ZERO, config(), cx)
            .unwrap();
        assert!(
            homogeneous
                .prescribed_adjoint(Default::default(), Default::default(), cx)
                .is_err()
        );
        assert!(homogeneous.adjoint(Default::default(), cx).is_ok());
    });
}
