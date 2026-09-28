//! G1/G3/G4: real bulk stress, independent loads and minimum-volume design.
use fs_ascent::projected_al::{ProjectedAlError, ProjectedAlOptions, ProjectedAlStop};
use fs_cutfem::elastic3::{CutElasticity3, ElasticityOptions3, adaptive::AdaptiveElasticity3};
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_topopt::pipeline::LoadCase;
use fs_topopt::sdf3::CutDensityStudy3;
use fs_topopt::sdf3::design::{StressDesignOptions3, StressDesignStudy3};
use fs_topopt::sdf3::stress::{StressError3, StressOptions3};
use fs_topopt::{EvaluationStop, SimpParams, SolveBudget, SolveControl, SolveProgress};
use std::{cell::Cell, ops::ControlFlow};

struct Slab;
impl CutSdf3 for Slab {
    fn value(&self, p: [f64; 3]) -> f64 {
        p[2] - 0.73
    }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        Interval::new(lo[2], hi[2]) - Interval::new(0.73, 0.73)
    }
    fn derivative_enclose(&self, _: [f64; 3], _: [f64; 3], a: HeightAxis) -> Interval {
        let d = if a == HeightAxis::Z { 1.0 } else { 0.0 };
        Interval::new(d, d)
    }
}
fn uniform(count: usize, beta: f64) -> (CutDensityStudy3, Vec<f64>) {
    uniform_with_params(
        count,
        SimpParams {
            beta,
            ..Default::default()
        },
    )
}
fn uniform_with_params(count: usize, params: SimpParams) -> (CutDensityStudy3, Vec<f64>) {
    let mut callback = |_| ControlFlow::Continue(());
    let mut quadrature = QuadratureControl3::new(
        QuadratureOptions3 {
            depth: 1,
            ..Default::default()
        },
        &mut callback,
    )
    .unwrap();
    let op = CutElasticity3::build(
        HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(),
        [count; 3],
        &Slab,
        &IsotropicElastic::new(3.0, 0.27, 1.0).unwrap(),
        &|p| p[0] == 0.0,
        ElasticityOptions3::default(),
        &mut quadrature,
    )
    .unwrap();
    let force = op
        .body_load(&|_| [0.0, 0.0, -1.0], || ControlFlow::Continue(()))
        .unwrap();
    (CutDensityStudy3::new(op, 0.15, params), force)
}
fn adaptive(beta: f64) -> (CutDensityStudy3<AdaptiveElasticity3>, Vec<f64>, Vec<f64>) {
    let tree = Octree3::uniform(1, 4, 1000).unwrap();
    let mark = *tree
        .leaves()
        .iter()
        .find(|c| c.index() == [0, 0, 1])
        .unwrap();
    let tree = tree.refined(&[mark], || ControlFlow::Continue(())).unwrap();
    let mut callback = |_| ControlFlow::Continue(());
    let mut quadrature = QuadratureControl3::new(
        QuadratureOptions3 {
            depth: 1,
            ..Default::default()
        },
        &mut callback,
    )
    .unwrap();
    let op = AdaptiveElasticity3::build(
        HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(),
        &tree,
        &Slab,
        &IsotropicElastic::new(3.0, 0.27, 1.0).unwrap(),
        &|p| p[0] == 0.0,
        ElasticityOptions3::default(),
        &mut quadrature,
    )
    .unwrap();
    let y = op
        .body_load(&|p| [0.0, -(0.4 + p[0]), 0.0], || ControlFlow::Continue(()))
        .unwrap();
    let z = op
        .body_load(&|p| [0.1 * p[1], 0.0, -1.0], || ControlFlow::Continue(()))
        .unwrap();
    (
        CutDensityStudy3::new(
            op,
            0.15,
            SimpParams {
                beta,
                ..Default::default()
            },
        ),
        y,
        z,
    )
}

#[test]
fn g1_hanging_node_multiload_stress_gradient_has_the_full_filter_and_adjoint_chain() {
    for (beta, q) in [(2.0, 1.0), (5.0, 1.4)] {
        let (mut study, y, z) = adaptive(beta);
        assert!(study.operator().physical_nodes().len() > study.operator().nodes().len());
        let loads = [
            LoadCase {
                force: &y,
                weight: 2.0,
            },
            LoadCase {
                force: &z,
                weight: 3.0,
            },
        ];
        let rho: Vec<f64> = (0..study.cells())
            .map(|i| 0.37 + 0.019 * i as f64)
            .collect();
        let options = StressOptions3 {
            relaxation_power: q,
            ..Default::default()
        };
        let mut callback = |_| ControlFlow::Continue(());
        let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
        let scales = study.operator().scales().to_vec();
        let exact = study
            .evaluate_stress(&rho, &loads, options, &mut control)
            .unwrap();
        assert_eq!(study.operator().scales(), scales);
        for (actual, expected) in exact.normalized_load_weights.iter().zip([0.4, 0.6]) {
            assert!((actual - expected).abs() < 2e-16);
        }
        assert!(exact.aggregate > 0.0 && exact.aggregate < exact.sampled_relaxed_max);
        assert_eq!(exact.adjoints.len(), 2);
        for index in 0..rho.len() {
            let h = 2e-5;
            let mut plus = rho.clone();
            let mut minus = rho.clone();
            plus[index] += h;
            minus[index] -= h;
            let a = study
                .evaluate_stress(&plus, &loads, options, &mut control)
                .unwrap();
            let b = study
                .evaluate_stress(&minus, &loads, options, &mut control)
                .unwrap();
            let fd = (a.aggregate - b.aggregate) / (2.0 * h);
            let relative =
                (fd - exact.gradient[index]).abs() / fd.abs().max(exact.aggregate * 1e-7);
            assert!(
                relative < 3e-4,
                "beta={beta} q={q} cell={index} fd={fd} adjoint={} error={relative}",
                exact.gradient[index]
            );
        }
    }
}

#[test]
fn g3_opposite_loads_do_not_cancel_and_stress_weights_are_explicitly_normalized() {
    let (mut study, force) = uniform(2, 2.0);
    let opposite: Vec<f64> = force.iter().map(|f| -f).collect();
    let rho = vec![0.6; study.cells()];
    let mut callback = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let one = study
        .evaluate_stress(
            &rho,
            &[LoadCase {
                force: &force,
                weight: 1.0,
            }],
            Default::default(),
            &mut control,
        )
        .unwrap();
    let both = study
        .evaluate_stress(
            &rho,
            &[
                LoadCase {
                    force: &force,
                    weight: 3.0,
                },
                LoadCase {
                    force: &opposite,
                    weight: 7.0,
                },
            ],
            Default::default(),
            &mut control,
        )
        .unwrap();
    assert!(both.aggregate > 0.0);
    assert!((both.aggregate / one.aggregate - 1.0).abs() < 1e-10);
    for (a, b) in both.gradient.iter().zip(&one.gradient) {
        assert!((a - b).abs() < 1e-9 * a.abs().max(1.0));
    }
    assert_eq!(both.point_count, 2 * one.point_count);
    assert_eq!(both.case_relaxed_max[0], both.case_relaxed_max[1]);
}

#[test]
fn g1_void_floor_is_not_given_spurious_strength_by_double_scaling_the_stress() {
    let (mut study, force) = uniform(2, 0.0);
    let loads = [LoadCase {
        force: &force,
        weight: 1.0,
    }];
    let mut callback = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let dense = study
        .evaluate_stress(
            &vec![1.0; study.cells()],
            &loads,
            Default::default(),
            &mut control,
        )
        .unwrap();
    for rho in [0.5, 0.1, 0.001] {
        let result = study
            .evaluate_stress(
                &vec![rho; study.cells()],
                &loads,
                Default::default(),
                &mut control,
            )
            .unwrap();
        let scale = 1e-6 + (1.0 - 1e-6) * rho * rho * rho;
        let expected = dense.aggregate * rho / scale;
        assert!(
            (result.aggregate / expected - 1.0).abs() < 2e-8,
            "rho={rho} actual={} expected={expected}",
            result.aggregate
        );
        assert!(result.aggregate > 3.0 * dense.aggregate);
        assert!((result.sampled_physical_max / dense.sampled_physical_max - 1.0).abs() < 2e-8);
    }
}

fn design_options(cap: f64) -> StressDesignOptions3 {
    StressDesignOptions3 {
        stress_limit: cap,
        density_floor: 0.05,
        optimizer: ProjectedAlOptions {
            tolerance: 2e-6,
            max_evaluations: 2000,
            ..Default::default()
        },
        ..Default::default()
    }
}

#[test]
fn g0_design_refuses_ersatz_strength_foldback_before_any_physical_solve() {
    for (e_min, beta, floor) in [(1e-6, 2.0, 0.001), (0.05, 2.0, 0.05), (1e-6, 8.0, 0.05)] {
        let (mut study, force) = uniform_with_params(
            1,
            SimpParams {
                e_min,
                beta,
                ..Default::default()
            },
        );
        let loads = [LoadCase {
            force: &force,
            weight: 1.0,
        }];
        let mut callback = |_| ControlFlow::Continue(());
        let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
        let original = study.operator().scales().to_vec();
        let options = StressDesignOptions3 {
            density_floor: floor,
            ..Default::default()
        };
        assert!(matches!(
            StressDesignStudy3::new(&mut study, &loads, &[0.8], options, &mut control),
            Err(ProjectedAlError::Invalid(
                "density floor and projection permit ersatz stress foldback"
            ))
        ));
        assert_eq!(control.work().linear_solves, 0);
        assert_eq!(study.operator().scales(), original);
    }
}

#[test]
fn g1_minimum_volume_trajectory_reaches_a_binding_stress_cap() {
    // One real cut Q1 element gives a scalar density oracle: its entire
    // stiffness scales together, so the known active solution is rho=0.55.
    let (mut study, force) = uniform(1, 0.0);
    let loads = [LoadCase {
        force: &force,
        weight: 1.0,
    }];
    let mut callback = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let cap = study
        .evaluate_stress(&[0.55], &loads, Default::default(), &mut control)
        .unwrap()
        .aggregate;
    let options = design_options(cap);
    let mut design =
        StressDesignStudy3::new(&mut study, &loads, &[0.9], options, &mut control).unwrap();
    let result = design.run(400).unwrap();
    eprintln!("active 3D stress cap: {result:?}; rho={:?}", design.point());
    assert_eq!(result.stop, ProjectedAlStop::Converged);
    assert!(design.feasible());
    assert!((design.point()[0] - 0.55).abs() < 1e-5);
    assert!(design.accepted().volume_fraction < 0.7 * design.history()[0].volume_fraction);
    assert!((design.accepted().aggregate / cap - 1.0).abs() < 1e-5);
    assert!(result.multiplier > 0.0);
    assert!(design.best_feasible().unwrap().volume_fraction < 0.56);
    assert!(design.history().len() > 2);
    let accepted = design.accepted().clone();
    assert_eq!(design.study().operator().scales(), accepted.scales);
    drop(design);
    let replay = study
        .evaluate_stress(&accepted.rho, &loads, options.stress, &mut control)
        .unwrap();
    assert_eq!(replay.displacements, accepted.displacements);
    assert_eq!(replay.aggregate, accepted.aggregate);
}

#[test]
fn g4_cancelled_trial_and_postaccept_poll_preserve_matching_fields_and_resume() {
    for stage in ["sdf3-stress-adjoint", "sdf3-stress-design-accepted"] {
        let (mut study, force) = uniform(1, 0.0);
        let loads = [LoadCase {
            force: &force,
            weight: 1.0,
        }];
        let enabled = Cell::new(false);
        let mut callback = |p: SolveProgress| {
            if enabled.get() && p.stage == stage {
                ControlFlow::Break(())
            } else {
                ControlFlow::Continue(())
            }
        };
        let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
        let cap = study
            .evaluate_stress(&[0.55], &loads, Default::default(), &mut control)
            .unwrap()
            .aggregate;
        let mut design = StressDesignStudy3::new(
            &mut study,
            &loads,
            &[0.9],
            design_options(cap),
            &mut control,
        )
        .unwrap();
        enabled.set(true);
        assert!(design.run(1).is_err());
        enabled.set(false);
        if stage == "sdf3-stress-adjoint" {
            assert_eq!(design.point(), [0.9]);
        } else {
            assert_eq!(design.optimizer_work().iterations, 1);
        }
        assert_eq!(design.study().operator().scales(), design.accepted().scales);
        assert!(design.best_feasible().is_some());
        let before = design.work();
        design.run(2).unwrap();
        assert!(design.work().linear_iterations > before.linear_iterations);
        let accepted = design.accepted().clone();
        drop(design);
        let replay = study
            .evaluate_stress(&accepted.rho, &loads, Default::default(), &mut control)
            .unwrap();
        assert_eq!(accepted.displacements, replay.displacements);
        assert_eq!(accepted.gradient, replay.gradient);
    }
}

#[test]
fn g4_cumulative_points_and_linear_budget_refuse_partial_evidence_without_mutation() {
    let (mut study, force) = uniform(2, 2.0);
    let loads = [LoadCase {
        force: &force,
        weight: 1.0,
    }];
    let rho = vec![0.6; study.cells()];
    let original = study.operator().scales().to_vec();
    let mut callback = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let one = study
        .evaluate_stress(&rho, &loads, Default::default(), &mut control)
        .unwrap();
    let limited = StressOptions3 {
        max_points: one.point_count,
        ..Default::default()
    };
    let two = [loads[0], loads[0]];
    assert!(matches!(
        study.evaluate_stress(&rho, &two, limited, &mut control),
        Err(StressError3::PointBudget)
    ));
    assert_eq!(study.operator().scales(), original);
    let mut control = SolveControl::new(
        SolveBudget {
            total_iterations: 1,
            ..Default::default()
        },
        &mut callback,
    );
    assert!(matches!(
        study.evaluate_stress(&rho, &loads, Default::default(), &mut control),
        Err(StressError3::Evaluation(EvaluationStop::TotalBudget { .. }))
    ));
    assert_eq!(control.work().linear_iterations, 1);
    assert_eq!(study.operator().scales(), original);
}
