//! Held thermal fields through the original minimum-volume stress optimizer.
use std::{cell::Cell, ops::ControlFlow};

use fs_ascent::projected_al::{ProjectedAlError, ProjectedAlOptions};
use fs_blake3::ContentHash;
use fs_cutfem::elastic3::eigenstrain::PreparedEigenstrain3;
use fs_cutfem::elastic3::{CutElasticity3, ElasticMaterialField3};
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_material::elastic::acoustic::ElasticAcousticState;
use fs_solver::LinearOp;
use fs_topopt::pipeline::LoadCase;
use fs_topopt::sdf3::CutDensityStudy3;
use fs_topopt::sdf3::design::{StressDesignOptions3, StressDesignStudy3};
use fs_topopt::sdf3::eigenstrain::{CutEigenstrainCase3, CutEigenstrainPolicy3};
use fs_topopt::sdf3::stress::{StressError3, StressMeasure3};
use fs_topopt::{EvaluationStop, SimpParams, SolveBudget, SolveControl, SolveProgress};

const PEAK: StressMeasure3 = StressMeasure3::SampledPeakBound;

struct Slab;
impl CutSdf3 for Slab {
    fn value(&self, p: [f64; 3]) -> f64 {
        p[2] - 0.73
    }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        Interval::new(lo[2] - 0.73, hi[2] - 0.73)
    }
    fn derivative_enclose(&self, _: [f64; 3], _: [f64; 3], axis: HeightAxis) -> Interval {
        let value = if axis == HeightAxis::Z { 1.0 } else { 0.0 };
        Interval::new(value, value)
    }
}
fn proceed() -> ControlFlow<()> {
    ControlFlow::Continue(())
}
fn policy() -> CutEigenstrainPolicy3 {
    CutEigenstrainPolicy3 {
        maximum_kinematic_strain: 0.1,
        max_points_per_case: 100_000,
    }
}
fn options(limit: f64) -> StressDesignOptions3 {
    StressDesignOptions3 {
        stress_limit: limit,
        density_floor: 0.25,
        optimizer: ProjectedAlOptions {
            tolerance: 2e-6,
            max_evaluations: 1000,
            max_spectral_step: 0.1,
            ..Default::default()
        },
        ..Default::default()
    }
}
fn fixture() -> (
    CutDensityStudy3<CutElasticity3>,
    PreparedEigenstrain3,
    PreparedEigenstrain3,
    Vec<f64>,
) {
    let material = ElasticAcousticState::from_isotropic(
        &IsotropicElastic::new(20.0, 0.25, 0.1).unwrap(),
        3.0,
        ContentHash([18; 32]),
    )
    .unwrap();
    let mut callback = |_| proceed();
    let mut geometry = QuadratureControl3::new(
        QuadratureOptions3 {
            depth: 1,
            ..Default::default()
        },
        &mut callback,
    )
    .unwrap();
    let op = CutElasticity3::build_with_material_field_and_embedded_dirichlet(
        HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(),
        [1; 3],
        &Slab,
        ElasticMaterialField3::Uniform(&material),
        &|p| p[0] == 0.0,
        &|_, _| true,
        Default::default(),
        Default::default(),
        Default::default(),
        &mut geometry,
    )
    .unwrap();
    let free = |_: usize, p: [f64; 3]| {
        Ok([
            0.0002 * (1.0 + p[1]),
            -0.0001,
            0.00015,
            0.00005 * p[0],
            0.0,
            0.0,
        ])
    };
    let field = op.prepare_eigenstrain(&free, 100_000, proceed).unwrap();
    let changed = op
        .prepare_eigenstrain(
            &|cell, p| free(cell, p).map(|e| e.map(|v| 2.5 * v)),
            100_000,
            proceed,
        )
        .unwrap();
    let force = op
        .body_load(&|p| [0.001 * p[1], -0.01, 0.002], proceed)
        .unwrap();
    let study = CutDensityStudy3::new(
        op,
        0.0,
        SimpParams {
            penal: 3.0,
            beta: 2.0,
            eta: 0.5,
            e_min: 1e-3,
        },
    );
    (study, field, changed, force)
}
fn cases<'a>(field: &'a PreparedEigenstrain3, force: &'a [f64]) -> [CutEigenstrainCase3<'a>; 1] {
    [CutEigenstrainCase3 {
        force,
        weight: 1.0,
        eigenstrain: field,
    }]
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

#[test]
fn thermal_strength_design_reduces_volume_and_keeps_matching_stress_fields() {
    for measure in [StressMeasure3::NormalizedAverage, PEAK] {
        let (mut study, field, _, force) = fixture();
        let cases = cases(&field, &force);
        let mut callback = |_| proceed();
        let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
        let cap = study
            .evaluate_eigenstrain_stress(
                &[0.55],
                &cases,
                policy(),
                Default::default(),
                measure,
                &mut control,
            )
            .unwrap()
            .aggregate;
        let options = options(cap);
        let accepted = {
            let mut session = StressDesignStudy3::new_with_eigenstrain(
                &mut study,
                &cases,
                policy(),
                &[0.85],
                options,
                measure,
                &mut control,
            )
            .unwrap();
            assert!(session.feasible());
            let work = session.work();
            let optimizer_work = session.optimizer_work();
            session.run(0).unwrap();
            assert_eq!(session.work(), work);
            assert_eq!(session.optimizer_work(), optimizer_work);
            session.run(8).unwrap();
            assert!(session.history().len() > 1);
            let best = session.best_feasible().unwrap();
            assert!(best.volume_fraction < session.history()[0].volume_fraction);
            assert!(best.aggregate <= cap * (1.0 + options.optimizer.tolerance));
            if measure == PEAK {
                assert!(best.sampled_physical_max <= cap * (1.0 + options.optimizer.tolerance));
            }
            assert_eq!(session.checkpoint().eigenstrain_policy, Some(policy()));
            assert_eq!(
                session.study().operator().scales(),
                session.accepted().scales
            );
            session.accepted().clone()
        };
        let replay = study
            .evaluate_eigenstrain_stress(
                &accepted.rho,
                &cases,
                policy(),
                options.stress,
                measure,
                &mut control,
            )
            .unwrap();
        assert_eq!(replay.displacements, accepted.displacements);
        assert_eq!(replay.adjoints, accepted.adjoints);
        assert_eq!(replay.gradient, accepted.gradient);
        assert_eq!(replay.aggregate, accepted.aggregate);
        // This evidence remains external-force work, not thermal load work.
        assert_eq!(
            accepted.case_compliances[0],
            dot(&force, &accepted.displacements[0])
        );
        assert_eq!(study.operator().scales(), accepted.scales);
    }
}

#[test]
fn a_pure_thermal_strength_constraint_can_have_a_positive_density_gradient() {
    let (mut study, field, _, _) = fixture();
    let force = vec![0.0; study.operator().n()];
    let cases = cases(&field, &force);
    let mut callback = |_| proceed();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let base = study
        .evaluate_eigenstrain_stress(
            &[0.65],
            &cases,
            policy(),
            Default::default(),
            PEAK,
            &mut control,
        )
        .unwrap();
    assert!(base.aggregate > 0.0);
    assert!(base.gradient[0] > 0.0);
    assert_eq!(base.case_compliances, vec![0.0]);
    let mut session = StressDesignStudy3::new_with_eigenstrain(
        &mut study,
        &cases,
        policy(),
        &[0.65],
        options(1.1 * base.aggregate),
        PEAK,
        &mut control,
    )
    .unwrap();
    session.run(3).unwrap();
    assert!(session.accepted().volume_fraction < base.volume_fraction);
    assert!(session.accepted().aggregate < base.aggregate);
    assert!(session.feasible());
    assert!(session.best_feasible().is_some());
    assert_eq!(session.accepted().case_compliances, vec![0.0]);
}

#[test]
fn held_field_checkpoint_rebuilds_the_same_accepted_and_feasible_designs() {
    let (mut study, field, _, force) = fixture();
    let cases = cases(&field, &force);
    let mut callback = |_| proceed();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let cap = study
        .evaluate_eigenstrain_stress(
            &[0.55],
            &cases,
            policy(),
            Default::default(),
            PEAK,
            &mut control,
        )
        .unwrap()
        .aggregate;
    let options = options(cap);
    let (checkpoint, expected, best, history, work) = {
        let mut session = StressDesignStudy3::new_with_eigenstrain(
            &mut study,
            &cases,
            policy(),
            &[0.85],
            options,
            PEAK,
            &mut control,
        )
        .unwrap();
        session.run(2).unwrap();
        let checkpoint = session.checkpoint();
        session.run(2).unwrap();
        (
            checkpoint,
            session.accepted().clone(),
            session.best_feasible().unwrap().clone(),
            session.history().to_vec(),
            session.optimizer_work(),
        )
    };
    let cost = checkpoint.restoration_cost();
    let mut restored = StressDesignStudy3::restore_with_eigenstrain(
        &mut study,
        &cases,
        policy(),
        checkpoint,
        options,
        PEAK,
        &mut control,
    )
    .unwrap();
    restored.run(2).unwrap();
    assert_eq!(restored.accepted().rho, expected.rho);
    assert_eq!(restored.accepted().displacements, expected.displacements);
    assert_eq!(restored.accepted().gradient, expected.gradient);
    assert_eq!(restored.accepted().aggregate, expected.aggregate);
    assert_eq!(restored.best_feasible().unwrap().rho, best.rho);
    assert_eq!(restored.history(), history);
    assert_eq!(restored.optimizer_work().iterations, work.iterations);
    assert_eq!(
        restored.optimizer_work().evaluations,
        work.evaluations + cost
    );
    assert_eq!(restored.study().operator().scales(), expected.scales);
}

#[test]
fn cold_thermal_policy_and_field_changes_cannot_silently_resume_the_wrong_problem() {
    let (mut study, field, changed, force) = fixture();
    let (_, foreign, _, _) = fixture();
    let cases = cases(&field, &force);
    let changed_cases = [CutEigenstrainCase3 {
        eigenstrain: &changed,
        ..cases[0]
    }];
    let foreign_cases = [CutEigenstrainCase3 {
        eigenstrain: &foreign,
        ..cases[0]
    }];
    let loads = [LoadCase {
        force: &force,
        weight: 1.0,
    }];
    let polls = Cell::new(0usize);
    let mut callback = |_| {
        polls.set(polls.get() + 1);
        proceed()
    };
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let options = options(1.0);
    let cold = StressDesignStudy3::new_with_measure(
        &mut study,
        &loads,
        &[0.75],
        options,
        PEAK,
        &mut control,
    )
    .unwrap()
    .checkpoint();
    let hot = StressDesignStudy3::new_with_eigenstrain(
        &mut study,
        &cases,
        policy(),
        &[0.75],
        options,
        PEAK,
        &mut control,
    )
    .unwrap()
    .checkpoint();
    let original = study.operator().scales().to_vec();
    let work = control.work();
    let before_polls = polls.get();
    assert!(matches!(
        StressDesignStudy3::restore_with_measure(
            &mut study,
            &loads,
            hot.clone(),
            options,
            PEAK,
            &mut control
        ),
        Err(ProjectedAlError::Invalid(_))
    ));
    assert!(matches!(
        StressDesignStudy3::restore_with_eigenstrain(
            &mut study,
            &cases,
            policy(),
            cold,
            options,
            PEAK,
            &mut control
        ),
        Err(ProjectedAlError::Invalid(_))
    ));
    let changed_policy = CutEigenstrainPolicy3 {
        maximum_kinematic_strain: 0.05,
        ..policy()
    };
    assert!(matches!(
        StressDesignStudy3::restore_with_eigenstrain(
            &mut study,
            &cases,
            changed_policy,
            hot.clone(),
            options,
            PEAK,
            &mut control
        ),
        Err(ProjectedAlError::Invalid(_))
    ));
    assert_eq!(
        polls.get(),
        before_polls,
        "kind/policy mismatches must precede admission work"
    );
    assert_eq!(control.work(), work);
    assert_eq!(study.operator().scales(), original);
    assert!(
        StressDesignStudy3::restore_with_eigenstrain(
            &mut study,
            &foreign_cases,
            policy(),
            hot.clone(),
            options,
            PEAK,
            &mut control
        )
        .is_err()
    );
    assert_eq!(
        control.work(),
        work,
        "foreign owners must refuse before any solve"
    );
    assert_eq!(study.operator().scales(), original);
    assert!(matches!(
        StressDesignStudy3::restore_with_eigenstrain(
            &mut study,
            &changed_cases,
            policy(),
            hot,
            options,
            PEAK,
            &mut control
        ),
        Err(ProjectedAlError::Invalid(_))
    ));
    assert!(control.work().linear_solves > work.linear_solves);
    assert_eq!(study.operator().scales(), original);
}

#[test]
fn thermal_point_limits_and_interrupted_trials_preserve_the_accepted_state() {
    let (mut study, field, _, force) = fixture();
    let cases = cases(&field, &force);
    let armed = Cell::new(false);
    let mut callback = |progress: SolveProgress| {
        if armed.get() && progress.stage == "sdf3-eigenstrain-stress" {
            ControlFlow::Break(())
        } else {
            proceed()
        }
    };
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let original = study.operator().scales().to_vec();
    let mut too_small = options(1.0);
    too_small.stress.max_points = field.bulk_points() - 1;
    let work = control.work();
    assert!(matches!(
        StressDesignStudy3::new_with_eigenstrain(
            &mut study,
            &cases,
            policy(),
            &[0.75],
            too_small,
            PEAK,
            &mut control
        ),
        Err(ProjectedAlError::Evaluation(StressError3::PointBudget))
    ));
    assert_eq!(control.work(), work);
    assert_eq!(study.operator().scales(), original);
    let mut session = StressDesignStudy3::new_with_eigenstrain(
        &mut study,
        &cases,
        policy(),
        &[0.75],
        options(1.0),
        PEAK,
        &mut control,
    )
    .unwrap();
    let accepted = session.accepted().clone();
    let best = session.best_feasible().unwrap().clone();
    let history = session.history().to_vec();
    armed.set(true);
    assert!(matches!(
        session.run(1),
        Err(ProjectedAlError::Evaluation(StressError3::Evaluation(
            EvaluationStop::Cancelled
        )))
    ));
    assert_eq!(session.accepted().rho, accepted.rho);
    assert_eq!(session.accepted().displacements, accepted.displacements);
    assert_eq!(session.accepted().gradient, accepted.gradient);
    assert_eq!(session.best_feasible().unwrap().rho, best.rho);
    assert_eq!(session.history(), history);
    assert_eq!(session.study().operator().scales(), accepted.scales);
    armed.set(false);
    session.run(1).unwrap();
    assert!(session.history().len() > history.len());
}
