//! Physical material validity must govern stress-constrained density steps.
#![cfg(feature = "cutfem-marquee")]

use fs_ascent::projected_al::{ProjectedAlError, ProjectedAlOptions};
use fs_cutfem::elastic3::adaptive::AdaptiveElasticity3;
use fs_cutfem::elastic3::{CutElasticity3, ElasticityFormulation3, ElasticityOptions3};
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::QuadratureControl3;
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_solver::LinearOp;
use fs_topopt::pipeline::LoadCase;
use fs_topopt::sdf3::design::{StressDesignOptions3, StressDesignStudy3};
use fs_topopt::sdf3::stress::{StressError3, StressOptions3};
use fs_topopt::sdf3::{CutDensityStudy3, Sdf3Elasticity};
use fs_topopt::{EvaluationStop, SimpParams, SolveBudget, SolveControl, SolveProgress};
use std::{cell::Cell, ops::ControlFlow};

struct Full;
impl CutSdf3 for Full {
    fn value(&self, _: [f64; 3]) -> f64 {
        -1.0
    }
    fn enclose(&self, _: [f64; 3], _: [f64; 3]) -> Interval {
        Interval::new(-1.0, -1.0)
    }
    fn derivative_enclose(&self, _: [f64; 3], _: [f64; 3], _: HeightAxis) -> Interval {
        Interval::new(0.0, 0.0)
    }
}
fn poll() -> ControlFlow<()> {
    ControlFlow::Continue(())
}
fn options(formulation: ElasticityFormulation3) -> ElasticityOptions3 {
    ElasticityOptions3 {
        formulation,
        ghost_gamma: 0.0,
        ..Default::default()
    }
}
fn cartesian(limit: f64, formulation: ElasticityFormulation3) -> CutElasticity3 {
    let mut callback = |_| poll();
    let mut control = QuadratureControl3::new(Default::default(), &mut callback).unwrap();
    CutElasticity3::build(
        HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(),
        [1; 3],
        &Full,
        &IsotropicElastic::new(10.0, 0.0, limit).unwrap(),
        &|p| p[0] == 0.0,
        options(formulation),
        &mut control,
    )
    .unwrap()
}
fn adaptive(formulation: ElasticityFormulation3) -> AdaptiveElasticity3 {
    let tree = Octree3::uniform(1, 3, 128).unwrap();
    let mark = *tree
        .leaves()
        .iter()
        .find(|c| c.index() == [1, 0, 0])
        .unwrap();
    let tree = tree.refined(&[mark], poll).unwrap();
    let mut callback = |_| poll();
    let mut control = QuadratureControl3::new(Default::default(), &mut callback).unwrap();
    AdaptiveElasticity3::build(
        HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(),
        &tree,
        &Full,
        &IsotropicElastic::new(10.0, 0.0, 0.1).unwrap(),
        &|p| p[0] == 0.0,
        options(formulation),
        &mut control,
    )
    .unwrap()
}
fn params() -> SimpParams {
    SimpParams {
        e_min: 1e-3,
        penal: 3.0,
        beta: 0.0,
        eta: 0.5,
    }
}
fn field(nodes: &[[f64; 3]], f: impl Fn([f64; 3]) -> [f64; 3]) -> Vec<f64> {
    nodes.iter().flat_map(|p| f(*p)).collect()
}
fn force(op: &impl LinearOp, u: &[f64]) -> Vec<f64> {
    let mut f = vec![0.0; op.n()];
    op.apply(u, &mut f);
    f
}
fn design_options() -> StressDesignOptions3 {
    StressDesignOptions3 {
        stress_limit: 1.0,
        density_floor: 0.2,
        optimizer: ProjectedAlOptions {
            max_evaluations: 128,
            max_backtracks: 20,
            ..Default::default()
        },
        ..Default::default()
    }
}
fn axial_study(limit: f64) -> (CutDensityStudy3, Vec<f64>) {
    let op = cartesian(limit, ElasticityFormulation3::Displacement);
    let u = field(op.nodes(), |p| [0.005 * p[0], 0.0, 0.0]);
    let f = force(&op, &u);
    (CutDensityStudy3::new(op, 0.0, params()), f)
}

fn corner_refusal<O: Sdf3Elasticity>(mut op: O, u: Vec<f64>) {
    // At the corner, ||eps||_F=0.085*sqrt(1.5)>0.1. Every retained
    // integration-point sample is still admissible on these exact grids.
    op.bulk_stress(&u, 100_000, poll).unwrap();
    assert!(op.admit_material_domain(&u, poll).is_err());
    let f = force(&op, &u);
    let small: Vec<_> = f.iter().map(|v| 0.01 * v).collect();
    let incoming = vec![0.63; op.cells()];
    op.set_scales(&incoming).unwrap();
    let mut study = CutDensityStudy3::new(op, 0.0, params());
    let rho = vec![1.0; study.cells()];
    let mut callback = |_| poll();
    let mut control = SolveControl::new(Default::default(), &mut callback);
    for loads in [
        vec![LoadCase {
            force: &f,
            weight: 1.0,
        }],
        vec![
            LoadCase {
                force: &small,
                weight: 1.0,
            },
            LoadCase {
                force: &f,
                weight: 0.0,
            },
        ],
    ] {
        assert!(matches!(
            study.evaluate_stress(&rho, &loads, StressOptions3::default(), &mut control),
            Err(StressError3::MaterialDomain(_))
        ));
        assert_eq!(study.operator().scales(), incoming);
    }
    assert!(control.work().linear_solves > 0);
}

#[test]
fn stress_refuses_corner_domain_exits_on_cartesian_and_hanging_bbar_fields() {
    for formulation in [
        ElasticityFormulation3::Displacement,
        ElasticityFormulation3::BBar,
    ] {
        let op = cartesian(0.1, formulation);
        let u = field(op.nodes(), |p| [0.085 * p[0] * p[1], 0.0, 0.0]);
        corner_refusal(op, u);
        let op = adaptive(formulation);
        assert!(op.physical_nodes().len() > op.nodes().len());
        let u = field(op.nodes(), |p| [0.085 * p[0] * p[1], 0.0, 0.0]);
        corner_refusal(op, u);
    }
}

#[test]
fn constrained_design_backtracks_to_valid_elastic_fields_without_widening_the_limit() {
    let (mut study, f) = axial_study(0.01);
    let loads = [LoadCase {
        force: &f,
        weight: 1.0,
    }];
    let mut callback = |_| poll();
    let mut control = SolveControl::new(Default::default(), &mut callback);
    let mut session =
        StressDesignStudy3::new(&mut study, &loads, &[0.9], design_options(), &mut control)
            .unwrap();
    let initial_work = session.work();
    let report = session.run(3).unwrap();
    assert_eq!(report.work.iterations, 3);
    assert!(report.work.rejected_trials > 0);
    assert!(session.accepted().volume_fraction < 0.9);
    assert!(session.accepted().volume_fraction > 0.79);
    assert!(session.work().linear_iterations > initial_work.linear_iterations);
    assert!(session.work().linear_solves > initial_work.linear_solves + 3);
    assert_eq!(
        session.study().operator().scales(),
        session.accepted().scales
    );
    for u in &session.accepted().displacements {
        session
            .study()
            .operator()
            .admit_material_domain(u, poll)
            .unwrap();
        let strain = u.iter().step_by(3).copied().fold(0.0_f64, f64::max);
        assert!(strain <= 0.01 * (1.0 + 1e-10));
    }
    assert!(session.best_feasible().unwrap().volume_fraction < 0.9);
    assert!(
        session
            .history()
            .windows(2)
            .all(|rows| rows[1].volume_fraction < rows[0].volume_fraction)
    );
}

#[test]
fn initialization_and_restoration_still_refuse_invalid_physical_fields() {
    let (mut study, f) = axial_study(0.01);
    let loads = [LoadCase {
        force: &f,
        weight: 1.0,
    }];
    let mut callback = |_| poll();
    let mut control = SolveControl::new(Default::default(), &mut callback);
    let incoming = study.operator().scales().to_vec();
    assert!(matches!(
        StressDesignStudy3::new(&mut study, &loads, &[0.2], design_options(), &mut control),
        Err(ProjectedAlError::Evaluation(StressError3::MaterialDomain(
            _
        )))
    ));
    assert_eq!(study.operator().scales(), incoming);
    let checkpoint = {
        let session =
            StressDesignStudy3::new(&mut study, &loads, &[0.9], design_options(), &mut control)
                .unwrap();
        session.checkpoint()
    };
    // A saved endpoint is re-solved against the supplied material. Neither
    // prior optimization success nor trial-rejection policy excuses this exit.
    let (mut changed, _) = axial_study(0.002);
    let incoming = changed.operator().scales().to_vec();
    assert!(matches!(
        StressDesignStudy3::restore(
            &mut changed,
            &loads,
            checkpoint,
            design_options(),
            &mut control
        ),
        Err(ProjectedAlError::Evaluation(StressError3::MaterialDomain(
            _
        )))
    ));
    assert_eq!(changed.operator().scales(), incoming);
}

#[test]
fn cancellation_and_spent_linear_budget_are_not_material_trial_rejections() {
    let (mut study, f) = axial_study(0.01);
    let loads = [LoadCase {
        force: &f,
        weight: 1.0,
    }];
    let stop = Cell::new(false);
    let mut callback = |progress: SolveProgress| {
        if stop.get() && progress.stage == "sdf3-stress-material-domain" {
            ControlFlow::Break(())
        } else {
            poll()
        }
    };
    let mut control = SolveControl::new(Default::default(), &mut callback);
    let initialization_iterations;
    {
        let mut session =
            StressDesignStudy3::new(&mut study, &loads, &[0.9], design_options(), &mut control)
                .unwrap();
        initialization_iterations = session.work().linear_iterations;
        let before = session.point().to_vec();
        let scales = session.study().operator().scales().to_vec();
        stop.set(true);
        assert!(matches!(
            session.run(3),
            Err(ProjectedAlError::Evaluation(StressError3::Evaluation(
                EvaluationStop::Cancelled
            )))
        ));
        assert_eq!(session.point(), before);
        assert_eq!(session.study().operator().scales(), scales);
        assert_eq!(session.optimizer_work().rejected_trials, 0);
    }

    let (mut study, f) = axial_study(0.01);
    let loads = [LoadCase {
        force: &f,
        weight: 1.0,
    }];
    let mut callback = |_| poll();
    let mut control = SolveControl::new(
        SolveBudget {
            total_iterations: initialization_iterations + 1,
            ..Default::default()
        },
        &mut callback,
    );
    let mut session =
        StressDesignStudy3::new(&mut study, &loads, &[0.9], design_options(), &mut control)
            .unwrap();
    let before = session.point().to_vec();
    let scales = session.study().operator().scales().to_vec();
    assert!(matches!(
        session.run(3),
        Err(ProjectedAlError::Evaluation(StressError3::Evaluation(
            EvaluationStop::TotalBudget { .. }
        )))
    ));
    assert_eq!(session.point(), before);
    assert_eq!(session.study().operator().scales(), scales);
    assert_eq!(session.optimizer_work().rejected_trials, 0);
}
