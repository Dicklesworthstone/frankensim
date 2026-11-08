//! G1 analytic optimum and G3 crossing/transaction gates on actual continua.
use fs_ascent::projected_al::{ProjectedAlError, ProjectedAlOptions};
use fs_solid::SolidError;
use fs_solid::linear::{Formulation, LinearProblem, PlaneKind};
use fs_solid::mesh2::{Mesh2, Patch};
use fs_solid::stability::{buckling_loads, reduced_pencil_with_element_scales};
use fs_topopt::buckling::{
    AggregateBucklingDesignOptions, AggregateBucklingDesignStudy, BucklingDesignOptions,
};
use std::cell::Cell;

fn columns() -> Mesh2 {
    let mut mesh = Mesh2::quads(2.0, 0.7, 1, 1);
    let second = mesh.clone();
    let offset = mesh.nodes.len();
    mesh.nodes
        .extend(second.nodes.iter().map(|p| [p[0], p[1] + 2.0]));
    mesh.elems.extend(
        second
            .elems
            .iter()
            .map(|e| e.iter().map(|n| n + offset).collect()),
    );
    for ((_, edges), (_, extra)) in mesh.patches.iter_mut().zip(&second.patches) {
        edges.extend(extra.iter().map(|(a, b)| (a + offset, b + offset)));
    }
    mesh
}

fn problem(mesh: &Mesh2, formulation: Formulation) -> LinearProblem<'_> {
    LinearProblem {
        mesh,
        youngs: 1000.0,
        poisson: 0.0,
        plane: PlaneKind::Stress,
        formulation,
        body_force: None,
        dirichlet: vec![(Patch::Left, &|_, _| [0.0, 0.0])],
        traction: vec![(Patch::Right, &|_, _| [-1.0, 0.0])],
        symmetry: vec![],
    }
}

fn fresh(problem: &LinearProblem<'_>, scales: &[f64], sharpness: f64) -> (f64, f64) {
    let (k, g, map, _) = reduced_pencil_with_element_scales(problem, scales).unwrap();
    let modes = buckling_loads(&k, &g, &map, 3, 600).unwrap();
    let first = modes.loads[0];
    let aggregate = first - (1.0 + (-sharpness * (modes.loads[1] - first)).exp()).ln() / sharpness;
    (first, aggregate)
}

fn options(unit: f64, floor: f64, normalized_sharpness: f64) -> AggregateBucklingDesignOptions {
    AggregateBucklingDesignOptions {
        design: BucklingDesignOptions {
            minimum_load_factor: floor * unit,
            spectral_steps: 600,
            optimizer: ProjectedAlOptions {
                max_evaluations: 2048,
                tolerance: 1e-6,
                ..Default::default()
            },
        },
        modes: 2,
        sharpness: normalized_sharpness / unit,
    }
}

#[test]
fn a_repeated_mode_optimum_matches_the_independent_analytic_design() {
    let mesh = columns();
    for formulation in [Formulation::Standard, Formulation::BBar] {
        let p = problem(&mesh, formulation);
        let unit = fresh(&p, &[1.0, 1.0], 1.0).0;
        let policy = options(unit, 0.8, 8.0);
        let mut study = AggregateBucklingDesignStudy::new(
            &p, &[1.3; 2], &[0.2; 2], &[2.0; 2], &[0.5; 2], policy,
        )
        .unwrap();
        let initial_cost = study.accepted().sample().objective;
        let report = study.run(240).unwrap();
        // Each force-loaded column has lambda_e=unit*s_e. Equal positive costs
        // give the exact smooth-constrained optimum s0=s1=.8+ln(2)/8.
        let exact = 0.8 + 2.0_f64.ln() / 8.0;
        for &scale in study.element_scales() {
            assert!(
                (scale - exact).abs() < 4e-4,
                "{formulation:?}: scale={scale}, exact={exact}, report={report:?}"
            );
        }
        assert!(report.objective < 0.8 * initial_cost);
        assert_eq!(
            study.accepted().physics().element_scales(),
            study.element_scales()
        );
        let (first, value) = fresh(&p, study.element_scales(), policy.sharpness);
        assert!((study.accepted().physics().value() - value).abs() < 1e-8 * unit);
        assert!(
            (report.constraint - (1.0 - value / policy.design.minimum_load_factor)).abs() < 1e-9
        );
        assert!(first >= value);
        let best = study.best_feasible().unwrap();
        let (first, value) = fresh(&p, best.physics().element_scales(), policy.sharpness);
        assert!(best.meets_load_floor());
        assert!(best.physics().value() >= policy.design.minimum_load_factor);
        assert!(value >= policy.design.minimum_load_factor * (1.0 - 1e-10));
        assert!(first >= value);
        assert!(best.sample().objective < 0.8 * initial_cost);
        assert!(best.physics().adjoint_relative_residual() <= 1e-8);
        let work = study.work();
        study.run(0).unwrap();
        assert_eq!(study.work(), work);
    }
}

#[test]
fn an_accepted_step_can_cross_exactly_into_a_repeated_mode() {
    let mesh = columns();
    let p = problem(&mesh, Formulation::Standard);
    let unit = fresh(&p, &[1.0, 1.0], 1.0).0;
    let policy = options(unit, 0.5, 80.0);
    let mut study = AggregateBucklingDesignStudy::new(
        &p,
        &[1.0, 2.0],
        &[0.2; 2],
        &[3.0; 2],
        &[0.1, 1.1],
        policy,
    )
    .unwrap();
    let initial_cost = study.accepted().sample().objective;
    study.run(1).unwrap();
    // Native AL's initial cost step is [1,2]-[.1,1.1]=[.9,.9]. The
    // first-mode consumer refuses here; the complete packet has a smooth trace.
    assert!(
        study
            .element_scales()
            .iter()
            .all(|s| (*s - 0.9).abs() < 1e-12)
    );
    let physics = study.accepted().physics();
    assert!((physics.buckling().loads[1] / physics.buckling().loads[0] - 1.0).abs() < 1e-9);
    assert!((physics.gradient()[0] - physics.gradient()[1]).abs() < 1e-8 * unit);
    assert!(study.accepted().sample().objective < initial_cost);
    assert!(study.accepted().meets_load_floor());
    study.run(8).unwrap();
    assert!(study.work().iterations > 1);
    assert_eq!(
        study.accepted().physics().element_scales(),
        study.element_scales()
    );
}

#[test]
fn a_negative_initial_envelope_can_reach_a_positive_conservative_floor() {
    let mesh = columns();
    let p = problem(&mesh, Formulation::Standard);
    let unit = fresh(&p, &[1.0, 1.0], 1.0).0;
    let policy = options(unit, 0.2, 0.5);
    let mut study =
        AggregateBucklingDesignStudy::new(&p, &[1.0; 2], &[0.2; 2], &[2.5; 2], &[0.2; 2], policy)
            .unwrap();
    assert!(study.accepted().physics().value() < 0.0);
    assert!(study.accepted().physics().buckling().loads[0] > policy.design.minimum_load_factor);
    assert!(!study.accepted().meets_load_floor());
    assert!(study.best_feasible().is_none());
    study.run(240).unwrap();
    let exact = 0.2 + 2.0_f64.ln() / 0.5;
    assert!(
        study
            .element_scales()
            .iter()
            .all(|s| (*s - exact).abs() < 8e-4)
    );
    let best = study.best_feasible().unwrap();
    assert!(best.meets_load_floor());
    assert!(
        fresh(&p, best.physics().element_scales(), policy.sharpness).1
            >= policy.design.minimum_load_factor * (1.0 - 1e-10)
    );
}

#[test]
fn invalid_packets_and_late_cancellation_preserve_the_last_accepted_design() {
    let mesh = columns();
    let p = problem(&mesh, Formulation::Standard);
    let unit = fresh(&p, &[1.0, 1.0], 1.0).0;
    let policy = options(unit, 0.8, 8.0);
    for sharpness in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            AggregateBucklingDesignStudy::new(
                &p,
                &[1.3; 2],
                &[0.2; 2],
                &[2.0; 2],
                &[0.5; 2],
                AggregateBucklingDesignOptions {
                    sharpness,
                    ..policy
                },
            ),
            Err(ProjectedAlError::Invalid(_))
        ));
    }
    assert!(matches!(
        AggregateBucklingDesignStudy::new(
            &p,
            &[1.3; 2],
            &[0.2; 2],
            &[2.0; 2],
            &[0.5; 2],
            AggregateBucklingDesignOptions { modes: 0, ..policy },
        ),
        Err(ProjectedAlError::Invalid(_))
    ));
    assert!(matches!(AggregateBucklingDesignStudy::new(
        &p, &[1.3; 2], &[0.2; 2], &[2.0; 2], &[0.5; 2],
        AggregateBucklingDesignOptions { modes: 1, ..policy },
    ), Err(ProjectedAlError::Evaluation(SolidError::InvalidInput { what }))
        if what.contains("packet boundary")));
    let mut study =
        AggregateBucklingDesignStudy::new(&p, &[1.3; 2], &[0.2; 2], &[2.0; 2], &[0.5; 2], policy)
            .unwrap();
    study.run(1).unwrap();
    let scales = study.element_scales().to_vec();
    let value = study.accepted().physics().value();
    let best = study
        .best_feasible()
        .unwrap()
        .physics()
        .element_scales()
        .to_vec();
    let before = study.work();
    let polls = Cell::new(0);
    assert!(matches!(
        study.run_cancellable(20, &|| {
            polls.set(polls.get() + 1);
            polls.get() >= 10
        }),
        Err(ProjectedAlError::Cancelled) | Err(ProjectedAlError::Evaluation(SolidError::Cancelled))
    ));
    assert_eq!(study.element_scales(), scales);
    assert_eq!(study.accepted().physics().value(), value);
    assert_eq!(
        study.best_feasible().unwrap().physics().element_scales(),
        best
    );
    assert_eq!(study.work().iterations, before.iterations);
    assert_eq!(study.work().evaluations, before.evaluations + 1);
    study.run(1).unwrap();
    assert!(study.work().iterations > before.iterations);
}
