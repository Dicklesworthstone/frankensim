//! Actual continuum buckling, total gradients and constrained modulus allocation.
use fs_ascent::projected_al::{ProjectedAlError, ProjectedAlOptions};
use fs_solid::SolidError;
use fs_solid::linear::{Formulation, LinearProblem, PlaneKind};
use fs_solid::mesh2::{Mesh2, Patch};
use fs_solid::stability::sensitivity::element_scale_sensitivities;
use fs_solid::stability::{buckling_loads, reduced_pencil_with_element_scales};
use fs_topopt::buckling::{BucklingDesignOptions, BucklingDesignStudy};
use std::cell::Cell;

fn column(mesh: &Mesh2, formulation: Formulation) -> LinearProblem<'_> {
    LinearProblem {
        mesh,
        youngs: 1000.0,
        poisson: 0.28,
        plane: PlaneKind::Stress,
        formulation,
        body_force: None,
        dirichlet: vec![(Patch::Left, &|_, _| [0.0, 0.0])],
        traction: vec![(Patch::Right, &|_, y| [-1.2 - 0.6 * y, 0.08])],
        symmetry: vec![],
    }
}
fn first_factor(problem: &LinearProblem<'_>, scales: &[f64]) -> f64 {
    let (k, g, map, _) = reduced_pencil_with_element_scales(problem, scales).unwrap();
    buckling_loads(&k, &g, &map, 2, 600).unwrap().loads[0]
}
fn options(floor: f64) -> BucklingDesignOptions {
    BucklingDesignOptions {
        minimum_load_factor: floor,
        spectral_steps: 600,
        optimizer: ProjectedAlOptions {
            max_evaluations: 1024,
            tolerance: 1e-6,
            ..Default::default()
        },
    }
}

#[test]
fn scalar_force_loaded_design_reaches_its_independent_exact_optimum() {
    let mesh = Mesh2::quads(2.0, 0.7, 1, 1);
    let problem = column(&mesh, Formulation::BBar);
    let unit_factor = first_factor(&problem, &[1.0]);
    let floor = 0.8 * unit_factor;
    let mut study =
        BucklingDesignStudy::new(&problem, &[1.2], &[0.2], &[2.0], &[0.4], options(floor)).unwrap();
    let initial_cost = study.accepted().sample().objective;
    let report = study.run(160).unwrap();
    // For one force-loaded element with homogeneous supports, stress is
    // independent of the modulus scale: lambda(s)=s*lambda(1). The constrained
    // cost optimum is exactly s=0.8, independent of the optimizer's formulas.
    assert!(
        (study.element_scales()[0] - 0.8).abs() < 3e-4,
        "scale={}, stop={:?}, constraint={}",
        study.element_scales()[0],
        report.stop,
        report.constraint
    );
    assert!(report.work.iterations > 0);
    assert!(report.objective < initial_cost);
    let current = study.accepted();
    let remeasured = first_factor(&problem, study.element_scales());
    assert!((remeasured / (unit_factor * study.element_scales()[0]) - 1.0).abs() < 1e-8);
    assert!((current.physics().buckling().loads[0] / remeasured - 1.0).abs() < 1e-8);
    assert!((report.constraint - (1.0 - remeasured / floor)).abs() < 1e-10);
    assert_eq!(current.sample().objective, report.objective);
    let best = study.best_feasible().unwrap();
    assert!(best.meets_load_floor());
    assert!(best.sample().objective < 0.98 * initial_cost);
    assert!(best.physics().buckling().loads[0] >= floor);
    let fresh_best = first_factor(&problem, best.physics().element_scales());
    assert!(fresh_best >= floor * (1.0 - 1e-10));
    assert!(
        current
            .physics()
            .adjoint_relative_residuals()
            .iter()
            .all(|r| *r <= 1e-8)
    );
    let work = study.work();
    study.run(0).unwrap();
    assert_eq!(study.work(), work);
}

#[test]
fn normalized_constraint_uses_total_element_derivatives_on_each_native_formulation() {
    let mesh = Mesh2::mapped_quads(3, 1, &|s, t| [3.0 * s + 0.1 * s * t, 0.7 * t]);
    let scales = [0.65, 1.1, 1.6];
    for formulation in [Formulation::Standard, Formulation::BBar] {
        let problem = column(&mesh, formulation);
        let floor = 0.9 * first_factor(&problem, &scales);
        let study = BucklingDesignStudy::new(
            &problem,
            &scales,
            &[0.2; 3],
            &[2.0; 3],
            &[0.2, 0.3, 0.5],
            options(floor),
        )
        .unwrap();
        let sample = study.accepted().sample();
        assert_eq!(sample.gradient, vec![0.2, 0.3, 0.5]);
        for element in 0..3 {
            let h = 5e-5 * scales[element];
            let mut plus = scales;
            let mut minus = scales;
            plus[element] += h;
            minus[element] -= h;
            let fd = -(first_factor(&problem, &plus) - first_factor(&problem, &minus))
                / (2.0 * h * floor);
            assert!((sample.constraint_gradient[element] - fd).abs() < 3e-5 * fd.abs().max(1e-3));
        }
        assert_eq!(
            study.accepted().physics().element_scales(),
            study.element_scales()
        );
    }
}

#[test]
fn a_mode_crossing_stops_the_trial_without_freezing_or_installing_its_physics() {
    // Two identical disconnected columns have first loads proportional to their
    // separate element scales. The initial modes are separated; the first AL
    // trial makes both scales 0.9 and must refuse the exact repeated mode.
    let mut mesh = Mesh2::quads(2.0, 0.7, 1, 1);
    let duplicate = mesh.clone();
    let offset = mesh.nodes.len();
    mesh.nodes
        .extend(duplicate.nodes.iter().map(|p| [p[0], p[1] + 2.0]));
    mesh.elems.extend(
        duplicate
            .elems
            .iter()
            .map(|e| e.iter().map(|n| n + offset).collect()),
    );
    for ((_, edges), (_, additional)) in mesh.patches.iter_mut().zip(&duplicate.patches) {
        edges.extend(additional.iter().map(|(a, b)| (a + offset, b + offset)));
    }
    let mut problem = column(&mesh, Formulation::Standard);
    problem.poisson = 0.0;
    problem.traction = vec![(Patch::Right, &|_, _| [-1.0, 0.0])];
    let initial = [1.0, 2.0];
    let floor = 0.5 * first_factor(&problem, &initial);
    let mut study = BucklingDesignStudy::new(
        &problem,
        &initial,
        &[0.2; 2],
        &[3.0; 2],
        &[0.1, 1.1],
        options(floor),
    )
    .unwrap();
    let before = study.accepted().physics().buckling().loads.clone();
    assert!(matches!(study.run(1),
        Err(ProjectedAlError::Evaluation(SolidError::InvalidInput { what }))
        if what.contains("separated eigenvalues")));
    assert_eq!(study.element_scales(), initial);
    assert_eq!(study.accepted().physics().buckling().loads, before);
    assert_eq!(
        study.best_feasible().unwrap().physics().element_scales(),
        initial
    );
    assert_eq!(study.work().iterations, 0);
    assert_eq!(study.work().evaluations, 2);
    assert_eq!(study.work().rejected_trials, 0);
}

#[test]
fn invalid_boxes_weights_and_cancellation_cannot_publish_a_partial_design() {
    let mesh = Mesh2::quads(2.0, 0.7, 1, 1);
    let problem = column(&mesh, Formulation::Standard);
    let floor = 0.8 * first_factor(&problem, &[1.0]);
    for weight in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::MAX] {
        assert!(matches!(
            BucklingDesignStudy::new(&problem, &[1.2], &[0.2], &[2.0], &[weight], options(floor)),
            Err(ProjectedAlError::Invalid(_))
        ));
    }
    for bounds in [
        (0.0, 2.0),
        (1.5, 2.0),
        (0.2, 1.0),
        (2.0, 1.0),
        (0.2, f64::INFINITY),
    ] {
        assert!(matches!(
            BucklingDesignStudy::new(
                &problem,
                &[1.2],
                &[bounds.0],
                &[bounds.1],
                &[0.4],
                options(floor)
            ),
            Err(ProjectedAlError::Invalid(_))
        ));
    }
    for bad_floor in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            BucklingDesignStudy::new(&problem, &[1.2], &[0.2], &[2.0], &[0.4], options(bad_floor)),
            Err(ProjectedAlError::Invalid(_))
        ));
    }
    assert!(matches!(
        BucklingDesignStudy::new(&problem, &[1.2], &[0.2], &[2.0], &[], options(floor)),
        Err(ProjectedAlError::Invalid(_))
    ));
    let enabled = Cell::new(false);
    let polls = Cell::new(0);
    let cancelled = || {
        if !enabled.get() {
            return false;
        }
        polls.set(polls.get() + 1);
        polls.get() >= 8
    };
    let mut study = BucklingDesignStudy::new_cancellable(
        &problem,
        &[1.2],
        &[0.2],
        &[2.0],
        &[0.4],
        options(floor),
        &cancelled,
    )
    .unwrap();
    let before = study.accepted().physics().buckling().loads.clone();
    enabled.set(true);
    assert!(matches!(
        study.run_cancellable(10, &cancelled),
        Err(ProjectedAlError::Cancelled) | Err(ProjectedAlError::Evaluation(SolidError::Cancelled))
    ));
    assert_eq!(study.element_scales(), [1.2]);
    assert_eq!(study.accepted().physics().buckling().loads, before);
    assert_eq!(
        study.best_feasible().unwrap().physics().element_scales(),
        [1.2]
    );
    assert_eq!(study.work().iterations, 0);
    assert_eq!(study.work().evaluations, 2);
    // Source total sensitivities remain directly available to other consumers.
    assert!(element_scale_sensitivities(&problem, &[1.0], 1, 600).is_ok());
}
