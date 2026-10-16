//! `frankensim-v6dn.2` graded vector-elasticity completion evidence.
//!
//! - cte-006 G1: fixed-pattern 2:1-graded Q1 MMS convergence;
//! - cte-007 G0/G3: uniform/graded affine field and compliance equivalence.

use fs_cutfem::{CutElasticity, CutSdf, Quadtree};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use std::collections::BTreeSet;
use std::f64::consts::PI;

#[derive(Debug, Clone, Copy)]
struct AlwaysInside;

impl CutSdf for AlwaysInside {
    fn value(&self, _point: [f64; 2]) -> f64 {
        -1.0
    }

    fn gradient(&self, _point: [f64; 2]) -> [f64; 2] {
        [1.0, 0.0]
    }

    fn enclose(&self, _lo: [f64; 2], _hi: [f64; 2]) -> Interval {
        Interval::point(-1.0)
    }
}

fn material() -> IsotropicElastic {
    IsotropicElastic::new(1.0, 0.3, 10.0).expect("compressible fixture material")
}

fn graded_grid(base_level: u32) -> Quadtree {
    let mut grid = Quadtree::with_room(base_level, base_level + 1);
    grid.refine_where(base_level + 1, &|lo, _hi| lo[0] >= 0.5);
    let levels: BTreeSet<_> = grid.leaves().map(|cell| cell.0).collect();
    assert_eq!(levels, BTreeSet::from([base_level, base_level + 1]));
    grid
}

fn problem<'a>(
    grid: &'a Quadtree,
    material: &'a IsotropicElastic,
    clamp: Option<&'a dyn Fn(f64, f64) -> bool>,
    traction: Option<&'a dyn Fn(f64, f64) -> [f64; 2]>,
) -> CutElasticity<'a> {
    CutElasticity {
        grid,
        sdf: &AlwaysInside,
        material,
        nitsche_beta: 20.0,
        ghost_gamma: 0.5,
        stabilization_scaling: fs_cutfem::CutStabilizationScaling::MuScaled,
        quad_depth: 3,
        clamp,
        boundary_traction: traction,
        traction_free_interface: true,
        solver_tol: 1e-13,
        solver_max_iters: 60_000,
    }
}

fn mms_u(x: f64, y: f64) -> [f64; 2] {
    [
        (PI * x).sin() * (PI * y).sin(),
        (2.0 * PI * x).sin() * (PI * y).sin(),
    ]
}

fn mms_gradient(x: f64, y: f64) -> [[f64; 2]; 2] {
    [
        [
            PI * (PI * x).cos() * (PI * y).sin(),
            PI * (PI * x).sin() * (PI * y).cos(),
        ],
        [
            2.0 * PI * (2.0 * PI * x).cos() * (PI * y).sin(),
            PI * (2.0 * PI * x).sin() * (PI * y).cos(),
        ],
    ]
}

fn mms_body(lambda: f64, mu: f64, x: f64, y: f64) -> [f64; 2] {
    let u = mms_u(x, y);
    let pi2 = PI * PI;
    [
        (lambda + 3.0 * mu) * pi2 * u[0]
            - 2.0 * (lambda + mu) * pi2 * (2.0 * PI * x).cos() * (PI * y).cos(),
        (lambda + 6.0 * mu) * pi2 * u[1] - (lambda + mu) * pi2 * (PI * x).cos() * (PI * y).cos(),
    ]
}

#[test]
#[allow(
    clippy::float_cmp,
    reason = "the manufactured trace is exactly zero at dyadic box coordinates"
)]
fn cte_006_fixed_pattern_graded_mms_order_gate() {
    let material = material();
    let (lambda, mu) = material.lame();
    let body = |x: f64, y: f64| mms_body(lambda, mu, x, y);
    let clamp_box = |x: f64, y: f64| x == 0.0 || x == 1.0 || y == 0.0 || y == 1.0;
    let mut errors = Vec::new();

    for level in [3u32, 4, 5] {
        let grid = graded_grid(level);
        let cut = problem(&grid, &material, Some(&clamp_box), None);
        let solution = cut.solve(&body, &mms_u).expect("graded MMS solve");
        let active_levels: BTreeSet<_> =
            solution.active_cells().iter().map(|cell| cell.0).collect();
        assert_eq!(
            active_levels,
            BTreeSet::from([level, level + 1]),
            "the physical solve must retain the fixed mixed-level pattern"
        );
        let (l2, h1) = cut.l2_h1_error(&solution, &mms_u, &mms_gradient);
        assert!(
            l2.is_finite() && l2 > 0.0 && h1.is_finite() && h1 > 0.0,
            "level {level} errors must be finite and nondegenerate: L2={l2:e}, H1={h1:e}"
        );
        errors.push((l2, h1, solution.dof_count(), solution.iters));
    }

    assert!(
        errors
            .windows(2)
            .all(|pair| pair[1].0 < pair[0].0 && pair[1].1 < pair[0].1),
        "graded MMS errors must decrease strictly: {errors:?}"
    );
    let l2_fit_order = 0.5 * (errors[0].0 / errors[2].0).log2();
    let h1_fit_order = 0.5 * (errors[0].1 / errors[2].1).log2();
    assert!(
        (1.75..=2.25).contains(&l2_fit_order),
        "graded Q1 L2 fitted order {l2_fit_order:.6} is not approximately two"
    );
    assert!(
        (0.75..=1.25).contains(&h1_fit_order),
        "graded Q1 H1 fitted order {h1_fit_order:.6} is not approximately one"
    );
    println!(
        "{{\"suite\":\"fs-cutfem/graded-elasticity\",\"case\":\"cte-006\",\"levels\":[3,4,5],\"l2\":[{:.10e},{:.10e},{:.10e}],\"h1\":[{:.10e},{:.10e},{:.10e}],\"dofs\":[{},{},{}],\"iters\":[{},{},{}],\"l2_fit_order\":{l2_fit_order:.8},\"h1_fit_order\":{h1_fit_order:.8}}}",
        errors[0].0,
        errors[1].0,
        errors[2].0,
        errors[0].1,
        errors[1].1,
        errors[2].1,
        errors[0].2,
        errors[1].2,
        errors[2].2,
        errors[0].3,
        errors[1].3,
        errors[2].3,
    );
}

#[test]
#[allow(
    clippy::float_cmp,
    clippy::too_many_lines,
    reason = "box-face identity is exact on the dyadic grid and the test retains one complete equivalence narrative"
)]
fn cte_007_uniform_and_graded_affine_field_and_compliance_equivalence() {
    let material = material();
    let (lambda, mu) = material.lame();
    let displacement_scale = 0.01;
    let sigma_xx = (lambda + 2.0 * mu) * displacement_scale;
    let sigma_yy = lambda * displacement_scale;
    let exact = |x: f64, _y: f64| [displacement_scale * x, 0.0];
    let exact_gradient = |_: f64, _: f64| [[displacement_scale, 0.0], [0.0, 0.0]];
    let clamp_left = |x: f64, _: f64| x == 0.0;
    let traction = |x: f64, y: f64| {
        if x == 1.0 {
            [sigma_xx, 0.0]
        } else if y == 0.0 {
            [0.0, -sigma_yy]
        } else if y == 1.0 {
            [0.0, sigma_yy]
        } else {
            [0.0, 0.0]
        }
    };
    let zero = |_: f64, _: f64| [0.0, 0.0];
    let uniform_grid = Quadtree::with_room(3, 4);
    let graded_grid = graded_grid(3);
    let uniform_problem = problem(&uniform_grid, &material, Some(&clamp_left), Some(&traction));
    let graded_problem = problem(&graded_grid, &material, Some(&clamp_left), Some(&traction));
    let uniform = uniform_problem
        .solve(&zero, &zero)
        .expect("uniform affine solve");
    let graded = graded_problem
        .solve(&zero, &zero)
        .expect("graded affine solve");
    let replay = graded_problem
        .solve(&zero, &zero)
        .expect("graded affine deterministic replay");

    let (uniform_l2, uniform_h1) = uniform_problem.l2_h1_error(&uniform, &exact, &exact_gradient);
    let (graded_l2, graded_h1) = graded_problem.l2_h1_error(&graded, &exact, &exact_gradient);
    assert!(
        uniform_l2 < 2e-9 && uniform_h1 < 2e-8,
        "uniform affine field error L2={uniform_l2:e}, H1={uniform_h1:e}"
    );
    assert!(
        graded_l2 < 2e-9 && graded_h1 < 2e-8,
        "graded affine field error L2={graded_l2:e}, H1={graded_h1:e}"
    );

    for point in [[0.1875, 0.375], [0.4375, 0.625], [0.6875, 0.375]] {
        let uniform_cell = uniform_grid
            .find_leaf_at(point[0], point[1])
            .expect("uniform probe cell");
        let graded_cell = graded_grid
            .find_leaf_at(point[0], point[1])
            .expect("graded probe cell");
        let (uniform_value, uniform_gradient) = uniform
            .value_gradient(&uniform_grid, uniform_cell, point)
            .expect("uniform physical probe");
        let (graded_value, graded_gradient) = graded
            .value_gradient(&graded_grid, graded_cell, point)
            .expect("graded physical probe");
        for (left, right) in uniform_value
            .iter()
            .chain(uniform_gradient.iter().flatten())
            .zip(graded_value.iter().chain(graded_gradient.iter().flatten()))
        {
            assert!(
                (left - right).abs() <= 2e-8,
                "uniform/graded physical probe mismatch at {point:?}: {left:e} vs {right:e}"
            );
        }
    }

    let expected_compliance = sigma_xx * displacement_scale;
    for (label, compliance) in [
        ("uniform", uniform.compliance()),
        ("graded", graded.compliance()),
    ] {
        assert!(
            (compliance - expected_compliance).abs() <= 2e-9,
            "{label} compliance {compliance:.17e} differs from exact external work {expected_compliance:.17e}"
        );
    }
    assert!(
        (uniform.compliance() - graded.compliance()).abs() <= 2e-9,
        "uniform/graded compliance mismatch: {:.17e} vs {:.17e}",
        uniform.compliance(),
        graded.compliance()
    );

    assert_eq!(graded.active_cells(), replay.active_cells());
    assert_eq!(graded.coefficients().len(), replay.coefficients().len());
    for (left, right) in graded.coefficients().iter().zip(replay.coefficients()) {
        assert_eq!(left.to_bits(), right.to_bits(), "coefficient replay bits");
    }
    assert_eq!(graded.nodal().len(), replay.nodal().len());
    for ((left_node, left), (right_node, right)) in graded.nodal().iter().zip(replay.nodal()) {
        assert_eq!(left_node, right_node);
        assert_eq!(left[0].to_bits(), right[0].to_bits());
        assert_eq!(left[1].to_bits(), right[1].to_bits());
    }
    assert_eq!(graded.compliance().to_bits(), replay.compliance().to_bits());

    println!(
        "{{\"suite\":\"fs-cutfem/graded-elasticity\",\"case\":\"cte-007\",\"uniform_l2\":{uniform_l2:.10e},\"uniform_h1\":{uniform_h1:.10e},\"graded_l2\":{graded_l2:.10e},\"graded_h1\":{graded_h1:.10e},\"uniform_compliance\":{:.10e},\"graded_compliance\":{:.10e},\"expected_compliance\":{expected_compliance:.10e},\"uniform_dofs\":{},\"graded_dofs\":{}}}",
        uniform.compliance(),
        graded.compliance(),
        uniform.dof_count(),
        graded.dof_count(),
    );
}
