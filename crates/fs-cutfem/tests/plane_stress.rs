//! G1/G3 plane-stress elasticity through the original CutFEM operator.
use fs_cutfem::{
    BoundaryTraction, Circle, CutElasticity, CutSdf, CutStabilizationScaling, DesignBoxEdge,
    EdgeBand, ElasticPlane, HalfPlane, Quadtree,
};
use fs_material::IsotropicElastic;
use fs_solver::op::LinearOp;
use std::f64::consts::PI;

fn problem<'a>(
    grid: &'a Quadtree,
    sdf: &'a dyn CutSdf,
    material: &'a IsotropicElastic,
) -> CutElasticity<'a> {
    CutElasticity {
        grid,
        sdf,
        material,
        nitsche_beta: 80.0,
        ghost_gamma: 0.1,
        stabilization_scaling: CutStabilizationScaling::LongitudinalModulus,
        quad_depth: 4,
        clamp: None,
        boundary_traction: None,
        traction_free_interface: false,
        solver_tol: 1e-11,
        solver_max_iters: 20_000,
    }
}

fn zero(_: f64, _: f64) -> [f64; 2] {
    [0.0; 2]
}

#[test]
fn affine_cut_energy_is_the_three_dimensional_traction_free_condensation() {
    // All affine gradients survive hanging-node elimination exactly. The
    // independent 3-D compliance gives sigma_xx=E(ex+nu ey)/(1-nu^2),
    // sigma_yy=E(ey+nu ex)/(1-nu^2), tau_xy=G gamma_xy.
    for graded in [false, true] {
        let mut grid = Quadtree::with_room(2, 3);
        if graded {
            grid.refine_where(3, &|lo, _| lo[0] >= 0.5);
        }
        let sdf = HalfPlane {
            normal: [1.0, 0.0],
            offset: 0.83,
        };
        for nu in [-0.7, 0.3, 0.499_999_9] {
            let material = IsotropicElastic::new(7.3, nu, 0.1).unwrap();
            let mut cut = problem(&grid, &sdf, &material);
            cut.traction_free_interface = true;
            for scaling in [
                CutStabilizationScaling::MuScaled,
                CutStabilizationScaling::LongitudinalModulus,
            ] {
                cut.stabilization_scaling = scaling;
                let operator = cut.assemble_plane_stress(&zero, &zero).unwrap();
                assert_eq!(operator.plane(), ElasticPlane::Stress);
                assert!(!operator.ghost_faces().is_empty());
                let [ex, ey, gamma] = [0.013, -0.007, 0.011];
                let mut u = vec![0.0; operator.dof_count()];
                for (&node, &id) in operator.node_ids() {
                    let [x, y] = grid.node_pos(node);
                    u[2 * id] = ex * x + gamma * y;
                    u[2 * id + 1] = ey * y;
                }
                let mut ku = vec![0.0; u.len()];
                operator.apply(&u, &mut ku);
                let work: f64 = u.iter().zip(&ku).map(|(x, y)| x * y).sum();
                let c = material.youngs / (1.0 - nu * nu);
                let shear = material.youngs / (2.0 * (1.0 + nu));
                let expected =
                    0.83 * (c * (ex * ex + ey * ey + 2.0 * nu * ex * ey) + shear * gamma * gamma);
                assert!(
                    (work / expected - 1.0).abs() < 2e-12,
                    "graded={graded}, nu={nu}: {work} != {expected}"
                );
            }
            if nu > 0.49 {
                assert!(
                    cut.assemble(&zero, &zero).is_err(),
                    "plane stress must not weaken the plane-strain admission gate"
                );
            }
        }
    }
}

fn exact(x: f64, y: f64) -> [f64; 2] {
    let value = (PI * x).sin() * (PI * y).sin();
    [0.01 * value, -0.02 * value]
}

fn exact_gradient(x: f64, y: f64) -> [[f64; 2]; 2] {
    let dx = PI * (PI * x).cos() * (PI * y).sin();
    let dy = PI * (PI * x).sin() * (PI * y).cos();
    [[0.01 * dx, 0.01 * dy], [-0.02 * dx, -0.02 * dy]]
}

#[test]
fn nearly_incompressible_thin_body_has_measured_displacement_orders() {
    let material = IsotropicElastic::new(11.0, 0.499_999_9, 0.1).unwrap();
    let sdf = Circle {
        center: [0.5, 0.5],
        radius: 0.3,
    };
    // Independent divergence of the condensed plane-stress tensor.
    let c11 = material.youngs / (1.0 - material.poisson.powi(2));
    let c12 = material.poisson * c11;
    let shear = material.youngs / (2.0 * (1.0 + material.poisson));
    let force = |x: f64, y: f64| {
        let ss = (PI * x).sin() * (PI * y).sin();
        let cc = (PI * x).cos() * (PI * y).cos();
        [
            PI * PI * ((c11 + shear) * 0.01 * ss + (c12 + shear) * 0.02 * cc),
            PI * PI * (-(c11 + shear) * 0.02 * ss - (c12 + shear) * 0.01 * cc),
        ]
    };
    let mut errors = Vec::new();
    // Coarse curved cuts are pre-asymptotic; use the original elasticity
    // suite's fixed asymptotic ladder and a three-level log-log fit. Moving
    // the cut relative to the lattice can perturb individual adjacent slopes.
    for level in [5, 6, 7] {
        let grid = Quadtree::uniform(level);
        let cut = problem(&grid, &sdf, &material);
        let solution = cut.solve_plane_stress(&force, &exact).unwrap();
        assert_eq!(solution.plane(), ElasticPlane::Stress);
        assert!(solution.rel_residual < cut.solver_tol);
        errors.push(cut.l2_h1_error(&solution, &exact, &exact_gradient));
    }
    assert!(
        errors
            .iter()
            .all(|(l2, h1)| { *l2 > 0.0 && l2.is_finite() && *h1 > 0.0 && h1.is_finite() })
    );
    assert!(
        errors
            .windows(2)
            .all(|pair| pair[1].0 < pair[0].0 && pair[1].1 < pair[0].1),
        "non-monotone errors: {errors:?}"
    );
    let l2_order = 0.5 * (errors[0].0 / errors[2].0).log2();
    let h1_order = 0.5 * (errors[0].1 / errors[2].1).log2();
    assert!(
        (1.8..2.2).contains(&l2_order),
        "L2 order {l2_order}: {errors:?}"
    );
    assert!(
        (0.8..1.2).contains(&h1_order),
        "H1 order {h1_order}: {errors:?}"
    );
}

#[test]
fn checked_traction_and_controlled_solve_share_the_plane_stress_operator() {
    let grid = Quadtree::uniform(3);
    let sdf = HalfPlane {
        normal: [1.0, 0.0],
        offset: 2.0,
    };
    let material = IsotropicElastic::new(13.0, 0.49, 0.1).unwrap();
    let clamp = |x: f64, _: f64| x == 0.0;
    let mut cut = problem(&grid, &sdf, &material);
    cut.clamp = Some(&clamp);
    cut.traction_free_interface = true;
    let support = EdgeBand::new(DesignBoxEdge::Right, 0.2, 0.8).unwrap();
    let load = |_: f64, _: f64| [0.03, -0.01];
    let traction = BoundaryTraction::EdgeBand {
        support,
        value: &load,
    };
    let operator = cut
        .assemble_plane_stress_with_boundary_traction(&zero, &zero, traction)
        .unwrap();
    let solution = cut
        .solve_plane_stress_with_boundary_traction(&zero, &zero, traction)
        .unwrap();
    let std::ops::ControlFlow::Continue(controlled) = operator
        .solve_controlled(cut.solver_tol, cut.solver_max_iters, 7, |_| {
            std::ops::ControlFlow::<()>::Continue(())
        })
        .unwrap()
    else {
        panic!("unexpected cancellation")
    };
    assert_eq!(controlled.coefficients(), solution.coefficients());
    assert_eq!(controlled.compliance(), solution.compliance());
    let total_x: f64 = operator.rhs().iter().step_by(2).sum();
    let total_y: f64 = operator.rhs().iter().skip(1).step_by(2).sum();
    assert!((total_x - 0.018).abs() < 1e-14);
    assert!((total_y + 0.006).abs() < 1e-14);
    cut.boundary_traction = Some(&load);
    assert!(
        cut.assemble_plane_stress_with_boundary_traction(&zero, &zero, traction)
            .is_err()
    );
}

#[test]
fn zero_poisson_reductions_have_identical_operator_and_solution_bits() {
    let grid = Quadtree::uniform(3);
    let sdf = Circle {
        center: [0.5, 0.5],
        radius: 0.3,
    };
    let material = IsotropicElastic::new(1.0, 0.0, 0.1).unwrap();
    let cut = problem(&grid, &sdf, &material);
    let strain = cut.assemble(&zero, &exact).unwrap();
    let stress = cut.assemble_plane_stress(&zero, &exact).unwrap();
    assert_eq!(strain.plane(), ElasticPlane::Strain);
    for row in 0..strain.dof_count() {
        assert_eq!(strain.matrix().row(row), stress.matrix().row(row));
    }
    assert_eq!(strain.rhs(), stress.rhs());
    assert_eq!(
        cut.solve(&zero, &exact).unwrap().coefficients(),
        cut.solve_plane_stress(&zero, &exact)
            .unwrap()
            .coefficients()
    );
}
