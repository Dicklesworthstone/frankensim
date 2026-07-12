//! Feature-gated `frankensim-gmik` elasticity-apply adjoint acceptance.

use fs_adjoint::transpose::{Tape, VjpRegistry};
use fs_adjoint::verify_gradient;
use fs_cutfem::{
    Circle, CutElasticity, CutElasticityOperator, ELASTICITY_APPLY_VJP_OP, Quadtree,
    register_elasticity_apply_vjp,
};
use fs_material::IsotropicElastic;

fn problem<'a>(
    grid: &'a Quadtree,
    sdf: &'a Circle,
    material: &'a IsotropicElastic,
) -> CutElasticity<'a> {
    CutElasticity {
        grid,
        sdf,
        material,
        nitsche_beta: 20.0,
        ghost_gamma: 0.5,
        quad_depth: 3,
        clamp: None,
        boundary_traction: None,
        traction_free_interface: false,
        solver_tol: 1e-13,
        solver_max_iters: 60_000,
    }
}

fn assert_bit_symmetric(operator: &CutElasticityOperator) {
    let matrix = operator.matrix();
    for row in 0..matrix.nrows() {
        let (columns, values) = matrix.row(row);
        for (&column, &value) in columns.iter().zip(values) {
            assert!(
                value.is_finite(),
                "elasticity matrix contains a non-finite entry at ({row}, {column})"
            );
            assert_eq!(
                value.to_bits(),
                matrix.get(column, row).to_bits(),
                "elasticity matrix lost exact symmetry at ({row}, {column})"
            );
        }
    }
}

#[test]
fn cte_004_registered_vjp_matches_fd_gradient_gate() {
    let grid = Quadtree::uniform(3);
    let sdf = Circle {
        center: [0.49, 0.52],
        radius: 0.31,
    };
    let material = IsotropicElastic::new(1.0, 0.3, 10.0).expect("fixture material");
    let cut = problem(&grid, &sdf, &material);
    let operator = cut
        .assemble(&|_, _| [0.0, 0.0], &|_, _| [0.0, 0.0])
        .expect("adjoint operator");
    assert_bit_symmetric(&operator);
    let n = operator.dof_count();
    let point: Vec<f64> = (0..n)
        .map(|i| 1e-3 * ((i * 17 + 3) % 29) as f64 / 29.0)
        .collect();
    let target: Vec<f64> = (0..n)
        .map(|i| 1e-4 * ((i * 11 + 5) % 23) as f64 / 23.0)
        .collect();
    let output = operator.apply_vec(&point);
    let residual: Vec<f64> = output.iter().zip(&target).map(|(a, b)| a - b).collect();
    let expected = operator.apply_transpose_vec(&residual);

    let mut registry = VjpRegistry::new();
    let op_key = register_elasticity_apply_vjp(&mut registry, &operator);
    assert!(op_key.starts_with(ELASTICITY_APPLY_VJP_OP));

    // Register a second same-sized, materially different operator *after* the
    // first. A fixed global key would overwrite the first VJP and silently
    // corrupt the tape below; content-addressed keys keep both live.
    let other_material = IsotropicElastic::new(1.3, 0.22, 10.0).expect("other material");
    let other_cut = problem(&grid, &sdf, &other_material);
    let other_operator = other_cut
        .assemble(&|_, _| [0.0, 0.0], &|_, _| [0.0, 0.0])
        .expect("other adjoint operator");
    let other_key = register_elasticity_apply_vjp(&mut registry, &other_operator);
    assert_ne!(
        op_key, other_key,
        "different CSR values need different VJP keys"
    );

    let mut tape = Tape::new();
    let leaf = tape.leaf(point.clone());
    let finite_output = output.iter().all(|value| value.is_finite());
    let result = tape.apply(&op_key, &[leaf], output);
    let gradients = tape
        .transpose(&registry, result, &residual)
        .expect("registered VJP");
    let gradient = &gradients[&leaf];
    assert_eq!(gradient.len(), expected.len());
    assert!(
        gradient.iter().zip(&expected).all(|(actual, want)| {
            actual.is_finite() && want.is_finite() && actual.to_bits() == want.to_bits()
        }),
        "registry VJP must be the exact symmetric transpose apply"
    );

    let objective = |candidate: &[f64]| {
        let applied = operator.apply_vec(candidate);
        0.5 * applied
            .iter()
            .zip(&target)
            .map(|(a, b)| {
                let residual = a - b;
                residual * residual
            })
            .sum::<f64>()
    };
    let directions: Vec<Vec<f64>> = (0..4)
        .map(|probe| {
            (0..n)
                .map(|i| (((i + 3 * probe) * 13 + 7) % 31) as f64 / 31.0 - 0.5)
                .collect()
        })
        .collect();
    let check = verify_gradient(&objective, &point, gradient, &directions, 1e-6, 2e-6);
    let coverage = registry.coverage();
    let finite_primal = finite_output
        && residual
            .iter()
            .chain(&expected)
            .chain(gradient)
            .all(|value| value.is_finite());
    let finite_fd = check.max_rel_err.is_finite()
        && check
            .pairs
            .iter()
            .all(|(analytic, fd)| analytic.is_finite() && fd.is_finite());
    let pass = finite_primal
        && finite_fd
        && check.pass
        && coverage.0.len() == 2
        && coverage.0.contains(&op_key.as_str())
        && coverage.0.contains(&other_key.as_str());
    println!(
        "{{\"test\":\"cte-004\",\"verdict\":\"{}\",\
         \"detail\":\"content-addressed elasticity apply VJP vs central FD\",\
         \"max_rel_error\":{:.3e},\"directions\":{},\"registered_ops\":{}}}",
        if pass { "pass" } else { "fail" },
        check.max_rel_err,
        check.pairs.len(),
        coverage.0.len()
    );
    assert!(pass, "cte-004 failed: {check:?}");
}
