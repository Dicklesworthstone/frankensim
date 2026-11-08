//! Sparse buckling is checked against physical pencils and independent roots.
use fs_solid::{
    Mesh2, Patch, SolidError,
    linear::{Formulation, LinearProblem, PlaneKind},
    stability::{
        self,
        sparse::{self, SparseBucklingLimits},
    },
};
use fs_sparse::{Coo, Csr};
use std::cell::Cell;

fn map(n: usize) -> Vec<Option<usize>> {
    (0..n).map(Some).collect()
}
fn diagonal(values: &[f64]) -> Csr {
    let mut coo = Coo::new(values.len(), values.len());
    for (i, &value) in values.iter().enumerate() {
        coo.push(i, i, value);
    }
    coo.assemble()
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-7 * expected.abs(),
        "{actual:e} != {expected:e}"
    );
}

// K = Rᵀ R, G = -Rᵀ diag(theta) R, upper bidiagonal R.
// The generalized roots are exactly theta, while K G != G K.
fn pencil(theta: &[f64], unit: f64, reference: f64) -> (Csr, Csr) {
    let n = theta.len();
    let mut k = Coo::new(n, n);
    let mut g = Coo::new(n, n);
    for (i, &value) in theta.iter().enumerate() {
        let row = [
            (i, 1.0),
            ((i + 1).min(n - 1), if i + 1 < n { 0.37 } else { 0.0 }),
        ];
        for &(a, ra) in &row {
            for &(b, rb) in &row {
                if ra * rb != 0.0 {
                    k.push(a, b, unit * ra * rb);
                    g.push(a, b, -unit * reference * value * ra * rb);
                }
            }
        }
    }
    (k.assemble(), g.assemble())
}
fn spectrum() -> Vec<f64> {
    vec![
        0.5, 0.25, 0.125, 0.06, 0.02, -25., -20., -15., -10., -8., -6., -5., -4., -3., -2., -1.,
        -0.5, 0.,
    ]
}

#[test]
fn noncommuting_indefinite_geometric_pencil_selects_largest_algebraic_roots() {
    let (k, g) = pencil(&spectrum(), 1.0, 1.0);
    let n = k.nrows();
    let commutator = (0..n)
        .map(|j| k.get(0, j) * g.get(j, 1) - g.get(0, j) * k.get(j, 1))
        .sum::<f64>();
    assert!(commutator.abs() > 1e-3);
    let result =
        sparse::buckling_loads(&k, &g, &map(n), 3, 250, &SparseBucklingLimits::default()).unwrap();
    let dense = stability::buckling_loads(&k, &g, &map(n), 3, 250).unwrap();
    assert_eq!(result.loads.len(), 3);
    for (i, exact) in [2.0, 4.0, 8.0].into_iter().enumerate() {
        close(result.loads[i], exact);
        close(result.loads[i], dense.loads[i]);
        let mut kp = vec![0.0; n];
        let mut gp = vec![0.0; n];
        k.spmv(&result.modes[i], &mut kp);
        g.spmv(&result.modes[i], &mut gp);
        let r = kp
            .iter()
            .zip(&gp)
            .fold(0.0_f64, |sum, (&a, &b)| sum.hypot(a + result.loads[i] * b));
        let scale = kp.iter().fold(0.0_f64, |sum, &v| sum.hypot(v));
        assert!(r / scale < 2e-8);
        close(
            -result.modes[i]
                .iter()
                .zip(&gp)
                .map(|(a, b)| a * b)
                .sum::<f64>(),
            1.0,
        );
    }
}

#[test]
fn sparse_route_preserves_units_and_reference_load_scaling() {
    for unit in [1e-200, 1.0, 1e200] {
        for reference in [1e-20, 1.0, 1e20] {
            let (k, g) = pencil(
                &[0.5, 0.25, 0.125, 0.08, 0.05, 0.02, 0.01, 0.005, 0.001],
                unit,
                reference,
            );
            let result =
                sparse::buckling_loads(&k, &g, &map(9), 2, 200, &SparseBucklingLimits::default())
                    .unwrap();
            close(result.loads[0] * reference, 2.0);
            close(result.loads[1] * reference, 4.0);
            assert!(result.relative_residuals.iter().all(|r| *r <= 1e-8));
        }
    }
}

#[test]
fn automatic_route_solves_more_than_4096_dofs_with_sparse_fill() {
    let n = 5003;
    let mut diagonal_k: Vec<f64> = (0..n).map(|i| 30.0 + (i % 37) as f64).collect();
    diagonal_k[0] = 2.0;
    diagonal_k[1] = 3.0;
    let k = diagonal(&diagonal_k);
    let g = diagonal(&vec![-1.0; n]);
    let result = stability::buckling_loads(&k, &g, &map(n), 2, 100).unwrap();
    assert_eq!(result.loads.len(), 2);
    close(result.loads[0], 2.0);
    close(result.loads[1], 3.0);
    // The native preparation independently records linear fill and one-entry
    // fronts; a global n² payload is outside this explicit allowance.
    let mut limits = SparseBucklingLimits::default();
    limits.max_factor_bytes = 4 * 1024 * 1024;
    limits.factor.max_front_dim = 2;
    let symbolic = fs_sparse::SymbolicLdlt::analyze_bounded(
        &k,
        fs_sparse::DirectOrdering::Amd,
        &limits.factor,
        &|| false,
    )
    .unwrap();
    let factor = symbolic
        .factor_bounded(
            &k,
            &fs_sparse::LdltOptions {
                allow_2x2: false,
                max_factor_bytes: Some(limits.max_factor_bytes),
                ..fs_sparse::LdltOptions::default()
            },
            &limits.factor,
            &|| false,
        )
        .unwrap();
    assert_eq!(factor.stats().nnz_l, n);
    assert_eq!(factor.stats().max_front_dim, 1);
    assert!(factor.stats().peak_factor_bytes < limits.max_factor_bytes);
    assert!(n * n * 8 > limits.max_factor_bytes);
}

#[test]
fn physical_p1_q1_and_bbar_pencils_agree_with_dense_whitening() {
    for (mesh, formulation, poisson) in [
        (
            Mesh2::triangles(3.0, 0.8, 6, 2),
            Formulation::Standard,
            0.25,
        ),
        (Mesh2::quads(3.0, 0.8, 6, 2), Formulation::Standard, 0.25),
        (Mesh2::quads(3.0, 0.8, 6, 2), Formulation::BBar, 0.499),
    ] {
        let problem = LinearProblem {
            mesh: &mesh,
            youngs: 1000.0,
            poisson,
            plane: PlaneKind::Strain,
            formulation,
            body_force: None,
            dirichlet: vec![(Patch::Left, &|_, _| [0.0, 0.0])],
            traction: vec![(Patch::Right, &|_, _| [-1.0, 0.0])],
            symmetry: vec![],
        };
        let (k, g, map, _) = stability::reduced_pencil(&problem).unwrap();
        let dense = stability::buckling_loads(&k, &g, &map, 2, 400).unwrap();
        let result =
            sparse::buckling_loads(&k, &g, &map, 2, 400, &SparseBucklingLimits::default()).unwrap();
        assert_eq!(result.loads.len(), dense.loads.len());
        for (a, b) in result.loads.iter().zip(&dense.loads) {
            close(*a, *b);
        }
        assert!(result.relative_residuals.iter().all(|r| *r <= 1e-8));
    }
}

#[test]
fn noncompressive_reference_is_empty_but_indefinite_base_never_admits() {
    let (k, g) = pencil(&[-1.0; 12], 1.0, 1.0);
    assert!(
        sparse::buckling_loads(&k, &g, &map(12), 2, 150, &SparseBucklingLimits::default())
            .unwrap()
            .loads
            .is_empty()
    );
    let zero = Coo::new(12, 12).assemble();
    assert!(
        sparse::buckling_loads(&k, &zero, &map(12), 2, 10, &SparseBucklingLimits::default())
            .unwrap()
            .loads
            .is_empty()
    );
    let mut diagonal_k = vec![1.0; 12];
    diagonal_k[11] = -0.01;
    for geometric in [&g, &zero] {
        assert!(matches!(
            sparse::buckling_loads(
                &diagonal(&diagonal_k),
                geometric,
                &map(12),
                2,
                150,
                &SparseBucklingLimits::default()
            ),
            Err(SolidError::SolveFailed { .. })
        ));
    }
}

#[test]
fn each_inner_and_outer_budget_refuses_without_partial_modes() {
    let (k, g) = pencil(&spectrum(), 1.0, 1.0);
    let n = k.nrows();
    let defaults = SparseBucklingLimits::default();
    let mut cases = Vec::new();
    let mut l = defaults;
    l.factor.max_analysis_work = 1;
    cases.push(l);
    let mut l = defaults;
    l.factor.max_numeric_work = 1;
    cases.push(l);
    let mut l = defaults;
    l.factor.max_symbolic_entries = 1;
    cases.push(l);
    let mut l = defaults;
    l.factor.max_front_dim = 1;
    cases.push(l);
    let mut l = defaults;
    l.max_factor_bytes = 1;
    cases.push(l);
    let mut l = defaults;
    l.max_operator_work = 1;
    cases.push(l);
    let mut l = defaults;
    l.max_operator_applications = 1;
    cases.push(l);
    for limits in cases {
        assert!(matches!(
            sparse::buckling_loads(&k, &g, &map(n), 2, 250, &limits),
            Err(SolidError::InvalidInput { .. })
        ));
    }
    assert!(matches!(
        sparse::buckling_loads(&k, &g, &map(n), 2, 1, &defaults),
        Err(SolidError::SolveFailed { .. })
    ));
}

#[test]
fn cancellation_and_invalid_input_preserve_fail_closed_behavior() {
    let (k, g) = pencil(&spectrum(), 1.0, 1.0);
    let n = k.nrows();
    let limits = SparseBucklingLimits::default();
    let polls = Cell::new(0usize);
    sparse::buckling_loads_cancellable(&k, &g, &map(n), 2, 250, &limits, &|| {
        polls.set(polls.get() + 1);
        false
    })
    .unwrap();
    let total = polls.get();
    for stop in [1, 12, 80, total / 2, total - 1] {
        polls.set(0);
        let result = sparse::buckling_loads_cancellable(&k, &g, &map(n), 2, 250, &limits, &|| {
            polls.set(polls.get() + 1);
            polls.get() >= stop
        });
        assert!(
            matches!(result, Err(SolidError::Cancelled)),
            "stop {stop}/{total}: {result:?}"
        );
    }
    let mut wrong_map = map(n);
    wrong_map[n - 1] = Some(0);
    assert!(matches!(
        sparse::buckling_loads(&k, &g, &wrong_map, 2, 100, &limits),
        Err(SolidError::InvalidInput { .. })
    ));
    let mut bad = k.clone();
    bad.row_values_mut(0)[0] = f64::NAN;
    assert!(matches!(
        sparse::buckling_loads(&bad, &g, &map(n), 2, 100, &limits),
        Err(SolidError::InvalidInput { .. })
    ));
    let mut asymmetric = k.clone();
    asymmetric.row_values_mut(0)[1] *= 1.1;
    assert!(matches!(
        sparse::buckling_loads(&asymmetric, &g, &map(n), 2, 100, &limits),
        Err(SolidError::InvalidInput { .. })
    ));
}
