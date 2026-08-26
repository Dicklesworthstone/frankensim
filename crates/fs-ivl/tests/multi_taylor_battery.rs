//! Conformance, property, and boundary test battery for multivariate Taylor models
//! (bead frankensim-epic-bedrock-6ys.23.2.1).

#![deny(unsafe_code)]

use fs_ivl::{
    Interval, MAX_MULTIVARIATE_DIM, MAX_MULTIVARIATE_ORDER, TaylorModel, TaylorModel1,
    TaylorModelError, VariableInfo, binomial, generate_multi_indices, term_count,
};

#[test]
fn test_binomial_and_term_counts() {
    assert_eq!(binomial(0, 0), Some(1));
    assert_eq!(binomial(5, 0), Some(1));
    assert_eq!(binomial(5, 5), Some(1));
    assert_eq!(binomial(5, 2), Some(10));
    assert_eq!(binomial(4, 2), Some(6));

    // d=1, p=3 -> binom(1+3, 3) = 4 terms (1, x, x^2, x^3)
    assert_eq!(term_count(1, 3).unwrap(), 4);

    // d=2, p=2 -> binom(2+2, 2) = 6 terms (1, x, y, x^2, xy, y^2)
    assert_eq!(term_count(2, 2).unwrap(), 6);

    // d=3, p=2 -> binom(3+2, 2) = 10 terms
    assert_eq!(term_count(3, 2).unwrap(), 10);

    // Excessive dimension refusal
    let err_dim = term_count(MAX_MULTIVARIATE_DIM + 1, 1).expect_err("dimension over cap");
    assert!(matches!(err_dim, TaylorModelError::OrderTooLarge { .. }));

    // Excessive order refusal
    let err_ord = term_count(1, MAX_MULTIVARIATE_ORDER + 1).expect_err("order over cap");
    assert!(matches!(err_ord, TaylorModelError::OrderTooLarge { .. }));
}

#[test]
fn test_multi_index_graded_properties() {
    for dim in 1..=4 {
        for order in 0..=4 {
            let count = term_count(dim, order).unwrap();
            let indices = generate_multi_indices(dim, order);
            assert_eq!(
                indices.len(),
                count,
                "multi-indices count must match binom(dim+order, order)"
            );

            // Verify graded ordering: degrees must be non-decreasing
            let mut prev_deg = 0;
            for mi in &indices {
                assert_eq!(mi.len(), dim);
                let deg: usize = mi.iter().map(|&d| d as usize).sum();
                assert!(deg >= prev_deg, "degrees must be non-decreasing (grlex)");
                assert!(deg <= order, "degree must not exceed max order");
                prev_deg = deg;
            }

            // Verify uniqueness
            let mut sorted = indices.clone();
            sorted.sort();
            sorted.dedup();
            assert_eq!(
                sorted.len(),
                indices.len(),
                "all multi-indices must be unique"
            );
        }
    }
}

#[test]
fn test_variable_info_and_fixed_axes() {
    let domain = Interval::new(0.0, 2.0);
    let var = VariableInfo::new("x", domain).unwrap().with_unit("meter");
    assert_eq!(var.name, "x");
    assert_eq!(var.unit.as_deref(), Some("meter"));
    assert!((var.center() - 1.0).abs() < 1e-15);
    assert!((var.radius() - 1.0).abs() < 1e-15);
    assert!(!var.is_fixed());

    // Degenerate/fixed zero-width axis
    let fixed_dom = Interval::point(3.1415);
    let fixed_var = VariableInfo::new("param", fixed_dom).unwrap();
    assert!(fixed_var.is_fixed());
    assert_eq!(fixed_var.center(), 3.1415);
    assert_eq!(fixed_var.radius(), 0.0);

    // Non-finite domain refusal
    let err_inf =
        VariableInfo::new("inf", Interval::WHOLE).expect_err("infinite domain must refuse");
    assert_eq!(err_inf, TaylorModelError::NonFiniteDomain);
}

#[test]
fn test_multivariate_variable_construction_and_eval() {
    let x_dom = Interval::new(1.0, 3.0);
    let y_dom = Interval::new(-1.0, 1.0);

    let vars = vec![
        VariableInfo::new("x", x_dom).unwrap(),
        VariableInfo::new("y", y_dom).unwrap(),
    ];

    let tm_x = TaylorModel::variable(0, vars.clone(), 2).expect("variable x");
    let tm_y = TaylorModel::variable(1, vars.clone(), 2).expect("variable y");

    assert_eq!(tm_x.dim(), 2);
    assert_eq!(tm_x.order(), 2);
    assert_eq!(tm_x.term_count(), 6);

    // Evaluate x at box [ [2.0, 2.0], [0.0, 0.0] ] (the center)
    let box_center = vec![Interval::point(2.0), Interval::point(0.0)];
    let val_x = tm_x.eval_box(&box_center).expect("eval center");
    assert!(val_x.contains(2.0));
    assert!((val_x.midpoint() - 2.0).abs() < 1e-14);

    // Evaluate y at box [ [2.0, 2.0], [0.5, 0.5] ]
    let box_eval = vec![Interval::point(2.0), Interval::point(0.5)];
    let val_y = tm_y.eval_box(&box_eval).expect("eval y");
    assert!(val_y.contains(0.5));
    assert!((val_y.midpoint() - 0.5).abs() < 1e-14);

    // Arithmetic: z = x + 2*y
    let tm_2y = tm_y.scale(2.0).expect("scale 2*y");
    let tm_z = tm_x.add(&tm_2y).expect("add x + 2*y");

    // z at (2.0, 0.5) should be 2.0 + 2*0.5 = 3.0
    let val_z = tm_z.eval_box(&box_eval).expect("eval z");
    assert!(val_z.contains(3.0));
    assert!((val_z.midpoint() - 3.0).abs() < 1e-14);
    assert!(val_z.width() < 1e-14);
}

#[test]
fn test_tm1_conversion_round_trip() {
    let dom = Interval::new(-2.0, 2.0);
    let tm1 = TaylorModel1::variable(dom, 3).expect("tm1 variable");

    // Convert to multivariate TaylorModel
    let mtm = TaylorModel::from_tm1(&tm1, "t", Some("s".into())).expect("from tm1");
    assert_eq!(mtm.dim(), 1);
    assert_eq!(mtm.order(), 3);
    assert_eq!(mtm.term_count(), 4);
    assert_eq!(mtm.variables()[0].name, "t");
    assert_eq!(mtm.variables()[0].unit.as_deref(), Some("s"));

    // Evaluate both at an interval
    let test_iv = Interval::new(-1.0, 1.0);
    let tm1_eval = tm1.eval_interval(test_iv);
    let mtm_eval = mtm.eval_box(&[test_iv]).expect("mtm eval");

    assert_eq!(tm1_eval.lo(), mtm_eval.lo());
    assert_eq!(tm1_eval.hi(), mtm_eval.hi());

    // Round-trip back to TaylorModel1
    let tm1_back = mtm.to_tm1().expect("to tm1");
    assert_eq!(tm1_back.order(), tm1.order());
    assert_eq!(tm1_back.domain(), tm1.domain());
    assert_eq!(tm1_back.center(), tm1.center());
}

#[test]
fn test_multivariate_multiplication_and_truncation() {
    let x_dom = Interval::new(0.0, 1.0);
    let y_dom = Interval::new(0.0, 1.0);
    let vars = vec![
        VariableInfo::new("x", x_dom).unwrap(),
        VariableInfo::new("y", y_dom).unwrap(),
    ];

    // Degree 2 model: (x + y) * (x + y) = x^2 + 2xy + y^2 exactly
    let tm_x = TaylorModel::variable(0, vars.clone(), 2).unwrap();
    let tm_y = TaylorModel::variable(1, vars.clone(), 2).unwrap();
    let tm_sum = tm_x.add(&tm_y).unwrap();
    let tm_sq = tm_sum.mul(&tm_sum).unwrap();

    // At x=0.5, y=0.5: sum is 1.0, sq is 1.0
    let pt_eval = vec![Interval::point(0.5), Interval::point(0.5)];
    let val_sq = tm_sq.eval_box(&pt_eval).unwrap();
    assert!(val_sq.contains(1.0));
    assert!((val_sq.midpoint() - 1.0).abs() < 1e-14);

    // Degree 1 model multiplication: tail x^2, 2xy, y^2 is truncated into remainder
    let tm_x1 = TaylorModel::variable(0, vars.clone(), 1).unwrap();
    let tm_y1 = TaylorModel::variable(1, vars.clone(), 1).unwrap();
    let tm_sum1 = tm_x1.add(&tm_y1).unwrap();
    let tm_sq1 = tm_sum1.mul(&tm_sum1).unwrap();

    // Truncated remainder must be non-zero and must enclose true point evaluations
    assert!(tm_sq1.remainder().width() > 0.0);
    let val_sq1 = tm_sq1.eval_box(&pt_eval).unwrap();
    assert!(val_sq1.contains(1.0));
}

#[test]
fn test_multivariate_division_and_refusal() {
    let x_dom = Interval::new(-0.5, 0.5);
    let vars = vec![VariableInfo::new("x", x_dom).unwrap()];

    // Denominator 1 + 0.2*x is strictly positive over [-0.5, 0.5] ([0.9, 1.1])
    let one = TaylorModel::constant(1.0, vars.clone(), 2).unwrap();
    let x = TaylorModel::variable(0, vars.clone(), 2).unwrap();
    let denom = one.add(&x.scale(0.2).unwrap()).unwrap();

    let recip = denom.reciprocal().expect("reciprocal of positive denom");
    assert!(
        recip
            .range()
            .unwrap()
            .encloses(Interval::new(1.0 / 1.1, 1.0 / 0.9))
    );

    // Denominator containing zero must refuse
    let bad_denom = TaylorModel::variable(0, vars.clone(), 2).unwrap(); // x over [-0.5, 0.5] contains 0
    let err_div = one
        .div(&bad_denom)
        .expect_err("div by zero-straddling denom must refuse");
    assert_eq!(err_div, TaylorModelError::DenominatorContainsZero);
}

#[test]
fn test_truncation_order_reduction() {
    let x_dom = Interval::new(0.0, 2.0);
    let vars = vec![VariableInfo::new("x", x_dom).unwrap()];

    let x = TaylorModel::variable(0, vars.clone(), 3).unwrap();
    let x_cubed = x.mul(&x).unwrap().mul(&x).unwrap();
    assert_eq!(x_cubed.order(), 3);

    // Truncate to order 1
    let x_trunc = x_cubed.truncate(1).unwrap();
    assert_eq!(x_trunc.order(), 1);
    assert_eq!(x_trunc.term_count(), 2);

    // Truncation to higher order refuses
    let err_hi = x_trunc
        .truncate(5)
        .expect_err("higher order truncate must refuse");
    assert!(matches!(
        err_hi,
        TaylorModelError::TruncationOrderTooLarge { .. }
    ));

    // Range containment: truncated model must still enclose x^3 over [0, 2] ([0, 8])
    let full_range = x_trunc.range().unwrap();
    assert!(full_range.contains(0.0));
    assert!(full_range.contains(8.0));
}

#[test]
fn test_recentering_exact_invariance() {
    let x_dom = Interval::new(0.0, 2.0);
    let y_dom = Interval::new(0.0, 2.0);
    let vars = vec![
        VariableInfo::new("x", x_dom).unwrap(),
        VariableInfo::new("y", y_dom).unwrap(),
    ];

    // f(x, y) = x^2 + 2*y
    let x = TaylorModel::variable(0, vars.clone(), 2).unwrap();
    let y = TaylorModel::variable(1, vars.clone(), 2).unwrap();
    let x_sq = x.mul(&x).unwrap();
    let y_2 = y.scale(2.0).unwrap();
    let f = x_sq.add(&y_2).unwrap();

    // Recenter to sub-box [0.0, 1.0] x [0.0, 1.0]
    let sub_box = vec![Interval::new(0.0, 1.0), Interval::new(0.0, 1.0)];
    let f_recenter = f.recenter(&sub_box).expect("recenter sub-box");

    assert_eq!(f_recenter.centers(), &[0.5, 0.5]);

    // Test point (0.3, 0.7) in both models
    let pt = vec![Interval::point(0.3), Interval::point(0.7)];
    let val_orig = f.eval_box(&pt).unwrap();
    let val_recenter = f_recenter.eval_box(&pt).unwrap();

    let expected = 0.3 * 0.3 + 2.0 * 0.7; // 0.09 + 1.40 = 1.49
    assert!(val_orig.contains(expected));
    assert!(val_recenter.contains(expected));
    assert!((val_orig.midpoint() - val_recenter.midpoint()).abs() < 1e-13);
}

#[test]
fn test_subdivision_inclusion_preservation() {
    let x_dom = Interval::new(-1.0, 1.0);
    let vars = vec![VariableInfo::new("x", x_dom).unwrap()];

    // Nonlinear function f(x) = x^3
    let x = TaylorModel::variable(0, vars.clone(), 3).unwrap();
    let f = x.mul(&x).unwrap().mul(&x).unwrap();
    let parent_range = f.range().unwrap();

    let (left, right) = f.subdivide_axis(0).unwrap();
    let left_range = left.range().unwrap();
    let right_range = right.range().unwrap();

    assert!(parent_range.encloses(left_range));
    assert!(parent_range.encloses(right_range));

    // 2D subdivision into 4 quadrants
    let vars2 = vec![
        VariableInfo::new("x", Interval::new(0.0, 2.0)).unwrap(),
        VariableInfo::new("y", Interval::new(0.0, 2.0)).unwrap(),
    ];
    let x2 = TaylorModel::variable(0, vars2.clone(), 2).unwrap();
    let y2 = TaylorModel::variable(1, vars2.clone(), 2).unwrap();
    let f2 = x2.mul(&y2).unwrap();

    let quads = f2.subdivide_all_axes().unwrap();
    assert_eq!(quads.len(), 4);

    let f2_range = f2.range().unwrap();
    for q in &quads {
        let q_range = q.range().unwrap();
        assert!(f2_range.encloses(q_range));
    }
}

#[test]
fn test_multivariate_composition() {
    // Outer function: f(u) = u^2
    let u_dom = Interval::new(0.0, 4.0);
    let outer_vars = vec![VariableInfo::new("u", u_dom).unwrap()];
    let u = TaylorModel::variable(0, outer_vars.clone(), 2).unwrap();
    let f_outer = u.mul(&u).unwrap();

    // Inner function: g(x, y) = x + y on [0, 1] x [0, 1] (range [0, 2] ⊆ [0, 4])
    let in_vars = vec![
        VariableInfo::new("x", Interval::new(0.0, 1.0)).unwrap(),
        VariableInfo::new("y", Interval::new(0.0, 1.0)).unwrap(),
    ];
    let x = TaylorModel::variable(0, in_vars.clone(), 2).unwrap();
    let y = TaylorModel::variable(1, in_vars.clone(), 2).unwrap();
    let g_inner = x.add(&y).unwrap();

    // Composed model: (x + y)^2
    let f_composed = f_outer.compose(&[g_inner.clone()]).expect("compose");
    assert_eq!(f_composed.dim(), 2);

    // Direct multiplication model: (x + y) * (x + y)
    let f_direct = g_inner.mul(&g_inner).unwrap();

    // Compare evaluation at (0.3, 0.4): expected (0.7)^2 = 0.49
    let pt = vec![Interval::point(0.3), Interval::point(0.4)];
    let val_comp = f_composed.eval_box(&pt).unwrap();
    let val_dir = f_direct.eval_box(&pt).unwrap();

    assert!(val_comp.contains(0.49));
    assert!(val_dir.contains(0.49));
    assert!((val_comp.midpoint() - val_dir.midpoint()).abs() < 1e-14);
}

#[test]
fn test_differentiation_and_gradient() {
    let vars = vec![
        VariableInfo::new("x", Interval::new(0.0, 2.0)).unwrap(),
        VariableInfo::new("y", Interval::new(0.0, 2.0)).unwrap(),
    ];

    // f(x, y) = 3*x^2*y + 4*y^3
    let x = TaylorModel::variable(0, vars.clone(), 3).unwrap();
    let y = TaylorModel::variable(1, vars.clone(), 3).unwrap();

    let x2 = x.mul(&x).unwrap();
    let x2_y = x2.mul(&y).unwrap().scale(3.0).unwrap();

    let y2 = y.mul(&y).unwrap();
    let y3 = y2.mul(&y).unwrap().scale(4.0).unwrap();

    let f = x2_y.add(&y3).unwrap();

    // Partial derivatives
    let df_dx = f.diff(0).expect("diff x");
    let df_dy = f.diff(1).expect("diff y");

    // Evaluate at point (1.5, 0.5):
    // df/dx = 6*x*y = 6 * 1.5 * 0.5 = 4.5
    // df/dy = 3*x^2 + 12*y^2 = 3 * 2.25 + 12 * 0.25 = 6.75 + 3 = 9.75
    let pt = vec![Interval::point(1.5), Interval::point(0.5)];
    let val_df_dx = df_dx.eval_box(&pt).unwrap();
    let val_df_dy = df_dy.eval_box(&pt).unwrap();

    assert!(val_df_dx.contains(4.5));
    assert!((val_df_dx.midpoint() - 4.5).abs() < 1e-13);

    assert!(val_df_dy.contains(9.75));
    assert!((val_df_dy.midpoint() - 9.75).abs() < 1e-13);

    // Gradient
    let grad = f.gradient().unwrap();
    assert_eq!(grad.len(), 2);

    // Out of bounds axis refuses
    let err_axis = f.diff(5).expect_err("out of bounds axis must refuse");
    assert!(matches!(
        err_axis,
        TaylorModelError::AxisIndexOutOfBounds { .. }
    ));
}

#[test]
fn test_degenerate_fixed_axis_differentiation() {
    let vars = vec![
        VariableInfo::new("x", Interval::new(0.0, 2.0)).unwrap(),
        VariableInfo::new("param", Interval::point(42.0)).unwrap(),
    ];

    let x = TaylorModel::variable(0, vars.clone(), 2).unwrap();
    let p = TaylorModel::variable(1, vars.clone(), 2).unwrap();
    let f = x.mul(&p).unwrap();

    // Derivative along fixed axis 'param' must be exact zero
    let df_dp = f.diff(1).unwrap();
    let pt = vec![Interval::point(1.0), Interval::point(42.0)];
    let val = df_dp.eval_box(&pt).unwrap();
    assert!(val.contains(0.0));
    assert!(val.width() < 1e-15);
}

#[test]
fn test_lipschitz_bound_and_integration() {
    let vars = vec![VariableInfo::new("x", Interval::new(0.0, 2.0)).unwrap()];

    // f(x) = 3*x^2 on [0, 2]
    let x = TaylorModel::variable(0, vars.clone(), 2).unwrap();
    let f = x.mul(&x).unwrap().scale(3.0).unwrap();

    // Lipschitz bound: max |f'(x)| on [0, 2] is f'(2) = 6*(2) = 12
    let l_bound = f.lipschitz_bound().unwrap();
    assert!(l_bound >= 12.0);

    // Coordinate integration: integral of 3*x^2 on [0, 2] is x^3|_0^2 = 8.0
    let integral = f.integrate_axis(0, 0.0, 2.0).unwrap();
    assert!(integral.contains(8.0));
    assert!((integral.midpoint() - 8.0).abs() < 1e-13);
}
