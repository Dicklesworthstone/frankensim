use super::*;
use std::cell::Cell;

fn options(n: usize, nodes: usize) -> BeynOptions {
    BeynOptions {
        nodes,
        probe_columns: n,
        rank_relative_tolerance: 1e-9,
        count_tolerance: 0.1,
        quadrature_tolerance: 0.1,
        contour_min_reciprocal_condition: 1e-12,
        linear_residual_tolerance: 1e-12,
        svd_tolerance: 1e-13,
        svd_max_sweeps: 100,
        eigen_max_sweeps: 120,
        max_dimension: n,
        max_dense_bytes: 2_000_000,
        max_work: 1_000_000_000,
        refinement: PoleOptions {
            spectral_scale: 1.0,
            step_tolerance: 1e-13,
            residual_tolerance: 1e-13,
            max_iterations: 24,
            max_dimension: n,
            max_dense_bytes: 2_000_000,
            max_work: 100_000_000,
        },
    }
}

fn circle(center: C64, radius: f64) -> Ellipse {
    Ellipse {
        center,
        real_radius: radius,
        imaginary_radius: radius,
    }
}

// The nonnormal third direction is outside the target region. No symmetry
// assumption enters the nonlinear, LU, SVD or general QR owners.
fn nonnormal(s: C64, _: &mut dyn FnMut() -> bool) -> Result<PoleEvaluation, NonlinearEigError> {
    Ok(PoleEvaluation {
        matrix: vec![
            s - C64::new(-0.4, 0.1),
            C64::new(0.3, 0.2),
            C64::from_re(0.1),
            C64::ZERO,
            s - C64::new(0.25, -0.15),
            C64::new(-0.2, 0.1),
            C64::ZERO,
            C64::ZERO,
            s - C64::from_re(4.0),
        ],
        derivative: vec![
            C64::ONE,
            C64::ZERO,
            C64::ZERO,
            C64::ZERO,
            C64::ONE,
            C64::ZERO,
            C64::ZERO,
            C64::ZERO,
            C64::ONE,
        ],
        residual_scale: 6.0 + s.abs(),
    })
}

#[test]
fn g1_contour_discovers_nonnormal_modes_without_individual_seeds() {
    let result = beyn_with_cancel(
        3,
        circle(C64::ZERO, 0.9),
        &options(3, 64),
        nonnormal,
        || false,
    )
    .unwrap();
    assert_eq!(
        (result.argument_count, result.rank, result.poles.len()),
        (2, 2, 2)
    );
    for exact in [C64::new(-0.4, 0.1), C64::new(0.25, -0.15)] {
        let mode = result
            .poles
            .iter()
            .min_by(|a, b| (a.value - exact).abs().total_cmp(&(b.value - exact).abs()))
            .unwrap();
        assert!((mode.value - exact).abs() < 1e-12, "{result:?}");
        assert!(mode.residual < 1e-12);
        assert_eq!(mode.vector.len(), 3);
    }
    assert!((result.argument_integral - C64::from_re(2.0)).abs() < 1e-10);
}

#[test]
fn g1_quadratic_modes_match_independent_companion_spectrum() {
    let k = [5.0, -1.0, -1.0, 3.0];
    let c = [0.2, 0.05, 0.05, 0.1];
    let m = [2.0, 0.0, 0.0, 1.0];
    let operator = |s: C64, _: &mut dyn FnMut() -> bool| {
        Ok(PoleEvaluation {
            matrix: (0..4)
                .map(|i| C64::from_re(k[i]) + s.scale(c[i]) + (s * s).scale(m[i]))
                .collect(),
            derivative: (0..4)
                .map(|i| C64::from_re(c[i]) + s.scale(2.0 * m[i]))
                .collect(),
            residual_scale: 8.0 + s.abs() + 3.0 * s.abs() * s.abs(),
        })
    };
    let companion = vec![
        C64::ZERO,
        C64::ZERO,
        C64::ONE,
        C64::ZERO,
        C64::ZERO,
        C64::ZERO,
        C64::ZERO,
        C64::ONE,
        C64::from_re(-2.5),
        C64::from_re(0.5),
        C64::from_re(-0.1),
        C64::from_re(-0.025),
        C64::ONE,
        C64::from_re(-3.0),
        C64::from_re(-0.05),
        C64::from_re(-0.1),
    ];
    let truth = crate::eigen_complex::eig(&companion, 4).unwrap();
    let contour = Ellipse {
        center: C64::new(-0.05, 1.6),
        real_radius: 0.7,
        imaginary_radius: 0.8,
    };
    let result = beyn_with_cancel(2, contour, &options(2, 64), operator, || false).unwrap();
    assert_eq!(result.argument_count, 2);
    for exact in truth.into_iter().filter(|&z| contour.contains(z)) {
        assert!(
            result.poles.iter().any(|p| (p.value - exact).abs() < 1e-11),
            "{result:?}"
        );
    }
}

fn lambert(s: C64, _: &mut dyn FnMut() -> bool) -> Result<PoleEvaluation, NonlinearEigError> {
    let (sine, cosine) = det::sin_cos(-s.im);
    let exponential = C64::new(cosine, sine).scale(det::exp(-s.re));
    let a = [-0.45, 0.05, 0.42];
    let b = [0.08, 0.05, 0.06];
    let q = [0.6, -0.8, 0.0, 0.8, 0.6, 0.0, 0.0, 0.0, 1.0];
    let mut matrix = vec![C64::ZERO; 9];
    let mut derivative = vec![C64::ZERO; 9];
    for i in 0..3 {
        for j in 0..3 {
            for t in 0..3 {
                let weight = q[i * 3 + t] * q[j * 3 + t];
                matrix[i * 3 + j] = matrix[i * 3 + j]
                    + (C64::from_re(a[t]) - s + exponential.scale(b[t])).scale(weight);
                derivative[i * 3 + j] =
                    derivative[i * 3 + j] + (-C64::ONE - exponential.scale(b[t])).scale(weight);
            }
        }
    }
    Ok(PoleEvaluation {
        matrix,
        derivative,
        residual_scale: 3.0 + s.abs() + exponential.abs(),
    })
}

#[test]
fn g1_lambert_w_spectrum_has_a_live_exponential_quadrature_ladder() {
    // Independent scalar Lambert-W principal branch: w*exp(w)=b*exp(-a).
    let mut exact = Vec::new();
    for (a, b) in [(-0.45_f64, 0.08_f64), (0.05, 0.05), (0.42, 0.06)] {
        let target = b * (-a).exp();
        let mut w = target;
        for _ in 0..20 {
            let ew = w.exp();
            w -= (w * ew - target) / (ew * (1.0 + w));
        }
        assert!((w * w.exp() - target).abs() < 1e-15);
        exact.push(C64::from_re(a + w));
    }
    let mut errors = Vec::new();
    for nodes in [16, 32, 64] {
        // A broad contour brings the nonprincipal Lambert branches near
        // enough to make the analytic-background quadrature error measurable.
        // The same three principal roots remain the only enclosed roots.
        let result = beyn_with_cancel(
            3,
            circle(C64::ZERO, 3.2),
            &options(3, nodes),
            lambert,
            || false,
        )
        .unwrap();
        assert_eq!((result.argument_count, result.rank), (3, 3));
        let error = result
            .unrefined_values
            .iter()
            .map(|&z| {
                exact
                    .iter()
                    .map(|&v| (z - v).abs())
                    .fold(f64::INFINITY, f64::min)
            })
            .fold(0.0, f64::max);
        errors.push(error);
        for target in &exact {
            assert!(
                result
                    .poles
                    .iter()
                    .any(|p| (p.value - *target).abs() < 1e-11 && p.residual < 1e-12),
                "{result:?}"
            );
        }
        println!(
            "nodes={nodes}, raw_eigenvalue_error={error:e}, moment_change={:e}, count={:?}",
            result.moment_relative_change, result.argument_integral
        );
    }
    assert!(
        errors[0] > 1e-9,
        "ladder must exercise actual quadrature error: {errors:?}"
    );
    assert!(
        errors[1] < errors[0] * 0.05 && errors[2] < errors[1] * 0.05,
        "{errors:?}"
    );
}

#[test]
fn g1_singular_descriptor_keeps_its_finite_mode_and_algebraic_direction() {
    let operator = |s: C64, _: &mut dyn FnMut() -> bool| {
        Ok(PoleEvaluation {
            matrix: vec![
                s + C64::from_re(2.0),
                C64::from_re(-1.0),
                C64::from_re(-3.0),
                C64::from_re(4.0),
            ],
            derivative: vec![C64::ONE, C64::ZERO, C64::ZERO, C64::ZERO],
            residual_scale: 7.0 + s.abs(),
        })
    };
    let result = beyn_with_cancel(
        2,
        circle(C64::from_re(-1.3), 0.7),
        &options(2, 64),
        operator,
        || false,
    )
    .unwrap();
    assert_eq!(result.argument_count, 1);
    let pole = &result.poles[0];
    assert!((pole.value - C64::from_re(-1.25)).abs() < 1e-12);
    assert!((pole.vector[1] - pole.vector[0].scale(0.75)).abs() < 1e-12);
}

#[test]
fn g3_uniform_coefficient_scaling_does_not_change_the_discovered_spectrum() {
    for scale in [1e-250, 1.0, 1e250] {
        let result = beyn_with_cancel(
            3,
            circle(C64::ZERO, 0.9),
            &options(3, 64),
            |s, cancel| {
                let mut v = nonnormal(s, cancel)?;
                for z in v.matrix.iter_mut().chain(&mut v.derivative) {
                    *z = z.scale(scale);
                }
                v.residual_scale *= scale;
                Ok(v)
            },
            || false,
        )
        .unwrap();
        assert_eq!(result.argument_count, 2);
        assert!((result.poles[0].value - C64::new(-0.4, 0.1)).abs() < 1e-11);
    }
}

#[test]
fn g0_contour_and_probe_refusals_do_not_silently_drop_modes() {
    let mut opts = options(3, 64);
    opts.probe_columns = 2;
    assert!(matches!(
        beyn_with_cancel(3, circle(C64::ZERO, 0.9), &opts, nonnormal, || false),
        Err(BeynError::ProbeRankSaturated { count: 2, .. })
    ));
    opts.probe_columns = 1;
    assert!(matches!(
        beyn_with_cancel(3, circle(C64::ZERO, 0.9), &opts, nonnormal, || false),
        Err(BeynError::ProbeRankSaturated { count: 2, .. })
    ));
    let scalar = |z: C64, _: &mut dyn FnMut() -> bool| {
        Ok(PoleEvaluation {
            matrix: vec![z - C64::ONE],
            derivative: vec![C64::ONE],
            residual_scale: 1.0 + z.abs(),
        })
    };
    assert!(matches!(
        beyn_with_cancel(1, circle(C64::ZERO, 1.0), &options(1, 32), scalar, || false),
        Err(BeynError::ContourSingular { .. })
    ));
    let mut near = options(1, 32);
    near.contour_min_reciprocal_condition = 1e-8;
    assert!(matches!(
        beyn_with_cancel(1, circle(C64::ZERO, 1.0 + 1e-10), &near, scalar, || false),
        Err(BeynError::ContourIllConditioned { .. })
    ));
    let empty = beyn_with_cancel(
        3,
        circle(C64::from_re(10.0), 0.5),
        &options(3, 32),
        nonnormal,
        || false,
    )
    .unwrap();
    assert_eq!(empty.argument_count, 0);
    assert!(empty.poles.is_empty());
}

#[test]
fn g4_budget_and_cancellation_refuse_before_publication_and_replay_exactly() {
    let contour = circle(C64::ZERO, 0.9);
    let mut opts = options(3, 64);
    let evaluations = Cell::new(0);
    opts.max_dense_bytes = 0;
    assert!(matches!(
        beyn_with_cancel(
            3,
            contour,
            &opts,
            |z, c| {
                evaluations.set(evaluations.get() + 1);
                nonnormal(z, c)
            },
            || false
        ),
        Err(BeynError::Invalid { .. })
    ));
    assert_eq!(evaluations.get(), 0);
    opts = options(3, 64);
    let polls = Cell::new(0);
    assert!(matches!(
        beyn_with_cancel(3, contour, &opts, nonnormal, || {
            polls.set(polls.get() + 1);
            polls.get() > 100
        }),
        Err(BeynError::Operator(NonlinearEigError::Cancelled))
    ));
    let a = beyn_with_cancel(3, contour, &opts, nonnormal, || false).unwrap();
    let b = beyn_with_cancel(3, contour, &opts, nonnormal, || false).unwrap();
    assert_eq!(a.argument_integral, b.argument_integral);
    for (a, b) in a.poles.iter().zip(&b.poles) {
        assert_eq!(a.value, b.value);
        assert_eq!(a.vector, b.vector);
        assert_eq!(a.residual.to_bits(), b.residual.to_bits());
    }
}
