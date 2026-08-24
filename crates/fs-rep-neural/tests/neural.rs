//! Battery for neural implicit charts (fs-rep-neural). The load-bearing tests
//! are the CERTIFICATE ones: the spectral-norm upper bound never understates
//! the true norm, the certified Lipschitz constant is never violated by
//! sampling, interval bound propagation is sound over a box, the gradient is
//! bounded by L, and a sphere-trace step of |f|/L never tunnels through the
//! surface.

use fs_rep_neural::{
    EvaluationInput, InputDimensionError, Layer, MLP_ACTIVATION_SEMANTICS,
    MLP_ACTIVATION_SEMANTICS_VERSION, MLP_ACTIVATION_ULP_BUDGET, MLP_FIELD_IDENTITY_SCHEMA_VERSION,
    MlpSdf, SAFE_STEP_POLICY, SAFE_STEP_POLICY_VERSION, SafeStepStatus, TopologyHint,
    derive_safe_step, spectral_norm, spectral_norm_upper_bound, spectral_normalize,
};

// A deterministic pseudo-random point stream in [-1, 1)^2 (no rand crate).
fn points(n: usize, seed: u64) -> Vec<[f64; 2]> {
    let mut s = seed;
    let mut next = || {
        s = s
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        // `s >> 32` ∈ [0, 2^32); /2^31 ∈ [0, 2); −1 ∈ [-1, 1). (The former
        // `s >> 33` gave [0, 2^31)/2^31 − 1 = [-1, 0) — only the third
        // quadrant, so the certificate tests silently probed one octant of the
        // domain. `L` is a GLOBAL Lipschitz bound, so full-domain coverage is
        // still sound.)
        ((s >> 32) as f64 / f64::from(1u32 << 31)) - 1.0
    };
    (0..n).map(|_| [next(), next()]).collect()
}

fn sample_mlp() -> MlpSdf {
    let l1 = Layer::new(
        vec![vec![0.5, -0.3], vec![0.2, 0.7], vec![-0.4, 0.1]],
        vec![0.1, -0.2, 0.05],
    );
    let l2 = Layer::new(vec![vec![0.6, -0.5, 0.3]], vec![0.0]);
    MlpSdf::new(vec![l1, l2], 1.0)
}

#[test]
fn spectral_norm_matches_known_values() {
    // diagonal: σ_max is the largest |diagonal|.
    assert!((spectral_norm(&[vec![3.0, 0.0], vec![0.0, 1.0]]) - 3.0).abs() < 1e-9);
    // [[1,2],[3,4]] has σ_max = sqrt((30 + sqrt(884))/2) ≈ 5.464986.
    let s = spectral_norm(&[vec![1.0, 2.0], vec![3.0, 4.0]]);
    assert!((s - 5.464_985_704_219_04).abs() < 1e-6, "σ = {s}");
}

#[test]
fn spectral_normalization_certifies_at_most_the_bound() {
    let layer = Layer::new(vec![vec![1.0, 2.0], vec![3.0, 4.0]], vec![0.0, 0.0]);
    let normed = spectral_normalize(layer, 2.0);
    let certified = spectral_norm_upper_bound(&normed.weights);
    assert!(certified <= 2.0, "certified upper bound = {certified}");
    assert!(spectral_norm(&normed.weights) <= certified);
}

#[test]
fn certified_bound_covers_a_power_iteration_blind_direction() {
    // Recreate spectral_norm's normalized fixed starting direction, then make
    // a nonzero row orthogonal to it. The small scale keeps even a last-ULP
    // residual below power iteration's fixed cutoff; a safety certificate must
    // still see the matrix.
    let initial_norm = (1.0_f64 + 1.1_f64 * 1.1_f64).sqrt();
    let v0 = 1.0 / initial_norm;
    let v1 = 1.1 / initial_norm;
    let blind_scale = 1e-8;
    let weights = vec![vec![v1 * blind_scale, -v0 * blind_scale]];
    assert_eq!(spectral_norm(&weights).to_bits(), 0.0f64.to_bits());

    let exact = v0.hypot(v1) * blind_scale;
    let certified = spectral_norm_upper_bound(&weights);
    assert!(certified >= exact, "upper {certified} < exact {exact}");
    assert!(certified.is_finite());

    let normed = spectral_normalize(Layer::new(weights, vec![0.0]), 0.5);
    let normalized_exact = normed.weights[0][0].hypot(normed.weights[0][1]);
    let normalized_certified = spectral_norm_upper_bound(&normed.weights);
    assert!(normalized_exact <= normalized_certified);
    assert!(normalized_certified <= 0.5);
}

#[test]
fn certified_bound_preserves_sparse_axis_layer_scale() {
    // This is the hidden-layer structure used by NeuroShapeCert. Its exact
    // spectral norm is sqrt(18); the induced bound sees that structure while a
    // Frobenius-only bound would be 6 and unnecessarily alter the chart.
    let weights = vec![
        vec![3.0, 0.0],
        vec![-3.0, 0.0],
        vec![0.0, 3.0],
        vec![0.0, -3.0],
    ];
    let exact = (18.0_f64).sqrt();
    let certified = spectral_norm_upper_bound(&weights);
    assert!(certified >= exact);
    assert!(
        certified - exact < 1e-12,
        "overly conservative: {certified}"
    );
}

#[test]
fn certificate_handles_zero_rectangular_and_extreme_finite_matrices() {
    let zero = vec![vec![0.0; 3], vec![0.0; 3]];
    assert_eq!(spectral_norm_upper_bound(&zero).to_bits(), 0.0f64.to_bits());
    assert_eq!(
        spectral_normalize(Layer::new(zero.clone(), vec![0.0; 2]), 1.0).weights,
        zero
    );
    let collapsed = spectral_normalize(Layer::new(vec![vec![1.0, -2.0]], vec![0.0]), 0.0);
    assert_eq!(spectral_norm_upper_bound(&collapsed.weights).to_bits(), 0.0f64.to_bits());

    // A single maximum-finite entry has an exactly representable finite norm.
    assert_eq!(
        spectral_norm_upper_bound(&[vec![f64::MAX]]).to_bits(),
        f64::MAX.to_bits()
    );

    // Two maximum-finite entries have a mathematical norm greater than
    // f64::MAX, so a direct finite certificate is impossible and must fail
    // closed. Normalization operates on scaled entries and remains useful.
    assert!(
        std::panic::catch_unwind(|| { spectral_norm_upper_bound(&[vec![f64::MAX, f64::MAX]]) })
            .is_err()
    );
    let huge = spectral_normalize(Layer::new(vec![vec![f64::MAX, f64::MAX]], vec![0.0]), 1.0);
    assert!(huge.weights.iter().flatten().all(|w| w.is_finite()));
    assert!(huge.weights.iter().flatten().any(|w| *w != 0.0));
    assert!(spectral_norm_upper_bound(&huge.weights) <= 1.0);

    // Relative scaling also avoids overflow when growing subnormal weights.
    let tiny = spectral_normalize(
        Layer::new(vec![vec![f64::from_bits(1), f64::from_bits(1)]], vec![0.0]),
        1.0,
    );
    assert!(tiny.weights.iter().flatten().all(|w| w.is_finite()));
    assert!(tiny.weights.iter().flatten().any(|w| *w != 0.0));
    assert!(spectral_norm_upper_bound(&tiny.weights) <= 1.0);
}

#[test]
fn certificate_rejects_nonfinite_inputs_and_unrepresentable_products() {
    assert!(std::panic::catch_unwind(|| Layer::new(vec![vec![f64::NAN]], vec![0.0])).is_err());
    assert!(std::panic::catch_unwind(|| Layer::new(vec![vec![1.0]], vec![f64::INFINITY])).is_err());
    assert!(std::panic::catch_unwind(|| spectral_norm_upper_bound(&[vec![f64::NAN]])).is_err());
    assert!(
        std::panic::catch_unwind(|| {
            spectral_normalize(Layer::new(vec![vec![1.0]], vec![0.0]), f64::INFINITY)
        })
        .is_err()
    );
    assert!(
        std::panic::catch_unwind(|| {
            spectral_normalize(Layer::new(vec![vec![1.0]], vec![0.0]), -1.0)
        })
        .is_err()
    );

    let layers = vec![
        Layer::new(vec![vec![1.0]], vec![0.0]),
        Layer::new(vec![vec![1.0]], vec![0.0]),
    ];
    assert!(std::panic::catch_unwind(|| MlpSdf::new(layers, f64::MAX)).is_err());
}

#[test]
fn interval_sign_margin_step_rounds_toward_the_conservative_side() {
    // The nearest f64 representation of exact 1/10 lies above the rational
    // value. A safety radius must therefore step down one ulp.
    let nearest = 1.0_f64 / 10.0;
    let positive = derive_safe_step((1.0, 1.0), 10.0);
    assert_eq!(positive.status(), SafeStepStatus::SignSeparated);
    assert_eq!(positive.magnitude_lower_bound().to_bits(), 1.0f64.to_bits());
    assert_eq!(positive.radius().to_bits() + 1, nearest.to_bits());
    assert_eq!(positive.policy_version(), SAFE_STEP_POLICY_VERSION);
    assert_eq!(SAFE_STEP_POLICY_VERSION, 1);
    assert_eq!(SAFE_STEP_POLICY, "fs-rep-neural-interval-sign-margin-v1");
    assert_eq!(MLP_ACTIVATION_ULP_BUDGET, fs_math::det::TANH_ULP_BUDGET);

    let negative = derive_safe_step((-2.0, -1.0), 10.0);
    assert_eq!(negative.status(), SafeStepStatus::SignSeparated);
    assert_eq!(negative.magnitude_lower_bound().to_bits(), 1.0f64.to_bits());
    assert_eq!(negative.radius().to_bits(), positive.radius().to_bits());

    // A rounded nominal value cannot create authority when the sound
    // enclosure still straddles zero.
    let unresolved = derive_safe_step((-f64::from_bits(1), f64::from_bits(1)), 1.0);
    assert_eq!(unresolved.status(), SafeStepStatus::NoFiniteSignMargin);
    assert_eq!(unresolved.radius().to_bits(), 0.0f64.to_bits());

    // Invalid claimed bounds cannot create motion authority. A genuinely
    // nonzero constant field retains mathematically unbounded clearance, while
    // the identically-zero field has zero clearance.
    for enclosure in [(f64::NAN, 1.0), (1.0, f64::NAN), (1.0, -1.0)] {
        let invalid = derive_safe_step(enclosure, 1.0);
        assert_eq!(invalid.status(), SafeStepStatus::InvalidEnclosure);
        assert_eq!(invalid.radius().to_bits(), 0.0f64.to_bits());
    }
    for invalid_lipschitz in [f64::NAN, f64::INFINITY, -1.0] {
        let invalid = derive_safe_step((1.0, 1.0), invalid_lipschitz);
        assert_eq!(invalid.status(), SafeStepStatus::InvalidLipschitz);
        assert_eq!(invalid.radius().to_bits(), 0.0f64.to_bits());
    }
    assert_eq!(
        derive_safe_step((1.0, 1.0), 0.0)
            .radius()
            .to_bits(),
        f64::INFINITY.to_bits()
    );
    assert_eq!(
        derive_safe_step((-1.0, -1.0), -0.0)
            .radius()
            .to_bits(),
        f64::INFINITY.to_bits()
    );
    assert_eq!(
        derive_safe_step((0.0, 0.0), 0.0)
            .radius()
            .to_bits(),
        0.0f64.to_bits()
    );
    assert_eq!(
        derive_safe_step((-0.0, 0.0), -0.0)
            .radius()
            .to_bits(),
        0.0f64.to_bits()
    );
}

#[test]
fn interval_sign_margin_handles_extremes_scaling_and_one_sided_infinity() {
    let positive_tail = derive_safe_step((f64::MAX, f64::INFINITY), 1.0);
    assert_eq!(positive_tail.status(), SafeStepStatus::SignSeparated);
    assert_eq!(positive_tail.magnitude_lower_bound().to_bits(), f64::MAX.to_bits());
    assert_eq!(positive_tail.radius().to_bits() + 1, f64::MAX.to_bits());

    let negative_tail = derive_safe_step((f64::NEG_INFINITY, -f64::MAX), 1.0);
    assert_eq!(negative_tail.status(), SafeStepStatus::SignSeparated);
    assert_eq!(negative_tail.magnitude_lower_bound().to_bits(), f64::MAX.to_bits());
    assert_eq!(
        negative_tail.radius().to_bits(),
        positive_tail.radius().to_bits()
    );

    let whole = derive_safe_step((f64::NEG_INFINITY, f64::INFINITY), 1.0);
    assert_eq!(whole.status(), SafeStepStatus::NoFiniteSignMargin);
    assert_eq!(whole.radius().to_bits(), 0.0f64.to_bits());

    let underflow = derive_safe_step((f64::from_bits(1), f64::from_bits(2)), f64::MAX);
    assert_eq!(underflow.status(), SafeStepStatus::SignSeparated);
    assert_eq!(underflow.radius().to_bits(), 0.0f64.to_bits());

    let overflow = derive_safe_step((f64::MAX, f64::MAX), f64::MIN_POSITIVE);
    assert_eq!(overflow.radius().to_bits(), f64::MAX.to_bits());

    let baseline = derive_safe_step((0.75, 1.25), 3.0);
    let scaled = derive_safe_step((1.5, 2.5), 6.0);
    assert_eq!(baseline.radius().to_bits(), scaled.radius().to_bits());
    for _ in 0..32 {
        let replay = derive_safe_step(baseline.enclosure(), baseline.lipschitz_upper_bound());
        assert_eq!(replay.status(), baseline.status());
        assert_eq!(replay.radius().to_bits(), baseline.radius().to_bits());
    }
}

#[test]
fn the_lipschitz_certificate_is_never_violated() {
    let net = sample_mlp();
    let l = net.lipschitz();
    assert!(l > 0.0 && l <= 1.0 + 1e-6, "L = {l}"); // bound 1, 1-Lipschitz tanh
    let pts = points(300, 0x1234_5678);
    let mut max_ratio = 0.0_f64;
    for w in pts.windows(2) {
        let (a, b) = (w[0], w[1]);
        let dist = ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2)).sqrt();
        if dist < 1e-9 {
            continue;
        }
        let ratio = (net.eval(&a) - net.eval(&b)).abs() / dist;
        max_ratio = max_ratio.max(ratio);
    }
    // NO sampled violation of the certified bound.
    assert!(
        max_ratio <= l + 1e-6,
        "sampled ratio {max_ratio} exceeds L {l}"
    );
}

#[test]
fn interval_bound_propagation_is_sound() {
    let net = sample_mlp();
    let (lo, hi) = ([-0.5, -0.5], [0.5, 0.5]);
    let (blo, bhi) = net.eval_interval(&lo, &hi);
    assert!(blo <= bhi);
    // every point in the box has f(x) inside the guaranteed enclosure.
    for p in points(400, 0xBEEF) {
        let x = [p[0] * 0.5, p[1] * 0.5]; // remap into the box
        let v = net.eval(&x);
        assert!(v >= blo && v <= bhi, "f={v} escaped [{blo}, {bhi}]");
    }
    // A degenerate input box still accumulates outward-rounding slack. It need
    // not equal the separately grouped point evaluator, but it must contain it
    // without an arbitrary comparison epsilon.
    let x = [0.2, -0.1];
    let (dlo, dhi) = net.eval_interval(&x, &x);
    let point = net.eval(&x);
    assert!(
        dlo <= point && point <= dhi,
        "f={point} escaped [{dlo}, {dhi}]"
    );
}

#[test]
fn interval_bound_propagation_refuses_malformed_boxes() {
    let net = sample_mlp();

    for (lo, hi) in [
        ([f64::NAN, 0.0], [1.0, 1.0]),
        ([f64::NEG_INFINITY, 0.0], [1.0, 1.0]),
        ([0.0, 0.0], [f64::INFINITY, 1.0]),
        ([1.0, 0.0], [0.0, 1.0]),
    ] {
        assert_eq!(
            net.eval_interval(&lo, &hi),
            (f64::NEG_INFINITY, f64::INFINITY)
        );
    }

    assert!(std::panic::catch_unwind(|| net.eval_interval(&[0.0], &[0.0])).is_err());
    assert!(std::panic::catch_unwind(|| net.eval_interval(&[0.0, 0.0], &[0.0])).is_err());
}

#[test]
fn every_evaluator_requires_the_exact_declared_input_dimension() {
    let net = sample_mlp();
    assert_eq!(net.input_dim(), 2);

    for (point, actual) in [(&[][..], 0), (&[0.0][..], 1), (&[0.0, 0.0, 7.0][..], 3)] {
        assert_eq!(
            net.try_eval(point),
            Err(InputDimensionError {
                input: EvaluationInput::Point,
                expected: 2,
                actual,
            })
        );
        assert_eq!(
            net.try_eval_grad(point),
            Err(InputDimensionError {
                input: EvaluationInput::Gradient,
                expected: 2,
                actual,
            })
        );
        assert!(std::panic::catch_unwind(|| net.eval(point)).is_err());
        assert!(std::panic::catch_unwind(|| net.eval_grad(point)).is_err());
    }

    assert_eq!(
        net.try_eval_interval(&[0.0], &[0.0, 0.0]),
        Err(InputDimensionError {
            input: EvaluationInput::IntervalLower,
            expected: 2,
            actual: 1,
        })
    );
    assert_eq!(
        net.try_eval_interval(&[0.0, 0.0], &[0.0, 0.0, 0.0]),
        Err(InputDimensionError {
            input: EvaluationInput::IntervalUpper,
            expected: 2,
            actual: 3,
        })
    );

    let point = [0.25, -0.5];
    assert_eq!(
        net.try_eval(&point).unwrap().to_bits(),
        net.eval(&point).to_bits()
    );
    assert_eq!(net.try_eval_grad(&point).unwrap(), net.eval_grad(&point));
    assert_eq!(
        net.try_eval_interval(&point, &point).unwrap(),
        net.eval_interval(&point, &point)
    );

    let diagnostic = InputDimensionError {
        input: EvaluationInput::IntervalUpper,
        expected: 2,
        actual: 3,
    };
    assert_eq!(
        diagnostic.to_string(),
        "interval upper endpoint input dimension mismatch: expected 2, got 3"
    );
}

#[test]
fn interval_bound_propagation_is_bit_deterministic() {
    let net = sample_mlp();
    let lo = [-0.375, -0.625];
    let hi = [0.875, 0.125];
    let first = net.eval_interval(&lo, &hi);
    for _ in 0..32 {
        let replay = net.eval_interval(&lo, &hi);
        assert_eq!(replay.0.to_bits(), first.0.to_bits());
        assert_eq!(replay.1.to_bits(), first.1.to_bits());
    }
}

#[test]
fn neural_affine_overflow_remains_a_sign_separating_enclosure() {
    let net = MlpSdf::new(vec![Layer::new(vec![vec![1.0]], vec![f64::MAX])], 1.0);
    let enclosure = net.eval_interval(&[f64::MAX], &[f64::MAX]);
    assert!(enclosure.0.is_finite() && enclosure.0 > 0.0);
    assert_eq!(enclosure.1.to_bits(), f64::INFINITY.to_bits());
    let step = derive_safe_step(enclosure, net.lipschitz());
    assert_eq!(step.status(), SafeStepStatus::SignSeparated);
    assert!(step.radius() > 0.0);
}

#[test]
fn point_evaluation_uses_the_interval_certifiers_deterministic_tanh() {
    assert_eq!(MLP_ACTIVATION_SEMANTICS_VERSION, 1);
    assert_eq!(MLP_ACTIVATION_SEMANTICS, "fs-rep-neural-det-tanh-v1");

    // A zero hidden weight makes the hidden preactivation exactly the bias.
    // A scalar output layer is normalized to the same deterministic stored
    // weight whether normalized here or by MlpSdf::new, so this independently
    // pins the end-to-end point evaluator to fs_math::det::tanh rather than a
    // platform elementary function outside fs-ivl's declared ULP budget.
    let normalized_output = spectral_normalize(Layer::new(vec![vec![1.0]], vec![0.0]), 1.0);
    let output_weight = normalized_output.weights[0][0];
    assert!(output_weight > 0.0);

    for bias in [
        -19.0,
        -1.234_567_890_123_456,
        -0.73,
        -1.0e-8,
        1.0e-8,
        0.73,
        1.234_567_890_123_456,
        19.0,
    ] {
        let net = MlpSdf::new(
            vec![
                Layer::new(vec![vec![0.0]], vec![bias]),
                Layer::new(vec![vec![1.0]], vec![0.0]),
            ],
            1.0,
        );
        let expected = output_weight * fs_math::det::tanh(bias);
        let observed = net.eval(&[0.0]);
        assert_eq!(
            observed.to_bits(),
            expected.to_bits(),
            "point activation drift at bias {bias:.17e}"
        );

        let (lo, hi) = net.eval_interval(&[0.0], &[0.0]);
        assert!(
            lo <= observed && observed <= hi,
            "deterministic point {observed:.17e} escaped [{lo:.17e}, {hi:.17e}]"
        );
    }
}

#[test]
fn finite_difference_gradient_is_a_bounded_deterministic_diagnostic() {
    let net = sample_mlp();
    let l = net.lipschitz();
    for p in points(200, 0xC0FFEE) {
        let g = net.eval_grad(&p);
        let replay = net.eval_grad(&p);
        assert_eq!(
            g.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
            replay.iter().map(|x| x.to_bits()).collect::<Vec<_>>()
        );
        // Each centered coordinate secant is individually bounded by L in
        // exact arithmetic. Their vector does not represent one common-point
        // analytic gradient, so only the generic sqrt(d)*L diagnostic bound
        // is checked here, with allowance for rounded point evaluations.
        assert!(g.iter().all(|component| component.abs() <= l + 1e-4));
        let gnorm = (g[0] * g[0] + g[1] * g[1]).sqrt();
        assert!(
            gnorm <= 2.0_f64.sqrt() * l + 2e-4,
            "finite-difference diagnostic norm {gnorm} exceeds sqrt(d)*L for L={l}"
        );
    }
}

#[test]
fn a_certified_sphere_trace_step_never_tunnels() {
    let net = sample_mlp();
    let l = net.lipschitz();
    let x0 = [0.3, 0.2];
    let v = net.eval(&x0);
    let enclosure = net.eval_interval(&x0, &x0);
    let r = derive_safe_step(enclosure, l).radius();
    // f cannot change sign anywhere within radius r of x0 -> no tunneling.
    for k in 0..64 {
        let theta = f64::from(k) * std::f64::consts::TAU / 64.0;
        let y = [
            x0[0] + 0.999 * r * theta.cos(),
            x0[1] + 0.999 * r * theta.sin(),
        ];
        assert!(
            net.eval(&y).signum() == v.signum(),
            "tunneled at angle {theta}"
        );
    }
}

#[test]
fn topology_is_honestly_unknown() {
    assert_eq!(sample_mlp().topology_hint(), TopologyHint::Unknown);
}

#[test]
fn evaluation_is_deterministic() {
    let a = sample_mlp();
    let b = sample_mlp();
    let x = [0.4, -0.6];
    assert_eq!(a.eval(&x).to_bits(), b.eval(&x).to_bits());
    assert_eq!(a.lipschitz().to_bits(), b.lipschitz().to_bits());
}

#[test]
fn normalized_field_identity_binds_parameters_and_certificate_semantics() {
    assert_eq!(MLP_FIELD_IDENTITY_SCHEMA_VERSION, 1);
    let a = sample_mlp();
    let b = sample_mlp();
    assert_eq!(a.identity(), b.identity());

    let changed_bias = MlpSdf::new(
        vec![
            Layer::new(
                vec![vec![0.5, -0.3], vec![0.2, 0.7], vec![-0.4, 0.1]],
                vec![0.1, -0.2, 0.050_000_000_000_000_01],
            ),
            Layer::new(vec![vec![0.6, -0.5, 0.3]], vec![0.0]),
        ],
        1.0,
    );
    assert_ne!(a.identity(), changed_bias.identity());
    assert_eq!(fs_ivl::INTERVAL_SEMANTICS_VERSION, 1);
    assert_eq!(fs_math::STRICT_CORE_SEMANTICS_VERSION, 1);
    assert_eq!(fs_math::STRICT_CORE_GOLDEN_HASH, 0xeb79_cab7_a016_43e5);
}
