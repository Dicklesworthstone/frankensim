//! Taylor-model + certified-root battery: functional containment,
//! O(wⁿ⁺¹) subdivision convergence (the reason TMs exist), interval
//! Newton/Krawczyk certification incl. the double-root honesty case,
//! Lipschitz extraction, and the cross-ISA golden hash.

use fs_ivl::{
    Interval, RootBox, RootSearchConfig, RootSearchError, TaylorModel1, krawczyk_step,
    lipschitz_bound, newton_roots, newton_roots_bounded,
};

/// f(x) = x·sin(x) + exp(x·0.3) as a Taylor model (order 5).
fn model_f(domain: Interval) -> TaylorModel1 {
    let x = TaylorModel1::variable(domain, 5);
    let xs = &x * &x.sin();
    &xs + &x.scale(0.3).exp()
}

fn point_f(x: f64) -> f64 {
    x * fs_math::det::sin(x) + fs_math::det::exp(0.3 * x)
}

#[test]
fn functional_containment() {
    // The containment law extended to functions: model bound over any
    // subbox must contain point evaluations inside it.
    let domain = Interval::new(-0.8, 1.2);
    let tm = model_f(domain);
    for k in 0..=40 {
        let x = -0.8 + 2.0 * f64::from(k) / 40.0;
        let enc = tm.eval_interval(Interval::point(x));
        let truth = point_f(x);
        assert!(
            enc.contains(truth)
                || enc.contains(fs_math::next_down(truth))
                || enc.contains(fs_math::next_up(truth)),
            "containment violated at {x}: {enc:?} vs {truth}"
        );
    }
    println!(
        "{{\"suite\":\"fs-ivl\",\"case\":\"tm-containment\",\"verdict\":\"pass\",\"detail\":\"41 point checks inside model enclosures\"}}"
    );
}

#[test]
fn subdivision_convergence_beats_intervals() {
    // The right metric is the REMAINDER (the model's overestimation):
    // total bound width is dominated by the function's genuine O(w)
    // variation over the box, for TMs and intervals alike. The TM claim
    // is that the EXCESS shrinks like O(wⁿ⁺¹) where interval
    // arithmetic's excess shrinks like O(w).
    let widths = [0.4f64, 0.2, 0.1, 0.05];
    let mut rem_w = Vec::new();
    let mut ia_excess = Vec::new();
    for &w in &widths {
        let box_ = Interval::new(0.4 - w / 2.0, 0.4 + w / 2.0);
        let tm = model_f(box_);
        rem_w.push(tm.remainder().width());
        // Interval-arithmetic excess: enclosure width minus the true
        // range width (sampled densely — adequate for a smooth function).
        let x = box_;
        let ia = x * x.sin() + (x * Interval::point(0.3)).exp();
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        for k in 0..=100 {
            let v = point_f(box_.lo() + w * f64::from(k) / 100.0);
            lo = lo.min(v);
            hi = hi.max(v);
        }
        ia_excess.push(ia.width() - (hi - lo));
    }
    // Remainder must shrink SUPERLINEARLY: theoretical gain per halving
    // at order 5 is 2⁶ = 64; demand > 20 (slack for low-order cross
    // terms), and the final remainder is far below the IA excess.
    for i in 1..widths.len() {
        let gain = rem_w[i - 1] / rem_w[i].max(1e-300);
        assert!(
            gain > 20.0,
            "remainder gain per halving only {gain:.2} (rem widths {rem_w:?})"
        );
    }
    // NOTE: model_f is a POOR fixture for beating interval arithmetic —
    // every subterm is monotone and single-occurrence, so IA has no
    // dependency problem there (its excess measured ~2e-15). The decisive
    // comparison needs a dependency-heavy fixture: x − x² around 0.4,
    // where IA overestimates by O(w) and an order-5 TM's remainder is
    // pure rounding.
    let _ = ia_excess;
    let w = 0.05;
    let box_ = Interval::new(0.4 - w / 2.0, 0.4 + w / 2.0);
    let x = TaylorModel1::variable(box_, 5);
    let p = &x - &(&x * &x);
    let tm_excess = p.remainder().width();
    let ia = box_ - box_ * box_;
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for k in 0..=100 {
        let t = box_.lo() + w * f64::from(k) / 100.0;
        let v = t - t * t;
        lo = lo.min(v);
        hi = hi.max(v);
    }
    let dep_excess = ia.width() - (hi - lo);
    assert!(
        dep_excess > 0.01,
        "IA must visibly overestimate the dependency fixture: {dep_excess}"
    );
    assert!(
        tm_excess < dep_excess / 1e6,
        "TM excess {tm_excess} not decisively below IA excess {dep_excess}"
    );
    println!(
        "{{\"suite\":\"fs-ivl\",\"case\":\"tm-convergence\",\"verdict\":\"pass\",\"detail\":\"remainders {rem_w:?}; dependency fixture: tm {tm_excess:.2e} vs ia {dep_excess:.2e}\"}}"
    );
}

#[test]
fn newton_certifies_simple_roots() {
    // x² − 2 on [0, 3]: one certified root at √2.
    let f = |x: Interval| x * x - Interval::point(2.0);
    let fp = |x: Interval| Interval::point(2.0) * x;
    let roots = newton_roots(&f, &fp, Interval::new(0.0, 3.0), 1e-10);
    assert_eq!(roots.len(), 1, "{roots:?}");
    assert!(
        roots[0].is_certified(),
        "sqrt2 must be certified: {roots:?}"
    );
    let bx = roots[0].interval();
    assert!(
        bx.contains(std::f64::consts::SQRT_2),
        "box must contain sqrt2"
    );
    assert!(bx.width() < 1e-9, "certified box should be tight: {bx:?}");
    // sin on [2, 7]: roots at π and 2π, both certified.
    let fs = |x: Interval| x.sin();
    let fc = |x: Interval| x.cos();
    let roots = newton_roots(&fs, &fc, Interval::new(2.0, 7.0), 1e-10);
    let certified: Vec<f64> = roots
        .iter()
        .filter(|r| r.is_certified())
        .map(|r| r.interval().midpoint())
        .collect();
    assert_eq!(certified.len(), 2, "{roots:?}");
    assert!((certified[0] - std::f64::consts::PI).abs() < 1e-8);
    assert!((certified[1] - 2.0 * std::f64::consts::PI).abs() < 1e-8);
    println!(
        "{{\"suite\":\"fs-ivl\",\"case\":\"newton-certify\",\"verdict\":\"pass\",\"detail\":\"sqrt2 + sin roots certified, boxes tight\"}}"
    );
}

#[test]
fn double_root_is_never_falsely_certified() {
    // x² at 0: a double root — no interval-Newton certificate EXISTS
    // (the derivative enclosure straddles zero); the honest answer is
    // Possible boxes around 0, never Certified.
    let f = |x: Interval| x * x;
    let fp = |x: Interval| Interval::point(2.0) * x;
    let roots = newton_roots(&f, &fp, Interval::new(-1.0, 1.0), 1e-6);
    assert!(!roots.is_empty(), "the root region must be reported");
    for r in &roots {
        assert!(
            !r.is_certified(),
            "double root falsely certified: {roots:?}"
        );
        assert!(r.interval().contains(0.0) || r.interval().width() <= 2e-6);
    }
    println!(
        "{{\"suite\":\"fs-ivl\",\"case\":\"newton-honesty\",\"verdict\":\"pass\",\"detail\":\"x^2 double root -> Possible only, never Certified\"}}"
    );
}

#[test]
fn krawczyk_agrees_with_newton() {
    // Krawczyk on a box around √2 must certify and contract.
    let f = |x: Interval| x * x - Interval::point(2.0);
    let fp = |x: Interval| Interval::point(2.0) * x;
    let x0 = Interval::new(1.3, 1.5);
    let (k1, strict) = krawczyk_step(&f, &fp, x0).expect("nonempty");
    assert!(strict, "Krawczyk must certify on this box");
    assert!(k1.contains(std::f64::consts::SQRT_2));
    assert!(k1.width() < x0.width(), "must contract");
}

#[test]
fn root_search_rejects_non_finite_or_unbounded_requests() {
    let zero = |_x: Interval| Interval::point(0.0);
    let config = RootSearchConfig {
        min_width: 1.0,
        max_boxes: 8,
    };
    assert_eq!(
        newton_roots_bounded(&zero, &zero, Interval::WHOLE, config),
        Err(RootSearchError::NonFiniteDomain)
    );
    assert_eq!(
        newton_roots_bounded(
            &zero,
            &zero,
            Interval::new(-1.0, 1.0),
            RootSearchConfig {
                min_width: f64::NAN,
                max_boxes: 8,
            },
        ),
        Err(RootSearchError::InvalidMinWidth)
    );
    assert_eq!(
        newton_roots_bounded(
            &zero,
            &zero,
            Interval::new(-1.0, 1.0),
            RootSearchConfig {
                min_width: 1.0,
                max_boxes: 0,
            },
        ),
        Err(RootSearchError::EmptyBudget)
    );
}

#[test]
fn root_search_budget_exhaustion_is_reported_as_possible() {
    let zero = |_x: Interval| Interval::point(0.0);
    let domain = Interval::new(-1.0, 1.0);
    let report = newton_roots_bounded(
        &zero,
        &zero,
        domain,
        RootSearchConfig {
            min_width: f64::MIN_POSITIVE,
            max_boxes: 8,
        },
    )
    .expect("valid bounded search");
    assert_eq!(report.boxes_examined, 8);
    assert!(!report.complete);
    assert!(report.roots.iter().all(|root| !root.is_certified()));
    assert!(report.roots.iter().any(
        |root| matches!(root, RootBox::Possible(iv) if iv.lo() <= domain.lo() && iv.hi() >= domain.hi())
    ));
}

#[test]
fn root_search_stops_when_midpoint_cannot_split_adjacent_floats() {
    let zero = |_x: Interval| Interval::point(0.0);
    let lo = 1.0;
    let hi = fs_math::next_up(lo);
    let report = newton_roots_bounded(
        &zero,
        &zero,
        Interval::new(lo, hi),
        RootSearchConfig {
            min_width: f64::MIN_POSITIVE,
            max_boxes: 8,
        },
    )
    .expect("valid bounded search");
    assert!(report.complete);
    assert_eq!(report.boxes_examined, 1);
    assert_eq!(report.roots, vec![RootBox::Possible(Interval::new(lo, hi))]);
}

#[test]
fn lipschitz_bounds_are_certified_and_tight() {
    // sin on [0, π/2]: sup|cos| = 1 (at 0). Bound ≥ 1, and tight.
    let fc = |x: Interval| x.cos();
    let l = lipschitz_bound(&fc, Interval::new(0.0, std::f64::consts::FRAC_PI_2));
    assert!(l >= 1.0, "must enclose the true Lipschitz constant 1: {l}");
    assert!(l < 1.0 + 1e-9, "should be tight on this monotone case: {l}");
    // Cubic x³ on [-2, 2]: sup|3x²| = 12.
    let fp = |x: Interval| Interval::point(3.0) * x * x;
    let l3 = lipschitz_bound(&fp, Interval::new(-2.0, 2.0));
    assert!((12.0..12.0 + 1e-9).contains(&l3), "cubic Lipschitz {l3}");
    println!(
        "{{\"suite\":\"fs-ivl\",\"case\":\"lipschitz\",\"verdict\":\"pass\",\"detail\":\"sin L={l:.3}, cubic L={l3:.3}\"}}"
    );
}

/// Recorded on aarch64-apple (M4 Pro); must match on x86-64 (trj).
const GOLDEN_HASH: u64 = 0xa16d_9759_a352_95d8;

#[test]
fn taylor_golden_hash() {
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    let mut feed = |v: f64| {
        for byte in v.to_bits().to_le_bytes() {
            acc ^= u64::from(byte);
            acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    let tm = model_f(Interval::new(-0.8, 1.2));
    let b = tm.bound();
    feed(b.lo());
    feed(b.hi());
    feed(tm.remainder().lo());
    feed(tm.remainder().hi());
    let f = |x: Interval| x * x - Interval::point(2.0);
    let fp = |x: Interval| Interval::point(2.0) * x;
    for r in newton_roots(&f, &fp, Interval::new(0.0, 3.0), 1e-10) {
        feed(r.interval().lo());
        feed(r.interval().hi());
        feed(if r.is_certified() { 1.0 } else { 0.0 });
    }
    feed(lipschitz_bound(
        &|x: Interval| x.cos(),
        Interval::new(0.0, std::f64::consts::FRAC_PI_2),
    ));
    println!(
        "{{\"suite\":\"fs-ivl\",\"case\":\"taylor-golden\",\"verdict\":\"info\",\"detail\":\"{acc:#018x}\"}}"
    );
    assert_eq!(
        acc, GOLDEN_HASH,
        "taylor/newton bits changed: {acc:#018x} vs {GOLDEN_HASH:#018x} — bump only with \
         semantic justification (golden-evidence policy)"
    );
}
