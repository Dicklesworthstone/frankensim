//! E4.0c battery (bead wf-root-guzez.5.1.3): coherent-draw coherence law
//! (same id = same surface bit-exact; the per-query-interval anti-pattern
//! is structurally impossible), strict OOD refusals over fitted boxes,
//! indicial kernel exact references (Wagner/Küssner values, sub-step
//! composition exactness, freeze-at-zero, reversed-flow refusal), caps at
//! cap AND cap+1, and a pinned Wagner-trace golden.
//! Repro: cargo test -p fs-airfoil --test uncertainty_indicial_battery

use fs_airfoil::fit::{BsplineAxis, FitSample, ResidualSurface};
use fs_airfoil::indicial::{
    IndicialKernel, IndicialState, KUSSNER_2POLE, MAX_DS, WAGNER_JONES, reduced_time_increment,
};
use fs_airfoil::table::{CoefficientTable, ConventionBlock, RegimePatch, SurfaceKind};
use fs_airfoil::uncertainty::{MAX_MODES, UncertainSurface};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-airfoil-e40c\",\"case\":\"{case}\",{payload}}}");
}

fn axis(name: &'static str, lo: f64, hi: f64, n: usize) -> BsplineAxis {
    BsplineAxis {
        name,
        lo,
        hi,
        n_coef: n,
    }
}

fn mean_surface() -> ResidualSurface {
    let axes = [
        axis("alpha_rad", -0.2, 0.3, 5),
        axis("log10_re", 5.0, 7.0, 4),
        axis("delta_rad", 0.0, 1.0, 1),
    ];
    let samples: Vec<FitSample> = (0..60)
        .map(|i| {
            let t = f64::from(i) / 59.0;
            let x = [
                -0.2 + 0.5 * t,
                5.0 + 2.0 * f64::from((i * 11) % 60) / 59.0,
                0.5,
            ];
            FitSample {
                x,
                y: 0.2 * x[0] + 0.01 * (x[1] - 6.0),
            }
        })
        .collect();
    ResidualSurface::fit(axes, vec![], &samples).unwrap()
}

#[test]
fn coherent_draw_law() {
    let mean = mean_surface();
    let n = mean.coef.len();
    // Two smooth modes (deterministic construction).
    let modes: Vec<Vec<f64>> = (0..2)
        .map(|k| {
            (0..n)
                .map(|i| 0.01 * f64::from(i as u32 % 7) * f64::from(k as u32 + 1))
                .collect()
        })
        .collect();
    let unc = UncertainSurface {
        mean: mean.clone(),
        modes,
    };
    // SAME realization id → bit-identical surface (coherence).
    let a = unc.realize("model-real-0001").unwrap();
    let b = unc.realize("model-real-0001").unwrap();
    assert_eq!(a, b, "one id must yield one surface, bit-exact");
    // DIFFERENT id → different surface.
    let c = unc.realize("model-real-0002").unwrap();
    assert_ne!(a.coef, c.coef, "distinct realizations must differ");
    // Coherence across queries: the draw shifts the WHOLE surface together —
    // the difference between two realizations at two query points has the
    // deterministic mode structure, not independent noise. Verify the shift
    // at two points matches the mode-combination prediction exactly.
    let q1 = [0.05, 6.0, 0.5];
    let q2 = [0.21, 5.4, 0.5];
    let w_a = UncertainSurface::weights("model-real-0001", 2);
    let mode_surface = |k: usize| ResidualSurface {
        axes: mean.axes.clone(),
        coef: unc.modes[k].clone(),
        constraints: vec![],
    };
    for q in [q1, q2] {
        let predicted =
            mean.eval(q) + w_a[0] * mode_surface(0).eval(q) + w_a[1] * mode_surface(1).eval(q);
        assert!(
            (a.eval(q) - predicted).abs() < 1e-12,
            "coherent structure broke at {q:?}"
        );
    }
    // Anonymous draws refuse; mode caps at MAX_MODES and MAX_MODES+1.
    assert_eq!(unc.realize("").unwrap_err().code, "realization-id-empty");
    let capped = UncertainSurface {
        mean: mean.clone(),
        modes: vec![vec![0.0; n]; MAX_MODES],
    };
    assert!(capped.realize("id").is_ok());
    let over = UncertainSurface {
        mean: mean.clone(),
        modes: vec![vec![0.0; n]; MAX_MODES + 1],
    };
    assert_eq!(
        over.realize("id").unwrap_err().code,
        "uncertainty-modes-invalid"
    );
    // Shape falsifier: wrong-length mode refuses.
    let bad = UncertainSurface {
        mean,
        modes: vec![vec![0.0; n - 1]],
    };
    assert_eq!(
        bad.realize("id").unwrap_err().code,
        "uncertainty-modes-invalid"
    );
    jlog(
        "coherent-draw",
        "\"law\":\"same-id-same-surface\",\"falsifiers\":\"executed\"",
    );
}

#[test]
fn strict_ood_refusal_states_the_fitted_box() {
    let table = CoefficientTable {
        kind: SurfaceKind::Canard,
        channel: "cl-residual",
        dossier_record: "a2-ames-fullscale-1999".into(),
        conventions: ConventionBlock {
            axes_id: "frd-body-v1".into(),
            moment_signs_id: "moment-signs-v1".into(),
            angles_id: "angles-v1".into(),
            reference_area_convention: "per-surface projected (geometry-conventions-v1)".into(),
            wind_reference: "not-applicable".into(),
        },
        patches: vec![RegimePatch {
            regime: "attached",
            surface: mean_surface(),
        }],
    };
    assert!(table.validate(1e-9).is_ok());
    // In-domain passes both paths identically.
    let inside = [0.1, 6.0, 0.5];
    assert_eq!(
        table.eval(inside).unwrap().to_bits(),
        table.eval_strict(inside).unwrap().to_bits()
    );
    // log Re at the box edge admitted; one float above refused with the box.
    assert!(table.eval_strict([0.1, 7.0, 0.5]).is_ok());
    let above = f64::from_bits(7.0f64.to_bits() + 1);
    let refusal = table.eval_strict([0.1, above, 0.5]).unwrap_err();
    assert_eq!(refusal.code, "query-outside-fitted-domain");
    assert!(
        refusal.message.contains("[5, 7]"),
        "box must be stated: {}",
        refusal.message
    );
    // The permissive path would have CLAMPED here — that hazard is why
    // eval_strict exists; α outside still refuses via the partition.
    assert_eq!(
        table.eval_strict([0.9, 6.0, 0.5]).unwrap_err().code,
        "alpha-outside-table"
    );
    jlog("strict-ood", "\"box_stated\":true");
}

#[test]
fn indicial_kernels_match_exact_references() {
    for k in [WAGNER_JONES, KUSSNER_2POLE] {
        k.admit().unwrap();
    }
    // Wagner: φ(0) = 1 − 0.165 − 0.335 = 0.5 EXACT; φ(∞) → 1.
    assert!((WAGNER_JONES.phi(0.0) - 0.5).abs() < 1e-15);
    assert!((WAGNER_JONES.phi(1.0e4) - 1.0).abs() < 1e-12);
    // Küssner: ψ(0) = 0 EXACT.
    assert!(KUSSNER_2POLE.phi(0.0).abs() < 1e-15);
    // Monotone non-decreasing on a grid (per-point oracle).
    let mut prev = -1.0;
    for i in 0..=400 {
        let s = 0.05 * f64::from(i);
        let v = WAGNER_JONES.phi(s);
        assert!(v >= prev - 1e-15, "Wagner must be monotone at s={s}");
        prev = v;
    }
    // State update matches the closed form: impulsive start advanced by Δs
    // equals φ(s) exactly (diagonal exponential = the true solution).
    let mut state = IndicialState::impulsive_start(&WAGNER_JONES);
    let mut s = 0.0;
    for _ in 0..48 {
        state.advance(&WAGNER_JONES, 0.25).unwrap();
        s += 0.25;
        assert!(
            (state.response() - WAGNER_JONES.phi(s)).abs() < 1e-13,
            "state update diverged from the closed form at s={s}"
        );
    }
    // Sub-step composition EXACTNESS: 1×0.4 vs 4×0.1 agree to 1e-14.
    let mut one = IndicialState::impulsive_start(&WAGNER_JONES);
    one.advance(&WAGNER_JONES, 0.4).unwrap();
    let mut four = IndicialState::impulsive_start(&WAGNER_JONES);
    for _ in 0..4 {
        four.advance(&WAGNER_JONES, 0.1).unwrap();
    }
    assert!(
        (one.response() - four.response()).abs() < 1e-14,
        "transition must be exact"
    );
    // Trim start has zero deficiency (fully developed memory).
    assert!((IndicialState::trim().response() - 1.0).abs() < 1e-15);
    // Kernel-parameter falsifier.
    let bad = IndicialKernel {
        kernel_id: "bad",
        a: [0.7, 0.5],
        b: [0.1, 1.0],
    };
    assert_eq!(bad.admit().unwrap_err().code, "kernel-params-invalid");
    jlog(
        "kernels",
        &format!("\"wagner_phi0\":{}", WAGNER_JONES.phi(0.0)),
    );
}

#[test]
fn chordwise_clock_freezes_and_refuses() {
    // ds = 2·U_conv·dt/c — the CHORDWISE clock (never the 3-D speed norm).
    let ds = reduced_time_increment(14.0, 1.981, 1.0 / 120.0).unwrap();
    assert!((ds - 2.0 * 14.0 * (1.0 / 120.0) / 1.981).abs() < 1e-15);
    // U_conv = 0: the clock FREEZES (a vertical gust advances nothing).
    assert_eq!(
        reduced_time_increment(0.0, 1.981, 1.0 / 120.0).unwrap().to_bits(),
        0.0f64.to_bits()
    );
    let mut state = IndicialState::impulsive_start(&WAGNER_JONES);
    let before = state.response();
    state.advance(&WAGNER_JONES, 0.0).unwrap();
    assert_eq!(
        state.response().to_bits(),
        before.to_bits(),
        "frozen clock must not move state"
    );
    // Reversed chordwise flow REFUSES — never |U|.
    let refusal = reduced_time_increment(-0.1, 1.981, 1.0 / 120.0).unwrap_err();
    assert_eq!(refusal.code, "indicial-flow-reversed");
    assert!(refusal.message.contains("reversed"));
    // Caps: ds at MAX_DS admitted, one float above refused.
    let mut s = IndicialState::impulsive_start(&WAGNER_JONES);
    assert!(s.advance(&WAGNER_JONES, MAX_DS).is_ok());
    let above = f64::from_bits(MAX_DS.to_bits() + 1);
    assert_eq!(
        s.advance(&WAGNER_JONES, above).unwrap_err().code,
        "reduced-time-increment-invalid"
    );
    // Negative dt and bad chord refuse with their own codes.
    assert_eq!(
        reduced_time_increment(1.0, 0.0, 0.01).unwrap_err().code,
        "chord-outside-domain"
    );
    assert_eq!(
        reduced_time_increment(1.0, 1.0, -0.01).unwrap_err().code,
        "timestep-invalid"
    );
    jlog("clock", &format!("\"ds_120hz\":{ds}"));
}

#[test]
fn wagner_trace_golden() {
    // Deterministic golden: 120 Hz Wagner build-up at the Dec-17 airspeed
    // over the 1903 chord (golden-bump protocol; measure-then-pin).
    let mut state = IndicialState::impulsive_start(&WAGNER_JONES);
    let ds = reduced_time_increment(13.86, 1.981, 1.0 / 120.0).unwrap();
    let mut payload = Vec::new();
    for _ in 0..240 {
        state.advance(&WAGNER_JONES, ds).unwrap();
        payload.extend_from_slice(&state.response().to_bits().to_le_bytes());
    }
    let digest =
        fs_blake3::hash_domain("org.frankensim.fs-airfoil.wagner-trace.v1", &payload).to_hex();
    jlog("golden", &format!("\"digest\":\"{digest}\""));
    assert_eq!(
        digest, "7897e5d780fb447d95375855996baeb326214b4a31f1e50c6146fb71b9e9e2f5",
        "Wagner trace moved — determinism regression or an intentional kernel \
         change requiring the golden-bump protocol"
    );
}
