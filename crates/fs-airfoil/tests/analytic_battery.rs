//! E4.0a analytic-limit battery (bead wf-root-guzez.5.1.1). Per-item
//! oracles against CLASSICAL EXACT results (never totals-only), refusals
//! at cap AND cap+1, falsifier-style negatives, a pinned polar golden,
//! and JSONL receipts. Repro: cargo test -p fs-airfoil --test analytic_battery

use fs_airfoil::{
    FLAT_PLATE_CD90, MAX_ABS_ALPHA_RAD, MAX_CAMBER_RATIO, MAX_LOG10_RE, MIN_LOG10_RE, NormalAxial,
    SectionGeometry, admit_query, body_to_wind, flat_plate_separated, thin_airfoil,
    transfer_moment, wind_to_body,
};
use std::f64::consts::PI;

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-airfoil-analytic\",\"case\":\"{case}\",{payload}}}");
}

const RE: f64 = 6.0; // log10(Re) = 1e6, mid-envelope

#[test]
fn thin_airfoil_exact_classical_results() {
    let f = 0.05; // the 1903 camber
    // Zero-lift angle: cl(α₀) = 0 at α₀ = −2f — EXACT.
    let at_a0 = thin_airfoil(-2.0 * f, f, RE).unwrap();
    assert!(at_a0.cl.abs() < 1e-15, "cl(α₀) = {} ≠ 0", at_a0.cl);
    // cl(0) = 4πf — EXACT.
    let at_zero = thin_airfoil(0.0, f, RE).unwrap();
    assert!((at_zero.cl - 4.0 * PI * f).abs() < 1e-15);
    // Lift slope dcl/dα = 2π — measured by central difference.
    let h = 1e-6;
    let slope =
        (thin_airfoil(h, f, RE).unwrap().cl - thin_airfoil(-h, f, RE).unwrap().cl) / (2.0 * h);
    assert!((slope - 2.0 * PI).abs() < 1e-6, "slope {slope} vs 2π");
    // Quarter-chord moment cm = −πf, independent of α (aerodynamic centre).
    for alpha in [-0.1, 0.0, 0.12] {
        let cm = thin_airfoil(alpha, f, RE).unwrap().cm_quarter;
        assert!((cm - (-PI * f)).abs() < 1e-15, "cm_c/4 {cm} at α={alpha}");
    }
    // Symmetric section: cl odd in α.
    let s = thin_airfoil(0.07, 0.0, RE).unwrap().cl;
    let n = thin_airfoil(-0.07, 0.0, RE).unwrap().cl;
    assert_eq!(
        s.to_bits(),
        (-n).to_bits(),
        "symmetric section must be exactly odd"
    );
    jlog(
        "thin-airfoil",
        &format!("\"slope\":{slope},\"cm_c4\":{}", at_zero.cm_quarter),
    );
}

#[test]
fn flat_plate_separated_shape() {
    // |cn| at ±90° is exactly CD90; cp at mid-chord.
    let up = flat_plate_separated(PI / 2.0, RE).unwrap();
    assert!((up.cn - FLAT_PLATE_CD90).abs() < 1e-12);
    let x_cp = 0.25 - up.cm_ref / up.cn;
    assert!(
        (x_cp - 0.5).abs() < 1e-12,
        "cp at 90° must be mid-chord, got {x_cp}"
    );
    // Odd symmetry in α (per-item, both branches).
    let a = flat_plate_separated(0.6, RE).unwrap();
    let b = flat_plate_separated(-0.6, RE).unwrap();
    assert_eq!(a.cn.to_bits(), (-b.cn).to_bits());
    assert_eq!(a.cm_ref.to_bits(), (-b.cm_ref).to_bits());
    // Small-α linearization: cn ≈ CD90·α (the separated branch is BELOW the
    // attached 2π slope — the baselines must not be conflated).
    let small = flat_plate_separated(0.01, RE).unwrap();
    assert!((small.cn - FLAT_PLATE_CD90 * 0.01).abs() < 1e-5);
    // No leading-edge suction when separated.
    assert!(a.ca.abs() == 0.0, "no suction when separated");
    jlog(
        "flat-plate",
        &format!("\"cn90\":{},\"x_cp90\":{x_cp}", up.cn),
    );
}

#[test]
fn decomposition_round_trips_and_moment_transfer_identity() {
    // wind→body→wind is the identity (exact rotation pair).
    for (cl, cd, alpha) in [(1.2, 0.08, 0.15), (-0.7, 0.3, -0.9), (0.0, 1.98, 1.4)] {
        let (cn, ca) = wind_to_body(cl, cd, alpha);
        let (cl2, cd2) = body_to_wind(cn, ca, alpha);
        assert!(
            (cl - cl2).abs() < 1e-14 && (cd - cd2).abs() < 1e-14,
            "round-trip at α={alpha}"
        );
    }
    // Moment transfer A→B→A returns exactly; A→B matches direct C-of-P form.
    let (cm_a, cn) = (-0.16, 1.1);
    let cm_b = transfer_moment(cm_a, cn, 0.25, 0.5);
    assert_eq!(
        transfer_moment(cm_b, cn, 0.5, 0.25).to_bits(),
        cm_a.to_bits(),
        "transfer must invert exactly"
    );
    // A moment that is pure cp offset: cm about cp is zero.
    let x_cp = 0.25 - cm_a / cn;
    let at_cp = transfer_moment(cm_a, cn, 0.25, x_cp);
    assert!(
        at_cp.abs() < 1e-15,
        "cm about the centre of pressure must vanish"
    );
    jlog("decomposition", &format!("\"cm_b\":{cm_b},\"x_cp\":{x_cp}"));
}

#[test]
fn refusals_at_cap_and_cap_plus_one() {
    // α: π admitted, next float above refused (cap AND cap+1 law).
    assert!(admit_query(MAX_ABS_ALPHA_RAD, RE).is_ok());
    let above = f64::from_bits(MAX_ABS_ALPHA_RAD.to_bits() + 1);
    assert_eq!(
        admit_query(above, RE).unwrap_err().code,
        "alpha-outside-domain"
    );
    // log10 Re: both edges admitted, just outside refused.
    assert!(admit_query(0.0, MIN_LOG10_RE).is_ok());
    assert!(admit_query(0.0, MAX_LOG10_RE).is_ok());
    let re_hi = f64::from_bits(MAX_LOG10_RE.to_bits() + 1);
    let re_lo = f64::from_bits(MIN_LOG10_RE.to_bits() - 1);
    assert_eq!(
        admit_query(0.0, re_hi).unwrap_err().code,
        "reynolds-outside-domain"
    );
    assert_eq!(
        admit_query(0.0, re_lo).unwrap_err().code,
        "reynolds-outside-domain"
    );
    // The refusal STATES the admitted domain (applicability-domain law).
    let msg = admit_query(0.0, 99.0).unwrap_err().message;
    assert!(msg.contains("[4, 8]"), "domain must be stated, got {msg}");
    // Camber cap and cap+1.
    assert!(thin_airfoil(0.0, MAX_CAMBER_RATIO, RE).is_ok());
    let camber_above = f64::from_bits(MAX_CAMBER_RATIO.to_bits() + 1);
    assert_eq!(
        thin_airfoil(0.0, camber_above, RE).unwrap_err().code,
        "camber-outside-domain"
    );
    // Non-finite is its own code, checked before domains.
    assert_eq!(
        admit_query(f64::NAN, RE).unwrap_err().code,
        "non-finite-input"
    );
    jlog("refusals", "\"caps\":\"alpha,re,camber at cap and cap+1\"");
}

#[test]
fn geometry_admission_and_provenance_gate() {
    let good = SectionGeometry {
        chord_m: 1.981,
        camber_ratio: 0.05,
        dossier_record: "a1-wright-1901-tunnel".into(),
        digitization_class: "drawings-pm-2-5-mm".into(),
    };
    assert!(good.admit().is_ok());
    // Falsifier: provenance stripped → refusal (the provenance gate is live).
    let mut anon = good.clone();
    anon.dossier_record.clear();
    assert_eq!(anon.admit().unwrap_err().code, "provenance-missing");
    let mut bad_chord = good.clone();
    bad_chord.chord_m = 0.0;
    assert_eq!(bad_chord.admit().unwrap_err().code, "chord-outside-domain");
    jlog("geometry", "\"provenance_gate\":\"live\"");
}

#[test]
fn polar_sweep_golden_digest() {
    // Deterministic golden over a 33-point α sweep of both baselines.
    // Golden-bump protocol applies to any change (measure-then-pin).
    let mut payload = Vec::new();
    for i in 0..33 {
        let alpha = -MAX_ABS_ALPHA_RAD + (2.0 * MAX_ABS_ALPHA_RAD) * (f64::from(i) / 32.0);
        let alpha = alpha.clamp(-MAX_ABS_ALPHA_RAD, MAX_ABS_ALPHA_RAD);
        let NormalAxial { cn, ca, cm_ref, .. } = flat_plate_separated(alpha, RE).unwrap();
        for v in [cn, ca, cm_ref] {
            payload.extend_from_slice(&v.to_bits().to_le_bytes());
        }
        if alpha.abs() <= 0.3 {
            let c = thin_airfoil(alpha, 0.05, RE).unwrap();
            for v in [c.cl, c.cd, c.cm_quarter] {
                payload.extend_from_slice(&v.to_bits().to_le_bytes());
            }
        }
    }
    let digest =
        fs_blake3::hash_domain("org.frankensim.fs-airfoil.analytic-polar.v1", &payload).to_hex();
    jlog("golden", &format!("\"digest\":\"{digest}\""));
    assert_eq!(
        digest, "51ad6b2e05e67a4c316aa5eee02c2688a90c3ccdf4d259384bf1b11da117f719",
        "analytic polar digest moved — determinism regression or an intentional \
         baseline change requiring the golden-bump protocol"
    );
}
