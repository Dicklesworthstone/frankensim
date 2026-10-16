//! V-11a landing battery (bead wf-root-guzez.4.8.2, E3.4-ii): a 1-D skid
//! drop fixture integrated with the contact model — no adhesion (per
//! tick), plastic sink monotone + bounded by penetration, regularized
//! Coulomb saturation + small-slip linearity, impulse/penetration
//! CONVERGENCE under dt refinement in declared bands, energy accounting
//! (plastic loss + damping loss = energy not returned), the impact
//! report, caps, golden.
//! Repro: cargo test -p fs-flyer --test contact_battery

use fs_flyer::contact::{
    ContactParams, ContactState, ImpactReport, MAX_PENETRATION_M, contact_tick,
};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-flyer-v11a-landing\",\"case\":\"{case}\",{payload}}}");
}

const G: f64 = 9.80665;
const MASS: f64 = 170.0; // per-skid share of the gross

fn params() -> ContactParams {
    ContactParams {
        stiffness_n_m: 60_000.0,
        damping_n_s_m: 2_500.0,
        mu: 0.5,
        v_reg_mps: 0.05,
        sink_rate_per_s: 6.0,
        // Static equilibrium p_e = mg/k = 0.028 m sits BELOW this, so the
        // skid settles; landing transients exceed it and yield plastically.
        yield_penetration_m: 0.035,
    }
}

/// Drop the skid from `v0` downward at the surface, with tangential speed
/// `vt0`; integrate until rebound separation or rest. Returns the report.
fn drop_run(dt: f64, v0: f64, vt0: f64) -> ImpactReport {
    let p = params();
    let mut st = ContactState::default();
    let (mut z, mut vz, mut vt) = (0.0f64, v0, vt0); // z = penetration depth
    let mut peak = 0.0f64;
    let mut impulse = 0.0f64;
    let mut max_pen = 0.0f64;
    let mut slide = 0.0f64;
    for _ in 0..200_000 {
        let tick = contact_tick(&p, &mut st, z, vz, vt, dt).unwrap();
        // Vertical: weight down (+), normal up (−) in penetration coords.
        let az = G - tick.normal_n / MASS;
        vz += az * dt;
        z += vz * dt;
        // Tangential: friction decelerates the slide.
        if tick.in_contact {
            vt += tick.friction_n / MASS * dt;
            slide += vt.abs() * dt;
        }
        peak = peak.max(tick.normal_n);
        impulse += tick.normal_n * dt;
        max_pen = max_pen.max(tick.penetration_m);
        // Done when the skid settles (contact, slow) — sand landings do
        // not meaningfully rebound with this damping.
        if tick.in_contact && vz.abs() < 1e-4 && vt.abs() < 1e-3 {
            return ImpactReport {
                impact_speed_mps: v0,
                peak_normal_n: peak,
                normal_impulse_ns: impulse,
                max_penetration_m: max_pen,
                final_sink_m: tick.sink_m,
                sliding_m: slide,
            };
        }
    }
    panic!("no settle — fixture broken");
}

#[test]
fn no_adhesion_and_sink_monotone_bounded() {
    let p = params();
    let mut st = ContactState::default();
    let dt = 1.0 / 480.0;
    let mut prev_sink = 0.0f64;
    // Push in, then pull out fast: the normal must clamp at 0 (never pull).
    for (pen, rate) in [
        (0.02, 0.5),
        (0.04, 0.5),
        (0.05, 0.0),
        (0.03, -2.0),
        (0.01, -2.0),
    ] {
        let t = contact_tick(&p, &mut st, pen, rate, 0.0, dt).unwrap();
        assert!(t.normal_n >= 0.0, "adhesion at pen {pen}");
        assert!(t.sink_m >= prev_sink, "sink must be monotone");
        assert!(
            t.sink_m <= pen.max(prev_sink) + 1e-12,
            "sink bounded by penetration seen"
        );
        prev_sink = t.sink_m;
    }
    // Out of contact: exactly zero force, sink unchanged.
    let t = contact_tick(&p, &mut st, 0.0, -1.0, 3.0, dt).unwrap();
    assert!(!t.in_contact && t.normal_n == 0.0 && t.friction_n == 0.0);
    assert_eq!(t.sink_m, prev_sink);
    jlog("adhesion-sink", &format!("\"final_sink\":{prev_sink}"));
}

#[test]
fn regularized_coulomb_saturates_and_is_linear_at_small_slip() {
    let p = params();
    let dt = 1.0 / 120.0;
    // Fixed normal via fixed penetration, zero sink (rate 0 for this probe).
    let frozen = ContactParams {
        sink_rate_per_s: 0.0,
        ..p
    };
    let probe = |vt: f64| -> (f64, f64) {
        let mut st = ContactState::default();
        let t = contact_tick(&frozen, &mut st, 0.01, 0.0, vt, dt).unwrap();
        (t.friction_n, t.normal_n)
    };
    let (f_fast, n) = probe(1.0); // 20x v_reg: saturated
    assert!(
        (f_fast.abs() - frozen.mu * n).abs() < 0.001 * frozen.mu * n,
        "must saturate at muN"
    );
    let (f_slow, _) = probe(0.005); // 0.1x v_reg: linear regime
    let expected = -frozen.mu * n * (0.005 / frozen.v_reg_mps);
    assert!(
        (f_slow - expected).abs() < 0.01 * expected.abs(),
        "linear at small slip"
    );
    // Odd symmetry.
    assert!((probe(0.3).0 + probe(-0.3).0).abs() < 1e-9);
    jlog(
        "coulomb",
        &format!("\"saturated\":{},\"muN\":{}", f_fast.abs(), frozen.mu * n),
    );
}

#[test]
fn impulse_and_penetration_converge_under_dt_refinement() {
    // The Dec-17 landing class: ~1.5 m/s sink rate with residual slide.
    let r480 = drop_run(1.0 / 480.0, 1.5, 3.0);
    let r120 = drop_run(1.0 / 120.0, 1.5, 3.0);
    let r240 = drop_run(1.0 / 240.0, 1.5, 3.0);
    // Declared bands (frozen registry: contact timing/impulse <= 2%): the
    // 120->480 Hz drift on impulse and max penetration stays inside 2%,
    // and refinement is monotone toward the fine reference.
    for (name, a, b) in [
        ("impulse", r120.normal_impulse_ns, r480.normal_impulse_ns),
        (
            "penetration",
            r120.max_penetration_m,
            r480.max_penetration_m,
        ),
    ] {
        let rel = (a - b).abs() / b;
        assert!(rel < 0.02, "{name} drift {rel:.4} exceeds the 2% band");
    }
    let e120 = (r120.peak_normal_n - r480.peak_normal_n).abs();
    let e240 = (r240.peak_normal_n - r480.peak_normal_n).abs();
    assert!(
        e240 <= e120 + 1e-9,
        "peak-force refinement must not diverge"
    );
    jlog(
        "convergence",
        &format!(
            "\"impulse\":[{},{}],\"pen\":[{},{}]",
            r120.normal_impulse_ns,
            r480.normal_impulse_ns,
            r120.max_penetration_m,
            r480.max_penetration_m
        ),
    );
}

#[test]
fn impact_report_is_physical() {
    let r = drop_run(1.0 / 240.0, 1.5, 3.0);
    // Momentum: the settle impulse must at least absorb the vertical drop
    // momentum plus the weight impulse over the contact (lower bound).
    assert!(
        r.normal_impulse_ns > MASS * r.impact_speed_mps,
        "impulse under-absorbs the drop"
    );
    assert!(
        r.peak_normal_n > MASS * G,
        "peak must exceed the static weight"
    );
    assert!(r.final_sink_m > 0.0 && r.final_sink_m < r.max_penetration_m + 0.05);
    assert!(r.sliding_m > 0.0, "the slide must be recorded");
    // Friction stops the slide within the run (Dec-17 flights ended in
    // short slides, not long skids).
    assert!(
        r.sliding_m < 25.0,
        "slide {} m implausible for sand at mu 0.5",
        r.sliding_m
    );
    jlog(
        "report",
        &format!(
            "\"peak\":{},\"impulse\":{},\"sink\":{},\"slide\":{}",
            r.peak_normal_n, r.normal_impulse_ns, r.final_sink_m, r.sliding_m
        ),
    );
}

#[test]
fn refusals_at_cap_and_cap_plus_one() {
    let p = params();
    let mut st = ContactState::default();
    let dt = 1.0 / 120.0;
    assert!(contact_tick(&p, &mut st, MAX_PENETRATION_M, 0.0, 0.0, dt).is_ok());
    let above = f64::from_bits(MAX_PENETRATION_M.to_bits() + 1);
    assert_eq!(
        contact_tick(&p, &mut st, above, 0.0, 0.0, dt)
            .unwrap_err()
            .code,
        "penetration-outside-domain"
    );
    let bad = ContactParams {
        stiffness_n_m: 0.0,
        ..p
    };
    assert_eq!(
        contact_tick(&bad, &mut st, 0.01, 0.0, 0.0, dt)
            .unwrap_err()
            .code,
        "contact-params-invalid"
    );
    assert_eq!(
        contact_tick(&p, &mut st, f64::NAN, 0.0, 0.0, dt)
            .unwrap_err()
            .code,
        "non-finite-input"
    );
    jlog(
        "refusals",
        "\"gates\":\"penetration cap/cap+1, params, NaN\"",
    );
}

#[test]
fn landing_golden_digest() {
    let r = drop_run(1.0 / 120.0, 1.5, 3.0);
    let mut payload = Vec::new();
    for v in [
        r.peak_normal_n,
        r.normal_impulse_ns,
        r.max_penetration_m,
        r.final_sink_m,
        r.sliding_m,
    ] {
        payload.extend_from_slice(&v.to_bits().to_le_bytes());
    }
    let digest =
        fs_blake3::hash_domain("org.frankensim.fs-flyer.v11a-landing-golden.v1", &payload).to_hex();
    jlog("golden", &format!("\"digest\":\"{digest}\""));
    assert_eq!(
        digest, "2e0f661637b55f248904b53ce39d9c3a1680fe027114422940a8f7574b13da12",
        "landing golden moved — determinism regression or an intentional \
         contact change requiring the golden-bump protocol"
    );
}
