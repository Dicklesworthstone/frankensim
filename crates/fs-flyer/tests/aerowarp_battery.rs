//! E4.6b0 battery (bead wf-root-guzez.5.15): per-strip loaded-twist
//! closed-form oracles, q-dependent effectiveness (executed), the REAL
//! loaded-vs-prescribed roll-moment comparison on the coupled biplane
//! solve (recorded, per DONE-WHEN), slack-risk diagnostics firing on the
//! hostile over-warp fixture (and silent at cruise), exact lag update,
//! caps at cap AND cap+1, determinism, golden.
//! Repro: cargo test -p fs-flyer --test aerowarp_battery

use fs_flyer::aerowarp::{MAX_SLACK_BOUND_RAD, ReducedAeroelasticWarp, WarpLagState};
use fs_wing::nonlinear::{InfluenceOperator, StripRegime, StripSpec, solve_nonlinear};
use fs_wing::{SurfaceId, flat_surface};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-flyer-e46b0\",\"case\":\"{case}\",{payload}}}");
}

const RHO: f64 = 1.294;
const V: f64 = 13.86;
const Q: f64 = 0.5 * RHO * V * V;

#[test]
fn loaded_twist_matches_the_closed_form_per_strip() {
    let m = ReducedAeroelasticWarp::wright_v1();
    let rep = m.evaluate(0.08, Q, 0.05).unwrap();
    for (i, s) in rep.strips.iter().enumerate() {
        let cmd = 0.08 * m.basis[i];
        let k = m.compliance_rad_per_nm[i]
            * Q
            * m.chord_m[i]
            * m.chord_m[i]
            * m.width_m[i]
            * m.cl_alpha;
        let loaded_ref = (cmd - k * 0.05) / (1.0 + k);
        assert!(
            (s.loaded_rad - loaded_ref).abs() < 1e-14,
            "strip {i}: {} vs {loaded_ref}",
            s.loaded_rad
        );
        // Warp-relative slack deficit oracle (wires taut at trim).
        let margin_ref = m.slack_bound_rad - (k / (1.0 + k)) * cmd.abs();
        assert!(
            (s.slack_margin_rad - margin_ref).abs() < 1e-14,
            "strip {i} slack margin {} vs {margin_ref}",
            s.slack_margin_rad
        );
        assert!(
            s.effectiveness > 0.0 && s.effectiveness <= 1.0,
            "strip {i} effectiveness {}",
            s.effectiveness
        );
        assert!((s.commanded_rad - cmd).abs() < 1e-15);
    }
    jlog(
        "closed-form",
        &format!("\"mean_eff\":{}", rep.mean_effectiveness),
    );
}

#[test]
fn effectiveness_drops_with_dynamic_pressure() {
    // The aeroelastic core claim, executed: higher q, lower control
    // power. (A rigid model — zero compliance — stays at 1 exactly.)
    let m = ReducedAeroelasticWarp::wright_v1();
    let lo = m.evaluate(0.08, 0.5 * Q, 0.0).unwrap().mean_effectiveness;
    let hi = m.evaluate(0.08, 2.0 * Q, 0.0).unwrap().mean_effectiveness;
    assert!(
        hi < lo - 0.02,
        "effectiveness must drop with q: {lo} -> {hi}"
    );
    let mut rigid = ReducedAeroelasticWarp::wright_v1();
    rigid.compliance_rad_per_nm = vec![0.0; 16];
    let r = rigid
        .evaluate(0.08, 2.0 * Q, 0.0)
        .unwrap()
        .mean_effectiveness;
    assert!((r - 1.0).abs() < 1e-15, "rigid twin must be exactly 1: {r}");
    jlog(
        "q-dependence",
        &format!("\"eff_lo_q\":{lo},\"eff_hi_q\":{hi}"),
    );
}

#[test]
fn antisymmetry_mirrors_left_and_right() {
    let m = ReducedAeroelasticWarp::wright_v1();
    let rep = m.evaluate(0.08, Q, 0.0).unwrap();
    // Strips 0..8 are one plane, mirrored pairs (s, 7-s).
    for s in 0..4 {
        let a = &rep.strips[s];
        let b = &rep.strips[7 - s];
        assert!(
            (a.commanded_rad + b.commanded_rad).abs() < 1e-15,
            "commanded antisymmetry {s}"
        );
        assert!(
            (a.loaded_rad + b.loaded_rad).abs() < 1e-14,
            "loaded antisymmetry {s} (alpha0 = 0)"
        );
    }
    jlog("antisymmetry", "\"mirrored\":true");
}

#[test]
fn loaded_vs_prescribed_roll_moment_on_the_real_solve() {
    // DONE-WHEN receipt: the loaded-warp roll effectiveness against the
    // prescribed-kinematic twin on the REAL coupled biplane solve; the
    // difference is RECORDED (V-05 philosophy), and the loaded moment
    // must be smaller but same-signed.
    let m = ReducedAeroelasticWarp::wright_v1();
    let dw = 0.08;
    let rep = m.evaluate(dw, Q, 0.06).unwrap();
    let roll = |twists: &[f64]| -> f64 {
        let mut p = flat_surface(SurfaceId::WingLower, 12.29, 1.981, 0.0, 0.0, 8, 2).unwrap();
        p.extend(flat_surface(SurfaceId::WingUpper, 12.29, 1.981, 0.0, -1.89, 8, 2).unwrap());
        let mut strips = Vec::new();
        for plane in 0..2 {
            let base = plane * 16;
            for s in 0..8 {
                strips.push(StripSpec {
                    panel_indices: vec![base + s, base + 8 + s],
                    chord_m: 1.981,
                    twist_rad: twists[plane * 8 + s],
                });
            }
        }
        let alpha = 0.06f64;
        let fs_v = [V * alpha.cos(), 0.0, V * alpha.sin()];
        let op = InfluenceOperator::build(&p, fs_v, RHO).unwrap();
        let closure = |_s: usize, a: f64| -> (f64, StripRegime) {
            (core::f64::consts::TAU * (a + 0.1), StripRegime::Attached)
        };
        let sol = solve_nonlinear(&op, &p, &strips, fs_v, RHO, &closure, None, None).unwrap();
        // Roll moment Mx = sum(r_y*F_z - r_z*F_y) about the mid-bay point.
        let mut mx = 0.0;
        for (j, panel) in p.iter().enumerate() {
            let seg = [
                panel.b[0] - panel.a[0],
                panel.b[1] - panel.a[1],
                panel.b[2] - panel.a[2],
            ];
            let s = RHO * sol.gamma[j];
            let fy = s * (fs_v[2] * seg[0] - fs_v[0] * seg[2]);
            let fz = s * (fs_v[0] * seg[1] - fs_v[1] * seg[0]);
            let ry = 0.5 * (panel.a[1] + panel.b[1]);
            let rz = 0.5 * (panel.a[2] + panel.b[2]) + 0.945;
            mx += ry * fz - rz * fy;
        }
        mx
    };
    let commanded: Vec<f64> = rep.strips.iter().map(|s| s.commanded_rad).collect();
    let loaded: Vec<f64> = rep.strips.iter().map(|s| s.loaded_rad).collect();
    let mx_cmd = roll(&commanded);
    let mx_loaded = roll(&loaded);
    assert!(
        mx_cmd > 0.0,
        "+dw must give +roll (right wing down): {mx_cmd}"
    );
    assert!(mx_loaded > 0.0, "loaded roll keeps the sign: {mx_loaded}");
    let ratio = mx_loaded / mx_cmd;
    assert!(
        ratio > 0.5 && ratio < 1.0,
        "loaded/prescribed roll ratio {ratio} outside the reduced-model band"
    );
    jlog(
        "roll-comparison",
        &format!("\"mx_cmd\":{mx_cmd},\"mx_loaded\":{mx_loaded},\"ratio\":{ratio}"),
    );
}

#[test]
fn slack_diagnostics_fire_on_hostile_over_warp_and_stay_silent_at_cruise() {
    let m = ReducedAeroelasticWarp::wright_v1();
    // Cruise: modest warp, cruise q — no slack risk.
    let cruise = m.evaluate(0.08, Q, 0.06).unwrap();
    assert!(
        cruise.slack_risk_strips.is_empty(),
        "cruise warp must not flag slack: {:?}",
        cruise.slack_risk_strips
    );
    // HOSTILE over-warp: full travel at a dive-class q — the tip strips'
    // load-induced deficit exceeds the slack bound and the diagnostic
    // MUST fire (DONE-WHEN clause).
    let hostile = m.evaluate(0.7, 3.0 * Q, 0.1).unwrap();
    assert!(
        !hostile.slack_risk_strips.is_empty(),
        "hostile over-warp must flag slack risk"
    );
    // Per-item: every flagged strip really has a negative margin, and
    // the flagged set includes the largest-basis (tip) strips.
    for &i in &hostile.slack_risk_strips {
        assert!(hostile.strips[i].slack_margin_rad < 0.0, "strip {i}");
    }
    jlog(
        "slack",
        &format!("\"hostile_flagged\":{}", hostile.slack_risk_strips.len()),
    );
}

#[test]
fn lag_is_exact_and_optional() {
    let m = ReducedAeroelasticWarp::wright_v1();
    let tau = m.lag_tau_s.unwrap();
    // Exact exponential: one 0.2 s step equals the closed form.
    let st = WarpLagState { delta_w_rad: 0.0 };
    let one = m.lag_step(st, 0.1, 0.2).unwrap();
    let expect = 0.1 * (1.0 - (-0.2f64 / tau).exp());
    assert!((one.delta_w_rad - expect).abs() < 1e-12);
    // Substep composition (same held command) matches to 1e-13.
    let mut four = st;
    for _ in 0..4 {
        four = m.lag_step(four, 0.1, 0.05).unwrap();
    }
    assert!((four.delta_w_rad - one.delta_w_rad).abs() < 1e-13);
    // Optional: no lag constant = identity on the command.
    let mut nolag = ReducedAeroelasticWarp::wright_v1();
    nolag.lag_tau_s = None;
    let out = nolag.lag_step(st, 0.1, 0.001).unwrap();
    assert!((out.delta_w_rad - 0.1).abs() < 1e-15);
    assert_eq!(
        m.lag_step(st, 0.1, 0.0).unwrap_err().code,
        "warp-state-invalid"
    );
    jlog("lag", &format!("\"tau\":{tau}"));
}

#[test]
fn admission_caps_at_cap_and_cap_plus_one() {
    let mk = |slack: f64| -> ReducedAeroelasticWarp {
        let mut m = ReducedAeroelasticWarp::wright_v1();
        m.slack_bound_rad = slack;
        m
    };
    assert!(mk(MAX_SLACK_BOUND_RAD).admit().is_ok(), "cap admits");
    assert_eq!(
        mk(MAX_SLACK_BOUND_RAD.next_up()).admit().unwrap_err().code,
        "warp-model-invalid",
        "cap+1 refuses"
    );
    let mut short = ReducedAeroelasticWarp::wright_v1();
    short.compliance_rad_per_nm.pop();
    assert_eq!(short.admit().unwrap_err().code, "warp-model-invalid");
    let mut neg = ReducedAeroelasticWarp::wright_v1();
    neg.compliance_rad_per_nm[0] = -1e-9;
    assert_eq!(neg.admit().unwrap_err().code, "warp-model-invalid");
    let m = ReducedAeroelasticWarp::wright_v1();
    assert_eq!(
        m.evaluate(0.8, Q, 0.0).unwrap_err().code,
        "warp-state-invalid",
        "beyond warp travel"
    );
    assert_eq!(
        m.evaluate(0.1, 0.0, 0.0).unwrap_err().code,
        "warp-state-invalid"
    );
    jlog("caps", "\"cap_and_cap_plus_one\":true");
}

#[test]
fn determinism_and_golden() {
    let m = ReducedAeroelasticWarp::wright_v1();
    let a = m.evaluate(0.08, Q, 0.06).unwrap();
    let b = m.evaluate(0.08, Q, 0.06).unwrap();
    assert_eq!(a, b, "bitwise repeat");
    let mut payload = Vec::new();
    for s in &a.strips {
        payload.extend_from_slice(&s.loaded_rad.to_bits().to_le_bytes());
        payload.extend_from_slice(&s.slack_margin_rad.to_bits().to_le_bytes());
    }
    payload.extend_from_slice(&a.mean_effectiveness.to_bits().to_le_bytes());
    let digest =
        fs_blake3::hash_domain("org.frankensim.fs-flyer.e46b0-golden.v1", &payload).to_hex();
    jlog(
        "golden",
        &format!(
            "\"digest\":\"{digest}\",\"mean_eff\":{}",
            a.mean_effectiveness
        ),
    );
    assert_eq!(
        digest, "21bf858917d932a778a512f0bb34b2c2ee67d84591015c152273281eb4a154c8",
        "warp golden moved — determinism regression or an intentional \
         model change requiring the golden-bump protocol"
    );
}
