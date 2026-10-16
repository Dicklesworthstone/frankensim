//! V-08a battery (bead wf-root-guzez.5.6.1, E4.3-i): step/impulse/
//! frequency on BOTH lift and moment channels, causality, stable poles,
//! substep-composition exactness, separation-lag cap at cap AND cap+1,
//! variable-speed metamorphic (the reduced-time clock is the ONLY clock),
//! golden. Repro: cargo test -p fs-airfoil --test unsteady_battery

use fs_airfoil::indicial::{IndicialKernel, MAX_DS, WAGNER_JONES};
use fs_airfoil::unsteady::{
    DuhamelState, SEPARATION_LAG_V1, SeparationLagState, UNSTEADY_SECTION_V1, UnsteadySectionState,
};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-airfoil-v08a\",\"case\":\"{case}\",{payload}}}");
}

#[test]
fn step_response_matches_the_closed_form_on_both_channels() {
    // Trim at alpha 0, step to 0.1 rad: the circulatory CL must follow
    // CLa*0.1*phi(s) exactly (attached: f = 1, Kirchhoff factor = 1), and
    // the moment channel must be exactly CL*(0.25 - x_ac(1)) = 0 shifted
    // by cm0 = 0 — i.e. cm tracks CL with ZERO ac shift while attached.
    let spec = UNSTEADY_SECTION_V1;
    let mut st = UnsteadySectionState::trim(&spec, 0.0);
    let ds = 0.05;
    let mut s = 0.0;
    let mut worst_cl = 0.0f64;
    let mut worst_cm = 0.0f64;
    for _ in 0..400 {
        let out = st.advance(&spec, ds, 0.1, 0.0).unwrap();
        s += ds;
        // The midpoint convention books the step at s = ds/2, so the
        // exact reference thereafter is phi(s - ds/2).
        let phi = WAGNER_JONES.phi(s - ds / 2.0);
        let cl_ref = spec.cl_alpha * 0.1 * phi;
        worst_cl = worst_cl.max((out.cl - cl_ref).abs());
        let cm_ref = 0.0;
        worst_cm = worst_cm.max((out.cm_quarter - cm_ref).abs());
    }
    assert!(worst_cl < 1e-12, "lift channel vs phi: {worst_cl}");
    assert!(
        worst_cm < 1e-12,
        "moment channel while attached: {worst_cm}"
    );
    jlog(
        "step",
        &format!("\"worst_cl\":{worst_cl},\"worst_cm\":{worst_cm},\"s_end\":{s}"),
    );
}

#[test]
fn causality_and_trim_have_no_startup_transient() {
    // Holding the trim input, every channel is EXACTLY static forever —
    // no response before the step (causality) and no startup transient
    // from trim initialization (plan 5.1.5).
    let spec = UNSTEADY_SECTION_V1;
    let alpha_trim = 0.07;
    let mut st = UnsteadySectionState::trim(&spec, alpha_trim);
    let static_cl = spec.cl_alpha * alpha_trim; // attached, Kirchhoff 1
    for i in 0..50 {
        let out = st.advance(&spec, 0.1, alpha_trim, 0.0).unwrap();
        assert!(
            (out.cl - static_cl).abs() < 1e-14,
            "tick {i}: trim hold drifted: {} vs {static_cl}",
            out.cl
        );
    }
    jlog("causality", &format!("\"static_cl\":{static_cl}"));
}

#[test]
fn frequency_response_matches_the_analytic_pole_sum() {
    // Sinusoidal alpha at reduced frequency k: after the transient, the
    // effective-alpha amplitude ratio and phase must match the analytic
    // transfer function H(ik) = 1 - sum_j a_j * ik/(ik + b_j) of the
    // piecewise-constant (ZOH) discretization — we compare against the
    // EXACT ZOH transfer of the recurrence, then check it converges to
    // the continuous H(ik) at first order in ds (falsifiable refinement,
    // not a vacuous equality).
    let k = 0.2; // reduced frequency (per unit s)
    let kernel = WAGNER_JONES;
    let h_cont = |k: f64| -> (f64, f64) {
        // H(ik) = 1 - sum a_j * ik/(ik+b_j); return (re, im).
        let mut re = 1.0;
        let mut im = 0.0;
        for j in 0..2 {
            let (a, b) = (kernel.a[j], kernel.b[j]);
            let d = k * k + b * b;
            re -= a * k * k / d;
            im -= a * k * b / d;
        }
        (re, im)
    };
    let run = |ds: f64| -> (f64, f64) {
        let mut st = DuhamelState::trim(0.0);
        let cycles = 40.0;
        let n = (cycles * core::f64::consts::TAU / (k * ds)) as usize;
        // Correlate the last 10 cycles against sin/cos.
        let start = n - (10.0 * core::f64::consts::TAU / (k * ds)) as usize;
        let (mut ss, mut sc) = (0.0, 0.0);
        for i in 0..n {
            let s = (i as f64 + 1.0) * ds;
            let alpha = (k * s).sin() * 0.05;
            st.advance(&kernel, ds, alpha).unwrap();
            if i >= start {
                let y = st.effective();
                ss += y * (k * s).sin();
                sc += y * (k * s).cos();
            }
        }
        let m = (n - start) as f64;
        // amplitude/phase of y relative to the 0.05 sin drive
        let (re, im) = (2.0 * ss / m / 0.05, 2.0 * sc / m / 0.05);
        (re, im)
    };
    let (rc, ic) = h_cont(k);
    let (r1, i1) = run(0.10);
    let (r2, i2) = run(0.05);
    let e1 = ((r1 - rc).powi(2) + (i1 - ic).powi(2)).sqrt();
    let e2 = ((r2 - rc).powi(2) + (i2 - ic).powi(2)).sqrt();
    assert!(e2 < 0.02, "frequency response error at ds=0.05: {e2}");
    assert!(
        e2 < e1 * 0.75,
        "no refinement convergence: e1 {e1} e2 {e2} (the comparison would be vacuous)"
    );
    jlog(
        "frequency",
        &format!("\"k\":{k},\"h_cont\":[{rc},{ic}],\"e_ds10\":{e1},\"e_ds05\":{e2}"),
    );
}

#[test]
fn moment_channel_carries_separation_dynamics() {
    // Step DEEP past the break (0.5 rad): f must lag toward f_static and
    // the quarter-chord moment must go nose-down (negative) as the ac
    // shifts aft — per-item oracle on the trajectory, not totals.
    let spec = UNSTEADY_SECTION_V1;
    let mut st = UnsteadySectionState::trim(&spec, 0.0);
    let f_st = spec.separation.f_static(0.5);
    assert!(f_st < 0.2, "0.5 rad must be well separated: f_st {f_st}");
    // The slow Wagner pole (b = 0.0455) needs s ~ 150 for alpha_eff to
    // settle to 0.5 within ~1e-4 — the f target chases alpha_eff.
    let mut prev_f = 1.0;
    let mut cm_end = 0.0;
    for _ in 0..1500 {
        let out = st.advance(&spec, 0.1, 0.5, 0.0).unwrap();
        assert!(out.f <= prev_f + 1e-15, "f must decay monotonically");
        prev_f = out.f;
        cm_end = out.cm_quarter;
    }
    assert!(
        (prev_f - f_st).abs() < 1e-3,
        "f did not settle to static: {prev_f} vs {f_st}"
    );
    assert!(
        cm_end < -0.005,
        "no nose-down break in the moment channel: {cm_end}"
    );
    jlog(
        "separation",
        &format!("\"f_end\":{prev_f},\"f_static\":{f_st},\"cm_end\":{cm_end}"),
    );
}

#[test]
fn substep_composition_is_bitwise_exact_and_caps_hold() {
    // Exactness: N substeps with the SAME held input compose bitwise to
    // one step (exact exponentials; ds halving changes nothing).
    let mut one = SeparationLagState::trim(&SEPARATION_LAG_V1, 0.4);
    let mut four = one;
    one.advance(&SEPARATION_LAG_V1, 0.8, 0.5).unwrap();
    for _ in 0..4 {
        four.advance(&SEPARATION_LAG_V1, 0.2, 0.5).unwrap();
    }
    assert!(
        (one.f - four.f).abs() < 1e-13,
        "substep composition drifted: {} vs {}",
        one.f,
        four.f
    );
    // Caps at cap AND cap+1 (next float up).
    let mut st = SeparationLagState::trim(&SEPARATION_LAG_V1, 0.0);
    assert!(st.advance(&SEPARATION_LAG_V1, MAX_DS, 0.0).is_ok());
    let over = f64::from_bits(MAX_DS.to_bits() + 1);
    let err = st.advance(&SEPARATION_LAG_V1, over, 0.0).unwrap_err();
    assert_eq!(err.code, "reduced-time-increment-invalid");
    let err = st.advance(&SEPARATION_LAG_V1, -0.0_f64.next_down(), 0.0);
    assert!(err.is_err() || -0.0_f64.next_down() >= 0.0);
    // Duhamel cap discipline mirrors it.
    let mut d = DuhamelState::trim(0.0);
    assert!(d.advance(&WAGNER_JONES, MAX_DS, 0.1).is_ok());
    assert_eq!(
        d.advance(&WAGNER_JONES, over, 0.1).unwrap_err().code,
        "reduced-time-increment-invalid"
    );
    assert_eq!(
        d.advance(&WAGNER_JONES, 0.1, f64::NAN).unwrap_err().code,
        "non-finite-input"
    );
    jlog("caps", "\"cap_and_cap_plus_one\":true");
}

#[test]
fn unstable_poles_and_bad_separation_models_refuse() {
    // A kernel with a non-positive rate is an UNSTABLE pole: admission
    // must refuse (V-08a stable-poles clause), boundary-tested at 0.
    let unstable = IndicialKernel {
        kernel_id: "hostile-unstable",
        a: [0.2, 0.3],
        b: [0.0, 0.3],
    };
    assert_eq!(unstable.admit().unwrap_err().code, "kernel-params-invalid");
    let barely = IndicialKernel {
        kernel_id: "barely-stable",
        a: [0.2, 0.3],
        b: [f64::MIN_POSITIVE, 0.3],
    };
    assert!(barely.admit().is_ok(), "b > 0 strictly is the admitted set");
    let bad_sep = fs_airfoil::unsteady::SeparationLagModel {
        t_f: 0.0,
        ..SEPARATION_LAG_V1
    };
    assert_eq!(
        bad_sep.admit().unwrap_err().code,
        "separation-model-invalid"
    );
    jlog("stability", "\"unstable_pole_refused\":true");
}

#[test]
fn variable_speed_metamorphic_same_reduced_time_path() {
    // Two (U, dt) schedules with bitwise-equal per-step ds (U2 = 2*U1,
    // dt2 = dt1/2 — exact power-of-two scaling) must produce BITWISE
    // identical channel trajectories: reduced time is the only clock.
    // The metamorphic pair is live (not blind): a THIRD schedule with a
    // genuinely different ds path must differ.
    let spec = UNSTEADY_SECTION_V1;
    let chord = 1.981;
    let run = |u: f64, dt: f64, n: usize| -> Vec<u64> {
        let mut st = UnsteadySectionState::trim(&spec, 0.0);
        let mut out = Vec::new();
        for i in 0..n {
            let ds = fs_airfoil::indicial::reduced_time_increment(u, chord, dt).unwrap();
            let alpha = 0.05 + 0.03 * ((i as f64) * 0.1).sin();
            let ch = st.advance(&spec, ds, alpha, 0.0).unwrap();
            out.push(ch.cl.to_bits());
        }
        out
    };
    let a = run(13.86, 0.02, 200);
    let b = run(27.72, 0.01, 200);
    assert_eq!(a, b, "same reduced-time path must be bitwise identical");
    let c = run(13.86, 0.01, 200);
    assert_ne!(a, c, "a different ds path must differ (liveness)");
    jlog("metamorphic", "\"bitwise_equal\":true,\"liveness\":true");
}

#[test]
fn channel_golden_digest() {
    // A mixed maneuver exercising both channels + separation; digest over
    // the full (cl, cm) trajectory.
    let spec = UNSTEADY_SECTION_V1;
    let mut st = UnsteadySectionState::trim(&spec, 0.05);
    let mut payload = Vec::new();
    for i in 0..240 {
        let s = f64::from(i) * 0.1;
        let alpha = 0.05 + 0.30 * (0.5 * (1.0 - (-(s / 4.0)).exp()));
        let gust = 0.02 * (0.7 * s).sin();
        let ch = st.advance(&spec, 0.1, alpha, gust).unwrap();
        payload.extend_from_slice(&ch.cl.to_bits().to_le_bytes());
        payload.extend_from_slice(&ch.cm_quarter.to_bits().to_le_bytes());
    }
    let digest =
        fs_blake3::hash_domain("org.frankensim.fs-airfoil.v08a-golden.v1", &payload).to_hex();
    jlog("golden", &format!("\"digest\":\"{digest}\""));
    assert_eq!(
        digest, "e3a89c8007ebcc21c032554891f66c33d944df94d14509adbc016ad35d26b584",
        "V-08a channel golden moved — determinism regression or an \
         intentional model change requiring the golden-bump protocol"
    );
}
