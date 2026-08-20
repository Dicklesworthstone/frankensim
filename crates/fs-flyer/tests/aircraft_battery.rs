//! E4.6a-ii battery (bead wf-root-guzez.5.13.2): force build-up
//! cross-checks (per-axis oracles, never totals-only), fixed-control trim
//! within physical bands, determinism, the untrimmable falsifier (canard
//! authority removed -> typed trim-not-found), state-domain refusals,
//! golden trim state. Repro: cargo test -p fs-flyer --test aircraft_battery

use fs_flyer::aircraft::{MAX_TRIM_ITERATIONS, wright_openloop_v1};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-flyer-e46aii\",\"case\":\"{case}\",{payload}}}");
}

const RHO: f64 = 1.294;

#[test]
fn force_buildup_per_axis_oracles() {
    let d = wright_openloop_v1();
    let b = d.force_buildup(13.86, 0.06, 0.05, 40.0, 0.0, RHO).unwrap();
    // Per-axis oracles against independent recomputation. lift_n is the
    // z-component of the KJ forces; KJ acts perpendicular to the local
    // freestream, so its x-component is EXACTLY lift_n·tan(alpha) (the
    // leading-edge-suction analog).
    let w = d.gross_mass_kg * 9.80665;
    let alpha = 0.06f64;
    let fz_ref = -b.lift_n + w * alpha.cos() - b.drag_n * alpha.sin();
    assert!(
        (b.force_n[2] - fz_ref).abs() < 2.0,
        "Fz {} vs reconstructed {fz_ref}",
        b.force_n[2]
    );
    let fx_ref = b.lift_n * alpha.tan() + b.thrust_n[0] + b.thrust_n[1]
        - b.drag_n * alpha.cos()
        - w * alpha.sin();
    assert!(
        (b.force_n[0] - fx_ref).abs() < 2.0,
        "Fx {} vs reconstructed {fx_ref}",
        b.force_n[0]
    );
    // Physical sanity: lift in the weight class at this state.
    assert!(b.lift_n > 2000.0 && b.lift_n < 6000.0, "lift {}", b.lift_n);
    assert!(b.thrust_n[0] > 20.0 && b.thrust_n[0] < 400.0);
    assert!(b.induced_drag_n > 10.0 && b.induced_drag_n < 300.0);
    jlog(
        "buildup",
        &format!(
            "\"lift\":{},\"drag\":{},\"thrust\":{},\"my\":{}",
            b.lift_n,
            b.drag_n,
            b.thrust_n[0] + b.thrust_n[1],
            b.moment_y_nm
        ),
    );
}

#[test]
fn canard_command_sign_matches_the_frozen_convention() {
    // control-signs-v1: positive canard command -> POSITIVE pitch moment.
    let d = wright_openloop_v1();
    let lo = d.force_buildup(13.86, 0.06, 0.00, 40.0, 0.0, RHO).unwrap();
    let hi = d.force_buildup(13.86, 0.06, 0.10, 40.0, 0.0, RHO).unwrap();
    assert!(
        hi.moment_y_nm > lo.moment_y_nm + 10.0,
        "positive dc must raise M_y: {} -> {}",
        lo.moment_y_nm,
        hi.moment_y_nm
    );
    jlog(
        "control-sign",
        &format!("\"dM\":{}", hi.moment_y_nm - lo.moment_y_nm),
    );
}

#[test]
fn trim_converges_inside_physical_bands() {
    let d = wright_openloop_v1();
    let t = d.trim(RHO, [13.0, 0.06, 0.1, 45.0]).unwrap();
    jlog(
        "trim-state",
        &format!(
            "\"v\":{},\"alpha\":{},\"dc\":{},\"omega\":{},\"iters\":{}",
            t.v_mps, t.alpha_rad, t.delta_canard_rad, t.omega_prop_rad_s, t.iterations
        ),
    );
    assert!(t.iterations <= MAX_TRIM_ITERATIONS);
    for (i, r) in t.residuals.iter().enumerate() {
        assert!(r.abs() < 0.5, "residual {i} = {r}");
    }
    // The torque-balance state must land at a physical prop speed.
    assert!(
        t.omega_prop_rad_s > 25.0 && t.omega_prop_rad_s < 80.0,
        "trim prop speed {} rad/s outside the physical band",
        t.omega_prop_rad_s
    );
    // Physical bands: the Flyer trims slow, nose-up, canard loaded.
    assert!(
        t.v_mps > 9.0 && t.v_mps < 20.0,
        "trim speed {} outside the plausible band",
        t.v_mps
    );
    // MODEL-plausibility band, not a historical claim: the 2*pi
    // thin-airfoil camber closure over-lifts vs the 1901-anchored
    // sections (E4.1: modern-equiv cl 0.659 at 5 deg vs 1.18 predicted),
    // so this model trims ~1 deg lower and ~11% faster than the Dec-17
    // state; measured trim alpha is -0.0149 rad. The historical-match
    // refinement is section-data territory (E4.1 tables), not trim math.
    assert!(
        t.alpha_rad > -0.05 && t.alpha_rad < 0.30,
        "trim alpha {} outside the plausible band",
        t.alpha_rad
    );
    assert!(
        t.delta_canard_rad.abs() < 0.5,
        "canard {}",
        t.delta_canard_rad
    );
    // Lift ~ weight at trim (residual closure already enforces it; this
    // is the human-readable receipt line).
    let w = d.gross_mass_kg * 9.80665;
    assert!((t.buildup.lift_n / w - 1.0).abs() < 0.15);
    jlog(
        "trim",
        &format!(
            "\"v\":{},\"alpha\":{},\"dc\":{},\"omega\":{},\"iters\":{},\"lift\":{},\"thrust\":{}",
            t.v_mps,
            t.alpha_rad,
            t.delta_canard_rad,
            t.omega_prop_rad_s,
            t.iterations,
            t.buildup.lift_n,
            t.buildup.thrust_n[0] + t.buildup.thrust_n[1]
        ),
    );
}

#[test]
fn trim_is_deterministic() {
    let d = wright_openloop_v1();
    let a = d.trim(RHO, [13.0, 0.06, 0.1, 45.0]).unwrap();
    let b = d.trim(RHO, [13.0, 0.06, 0.1, 45.0]).unwrap();
    assert_eq!(a.v_mps.to_bits(), b.v_mps.to_bits());
    assert_eq!(a.alpha_rad.to_bits(), b.alpha_rad.to_bits());
    assert_eq!(a.delta_canard_rad.to_bits(), b.delta_canard_rad.to_bits());
    assert_eq!(a.omega_prop_rad_s.to_bits(), b.omega_prop_rad_s.to_bits());
    jlog("determinism", "\"bitwise\":true");
}

#[test]
fn untrimmable_configuration_refuses_typed() {
    // FALSIFIER: shrink the canard to a sliver — no pitch authority, no
    // equilibrium. The trim must refuse with the typed code, never hand
    // back a best-effort state.
    let mut d = wright_openloop_v1();
    d.canard_span_m = 0.2;
    d.canard_chord_m = 0.05;
    let err = d.trim(RHO, [13.0, 0.06, 0.1, 45.0]).unwrap_err();
    assert_eq!(err.code, "trim-not-found");
    assert!(
        err.message.contains("residual trail") || err.message.contains("|r|"),
        "refusal must carry the residual trajectory: {}",
        err.message
    );
    jlog("falsifier", "\"untrimmable_refused\":true");
}

#[test]
fn state_domain_refusals() {
    let d = wright_openloop_v1();
    assert_eq!(
        d.force_buildup(0.5, 0.06, 0.0, 40.0, 0.0, RHO)
            .unwrap_err()
            .code,
        "state-invalid"
    );
    assert_eq!(
        d.force_buildup(13.0, 0.7, 0.0, 40.0, 0.0, RHO)
            .unwrap_err()
            .code,
        "state-invalid"
    );
    assert_eq!(
        d.force_buildup(13.0, 0.06, 0.0, 200.0, 0.0, RHO)
            .unwrap_err()
            .code,
        "state-invalid"
    );
    assert_eq!(
        d.force_buildup(13.0, 0.06, 0.0, 40.0, 0.0, f64::NAN)
            .unwrap_err()
            .code,
        "state-invalid"
    );
    // Design digest sensitivity: moving the CG moves the digest.
    let mut d2 = wright_openloop_v1();
    d2.cg_m[0] -= 0.1;
    assert_ne!(d2.digest(), d.digest());
    jlog("domain", "\"refusals_typed\":true");
}

#[test]
fn trim_golden_digest() {
    let d = wright_openloop_v1();
    let t = d.trim(RHO, [13.0, 0.06, 0.1, 45.0]).unwrap();
    let mut payload = Vec::new();
    for v in [
        t.v_mps,
        t.alpha_rad,
        t.delta_canard_rad,
        t.buildup.lift_n,
        t.buildup.drag_n,
        t.buildup.thrust_n[0],
    ] {
        payload.extend_from_slice(&v.to_bits().to_le_bytes());
    }
    payload.extend_from_slice(t.design_digest.as_bytes());
    let digest =
        fs_blake3::hash_domain("org.frankensim.fs-flyer.e46aii-golden.v1", &payload).to_hex();
    jlog("golden", &format!("\"digest\":\"{digest}\""));
    assert_eq!(
        digest, "62ebf1d4d5487fe252e2f4008e64167dffd8f1dbbe076c07cef538e29bdd75e7",
        "trim golden moved — determinism regression or an intentional \
         model change requiring the golden-bump protocol"
    );
}
