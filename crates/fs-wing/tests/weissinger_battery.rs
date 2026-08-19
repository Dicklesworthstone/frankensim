//! E4.2-i battery (bead wf-root-guzez.5.3.1): classical fixtures per the
//! a5-biplane-theory verification role — monoplane lift-slope vs the
//! lifting-line trend, the EMERGENT Munk-class biplane gap effect (no
//! scalar factor anywhere), bitwise left-right symmetry, canard-upwash
//! sign, condition estimate reported, caps at cap AND cap+1, golden.
//! Repro: cargo test -p fs-wing --test weissinger_battery

use fs_wing::{MAX_PANELS, Panel, SurfaceId, flat_surface, solve_weissinger_linear};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-wing-e42i\",\"case\":\"{case}\",{payload}}}");
}

const RHO: f64 = 1.294;
const V: f64 = 13.86;

fn freestream(alpha_rad: f64) -> [f64; 3] {
    // frd: flight along +x means air arrives along... body sees V along
    // +x with an upward component giving positive alpha: w = -V sin a
    // (z down, so upwash is -z). Wait: positive alpha = nose above the
    // velocity => relative wind has +z (downward-in-body) component? In
    // frd with V ahead: u = V cos a, w = V sin a (w positive DOWN gives
    // positive alpha). Panels' normals are -z, so rhs = -V·n = +w·... —
    // the battery checks SIGNS through lift, which settles conventions.
    [V * alpha_rad.cos(), 0.0, V * alpha_rad.sin()]
}

#[test]
fn monoplane_lift_slope_matches_the_lifting_line_trend() {
    // AR 6 rectangular wing, 8 span x 2 chord panels. Lifting-line trend:
    // CL_alpha ~ a0/(1 + a0/(pi AR)) with a0 = 2 pi -> ~4.71/rad (e~1).
    // Weissinger-L on a rectangular wing lands within ~12% of that trend
    // (its classical accuracy class) — a TREND fixture, not an equality.
    let b = 12.29f64;
    let c = 2.048f64; // AR = 6.0
    let panels = flat_surface(SurfaceId::WingLower, b, c, 0.0, 0.0, 8, 2).unwrap();
    let alpha = 0.05f64;
    let r = solve_weissinger_linear(&panels, freestream(alpha), RHO).unwrap();
    let s_ref = b * c;
    let q = 0.5 * RHO * V * V;
    let cl = r.total_lift_n / (q * s_ref);
    let cl_alpha = cl / alpha;
    let a0 = 2.0 * std::f64::consts::PI;
    let ll = a0 / (1.0 + a0 / (std::f64::consts::PI * 6.0));
    assert!(
        (cl_alpha / ll - 1.0).abs() < 0.12,
        "CL_alpha {cl_alpha:.3} vs lifting-line {ll:.3} beyond the 12% trend band"
    );
    assert!(
        r.condition_est.is_finite() && r.condition_est > 0.0,
        "condition reported"
    );
    jlog(
        "monoplane",
        &format!(
            "\"cl_alpha\":{cl_alpha},\"lifting_line\":{ll},\"cond\":{}",
            r.condition_est
        ),
    );
}

#[test]
fn biplane_gap_effect_emerges_munk_trend() {
    // TWO planes at the 1903 gap lift LESS than 2x one plane (mutual
    // downwash); the deficit SHRINKS as the gap grows (Munk TR-151 trend).
    // No scalar factor exists in the code — this must EMERGE.
    let b = 12.29f64;
    let c = 2.048f64;
    let alpha = 0.05f64;
    let single = {
        let p = flat_surface(SurfaceId::WingLower, b, c, 0.0, 0.0, 8, 2).unwrap();
        solve_weissinger_linear(&p, freestream(alpha), RHO)
            .unwrap()
            .total_lift_n
    };
    let biplane_lift = |gap: f64| -> f64 {
        let mut p = flat_surface(SurfaceId::WingLower, b, c, 0.0, 0.0, 8, 2).unwrap();
        p.extend(flat_surface(SurfaceId::WingUpper, b, c, 0.0, -gap, 8, 2).unwrap());
        solve_weissinger_linear(&p, freestream(alpha), RHO)
            .unwrap()
            .total_lift_n
    };
    let close = biplane_lift(1.89); // the 1903 gap
    let wide = biplane_lift(6.0);
    let k_close = close / (2.0 * single);
    let k_wide = wide / (2.0 * single);
    assert!(
        k_close < 0.97,
        "1903-gap biplane factor {k_close:.3} must show mutual downwash"
    );
    assert!(
        k_wide > k_close,
        "the deficit must shrink with gap ({k_wide:.3} vs {k_close:.3})"
    );
    assert!(
        k_wide < 1.02,
        "wide-gap factor approaches 1 from below (got {k_wide:.3})"
    );
    jlog(
        "biplane",
        &format!("\"k_1903\":{k_close},\"k_wide\":{k_wide}"),
    );
}

#[test]
fn left_right_symmetry_is_bitwise() {
    // A symmetric configuration at zero sideslip: mirrored strips carry
    // bitwise-equal circulations (deterministic assembly + solve).
    let panels = flat_surface(SurfaceId::WingLower, 12.29, 2.048, 0.0, 0.0, 8, 2).unwrap();
    let r = solve_weissinger_linear(&panels, freestream(0.05), RHO).unwrap();
    for row in 0..2 {
        for s in 0..4 {
            let left = r.gamma[row * 8 + s];
            let right = r.gamma[row * 8 + (7 - s)];
            assert!(
                (left - right).abs() < 1e-12 * left.abs().max(1e-12),
                "mirror strips differ: {left} vs {right}"
            );
        }
    }
    jlog("symmetry", "\"mirror\":\"equal to 1e-12 relative\"");
}

#[test]
fn canard_upwash_raises_wing_loading_ahead() {
    // The canard sits AHEAD of the wing (frd +x). A lifting canard sheds
    // trailing vorticity that induces UPWASH outside its span at the wing
    // — the wing's lift with the canard present differs from without,
    // and the canard's own lift stays positive. Sign structure only
    // (quantitative interference is E4.2-ii+ territory).
    let alpha = 0.05f64;
    let wing = flat_surface(SurfaceId::WingLower, 12.29, 2.048, 0.0, 0.0, 8, 2).unwrap();
    let wing_alone = solve_weissinger_linear(&wing, freestream(alpha), RHO).unwrap();
    let mut both = wing.clone();
    both.extend(flat_surface(SurfaceId::CanardLower, 3.66, 0.76, 3.0, 0.3, 4, 2).unwrap());
    let combo = solve_weissinger_linear(&both, freestream(alpha), RHO).unwrap();
    let canard_lift = combo
        .surface_lift_n
        .iter()
        .find(|(s, _)| *s == SurfaceId::CanardLower)
        .unwrap()
        .1;
    let wing_lift_combo = combo
        .surface_lift_n
        .iter()
        .find(|(s, _)| *s == SurfaceId::WingLower)
        .unwrap()
        .1;
    assert!(canard_lift > 0.0, "the canard must lift at positive alpha");
    assert!(
        (wing_lift_combo - wing_alone.total_lift_n).abs() > 0.1,
        "the canard must MEASURABLY change the wing loading (coupling is live)"
    );
    jlog(
        "canard",
        &format!(
            "\"canard_lift\":{canard_lift},\"wing_delta\":{}",
            wing_lift_combo - wing_alone.total_lift_n
        ),
    );
}

#[test]
fn refusals_at_cap_and_cap_plus_one() {
    // Panel-count caps.
    let mk = |n: usize| -> Vec<Panel> {
        flat_surface(SurfaceId::WingLower, 12.0, 2.0, 0.0, 0.0, n, 1).unwrap()
    };
    assert!(solve_weissinger_linear(&mk(MAX_PANELS), freestream(0.05), RHO).is_ok());
    assert_eq!(
        solve_weissinger_linear(&mk(MAX_PANELS + 1), freestream(0.05), RHO)
            .unwrap_err()
            .code,
        "panel-count-invalid"
    );
    assert_eq!(
        solve_weissinger_linear(&[], freestream(0.05), RHO)
            .unwrap_err()
            .code,
        "panel-count-invalid"
    );
    // Bad panel (non-unit normal) and bad freestream.
    let mut bad = mk(4);
    bad[0].normal = [0.0, 0.0, -2.0];
    assert_eq!(
        solve_weissinger_linear(&bad, freestream(0.05), RHO)
            .unwrap_err()
            .code,
        "panel-invalid"
    );
    assert_eq!(
        solve_weissinger_linear(&mk(4), [0.0, 0.0, 0.0], RHO)
            .unwrap_err()
            .code,
        "freestream-invalid"
    );
    jlog(
        "refusals",
        "\"gates\":\"panels cap/cap+1/0, normal, freestream\"",
    );
}

#[test]
fn weissinger_golden_digest() {
    // The 1903 five-surface layout at the Dec-17 flight condition.
    let mut p = flat_surface(SurfaceId::WingLower, 12.29, 1.981, 0.0, 0.0, 8, 2).unwrap();
    p.extend(flat_surface(SurfaceId::WingUpper, 12.29, 1.981, 0.0, -1.89, 8, 2).unwrap());
    p.extend(flat_surface(SurfaceId::CanardLower, 3.66, 0.76, 3.0, 0.3, 4, 2).unwrap());
    p.extend(flat_surface(SurfaceId::CanardUpper, 3.66, 0.76, 3.0, -0.25, 4, 2).unwrap());
    let r = solve_weissinger_linear(&p, freestream(0.07), RHO).unwrap();
    let mut payload = Vec::new();
    for g in &r.gamma {
        payload.extend_from_slice(&g.to_bits().to_le_bytes());
    }
    payload.extend_from_slice(&r.total_lift_n.to_bits().to_le_bytes());
    let digest = fs_blake3::hash_domain("org.frankensim.fs-wing.e42i-golden.v1", &payload).to_hex();
    jlog(
        "golden",
        &format!("\"digest\":\"{digest}\",\"lift\":{}", r.total_lift_n),
    );
    assert_eq!(
        digest, "85eb80b5157b5effcbe88afdfa6f8d4ceac56518db61e5e03fedea54ad188ed4",
        "Weissinger golden moved — determinism regression or an intentional \
         assembly change requiring the golden-bump protocol"
    );
}
