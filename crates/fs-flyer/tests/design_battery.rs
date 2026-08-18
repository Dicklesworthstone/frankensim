//! E3.1 design battery (bead wf-root-guzez.4.1). Per-item oracles: the
//! reference config pins mass/CG against dossier values and reproduces the
//! published single-lineage inertias within the documented band; admission
//! refuses non-physical inputs at cap AND cap+1; the derived panel matches
//! hand calculations; design digest golden (measure-then-pin).
//! Repro: cargo test -p fs-flyer --test design_battery

use fs_flyer::{
    FlyerDesign, LateralControlTopology, MAX_PILOT_KG, MIN_PILOT_KG, PUBLISHED_IXX_KGM2,
    PUBLISHED_IYY_KGM2, PUBLISHED_IZZ_KGM2,
};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-flyer-design\",\"case\":\"{case}\",{payload}}}");
}

#[test]
fn reference_mass_cg_pinned_to_dossier() {
    let d = FlyerDesign::reference_1903();
    d.admit().expect("the reference design must admit");
    let b = d.mass_build_up().unwrap();
    // Gross 340.2 kg (flyer-reference: 605 lb empty + 145 lb pilot).
    assert!((b.gross_kg - 340.17).abs() < 0.5, "gross {}", b.gross_kg);
    // Gross CG at 29.7% chord aft of the wing LE (canard-mechanics dossier):
    // x_cg = −0.297·1.981 = −0.588 m.
    let cg_xc = -b.cg_m[0] / 1.981;
    assert!(
        (cg_xc - 0.297).abs() < 0.01,
        "CG at {cg_xc} x/c vs dossier 0.297"
    );
    // Lateral symmetry: |y_cg| small (engine right, pilot left).
    assert!(b.cg_m[1].abs() < 0.06, "y_cg {}", b.cg_m[1]);
    jlog(
        "mass-cg",
        &format!("\"gross\":{},\"cg_xc\":{cg_xc}", b.gross_kg),
    );
}

#[test]
fn inertia_reproduces_published_lineage_within_band() {
    // The published inertias are a SINGLE-LINEAGE reconstruction (Jex &
    // Culick 85-1804 via UIUC); our component build-up is an independent
    // reconstruction CALIBRATED to land within a documented ±15% band —
    // agreement inside the band is a cross-check, not a validation claim.
    let b = FlyerDesign::reference_1903().mass_build_up().unwrap();
    let [ixx, iyy, izz] = b.inertia_kgm2;
    for (got, want, name) in [
        (ixx, PUBLISHED_IXX_KGM2, "Ixx"),
        (iyy, PUBLISHED_IYY_KGM2, "Iyy"),
        (izz, PUBLISHED_IZZ_KGM2, "Izz"),
    ] {
        let ratio = got / want;
        assert!(
            (0.85..=1.15).contains(&ratio),
            "{name} {got:.0} vs published {want:.0} (ratio {ratio:.3}) outside ±15%"
        );
        jlog(
            "inertia",
            &format!("\"axis\":\"{name}\",\"got\":{got:.1},\"published\":{want:.1}"),
        );
    }
    // Slender-biplane ordering must hold (per-item, not totals).
    assert!(iyy < ixx && ixx < izz, "ordering Iyy < Ixx < Izz violated");
}

#[test]
fn admission_refuses_at_cap_and_cap_plus_one() {
    let base = FlyerDesign::reference_1903;
    // Pilot mass: both edges admitted, one float outside refused.
    let mut d = base();
    d.pilot_mass_kg = MIN_PILOT_KG;
    assert!(d.admit().is_ok());
    d.pilot_mass_kg = MAX_PILOT_KG;
    assert!(d.admit().is_ok());
    d.pilot_mass_kg = f64::from_bits(MAX_PILOT_KG.to_bits() + 1);
    assert_eq!(d.admit().unwrap_err().code, "pilot-mass-outside-domain");
    d.pilot_mass_kg = f64::from_bits(MIN_PILOT_KG.to_bits() - 1);
    assert_eq!(d.admit().unwrap_err().code, "pilot-mass-outside-domain");
    // Area consistency: the real planform admits (2.7% off rectangular);
    // a fabricated 20%-off area refuses with the deviation stated.
    let mut bad_area = base();
    bad_area.wing.area_both_m2 = 2.0 * 12.29 * 1.981 * 1.25;
    let refusal = bad_area.admit().unwrap_err();
    assert_eq!(refusal.code, "area-inconsistent");
    assert!(refusal.message.contains('%'), "deviation must be stated");
    // Hinge axis outside the E1.5 prior refuses (the prior IS the domain).
    let mut bad_hinge = base();
    bad_hinge.canard.hinge_axis_xc = 0.24;
    assert_eq!(
        bad_hinge.admit().unwrap_err().code,
        "hinge-axis-outside-prior"
    );
    // Mass-spec falsifier: silently dropping a component turns red.
    let mut hidden = base();
    hidden.components.pop();
    assert_eq!(hidden.admit().unwrap_err().code, "mass-spec-mismatch");
    // Non-finite and negative-span refusals.
    let mut nan = base();
    nan.wing.span_m = f64::NAN;
    assert_eq!(nan.admit().unwrap_err().code, "non-finite-input");
    let mut neg = base();
    neg.wing.span_m = 0.0;
    assert_eq!(neg.admit().unwrap_err().code, "span-outside-domain");
    jlog(
        "admission",
        "\"caps\":\"pilot both edges, area, hinge prior, mass-spec\"",
    );
}

#[test]
fn derived_panel_matches_hand_calculations() {
    let d = FlyerDesign::reference_1903();
    let p = d.derived_panel(0.3).unwrap();
    let b = d.mass_build_up().unwrap();
    // HAND CALCULATION (same documented model, computed independently):
    // canard AC x = 2.231 + 0.75·0.762 = 2.8025 m; wing AC x = −0.49525 m.
    let x_ac_c = 2.231 + 0.75 * 0.762;
    let x_ac_w = -0.25 * 1.981;
    let l_c = x_ac_c - b.cg_m[0];
    let vc_hand = (4.46 * l_c) / (47.38 * 1.981);
    assert!(
        (p.canard_volume - vc_hand).abs() < 1e-12,
        "V_c {} vs {vc_hand}",
        p.canard_volume
    );
    assert!(
        (p.canard_volume - 0.161).abs() < 0.005,
        "V_c magnitude sanity"
    );
    // Naive NP by hand: (S_w·x_w + S_c·x_c)/(S_w+S_c), reported x/c aft.
    let x_np = (47.38 * x_ac_w + 4.46 * x_ac_c) / (47.38 + 4.46);
    let np_hand = -x_np / 1.981;
    assert!((p.neutral_point_xc_naive - np_hand).abs() < 1e-12);
    // Margins: fixed = NP − CG (aft-positive); free more negative under
    // overbalance (the floating canard amplifies the destabilizer).
    let cg_xc = -b.cg_m[0] / 1.981;
    assert!((p.static_margin_fixed - (np_hand - cg_xc)).abs() < 1e-12);
    assert!(
        p.static_margin_fixed < 0.0,
        "the Flyer must come out unstable"
    );
    assert!(
        p.static_margin_free < p.static_margin_fixed,
        "free-control must be worse"
    );
    // The naive two-surface NP sits AFT of the Culick-lineage 3.9%c value —
    // interference/flex effects (absent from the naive model) destabilize.
    // Recorded as a documented comparison, not reconciled.
    assert!(p.neutral_point_xc_naive > 0.039);
    // Hinge gradient: zero exactly at x_h = 1/4; hand value at 0.375.
    let mut at_quarter = d.clone();
    at_quarter.canard.hinge_axis_xc = 0.25;
    let p0 = at_quarter.derived_panel(0.0).unwrap();
    assert!(
        p0.hinge_moment_gradient.abs() < 1e-15,
        "balanced hinge ⇒ zero gradient"
    );
    let hand = 2.0 * std::f64::consts::PI * 0.125;
    assert!((p.hinge_moment_gradient - hand).abs() < 1e-12);
    assert!(
        p.hinge_moment_gradient > 0.0,
        "self-driving sign per the Orville evidence"
    );
    // Hinge-ratio domain refusal.
    assert_eq!(
        d.derived_panel(1.5).unwrap_err().code,
        "hinge-ratio-invalid"
    );
    jlog(
        "panel",
        &format!(
            "\"vc\":{},\"np_naive_xc\":{},\"margin_fixed\":{},\"margin_free\":{},\"ch_a\":{}",
            p.canard_volume,
            p.neutral_point_xc_naive,
            p.static_margin_fixed,
            p.static_margin_free,
            p.hinge_moment_gradient
        ),
    );
}

#[test]
fn design_digest_golden_and_sensitivity() {
    let d = FlyerDesign::reference_1903();
    let base = d.digest();
    assert_eq!(base, d.digest(), "digest must be deterministic");
    // Any field change moves the digest (identity sensitivity).
    let mut tweaked = d.clone();
    tweaked.canard.hinge_axis_xc += 1e-12;
    assert_ne!(base, tweaked.digest());
    let mut retopo = d.clone();
    retopo.lateral = LateralControlTopology::WarpIndependentRudder;
    assert_ne!(base, retopo.digest());
    jlog("golden", &format!("\"digest\":\"{base}\""));
    assert_eq!(
        base, "27502621e2bd2394be18b5c189b06e4847e1ff478c5d84956d6bbbf218780c7d",
        "reference design digest moved — identity regression or an intentional \
         reference change requiring the golden-bump protocol"
    );
}
