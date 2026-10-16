//! E4.9a battery (bead wf-root-guzez.5.21): the FIRST discrepancy
//! receipts recorded in the E10.1 harness format — attached monoplane
//! vs classical lifting line, biplane vs the Prandtl/Munk correction,
//! BEMT vs actuator-disk momentum, and the coupled-prop interference
//! delta. Differences REPORTED, never forced. Receipt caps at cap AND
//! cap+1; format stability golden.
//! Repro: cargo test -p fs-flyer --test referee_battery

use fs_flyer::referee::{
    DiscrepancyReceipt, MAX_REL_DISCREPANCY, RECEIPT_SCHEMA, referee_biplane_efficiency,
    referee_liftingline_cl, referee_momentum_thrust,
};
use fs_wing::nonlinear::{InfluenceOperator, StripRegime, StripSpec, solve_nonlinear};
use fs_wing::{SurfaceId, flat_surface};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-flyer-e49a\",\"case\":\"{case}\",{payload}}}");
}

const RHO: f64 = 1.294;
const V: f64 = 13.86;
const ALPHA: f64 = 0.05;
const CAMBER: f64 = 0.05;

fn closure(_s: usize, a: f64) -> (f64, StripRegime) {
    (
        2.0 * std::f64::consts::PI * (a + 2.0 * CAMBER),
        StripRegime::Attached,
    )
}

fn solve_lift(panels: &[fs_wing::Panel], strips: &[StripSpec]) -> f64 {
    let fs_v = [V * ALPHA.cos(), 0.0, V * ALPHA.sin()];
    let op = InfluenceOperator::build(panels, fs_v, RHO).unwrap();
    solve_nonlinear(&op, panels, strips, fs_v, RHO, &closure, None, None)
        .unwrap()
        .total_lift_n
}

fn monoplane() -> (Vec<fs_wing::Panel>, Vec<StripSpec>) {
    let p = flat_surface(SurfaceId::WingLower, 12.29, 1.981, 0.0, 0.0, 8, 2).unwrap();
    let strips = (0..8)
        .map(|s| StripSpec {
            panel_indices: vec![s, 8 + s],
            chord_m: 1.981,
            twist_rad: 0.0,
        })
        .collect();
    (p, strips)
}

#[test]
fn monoplane_attached_case_agrees_with_lifting_line() {
    let (p, strips) = monoplane();
    let production = solve_lift(&p, &strips);
    let q = 0.5 * RHO * V * V;
    let s_ref = 12.29 * 1.981;
    let ar = 12.29 / 1.981;
    let cl_ref = referee_liftingline_cl(ALPHA, CAMBER, ar, 1.0);
    let referee = q * s_ref * cl_ref;
    let r = DiscrepancyReceipt::new(
        "e49a-mono-attached",
        "wing_lift",
        "N",
        production,
        referee,
        "closed-form lifting line, no shared code with the panel solver",
        "formulation-band-0.15",
    )
    .unwrap();
    println!("{}", r.to_jsonl());
    // The declared band: Weissinger vs lifting line at AR 6 within 15%.
    assert!(
        r.rel_discrepancy.abs() < 0.15,
        "attached monoplane outside the declared band: {}",
        r.rel_discrepancy
    );
    jlog(
        "mono",
        &format!(
            "\"production\":{production},\"referee\":{referee},\"rel\":{}",
            r.rel_discrepancy
        ),
    );
}

#[test]
fn biplane_case_records_the_interference_discrepancy() {
    // Biplane: production coupled solve vs the classical Prandtl/Munk
    // correction applied to the monoplane referee. The delta is
    // RECORDED — the two formulations genuinely differ (bound-bound vs
    // classical trailing-only interference), and forcing agreement
    // would falsify one of them.
    let (mp, ms) = monoplane();
    let mono_lift = solve_lift(&mp, &ms);
    let mut p = flat_surface(SurfaceId::WingLower, 12.29, 1.981, 0.0, 0.0, 8, 2).unwrap();
    p.extend(flat_surface(SurfaceId::WingUpper, 12.29, 1.981, 0.0, -1.89, 8, 2).unwrap());
    let mut strips = Vec::new();
    for plane in 0..2 {
        let base = plane * 16;
        for s in 0..8 {
            strips.push(StripSpec {
                panel_indices: vec![base + s, base + 8 + s],
                chord_m: 1.981,
                twist_rad: 0.0,
            });
        }
    }
    let production = solve_lift(&p, &strips);
    // Referee: 2x monoplane scaled by the biplane efficiency at
    // g/b = 1.89/12.29.
    let eff = referee_biplane_efficiency(1.89 / 12.29);
    let referee = 2.0 * mono_lift * eff;
    let r = DiscrepancyReceipt::new(
        "e49a-biplane-interference",
        "biplane_lift",
        "N",
        production,
        referee,
        "Prandtl/Munk biplane correction on the closed-form monoplane",
        "reported-only",
    )
    .unwrap();
    println!("{}", r.to_jsonl());
    // Reported-only class: the receipt must EXIST and be sane, and the
    // discrepancy is data (measured on first run and logged).
    assert!(
        r.rel_discrepancy.abs() < 0.5,
        "biplane discrepancy implausibly large: {}",
        r.rel_discrepancy
    );
    jlog(
        "biplane",
        &format!(
            "\"production\":{production},\"referee\":{referee},\"rel\":{},\"munk_eff\":{eff}",
            r.rel_discrepancy
        ),
    );
}

#[test]
fn bemt_thrust_agrees_with_momentum_at_matched_power() {
    // Production: BEMT at the aircraft trim class (fs-airscrew).
    // Referee: actuator-disk momentum at the SAME absorbed power —
    // BEMT thrust must sit BELOW the momentum ideal (losses) but within
    // the declared class band.
    let rotor = fs_flyer::aircraft::wright_rotor_v1();
    let omega = 52.0;
    let sol = fs_airscrew::bemt_solve(&rotor, RHO, V, omega).unwrap();
    let power = sol.torque_nm * omega;
    let area = core::f64::consts::PI * rotor.radius_m * rotor.radius_m;
    let ideal = referee_momentum_thrust(power, V, RHO, area).unwrap();
    let r = DiscrepancyReceipt::new(
        "e49a-prop-momentum",
        "thrust",
        "N",
        sol.thrust_n,
        ideal,
        "actuator-disk momentum bisection at matched power, no BEMT code",
        "below-ideal-within-0.5",
    )
    .unwrap();
    println!("{}", r.to_jsonl());
    assert!(
        sol.thrust_n < ideal,
        "BEMT must not beat the momentum ideal: {} vs {ideal}",
        sol.thrust_n
    );
    assert!(
        r.rel_discrepancy > -0.5,
        "BEMT more than 50% below ideal: {}",
        r.rel_discrepancy
    );
    jlog(
        "prop",
        &format!(
            "\"bemt_n\":{},\"ideal_n\":{ideal},\"rel\":{},\"power_w\":{power}",
            sol.thrust_n, r.rel_discrepancy
        ),
    );
}

#[test]
fn receipt_caps_at_cap_and_cap_plus_one() {
    // At the discrepancy cap: admitted.
    let at = DiscrepancyReceipt::new("cap", "x", "-", 1.0 + MAX_REL_DISCREPANCY, 1.0, "n", "c");
    assert!(at.is_ok(), "cap admits");
    // One ulp past: refused.
    let past = DiscrepancyReceipt::new(
        "cap",
        "x",
        "-",
        (1.0 + MAX_REL_DISCREPANCY).next_up(),
        1.0,
        "n",
        "c",
    );
    assert_eq!(past.unwrap_err().code, "referee-receipt-invalid");
    assert_eq!(
        DiscrepancyReceipt::new("z", "x", "-", 1.0, 0.0, "n", "c")
            .unwrap_err()
            .code,
        "referee-receipt-invalid",
        "zero referee refuses"
    );
    assert_eq!(
        DiscrepancyReceipt::new("", "x", "-", 1.0, 1.0, "n", "c")
            .unwrap_err()
            .code,
        "referee-receipt-invalid",
        "empty case id refuses"
    );
    assert_eq!(
        referee_momentum_thrust(-1.0, 10.0, RHO, 5.0)
            .unwrap_err()
            .code,
        "referee-momentum-invalid"
    );
    jlog("caps", "\"cap_and_cap_plus_one\":true");
}

#[test]
fn receipt_format_is_stable_and_golden() {
    let r = DiscrepancyReceipt::new(
        "e49a-format-fixture",
        "lift",
        "N",
        1234.5,
        1200.0,
        "fixture",
        "reported-only",
    )
    .unwrap();
    let row = r.to_jsonl();
    assert!(row.starts_with(&format!("{{\"schema\":\"{RECEIPT_SCHEMA}\"")));
    for field in [
        "\"case_id\":",
        "\"quantity\":",
        "\"units\":",
        "\"production\":",
        "\"referee\":",
        "\"rel_discrepancy\":",
        "\"independence\":",
        "\"comparison_class\":",
    ] {
        assert!(row.contains(field), "harness field missing: {field}");
    }
    assert_eq!(r.digest(), r.digest(), "digest stable");
    let digest =
        fs_blake3::hash_domain("org.frankensim.fs-flyer.e49a-golden.v1", row.as_bytes()).to_hex();
    jlog("golden", &format!("\"digest\":\"{digest}\""));
    assert_eq!(
        digest, "1bfcedf0c824aef017b375b08e90b0714108b2aca41ef43216727732785735e6",
        "receipt-format golden moved — the E10.1 harness format is a \
         frozen contract; changes ride the golden-bump protocol"
    );
}
