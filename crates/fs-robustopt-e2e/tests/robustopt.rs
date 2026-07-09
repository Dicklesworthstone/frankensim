//! End-to-end battery: every family's global optimum is SOS-PROVEN, robustness
//! reorders the ranking, and the headline claim is honestly Estimated.

use fs_evidence::{Color, ColorRank};
use fs_robustopt_e2e::{Family, demo_families, run_campaign};

#[test]
fn global_optima_are_sos_proven_and_robustness_reorders() {
    let report = run_campaign(&demo_families(), 0.9, 2.0, 41);
    // every convex family's global optimum carries an SOS proof.
    assert_eq!(report.certified_count, 3);
    for v in &report.families {
        assert!(
            matches!(v.nominal_color, Color::Verified { .. }),
            "{} unproven",
            v.name
        );
        assert!((v.x_star - 2.0).abs() < 1e-9, "{} x*={}", v.name, v.x_star);
    }
    // nominal optima match the closed form c − b²/4a.
    let cost = |n: &str| {
        report
            .families
            .iter()
            .find(|v| v.name == n)
            .unwrap()
            .nominal_cost
    };
    assert!((cost("champion") - 1.2).abs() < 1e-9);
    assert!((cost("flat") - 2.0).abs() < 1e-9);
    assert!((cost("sharp") - 2.0).abs() < 1e-9);
    // the LOWEST-nominal family is champion, but ROBUSTNESS crowns the flat one.
    assert_eq!(report.nominal_winner, "champion");
    assert_eq!(report.robust_winner, "flat");
    assert!(report.robustness_reorders);
    // the robust winner's worst case beats the nominal winner's worst case.
    let robust = |n: &str| {
        report
            .families
            .iter()
            .find(|v| v.name == n)
            .unwrap()
            .robust_cost
    };
    assert!(
        robust("flat") < robust("champion"),
        "flat {} !< champion {}",
        robust("flat"),
        robust("champion")
    );
    // NO LAUNDERING: a CVaR robust ranking is a sample statistic → Estimated.
    assert_eq!(report.headline_rank, ColorRank::Estimated);
    println!(
        "{{\"campaign\":\"proofrobust\",\"proven\":{},\"nominal_winner\":\"{}\",\
         \"robust_winner\":\"{}\",\"champion_robust\":{:.3},\"flat_robust\":{:.3}}}",
        report.certified_count,
        report.nominal_winner,
        report.robust_winner,
        robust("champion"),
        robust("flat"),
    );
}

#[test]
fn a_downward_family_is_not_certified() {
    // an unbounded-below "family" (a < 0) has no global minimum: no SOS proof.
    let report = run_campaign(
        &[
            Family::new("bad", -1.0, 0.0, 0.0),
            Family::new("ok", 1.0, 0.0, 1.0),
        ],
        0.9,
        1.0,
        21,
    );
    assert_eq!(report.certified_count, 1);
    let bad = report.families.iter().find(|v| v.name == "bad").unwrap();
    assert!(matches!(bad.nominal_color, Color::Estimated { .. }));
}

#[test]
fn the_campaign_is_deterministic() {
    let a = run_campaign(&demo_families(), 0.9, 2.0, 41);
    let b = run_campaign(&demo_families(), 0.9, 2.0, 41);
    assert_eq!(a.robust_winner, b.robust_winner);
    assert_eq!(
        a.families[0].robust_cost.to_bits(),
        b.families[0].robust_cost.to_bits()
    );
}
