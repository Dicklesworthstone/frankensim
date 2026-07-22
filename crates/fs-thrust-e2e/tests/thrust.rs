//! End-to-end battery for the CertQD-Thrust campaign: the physics self-propels,
//! the campaign illuminates a diverse family under certify-or-escalate inside
//! the calibration set's per-axis support, every published drift is honestly
//! `Estimated`, and the whole run is deterministic.

use fs_evidence::ColorRank;
use fs_thrust_e2e::{
    CampaignBudget, Design, calibration_designs, calibration_support, design_grid, run_campaign,
    simulate_thrust,
};

#[test]
fn a_four_vortex_thruster_self_propels_and_conserves_impulse() {
    // physics sanity: a leading dipole with a weak trailing pair drifts +x, and
    // a converged inviscid sim conserves the exact linear impulse.
    let d = Design {
        gamma: 1.0,
        d: 0.6,
        l: 1.5,
        ratio: 0.3,
    };
    let res = simulate_thrust(&d, 300, 0.02, 0.05);
    assert!(
        res.drift > 0.1,
        "thruster should self-propel +x, drift={}",
        res.drift
    );
    assert!(
        res.impulse_error < 1e-2,
        "impulse leaked: {}",
        res.impulse_error
    );
}

#[test]
fn the_campaign_illuminates_a_screened_diverse_family() {
    let report = run_campaign(&CampaignBudget::default());
    // surface the headline numbers FIRST (structured, for the run log).
    println!(
        "{{\"campaign\":\"certqd-thrust\",\"niches\":{},\"coverage\":{:.3},\"qd_score\":{:.3},\
         \"best_drift\":{:.3},\"best\":{:?},\"screened\":{},\"unscreened\":{},\"full_sims\":{},\
         \"short_sims\":{},\"steps_spent\":{},\"steps_all_full\":{},\"step_savings\":{:.3},\
         \"band\":{:.4},\"screened_drift_hull\":{:?}}}",
        report.num_elites,
        report.coverage,
        report.qd_score,
        report.best_drift,
        report.best,
        report.conservation_screened_elites,
        report.unscreened_elites,
        report.full_sims,
        report.short_sims,
        report.steps_spent,
        report.steps_all_full,
        1.0 - report.steps_spent as f64 / report.steps_all_full as f64,
        report.band_half_width,
        report.conservation_screened_drift_hull,
    );
    for e in &report.atlas {
        println!(
            "ATLAS {:.4} {:.4} {:.5} {} {:?}",
            e.budget, e.length, e.drift, e.conservation_screened, e.rank
        );
    }
    // ILLUMINATION: a diverse archive of niches, not a single optimum.
    assert!(
        report.num_elites >= 5,
        "too few niches: {}",
        report.num_elites
    );
    assert!(report.coverage > 0.0 && report.qd_score > 0.0);
    // the best thruster genuinely self-propels.
    assert!(report.best_drift > 0.0, "best drift {}", report.best_drift);
    // CERTIFY-OR-ESCALATE: both fidelities were used, and it beat all-full cost.
    assert!(report.short_sims > 0, "surrogate never used");
    assert!(report.full_sims > 0, "nothing escalated");
    assert!(
        report.steps_spent < report.steps_all_full,
        "no savings: {} vs {}",
        report.steps_spent,
        report.steps_all_full
    );
    // …and the saving is not a promise, it is arithmetic: the served designs
    // must repay the eight paired calibration sims. Check the identity that
    // actually governs it (bead .2.37).
    let budget = CampaignBudget::default();
    let calibration = calibration_designs().len();
    let repaid = report.short_sims * (budget.full_steps - budget.short_steps);
    let overhead = calibration * (budget.short_steps + budget.full_steps);
    assert!(
        repaid > overhead,
        "served designs must repay calibration: {repaid} vs {overhead}"
    );
    assert_eq!(
        report.steps_spent,
        overhead + report.short_sims * budget.short_steps + report.full_sims * budget.full_steps
    );
    // SCREEN: escalated sims that conserved impulse passed the screen.
    assert!(
        report.conservation_screened_elites > 0,
        "no screened elites"
    );
    assert_eq!(
        report.conservation_screened_elites + report.unscreened_elites,
        report.num_elites
    );
    let (lo, hi) = report
        .conservation_screened_drift_hull
        .expect("a screened elite exists");
    assert!(lo <= hi, "screened drift hull [{lo}, {hi}]");
    // Once the validity domain covers every calibrated axis, the surrogate's
    // 18 in-domain designs (all at d = 0.7) never out-drift the tighter-spaced
    // ones, so no surrogate estimate reaches the atlas: every elite here is a
    // screened full sim. That is the honest reading of a narrow calibration
    // set, and it is pinned so a wider domain has to move it.
    assert_eq!(
        report.unscreened_elites, 0,
        "every default-budget elite is a screened full sim"
    );
    // NO LAUNDERING: the campaign claim is Estimated — the impulse screen is a
    // diagnostic on a different functional, so it cannot certify a drift.
    assert_eq!(report.campaign_rank, ColorRank::Estimated);
    // the reproducible notebook carries the story and is content-addressed.
    assert!(report.notebook_markdown.contains("CertQD-Thrust"));
    assert!(report.notebook_markdown.contains("best_drift"));
    assert_ne!(report.content_hash, 0);
}

/// REGRESSION (bead `frankensim-extreal-program-f85xj.2.30`): the campaign's
/// public atlas must not publish a single interval-certified elite. The drift
/// comes from unchecked RK4 and the only trust signal is an impulse-conservation
/// residual on a DIFFERENT functional, so `Verified` is unearnable here. Before
/// the fix this run published 19 `Verified` elites and a "certified drift
/// envelope" composed from their invented ±1e-9 bands.
#[test]
fn no_elite_claims_an_interval_certificate_from_the_impulse_screen() {
    let report = run_campaign(&CampaignBudget::default());
    for e in &report.atlas {
        assert_eq!(
            e.rank,
            ColorRank::Estimated,
            "elite at (budget {:.3}, length {:.3}) claims {:?}",
            e.budget,
            e.length,
            e.rank
        );
    }
    assert_eq!(report.campaign_rank, ColorRank::Estimated);
    // This is a DOWNGRADE of the claim, not a deletion of the evidence: the
    // screen still runs, still admits, and still publishes its residual hull.
    // (That it still REFUSES is pinned by
    // `malformed_conservation_tolerance_cannot_mint_screened_elites`.)
    assert_eq!(report.conservation_screened_elites, report.num_elites);
    let (lo, hi) = report
        .conservation_screened_drift_hull
        .expect("screened elites publish their residual hull");
    assert!(lo.is_finite() && hi.is_finite() && lo < hi);
}

/// REGRESSION (bead `frankensim-extreal-program-f85xj.2.29`): no design may be
/// served by the surrogate at a gene value the calibration residuals never saw.
/// Before the fix the validity test was the 2-D descriptor hull only, so 63 of
/// 84 surrogate-served designs sat at transverse spacings `d ∈ {0.4, 1.0, 1.3}`
/// — off-calibration on the axis that drives dipole self-advection — and each
/// received the `d = 0.7`-calibrated half-width as its uncertainty.
#[test]
fn every_surrogate_served_design_is_inside_the_calibration_support() {
    let support = calibration_support();
    let served: Vec<Design> = design_grid()
        .into_iter()
        .filter(|d| support.contains(d))
        .collect();
    for design in &served {
        assert!(
            (design.d - 0.7).abs() < 1e-12,
            "off-calibration spacing served: {design:?} (calibration support {:?})",
            support.d
        );
    }
    // The campaign's own accounting agrees with the predicate.
    let report = run_campaign(&CampaignBudget::default());
    assert_eq!(report.short_sims, served.len());
    assert_eq!(report.full_sims, design_grid().len() - served.len());
}

#[test]
fn the_campaign_is_deterministic() {
    let a = run_campaign(&CampaignBudget::default());
    let b = run_campaign(&CampaignBudget::default());
    assert_eq!(a.content_hash, b.content_hash);
    assert_eq!(a.num_elites, b.num_elites);
    assert!((a.qd_score - b.qd_score).abs() < 1e-12);
    assert_eq!(a.best.gamma.to_bits(), b.best.gamma.to_bits());
}

#[test]
fn malformed_conservation_tolerance_cannot_mint_screened_elites() {
    // Keep the public-path witness cheap while forcing every design through a
    // full simulation. Before the fail-closed guard, +infinity authorized every
    // finite conservation ratio and produced false screened archive entries.
    let budget = CampaignBudget {
        full_steps: 2,
        short_steps: 1,
        bins: 2,
        alpha: 0.1,
        decision_tol: 0.0,
        conserve_tol: f64::INFINITY,
        ..CampaignBudget::default()
    };
    let report = run_campaign(&budget);
    assert!(report.full_sims > 0);
    assert_eq!(report.conservation_screened_elites, 0);
    assert_eq!(report.unscreened_elites, report.num_elites);
    assert_eq!(report.conservation_screened_drift_hull, None);
    assert_eq!(report.campaign_rank, ColorRank::Estimated);
}
