//! Battery for discrepancy probes + the budget pie (addendum Proposal 3).
//! Covers: adjacent-rung probes produce localized estimated-color fields;
//! a near-zero gap manufactures no error; dim-mismatch and top-rung errors;
//! the budget pie attributes a model-form-dominated case to the closure (not
//! the mesh) and vice versa; the probe-budget cap is a HARD ceiling at the
//! exact limit; and determinism.

use fs_evidence::{Color, ColorRank, ValidityDomain};
use fs_ladder::{Ladder, Refine1d};
use fs_probe::{BudgetPie, ErrorContribution, ProbeBudget, ProbeError, probe_adjacent};

fn ladder() -> Ladder {
    Ladder::new("cht", "correlation", 1.0, "bottom").then(
        Box::new(Refine1d),
        "RANS",
        40.0,
        "steady RANS",
    )
}

fn dispersion_of(c: &Color) -> f64 {
    match c {
        Color::Estimated { dispersion, .. } => *dispersion,
        other => panic!("expected Estimated color, got {other:?}"),
    }
}

#[test]
fn near_zero_gap_manufactures_no_error() {
    // fine == prolongate(coarse) → zero discrepancy, estimated-color, zero dispersion.
    let l = ladder();
    let coarse = vec![0.0, 2.0, 4.0];
    let fine = vec![0.0, 1.0, 2.0, 3.0, 4.0]; // == Refine1d prolongation of coarse
    let field = probe_adjacent(&l, 0, &coarse, &fine).unwrap();
    assert_eq!(field.from_rung, 0);
    assert_eq!(field.to_rung, 1);
    assert_eq!(field.l_inf.to_bits(), 0.0f64.to_bits());
    assert_eq!(dispersion_of(&field.color).to_bits(), 0.0f64.to_bits());
    assert_eq!(field.color.rank(), ColorRank::Estimated);
}

#[test]
fn localized_discrepancy_is_measured_and_estimated_color() {
    let l = ladder();
    let coarse = vec![0.0, 2.0, 4.0];
    // fine deviates only at the last point by 6.0.
    let fine = vec![0.0, 1.0, 2.0, 3.0, 10.0];
    let field = probe_adjacent(&l, 0, &coarse, &fine).unwrap();
    assert_eq!(field.per_subdomain.len(), 5);
    assert_eq!(field.l_inf.to_bits(), 6.0f64.to_bits());
    // localized: only the last subdomain carries the error.
    assert_eq!(field.per_subdomain[4].to_bits(), 6.0f64.to_bits());
    assert!(field.per_subdomain[..4].iter().all(|&d| d == 0.0));
    // estimated color names the probe.
    match &field.color {
        Color::Estimated {
            estimator,
            dispersion,
        } => {
            assert!(
                estimator.contains("adjacent-rung-probe:cht:0->1"),
                "{estimator}"
            );
            assert_eq!(dispersion.to_bits(), 6.0f64.to_bits());
        }
        other => panic!("expected Estimated, got {other:?}"),
    }
}

#[test]
fn dim_mismatch_is_a_structured_error() {
    let l = ladder();
    let coarse = vec![0.0, 2.0, 4.0]; // prolongates to 5 points
    let fine = vec![0.0, 1.0, 2.0]; // wrong length
    match probe_adjacent(&l, 0, &coarse, &fine) {
        Err(ProbeError::DimMismatch {
            prolongated_len,
            fine_len,
        }) => {
            assert_eq!(prolongated_len, 5);
            assert_eq!(fine_len, 3);
        }
        other => panic!("expected DimMismatch, got {other:?}"),
    }
}

#[test]
fn probe_at_the_top_rung_bubbles_a_ladder_error() {
    let l = ladder(); // rungs 0,1; rung 1 is the top
    let err = probe_adjacent(&l, 1, &[1.0, 2.0], &[1.0]).unwrap_err();
    assert!(matches!(err, ProbeError::Ladder(_)));
}

#[test]
fn probe_is_deterministic() {
    let l = ladder();
    let coarse = vec![0.5, 1.5, 2.5];
    let fine = vec![0.0, 1.0, 2.0, 3.0, 4.0];
    assert_eq!(
        probe_adjacent(&l, 0, &coarse, &fine).unwrap(),
        probe_adjacent(&l, 0, &coarse, &fine).unwrap()
    );
}

// ---- budget pie -----------------------------------------------------------

fn verified(mag: f64) -> ErrorContribution {
    ErrorContribution::new("mesh", Color::Verified { lo: -mag, hi: mag }, mag)
}
fn estimated(mag: f64) -> ErrorContribution {
    ErrorContribution::new(
        "turbulence-closure",
        Color::Estimated {
            estimator: "rung-probe".into(),
            dispersion: mag,
        },
        mag,
    )
}
fn validated(mag: f64) -> ErrorContribution {
    ErrorContribution::new(
        "wind-tunnel-anchor",
        Color::Validated {
            regime: ValidityDomain::default(),
            dataset: "wt-2026".into(),
        },
        mag,
    )
}

#[test]
fn budget_pie_points_at_the_closure_when_model_form_dominates() {
    // small numerical error, LARGE model-form error → the pie must say
    // "fix the closure", not "refine the mesh".
    let pie = BudgetPie::of(&[verified(0.7), estimated(12.0), validated(1.0)]);
    assert_eq!(pie.dominant(), Some(ColorRank::Estimated));
    assert!(pie.verdict().contains("MODEL-FORM"));
    assert!(!pie.verdict().contains("refine the mesh or raise"));
    // fractions sum to ~1.
    let sum = pie.fraction(ColorRank::Verified)
        + pie.fraction(ColorRank::Validated)
        + pie.fraction(ColorRank::Estimated);
    assert!((sum - 1.0).abs() < 1e-12);
}

#[test]
fn budget_pie_points_at_the_mesh_when_numerical_dominates() {
    let pie = BudgetPie::of(&[verified(20.0), estimated(0.5)]);
    assert_eq!(pie.dominant(), Some(ColorRank::Verified));
    assert!(pie.verdict().contains("refine the mesh"));
}

#[test]
fn empty_budget_has_no_dominant() {
    let pie = BudgetPie::of(&[]);
    assert_eq!(pie.total.to_bits(), 0.0f64.to_bits());
    assert_eq!(pie.dominant(), None);
    assert_eq!(pie.verdict(), "no error budget recorded");
}

#[test]
fn ties_resolve_conservatively_to_the_weaker_color() {
    // equal verified and estimated magnitudes → report the WEAKER (estimated).
    let pie = BudgetPie::of(&[verified(5.0), estimated(5.0)]);
    assert_eq!(pie.dominant(), Some(ColorRank::Estimated));
}

// ---- probe budget cap -----------------------------------------------------

#[test]
fn probe_budget_cap_is_a_hard_ceiling_at_the_exact_limit() {
    // fleet 100, cap 10% → cap = 10.0.
    let mut b = ProbeBudget::new(100.0, 0.10);
    assert_eq!(b.cap().to_bits(), 10.0f64.to_bits());
    // spend up to EXACTLY the cap: allowed.
    assert!(b.try_spend(7.0).is_ok());
    assert!(b.try_spend(3.0).is_ok()); // now at exactly 10.0
    assert_eq!(b.remaining().to_bits(), 0.0f64.to_bits());
    // one iota beyond the cap: refused.
    match b.try_spend(0.001) {
        Err(ProbeError::BudgetExceeded { .. }) => {}
        other => panic!("expected BudgetExceeded, got {other:?}"),
    }
    // and the refused spend did not mutate the budget.
    assert_eq!(b.spent().to_bits(), 10.0f64.to_bits());
}

#[test]
fn probe_budget_rejects_bad_costs() {
    let mut b = ProbeBudget::new(100.0, 0.5);
    assert!(matches!(b.try_spend(-1.0), Err(ProbeError::BadCost { .. })));
    assert!(matches!(
        b.try_spend(f64::NAN),
        Err(ProbeError::BadCost { .. })
    ));
    assert!(matches!(
        b.try_spend(f64::INFINITY),
        Err(ProbeError::BadCost { .. })
    ));
    // clamping: a >1 fraction is clamped to the whole fleet budget.
    let over = ProbeBudget::new(100.0, 2.0);
    assert_eq!(over.cap().to_bits(), 100.0f64.to_bits());
}

#[test]
fn probe_budget_nonfinite_inputs_fail_closed() {
    for mut budget in [
        ProbeBudget::new(f64::NAN, 0.5),
        ProbeBudget::new(f64::INFINITY, 0.5),
        ProbeBudget::new(100.0, f64::NAN),
        ProbeBudget::new(100.0, f64::INFINITY),
    ] {
        assert_eq!(budget.cap().to_bits(), 0.0f64.to_bits());
        assert_eq!(budget.remaining().to_bits(), 0.0f64.to_bits());
        assert!(matches!(
            budget.try_spend(1.0),
            Err(ProbeError::BudgetExceeded { .. })
        ));
        assert_eq!(budget.spent().to_bits(), 0.0f64.to_bits());
    }
}
