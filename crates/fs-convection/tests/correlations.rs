//! G0/G3 correlation-card, refusal, typing, and metamorphic coverage.

use fs_convection::{
    CorrelationError, CorrelationId, CorrelationInputs, DiscrepancyBasis, HeatTransferCoefficient,
    ThermalConductivity, ThermalDirection, correlation_catalog, evaluate,
};
use fs_evidence::{CertifyError, NumericalKind};
use fs_math::det;
use fs_qty::Length;
use fs_vvreg::thermal_level_a::{
    ThermalLevelAAcceptance, ThermalLevelACase, ThermalLevelAFamily, ThermalLevelAKind,
    thermal_level_a_cases,
};

const LEVEL_A_CONVECTION_BINDINGS: [(&str, &str); 2] = [
    (
        "thermal-a-duct-nu-cwt",
        "tests/correlations.rs::level_a_fully_developed_limits_and_rectangular_square_values_are_frozen",
    ),
    (
        "thermal-a-duct-nu-chf",
        "tests/correlations.rs::level_a_fully_developed_limits_and_rectangular_square_values_are_frozen",
    ),
];

fn close(observed: f64, expected: f64, relative: f64) {
    let scale = expected.abs().max(1.0);
    assert!(
        (observed - expected).abs() <= relative * scale,
        "observed={observed:.17e} expected={expected:.17e} tolerance={relative:.3e}"
    );
}

fn duct_inputs() -> CorrelationInputs {
    CorrelationInputs::forced(1_000.0, 0.7).with_length_ratio(100.0)
}

fn level_a_convection_case(case_id: &str) -> &'static ThermalLevelACase {
    let case = thermal_level_a_cases()
        .iter()
        .find(|case| case.id == case_id)
        .unwrap_or_else(|| panic!("missing Level-A case {case_id}"));
    assert_eq!(case.family, ThermalLevelAFamily::ConvectionLimit);
    assert_eq!(case.kind, ThermalLevelAKind::AnalyticReference);
    assert_eq!(case.metric, "nusselt-number");
    assert!(
        LEVEL_A_CONVECTION_BINDINGS
            .iter()
            .any(|(id, _)| *id == case_id),
        "{case_id} is not declared as an executing fs-convection binding"
    );
    let reynolds = case
        .context
        .iter()
        .find(|axis| axis.name == "reynolds-number")
        .unwrap_or_else(|| panic!("{case_id} must declare its Reynolds context"));
    assert!(reynolds.lo <= 1_000.0 && 1_000.0 <= reynolds.hi);
    case
}

fn assert_level_a_convection_limit(case_id: &str, observed: f64) {
    let case = level_a_convection_case(case_id);
    let ThermalLevelAAcceptance::Tolerance { atol, rtol } = case.acceptance else {
        panic!("{case_id}: analytic Level-A row must carry a scalar tolerance");
    };
    let absolute_error = (observed - case.reference_value_si).abs();
    let envelope = atol + rtol * case.reference_value_si.abs();
    assert!(
        absolute_error <= envelope,
        "{case_id}: Nu={observed:.17e}, reference={:.17e}, error={absolute_error:.3e}, \
         envelope={envelope:.3e}",
        case.reference_value_si
    );
    assert_eq!(observed.to_bits(), case.reference_value_si.to_bits());
    println!(
        "{{\"suite\":\"fs-convection/level-a\",\"case_id\":\"{case_id}\",\
         \"computed\":{observed:.17e},\"reference\":{:.17e},\
         \"absolute_error\":{absolute_error:.17e},\"envelope\":{envelope:.17e},\
         \"authority\":\"executed-formula-limit-not-registry-receipt\",\
         \"verdict\":\"pass\"}}",
        case.reference_value_si
    );
}

#[test]
fn catalog_has_eleven_sourced_cards_and_no_unlabeled_discrepancy() {
    let catalog = correlation_catalog();
    assert_eq!(catalog.len(), 11);
    assert_eq!(catalog.len(), CorrelationId::ALL.len());

    let mut names = std::collections::BTreeSet::new();
    for card in catalog {
        assert!(
            names.insert(card.id.name()),
            "duplicate card {}",
            card.id.name()
        );
        assert_eq!(card.model.name, card.id.name());
        assert!(!card.model.validity.bounds().is_empty());
        assert!(!card.source.citation.trim().is_empty());
        assert!(!card.source.identifier.trim().is_empty());
        assert!(!card.model.assumptions.is_empty());
        assert!(!card.model.known_failures.is_empty());
        match card.discrepancy_basis {
            DiscrepancyBasis::AnalyticIdealLimit => {
                assert_eq!(card.model.discrepancy_rel.to_bits(), 0);
            }
            DiscrepancyBasis::EngineeringAllowance => {
                assert!(card.model.discrepancy_rel >= 0.10 && card.model.discrepancy_rel <= 0.25);
            }
        }
    }
}

#[test]
fn flat_plate_cards_retain_distinct_formula_authorities() {
    let catalog = correlation_catalog();
    let laminar = catalog
        .iter()
        .find(|card| card.id == CorrelationId::FlatPlateLaminarAverage)
        .expect("laminar flat-plate card");
    let mixed = catalog
        .iter()
        .find(|card| card.id == CorrelationId::FlatPlateTurbulentAverage)
        .expect("mixed flat-plate card");

    assert_ne!(laminar.source, mixed.source);
    assert!(laminar.source.citation.contains("Pohlhausen"));
    assert_eq!(
        mixed.source.identifier,
        "doi:10.1002/9781119686040; ISBN 978-1-119-68597-5"
    );
    assert!(mixed.source.citation.contains("Eq. (1.62)"));
    assert!(
        mixed
            .model
            .assumptions
            .iter()
            .any(|assumption| assumption.contains("laminar leading edge"))
    );
    assert!(
        mixed
            .model
            .known_failures
            .iter()
            .any(|failure| failure.contains("transition Reynolds number"))
    );
}

#[test]
fn level_a_fully_developed_limits_and_rectangular_square_values_are_frozen() {
    let cwt = evaluate(CorrelationId::CircularDuctLaminarCwt, duct_inputs()).expect("CWT");
    let chf = evaluate(CorrelationId::CircularDuctLaminarChf, duct_inputs()).expect("CHF");
    assert_level_a_convection_limit("thermal-a-duct-nu-cwt", cwt.evidence().value);
    assert_level_a_convection_limit("thermal-a-duct-nu-chf", chf.evidence().value);

    let square = CorrelationInputs::forced(1_000.0, 0.7)
        .with_length_ratio(100.0)
        .with_aspect_ratio(1.0);
    close(
        evaluate(CorrelationId::RectangularDuctLaminarCwt, square)
            .expect("square CWT")
            .evidence()
            .value,
        2.978_695,
        2.0e-15,
    );
    close(
        evaluate(CorrelationId::RectangularDuctLaminarChf, square)
            .expect("square CHF")
            .evidence()
            .value,
        3.610_224,
        2.0e-15,
    );
}

#[test]
fn level_a_convection_binding_partition_is_complete() {
    let catalog_ids = thermal_level_a_cases()
        .iter()
        .filter(|case| case.family == ThermalLevelAFamily::ConvectionLimit)
        .map(|case| case.id)
        .collect::<std::collections::BTreeSet<_>>();
    let binding_ids = LEVEL_A_CONVECTION_BINDINGS
        .iter()
        .map(|(id, _)| *id)
        .collect::<std::collections::BTreeSet<_>>();
    assert_eq!(binding_ids, catalog_ids);
    assert_eq!(binding_ids.len(), 2);
    for (id, test) in LEVEL_A_CONVECTION_BINDINGS {
        assert!(
            test.starts_with("tests/correlations.rs::"),
            "{id}: executing test path must be stable"
        );
    }
}

#[test]
fn every_nonconstant_formula_has_a_frozen_source_formula_spot_value() {
    let cases = [
        (
            CorrelationId::CircularDuctHausen,
            CorrelationInputs::forced(1_000.0, 7.0).with_length_ratio(20.0),
            11.488_360_610_697_356,
        ),
        (
            CorrelationId::RectangularDuctLaminarCwt,
            CorrelationInputs::forced(1_000.0, 0.7)
                .with_length_ratio(100.0)
                .with_aspect_ratio(0.5),
            3.388_736_875_000_000_6,
        ),
        (
            CorrelationId::RectangularDuctLaminarChf,
            CorrelationInputs::forced(1_000.0, 0.7)
                .with_length_ratio(100.0)
                .with_aspect_ratio(0.5),
            4.125_812_203_124_999,
        ),
        (
            CorrelationId::DittusBoelter,
            CorrelationInputs::forced(100_000.0, 0.7).with_length_ratio(100.0),
            199.419_237_807_658_48,
        ),
        (
            CorrelationId::Gnielinski,
            CorrelationInputs::forced(100_000.0, 0.7).with_length_ratio(100.0),
            178.622_951_779_291_2,
        ),
        (
            CorrelationId::FlatPlateLaminarAverage,
            CorrelationInputs::forced(100_000.0, 0.7),
            186.437_852_875_226_2,
        ),
        (
            CorrelationId::FlatPlateTurbulentAverage,
            CorrelationInputs::forced(1_000_000.0, 0.7),
            1_299.484_953_525_734_2,
        ),
        (
            CorrelationId::ChurchillBernsteinCylinder,
            CorrelationInputs::forced(10_000.0, 0.7),
            53.327_788_670_209_97,
        ),
        (
            CorrelationId::ChurchillChuVerticalPlate,
            CorrelationInputs::natural(1.0e6, 0.7),
            16.530_366_876_407_225,
        ),
    ];

    for (id, inputs, expected) in cases {
        let observed = evaluate(id, inputs).unwrap_or_else(|error| panic!("{id:?}: {error}"));
        close(observed.evidence().value, expected, 3.0e-13);
    }
}

#[test]
fn validity_edges_are_inclusive_and_missing_or_outside_axes_refuse() {
    evaluate(
        CorrelationId::DittusBoelter,
        CorrelationInputs::forced(10_000.0, 0.7).with_length_ratio(10.0),
    )
    .expect("inclusive lower edge");
    evaluate(
        CorrelationId::DittusBoelter,
        CorrelationInputs::forced(120_000.0, 120.0).with_length_ratio(1.0e6),
    )
    .expect("inclusive upper edge");

    let missing = evaluate(
        CorrelationId::DittusBoelter,
        CorrelationInputs::forced(100_000.0, 0.7),
    )
    .expect_err("L/Dh is mandatory");
    match missing {
        CorrelationError::OutOfDomain { violations, .. } => {
            assert_eq!(violations.len(), 1);
            assert_eq!(violations[0].axis, "L_over_Dh");
            assert_eq!(violations[0].value, None);
        }
        other => panic!("unexpected refusal: {other}"),
    }

    let outside = evaluate(
        CorrelationId::DittusBoelter,
        CorrelationInputs::forced(9_999.0, 0.69).with_length_ratio(9.0),
    )
    .expect_err("three axes are outside");
    match outside {
        CorrelationError::OutOfDomain { violations, .. } => {
            let axes = violations
                .iter()
                .map(|violation| violation.axis.as_str())
                .collect::<Vec<_>>();
            assert_eq!(axes, ["L_over_Dh", "Pr", "Re"]);
        }
        other => panic!("unexpected refusal: {other}"),
    }

    assert!(matches!(
        evaluate(
            CorrelationId::FlatPlateLaminarAverage,
            CorrelationInputs::forced(f64::NAN, 0.7)
        ),
        Err(CorrelationError::InvalidGroup { axis: "Re", .. })
    ));
}

#[test]
fn cylinder_card_checks_the_product_constraint_through_a_named_peclet_axis() {
    let accepted = evaluate(
        CorrelationId::ChurchillBernsteinCylinder,
        CorrelationInputs::forced(1.0, 0.2),
    )
    .expect("Pe=0.2 is inclusive");
    assert_eq!(accepted.groups().get("Pe"), Some(&0.2));

    let refused = evaluate(
        CorrelationId::ChurchillBernsteinCylinder,
        CorrelationInputs::forced(1.0, 0.199),
    )
    .expect_err("Pe and Pr are outside");
    let CorrelationError::OutOfDomain { violations, .. } = refused else {
        panic!("expected domain refusal");
    };
    assert!(violations.iter().any(|violation| violation.axis == "Pe"));
    assert!(violations.iter().any(|violation| violation.axis == "Pr"));
}

#[test]
fn typed_h_retains_model_evidence_and_empirical_result_cannot_certify() {
    fn require_htc(_: &fs_evidence::Evidence<HeatTransferCoefficient>) {}

    let evaluated = evaluate(
        CorrelationId::Gnielinski,
        CorrelationInputs::forced(100_000.0, 0.7).with_length_ratio(100.0),
    )
    .expect("in domain");
    let coefficient = evaluated
        .heat_transfer_coefficient(ThermalConductivity::new(0.026), Length::new(0.010))
        .expect("typed h");
    require_htc(&coefficient);
    close(
        coefficient.value.value(),
        evaluated.evidence().value * 0.026 / 0.010,
        1.0e-15,
    );
    assert_eq!(coefficient.numerical.kind, NumericalKind::Estimate);
    assert_eq!(coefficient.model.cards, [CorrelationId::Gnielinski.name()]);
    assert!(coefficient.model.in_domain);
    assert_eq!(coefficient.model.discrepancy_rel, 0.10);
    assert!(matches!(
        coefficient.certified(),
        Err(CertifyError::NotRigorous {
            kind: NumericalKind::Estimate
        })
    ));
}

#[test]
fn g3_coherent_unit_rescaling_leaves_nu_and_h_invariant() {
    let evaluated = evaluate(
        CorrelationId::FlatPlateLaminarAverage,
        CorrelationInputs::forced(100_000.0, 0.7),
    )
    .expect("in domain");
    let si = evaluated
        .heat_transfer_coefficient(ThermalConductivity::new(0.026), Length::new(0.1))
        .expect("SI");
    // The same k supplied as 26 mW/(m K), and L as 100 mm, normalized
    // explicitly at the boundary before entering the coherent-SI types.
    let rescaled = evaluated
        .heat_transfer_coefficient(
            ThermalConductivity::new(26.0 * 1.0e-3),
            Length::new(100.0 * 1.0e-3),
        )
        .expect("rescaled");
    close(si.value.value(), rescaled.value.value(), 1.0e-15);
    assert_eq!(si.model, rescaled.model);
}

#[test]
fn dittus_boelter_direction_is_semantic_and_provenance_bearing() {
    let heating = evaluate(
        CorrelationId::DittusBoelter,
        CorrelationInputs::forced(100_000.0, 0.7)
            .with_length_ratio(100.0)
            .with_direction(ThermalDirection::HeatingFluid),
    )
    .expect("heating");
    let cooling = evaluate(
        CorrelationId::DittusBoelter,
        CorrelationInputs::forced(100_000.0, 0.7)
            .with_length_ratio(100.0)
            .with_direction(ThermalDirection::CoolingFluid),
    )
    .expect("cooling");
    close(cooling.evidence().value, 206.660_391_611_847_25, 3.0e-13);
    assert_ne!(
        heating.evidence().value.to_bits(),
        cooling.evidence().value.to_bits()
    );
    assert_ne!(heating.evidence().provenance, cooling.evidence().provenance);
}

#[test]
fn dimensional_inputs_refuse_zero_negative_nan_and_overflow() {
    let evaluated = evaluate(
        CorrelationId::FlatPlateLaminarAverage,
        CorrelationInputs::forced(100_000.0, 0.7),
    )
    .expect("in domain");
    for conductivity in [0.0, -1.0, f64::NAN] {
        assert!(matches!(
            evaluated.heat_transfer_coefficient(
                ThermalConductivity::new(conductivity),
                Length::new(0.1)
            ),
            Err(CorrelationError::InvalidDimensionalInput {
                field: "fluid thermal conductivity",
                ..
            })
        ));
    }
    for length in [0.0, -1.0, f64::INFINITY] {
        assert!(matches!(
            evaluated
                .heat_transfer_coefficient(ThermalConductivity::new(0.026), Length::new(length)),
            Err(CorrelationError::InvalidDimensionalInput {
                field: "characteristic length",
                ..
            })
        ));
    }
    assert!(matches!(
        evaluated.heat_transfer_coefficient(
            ThermalConductivity::new(f64::MAX),
            Length::new(f64::MIN_POSITIVE)
        ),
        Err(CorrelationError::NonFiniteResult {
            stage: "Nu-to-h conversion",
            ..
        })
    ));
}

// ---------------------------------------------------------------------------
// Limiting-behavior checks for the NONCONSTANT cards.
//
// The Level-A rows above freeze the two constant-Nu duct values. Those cards
// return their literal unconditionally, so a value comparison against the
// registry constant exercises no arithmetic: it proves the domain gate admits
// the point and that nobody typo'd a constant, and nothing more. The tests
// below are the limiting-behavior half. Each one asserts a relation that is
// derivable from the cards' own published coefficients rather than from this
// implementation's output, so a drifted coefficient fails them.
// ---------------------------------------------------------------------------

/// Graetz-number ladder shared by [`CorrelationId::CircularDuctHausen`] and
/// [`CorrelationId::CircularDuctLaminarCwt`].
///
/// Hausen admits `L_over_Dh` in `[0.05, 1000]`; the fully-developed card admits
/// `[50, 1e6]`. Their Reynolds and Prandtl boxes are identical, so on the
/// intersection `[50, 1000]` both cards evaluate the SAME point with no
/// refusal. The ratios below are exact successive halvings of `Gz`.
const SHARED_DUCT_LENGTH_RATIOS: [f64; 5] = [62.5, 125.0, 250.0, 500.0, 1_000.0];

#[test]
fn hausen_converges_to_the_fully_developed_duct_limit_as_graetz_vanishes() {
    // Hausen is the thermally developing relation
    //     Nu = 3.66 + 0.0668 Gz / (1 + 0.04 Gz^(2/3)),   Gz = Re Pr / (L/Dh),
    // whose Gz -> 0 limit is precisely the fully-developed constant-wall-
    // temperature value returned by CircularDuctLaminarCwt. That makes the
    // constant an actual LIMIT of a nonconstant card rather than a frozen
    // literal.
    //
    // Gz = 0 is not admissible -- insert_group refuses non-positive groups and
    // L_over_Dh is bounded above -- so the limit is established by the
    // approach: a strictly positive gap that shrinks under the analytic
    // bracket and halves when Gz halves.
    const REYNOLDS: f64 = 1.0;
    const PRANDTL: f64 = 0.6;

    let mut previous_gap = f64::INFINITY;
    let mut previous_ratio_checked = false;
    for length_ratio in SHARED_DUCT_LENGTH_RATIOS {
        let point = CorrelationInputs::forced(REYNOLDS, PRANDTL).with_length_ratio(length_ratio);
        let developed = evaluate(CorrelationId::CircularDuctLaminarCwt, point)
            .expect("fully-developed card admits the shared point")
            .evidence()
            .value;
        let developing = evaluate(CorrelationId::CircularDuctHausen, point)
            .expect("Hausen card admits the shared point")
            .evidence()
            .value;

        let graetz = REYNOLDS * PRANDTL / length_ratio;
        let gap = developing - developed;

        // Entry-length enhancement is strictly positive: the developing
        // solution never falls below its own fully-developed limit.
        assert!(
            gap > 0.0,
            "L/Dh={length_ratio}: developing Nu={developing:.17e} must exceed \
             the fully-developed limit {developed:.17e}"
        );
        // Since 1 + 0.04 Gz^(2/3) >= 1 for every admissible Gz, the enhancement
        // is bounded above by the bare numerator. This bracket comes from the
        // published form, not from this implementation's output.
        assert!(
            gap <= 0.0668 * graetz,
            "L/Dh={length_ratio}: gap={gap:.9e} exceeds the analytic bracket \
             {:.9e}",
            0.0668 * graetz
        );
        // Monotone approach to the limit.
        assert!(
            gap < previous_gap,
            "L/Dh={length_ratio}: gap={gap:.9e} did not shrink below \
             {previous_gap:.9e}"
        );
        // First-order approach: halving Gz halves the residual gap. This is the
        // quantitative content of "the constant is the limit" -- a card that
        // merely stayed near 3.66 without converging would fail here.
        if previous_gap.is_finite() {
            let ratio = previous_gap / gap;
            assert!(
                (1.99..=2.01).contains(&ratio),
                "L/Dh={length_ratio}: gap ratio {ratio:.6} is not the \
                 first-order value 2 for a halved Graetz number"
            );
            previous_ratio_checked = true;
        }
        previous_gap = gap;
    }
    assert!(
        previous_ratio_checked,
        "the ladder must contain at least one halving"
    );
    // At the smallest admissible Graetz number the developing card has settled
    // onto the fully-developed constant to ~1.1e-5 relative.
    assert!(
        previous_gap / 3.66 < 2.0e-5,
        "residual relative gap {:.6e} at the smallest admissible Graetz number",
        previous_gap / 3.66
    );
}

#[test]
fn flat_plate_cards_are_continuous_across_their_shared_transition_reynolds_number() {
    // Re = 5e5 is simultaneously the INCLUSIVE upper bound of the laminar
    // average card and the INCLUSIVE lower bound of the mixed
    // laminar-to-turbulent card, so both admit exactly that point.
    //
    // The mixed card's leading-edge constant is not free: it is the offset that
    // makes the turbulent form reproduce the laminar average at transition,
    //     0.037 Re_c^0.8 - 0.664 Re_c^0.5 = 871.3234750958699,
    // rounded to 871 in the published relation. So the two cards must agree at
    // Re_c to within that rounding, and -- because both carry the same
    // Pr^(1/3) factor -- the relative disagreement must be EXACTLY independent
    // of Pr.
    const TRANSITION_REYNOLDS: f64 = 5.0e5;
    // Prandtl window shared by both cards: laminar [0.6, 50], mixed [0.6, 60].
    const SHARED_PRANDTL: [f64; 5] = [0.6, 0.7, 1.0, 10.0, 50.0];

    let mut relative_gaps = Vec::new();
    for prandtl in SHARED_PRANDTL {
        let point = CorrelationInputs::forced(TRANSITION_REYNOLDS, prandtl);
        let laminar = evaluate(CorrelationId::FlatPlateLaminarAverage, point)
            .expect("laminar card admits its inclusive upper Reynolds bound")
            .evidence()
            .value;
        let mixed = evaluate(CorrelationId::FlatPlateTurbulentAverage, point)
            .expect("mixed card admits its inclusive lower Reynolds bound")
            .evidence()
            .value;
        // The mixed card subtracts a constant; positivity at the transition
        // point is the boundary of where that subtraction stays physical. Its
        // zero crossing, Re = (871/0.037)^1.25 = 291588.6, sits below the
        // card's floor, so it can never return a non-positive Nu in-domain.
        assert!(
            mixed > 0.0,
            "Pr={prandtl}: mixed card returned a non-positive Nu at transition"
        );
        relative_gaps.push((mixed - laminar) / laminar);
    }

    // The two independent relations agree at transition to better than 0.1%.
    for (prandtl, gap) in SHARED_PRANDTL.iter().zip(&relative_gaps) {
        assert!(
            gap.abs() < 1.0e-3,
            "Pr={prandtl}: transition discontinuity {gap:.6e} exceeds 0.1%"
        );
    }
    // And the disagreement is Pr-independent, because the shared Pr^(1/3)
    // factor cancels in the ratio. A Prandtl exponent that drifted apart
    // between the two cards would break this even if each card stayed
    // plausible on its own.
    for (prandtl, gap) in SHARED_PRANDTL.iter().zip(&relative_gaps) {
        close(*gap, relative_gaps[0], 1.0e-12);
        assert!(gap.is_finite(), "Pr={prandtl}: non-finite relative gap");
    }

    // Pin the rounded constant itself. At Pr = 1 the Pr^(1/3) factors are
    // exactly one, so the raw difference between the cards is
    //     (0.037 Re_c^0.8 - 871) - 0.664 Re_c^0.5 = 871.3234750958699 - 871,
    // i.e. it measures the rounding residue of the published 871 directly.
    // Substituting 872 would move this to -0.6765: a sign flip.
    let unit_prandtl = CorrelationInputs::forced(TRANSITION_REYNOLDS, 1.0);
    let laminar = evaluate(CorrelationId::FlatPlateLaminarAverage, unit_prandtl)
        .expect("laminar at unit Prandtl")
        .evidence()
        .value;
    let mixed = evaluate(CorrelationId::FlatPlateTurbulentAverage, unit_prandtl)
        .expect("mixed at unit Prandtl")
        .evidence()
        .value;
    close(mixed - laminar, 0.323_475_095_869_866_9, 1.0e-9);
}

#[test]
fn dittus_boelter_direction_exponents_degenerate_at_unit_prandtl() {
    // Pr = 1 lies inside the card's [0.7, 120] window. The heating and cooling
    // forms differ ONLY in the Prandtl exponent (0.4 versus 0.3), and
    // 1^0.4 = 1^0.3 = 1 exactly, so at that point the two directions must
    // return bit-identical Nusselt numbers equal to the bare 0.023 Re^0.8
    // factor. The existing direction test evaluates at Pr = 0.7, where it can
    // only assert that the two differ; this pins the point where they must not.
    for reynolds in [10_000.0, 50_000.0, 120_000.0] {
        let heating = evaluate(
            CorrelationId::DittusBoelter,
            CorrelationInputs::forced(reynolds, 1.0)
                .with_length_ratio(100.0)
                .with_direction(ThermalDirection::HeatingFluid),
        )
        .expect("heating at unit Prandtl");
        let cooling = evaluate(
            CorrelationId::DittusBoelter,
            CorrelationInputs::forced(reynolds, 1.0)
                .with_length_ratio(100.0)
                .with_direction(ThermalDirection::CoolingFluid),
        )
        .expect("cooling at unit Prandtl");

        assert_eq!(
            heating.evidence().value.to_bits(),
            cooling.evidence().value.to_bits(),
            "Re={reynolds}: the direction exponents must degenerate at Pr = 1"
        );
        // The declared direction is still recorded even where it cannot move
        // the number: provenance tracks the input, not merely the output.
        assert_ne!(
            heating.evidence().provenance,
            cooling.evidence().provenance,
            "Re={reynolds}: provenance must retain the declared direction"
        );
    }
}

#[test]
fn gnielinski_correction_denominator_collapses_at_unit_prandtl() {
    // At Pr = 1 the correction denominator
    //     1 + 12.7 sqrt(f/8) (Pr^(2/3) - 1)
    // has a vanishing second term, so the relation reduces exactly to
    //     Nu = (f/8)(Re - 1000),   f = (0.79 ln Re - 1.64)^-2,
    // with no Prandtl dependence at all. Pr = 1 is inside the card's
    // [0.5, 2000] window, so this reduction is reachable through evaluate().
    // The reduced form is rebuilt here from the published coefficients using
    // the same deterministic primitives the crate evaluates with, so the
    // comparison isolates the correction structure rather than fs-math's
    // rounding.
    for reynolds in [3_000.0, 10_000.0, 100_000.0, 5.0e6] {
        let observed = evaluate(
            CorrelationId::Gnielinski,
            CorrelationInputs::forced(reynolds, 1.0).with_length_ratio(100.0),
        )
        .expect("Gnielinski at unit Prandtl")
        .evidence()
        .value;

        let friction = 1.0 / det::pow(0.79_f64.mul_add(det::ln(reynolds), -1.64), 2.0);
        let reduced = (friction / 8.0) * (reynolds - 1_000.0);

        assert_eq!(
            observed.to_bits(),
            reduced.to_bits(),
            "Re={reynolds}: Nu={observed:.17e} must collapse to the \
             Prandtl-free reduced form {reduced:.17e} at Pr = 1"
        );
    }
}

#[test]
fn churchill_chu_recovers_its_inadmissible_zero_rayleigh_intercept() {
    // The relation is
    //     Nu = (0.825 + 0.387 Ra^(1/6) / D(Pr))^2,
    // so sqrt(Nu) is AFFINE in x = Ra^(1/6) with intercept exactly 0.825 --
    // the Ra -> 0 pure-conduction limit.
    //
    // That endpoint is doubly inadmissible: insert_group refuses non-positive
    // groups, and the card floor is Ra >= 0.1. It therefore cannot be reached
    // through evaluate() at any tolerance. But two IN-DOMAIN points determine
    // the line, and its intercept recovers the endpoint exactly with no
    // refusal. This is the technique for the crate's several inadmissible
    // endpoints; without it they can only be read off the source.
    //
    // The Rayleigh values are exact sixth powers, so the abscissae 1 and 100
    // are exact in the test and introduce no rounding of their own.
    const PRANDTL: f64 = 0.7;
    const LOW_RAYLEIGH: f64 = 1.0;
    const HIGH_RAYLEIGH: f64 = 1.0e12;
    const LOW_ABSCISSA: f64 = 1.0; // LOW_RAYLEIGH^(1/6)
    const HIGH_ABSCISSA: f64 = 100.0; // HIGH_RAYLEIGH^(1/6)

    let root_nusselt = |rayleigh: f64| {
        det::sqrt(
            evaluate(
                CorrelationId::ChurchillChuVerticalPlate,
                CorrelationInputs::natural(rayleigh, PRANDTL),
            )
            .expect("Churchill-Chu admits the in-domain Rayleigh point")
            .evidence()
            .value,
        )
    };

    let high = root_nusselt(HIGH_RAYLEIGH);
    let low = root_nusselt(LOW_RAYLEIGH);
    let slope = (high - low) / (HIGH_ABSCISSA - LOW_ABSCISSA);
    let intercept = low - slope * LOW_ABSCISSA;

    // A drifted leading constant moves the intercept one-for-one, so this
    // recovers 0.825 to far better than any plausible perturbation.
    close(intercept, 0.825, 1.0e-9);
    // The fitted slope is 0.387/D(Pr) > 0: the relation is genuinely
    // increasing in Ra^(1/6) rather than accidentally flat.
    assert!(
        slope > 0.0,
        "fitted slope {slope:.17e} must be positive over 12 decades of Rayleigh"
    );
}
