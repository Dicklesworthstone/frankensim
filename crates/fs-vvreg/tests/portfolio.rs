//! G0/G3 tests for claim-scoped external-evidence portfolios.

use fs_vvreg::ContentHash;
use fs_vvreg::corpus::EvidenceLevel;
use fs_vvreg::portfolio::{
    EvidenceAxis, EvidencePortfolio, PortfolioClaimClass, PortfolioObservation, PortfolioRefusal,
    axes_for_level,
};

const QOI: &str = "component-peak-temperature";
const REGIME: &str = "cooling-vertical-v1";
const DOCTRINE_DOC: &str = include_str!("../../../docs/EVIDENCE_PORTFOLIO_DOCTRINE.md");
const VVREG_CONTRACT: &str = include_str!("../CONTRACT.md");
const EVIDENCE_CONTRACT: &str = include_str!("../../fs-evidence/CONTRACT.md");

fn hash(byte: u8) -> ContentHash {
    ContentHash([byte; 32])
}

fn observation(axis: EvidenceAxis, source: u8, independence_group: u8) -> PortfolioObservation {
    PortfolioObservation::try_new(axis, QOI, REGIME, hash(source), hash(independence_group))
        .expect("valid portfolio observation")
}

#[test]
fn doctrine_is_present_and_contract_referenced() {
    for heading in [
        "# External Evidence Is a Portfolio, Not a Pyramid",
        "## The seven coordinates",
        "## A-E corpus migration",
        "## Claim admission",
        "## Anti-laundering rules",
        "## Scorecard behavior",
        "### Worked refusal: field data are not a lab surrogate",
        "## No-claim boundaries",
    ] {
        assert!(
            DOCTRINE_DOC.contains(heading),
            "doctrine is missing required heading {heading}"
        );
    }
    for (axis, doctrine_name) in [
        (
            EvidenceAxis::NumericalVerification,
            "Numerical verification",
        ),
        (EvidenceAxis::CrossCodeAgreement, "Cross-code agreement"),
        (
            EvidenceAxis::ControlledExperimentalValidation,
            "Controlled experimental validation",
        ),
        (
            EvidenceAxis::BlindPredictiveValidation,
            "Blind predictive validation",
        ),
        (EvidenceAxis::FieldMonitoring, "Field monitoring"),
        (
            EvidenceAxis::TransferabilityAcrossRegimes,
            "Transferability across regimes",
        ),
        (
            EvidenceAxis::IndependentReproduction,
            "Independent reproduction",
        ),
    ] {
        assert!(
            DOCTRINE_DOC.contains(doctrine_name),
            "doctrine is missing axis {}",
            axis.slug()
        );
    }
    assert!(DOCTRINE_DOC.contains("Field monitoring alone cannot support `Validated`"));
    assert!(DOCTRINE_DOC.contains("missing axis controlled-experimental-validation"));
    assert!(VVREG_CONTRACT.contains("portfolio::{EvidenceAxis,EvidencePortfolio"));
    assert!(EVIDENCE_CONTRACT.contains("`fs-vvreg::portfolio`"));
}

#[test]
fn legacy_a_to_e_tags_map_to_coordinates_not_ranks() {
    assert_eq!(
        axes_for_level(EvidenceLevel::Analytic),
        &[EvidenceAxis::NumericalVerification]
    );
    assert_eq!(
        axes_for_level(EvidenceLevel::CrossCode),
        &[EvidenceAxis::CrossCodeAgreement]
    );
    assert_eq!(
        axes_for_level(EvidenceLevel::PublishedExperiment),
        &[EvidenceAxis::ControlledExperimentalValidation]
    );
    assert_eq!(
        axes_for_level(EvidenceLevel::Blind),
        &[
            EvidenceAxis::ControlledExperimentalValidation,
            EvidenceAxis::BlindPredictiveValidation,
        ]
    );
    assert_eq!(
        axes_for_level(EvidenceLevel::Field),
        &[EvidenceAxis::FieldMonitoring]
    );
}

#[test]
fn every_claim_class_names_each_missing_axis() {
    let cases = [
        PortfolioClaimClass::NumericallyVerified,
        PortfolioClaimClass::CrossCodeConsistent,
        PortfolioClaimClass::ValidatedPrediction,
        PortfolioClaimClass::BlindValidatedPrediction,
        PortfolioClaimClass::FieldSupportedPrediction,
        PortfolioClaimClass::TransferablePrediction,
        PortfolioClaimClass::IndependentlyReproducedPrediction,
    ];

    for claim_class in cases {
        for &missing in claim_class.required_axes() {
            let observations = claim_class
                .required_axes()
                .iter()
                .copied()
                .filter(|axis| *axis != missing)
                .enumerate()
                .map(|(index, axis)| observation(axis, index as u8 + 1, index as u8 + 20))
                .collect();
            let portfolio = EvidencePortfolio::try_new(observations).unwrap();
            assert!(matches!(
                portfolio.admit(claim_class, QOI, REGIME),
                Err(PortfolioRefusal::MissingAxis {
                    claim_class: got_class,
                    axis: got_axis,
                    ref qoi,
                    ref regime,
                }) if got_class == claim_class
                    && got_axis == missing
                    && qoi == QOI
                    && regime == REGIME
            ));
        }
    }
}

#[test]
fn field_evidence_alone_cannot_launder_a_validated_prediction() {
    let repeated_field = observation(EvidenceAxis::FieldMonitoring, 1, 7);
    let portfolio = EvidencePortfolio::try_new(vec![
        repeated_field.clone(),
        repeated_field.clone(),
        repeated_field,
    ])
    .unwrap();
    assert_eq!(portfolio.observations().len(), 1);
    assert!(matches!(
        portfolio.admit(PortfolioClaimClass::ValidatedPrediction, QOI, REGIME),
        Err(PortfolioRefusal::MissingAxis {
            axis: EvidenceAxis::ControlledExperimentalValidation,
            ..
        })
    ));
    assert!(matches!(
        portfolio.admit(PortfolioClaimClass::FieldSupportedPrediction, QOI, REGIME),
        Err(PortfolioRefusal::MissingAxis {
            axis: EvidenceAxis::ControlledExperimentalValidation,
            ..
        })
    ));
}

#[test]
fn independent_reproduction_requires_a_distinct_group() {
    let experiment = observation(EvidenceAxis::ControlledExperimentalValidation, 1, 9);
    let same_group = observation(EvidenceAxis::IndependentReproduction, 2, 9);
    let portfolio = EvidencePortfolio::try_new(vec![experiment.clone(), same_group]).unwrap();
    assert!(matches!(
        portfolio.admit(
            PortfolioClaimClass::IndependentlyReproducedPrediction,
            QOI,
            REGIME
        ),
        Err(PortfolioRefusal::IndependenceNotEstablished {
            experiment_group,
            ..
        }) if experiment_group == hash(9)
    ));

    let same_source = observation(EvidenceAxis::IndependentReproduction, 1, 10);
    let portfolio = EvidencePortfolio::try_new(vec![experiment.clone(), same_source]).unwrap();
    assert!(matches!(
        portfolio.admit(
            PortfolioClaimClass::IndependentlyReproducedPrediction,
            QOI,
            REGIME
        ),
        Err(PortfolioRefusal::IndependenceNotEstablished { .. })
    ));

    let distinct_group = observation(EvidenceAxis::IndependentReproduction, 3, 10);
    let portfolio =
        EvidencePortfolio::try_new(vec![experiment.clone(), distinct_group.clone()]).unwrap();
    let admission = portfolio
        .admit(
            PortfolioClaimClass::IndependentlyReproducedPrediction,
            QOI,
            REGIME,
        )
        .unwrap();
    assert_eq!(admission.support().len(), 2);
    assert_eq!(
        admission.support()[1].independence_group(),
        distinct_group.independence_group()
    );

    let alternate_experiment = observation(EvidenceAxis::ControlledExperimentalValidation, 4, 10);
    let reproduction = observation(EvidenceAxis::IndependentReproduction, 5, 9);
    let portfolio = EvidencePortfolio::try_new(vec![
        experiment,
        alternate_experiment.clone(),
        reproduction.clone(),
    ])
    .unwrap();
    let admission = portfolio
        .admit(
            PortfolioClaimClass::IndependentlyReproducedPrediction,
            QOI,
            REGIME,
        )
        .unwrap();
    assert_eq!(
        admission.support()[0].source(),
        alternate_experiment.source()
    );
    assert_eq!(admission.support()[1].source(), reproduction.source());
}

#[test]
fn portfolio_and_admission_identities_are_order_stable_and_mutation_sensitive() {
    let numerical = observation(EvidenceAxis::NumericalVerification, 1, 11);
    let cross_code = observation(EvidenceAxis::CrossCodeAgreement, 2, 12);
    let forward = EvidencePortfolio::try_new(vec![numerical.clone(), cross_code.clone()]).unwrap();
    let reverse = EvidencePortfolio::try_new(vec![cross_code, numerical]).unwrap();
    assert_eq!(forward.identity(), reverse.identity());

    let admission = forward
        .admit(PortfolioClaimClass::CrossCodeConsistent, QOI, REGIME)
        .unwrap();
    assert_eq!(
        admission.identity(),
        reverse
            .admit(PortfolioClaimClass::CrossCodeConsistent, QOI, REGIME)
            .unwrap()
            .identity()
    );
    assert!(
        admission
            .render_log()
            .contains("axes=numerical-verification,cross-code-agreement")
    );

    let changed_source = EvidencePortfolio::try_new(vec![
        observation(EvidenceAxis::NumericalVerification, 3, 11),
        observation(EvidenceAxis::CrossCodeAgreement, 2, 12),
    ])
    .unwrap();
    assert_ne!(forward.identity(), changed_source.identity());

    let changed_regime = EvidencePortfolio::try_new(vec![
        PortfolioObservation::try_new(
            EvidenceAxis::NumericalVerification,
            QOI,
            "cooling-vertical-v2",
            hash(1),
            hash(11),
        )
        .unwrap(),
        observation(EvidenceAxis::CrossCodeAgreement, 2, 12),
    ])
    .unwrap();
    assert_ne!(forward.identity(), changed_regime.identity());
}

#[test]
fn support_is_qoi_and_regime_exact() {
    let wrong_qoi = PortfolioObservation::try_new(
        EvidenceAxis::ControlledExperimentalValidation,
        "pressure-drop",
        REGIME,
        hash(1),
        hash(2),
    )
    .unwrap();
    let wrong_regime = PortfolioObservation::try_new(
        EvidenceAxis::ControlledExperimentalValidation,
        QOI,
        "different-regime",
        hash(3),
        hash(4),
    )
    .unwrap();
    let portfolio = EvidencePortfolio::try_new(vec![wrong_qoi, wrong_regime]).unwrap();
    assert!(matches!(
        portfolio.admit(PortfolioClaimClass::ValidatedPrediction, QOI, REGIME),
        Err(PortfolioRefusal::MissingAxis {
            axis: EvidenceAxis::ControlledExperimentalValidation,
            ..
        })
    ));
}

#[test]
fn malformed_observations_refuse_before_portfolio_publication() {
    assert!(matches!(
        PortfolioObservation::try_new(
            EvidenceAxis::NumericalVerification,
            " ",
            REGIME,
            hash(1),
            hash(2),
        ),
        Err(PortfolioRefusal::InvalidIdentifier { field: "qoi", .. })
    ));
    assert!(matches!(
        PortfolioObservation::try_new(
            EvidenceAxis::NumericalVerification,
            QOI,
            REGIME,
            ContentHash([0; 32]),
            hash(2),
        ),
        Err(PortfolioRefusal::ZeroHash { field: "source" })
    ));
}
