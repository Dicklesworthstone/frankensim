//! G0/G3 known answers, exact-reference cancellation, and source-domain refusals.
use super::*;
use crate::{ChargeVector, ElementId, ElementalMatrix, StoichiometricMatrix};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, ObservationDataset, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, UncertaintyModel,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

const REFERENCE: &str = "source-elemental-reference-authored-gibbs-test";
const PROPERTY: &str = "standard-molar-formation-gibbs-energy";

fn id(value: &str) -> SpeciesId {
    SpeciesId::new(value).unwrap()
}

fn exchange(scale: i128, charge: i128) -> ConservedReaction {
    // Two authored states with identical atoms and charge; these are not
    // measured mineral data or a claim about a real isomerization.
    let species = vec![id("state-a"), id("state-b")];
    let reaction = ReactionId::new("authored-isomerization").unwrap();
    let elemental = ElementalMatrix::new(
        vec![ElementId::new("H").unwrap()],
        species.clone(),
        vec![vec![1, 1]],
    )
    .unwrap();
    let stoichiometric = StoichiometricMatrix::new(
        species.clone(),
        vec![reaction.clone()],
        vec![vec![-scale], vec![scale]],
    )
    .unwrap();
    let charges = ChargeVector::new(species, vec![charge, charge]).unwrap();
    ConservedReaction::new(&elemental, &stoichiometric, &charges, reaction).unwrap()
}

fn point(phase: &str, ionic: bool) -> QueryPoint {
    let mut point = QueryPoint::new()
        .with_quantity(
            "temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            298.15,
        )
        .unwrap()
        .with_quantity(
            "pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            100_000.0,
        )
        .unwrap();
    let mut flags = vec![REFERENCE, "source-standard-state", phase];
    if phase == "source-phase-aqueous" {
        flags.push("source-aqueous-unit-activity-infinite-dilution");
        if ionic {
            flags.extend([
                "source-ionic-reference-hydrogen-zero",
                "source-individual-ion",
            ]);
        }
    } else {
        flags.push("source-pure-species");
    }
    for flag in flags {
        point = point
            .with_quantity(flag, QuantitySpec::dimensional(Dims::NONE), 1.0)
            .unwrap();
    }
    point
}

fn pack(value: f64, point: &QueryPoint, omit: Option<&str>, property: &str) -> NormalizedPack {
    let mut validity = ValidityDomain::unconstrained();
    for (axis, value) in point.axes() {
        if omit != Some(axis.as_str()) {
            validity = validity.with_quantity(axis, point.axis_quantities()[axis], *value, *value);
        }
    }
    let dims = MolarEnergyQuantityV1::DIMS;
    let mut claim = PropertyClaim {
        key: PropertyKey::new(property, dims),
        value: PropertyValue::Scalar { value, dims },
        validity,
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: Provenance {
            source: "authored Gibbs arithmetic fixture".into(),
            license: "project fixture".into(),
            artifact: None,
        },
    };
    let hash = claim.content_hash();
    let mut claims = ClaimSet::new();
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: "authored exactly conserved states".into(),
            method: "synthetic numerical fixture".into(),
            artifact: hash,
            caveats: "not experimental thermochemistry".into(),
            provenance: claim.provenance.clone(),
        })
        .unwrap();
    claim.observations.push(observation);
    claims.insert_claim(claim).unwrap();
    NormalizedPack::new(
        "fixture",
        "test",
        hash,
        "project fixture",
        claims,
        vec![],
        vec![],
    )
    .unwrap()
}

fn sources<'a>(
    reaction: &'a ConservedReaction,
    packs: &'a [NormalizedPack],
    points: &'a [QueryPoint],
) -> Vec<FormationGibbsSource<'a>> {
    reaction
        .terms()
        .iter()
        .zip(packs)
        .zip(points)
        .map(|((term, pack), at)| FormationGibbsSource {
            species: &term.species,
            pack,
            at,
            reference_axis: REFERENCE,
        })
        .collect()
}

fn evaluate(scale: i128, values: [f64; 2]) -> StandardReactionGibbs {
    let reaction = exchange(scale, 0);
    let points = [
        point("source-phase-solid", false),
        point("source-phase-solid", false),
    ];
    let packs = [
        pack(values[0], &points[0], None, PROPERTY),
        pack(values[1], &points[1], None, PROPERTY),
    ];
    reaction
        .standard_gibbs_change(&sources(&reaction, &packs, &points))
        .unwrap()
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-14 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

#[test]
fn decimal_known_answers_and_activity_direction_change() {
    // Independent 60-digit Decimal oracle, R=8.31446261815324, T=298.15.
    let standard = evaluate(1, [0.0, -10_000.0]);
    close(standard.molar().value(), -10_000.0);
    close(
        standard.log_equilibrium_constant().value(),
        4.033_954_554_510_348_5,
    );
    close(
        standard.equilibrium_constant().unwrap().value(),
        56.483_838_586_092_66,
    );
    assert_eq!(
        standard.status(),
        EquilibriumConstantStatusV1::FinitePositive
    );
    assert_eq!(standard.temperature().value(), 298.15);
    assert_eq!(standard.reference_pressure().value(), 100_000.0);
    assert_eq!(standard.reference_axis(), REFERENCE);
    for (quotient, expected) in [
        (1.0, -10_000.0),
        (10.0, -4291.990_497_464_741),
        (100.0, 1416.019_005_070_518),
    ] {
        let result = standard
            .at_activities(&[
                (id("state-b"), Dimensionless::new(quotient)),
                (id("state-a"), Dimensionless::new(1.0)),
            ])
            .unwrap();
        close(result.molar.value(), expected);
        close(result.affinity().value(), -expected);
        assert_eq!(result.activities[0].0, id("state-a"));
    }
}

#[test]
fn equation_reversal_scaling_and_reference_shifts() {
    let forward = evaluate(1, [0.0, -10_000.0]);
    let reversed = evaluate(-1, [0.0, -10_000.0]);
    let doubled = evaluate(2, [0.0, -10_000.0]);
    let shifted = evaluate(1, [1_000_000.0, 990_000.0]);
    assert_eq!(reversed.molar().value(), -forward.molar().value());
    assert_eq!(doubled.molar().value(), 2.0 * forward.molar().value());
    assert_eq!(shifted.molar().value(), forward.molar().value());
    assert_eq!(
        shifted.log_equilibrium_constant(),
        forward.log_equilibrium_constant()
    );
    let activities = [
        (id("state-a"), Dimensionless::new(0.1)),
        (id("state-b"), Dimensionless::new(0.5)),
    ];
    let a = forward.at_activities(&activities).unwrap();
    let b = reversed.at_activities(&activities).unwrap();
    let c = doubled.at_activities(&activities).unwrap();
    assert_eq!(b.log_quotient.value(), -a.log_quotient.value());
    assert_eq!(b.molar.value(), -a.molar.value());
    close(c.molar.value(), 2.0 * a.molar.value());
}

#[test]
fn source_and_activity_permutations_preserve_bits_and_receipts() {
    let reaction = exchange(1, 0);
    let points = [
        point("source-phase-solid", false),
        point("source-phase-liquid", false),
    ];
    let packs = [
        pack(100.0, &points[0], None, PROPERTY),
        pack(-9900.0, &points[1], None, PROPERTY),
    ];
    let mut input = sources(&reaction, &packs, &points);
    let first = reaction.standard_gibbs_change(&input).unwrap();
    input.reverse();
    let second = reaction.standard_gibbs_change(&input).unwrap();
    assert_eq!(
        first.molar().value().to_bits(),
        second.molar().value().to_bits()
    );
    for (i, term) in first.terms().iter().enumerate() {
        packs[i]
            .claims()
            .verify_receipt(&term.answer.receipt)
            .unwrap();
        assert_eq!(term.answer.receipt, second.terms()[i].answer.receipt);
        assert_eq!(
            term.answer.evidence.value.uncertainty,
            UncertaintyModel::Unstated
        );
    }
    let mut activities = vec![
        (id("state-a"), Dimensionless::new(0.2)),
        (id("state-b"), Dimensionless::new(1.7)),
    ];
    let before = first.at_activities(&activities).unwrap();
    activities.reverse();
    let after = first.at_activities(&activities).unwrap();
    assert_eq!(
        before.molar.value().to_bits(),
        after.molar.value().to_bits()
    );
    assert_eq!(before.activities, after.activities);
}

#[test]
fn finite_logs_survive_direct_constant_and_quotient_overflow() {
    for (energy, status) in [
        (-10_000_000.0, EquilibriumConstantStatusV1::LogOnlyOverflow),
        (10_000_000.0, EquilibriumConstantStatusV1::LogOnlyUnderflow),
    ] {
        let result = evaluate(1, [0.0, energy]);
        assert!(result.log_equilibrium_constant().value().is_finite());
        assert_eq!(result.status(), status);
        assert!(result.equilibrium_constant().is_none());
        let at = result
            .at_activities(&[
                (id("state-a"), Dimensionless::new(1e-300)),
                (id("state-b"), Dimensionless::new(1e300)),
            ])
            .unwrap();
        close(at.log_quotient.value(), 1381.551_055_796_427_4);
        assert!(at.molar.value().is_finite());
    }
}

#[test]
fn invalid_or_incomplete_activities_refuse_without_defaults() {
    let standard = evaluate(1, [0.0, -10_000.0]);
    for value in [0.0, -0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(
            standard
                .at_activities(&[
                    (id("state-a"), Dimensionless::new(1.0)),
                    (id("state-b"), Dimensionless::new(value))
                ])
                .is_err()
        );
    }
    let a = (id("state-a"), Dimensionless::new(1.0));
    let b = (id("state-b"), Dimensionless::new(1.0));
    assert!(standard.at_activities(&[]).is_err());
    assert!(standard.at_activities(&[a.clone()]).is_err());
    assert!(standard.at_activities(&[a.clone(), a.clone()]).is_err());
    assert!(
        standard
            .at_activities(&[a.clone(), (id("other"), Dimensionless::new(1.0))])
            .is_err()
    );
    assert!(standard.at_activities(&[a.clone(), b, a]).is_err());
}

#[test]
fn aqueous_and_ionic_conventions_belong_to_selected_claims() {
    for (charge, required) in [
        (0, "source-aqueous-unit-activity-infinite-dilution"),
        (1, "source-ionic-reference-hydrogen-zero"),
        (1, "source-individual-ion"),
    ] {
        let reaction = exchange(1, charge);
        let points = [
            point("source-phase-aqueous", charge != 0),
            point("source-phase-aqueous", charge != 0),
        ];
        let mut packs = [
            pack(0.0, &points[0], None, PROPERTY),
            pack(-10_000.0, &points[1], None, PROPERTY),
        ];
        assert!(
            reaction
                .standard_gibbs_change(&sources(&reaction, &packs, &points))
                .is_ok()
        );
        packs[1] = pack(-10_000.0, &points[1], Some(required), PROPERTY);
        // The query still supplies the flag; that cannot repair absent source authority.
        assert!(
            reaction
                .standard_gibbs_change(&sources(&reaction, &packs, &points))
                .is_err()
        );
    }
    let charged = exchange(1, 1);
    let points = [
        point("source-phase-gas", false),
        point("source-phase-solid", false),
    ];
    let packs = [
        pack(0.0, &points[0], None, PROPERTY),
        pack(-10_000.0, &points[1], None, PROPERTY),
    ];
    assert!(
        charged
            .standard_gibbs_change(&sources(&charged, &packs, &points))
            .is_err()
    );
}

#[test]
fn source_coverage_property_reference_and_domain_refusals() {
    let reaction = exchange(1, 0);
    let mut points = [
        point("source-phase-solid", false),
        point("source-phase-solid", false),
    ];
    let packs = [
        pack(0.0, &points[0], None, PROPERTY),
        pack(-10_000.0, &points[1], None, PROPERTY),
    ];
    let mut input = sources(&reaction, &packs, &points);
    assert!(reaction.standard_gibbs_change(&input[..1]).is_err());
    input[1].species = input[0].species;
    assert!(reaction.standard_gibbs_change(&input).is_err());
    input[1].species = &reaction.terms()[1].species;
    input[1].reference_axis = "source-elemental-reference-other";
    assert!(reaction.standard_gibbs_change(&input).is_err());
    for omitted in [
        "temperature",
        "pressure",
        "source-standard-state",
        REFERENCE,
        "source-pure-species",
        "source-phase-solid",
    ] {
        let bad = [
            pack(0.0, &points[0], None, PROPERTY),
            pack(-10_000.0, &points[1], Some(omitted), PROPERTY),
        ];
        assert!(
            reaction
                .standard_gibbs_change(&sources(&reaction, &bad, &points))
                .is_err(),
            "{omitted}"
        );
    }
    let enthalpy = [
        pack(0.0, &points[0], None, PROPERTY),
        pack(
            -10_000.0,
            &points[1],
            None,
            "standard-molar-formation-enthalpy",
        ),
    ];
    assert!(
        reaction
            .standard_gibbs_change(&sources(&reaction, &enthalpy, &points))
            .is_err()
    );
    points[1] = points[1]
        .clone()
        .with_quantity(
            "temperature",
            points[1].axis_quantities()["temperature"],
            300.0,
        )
        .unwrap();
    assert!(
        reaction
            .standard_gibbs_change(&sources(&reaction, &packs, &points))
            .is_err()
    );
}
