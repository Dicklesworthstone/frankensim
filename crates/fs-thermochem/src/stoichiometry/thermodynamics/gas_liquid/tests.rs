//! G0/G3 partition laws with original owners and authored Gibbs fixtures.

use super::*;
use crate::{
    ChargeVector, ElementId, ElementalMatrix, MolarEnergyQuantityV1, ReactionId,
    StoichiometricMatrix,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, NormalizedPack, ObservationDataset, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

const REFERENCE: &str = "source-elemental-reference-authored-partition-test";
const PROPERTY: &str = "standard-molar-formation-gibbs-energy";
const TEMPERATURE: f64 = 298.15;
const P0: f64 = 100_000.0;

fn id(name: &str) -> SpeciesId {
    SpeciesId::new(name).unwrap()
}
fn reaction(scale: i128, charge: i128) -> ConservedReaction {
    let species = vec![id("X(g)"), id("X(aq)")];
    let rid = ReactionId::new("authored-gas-dissolution").unwrap();
    let elements = ElementalMatrix::new(
        vec![ElementId::new("H").unwrap()],
        species.clone(),
        vec![vec![2, 2]],
    )
    .unwrap();
    let matrix = StoichiometricMatrix::new(
        species.clone(),
        vec![rid.clone()],
        vec![vec![-scale], vec![scale]],
    )
    .unwrap();
    let charges = ChargeVector::new(species, vec![charge, charge]).unwrap();
    ConservedReaction::new(&elements, &matrix, &charges, rid).unwrap()
}
fn point(gas: bool) -> QueryPoint {
    let mut at = QueryPoint::new()
        .with_quantity(
            "temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            TEMPERATURE,
        )
        .unwrap()
        .with_quantity(
            "pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            P0,
        )
        .unwrap();
    for flag in [
        REFERENCE,
        "source-standard-state",
        if gas {
            "source-phase-gas"
        } else {
            "source-phase-aqueous"
        },
        if gas {
            "source-pure-species"
        } else {
            "source-aqueous-unit-activity-infinite-dilution"
        },
    ] {
        at = at
            .with_quantity(flag, QuantitySpec::dimensional(Dims::NONE), 1.0)
            .unwrap();
    }
    at
}
fn pack(name: &str, value: f64, at: &QueryPoint, omit: Option<&str>) -> NormalizedPack {
    let mut validity = ValidityDomain::unconstrained();
    for (axis, value) in at.axes() {
        if omit != Some(axis.as_str()) {
            validity = validity.with_quantity(axis, at.axis_quantities()[axis], *value, *value);
        }
    }
    let mut claim = PropertyClaim {
        key: PropertyKey::new(PROPERTY, MolarEnergyQuantityV1::DIMS),
        value: PropertyValue::Scalar {
            value,
            dims: MolarEnergyQuantityV1::DIMS,
        },
        validity,
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: Provenance {
            source: "Authored gas-dissolution algebra fixture, not measured solubility".into(),
            license: "project test fixture".into(),
            artifact: None,
        },
    };
    let artifact = claim.content_hash();
    let mut claims = ClaimSet::new();
    claim.observations.push(
        claims
            .register_observation(ObservationDataset {
                specimen: name.into(),
                method: "manufactured phase Gibbs values".into(),
                artifact,
                caveats: "No empirical Henry coefficient or transfer rate".into(),
                provenance: claim.provenance.clone(),
            })
            .unwrap(),
    );
    claims.insert_claim(claim).unwrap();
    NormalizedPack::new(
        name,
        "fixture-v1",
        artifact,
        "project test fixture",
        claims,
        vec![],
        vec![],
    )
    .unwrap()
}
fn model(log_k: f64) -> IdealGasLiquidPartition {
    let points = [point(true), point(false)];
    let gibbs = if log_k == 0.0 {
        0.0
    } else {
        -UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * TEMPERATURE * log_k
    };
    let packs = [
        pack("gas", 0.0, &points[0], None),
        pack("aqueous", gibbs, &points[1], None),
    ];
    admit(reaction(1, 0), &packs, &points).unwrap()
}
fn admit(
    reaction: ConservedReaction,
    packs: &[NormalizedPack; 2],
    points: &[QueryPoint; 2],
) -> Result<IdealGasLiquidPartition> {
    let ids = [id("X(g)"), id("X(aq)")];
    let sources = (0..2)
        .map(|i| FormationGibbsSource {
            species: &ids[i],
            pack: &packs[i],
            at: &points[i],
            reference_axis: REFERENCE,
        })
        .collect::<Vec<_>>();
    IdealGasLiquidPartition::from_sources(reaction, ids[0].clone(), ids[1].clone(), &sources)
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-12 * expected.abs().max(1e-280),
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn supplied_partial_pressure_uses_the_original_activity_basis_and_sources() {
    let m = model(2.0_f64.ln());
    let amount = m
        .equilibrium_aqueous_amount(Pressure::new(0.5 * P0), Amount::new(3.0))
        .unwrap();
    close(amount.value(), 3.0);
    let pressure = m
        .ideal_partial_pressure(
            Amount::new(2.0),
            Volume::new(2.0 * UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * TEMPERATURE / (0.5 * P0)),
        )
        .unwrap();
    close(pressure.value(), 0.5 * P0);
    assert_eq!(m.gas(), &id("X(g)"));
    assert_eq!(m.aqueous(), &id("X(aq)"));
    assert_eq!(m.standard().reference_axis(), REFERENCE);
    assert_eq!(m.standard().terms().len(), 2);
    assert!(m.standard().terms().iter().all(|t| !t.pack_id.is_empty()));
}

#[test]
fn closed_partition_matches_independent_linear_mass_balance_and_mass_action() {
    let m = model(0.25_f64.ln());
    let total = 7.0;
    let basis = 4.0;
    let volume = 0.1;
    let b = basis * 0.25 * UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * TEMPERATURE / (volume * P0);
    let expected_gas = total / (1.0 + b);
    let expected_aqueous = b * expected_gas;
    let answer = m
        .finite_equilibrium(Amount::new(total), Amount::new(basis), Volume::new(volume))
        .unwrap();
    close(answer.gas_amount.value(), expected_gas);
    close(answer.aqueous_amount.value(), expected_aqueous);
    close(
        answer.partial_pressure.value(),
        expected_gas * UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * TEMPERATURE / volume,
    );
    assert!(answer.component_balance_residual.value().abs() < 3e-15);
    assert!(answer.log_equilibrium_residual.unwrap().value().abs() < 1e-13);
    close(
        m.equilibrium_aqueous_amount(answer.partial_pressure, Amount::new(basis))
            .unwrap()
            .value(),
        answer.aqueous_amount.value(),
    );
    assert_eq!(answer.model().reaction().id(), m.reaction().id());
}

#[test]
fn extensive_scaling_preserves_pressure_and_activity_while_basis_remains_explicit() {
    let m = model(-4.0);
    let reference = m
        .finite_equilibrium(Amount::new(0.3), Amount::new(8.0), Volume::new(0.02))
        .unwrap();
    for factor in [1e-8, 0.5, 8.0, 1e12] {
        let scaled = m
            .finite_equilibrium(
                Amount::new(0.3 * factor),
                Amount::new(8.0 * factor),
                Volume::new(0.02 * factor),
            )
            .unwrap();
        close(
            scaled.gas_amount.value() / factor,
            reference.gas_amount.value(),
        );
        close(
            scaled.aqueous_amount.value() / factor,
            reference.aqueous_amount.value(),
        );
        close(
            scaled.partial_pressure.value(),
            reference.partial_pressure.value(),
        );
    }
    let doubled_basis = m
        .finite_equilibrium(Amount::new(0.3), Amount::new(16.0), Volume::new(0.02))
        .unwrap();
    assert!(doubled_basis.aqueous_amount.value() > reference.aqueous_amount.value());
    assert!(doubled_basis.partial_pressure.value() < reference.partial_pressure.value());
    let repeated = m
        .finite_equilibrium(Amount::new(0.3), Amount::new(8.0), Volume::new(0.02))
        .unwrap();
    assert_eq!(
        reference.gas_amount.value().to_bits(),
        repeated.gas_amount.value().to_bits()
    );
    assert_eq!(
        reference.aqueous_amount.value().to_bits(),
        repeated.aqueous_amount.value().to_bits()
    );
}

#[test]
fn logarithmic_partition_retains_both_positive_tails_when_ratio_underflows() {
    let total: f64 = 1e300;
    let gas_volume = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * TEMPERATURE / P0;
    // exp(-1000) itself underflows, but total*exp(-1000) is representable.
    let expected_minor = (total.ln() - 1000.0_f64).exp();
    for log_k in [-1000.0, 1000.0] {
        let m = model(log_k);
        assert!(m.standard().equilibrium_constant().is_none());
        let answer = m
            .finite_equilibrium(
                Amount::new(total),
                Amount::new(1.0),
                Volume::new(gas_volume),
            )
            .unwrap();
        let (minor, major) = if log_k < 0.0 {
            (answer.aqueous_amount, answer.gas_amount)
        } else {
            (answer.gas_amount, answer.aqueous_amount)
        };
        close(minor.value(), expected_minor);
        assert_eq!(major.value(), total);
        assert!(minor.value() > 0.0);
        // The minor stock remains explicit even though it is below the ULP of total.
        assert_eq!(answer.component_balance_residual, minor);
        assert!(answer.log_equilibrium_residual.unwrap().value().abs() < 1e-12);
    }
    close(
        model(-1000.0)
            .equilibrium_aqueous_amount(Pressure::new(P0), Amount::new(total))
            .unwrap()
            .value(),
        expected_minor,
    );
}

#[test]
fn empty_component_is_distinct_and_nonphysical_or_unrepresentable_inputs_refuse() {
    let m = model(0.0);
    let zero = m
        .finite_equilibrium(Amount::new(0.0), Amount::new(1.0), Volume::new(1.0))
        .unwrap();
    assert_eq!(zero.gas_amount.value(), 0.0);
    assert_eq!(zero.aqueous_amount.value(), 0.0);
    assert_eq!(zero.partial_pressure.value(), 0.0);
    assert!(zero.log_equilibrium_residual.is_none());
    assert_eq!(
        m.equilibrium_aqueous_amount(Pressure::new(0.0), Amount::new(2.0))
            .unwrap()
            .value(),
        0.0
    );
    for bad in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            m.finite_equilibrium(Amount::new(bad), Amount::new(1.0), Volume::new(1.0))
                .is_err()
        );
        assert!(
            m.equilibrium_aqueous_amount(Pressure::new(bad), Amount::new(1.0))
                .is_err()
        );
        assert!(
            m.ideal_partial_pressure(Amount::new(bad), Volume::new(1.0))
                .is_err()
        );
    }
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            m.finite_equilibrium(Amount::new(0.0), Amount::new(bad), Volume::new(1.0))
                .is_err()
        );
        assert!(
            m.finite_equilibrium(Amount::new(0.0), Amount::new(1.0), Volume::new(bad))
                .is_err()
        );
    }
    assert!(matches!(
        model(-1000.0).finite_equilibrium(Amount::new(1.0), Amount::new(1.0), Volume::new(1.0)),
        Err(GasLiquidError::Arithmetic(_))
    ));
    assert!(matches!(
        model(1000.0).equilibrium_aqueous_amount(Pressure::new(P0), Amount::new(1.0)),
        Err(GasLiquidError::Arithmetic(_))
    ));
    assert!(
        m.ideal_partial_pressure(Amount::new(f64::MAX), Volume::new(f64::MIN_POSITIVE))
            .is_err()
    );
}

#[test]
fn scaled_ionic_reversed_and_source_phase_substitutions_refuse() {
    let points = [point(true), point(false)];
    let packs = [
        pack("gas", 0.0, &points[0], None),
        pack("aqueous", 0.0, &points[1], None),
    ];
    for equation in [reaction(2, 0), reaction(-1, 0), reaction(1, 1)] {
        assert!(matches!(
            admit(equation, &packs, &points),
            Err(GasLiquidError::Invalid(_))
        ));
    }
    let phases_swapped = [point(false), point(true)];
    let swapped_packs = [
        pack("gas", 0.0, &phases_swapped[0], None),
        pack("aqueous", 0.0, &phases_swapped[1], None),
    ];
    assert!(admit(reaction(1, 0), &swapped_packs, &phases_swapped).is_err());
    let missing_owned_phase = [
        pack("gas", 0.0, &points[0], Some("source-phase-gas")),
        packs[1].clone(),
    ];
    assert!(admit(reaction(1, 0), &missing_owned_phase, &points).is_err());
    let missing_owned_aqueous = [
        packs[0].clone(),
        pack(
            "aqueous",
            0.0,
            &points[1],
            Some("source-aqueous-unit-activity-infinite-dilution"),
        ),
    ];
    assert!(admit(reaction(1, 0), &missing_owned_aqueous, &points).is_err());
    let wrong_t = point(false)
        .with_quantity(
            "temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            310.0,
        )
        .unwrap();
    let wrong_t_pack = pack("aqueous", 0.0, &wrong_t, None);
    assert!(
        admit(
            reaction(1, 0),
            &[packs[0].clone(), wrong_t_pack],
            &[points[0].clone(), wrong_t]
        )
        .is_err()
    );
}
