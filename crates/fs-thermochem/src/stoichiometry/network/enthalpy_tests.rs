use super::*;
use crate::stoichiometry::minerals::MineralReaction;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, NormalizedPack, ObservationDataset, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

const REFERENCE: &str = "source-elemental-reference-synthetic";
const OTHER: &str = "source-elemental-reference-other-synthetic";

fn reference_point(reference: &str, temperature: f64, species: &SpeciesId) -> QueryPoint {
    let phase = if species.as_str().ends_with("(aq)") {
        "source-phase-aqueous"
    } else if species.as_str().ends_with("(g)") {
        "source-phase-gas"
    } else if species.as_str().ends_with("(l)") {
        "source-phase-liquid"
    } else {
        "source-phase-solid"
    };
    QueryPoint::new()
        .with_quantity(
            "temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            temperature,
        )
        .unwrap()
        .with_quantity(
            "pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            100_000.0,
        )
        .unwrap()
        .with_quantity(reference, QuantitySpec::dimensional(Dims::NONE), 1.0)
        .unwrap()
        .with_quantity(
            "source-standard-state",
            QuantitySpec::dimensional(Dims::NONE),
            1.0,
        )
        .unwrap()
        .with_quantity(phase, QuantitySpec::dimensional(Dims::NONE), 1.0)
        .unwrap()
}

fn pack(species: &SpeciesId, point: &QueryPoint) -> NormalizedPack {
    // Deliberately authored small numbers test the operation tree. They are
    // not mineral thermochemical values; real USGS packs have integration tests.
    let value = match species.as_str() {
        "CaCO3(s)" | "Cu2CO3(OH)2(s)" => -10.0,
        "H+(aq)" => 0.0,
        "Ca[2+](aq)" => -5.0,
        "Cu[2+](aq)" => 2.0,
        "CO2(aq)" => -3.0,
        "H2O(l)" => -4.0,
        _ => -1.0,
    };
    let mut validity = ValidityDomain::unconstrained();
    for (axis, value) in point.axes() {
        validity = validity.with_quantity(axis, point.axis_quantities()[axis], *value, *value);
    }
    let dims = Dims([2, 1, -2, 0, 0, -1]);
    let mut claim = PropertyClaim {
        key: PropertyKey::new("standard-molar-formation-enthalpy", dims),
        value: PropertyValue::Scalar { value, dims },
        validity,
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: Provenance {
            source: "authored network heat fixture, not measured mineral data".into(),
            license: "project fixture".into(),
            artifact: None,
        },
    };
    let hash = claim.content_hash();
    let mut claims = ClaimSet::new();
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: species.as_str().into(),
            method: "authored exact-arithmetic fixture".into(),
            artifact: hash,
            caveats: "synthetic values, not physical observations".into(),
            provenance: claim.provenance.clone(),
        })
        .unwrap();
    claim.observations.push(observation);
    claims.insert_claim(claim).unwrap();
    NormalizedPack::new(
        species.as_str(),
        "fixture",
        hash,
        "project fixture",
        claims,
        vec![],
        vec![],
    )
    .unwrap()
}

fn inventory(network: &ReactionNetwork, acid: f64) -> Vec<(SpeciesId, Amount)> {
    network
        .species()
        .iter()
        .map(|id| {
            let value = match id.as_str() {
                "CaCO3(s)" | "Cu2CO3(OH)2(s)" => 10.0,
                "H+(aq)" => acid,
                _ => 0.0,
            };
            (id.clone(), Amount::new(value))
        })
        .collect()
}

#[test]
fn shared_reagent_stage_heat_keeps_original_sources_and_admitted_extents() {
    let network = ReactionNetwork::new(vec![
        MineralReaction::MalachiteAcidLeach.build().unwrap(),
        MineralReaction::CalciteAcidNeutralization.build().unwrap(),
    ])
    .unwrap();
    let points: Vec<_> = network
        .species()
        .iter()
        .map(|id| reference_point(REFERENCE, 298.15, id))
        .collect();
    let packs: Vec<_> = network
        .species()
        .iter()
        .zip(&points)
        .map(|(id, point)| pack(id, point))
        .collect();
    let sources: Vec<_> = network
        .species()
        .iter()
        .zip(&points)
        .zip(&packs)
        .map(|((species, at), pack)| FormationEnthalpySource {
            species,
            pack,
            at,
            reference_axis: REFERENCE,
        })
        .collect();
    let requests: Vec<_> = network
        .reactions()
        .iter()
        .map(|reaction| {
            (
                reaction.id().clone(),
                Amount::new(if reaction.id().as_str().starts_with("calcite") {
                    3.0
                } else {
                    7.0
                }),
            )
        })
        .collect();
    let initial = inventory(&network, 34.0);
    let result = network
        .advance_with_standard_enthalpy(&initial, &requests, &sources)
        .unwrap();
    // Synthetic Hess sums: calcite -2 J/mol, malachite -1 J/mol.
    assert_eq!(result.enthalpy.value(), -13.0);
    assert_eq!(result.reactions[0].enthalpy.value(), -6.0);
    assert_eq!(result.reactions[1].enthalpy.value(), -7.0);
    assert_eq!(
        result
            .inventory
            .final_inventory
            .iter()
            .find(|(id, _)| id.as_str() == "H+(aq)")
            .unwrap()
            .1
            .value(),
        0.0
    );
    for contribution in &result.reactions {
        for term in contribution.standard.terms() {
            let source = sources
                .iter()
                .find(|source| source.species == &term.term.species)
                .unwrap();
            assert_eq!(term.pack_id, source.pack.pack_id());
            source
                .pack
                .claims()
                .verify_receipt(&term.answer.receipt)
                .unwrap();
        }
    }
    assert!(
        network
            .advance_with_standard_enthalpy(&inventory(&network, 33.0), &requests, &sources)
            .is_err()
    );
    assert_eq!(
        initial,
        inventory(&network, 34.0),
        "failure must not alter inputs"
    );
    assert!(
        network
            .advance_with_standard_enthalpy(&initial, &requests, &sources[..sources.len() - 1])
            .is_err()
    );
    let mut duplicate: Vec<_> = sources
        .iter()
        .map(|source| FormationEnthalpySource {
            species: source.species,
            pack: source.pack,
            at: source.at,
            reference_axis: source.reference_axis,
        })
        .collect();
    duplicate[0] = FormationEnthalpySource {
        species: sources[1].species,
        pack: sources[1].pack,
        at: sources[1].at,
        reference_axis: REFERENCE,
    };
    assert!(
        network
            .advance_with_standard_enthalpy(&initial, &requests, &duplicate)
            .is_err()
    );
}

#[test]
fn disconnected_zero_extent_reactions_still_require_one_reference_state() {
    let network = ReactionNetwork::new(vec![
        MineralReaction::LimeHydration.build().unwrap(),
        MineralReaction::CupricOxideCarbonMonoxideReduction
            .build()
            .unwrap(),
    ])
    .unwrap();
    let initial = inventory(&network, 0.0);
    let requests: Vec<_> = network
        .reactions()
        .iter()
        .map(|reaction| (reaction.id().clone(), Amount::new(0.0)))
        .collect();
    for (second_reference, second_temperature) in
        [(OTHER, 298.15), (REFERENCE, 310.0), (REFERENCE, 298.15)]
    {
        let references: Vec<_> = network
            .species()
            .iter()
            .map(|id| {
                if id.as_str().starts_with("Cu") || id.as_str().ends_with("(g)") {
                    (second_reference, second_temperature)
                } else {
                    (REFERENCE, 298.15)
                }
            })
            .collect();
        let points: Vec<_> = network
            .species()
            .iter()
            .zip(&references)
            .map(|(id, (reference, temperature))| reference_point(reference, *temperature, id))
            .collect();
        let packs: Vec<_> = network
            .species()
            .iter()
            .zip(&points)
            .map(|(id, point)| pack(id, point))
            .collect();
        let sources: Vec<_> = network
            .species()
            .iter()
            .zip(&points)
            .zip(&packs)
            .zip(&references)
            .map(
                |(((species, at), pack), (reference, _))| FormationEnthalpySource {
                    species,
                    pack,
                    at,
                    reference_axis: reference,
                },
            )
            .collect();
        let result = network.advance_with_standard_enthalpy(&initial, &requests, &sources);
        if second_reference == REFERENCE && second_temperature == 298.15 {
            assert_eq!(result.unwrap().enthalpy.value(), 0.0);
        } else {
            assert!(matches!(
                result,
                Err(ReactionNetworkError::Balance(
                    ReactionBalanceError::Invalid("network formation reference states differ")
                ))
            ));
        }
    }
}
