//! Actual admitted source-pack identities survive chemistry and active filtering.
use super::*;
use crate::state_point::MaterialPropertySelection;
use crate::thermochemical_gas::packs::{SpeciesPackRef, from_packs};
use fs_blake3::hash_bytes;
use fs_matdb::{
    MODEL_PACK_TARGET_BASIS, ModelNormalizationReceipt, ModelNormalizationTarget,
    NormalizedMaterialCardPack, NormalizedModelPack, NormalizedPack, NormalizedSpeciesPack,
    SPECIES_PACK_TARGET_BASIS, SpeciesAssociation, SpeciesNormalizationReceipt,
    SpeciesNormalizationTarget, ValidityBoundSide,
};
use fs_qty::{MolarMass, Pressure};

fn packs(
    model: Nasa9StandardStateModelV1,
) -> (
    NormalizedSpeciesPack,
    NormalizedModelPack,
    NormalizedMaterialCardPack,
) {
    let name = model.species().as_str();
    let artifact = hash_bytes(format!("manufactured chamber source for {name}").as_bytes());
    let mut card = model.regions()[0].card().clone();
    card.provenance.artifact = Some(artifact);
    card.sources = vec![artifact];
    let association = SpeciesAssociation::new(
        model.species().clone(),
        model.molar_mass().value(),
        "gas",
        "ideal-gas",
        100000.0,
        "synthetic-reference",
        vec![artifact],
        card.provenance.clone(),
    )
    .unwrap();
    let species = NormalizedSpeciesPack::new(
        name,
        "authored-chamber-fixture-v1",
        artifact,
        "CC0-1.0",
        association,
        vec![
            SpeciesNormalizationReceipt::new(
                SpeciesNormalizationTarget::MolarMass,
                hash_bytes(format!("{} kg/mol", model.molar_mass().value()).as_bytes()),
                MolarMass::DIMS,
                1.0,
                0.0,
                "kg/mol",
                SPECIES_PACK_TARGET_BASIS,
            ),
            SpeciesNormalizationReceipt::new(
                SpeciesNormalizationTarget::ReferencePressure,
                hash_bytes(b"100000 Pa"),
                Pressure::DIMS,
                1.0,
                0.0,
                "Pa",
                SPECIES_PACK_TARGET_BASIS,
            ),
        ],
    )
    .unwrap();
    let mut normalizations: Vec<_> = card
        .parameters
        .iter()
        .map(|(name, parameter)| {
            ModelNormalizationReceipt::new(
                ModelNormalizationTarget::Parameter {
                    model: card.content_hash(),
                    parameter: name.clone(),
                },
                hash_bytes(format!("{name}={} SI", parameter.value).as_bytes()),
                parameter.dims,
                1.0,
                0.0,
                "SI",
                MODEL_PACK_TARGET_BASIS,
                None,
                None,
            )
        })
        .collect();
    for (side, t) in [
        (ValidityBoundSide::Lower, 200.0),
        (ValidityBoundSide::Upper, 2000.0),
    ] {
        normalizations.push(ModelNormalizationReceipt::new(
            ModelNormalizationTarget::ValidityBound {
                model: card.content_hash(),
                axis: "T".into(),
                side,
            },
            hash_bytes(format!("{t} K").as_bytes()),
            Temperature::DIMS,
            1.0,
            0.0,
            "K",
            MODEL_PACK_TARGET_BASIS,
            None,
            None,
        ));
    }
    let nasa9 = NormalizedModelPack::new(
        name,
        "authored-chamber-fixture-v1",
        artifact,
        "CC0-1.0",
        vec![card],
        normalizations,
    )
    .unwrap();
    let material = transport(name, 1.7e-5, 0.025).material;
    let mut claims = fs_matdb::ClaimSet::new();
    let provenance = fs_matdb::Provenance {
        source: format!("authored {name} transport curves"),
        license: "CC0-1.0".into(),
        artifact: Some(artifact),
    };
    let observation = claims
        .register_observation(fs_matdb::ObservationDataset {
            specimen: format!("manufactured {name}"),
            method: "analytic test fixture".into(),
            artifact,
            caveats: "no measured physical evidence".into(),
            provenance: provenance.clone(),
        })
        .unwrap();
    for (_, claim) in material.claims().claims_ordered() {
        let mut claim = claim.clone();
        claim.observations = vec![observation];
        claim.provenance = provenance.clone();
        claims.insert_claim(claim).unwrap();
    }
    let transport = NormalizedMaterialCardPack::new(
        material.id().clone(),
        NormalizedPack::new(
            name,
            "authored-chamber-fixture-v1",
            artifact,
            "CC0-1.0",
            claims,
            Vec::new(),
            Vec::new(),
        )
        .unwrap(),
    )
    .unwrap();
    // Consumers see decoded canonical artifacts, not directly assigned hashes.
    (
        NormalizedSpeciesPack::from_bytes(&species.to_bytes()).unwrap(),
        NormalizedModelPack::from_bytes(&nasa9.to_bytes()).unwrap(),
        NormalizedMaterialCardPack::from_bytes(&transport.to_bytes()).unwrap(),
    )
}

#[test]
fn canonical_species_model_and_transport_pack_receipts_survive_equilibrium_freeze() {
    let packs = [
        packs(species("A", 0.01, 4.0, 0.0)),
        packs(species_with_energy_reference("B", 0.01, 4.0, 0.0, -2.4e6)),
        packs(species("I", 0.02, 3.5, 0.0)),
    ];
    let sources: Vec<_> = packs
        .iter()
        .map(|(species, nasa9, transport)| SpeciesPackRef {
            species,
            nasa9: Some(nasa9),
            transport: Some((transport, MaterialPropertySelection::SingleClaimOnly)),
        })
        .collect();
    let bank = from_packs(
        &sources,
        Composition::new(CompositionBasis::MoleFraction, vec![0.25, 0.25, 0.5]).unwrap(),
        None,
    )
    .unwrap();
    let ids = vec![id("A"), id("B"), id("I")];
    let chemistry = IdealGasEquilibriumV1::new(
        ElementalMatrix::new(
            vec![ElementId::new("H").unwrap(), ElementId::new("O").unwrap()],
            ids.clone(),
            vec![vec![1, 1, 0], vec![0, 0, 1]],
        )
        .unwrap(),
        ChargeVector::new(ids, vec![0, 0, 0]).unwrap(),
        bank.thermodynamics().components().to_vec(),
    )
    .unwrap();
    let result = bank
        .equilibrium_chamber_hp(
            &chemistry,
            &point(300.0, 3e6),
            &feed(),
            Temperature::new(300.0),
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(bank.source_packs().len(), 3);
    assert_eq!(
        result.parcel().model().source_packs(),
        &bank.source_packs()[..2]
    );
    assert_eq!(
        result.parcel().state().source_packs(),
        &bank.source_packs()[..2]
    );
    assert_eq!(result.equilibrium().species()[2].log_mole_fraction, None);
    for (receipt, (species, nasa9, transport)) in
        result.parcel().state().source_packs().iter().zip(&packs)
    {
        assert_eq!(receipt.species_pack, species.content_hash());
        assert_eq!(receipt.nasa9_pack, nasa9.content_hash());
        assert_eq!(receipt.transport_pack, transport.content_hash());
    }
}
