//! Bind existing verified artifacts without interpreting a species name as data.

use fs_blake3::ContentHash;
use fs_matdb::{
    ConstitutiveModelCard, NormalizedMaterialCardPack, NormalizedModelPack, NormalizedSpeciesPack,
};
use fs_qty::{MolarMass, Pressure};
use fs_thermochem::{
    Composition, ElementalReferenceIdV1, Nasa9RegionV1, Nasa9StandardStateModelV1,
    ReferenceEquationOfStateV1, SpeciesId, StandardStateConventionV1, StandardStatePhaseV1,
    mixture::{FrozenIdealGasMixtureModelV1, MAX_FROZEN_IDEAL_GAS_COMPONENTS_V1},
};

use super::{MaterialError, SpeciesTransportMaterial, ThermochemicalGasModel, failure, invalid};
use crate::state_point::MaterialPropertySelection;

/// One caller-declared association. Optional artifacts allow a missing required
/// source to produce a named refusal rather than a fallback material.
#[derive(Clone, Debug)]
pub struct SpeciesPackRef<'a> {
    /// Molar mass and explicit standard-state convention.
    pub species: &'a NormalizedSpeciesPack,
    /// NASA-9 regions whose canonical pack id must equal that species id.
    pub nasa9: Option<&'a NormalizedModelPack>,
    /// Pure-species transport and its explicit claim-selection policy.
    pub transport: Option<(&'a NormalizedMaterialCardPack, MaterialPropertySelection)>,
}

/// Immutable artifact join retained by the model and every resolved state.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ThermochemicalPackReceipt {
    /// Canonical identity, in the mixture's canonical component order.
    pub species: SpeciesId,
    /// Verified species association pack identity.
    pub species_pack: ContentHash,
    /// Verified NASA-9 model pack identity.
    pub nasa9_pack: ContentHash,
    /// Verified pure-transport material-card pack identity.
    pub transport_pack: ContentHash,
}

/// Construct one standard-state species from the existing pack formats.
/// No molar mass, reference pressure, coefficient or species alias is inferred.
pub fn standard_state_from_packs(
    species: &NormalizedSpeciesPack,
    nasa9: &NormalizedModelPack,
) -> Result<Nasa9StandardStateModelV1, MaterialError> {
    let association = species.association();
    if nasa9.pack_id() != association.species().as_str() {
        return Err(invalid(
            "NASA-9 pack id differs from the explicit species association",
        ));
    }
    if association.standard_state_phase() != "gas" || association.reference_eos() != "ideal-gas" {
        return Err(invalid(
            "thermochemical gas needs an explicit gas/ideal-gas standard state",
        ));
    }
    let reference = ElementalReferenceIdV1::new(association.elemental_reference())
        .map_err(|error| failure("elemental reference", error))?;
    let convention = StandardStateConventionV1::new(
        StandardStatePhaseV1::Gas,
        ReferenceEquationOfStateV1::IdealGas,
        Pressure::new(association.reference_pressure()),
        reference,
    )
    .map_err(|error| failure("standard state", error))?;
    Nasa9StandardStateModelV1::new(
        association.species().clone(),
        MolarMass::new(association.molar_mass()),
        convention,
        ordered_regions(nasa9.models())?,
    )
    .map_err(|error| failure("NASA-9 species", error))
}

fn ordered_regions(cards: &[ConstitutiveModelCard]) -> Result<Vec<Nasa9RegionV1>, MaterialError> {
    let mut regions = cards
        .iter()
        .cloned()
        .map(Nasa9RegionV1::from_card)
        .collect::<Result<Vec<_>, _>>()
        .map_err(|error| failure("NASA-9 regions", error))?;
    // Portable packs are canonicalized by card hash, not by temperature.
    // Recover physical interval order without changing a source coefficient;
    // the thermochemistry owner still refuses overlapping intervals.
    regions.sort_by(|a, b| {
        a.temperature_min()
            .value()
            .total_cmp(&b.temperature_min().value())
    });
    Ok(regions)
}

/// Fractions refer to `species` in the supplied order. The existing mixture
/// owner sorts each fraction together with its species, checks conventions,
/// converts a declared mass basis, and enforces the 128-component bound.
pub fn from_packs(
    species: &[SpeciesPackRef<'_>],
    composition: Composition,
    water: Option<&SpeciesId>,
) -> Result<ThermochemicalGasModel, MaterialError> {
    if species.len() > MAX_FROZEN_IDEAL_GAS_COMPONENTS_V1 {
        return Err(invalid(format!(
            "thermochemical component count exceeds MAX_FROZEN_IDEAL_GAS_COMPONENTS_V1 ({MAX_FROZEN_IDEAL_GAS_COMPONENTS_V1})"
        )));
    }
    let mut models = Vec::with_capacity(species.len());
    let mut transport = Vec::with_capacity(species.len());
    let mut receipts = Vec::with_capacity(species.len());
    for source in species {
        let id = source.species.association().species();
        let nasa9 = source
            .nasa9
            .ok_or_else(|| invalid(format!("species {id} has no NASA-9 pack")))?;
        let (pack, selection) = source
            .transport
            .as_ref()
            .ok_or_else(|| invalid(format!("species {id} has no transport card")))?;
        models.push(standard_state_from_packs(source.species, nasa9)?);
        transport.push(SpeciesTransportMaterial {
            species: id.clone(),
            material: pack.card().clone(),
            selection: selection.clone(),
        });
        receipts.push(ThermochemicalPackReceipt {
            species: id.clone(),
            species_pack: source.species.content_hash(),
            nasa9_pack: nasa9.content_hash(),
            transport_pack: pack.content_hash(),
        });
    }
    let thermodynamics = FrozenIdealGasMixtureModelV1::new(models, composition)
        .map_err(|error| failure("mixture construction", error))?;
    let mut model = ThermochemicalGasModel::new(thermodynamics, transport, water)?;
    receipts.sort_by(|a, b| a.species.cmp(&b.species));
    model.source_packs = receipts;
    Ok(model)
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_evidence::ValidityDomain;

    #[test]
    fn g0_hash_order_is_not_temperature_order_and_overlap_still_refuses() {
        let fixture = super::super::tests::species("synthetic-region-order", 0.028, 3.5, 0.0);
        let mut low = fixture.regions()[0].card().clone();
        let mut high = low.clone();
        low.validity = ValidityDomain::unconstrained().with("T", 200.0, 1000.0);
        high.validity = ValidityDomain::unconstrained().with("T", 1000.0, 6000.0);
        high.parameters.get_mut("a2").unwrap().value = 4.0;
        for cards in [
            vec![high.clone(), low.clone()],
            vec![low.clone(), high.clone()],
        ] {
            let regions = ordered_regions(&cards).unwrap();
            assert_eq!(
                regions
                    .iter()
                    .map(|r| (r.temperature_min().value(), r.temperature_max().value()))
                    .collect::<Vec<_>>(),
                vec![(200.0, 1000.0), (1000.0, 6000.0)]
            );
            for (region, source) in regions.iter().zip([&low, &high]) {
                assert_eq!(region.card_identity(), source.content_hash().as_bytes());
                for (i, name) in ["a0", "a1", "a2", "a3", "a4", "a5", "a6", "a7", "a8"]
                    .iter()
                    .enumerate()
                {
                    assert_eq!(
                        region.coefficients()[i].to_bits(),
                        source.parameters[*name].value.to_bits()
                    );
                }
            }
            assert!(
                Nasa9StandardStateModelV1::new(
                    fixture.species().clone(),
                    fixture.molar_mass(),
                    fixture.convention().clone(),
                    regions
                )
                .is_ok()
            );
        }
        low.validity = ValidityDomain::unconstrained().with("T", 200.0, 1200.0);
        let overlap = ordered_regions(&[high, low]).unwrap();
        assert!(
            Nasa9StandardStateModelV1::new(
                fixture.species().clone(),
                fixture.molar_mass(),
                fixture.convention().clone(),
                overlap
            )
            .is_err()
        );
    }
}
