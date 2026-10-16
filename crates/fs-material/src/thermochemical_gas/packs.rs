//! Bind existing verified artifacts without interpreting a species name as data.

use fs_blake3::ContentHash;
use fs_matdb::{NormalizedMaterialCardPack, NormalizedModelPack, NormalizedSpeciesPack};
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
    let regions = nasa9
        .models()
        .iter()
        .cloned()
        .map(Nasa9RegionV1::from_card)
        .collect::<Result<Vec<_>, _>>()
        .map_err(|error| failure("NASA-9 regions", error))?;
    Nasa9StandardStateModelV1::new(
        association.species().clone(),
        MolarMass::new(association.molar_mass()),
        convention,
        regions,
    )
    .map_err(|error| failure("NASA-9 species", error))
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
