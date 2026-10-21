//! Cold admission of canonical species, NASA-9 and transport material packs.
//! The explicit triple assigns transport to a species; no name selects data.
//! Numerical evaluation remains in `read_gas_record` and its original owners.

use super::{GasInputError, GasRecord, ThermochemicalPhase, ThermochemicalSources, read_gas_record};
use crate::state_point::MaterialPropertySelection;
use crate::thermochemical_gas::packs::SpeciesPackRef;
use fs_matdb::{NormalizedMaterialCardPack, NormalizedModelPack, NormalizedSpeciesPack};
use fs_thermochem::{SpeciesId, mixture::MAX_FROZEN_IDEAL_GAS_COMPONENTS_V1};

pub const MAX_SOURCE_PACK_BYTES: usize = 16 * 1024 * 1024;
pub const MAX_SOURCE_BYTES: usize = 64 * 1024 * 1024;
pub const MAX_RECORD_BYTES: usize = 8192;

/// Three unchanged canonical pack encodings, explicitly paired by the caller.
#[derive(Clone, Copy, Debug)]
pub struct SourcePackBytes<'a> {
    pub species: &'a [u8],
    pub nasa9: &'a [u8],
    pub transport: &'a [u8],
}

/// Resolve a complete thermochemical record without filesystem or store access.
/// The returned record owns the evaluated mixture, source models, transport
/// claims and pack receipts after these borrowed bytes have gone away.
///
/// All transport properties use the existing `SingleClaimOnly` policy. Import
/// selected claims as a canonical card when a source contains alternatives, or
/// use the original borrowed registry API with an explicit selection policy.
/// The supplied registry may contain unused species, but an atmosphere with
/// supplied sources must use the thermochemical grammar; sources are never
/// silently ignored in favor of a legacy gas. Phase and water identity are
/// explicit, independent decisions, not inferred from chemical spelling.
pub fn read_source_packs(
    text: &str,
    sources: &[SourcePackBytes<'_>],
    water_species: Option<&str>,
    phase: ThermochemicalPhase,
) -> Result<GasRecord, GasInputError> {
    if text.len() > MAX_RECORD_BYTES
        || text.split_ascii_whitespace().next() != Some("thermochemical") {
        return Err(GasInputError::Input("source packs require a thermochemical record of at most 8 KiB"));
    }
    if sources.is_empty() || sources.len() > MAX_FROZEN_IDEAL_GAS_COMPONENTS_V1 {
        return Err(GasInputError::Input("gas source registry requires 1..128 explicit pack triples"));
    }
    let mut size = 0usize;
    for source in sources {
        for bytes in [source.species, source.nasa9, source.transport] {
            size = size.checked_add(bytes.len()).ok_or(GasInputError::Input("gas source byte overflow"))?;
            if bytes.is_empty() || bytes.len() > MAX_SOURCE_PACK_BYTES || size > MAX_SOURCE_BYTES {
                return Err(GasInputError::Input("nonempty gas packs must fit 16 MiB each and 64 MiB in total"));
            }
        }
    }
    let water = water_species.map(|name| SpeciesId::new(name)
        .map_err(|_| GasInputError::Input("invalid explicit water species identity"))).transpose()?;
    if phase == ThermochemicalPhase::LiquidWaterLimit && water.is_none() {
        return Err(GasInputError::Input("liquid-water admission requires an explicit water identity"));
    }
    let mut loaded = Vec::with_capacity(sources.len());
    for source in sources {
        let species = NormalizedSpeciesPack::from_bytes(source.species)
            .map_err(|_| GasInputError::Input("invalid canonical species pack"))?;
        let nasa9 = NormalizedModelPack::from_bytes(source.nasa9)
            .map_err(|_| GasInputError::Input("invalid canonical NASA-9 model pack"))?;
        let transport = NormalizedMaterialCardPack::from_bytes(source.transport)
            .map_err(|_| GasInputError::Input("invalid canonical transport material-card pack"))?;
        loaded.push((species, nasa9, transport));
    }
    let bindings = loaded.iter().map(|(species, nasa9, transport)| SpeciesPackRef {
        species, nasa9: Some(nasa9),
        transport: Some((transport, MaterialPropertySelection::SingleClaimOnly)),
    }).collect::<Vec<_>>();
    let registry = ThermochemicalSources { species: &bindings, water_species: water.as_ref(), phase };
    let mut fields = text.split_ascii_whitespace();
    let record = read_gas_record(&mut fields, None, Some(&registry))?;
    if fields.next().is_some() {
        return Err(GasInputError::Input("unexpected fields after the gas declaration"));
    }
    Ok(record)
}

/// Length framing for WASM: each consecutive species/model/transport triple
/// refers to the unchanged canonical bytes in `bytes`. Every byte is consumed.
/// This does not define an alternate material serialization or gas law.
pub fn read_source_packet(
    text: &str,
    bytes: &[u8],
    lengths: &[u32],
    water_species: Option<&str>,
    phase: ThermochemicalPhase,
) -> Result<GasRecord, GasInputError> {
    if bytes.len() > MAX_SOURCE_BYTES || lengths.is_empty()
        || lengths.len() % 3 != 0 || lengths.len() > 3 * MAX_FROZEN_IDEAL_GAS_COMPONENTS_V1 {
        return Err(GasInputError::Input("gas packet requires complete source triples within the byte/count budget"));
    }
    let mut offset = 0usize;
    let mut take = |length: u32| -> Result<&[u8], GasInputError> {
        if length == 0 || length as usize > MAX_SOURCE_PACK_BYTES {
            return Err(GasInputError::Input("gas packet pack must contain 1 byte..16 MiB"));
        }
        let end = offset.checked_add(length as usize)
            .ok_or(GasInputError::Input("gas packet framing overflow"))?;
        let part = bytes.get(offset..end).ok_or(GasInputError::Input("truncated gas source packet"))?;
        offset = end;
        Ok(part)
    };
    let sources = lengths.chunks_exact(3).map(|triple| Ok(SourcePackBytes {
        species: take(triple[0])?, nasa9: take(triple[1])?, transport: take(triple[2])?,
    })).collect::<Result<Vec<_>, GasInputError>>()?;
    if offset != bytes.len() {
        return Err(GasInputError::Input("trailing bytes after gas source packs"));
    }
    read_source_packs(text, &sources, water_species, phase)
}
