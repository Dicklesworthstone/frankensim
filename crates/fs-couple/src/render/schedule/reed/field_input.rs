//! Authored primitive gas fields over the existing pack and field owners.
//!
//! This module parses data and associates explicit source IDs. It introduces no
//! thermodynamic, transport, interpolation, duct or pressure-observation law.

use super::{GasFieldResolution, GasRecord, Reader, ReedPerformanceError, ThermochemicalSources, bad};
use fs_material::gas_field::{AxialGasField, MAX_FIELD_STATIONS, PrimitiveGasState};
use fs_material::thermochemical_gas::packs::{SpeciesPackRef, from_packs};
use fs_matdb::QueryPoint;
use fs_qty::semantic::{Composition, CompositionBasis};

// Cold input budget, independently rechecked by the existing mixture owner.
const MAX_FIELD_SPECIES: usize = 128;
// Bound retained per-station thermochemical receipts, not just raw scalars.
const MAX_FIELD_STATE_COMPONENTS: usize = 16_384;

/// gas_field continuous STEPS_PER_SPAN STATIONS mole SPECIES_COUNT ID...
/// gas_field midpoint CELLS_PER_SEGMENT STATIONS mole SPECIES_COUNT ID...
/// gas_field CELLS_PER_SEGMENT STATIONS mole SPECIES_COUNT ID... (midpoint)
/// gas_station X_M T_K FRACTION... (one row per station, in declared ID order)
///
/// Ambient supplies the one common pressure. Species are sorted together with
/// their columns before entering the original canonical-composition owner.
pub(super) fn read(
    reader: &mut Reader<'_>,
    sources: Option<&ThermochemicalSources<'_>>,
    pressure_pa: f64,
    records: &mut Vec<GasRecord>,
) -> Result<(AxialGasField, GasFieldResolution), ReedPerformanceError> {
    let mut row = reader.row("gas_field")?;
    let line = row.line;
    let resolution = match row.words.clone().next() {
        Some("continuous") => {
            row.word()?;
            GasFieldResolution::Continuous(fs_duct::RiccatiSteps(
                row.count(fs_duct::profile::MAX_RICCATI_STEPS)?,
            ))
        }
        mode => {
            if mode == Some("midpoint") { row.word()?; }
            GasFieldResolution::Midpoint(
                row.count(crate::driving_point::regional::MAX_REGIONAL_SEGMENTS)?,
            )
        }
    };
    let stations = row.count(MAX_FIELD_STATIONS)?;
    if matches!(resolution, GasFieldResolution::Midpoint(0)
        | GasFieldResolution::Continuous(fs_duct::RiccatiSteps(0))) || stations < 2 {
        return Err(bad(line, "gas_field needs positive resolution and at least two stations"));
    }
    if row.word()? != "mole" {
        return Err(bad(line, "gas_field station compositions require explicit mole fractions"));
    }
    let count = row.count(MAX_FIELD_SPECIES)?;
    if count == 0 {
        return Err(bad(line, "gas_field needs at least one explicitly sourced species"));
    }
    if stations.checked_mul(count).is_none_or(|terms| terms > MAX_FIELD_STATE_COMPONENTS) {
        return Err(bad(line, "gas_field station/species product exceeds the retained source budget"));
    }
    let sources = sources.ok_or_else(|| bad(line, "gas_field requires an explicit thermochemical source registry"))?;
    let mut selected: Vec<(usize, SpeciesPackRef<'_>)> = Vec::with_capacity(count);
    for column in 0..count {
        let name = row.word()?;
        let mut candidates = sources.species.iter()
            .filter(|source| source.species.association().species().as_str() == name);
        let source = candidates.next()
            .ok_or_else(|| bad(line, "gas_field species has no source registry entry"))?;
        if candidates.next().is_some() {
            return Err(bad(line, "gas_field species has an ambiguous source registry entry"));
        }
        if selected.iter().any(|(_, entry)| entry.species.association().species().as_str() == name) {
            return Err(bad(line, "gas_field contains a duplicate species column"));
        }
        selected.push((column, (*source).clone()));
    }
    row.finish()?;
    selected.sort_by(|(_, a), (_, b)| {
        a.species.association().species().cmp(b.species.association().species())
    });

    let mut positions = Vec::with_capacity(stations);
    let mut states = Vec::with_capacity(stations);
    let mut station_lines = Vec::with_capacity(stations);
    for _ in 0..stations {
        let mut row = reader.row("gas_station")?;
        let station_line = row.line;
        let position = row.scalar()?;
        let temperature_k = row.positive()?;
        let mut supplied = Vec::with_capacity(count);
        for _ in 0..count {
            supplied.push(row.scalar()?);
        }
        row.finish()?;
        let fractions = selected.iter().map(|(column, _)| supplied[*column]).collect();
        let composition = Composition::new(CompositionBasis::MoleFraction, fractions)
            .map_err(|error| ReedPerformanceError::Gas {
                line: station_line,
                source: fs_material::MaterialError::Parameters {
                    what: format!("gas_field station composition: {error}"),
                },
            })?;
        positions.push(position);
        states.push(PrimitiveGasState { temperature_k, composition });
        station_lines.push(station_line);
    }

    // from_packs binds an active mixture. A positive, exactly summing dyadic
    // composition retains the full DECLARED species catalogue for later
    // with_composition calls. This seed is never evaluated at any T, p or x:
    // every physical sample uses the untouched authored station composition.
    // Exact zeros, including a species absent everywhere, remain exact zeros.
    let weight = 1.0 / count.next_power_of_two() as f64;
    let mut seed = vec![weight; count];
    seed[0] = 1.0 - weight * (count - 1) as f64;
    let composition = Composition::new(CompositionBasis::MoleFraction, seed)
        .map_err(|error| ReedPerformanceError::Gas {
            line,
            source: fs_material::MaterialError::Parameters {
                what: format!("gas_field catalogue composition: {error}"),
            },
        })?;
    let selected: Vec<_> = selected.into_iter().map(|(_, source)| source).collect();
    let water = sources.water_species.filter(|id| selected.iter()
        .any(|source| source.species.association().species() == *id));
    let model = from_packs(&selected, composition, water)
        .map_err(|source| ReedPerformanceError::Gas { line, source })?;
    let field = AxialGasField::new(
        positions, pressure_pa, states, model, QueryPoint::new(), sources.phase,
    ).map_err(|source| ReedPerformanceError::Gas { line, source })?;

    // Retain the authored knot receipts on the same gas-record surface used by
    // CLI/ensemble provenance. The field also retains the complete source set,
    // including a species which is inactive at one or every station.
    for (&position, &line) in field.stations_m().iter().zip(&station_lines) {
        let state = field.state_at(position)
            .map_err(|source| ReedPerformanceError::Gas { line, source })?;
        records.push(GasRecord::Thermochemical(state));
    }
    Ok((field, resolution))
}

#[cfg(test)]
mod tests;
