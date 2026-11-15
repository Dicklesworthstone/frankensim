//! `gas_cavity`: a closed, uniform gas of cited thermochemistry, `fs-gas`'s
//! [`GasCavity`], with a pressure-work port and, unless insulated, a heat
//! port.
//!
//! The gas is a thermochemical record resolved through explicit canonical
//! species, NASA-9 and transport packs, by the same loader the cht command
//! uses. The cavity owns the calorics, the pressure and the discrete-gradient
//! heat and work solve; nothing here evaluates the gas or integrates anything.

use std::path::Path;

use fs_gas::component::{ClosedGasParcel, GasCavity, GasEnergySolve};
use fs_material::gas_species::input::ThermochemicalPhase;
use fs_material::gas_species::input::packs::{
    MAX_SOURCE_BYTES, MAX_SOURCE_PACK_BYTES, SourcePackBytes, read_source_parcel,
};
use fs_ports::component::PortSpec;
use fs_system::{NodeId, System, SystemBuilder};

use super::materials::read_pack;
use super::{J, Refusal, add, fields, number_or, positive};

/// The most species a record may name: the frozen mixture's own bound.
const MAX_SOURCES: usize = 128;

/// The owner's caloric solve budget, as the `fs-gas` component tests declare
/// it: the energy residual per kilogram and the safeguarded iterations.
const ENERGY_TOLERANCE_J_KG: f64 = 1.0e-7;
const MAX_ITERATIONS: f64 = 128.0;

/// A `gas_cavity` element and its initial temperature [K], the starting
/// value of the temperatures joined to its heat port.
pub(super) fn cavity(
    builder: &mut SystemBuilder,
    name: &str,
    definition: &J,
    scene: &Path,
    what: &str,
) -> Result<((NodeId, Vec<PortSpec>), f64), Refusal> {
    fields(
        definition,
        &[
            "kind",
            "gas",
            "sources",
            "volume_m3",
            "temperature_min_k",
            "temperature_max_k",
            "energy_tolerance_j_kg",
            "max_iterations",
            "insulated",
        ],
        what,
    )?;
    let record = definition
        .str_field("gas")
        .ok_or_else(|| format!("{what}: gas is a thermochemical record"))?;
    let rows = definition
        .get("sources")
        .and_then(J::as_array)
        .filter(|rows| !rows.is_empty() && rows.len() <= MAX_SOURCES)
        .ok_or_else(|| {
            format!(
                "{what}: sources is an array of 1 to {MAX_SOURCES} species, nasa9 and \
                 transport pack triples"
            )
        })?;
    let base = scene.parent().unwrap_or(Path::new("."));
    let mut total = 0;
    let mut loaded = Vec::with_capacity(rows.len());
    for (index, row) in rows.iter().enumerate() {
        let inner = format!("{what}.sources[{index}]");
        fields(row, &["species", "nasa9", "transport"], &inner)?;
        let mut triple = Vec::with_capacity(3);
        for key in ["species", "nasa9", "transport"] {
            let path = row
                .str_field(key)
                .filter(|path| !path.trim().is_empty())
                .ok_or_else(|| format!("{inner}.{key}: a pack path, relative to the scene"))?;
            let limit = MAX_SOURCE_PACK_BYTES.min(MAX_SOURCE_BYTES.saturating_sub(total));
            let bytes = read_pack(&base.join(path), limit as u64, &format!("{inner}.{key}"))?;
            total += bytes.len();
            triple.push(bytes);
        }
        loaded.push(triple);
    }
    let sources: Vec<SourcePackBytes<'_>> = loaded
        .iter()
        .map(|triple| SourcePackBytes {
            species: &triple[0],
            nasa9: &triple[1],
            transport: &triple[2],
        })
        .collect();
    let volume_m3 = positive(definition, "volume_m3", what)?;
    let parcel = read_source_parcel(
        record,
        &sources,
        None,
        ThermochemicalPhase::FrozenGas,
        volume_m3,
    )
    .map_err(|error| format!("{what}.gas: {error}"))?;
    let iterations = number_or(definition, "max_iterations", MAX_ITERATIONS, what)?;
    if !(1.0..=MAX_ITERATIONS).contains(&iterations) || iterations.fract() != 0.0 {
        return Err(format!("{what}: max_iterations is a whole number from 1 to 128").into());
    }
    // A whole number from 1 to 128, checked above.
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let maximum_iterations = iterations as usize;
    let controls = GasEnergySolve {
        minimum_temperature_k: positive(definition, "temperature_min_k", what)?,
        maximum_temperature_k: positive(definition, "temperature_max_k", what)?,
        energy_tolerance_j_kg: number_or(
            definition,
            "energy_tolerance_j_kg",
            ENERGY_TOLERANCE_J_KG,
            what,
        )?,
        maximum_iterations,
    };
    let insulated = match definition.get("insulated") {
        None => false,
        Some(J::Bool(insulated)) => *insulated,
        Some(_) => return Err(format!("{what}: insulated is true or false").into()),
    };
    let reference = ClosedGasParcel::from_isochoric(parcel);
    let temperature_k = reference.state().state().temperature;
    let mut cavity =
        GasCavity::new(reference.clone(), controls).map_err(|error| format!("{what}: {error}"))?;
    if insulated {
        cavity = cavity
            .insulated()
            .map_err(|error| format!("{what}: {error}"))?;
    }
    Ok((add(builder, name, cavity, reference)?, temperature_k))
}

/// The report of a gas cavity: the state of the gas it holds.
pub(super) fn values(system: &System, node: NodeId) -> Vec<(String, f64)> {
    let Some(gas) = system.state::<ClosedGasParcel>(node) else {
        return Vec::new();
    };
    let state = gas.state();
    [
        ("temperature_k", state.state().temperature),
        ("pressure_pa", state.state().pressure),
        ("density_kg_m3", state.state().density),
        ("volume_m3", gas.volume_m3()),
        ("mass_kg", gas.mass_kg()),
        (
            "specific_internal_energy_j_kg",
            state.specific_internal_energy_j_kg(),
        ),
    ]
    .into_iter()
    .map(|(key, value)| (key.to_owned(), value))
    .collect()
}
