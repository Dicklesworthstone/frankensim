//! Cold material admission for the existing CHT solver. The canonical gas,
//! liquid and material-state owners supply every coefficient; this module only
//! reads files and retains the exact reference states used by the CLI.

use std::path::Path;

use fs_lbm::material::conjugate::MaterialChtFluid;
use fs_matdb::{NormalizedMaterialCardPack, QueryPoint};
use fs_material::gas_species::input::packs::{
    MAX_RECORD_BYTES, MAX_SOURCE_BYTES, MAX_SOURCE_PACK_BYTES, SourcePackBytes, read_source_packs,
};
use fs_material::gas_species::input::{GasRecord, ThermochemicalPhase, read_gas_record};
use fs_material::liquid::{
    ResolvedLiquidState, resolve_liquid_state, resolve_liquid_state_from_kinematic,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, ResolvedMaterialStatePoint,
    SPECIFIC_HEAT_CAPACITY_DIMS, SPECIFIC_HEAT_CAPACITY_PROPERTY, ScalarAdmissibility,
    ScalarPropertyRequirement, THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY,
    resolve_material_state_point,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, QuantitySpec};

use super::{FluidProperties, J, Result, SolidMaterial, bad, field, quote};

pub(super) enum FluidSource {
    Gas {
        record: GasRecord,
        declaration: String,
    },
    Liquid {
        state: ResolvedLiquidState,
        pack: String,
    },
}

pub(super) struct SolidSource {
    pub index: usize,
    state: ResolvedMaterialStatePoint,
    pack: String,
}

fn fields(value: &J, allowed: &[&str]) -> Result<()> {
    let object = value
        .as_object()
        .ok_or_else(|| bad("material source must be an object"))?;
    if let Some((key, _)) = object
        .iter()
        .find(|(key, _)| !allowed.contains(&key.as_str()))
    {
        return Err(bad(format!("unsupported material-source field {key:?}")));
    }
    Ok(())
}

fn string<'a>(value: &'a J, key: &str) -> Result<&'a str> {
    field(value, key, "material source")?
        .as_str()
        .filter(|s| !s.trim().is_empty())
        .ok_or_else(|| bad(format!("material source {key} must be a nonempty string")))
}

fn frozen(value: &J) -> Result<()> {
    if string(value, "mode")? != "frozen" {
        return Err(bad(
            "cooling-cht material sources require explicit mode frozen; evolving properties need the coupled material solver",
        ));
    }
    Ok(())
}

fn bytes(base: &Path, path: &str, total: &mut usize) -> Result<Vec<u8>> {
    use std::io::Read as _;
    let limit = MAX_SOURCE_PACK_BYTES.min(MAX_SOURCE_BYTES.saturating_sub(*total));
    let mut bytes = Vec::new();
    std::fs::File::open(base.join(path))
        .and_then(|file| file.take(limit as u64 + 1).read_to_end(&mut bytes))
        .map_err(|error| bad(format!("material source {path:?}: {error}")))?;
    if bytes.is_empty() || bytes.len() > limit {
        return Err(bad(
            "material sources require nonempty packs within 16 MiB per pack and 64 MiB per scene",
        ));
    }
    *total += bytes.len();
    Ok(bytes)
}

fn load_card(
    value: &J,
    base: &Path,
    total: &mut usize,
) -> Result<(NormalizedMaterialCardPack, QueryPoint)> {
    let pack = NormalizedMaterialCardPack::from_bytes(&bytes(base, string(value, "pack")?, total)?)
        .map_err(|error| bad(format!("material-card pack: {error}")))?;
    // Every physics command preserves the same explicit source coordinates;
    // the cold loader does not serialize/reparse a second command's JSON type.
    let point = crate::material_state_input::state_point(value)
        .map_err(bad)?;
    Ok((pack, point))
}

pub(super) fn fluid(
    value: &J,
    base: &Path,
    total: &mut usize,
) -> Result<(FluidProperties, FluidSource)> {
    frozen(value)?;
    if value.get("gas").is_some() {
        fields(value, &["mode", "gas", "sources", "phase", "water_species"])?;
        let declaration = string(value, "gas")?;
        if declaration.len() > MAX_RECORD_BYTES {
            return Err(bad("gas declaration exceeds 8 KiB"));
        }
        let record = match value.get("sources") {
            Some(sources) => {
                let sources = sources
                    .as_array()
                    .filter(|s| !s.is_empty() && s.len() <= 128)
                    .ok_or_else(|| {
                        bad("gas sources require 1..128 species/model/transport triples")
                    })?;
                let phase = match string(value, "phase")? {
                    "frozen-gas" => ThermochemicalPhase::FrozenGas,
                    "liquid-water-limit" => ThermochemicalPhase::LiquidWaterLimit,
                    _ => return Err(bad("gas phase must be frozen-gas or liquid-water-limit")),
                };
                let water = value
                    .get("water_species")
                    .map(|_| string(value, "water_species"))
                    .transpose()?;
                let mut loaded = Vec::with_capacity(sources.len());
                for source in sources {
                    fields(source, &["species", "nasa9", "transport"])?;
                    loaded.push([
                        bytes(base, string(source, "species")?, total)?,
                        bytes(base, string(source, "nasa9")?, total)?,
                        bytes(base, string(source, "transport")?, total)?,
                    ]);
                }
                let references = loaded
                    .iter()
                    .map(|s| SourcePackBytes {
                        species: &s[0],
                        nasa9: &s[1],
                        transport: &s[2],
                    })
                    .collect::<Vec<_>>();
                read_source_packs(declaration, &references, water, phase)
                    .map_err(|error| bad(error.to_string()))?
            }
            None => {
                if value.get("phase").is_some() || value.get("water_species").is_some() {
                    return Err(bad(
                        "phase and water_species belong to explicit thermochemical source packs",
                    ));
                }
                let mut fields = declaration.split_ascii_whitespace();
                let record = read_gas_record(&mut fields, None, None)
                    .map_err(|error| bad(error.to_string()))?;
                if fields.next().is_some() {
                    return Err(bad("unexpected fields after gas declaration"));
                }
                record
            }
        };
        let properties = MaterialChtFluid::try_new(record.state())
            .map_err(|error| bad(error.to_string()))?
            .properties()
            .clone();
        Ok((
            properties,
            FluidSource::Gas {
                record,
                declaration: declaration.into(),
            },
        ))
    } else {
        fields(value, &["mode", "pack", "state", "viscosity"])?;
        let (pack, point) = load_card(value, base, total)?;
        let state = match string(value, "viscosity")? {
            "dynamic" => resolve_liquid_state(
                pack.card(),
                &point,
                MaterialPropertySelection::SingleClaimOnly,
            ),
            "kinematic" => resolve_liquid_state_from_kinematic(
                pack.card(),
                &point,
                MaterialPropertySelection::SingleClaimOnly,
            ),
            _ => {
                return Err(bad(
                    "liquid source viscosity must explicitly be dynamic or kinematic",
                ));
            }
        }
        .map_err(|error| bad(error.to_string()))?;
        admit_static(state.material())?;
        let properties = MaterialChtFluid::try_new(&state)
            .map_err(|error| bad(error.to_string()))?
            .properties()
            .clone();
        Ok((
            properties,
            FluidSource::Liquid {
                state,
                pack: pack.content_hash().to_hex(),
            },
        ))
    }
}

pub(super) fn solid(
    value: &J,
    base: &Path,
    total: &mut usize,
    name: &str,
    index: usize,
) -> Result<(SolidMaterial, SolidSource)> {
    fields(value, &["mode", "pack", "state"])?;
    frozen(value)?;
    let (pack, point) = load_card(value, base, total)?;
    let positive = ScalarAdmissibility::StrictlyPositive;
    let conductivity_quantity = if pack
        .card()
        .claims()
        .registered_quantity(THERMAL_CONDUCTIVITY_PROPERTY)
        .is_some_and(|q| q.semantic_type().is_some())
    {
        QuantitySpec::semantic(SemanticType::new(
            QuantityKind::ThermalConductivity,
            ValueForm::Static,
        ))
    } else {
        QuantitySpec::dimensional(THERMAL_CONDUCTIVITY_DIMS)
    };
    let requirements = [
        ScalarPropertyRequirement::try_new(DENSITY_PROPERTY, Density::DIMS, positive),
        ScalarPropertyRequirement::try_new(
            SPECIFIC_HEAT_CAPACITY_PROPERTY,
            SPECIFIC_HEAT_CAPACITY_DIMS,
            positive,
        ),
        ScalarPropertyRequirement::try_with_quantity(
            THERMAL_CONDUCTIVITY_PROPERTY,
            conductivity_quantity,
            positive,
        ),
    ]
    .into_iter()
    .collect::<std::result::Result<Vec<_>, _>>()
    .map_err(|error| bad(error.to_string()))?;
    let state = resolve_material_state_point(
        pack.card(),
        &point,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .map_err(|error| bad(error.to_string()))?;
    admit_static(&state)?;
    let coefficient = |name| {
        state
            .property(name)
            .expect("required source property resolved")
            .value_si()
    };
    let capacity = coefficient(DENSITY_PROPERTY) * coefficient(SPECIFIC_HEAT_CAPACITY_PROPERTY);
    if !capacity.is_finite() || capacity <= 0.0 {
        return Err(bad("sourced volumetric heat capacity is unrepresentable"));
    }
    let material = SolidMaterial::new(name, coefficient(THERMAL_CONDUCTIVITY_PROPERTY))
        .with_heat_capacity(capacity);
    Ok((
        material,
        SolidSource {
            index,
            state,
            pack: pack.content_hash().to_hex(),
        },
    ))
}

fn admit_static(state: &ResolvedMaterialStatePoint) -> Result<()> {
    if state
        .properties()
        .iter()
        .any(|property| !property.is_static_at_fixed_external_conditions())
    {
        return Err(bad(
            "frozen CHT requires static material coefficients at fixed thermal, pressure and composition conditions; frequency, field, strain and shear-rate sources need their physical owner",
        ));
    }
    Ok(())
}

impl FluidSource {
    /// A liquid must never inherit the legacy ideal-gas expansion coefficient.
    pub(super) fn requires_expansion(&self) -> bool {
        matches!(self, Self::Liquid { .. })
    }

    pub(super) fn gas_reference_temperature(&self) -> Option<f64> {
        match self {
            Self::Gas { record, .. } => Some(record.state().temperature),
            Self::Liquid { .. } => None,
        }
    }

    pub(super) fn json(&self) -> String {
        match self {
            Self::Gas {
                record,
                declaration,
            } => {
                let gas = record.state();
                let packs = match record {
                    GasRecord::Thermochemical(state) => state.source_packs().iter().map(|source|
                        format!("{{\"species\":{},\"species_pack\":{},\"nasa9_pack\":{},\"transport_pack\":{}}}",
                            quote(source.species.as_str()), quote(&source.species_pack.to_hex()),
                            quote(&source.nasa9_pack.to_hex()), quote(&source.transport_pack.to_hex())))
                        .collect::<Vec<_>>().join(","),
                    _ => String::new(),
                };
                format!(
                    "{{\"kind\":\"gas\",\"mode\":\"frozen\",\"declaration\":{},\"temperature_k\":{},\"pressure_pa\":{},\"sound_speed_m_s\":{},\"packs\":[{}]}}",
                    quote(declaration),
                    gas.temperature,
                    gas.pressure,
                    gas.sound_speed,
                    packs
                )
            }
            Self::Liquid { state, pack } => format!(
                "{{\"kind\":\"liquid\",\"mode\":\"frozen\",\"pack_identity\":{},\"material_state\":{}}}",
                quote(pack),
                quote(&state.material().identity().to_hex())
            ),
        }
    }
}

impl SolidSource {
    pub(super) fn json(&self) -> String {
        let value = |name| {
            self.state
                .property(name)
                .expect("required source property resolved")
                .value_si()
        };
        format!(
            "{{\"index\":{},\"mode\":\"frozen\",\"pack_identity\":{},\"material_state\":{},\"conductivity_w_m_k\":{},\"density_kg_m3\":{},\"specific_heat_j_kg_k\":{}}}",
            self.index,
            quote(&self.pack),
            quote(&self.state.identity().to_hex()),
            value(THERMAL_CONDUCTIVITY_PROPERTY),
            value(DENSITY_PROPERTY),
            value(SPECIFIC_HEAT_CAPACITY_PROPERTY)
        )
    }
}
