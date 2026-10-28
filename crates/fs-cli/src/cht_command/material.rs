//! Material admission and source-preserving CHT feedback. The canonical gas,
//! liquid and material-state owners supply every coefficient; this module
//! retains the original cards and invokes the existing coupled solver.

use std::path::Path;

use fs_lbm::material::conjugate::MaterialChtFluid;
use fs_lbm::material::conjugate::solid::SolidConductivityLaw;
use fs_matdb::{NormalizedMaterialCardPack, QueryPoint};
use fs_material::gas_species::input::packs::{
    MAX_RECORD_BYTES, MAX_SOURCE_BYTES, MAX_SOURCE_PACK_BYTES, SourcePackBytes, read_source_packs,
};
use fs_material::gas_species::input::{GasRecord, ThermochemicalPhase, read_gas_record};
use fs_material::liquid::{
    ResolvedLiquidState, resolve_liquid_state, resolve_liquid_state_from_kinematic,
};
use fs_material::liquid_temperature::{LiquidTemperatureLaw, LiquidViscosityBasis};
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
        pack: NormalizedMaterialCardPack,
        temperature_axis: Option<String>,
        viscosity_basis: LiquidViscosityBasis,
    },
}

pub(super) struct SolidSource {
    pub index: usize,
    state: ResolvedMaterialStatePoint,
    pack: NormalizedMaterialCardPack,
    temperature_axis: Option<String>,
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

fn temperature_axis(value: &J) -> Result<Option<String>> {
    match string(value, "mode")? {
        "frozen" if value.get("temperature_axis").is_none() => Ok(None),
        "temperature-dependent" => Ok(Some(string(value, "temperature_axis")?.into())),
        "frozen" => Err(bad("temperature_axis belongs to temperature-dependent material mode")),
        _ => Err(bad("material mode must explicitly be frozen or temperature-dependent")),
    }
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
    let temperature_axis = temperature_axis(value)?;
    if value.get("gas").is_some() {
        if temperature_axis.is_some() {
            return Err(bad("temperature-dependent cooling-cht requires a source-backed liquid; gas sources currently support frozen mode"));
        }
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
        fields(value, &["mode", "pack", "state", "viscosity", "temperature_axis"])?;
        let (pack, point) = load_card(value, base, total)?;
        let viscosity_basis = match string(value, "viscosity")? {
            "dynamic" => LiquidViscosityBasis::Dynamic,
            "kinematic" => LiquidViscosityBasis::Kinematic,
            _ => return Err(bad("liquid source viscosity must explicitly be dynamic or kinematic")),
        };
        let state = match viscosity_basis {
            LiquidViscosityBasis::Dynamic => resolve_liquid_state(
                pack.card(),
                &point,
                MaterialPropertySelection::SingleClaimOnly,
            ),
            LiquidViscosityBasis::Kinematic => resolve_liquid_state_from_kinematic(
                pack.card(),
                &point,
                MaterialPropertySelection::SingleClaimOnly,
            ),
        }
        .map_err(|error| bad(error.to_string()))?;
        admit_static(state.material())?;
        if let Some(axis) = &temperature_axis {
            LiquidTemperatureLaw::try_from_material(
                pack.card(), state.material().clone(), axis, viscosity_basis,
            ).map_err(|error| bad(error.to_string()))?;
        }
        let properties = MaterialChtFluid::try_new(&state)
            .map_err(|error| bad(error.to_string()))?
            .properties()
            .clone();
        Ok((
            properties,
            FluidSource::Liquid {
                state,
                pack,
                temperature_axis,
                viscosity_basis,
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
    fields(value, &["mode", "pack", "state", "temperature_axis"])?;
    let temperature_axis = temperature_axis(value)?;
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
    if let Some(axis) = &temperature_axis {
        SolidConductivityLaw::try_new(
            pack.card(), state.clone(), axis, [THERMAL_CONDUCTIVITY_PROPERTY; 3],
        ).map_err(|error| bad(error.to_string()))?;
    }
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
            pack,
            temperature_axis,
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
            "CHT requires static material coefficients at fixed thermal, pressure and composition conditions; frequency, field, strain and shear-rate sources need their physical owner",
        ));
    }
    Ok(())
}

impl FluidSource {
    /// A liquid must never inherit the legacy ideal-gas expansion coefficient.
    pub(super) fn requires_expansion(&self) -> bool {
        matches!(self, Self::Liquid { temperature_axis: None, .. })
    }

    pub(super) fn temperature_dependent(&self) -> bool {
        matches!(self, Self::Liquid { temperature_axis: Some(_), .. })
    }

    pub(super) fn evolving_reference_temperature(&self) -> Option<f64> {
        match self {
            Self::Liquid { state, temperature_axis: Some(axis), .. } => {
                Some(state.material().query().axes()[axis])
            }
            _ => None,
        }
    }

    fn law(&self) -> Result<LiquidTemperatureLaw<'_>> {
        let Self::Liquid { state, pack, temperature_axis: Some(axis), viscosity_basis } = self else {
            return Err(bad("temperature-dependent CHT needs an explicitly evolving liquid source"));
        };
        LiquidTemperatureLaw::try_from_material(
            pack.card(), state.material().clone(), axis, *viscosity_basis,
        ).map_err(|error| bad(error.to_string()))
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
            Self::Liquid { state, pack, temperature_axis, .. } => format!(
                "{{\"kind\":\"liquid\",{},\"pack_identity\":{},\"material_state\":{}}}",
                mode_json(temperature_axis.as_deref()),
                quote(&pack.content_hash().to_hex()),
                quote(&state.material().identity().to_hex()),
            ),
        }
    }
}

impl SolidSource {
    pub(super) fn temperature_dependent(&self) -> bool {
        self.temperature_axis.is_some()
    }

    fn law(&self) -> Result<SolidConductivityLaw<'_>> {
        let axis = self.temperature_axis.as_deref()
            .ok_or_else(|| bad("temperature-dependent CHT needs every solid source in the same mode"))?;
        SolidConductivityLaw::try_new(
            self.pack.card(), self.state.clone(), axis, [THERMAL_CONDUCTIVITY_PROPERTY; 3],
        ).map_err(|error| bad(error.to_string()))
    }

    pub(super) fn json(&self) -> String {
        let value = |name| {
            self.state
                .property(name)
                .expect("required source property resolved")
                .value_si()
        };
        format!(
            "{{\"index\":{},{},\"pack_identity\":{},\"material_state\":{},\"conductivity_w_m_k\":{},\"density_kg_m3\":{},\"specific_heat_j_kg_k\":{}}}",
            self.index,
            mode_json(self.temperature_axis.as_deref()),
            quote(&self.pack.content_hash().to_hex()),
            quote(&self.state.identity().to_hex()),
            value(THERMAL_CONDUCTIVITY_PROPERTY),
            value(DENSITY_PROPERTY),
            value(SPECIFIC_HEAT_CAPACITY_PROPERTY)
        )
    }
}

fn mode_json(axis: Option<&str>) -> String {
    axis.map_or_else(|| "\"mode\":\"frozen\"".into(), |axis| {
        format!("\"mode\":\"temperature-dependent\",\"temperature_axis\":{}", quote(axis))
    })
}

pub(super) const EVOLVING_NO_CLAIM: &str = "steady laminar fixed-geometry Boussinesq flow and conservative conjugate heat transfer; original liquid density, viscosity and conductivity and solid conductivity reevaluated at local temperature with source claims pinned; reference liquid density and heat capacity retained in inertia, continuity and heat advection with at most 5 percent departures; isotropic source solids; no compressibility, solid storage, melting, deformation, radiation, turbulence closure or phase-stability claim; staircase geometry and one declared resolution unless grid_convergence is requested; Estimated numerical evidence, not experimental validation";

/// Existing CLI geometry and boundaries, existing coupled momentum/heat owner.
/// Every runtime query borrows the same owned cards admitted by the loader.
pub(super) fn solve(
    scene: &super::Scene,
    domain: &fs_lbm::conjugate::VoxelDomain,
    flow: &fs_lbm::conjugate::SimpleConfig,
    thermal: &fs_lbm::conjugate::ThermalSetup,
    gate: &fs_exec::CancelGate,
) -> Result<(fs_lbm::conjugate::FvFlow, fs_lbm::conjugate::EnergySolution, usize, String)> {
    use fs_lbm::material::conjugate::buoyant::{MaterialBuoyancyConfig, MaterialBuoyancyError};
    use fs_lbm::material::transport::boussinesq::BoussinesqConfig;
    use fs_material::MaterialError;

    let liquid = scene.fluid_source.as_ref()
        .ok_or_else(|| bad("temperature-dependent CHT requires a liquid source"))?.law()?;
    let solids = scene.solid_sources.iter().map(SolidSource::law).collect::<Result<Vec<_>>>()?;
    let reference_temperature = liquid.reference_point().axes()[liquid.temperature_axis()];
    let solid_initial: Vec<_> = scene.solid_sources.iter().map(|source| {
        let axis = source.temperature_axis.as_ref().expect("evolving scene admitted");
        source.state.query().axes()[axis]
    }).collect();
    let mut initial = vec![reference_temperature; domain.cell_count()];
    for (cell, temperature) in initial.iter_mut().enumerate() {
        if gate.is_requested() {
            return Err(super::solver_failure(&fs_lbm::conjugate::ChtError::Cancelled));
        }
        if let fs_lbm::conjugate::Voxel::Solid(index) = domain.voxel_at(cell) {
            *temperature = solid_initial[usize::from(index)];
        }
    }
    let bound = MaterialChtFluid::try_new(liquid.reference_state())
        .map_err(|error| bad(error.to_string()))?;
    let relative = scene.tolerance.min(1e-6);
    let gravity = scene.gravity.unwrap_or([0.0; 3]);
    let force_scale = scene.fluid.density_kg_m3
        * gravity.into_iter().map(f64::abs).fold(1.0, f64::max);
    let controls = MaterialBuoyancyConfig {
        boussinesq: BoussinesqConfig {
            gravity_m_s2: gravity,
            max_relative_density_departure: 0.05,
            max_relative_specific_heat_departure: 0.05,
        },
        maximum_passes: scene.max_iterations,
        temperature_tolerance_k: relative * reference_temperature,
        transport_relative_tolerance: relative,
        force_tolerance_n_m3: relative * force_scale,
        temperature_relaxation: 0.5,
    };
    let run = bound.solve_natural_convection_with_material_solids(
        domain, &solids, flow, thermal, &fs_lbm::conjugate::EnergyConfig::default(),
        &initial, controls, gate, |_, temperature, _| {
            liquid.at_temperature(fs_qty::Temperature::new(temperature))
                .map_err(|error| MaterialError::State { what: error.to_string() })
        },
    ).map_err(|error| match &error {
        MaterialBuoyancyError::Solver(source) => super::solver_failure(source),
        MaterialBuoyancyError::NotConverged(_) => super::Failure {
            code: "cooling-cht-not-steady", message: error.to_string(),
        },
        _ => bad(error.to_string()),
    })?;
    let fluid = endpoint_json(run.endpoint.iter().map(|sample| (
        sample.cell, sample.temperature_k,
        sample.material.resolved_transport().thermal_conductivity_w_m_k(),
        sample.material.material().material(),
    )), gate)?;
    let solid = endpoint_json(run.solid_endpoint.iter().map(|sample| (
        sample.cell, sample.temperature_k, sample.conductivity_w_m_k[0], &sample.material,
    )), gate)?;
    let report = run.report;
    let evidence = format!(concat!(
        "{{\"mode\":\"temperature-dependent\",\"passes\":{},",
        "\"temperature_defect_k\":{},\"transport_relative_defect\":{},\"force_defect_n_m3\":{},",
        "\"temperature_tolerance_k\":{},\"transport_relative_tolerance\":{},\"force_tolerance_n_m3\":{},",
        "\"max_relative_density_departure\":0.05,\"max_relative_specific_heat_departure\":0.05,\"solid_to_fluid_heat_w\":{},",
        "\"fluid_endpoint\":{},\"solid_endpoint\":{}}}"),
        report.passes, super::num(report.temperature_defect_k)?,
        super::num(report.transport_relative_defect)?, super::num(report.force_defect_n_m3)?,
        super::num(controls.temperature_tolerance_k)?, super::num(controls.transport_relative_tolerance)?,
        super::num(controls.force_tolerance_n_m3)?,
        super::num(run.energy.solid_to_fluid_heat_w(domain))?, fluid, solid);
    Ok((run.flow, run.energy, report.passes, evidence))
}

/// Compact identities of the actual endpoint states, in original cell order.
/// The digest includes each complete state identity (source choices and query),
/// rather than claiming that one representative temperature describes a body.
fn endpoint_json<'a>(
    samples: impl Iterator<Item = (usize, f64, f64, &'a ResolvedMaterialStatePoint)>,
    gate: &fs_exec::CancelGate,
) -> Result<String> {
    let mut digest = fs_blake3::DomainHasher::new("frankensim.cooling-cht.material-endpoint.v1");
    let mut count = 0;
    let (mut t_min, mut t_max, mut k_min, mut k_max) =
        (f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY);
    for (cell, temperature, conductivity, state) in samples {
        if gate.is_requested() {
            return Err(super::solver_failure(&fs_lbm::conjugate::ChtError::Cancelled));
        }
        count += 1;
        digest.update(&(cell as u64).to_le_bytes());
        digest.update(state.identity().as_bytes());
        t_min = t_min.min(temperature); t_max = t_max.max(temperature);
        k_min = k_min.min(conductivity); k_max = k_max.max(conductivity);
    }
    if count == 0 { return Ok("{\"cells\":0}".into()); }
    Ok(format!(concat!(
        "{{\"cells\":{},\"states_identity\":{},\"minimum_temperature_k\":{},\"maximum_temperature_k\":{},",
        "\"minimum_conductivity_w_m_k\":{},\"maximum_conductivity_w_m_k\":{}}}"),
        count, quote(&digest.finalize().to_hex()), super::num(t_min)?, super::num(t_max)?,
        super::num(k_min)?, super::num(k_max)?))
}
