//! Material admission and source-preserving CHT feedback. The canonical gas,
//! liquid and material-state owners supply every coefficient; this module
//! retains the original cards and invokes the existing coupled solver.

pub(super) mod enthalpy;

use std::path::Path;

use fs_lbm::material::conjugate::MaterialChtFluid;
use fs_lbm::material::conjugate::solid::SolidConductivityLaw;
use fs_matdb::{NormalizedMaterialCardPack, QueryPoint};
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_material::gas_species::input::packs::{
    MAX_RECORD_BYTES, MAX_SOURCE_BYTES, MAX_SOURCE_PACK_BYTES, SourcePackBytes, read_source_packs,
    read_source_temperature_law,
};
use fs_material::gas_species::input::temperature::IsobaricGasLaw;
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
        temperature_law: Option<IsobaricGasLaw>,
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
    storage_temperature_range_k: Option<[f64; 2]>,
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
    if value.get("gas").is_some() {
        fields(value, &["mode", "gas", "sources", "phase", "water_species"])?;
        let evolving = match string(value, "mode")? {
            "frozen" => false,
            "temperature-dependent" => true,
            _ => return Err(bad("material mode must explicitly be frozen or temperature-dependent")),
        };
        let declaration = string(value, "gas")?;
        if declaration.len() > MAX_RECORD_BYTES {
            return Err(bad("gas declaration exceeds 8 KiB"));
        }
        let mut temperature_law = None;
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
                if evolving {
                    let law = read_source_temperature_law(declaration, &references, water, phase)
                        .map_err(|error| bad(error.to_string()))?;
                    let record = law.reference().clone();
                    temperature_law = Some(law);
                    record
                } else {
                    read_source_packs(declaration, &references, water, phase)
                        .map_err(|error| bad(error.to_string()))?
                }
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
                if evolving {
                    temperature_law = Some(IsobaricGasLaw::try_from_record(record.clone())
                        .map_err(|error| bad(error.to_string()))?);
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
                temperature_law,
            },
        ))
    } else {
        let temperature_axis = temperature_axis(value)?;
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
    fields(value, &[
        "mode", "pack", "state", "temperature_axis", "storage_temperature_range_k",
    ])?;
    let temperature_axis = temperature_axis(value)?;
    let storage_temperature_range_k = value.get("storage_temperature_range_k")
        .map(|range| {
            if temperature_axis.is_none() {
                return Err(bad("storage_temperature_range_k belongs to temperature-dependent solid sources"));
            }
            let pair = range.as_array().filter(|pair| pair.len() == 2)
                .ok_or_else(|| bad("storage_temperature_range_k must be [lower_kelvin, upper_kelvin]"))?;
            let lo = pair[0].as_f64();
            let hi = pair[1].as_f64();
            match (lo, hi) {
                (Some(lo), Some(hi)) if lo.is_finite() && hi.is_finite() && lo > 0.0 && lo < hi => {
                    Ok([lo, hi])
                }
                _ => Err(bad("storage_temperature_range_k needs positive, finite, increasing bounds")),
            }
        }).transpose()?;
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
            storage_temperature_range_k,
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
        matches!(
            self,
            Self::Liquid {
                temperature_axis: Some(_),
                ..
            } | Self::Gas {
                temperature_law: Some(_),
                ..
            }
        )
    }

    pub(super) fn evolving_reference_temperature(&self) -> Option<f64> {
        match self {
            Self::Gas { record, temperature_law: Some(_), .. } => Some(record.state().temperature),
            Self::Liquid { state, temperature_axis: Some(axis), .. } => {
                Some(state.material().query().axes()[axis])
            }
            _ => None,
        }
    }

    fn law(&self) -> Result<FluidLaw<'_>> {
        match self {
            Self::Gas {
                temperature_law: Some(law),
                ..
            } => Ok(FluidLaw::Gas(law)),
            Self::Liquid {
                state,
                pack,
                temperature_axis: Some(axis),
                viscosity_basis,
            } => LiquidTemperatureLaw::try_from_material(
                pack.card(),
                state.material().clone(),
                axis,
                *viscosity_basis,
            )
            .map(FluidLaw::Liquid)
            .map_err(|error| bad(error.to_string())),
            _ => Err(bad(
                "temperature-dependent CHT needs an explicitly evolving fluid source",
            )),
        }
    }

    pub(super) fn gas_reference_temperature(&self) -> Option<f64> {
        match self {
            Self::Gas { record, temperature_law: None, .. } => Some(record.state().temperature),
            _ => None,
        }
    }

    pub(super) fn json(&self) -> String {
        match self {
            Self::Gas {
                record,
                declaration,
                temperature_law,
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
                    "{{\"kind\":\"gas\",\"mode\":{},\"declaration\":{},\"temperature_k\":{},\"pressure_pa\":{},\"sound_speed_m_s\":{},\"packs\":[{}]}}",
                    quote(if temperature_law.is_some() { "temperature-dependent" } else { "frozen" }),
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

enum FluidLaw<'a> {
    Gas(&'a IsobaricGasLaw),
    Liquid(LiquidTemperatureLaw<'a>),
}

enum FluidState {
    Gas(GasRecord),
    Liquid(ResolvedLiquidState),
}

impl FluidTransport for FluidState {
    fn transport_state(
        &self,
    ) -> std::result::Result<FluidTransportState, fs_material::MaterialError> {
        match self {
            Self::Gas(gas) => gas.transport_state(),
            Self::Liquid(liquid) => liquid.transport_state(),
        }
    }
}

impl FluidLaw<'_> {
    fn reference(&self) -> &dyn FluidTransport {
        match self {
            Self::Gas(law) => law.reference(),
            Self::Liquid(law) => law.reference_state(),
        }
    }

    fn at_temperature(
        &self,
        temperature_k: f64,
    ) -> std::result::Result<FluidState, fs_material::MaterialError> {
        let temperature = fs_qty::Temperature::new(temperature_k);
        match self {
            Self::Gas(law) => law.at_temperature(temperature).map(FluidState::Gas),
            Self::Liquid(law) => law
                .at_temperature(temperature)
                .map(FluidState::Liquid)
                .map_err(|error| fs_material::MaterialError::State {
                    what: error.to_string(),
                }),
        }
    }
}

impl SolidSource {
    pub(super) fn temperature_dependent(&self) -> bool {
        self.temperature_axis.is_some()
    }

    pub(super) fn has_storage_range(&self) -> bool {
        self.storage_temperature_range_k.is_some()
    }

    fn law(&self) -> Result<SolidConductivityLaw<'_>> {
        let axis = self.temperature_axis.as_deref()
            .ok_or_else(|| bad("temperature-dependent CHT needs every solid source in the same mode"))?;
        SolidConductivityLaw::try_new(
            self.pack.card(), self.state.clone(), axis, [THERMAL_CONDUCTIVITY_PROPERTY; 3],
        ).map_err(|error| bad(error.to_string()))
    }

    fn thermal_law(&self) -> Result<fs_lbm::material::conjugate::solid::SolidThermalLaw<'_>> {
        let span = self.storage_temperature_range_k.ok_or_else(|| {
            bad("temperature-dependent unsteady solids require storage_temperature_range_k")
        })?;
        self.law()?.with_heat_storage(
            DENSITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_PROPERTY, span,
        ).map_err(|error| bad(error.to_string()))
    }

    pub(super) fn json(&self) -> String {
        let value = |name| {
            self.state
                .property(name)
                .expect("required source property resolved")
                .value_si()
        };
        let storage = self.storage_temperature_range_k.map_or_else(String::new, |[lo, hi]| {
            format!(",\"storage_temperature_range_k\":[{lo},{hi}]")
        });
        format!(
            "{{\"index\":{},{},\"pack_identity\":{},\"material_state\":{},\"conductivity_w_m_k\":{},\"density_kg_m3\":{},\"specific_heat_j_kg_k\":{}{storage}}}",
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

pub(super) const EVOLVING_NO_CLAIM: &str = "steady laminar fixed-geometry Boussinesq flow and conservative conjugate heat transfer; original gas or liquid density, viscosity and conductivity and solid conductivity reevaluated at local temperature with source claims pinned; reference fluid density retained in inertia and continuity; heat advection uses reference Cp unless fluid_enthalpy selects original source h(T) with conservative first-order upwind; source density and Cp departures remain at most 5 percent; fixed thermodynamic pressure and composition; isotropic source solids; no compressibility, solid storage, melting, deformation, radiation, turbulence closure or phase-stability claim; staircase geometry and one declared resolution unless grid_convergence is requested; Estimated numerical evidence, not experimental validation";

/// Independent units and limits for the nonlinear solid heat-storage solve.
pub(super) fn storage_controls(
    value: &J,
) -> Result<fs_lbm::material::conjugate::buoyant::SolidThermalConfig> {
    use fs_lbm::material::conjugate::buoyant::SolidThermalConfig;
    fields(value, &[
        "maximum_iterations", "temperature_tolerance_k", "conductivity_relative_tolerance",
        "storage_tolerance_j", "temperature_relaxation",
    ])?;
    let defaults = SolidThermalConfig::default();
    let at = "transient.solid_storage";
    let maximum_iterations = super::optional_number(value, "maximum_iterations", at)?
        .unwrap_or(defaults.maximum_iterations as f64);
    if !maximum_iterations.is_finite()
        || !(1.0..=100_000.0).contains(&maximum_iterations)
        || maximum_iterations.fract() != 0.0
    {
        return Err(bad("transient.solid_storage.maximum_iterations must be an integer in 1..=100000"));
    }
    let positive = |name, default| {
        let number = super::optional_number(value, name, at)?.unwrap_or(default);
        if !number.is_finite() || number <= 0.0 {
            return Err(bad(format!("{at}.{name} must be finite and positive")));
        }
        Ok(number)
    };
    let controls = SolidThermalConfig {
        maximum_iterations: maximum_iterations as usize,
        temperature_tolerance_k: positive("temperature_tolerance_k", defaults.temperature_tolerance_k)?,
        conductivity_relative_tolerance: positive(
            "conductivity_relative_tolerance", defaults.conductivity_relative_tolerance,
        )?,
        storage_tolerance_j: positive("storage_tolerance_j", defaults.storage_tolerance_j)?,
        temperature_relaxation: positive("temperature_relaxation", defaults.temperature_relaxation)?,
    };
    if controls.temperature_relaxation > 1.0 {
        return Err(bad("transient.solid_storage.temperature_relaxation must lie in (0,1]"));
    }
    Ok(controls)
}

/// Existing CLI geometry and boundaries, existing coupled momentum/heat owner.
/// Every runtime query borrows the same owned cards admitted by the loader.
pub(super) fn solve(
    scene: &super::Scene,
    domain: &fs_lbm::conjugate::VoxelDomain,
    flow: &fs_lbm::conjugate::SimpleConfig,
    thermal: &fs_lbm::conjugate::ThermalSetup,
    gate: &fs_exec::CancelGate,
) -> Result<(
    fs_lbm::conjugate::FvFlow,
    fs_lbm::conjugate::EnergySolution,
    usize,
    String,
)> {
    use fs_lbm::material::conjugate::buoyant::{MaterialBuoyancyConfig, MaterialBuoyancyError};
    use fs_lbm::material::transport::boussinesq::BoussinesqConfig;
    let source = scene
        .fluid_source
        .as_ref()
        .ok_or_else(|| bad("temperature-dependent CHT requires a fluid source"))?;
    let fluid = source.law()?;
    let solids = scene.solid_sources.iter().map(SolidSource::law).collect::<Result<Vec<_>>>()?;
    let reference_temperature = source
        .evolving_reference_temperature()
        .expect("evolving fluid admitted");
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
    let bound = MaterialChtFluid::try_new(fluid.reference())
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
    let failure = |error: MaterialBuoyancyError| match &error {
        MaterialBuoyancyError::Solver(source) => super::solver_failure(source),
        MaterialBuoyancyError::NotConverged(_)
            | MaterialBuoyancyError::EnthalpyNotConverged { .. } => super::Failure {
            code: "cooling-cht-not-steady", message: error.to_string(),
        },
        _ => bad(error.to_string()),
    };
    let (run, caloric_report) = if let Some(input) = scene.fluid_enthalpy {
        let caloric = input.bind(source)?;
        let energy = fs_lbm::conjugate::EnergyConfig {
            scheme: fs_lbm::conjugate::ConvectionScheme::Upwind,
            ..fs_lbm::conjugate::EnergyConfig::default()
        };
        let run = bound.solve_natural_convection_with_material_solids_and_enthalpy(
            domain, &solids, flow, thermal, &energy, &initial, controls,
            &caloric, input.controls(reference_temperature), gate,
            |_, temperature, _| fluid.at_temperature(temperature),
        ).map_err(failure)?;
        let report = caloric.report(reference_temperature, &run.enthalpy)?;
        (run.convection, report)
    } else {
        let run = bound.solve_natural_convection_with_material_solids(
            domain, &solids, flow, thermal, &fs_lbm::conjugate::EnergyConfig::default(),
            &initial, controls, gate, |_, temperature, _| fluid.at_temperature(temperature),
        ).map_err(failure)?;
        (run, "null".into())
    };
    let fluid = fluid_endpoint_json(&run.endpoint, gate)?;
    let solid = endpoint_json(run.solid_endpoint.iter().map(|sample| (
        sample.cell, sample.temperature_k, sample.conductivity_w_m_k[0], &sample.material,
    )), gate)?;
    let report = run.report;
    let evidence = format!(concat!(
        "{{\"mode\":\"temperature-dependent\",\"passes\":{},",
        "\"temperature_defect_k\":{},\"transport_relative_defect\":{},\"force_defect_n_m3\":{},",
        "\"temperature_tolerance_k\":{},\"transport_relative_tolerance\":{},\"force_tolerance_n_m3\":{},",
        "\"max_relative_density_departure\":0.05,\"max_relative_specific_heat_departure\":0.05,\"solid_to_fluid_heat_w\":{},",
        "\"fluid_endpoint\":{},\"solid_endpoint\":{},\"fluid_enthalpy\":{}}}"),
        report.passes, super::num(report.temperature_defect_k)?,
        super::num(report.transport_relative_defect)?, super::num(report.force_defect_n_m3)?,
        super::num(controls.temperature_tolerance_k)?, super::num(controls.transport_relative_tolerance)?,
        super::num(controls.force_tolerance_n_m3)?,
        super::num(run.energy.solid_to_fluid_heat_w(domain))?, fluid, solid, caloric_report);
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

/// Report the actual gas EOS and transport samples, not a frozen reference.
/// Source identities remain on `material_inputs`; this digest covers the
/// numerical endpoint fields in original cell order and their source receipts.
fn gas_endpoint_json<'a>(
    samples: impl Iterator<Item = (usize, &'a GasRecord)>,
    gate: &fs_exec::CancelGate,
) -> Result<String> {
    let mut digest = fs_blake3::DomainHasher::new("frankensim.cooling-cht.gas-endpoint.v1");
    let mut count = 0;
    let mut minima = [f64::INFINITY; 5];
    let mut maxima = [f64::NEG_INFINITY; 5];
    for (cell, record) in samples {
        if gate.is_requested() {
            return Err(super::solver_failure(
                &fs_lbm::conjugate::ChtError::Cancelled,
            ));
        }
        count += 1;
        digest.update(&(cell as u64).to_le_bytes());
        let gas = record.state();
        for (i, value) in [
            gas.temperature,
            gas.density,
            gas.specific_heat_cp,
            gas.dynamic_viscosity,
            gas.thermal_conductivity,
        ]
        .into_iter()
        .enumerate()
        {
            minima[i] = minima[i].min(value);
            maxima[i] = maxima[i].max(value);
            digest.update(&value.to_le_bytes());
        }
        for value in [
            gas.pressure,
            gas.specific_gas_constant,
            gas.water_mole_fraction,
        ] {
            digest.update(&value.to_le_bytes());
        }
        if let GasRecord::Thermochemical(state) = record {
            for component in state.components() {
                digest.update(component.identity().as_bytes());
            }
            for source in state.source_packs() {
                digest.update(source.species_pack.as_bytes());
                digest.update(source.nasa9_pack.as_bytes());
                digest.update(source.transport_pack.as_bytes());
            }
        }
    }
    if count == 0 {
        return Ok("{\"cells\":0}".into());
    }
    let mut fields = format!(
        "{{\"kind\":\"gas\",\"cells\":{count},\"samples_identity\":{}",
        quote(&digest.finalize().to_hex())
    );
    for (i, name) in [
        "temperature_k",
        "density_kg_m3",
        "specific_heat_j_kg_k",
        "dynamic_viscosity_pa_s",
        "conductivity_w_m_k",
    ]
    .into_iter()
    .enumerate()
    {
        fields.push_str(&format!(
            ",\"minimum_{name}\":{},\"maximum_{name}\":{}",
            super::num(minima[i])?,
            super::num(maxima[i])?
        ));
    }
    fields.push('}');
    Ok(fields)
}


/// Advance the original unsteady momentum and storage owners with actual
/// fluid coefficients, retaining the same source law as the steady adapter.
#[allow(clippy::too_many_arguments)]
pub(super) fn march(
    scene: &super::Scene,
    domain: &fs_lbm::conjugate::VoxelDomain,
    flow: &fs_lbm::conjugate::SimpleConfig,
    thermal: &fs_lbm::conjugate::ThermalSetup,
    initial_temperature: &[f64],
    unsteady: &fs_lbm::conjugate::UnsteadyConfig,
    inlet_schedule: impl Fn(f64) -> f64,
    power_schedule: impl Fn(f64) -> f64,
    radiation: Option<&fs_lbm::conjugate::RadiationConfig>,
    gate: &fs_exec::CancelGate,
) -> Result<(fs_lbm::conjugate::ConjugateMarch, String)> {
    use fs_lbm::material::conjugate::buoyant::MaterialBuoyancyError;
    use fs_lbm::material::transport::boussinesq::BoussinesqConfig;

    let source = scene.fluid_source.as_ref()
        .ok_or_else(|| bad("unsteady material CHT requires a fluid source"))?;
    let fluid = source.law()?;
    let bound = MaterialChtFluid::try_new(fluid.reference())
        .map_err(|error| bad(error.to_string()))?;
    let bounds = BoussinesqConfig {
        gravity_m_s2: scene.gravity.unwrap_or([0.0; 3]),
        max_relative_density_departure: 0.05,
        max_relative_specific_heat_departure: 0.05,
    };
    let failure = |error: MaterialBuoyancyError| match &error {
        MaterialBuoyancyError::Solver(fs_lbm::conjugate::ChtError::InvalidInput {
            field: "unsteady.material_energy_solves", ..
        }) => super::Failure { code: "cooling-cht-budget", message: error.to_string() },
        MaterialBuoyancyError::Solver(source) => super::solver_failure(source),
        MaterialBuoyancyError::EnthalpyStepNotConverged { .. } => super::Failure {
            code: "cooling-cht-not-steady", message: error.to_string(),
        },
        _ => bad(error.to_string()),
    };
    let (run, caloric_report, solid_feedback) = if let Some(input) = scene.fluid_enthalpy {
        let caloric = input.bind(source)?;
        let reference_temperature = source.evolving_reference_temperature()
            .expect("evolving fluid admitted");
        let energy = fs_lbm::conjugate::EnergyConfig {
            scheme: fs_lbm::conjugate::ConvectionScheme::Upwind,
            ..fs_lbm::conjugate::EnergyConfig::default()
        };
        if scene.solid_sources.iter().any(SolidSource::temperature_dependent) {
            let solids = scene.solid_sources.iter().map(SolidSource::thermal_law)
                .collect::<Result<Vec<_>>>()?;
            let solid_controls = scene.transient.as_ref().and_then(|t| t.solid_storage)
                .unwrap_or_default();
            let result = bound.march_natural_convection_with_enthalpy_and_material_solids(
                domain, &solids, flow, thermal, &energy, initial_temperature,
                bounds, &caloric, input.transient_controls(reference_temperature)?,
                solid_controls, unsteady, inlet_schedule, power_schedule, radiation, gate,
                |_, temperature, _| fluid.at_temperature(temperature),
            ).map_err(failure)?;
            let report = caloric.transient_report(
                reference_temperature, &result.fluid.enthalpy, &result.fluid.records,
            )?;
            let endpoint = endpoint_json(result.solid_endpoint.iter().map(|sample| (
                sample.cell, sample.temperature_k, sample.conductivity_w_m_k[0], &sample.material,
            )), gate)?;
            let mut storage = format!(concat!(
                ",\"solid_coupling\":\"joint-implicit-endpoint\",\"solid_endpoint\":{},",
                "\"solid_storage\":{{\"model\":\"fixed-reference-mass-integral-cp\",",
                "\"maximum_iterations\":{},\"temperature_tolerance_k\":{},",
                "\"conductivity_relative_tolerance\":{},\"storage_tolerance_j\":{},",
                "\"temperature_relaxation\":{},\"sources\":["),
                endpoint, solid_controls.maximum_iterations,
                super::num(solid_controls.temperature_tolerance_k)?,
                super::num(solid_controls.conductivity_relative_tolerance)?,
                super::num(solid_controls.storage_tolerance_j)?,
                super::num(solid_controls.temperature_relaxation)?);
            for (index, solid) in solids.iter().enumerate() {
                if gate.is_requested() {
                    return Err(super::solver_failure(&fs_lbm::conjugate::ChtError::Cancelled));
                }
                if index > 0 { storage.push(','); }
                let [lo, hi] = solid.storage().temperature_span_k();
                storage.push_str(&format!(concat!(
                    "{{\"index\":{},\"storage_identity\":{},\"reference_density_kg_m3\":{},",
                    "\"temperature_range_k\":[{},{}]}}"),
                    index, quote(&solid.storage().identity().to_hex()),
                    super::num(solid.reference_density_kg_m3())?,
                    super::num(lo)?, super::num(hi)?));
            }
            storage.push_str("],\"records\":[");
            for (index, record) in result.solid_records.iter().enumerate() {
                if gate.is_requested() {
                    return Err(super::solver_failure(&fs_lbm::conjugate::ChtError::Cancelled));
                }
                if index > 0 { storage.push(','); }
                storage.push_str(&format!(concat!(
                    "{{\"time_s\":{},\"iterations\":{},\"temperature_defect_k\":{},",
                    "\"conductivity_relative_defect\":{},\"storage_defect_j\":{},",
                    "\"solid_stored_energy_change_j\":{}}}"),
                    super::num(record.time_s)?, record.iterations,
                    super::num(record.temperature_defect_k)?,
                    super::num(record.conductivity_relative_defect)?,
                    super::num(record.storage_defect_j)?,
                    super::num(record.solid_stored_energy_change_j)?));
            }
            storage.push_str("]}");
            (result.fluid.convection, report, storage)
        } else {
            let run = bound.march_natural_convection_with_enthalpy(
                domain, &scene.materials, flow, thermal, &energy, initial_temperature,
                bounds, &caloric, input.transient_controls(reference_temperature)?,
                unsteady, inlet_schedule, power_schedule, radiation, gate,
                |_, temperature, _| fluid.at_temperature(temperature),
            ).map_err(failure)?;
            let report = caloric.transient_report(
                reference_temperature, &run.enthalpy, &run.records,
            )?;
            (run.convection, report, String::new())
        }
    } else if scene.solid_sources.iter().any(SolidSource::temperature_dependent) {
        let solids = scene.solid_sources.iter().map(SolidSource::thermal_law)
            .collect::<Result<Vec<_>>>()?;
        let controls = scene.transient.as_ref().and_then(|t| t.solid_storage)
            .unwrap_or_default();
        let result = bound.march_natural_convection_with_material_solids(
            domain, &solids, flow, thermal,
            &fs_lbm::conjugate::EnergyConfig::default(), initial_temperature,
            bounds, unsteady, controls, inlet_schedule, power_schedule, radiation, gate,
            |_, temperature, _| fluid.at_temperature(temperature),
        ).map_err(failure)?;
        let endpoint = endpoint_json(result.solid_endpoint.iter().map(|sample| (
            sample.cell, sample.temperature_k, sample.conductivity_w_m_k[0], &sample.material,
        )), gate)?;
        let mut storage = format!(concat!(
            ",\"solid_coupling\":\"implicit-endpoint\",\"solid_endpoint\":{},",
            "\"solid_storage\":{{\"model\":\"fixed-reference-mass-integral-cp\",",
            "\"maximum_iterations\":{},\"temperature_tolerance_k\":{},",
            "\"conductivity_relative_tolerance\":{},\"storage_tolerance_j\":{},",
            "\"temperature_relaxation\":{},\"sources\":["),
            endpoint, controls.maximum_iterations, super::num(controls.temperature_tolerance_k)?,
            super::num(controls.conductivity_relative_tolerance)?,
            super::num(controls.storage_tolerance_j)?, super::num(controls.temperature_relaxation)?);
        for (index, solid) in solids.iter().enumerate() {
            if gate.is_requested() {
                return Err(super::solver_failure(&fs_lbm::conjugate::ChtError::Cancelled));
            }
            if index > 0 { storage.push(','); }
            let [lo, hi] = solid.storage().temperature_span_k();
            storage.push_str(&format!(concat!(
                "{{\"index\":{},\"storage_identity\":{},\"reference_density_kg_m3\":{},",
                "\"temperature_range_k\":[{},{}]}}"),
                index, quote(&solid.storage().identity().to_hex()),
                super::num(solid.reference_density_kg_m3())?, super::num(lo)?, super::num(hi)?));
        }
        storage.push_str("],\"records\":[");
        for (index, record) in result.solid_records.iter().enumerate() {
            if gate.is_requested() {
                return Err(super::solver_failure(&fs_lbm::conjugate::ChtError::Cancelled));
            }
            if index > 0 { storage.push(','); }
            storage.push_str(&format!(concat!(
                "{{\"time_s\":{},\"iterations\":{},\"temperature_defect_k\":{},",
                "\"conductivity_relative_defect\":{},\"storage_defect_j\":{},",
                "\"solid_stored_energy_change_j\":{}}}"),
                super::num(record.time_s)?, record.iterations,
                super::num(record.temperature_defect_k)?,
                super::num(record.conductivity_relative_defect)?,
                super::num(record.storage_defect_j)?,
                super::num(record.solid_stored_energy_change_j)?));
        }
        storage.push_str("]}");
        (result.fluid, "null".into(), storage)
    } else {
        let run = bound.march_natural_convection(
            domain, &scene.materials, flow, thermal,
            &fs_lbm::conjugate::EnergyConfig::default(), initial_temperature,
            bounds, unsteady, inlet_schedule, power_schedule, radiation, gate,
            |_, temperature, _| fluid.at_temperature(temperature),
        ).map_err(failure)?;
        (run, "null".into(), String::new())
    };
    let applied = fluid_endpoint_json(&run.applied, gate)?;
    let endpoint = fluid_endpoint_json(&run.endpoint, gate)?;
    let feedback = format!(concat!(
        "{{\"mode\":\"temperature-dependent\",\"coupling\":\"explicit-first-order\",",
        "\"max_relative_density_departure\":0.05,\"max_relative_specific_heat_departure\":0.05,",
        "\"applied_time_s\":{},\"endpoint_time_s\":{},",
        "\"fluid_applied\":{},\"fluid_endpoint\":{},\"fluid_enthalpy\":{}{}}}"),
        super::num(run.applied_time_s)?, super::num(run.endpoint_time_s)?,
        applied, endpoint, caloric_report, solid_feedback);
    Ok((run.march, feedback))
}

/// Source-caloric storage shares the same fixed-density physical approximation.
pub(super) const UNSTEADY_NO_CLAIM: &str = "unsteady laminar fixed-geometry Boussinesq flow and conservative conjugate heat storage; original fluid density, viscosity and conductivity reevaluated at each accepted temperature with fixed pressure, composition and selected source claims; reference fluid density retained in inertia and continuity; storage and advection use reference Cp unless fluid_enthalpy selects original source h(T), backward-Euler source storage and conservative first-order upwind; source density and Cp departures remain at most 5 percent; explicit first-order fluid-material coupling even with BDF2 momentum; solids use authored/frozen storage or implicit endpoint conductivity and source-integrated heat capacity on a declared single-phase span at fixed reference mass; source-caloric fluid and source-integral solid storage may close in one joint implicit thermal endpoint; original optional surface-radiation model; no compressibility, latent phase change, deformation, turbulence closure or time/grid convergence claim; Estimated numerical evidence, not experimental validation";

fn fluid_endpoint_json(
    samples: &[fs_lbm::material::conjugate::buoyant::MaterialBuoyancySample<FluidState>],
    gate: &fs_exec::CancelGate,
) -> Result<String> {
    match samples.first().map(|sample| sample.material.material()) {
        Some(FluidState::Liquid(_)) => endpoint_json(
            samples.iter().map(|sample| {
                let FluidState::Liquid(state) = sample.material.material() else {
                    unreachable!("one fluid source owns every cell")
                };
                (
                    sample.cell, sample.temperature_k,
                    sample.material.resolved_transport().thermal_conductivity_w_m_k(),
                    state.material(),
                )
            }),
            gate,
        ),
        Some(FluidState::Gas(_)) => gas_endpoint_json(
            samples.iter().map(|sample| {
                let FluidState::Gas(state) = sample.material.material() else {
                    unreachable!("one fluid source owns every cell")
                };
                (sample.cell, state)
            }),
            gate,
        ),
        None => Ok("{\"cells\":0}".into()),
    }
}
