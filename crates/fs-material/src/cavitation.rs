//! Shared liquid/vapor properties for cavitating flow and bubble dynamics.
//!
//! Saturation pressure is the equilibrium vapor PARTIAL pressure at the liquid
//! temperature, never a gauge pressure or the pressure of trapped gas. A finite
//! nucleus still has capillary and noncondensable-gas pressures; `p < p_sat`
//! alone is not a complete inception criterion. These immutable coefficients
//! are frozen during a consumer's isothermal advance. Thermal feedback must
//! refresh the original source before the next accepted advance.
//!
//! The property view below owns no flow, bubble history, phase-transfer rate or thermal
//! integrator. Saturated vapor density and latent enthalpy are explicitly
//! optional; their absence must never be replaced by an ideal-gas assumption.

/// Finite dissolved-gas exchange with an original gas inventory and calorics.
pub mod dissolved_gas;

use crate::MaterialError;
use crate::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use crate::liquid::surface::{SURFACE_TENSION_DIMS, require_frozen_coefficient_support};
use crate::liquid::{DYNAMIC_VISCOSITY_DIMS, KINEMATIC_VISCOSITY_DIMS};
use crate::liquid_temperature::LiquidViscosityBasis;
use crate::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_matdb::{MaterialCard, QueryPoint};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature, Velocity};
use fs_thermochem::water::if97::{If97SaturationPair, regions::If97FluidState};
use fs_thermochem::water::transport::If97TransportState;

/// Explicit saturated vapor data, at the same temperature as the liquid.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SaturatedVaporState {
    density_kg_m3: f64,
    latent_enthalpy_j_kg: f64,
}

impl SaturatedVaporState {
    /// Admit vapor density and positive vapor-minus-liquid enthalpy [J/kg].
    /// The caller owns phase equilibrium and the shared caloric reference.
    pub fn try_new(density_kg_m3: f64, latent_enthalpy_j_kg: f64) -> Result<Self, MaterialError> {
        if !positive(density_kg_m3) || !positive(latent_enthalpy_j_kg) {
            return Err(invalid(
                "saturated vapor density and latent enthalpy must be positive finite",
            ));
        }
        Ok(Self {
            density_kg_m3,
            latent_enthalpy_j_kg,
        })
    }
    /// Saturated vapor density [kg/m³].
    pub const fn density_kg_m3(self) -> f64 {
        self.density_kg_m3
    }
    /// Vapor-minus-liquid specific enthalpy [J/kg].
    pub const fn latent_enthalpy_j_kg(self) -> f64 {
        self.latent_enthalpy_j_kg
    }
}

/// Checked numerical liquid coefficients, without inferred source authority.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CavitationState {
    temperature_k: f64,
    liquid_density_kg_m3: f64,
    dynamic_viscosity_pa_s: f64,
    surface_tension_n_m: f64,
    saturation_pressure_pa: f64,
    liquid_sound_speed_m_s: f64,
    saturated_vapor: Option<SaturatedVaporState>,
}

impl CavitationState {
    /// Admit SI coefficients at one explicit temperature. Zero viscosity,
    /// surface tension or vapor pressure admits a caller-declared ideal limit;
    /// missing source data is never interpreted as such a limit.
    pub fn try_new(
        temperature_k: f64,
        liquid_density_kg_m3: f64,
        dynamic_viscosity_pa_s: f64,
        surface_tension_n_m: f64,
        saturation_pressure_pa: f64,
        liquid_sound_speed_m_s: f64,
    ) -> Result<Self, MaterialError> {
        if ![temperature_k, liquid_density_kg_m3, liquid_sound_speed_m_s]
            .into_iter()
            .all(positive)
            || ![
                dynamic_viscosity_pa_s,
                surface_tension_n_m,
                saturation_pressure_pa,
            ]
            .into_iter()
            .all(|v| v.is_finite() && v >= 0.0)
        {
            return Err(invalid(
                "cavitation requires positive finite T, liquid density and sound speed, and nonnegative finite viscosity, surface tension and saturation pressure",
            ));
        }
        Ok(Self {
            temperature_k,
            liquid_density_kg_m3,
            dynamic_viscosity_pa_s,
            surface_tension_n_m,
            saturation_pressure_pa,
            liquid_sound_speed_m_s,
            saturated_vapor: None,
        })
    }

    /// Attach explicit equilibrium vapor data. Critical collapse and a vapor
    /// denser than its liquid refuse; no gas composition or EOS is inferred.
    pub fn with_saturated_vapor(
        mut self,
        vapor: SaturatedVaporState,
    ) -> Result<Self, MaterialError> {
        if vapor.density_kg_m3 >= self.liquid_density_kg_m3 {
            return Err(invalid(
                "distinct saturated vapor must be less dense than its liquid",
            ));
        }
        let volume_gap = 1.0 / vapor.density_kg_m3 - 1.0 / self.liquid_density_kg_m3;
        let pressure_work = self.saturation_pressure_pa * volume_gap;
        if !positive(volume_gap)
            || !pressure_work.is_finite()
            || !(vapor.latent_enthalpy_j_kg - pressure_work).is_finite()
        {
            return Err(invalid(
                "liquid/vapor volume gap and pressure work must be representable",
            ));
        }
        self.saturated_vapor = Some(vapor);
        Ok(self)
    }

    /// Common equilibrium/interface temperature [K].
    pub const fn temperature_k(self) -> f64 {
        self.temperature_k
    }
    /// Surrounding liquid density [kg/m³].
    pub const fn liquid_density_kg_m3(self) -> f64 {
        self.liquid_density_kg_m3
    }
    /// Liquid dynamic shear viscosity [Pa s].
    pub const fn dynamic_viscosity_pa_s(self) -> f64 {
        self.dynamic_viscosity_pa_s
    }
    /// Liquid/bubble surface tension [N/m] in the declared interface environment.
    pub const fn surface_tension_n_m(self) -> f64 {
        self.surface_tension_n_m
    }
    /// Equilibrium vapor partial pressure [Pa], on an absolute reference.
    pub const fn saturation_pressure_pa(self) -> f64 {
        self.saturation_pressure_pa
    }
    /// Surrounding liquid sound speed [m/s], before bubble-mixture effects.
    pub const fn liquid_sound_speed_m_s(self) -> f64 {
        self.liquid_sound_speed_m_s
    }
    /// Explicit vapor/caloric data, or `None` when unavailable.
    pub const fn saturated_vapor(self) -> Option<SaturatedVaporState> {
        self.saturated_vapor
    }
}

/// Frozen coefficients from one original material or explicitly authored state.
pub trait CavitationProperties {
    /// Return the original owner's checked frozen cavitation coefficients.
    fn cavitation_state(&self) -> Result<CavitationState, MaterialError>;
}
impl CavitationProperties for CavitationState {
    fn cavitation_state(&self) -> Result<CavitationState, MaterialError> {
        Ok(*self)
    }
}
impl FluidAcoustics for CavitationState {
    /// Acoustic response of the surrounding liquid, not of a bubbly mixture.
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        FluidAcousticState::try_new(self.liquid_density_kg_m3, self.liquid_sound_speed_m_s)
    }
}
impl<T: CavitationProperties + ?Sized> CavitationProperties for &T {
    fn cavitation_state(&self) -> Result<CavitationState, MaterialError> {
        T::cavitation_state(*self)
    }
}

/// Exact source roles, preserving the source's own property vocabulary.
#[derive(Clone, Debug, PartialEq)]
pub struct CavitationPropertyRequirements {
    // Density, viscosity, surface tension, saturation pressure, sound speed.
    liquid: [ScalarPropertyRequirement; 5],
    vapor: Option<[ScalarPropertyRequirement; 2]>,
    viscosity_basis: LiquidViscosityBasis,
}

impl CavitationPropertyRequirements {
    /// Assign distinct scalar keys to the five liquid roles listed above and,
    /// optionally, [saturated vapor density, latent vaporization enthalpy].
    /// Latent enthalpy has units J/kg, not J/mol or a heat capacity.
    pub fn try_new(
        liquid: [ScalarPropertyRequirement; 5],
        vapor: Option<[ScalarPropertyRequirement; 2]>,
        viscosity_basis: LiquidViscosityBasis,
    ) -> Result<Self, MaterialStatePointError> {
        let viscosity_dims = match viscosity_basis {
            LiquidViscosityBasis::Dynamic => DYNAMIC_VISCOSITY_DIMS,
            LiquidViscosityBasis::Kinematic => KINEMATIC_VISCOSITY_DIMS,
        };
        let dims = [
            Density::DIMS,
            viscosity_dims,
            SURFACE_TENSION_DIMS,
            Pressure::DIMS,
            Velocity::DIMS,
            Density::DIMS,
            Dims([2, 0, -2, 0, 0, 0]),
        ];
        let all: Vec<_> = liquid.iter().chain(vapor.iter().flatten()).collect();
        for (i, (required, dims)) in all.iter().zip(dims).enumerate() {
            if required.quantity() != QuantitySpec::dimensional(dims)
                || required.hardness_test().is_some()
                || required.elastic_component().is_some()
                || required.strain_component().is_some()
                || required.stress_component().is_some()
                || all[..i].iter().any(|prior| prior.name() == required.name())
            {
                return Err(MaterialStatePointError::InvalidRequirement {
                    property: required.name().into(),
                    reason: "cavitation requires distinct exact scalar keys with the declared liquid/vapor SI quantities",
                });
            }
        }
        Ok(Self {
            liquid,
            vapor,
            viscosity_basis,
        })
    }

    /// Canonical names; the vapor/latent pair is requested only when needed.
    pub fn canonical(
        viscosity_basis: LiquidViscosityBasis,
        include_vapor_calorics: bool,
    ) -> Result<Self, MaterialStatePointError> {
        let req = |name, dims| {
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
        };
        let (viscosity_name, viscosity_dims) = match viscosity_basis {
            LiquidViscosityBasis::Dynamic => ("dynamic_viscosity", DYNAMIC_VISCOSITY_DIMS),
            LiquidViscosityBasis::Kinematic => ("kinematic_viscosity", KINEMATIC_VISCOSITY_DIMS),
        };
        Self::try_new(
            [
                req("density", Density::DIMS)?,
                req(viscosity_name, viscosity_dims)?,
                req("surface_tension", SURFACE_TENSION_DIMS)?,
                req("saturation_vapor_pressure", Pressure::DIMS)?,
                req("sound_speed", Velocity::DIMS)?,
            ],
            if include_vapor_calorics {
                Some([
                    req("saturated_vapor_density", Density::DIMS)?,
                    req("latent_vaporization_enthalpy", Dims([2, 0, -2, 0, 0, 0]))?,
                ])
            } else {
                None
            },
            viscosity_basis,
        )
    }
    /// Exact density, viscosity, surface tension, saturation pressure and sound roles.
    pub const fn liquid(&self) -> &[ScalarPropertyRequirement; 5] {
        &self.liquid
    }
    /// Optional exact saturated vapor density and latent enthalpy roles.
    pub const fn vapor(&self) -> Option<&[ScalarPropertyRequirement; 2]> {
        self.vapor.as_ref()
    }
    /// Original source viscosity basis; no alternate observation is selected.
    pub const fn viscosity_basis(&self) -> LiquidViscosityBasis {
        self.viscosity_basis
    }
}

/// Frozen cavitation coefficients retaining every selected source and query.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedCavitationState {
    material: ResolvedMaterialStatePoint,
    requirements: CavitationPropertyRequirements,
    temperature_axis: String,
    state: CavitationState,
}

impl ResolvedCavitationState {
    /// Adopt one already resolved bundle; no extra thermal properties are
    /// requested and no claim is replaced. Absolute temperature must be a
    /// semantic quantity, or the established legacy T-in-kelvin convention.
    pub fn try_from_material(
        material: ResolvedMaterialStatePoint,
        requirements: &CavitationPropertyRequirements,
        temperature_axis: &str,
    ) -> Result<Self, MaterialStatePointError> {
        let absolute_t = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        ));
        let quantity = material.axis_quantities().get(temperature_axis).copied();
        if quantity != Some(absolute_t)
            && !(temperature_axis == "T"
                && quantity.is_none_or(|q| q == QuantitySpec::dimensional(Temperature::DIMS)))
        {
            return Err(MaterialStatePointError::InvalidDerived {
                quantity: "cavitation temperature must be absolute kelvin",
            });
        }
        let temperature = material.state_coordinate(temperature_axis).ok_or(
            MaterialStatePointError::InvalidDerived {
                quantity: "cavitation requires its explicit temperature coordinate",
            },
        )?;
        let value = |required: &ScalarPropertyRequirement| {
            let p = material
                .property(required.name())
                .filter(|p| p.requirement() == required)
                .ok_or_else(|| MaterialStatePointError::InvalidRequirement {
                    property: required.name().into(),
                    reason: "cavitation bundle must retain the exact original requirement",
                })?;
            require_frozen_coefficient_support(p)?;
            if !p.has_only_fixed_external_conditions(Some(temperature_axis), false) {
                return Err(MaterialStatePointError::InvalidRequirement {
                    property: required.name().into(),
                    reason: "cavitation coefficients must share one temperature and fixed external conditions",
                });
            }
            Ok(p.value_si())
        };
        let [density, viscosity, sigma, pressure, sound] = &requirements.liquid;
        let rho = value(density)?;
        let source_viscosity = value(viscosity)?;
        let mu = match requirements.viscosity_basis {
            LiquidViscosityBasis::Dynamic => source_viscosity,
            LiquidViscosityBasis::Kinematic => rho * source_viscosity,
        };
        if source_viscosity > 0.0 && mu == 0.0 {
            return Err(MaterialStatePointError::InvalidDerived {
                quantity: "nonzero source viscosity cannot underflow to an inviscid cavitation limit",
            });
        }
        let map = |_| MaterialStatePointError::InvalidDerived {
            quantity: "coherent finite cavitation liquid/vapor coefficients",
        };
        let mut state = CavitationState::try_new(
            temperature,
            rho,
            mu,
            value(sigma)?,
            value(pressure)?,
            value(sound)?,
        )
        .map_err(map)?;
        if let Some([vapor_density, latent]) = &requirements.vapor {
            state = state
                .with_saturated_vapor(
                    SaturatedVaporState::try_new(value(vapor_density)?, value(latent)?)
                        .map_err(map)?,
                )
                .map_err(map)?;
        }
        Ok(Self {
            material,
            requirements: requirements.clone(),
            temperature_axis: temperature_axis.into(),
            state,
        })
    }
    /// Complete original material/query bundle and source usage receipts.
    pub const fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }
    /// Exact retained source roles, including optional vapor calorics.
    pub const fn requirements(&self) -> &CavitationPropertyRequirements {
        &self.requirements
    }
    /// Checked numerical snapshot derived from this source bundle.
    pub const fn state(&self) -> CavitationState {
        self.state
    }
    /// Atomically refresh all properties through the same selected source claims.
    pub fn at_query_point(
        &self,
        card: &MaterialCard,
        point: &QueryPoint,
    ) -> Result<Self, MaterialStatePointError> {
        Self::try_from_material(
            self.material.at_query_point(card, point)?,
            &self.requirements,
            &self.temperature_axis,
        )
    }
}
impl CavitationProperties for ResolvedCavitationState {
    fn cavitation_state(&self) -> Result<CavitationState, MaterialError> {
        Ok(self.state)
    }
}
impl FluidAcoustics for ResolvedCavitationState {
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        self.state.acoustic_state()
    }
}

/// Resolve only the cavitation consumer's explicitly selected source roles.
pub fn resolve_cavitation_state(
    card: &MaterialCard,
    point: &QueryPoint,
    requirements: &CavitationPropertyRequirements,
    temperature_axis: &str,
    selection: MaterialPropertySelection,
) -> Result<ResolvedCavitationState, MaterialStatePointError> {
    let properties: Vec<_> = requirements
        .liquid
        .iter()
        .chain(requirements.vapor.iter().flatten())
        .cloned()
        .collect();
    let material = resolve_material_state_point(card, point, &properties, selection)?;
    ResolvedCavitationState::try_from_material(material, requirements, temperature_axis)
}

/// Original IF97 saturation and industrial water transport, plus explicit sigma.
///
/// The incompressible surrounding liquid is evaluated on the saturated side at
/// the stated T, appropriate only when pressure-induced changes are negligible.
/// This is not a metastable EOS continuation or a critical-region bubble model.
/// Surface tension remains caller-supplied for the actual interface environment;
/// purity or surfactant effects are never inferred from the water source name.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SaturatedWaterCavitation {
    saturation: If97SaturationPair,
    liquid_transport: If97TransportState,
    state: CavitationState,
}
impl SaturatedWaterCavitation {
    /// Retain the original source pair and its unmodified enthalpy endpoints.
    pub fn try_new(
        saturation: If97SaturationPair,
        surface_tension_n_m: f64,
    ) -> Result<Self, MaterialError> {
        let liquid_transport =
            If97TransportState::from_state(If97FluidState::Gibbs(*saturation.liquid()))
                .map_err(|error| invalid(&format!("water cavitation transport: {error}")))?;
        let state = CavitationState::try_new(
            saturation.temperature().value(),
            saturation.liquid().density().value(),
            liquid_transport.dynamic_viscosity().value(),
            surface_tension_n_m,
            saturation.pressure().value(),
            saturation.liquid().sound_speed().value(),
        )?
        .with_saturated_vapor(SaturatedVaporState::try_new(
            saturation.vapor().density().value(),
            saturation.latent_enthalpy().value(),
        )?)?;
        Ok(Self {
            saturation,
            liquid_transport,
            state,
        })
    }
    /// Original saturated liquid/vapor thermodynamic pair.
    pub const fn saturation(&self) -> &If97SaturationPair {
        &self.saturation
    }
    /// Original industrial transport evaluated on the same saturated liquid.
    pub const fn liquid_transport(&self) -> &If97TransportState {
        &self.liquid_transport
    }
    /// Checked numerical coefficients, retaining the explicit interface tension.
    pub const fn state(&self) -> CavitationState {
        self.state
    }
}
impl CavitationProperties for SaturatedWaterCavitation {
    fn cavitation_state(&self) -> Result<CavitationState, MaterialError> {
        Ok(self.state)
    }
}
impl FluidAcoustics for SaturatedWaterCavitation {
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        self.state.acoustic_state()
    }
}

fn positive(v: f64) -> bool {
    v.is_finite() && v > 0.0
}
fn invalid(what: &str) -> MaterialError {
    MaterialError::Parameters { what: what.into() }
}
