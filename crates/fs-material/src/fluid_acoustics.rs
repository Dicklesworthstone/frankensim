//! One checked acoustic response for gases and liquids at their actual state.
//!
//! The original material or EOS owns density and adiabatic sound speed. This
//! view derives `K_s = rho c²`, `beta_s = 1/K_s` and `Z = rho c` without an
//! ideal-gas assumption. It describes local, small-signal, nondispersive bulk
//! response, not shear waves, attenuation or cavitation. Each original owner
//! explicitly selects its phase and relaxation closure; this view supplies none.
//! Source-card wrappers retain every query and usage receipt. The numerical
//! snapshot itself does not confer source or experimental authority.

/// Shared acoustic, caloric and transport state for viscous/thermal sound loss.
pub mod thermoviscous;

/// Explicit frozen two-phase compressibility shared by acoustic consumers.
pub mod mixture;
pub use mixture::FrozenHomogeneousMixture;

use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::{MaterialCard, QueryPoint};
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Velocity};

use crate::MaterialError;
use crate::fluid_transport::{
    FluidTransport, FluidTransportState, ThermalFluidTransport, ThermalFluidTransportState,
};
use crate::gas::{GasState, ResolvedGasState};
use crate::gas_mixture::ResolvedGasMixtureState;
use crate::liquid::{LiquidPropertyRequirements, ResolvedLiquidState};
use crate::liquid_temperature::LiquidViscosityBasis;
use crate::moist_air::ResolvedMoistAirState;
use crate::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
};
use crate::thermochemical_gas::ResolvedThermochemicalGasState;

/// Canonical adiabatic small-signal sound-speed property [m/s].
pub const SOUND_SPEED_PROPERTY: &str = "sound_speed";
/// Canonical adiabatic bulk-modulus property [Pa], distinct from isothermal K.
pub const ADIABATIC_BULK_MODULUS_PROPERTY: &str = "adiabatic_bulk_modulus";

/// Immutable, coherent numerical acoustic coefficients at one material state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FluidAcousticState {
    density: f64,
    sound_speed: f64,
    bulk_modulus: f64,
    compressibility: f64,
    impedance: f64,
}

impl FluidAcousticState {
    /// Admit explicitly supplied SI density and adiabatic sound speed.
    /// No material identity, source support or phase stability is inferred.
    pub fn try_new(density_kg_m3: f64, sound_speed_m_s: f64) -> Result<Self, MaterialError> {
        let impedance = density_kg_m3 * sound_speed_m_s;
        Self::admit(density_kg_m3, sound_speed_m_s, impedance * sound_speed_m_s)
    }

    /// Derive c from an explicitly adiabatic modulus, retaining the input K.
    /// Isothermal compressibility cannot be used here without caloric closure.
    pub fn from_adiabatic_bulk_modulus(
        density_kg_m3: f64,
        bulk_modulus_pa: f64,
    ) -> Result<Self, MaterialError> {
        Self::admit(
            density_kg_m3,
            (bulk_modulus_pa / density_kg_m3).sqrt(),
            bulk_modulus_pa,
        )
    }

    fn admit(density: f64, sound_speed: f64, bulk_modulus: f64) -> Result<Self, MaterialError> {
        let impedance = density * sound_speed;
        let compressibility = 1.0 / bulk_modulus;
        if [
            density,
            sound_speed,
            bulk_modulus,
            impedance,
            compressibility,
        ]
        .into_iter()
        .any(|v| !v.is_finite() || v <= 0.0)
        {
            return Err(MaterialError::Parameters {
                what: "fluid acoustic density, sound speed, adiabatic modulus, impedance and compressibility must be positive finite values".into(),
            });
        }
        Ok(Self {
            density,
            sound_speed,
            bulk_modulus,
            compressibility,
            impedance,
        })
    }

    /// Bulk density [kg/m³].
    pub const fn density_kg_m3(self) -> f64 {
        self.density
    }
    /// Adiabatic, small-signal sound speed [m/s].
    pub const fn sound_speed_m_s(self) -> f64 {
        self.sound_speed
    }
    /// Adiabatic bulk modulus `rho c²` [Pa].
    pub const fn adiabatic_bulk_modulus_pa(self) -> f64 {
        self.bulk_modulus
    }
    /// Adiabatic compressibility `1/K_s` [1/Pa].
    pub const fn adiabatic_compressibility_per_pa(self) -> f64 {
        self.compressibility
    }
    /// Plane traveling-wave specific acoustic impedance `rho c` [Pa s/m].
    pub const fn characteristic_impedance_pa_s_m(self) -> f64 {
        self.impedance
    }
}

/// Checked acoustic view of one original material/EOS evaluation.
///
/// Implementations must return density and sound speed from the same state.
/// Source identities, query coordinates and phase policy remain on that owner;
/// this trait never replaces them with a gas law or a material-name default.
pub trait FluidAcoustics {
    /// Return a coherent local small-signal acoustic snapshot.
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError>;
}

impl<T: FluidAcoustics + ?Sized> FluidAcoustics for &T {
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        T::acoustic_state(*self)
    }
}
impl FluidAcoustics for FluidAcousticState {
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        Ok(*self)
    }
}
impl FluidAcoustics for GasState {
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        // Public gas fields may be stale or independently mutated. Reuse the
        // original EOS/caloric/transport admission rather than accepting rho/c
        // as if they were an unrelated, explicitly authored numerical tuple.
        self.validate_coherence()?;
        FluidAcousticState::try_new(self.density, self.sound_speed)
    }
}

macro_rules! gas_owner {
    ($($ty:ty),+ $(,)?) => {$ (
        impl FluidAcoustics for $ty {
            fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
                self.state().acoustic_state()
            }
        }
    )+};
}
gas_owner!(
    ResolvedGasState,
    ResolvedGasMixtureState,
    ResolvedMoistAirState,
    ResolvedThermochemicalGasState
);

// The retained source().state() is an ideal reference. Acoustic consumers
// must use the checked real-gas root and its own caloric derivatives.
impl FluidAcoustics
    for crate::thermochemical_gas::nonideal::ResolvedThermochemicalRealGasState
{
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        FluidAcousticState::try_new(self.density_kg_m3(), self.sound_speed_m_s())
    }
}

macro_rules! water_owner {
    ($($ty:ty),+ $(,)?) => {$ (
        impl FluidAcoustics for $ty {
            fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
                // These immutable source states already admit phase/region and
                // derive c from the same Gibbs/Helmholtz derivatives as rho.
                FluidAcousticState::try_new(self.density().value(), self.sound_speed().value())
            }
        }
    )+};
}
water_owner!(
    fs_thermochem::water::if97::If97State,
    fs_thermochem::water::if97::region3::If97DenseState,
    fs_thermochem::water::if97::regions::If97FluidState,
);

/// Explicitly selected acoustic source parameter. The distinction between
/// adiabatic and isothermal response is a caller-declared physical role, not
/// something inferred from a property's name or its pressure dimensions.
#[derive(Clone, Debug, PartialEq)]
pub enum FluidAcousticParameter {
    /// Source adiabatic small-signal sound speed [m/s].
    SoundSpeed(ScalarPropertyRequirement),
    /// Source adiabatic bulk modulus [Pa].
    AdiabaticBulkModulus(ScalarPropertyRequirement),
}

impl FluidAcousticParameter {
    /// Select the canonical `sound_speed` dimension-only property.
    pub fn sound_speed() -> Result<Self, MaterialStatePointError> {
        Ok(Self::SoundSpeed(positive_requirement(
            SOUND_SPEED_PROPERTY,
            Velocity::DIMS,
        )?))
    }
    /// Select the canonical `adiabatic_bulk_modulus` dimension-only property.
    pub fn adiabatic_bulk_modulus() -> Result<Self, MaterialStatePointError> {
        Ok(Self::AdiabaticBulkModulus(positive_requirement(
            ADIABATIC_BULK_MODULUS_PROPERTY,
            Pressure::DIMS,
        )?))
    }
    /// Original exact requirement, including caller admissibility.
    pub const fn requirement(&self) -> &ScalarPropertyRequirement {
        match self {
            Self::SoundSpeed(r) | Self::AdiabaticBulkModulus(r) => r,
        }
    }
    fn tag(&self) -> u8 {
        match self {
            Self::SoundSpeed(_) => 0,
            Self::AdiabaticBulkModulus(_) => 1,
        }
    }
    fn dims(&self) -> Dims {
        match self {
            Self::SoundSpeed(_) => Velocity::DIMS,
            Self::AdiabaticBulkModulus(_) => Pressure::DIMS,
        }
    }
}

/// Acoustic response retaining the complete source bundle and parameter roles.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedFluidAcousticState {
    material: ResolvedMaterialStatePoint,
    density_requirement: ScalarPropertyRequirement,
    parameter: FluidAcousticParameter,
    state: FluidAcousticState,
    identity: ContentHash,
}

impl ResolvedFluidAcousticState {
    /// Adopt an existing complete source bundle without reselecting claims.
    /// Additional properties are retained and refreshed along with acoustics.
    pub fn try_from_material(
        material: ResolvedMaterialStatePoint,
        density: &ScalarPropertyRequirement,
        parameter: &FluidAcousticParameter,
    ) -> Result<Self, MaterialStatePointError> {
        require_scalar(density, Density::DIMS)?;
        require_scalar(parameter.requirement(), parameter.dims())?;
        let value = |required: &ScalarPropertyRequirement| {
            material.property(required.name())
                .filter(|p| p.requirement() == required)
                .map(|p| p.value_si())
                .ok_or_else(|| MaterialStatePointError::InvalidRequirement {
                    property: required.name().into(),
                    reason: "acoustic source bundle must contain the exact selected scalar requirement",
                })
        };
        let rho = value(density)?;
        let response = value(parameter.requirement())?;
        let state = match parameter {
            FluidAcousticParameter::SoundSpeed(_) => FluidAcousticState::try_new(rho, response),
            FluidAcousticParameter::AdiabaticBulkModulus(_) => {
                FluidAcousticState::from_adiabatic_bulk_modulus(rho, response)
            }
        }
        .map_err(|_| MaterialStatePointError::InvalidDerived {
            quantity: "finite positive fluid acoustic response",
        })?;
        let mut identity = DomainHasher::new("org.frankensim.fs-material.fluid-acoustics.v1");
        identity.update(material.identity().as_bytes());
        identity.update(&[parameter.tag()]);
        for r in [density, parameter.requirement()] {
            identity.update(&(r.name().len() as u64).to_le_bytes());
            identity.update(r.name().as_bytes());
        }
        for v in [
            state.density,
            state.sound_speed,
            state.bulk_modulus,
            state.compressibility,
            state.impedance,
        ] {
            identity.update(&v.to_bits().to_le_bytes());
        }
        Ok(Self {
            material,
            density_requirement: density.clone(),
            parameter: parameter.clone(),
            state,
            identity: identity.finalize(),
        })
    }
    /// Original material, query coordinates and all usage receipts.
    pub const fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }
    /// Explicit source parameter and physical role.
    pub const fn parameter(&self) -> &FluidAcousticParameter {
        &self.parameter
    }
    /// Selected density requirement, including its exact key and admissibility.
    pub const fn density_requirement(&self) -> &ScalarPropertyRequirement {
        &self.density_requirement
    }
    /// Checked numerical response; this borrows no alternate EOS or default.
    pub const fn state(&self) -> &FluidAcousticState {
        &self.state
    }
    /// Identity binds complete source/query receipts, roles and derived values.
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    /// Refresh every source claim at one new complete query. Unsupported states
    /// refuse atomically without changing this original snapshot or its sources.
    pub fn at_query_point(
        &self,
        card: &MaterialCard,
        point: &QueryPoint,
    ) -> Result<Self, MaterialStatePointError> {
        Self::try_from_material(
            self.material.at_query_point(card, point)?,
            &self.density_requirement,
            &self.parameter,
        )
    }
}
impl FluidAcoustics for ResolvedFluidAcousticState {
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        Ok(self.state)
    }
}

/// Resolve exactly density and the explicitly selected acoustic response.
/// No heat capacity or viscosity is required for a lossless acoustic consumer.
pub fn resolve_fluid_acoustic_state(
    card: &MaterialCard,
    point: &QueryPoint,
    density: &ScalarPropertyRequirement,
    parameter: &FluidAcousticParameter,
    selection: MaterialPropertySelection,
) -> Result<ResolvedFluidAcousticState, MaterialStatePointError> {
    require_scalar(density, Density::DIMS)?;
    require_scalar(parameter.requirement(), parameter.dims())?;
    let material = resolve_material_state_point(
        card,
        point,
        &[density.clone(), parameter.requirement().clone()],
        selection,
    )?;
    ResolvedFluidAcousticState::try_from_material(material, density, parameter)
}

/// Complete liquid transport and acoustic response from one common source query.
/// No cp/cv or ideal-gas identity is used to infer liquid compressibility.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedAcousticLiquidState {
    acoustic: ResolvedFluidAcousticState,
    liquid: ResolvedLiquidState,
    transport: FluidTransportState,
}
impl ResolvedAcousticLiquidState {
    fn admit(
        acoustic: ResolvedFluidAcousticState,
        requirements: &LiquidPropertyRequirements,
    ) -> Result<Self, MaterialStatePointError> {
        let liquid = ResolvedLiquidState::try_from_material_with_requirements(
            acoustic.material.clone(),
            requirements,
        )?;
        Self::from_views(acoustic, liquid)
    }

    /// Add acoustics to an existing complete liquid bundle without querying or
    /// replacing its source. The selected acoustic requirement must already be
    /// in that bundle. Retains any admitted density-expansion path, including
    /// its fixed pressure/composition context and analytic density tangent.
    ///
    /// All three views (sound, momentum/heat transport, and thermal expansion)
    /// refresh through the original liquid owner at one common query point.
    pub fn try_from_liquid(
        liquid: ResolvedLiquidState,
        parameter: &FluidAcousticParameter,
    ) -> Result<Self, MaterialStatePointError> {
        let acoustic = ResolvedFluidAcousticState::try_from_material(
            liquid.material().clone(),
            &liquid.requirements().properties()[0],
            parameter,
        )?;
        Self::from_views(acoustic, liquid)
    }

    fn from_views(
        acoustic: ResolvedFluidAcousticState,
        liquid: ResolvedLiquidState,
    ) -> Result<Self, MaterialStatePointError> {
        let transport =
            liquid
                .transport_state()
                .map_err(|_| MaterialStatePointError::InvalidDerived {
                    quantity: "finite positive acoustic liquid transport",
                })?;
        if transport.density_kg_m3().to_bits() != acoustic.state.density.to_bits() {
            return Err(MaterialStatePointError::InvalidDerived {
                quantity: "one density for liquid transport and acoustics",
            });
        }
        Ok(Self {
            acoustic,
            liquid,
            transport,
        })
    }
    /// Complete evidence-bearing source/query bundle.
    pub const fn material(&self) -> &ResolvedMaterialStatePoint {
        self.acoustic.material()
    }
    /// Source-bound acoustic state and its original parameter roles.
    pub const fn acoustic(&self) -> &ResolvedFluidAcousticState {
        &self.acoustic
    }
    /// Original liquid owner, including its expansion source and complete
    /// transport/acoustic bundle. No second set of claims is selected.
    pub const fn liquid(&self) -> &ResolvedLiquidState {
        &self.liquid
    }
    /// Basis actually reported by the source; no dynamic/kinematic fallback.
    pub const fn viscosity_basis(&self) -> LiquidViscosityBasis {
        self.liquid.requirements().viscosity_basis()
    }
    /// Exact original liquid property roles, without canonical-name substitution.
    pub const fn liquid_requirements(&self) -> &LiquidPropertyRequirements {
        self.liquid.requirements()
    }
    /// Recompute all transport and acoustic properties using the original claims.
    pub fn at_query_point(
        &self,
        card: &MaterialCard,
        point: &QueryPoint,
    ) -> Result<Self, MaterialStatePointError> {
        Self::try_from_liquid(
            self.liquid.at_query_point(card, point)?,
            self.acoustic.parameter(),
        )
    }
}
impl FluidAcoustics for ResolvedAcousticLiquidState {
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        self.acoustic.acoustic_state()
    }
}
impl FluidTransport for ResolvedAcousticLiquidState {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        Ok(self.transport)
    }
}
impl ThermalFluidTransport for ResolvedAcousticLiquidState {
    fn thermal_transport_state(&self) -> Result<ThermalFluidTransportState, MaterialError> {
        // Reuse the same temperature semantics and known/unknown/signed beta
        // as every other liquid consumer. In particular, no 1/T gas shortcut
        // and no silent zero-expansion assumption enters acoustic liquids.
        self.liquid.thermal_transport_state()
    }
}

/// Resolve liquid rho/cp/k/viscosity plus the selected adiabatic acoustic law
/// atomically at one query, for consumers needing both diffusion and sound.
pub fn resolve_acoustic_liquid_state(
    card: &MaterialCard,
    point: &QueryPoint,
    parameter: &FluidAcousticParameter,
    viscosity_basis: LiquidViscosityBasis,
    selection: MaterialPropertySelection,
) -> Result<ResolvedAcousticLiquidState, MaterialStatePointError> {
    resolve_acoustic_liquid_state_with_requirements(
        card,
        point,
        parameter,
        &LiquidPropertyRequirements::canonical(viscosity_basis)?,
        selection,
    )
}

/// Resolve original explicitly selected liquid and acoustic keys together.
/// This preserves sources whose vocabulary differs from the canonical names;
/// it does not rewrite their claims, receipts or physical quantity conventions.
pub fn resolve_acoustic_liquid_state_with_requirements(
    card: &MaterialCard,
    point: &QueryPoint,
    parameter: &FluidAcousticParameter,
    liquid_requirements: &LiquidPropertyRequirements,
    selection: MaterialPropertySelection,
) -> Result<ResolvedAcousticLiquidState, MaterialStatePointError> {
    require_scalar(parameter.requirement(), parameter.dims())?;
    let mut requirements = liquid_requirements.properties().to_vec();
    requirements.push(parameter.requirement().clone());
    let density = &liquid_requirements.properties()[0];
    let material = resolve_material_state_point(card, point, &requirements, selection)?;
    let acoustic = ResolvedFluidAcousticState::try_from_material(material, density, parameter)?;
    ResolvedAcousticLiquidState::admit(acoustic, liquid_requirements)
}

fn positive_requirement(
    name: &str,
    dims: Dims,
) -> Result<ScalarPropertyRequirement, MaterialStatePointError> {
    ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
}

fn require_scalar(
    requirement: &ScalarPropertyRequirement,
    dims: Dims,
) -> Result<(), MaterialStatePointError> {
    if requirement.quantity() != QuantitySpec::dimensional(dims)
        || requirement.hardness_test().is_some()
        || requirement.elastic_component().is_some()
        || requirement.strain_component().is_some()
        || requirement.stress_component().is_some()
    {
        return Err(MaterialStatePointError::InvalidRequirement {
            property: requirement.name().into(),
            reason: "fluid acoustic parameter requires an explicit dimension-only scalar of the correct units, not pressure, stress or a tensor component",
        });
    }
    Ok(())
}
