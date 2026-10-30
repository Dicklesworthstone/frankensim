//! Coherent small-signal fluid response for viscous and thermal boundary layers.
//!
//! A sound-loss consumer needs the same state's adiabatic sound speed, density,
//! viscosity, conductivity and Cp/Cv. Cp/Cv describes thermal coupling for a
//! liquid as well as a gas; it does not replace the original mechanical EOS.
//! No ideal-gas law, liquid Cv default, bulk-viscosity law or relaxation spectrum
//! is supplied here. Source-bearing owners retain their original phase and
//! query policy. A numerical snapshot alone carries no source authority.

use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::{MaterialCard, QueryPoint};

use super::{
    FluidAcousticParameter, FluidAcousticState, FluidAcoustics, ResolvedAcousticLiquidState,
    ResolvedFluidAcousticState, positive_requirement, require_scalar,
};
use crate::MaterialError;
use crate::fluid_transport::{
    FluidTransport, FluidTransportState, ThermalFluidTransport, ThermalFluidTransportState,
};
use crate::gas::{GasState, ResolvedGasState};
use crate::gas_mixture::ResolvedGasMixtureState;
use crate::liquid::{
    LiquidPropertyRequirements, VOLUMETRIC_EXPANSION_DIMS, VOLUMETRIC_EXPANSION_PROPERTY,
};
use crate::liquid_temperature::LiquidViscosityBasis;
use crate::moist_air::ResolvedMoistAirState;
use crate::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ResolvedMaterialStatePoint,
    SPECIFIC_HEAT_CAPACITY_DIMS, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use crate::thermochemical_gas::ResolvedThermochemicalGasState;
use crate::thermochemical_gas::nonideal::ResolvedThermochemicalRealGasState;

/// Canonical isochoric specific heat [J/(kg K)], distinct from the Cp role.
pub const SPECIFIC_HEAT_CAPACITY_CV_PROPERTY: &str = "specific_heat_capacity_cv";

/// One checked local Newtonian-fluid acoustic, caloric and transport tuple.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermoviscousFluidState {
    acoustic: FluidAcousticState,
    transport: FluidTransportState,
    specific_heat_cv: f64,
    heat_capacity_ratio: f64,
}

impl ThermoviscousFluidState {
    /// Bind independently declared Cv to acoustic and transport views of ONE
    /// state. Requires identical density, finite Cp >= Cv > 0, and finite Cp/Cv.
    /// Cv = Cp is allowed and gives exactly zero thermal acoustic coupling.
    pub fn try_new(
        acoustic: FluidAcousticState,
        transport: FluidTransportState,
        cv: f64,
    ) -> Result<Self, MaterialError> {
        Self::admit(
            acoustic,
            transport,
            cv,
            transport.specific_heat_j_kg_k() / cv,
        )
    }

    /// Preserve an original owner's admitted Cp/Cv, including its floating-point
    /// evaluation order. This is an explicitly supplied caloric ratio, never a
    /// relation used to infer density or sound speed from an ideal-gas law.
    pub fn from_heat_capacity_ratio(
        acoustic: FluidAcousticState,
        transport: FluidTransportState,
        ratio: f64,
    ) -> Result<Self, MaterialError> {
        Self::admit(
            acoustic,
            transport,
            transport.specific_heat_j_kg_k() / ratio,
            ratio,
        )
    }

    fn admit(
        acoustic: FluidAcousticState,
        transport: FluidTransportState,
        cv: f64,
        ratio: f64,
    ) -> Result<Self, MaterialError> {
        if acoustic.density_kg_m3().to_bits() != transport.density_kg_m3().to_bits()
            || !cv.is_finite()
            || cv <= 0.0
            || cv > transport.specific_heat_j_kg_k()
            || !ratio.is_finite()
            || ratio < 1.0
        {
            return Err(MaterialError::Parameters {
                what: "thermoviscous response requires one density and finite Cp >= Cv > 0 with representable Cp/Cv".into(),
            });
        }
        Ok(Self {
            acoustic,
            transport,
            specific_heat_cv: cv,
            heat_capacity_ratio: ratio,
        })
    }

    /// Original adiabatic acoustic coefficients, with no alternate EOS.
    pub const fn acoustic(&self) -> FluidAcousticState {
        self.acoustic
    }
    /// Original momentum/heat transport coefficients, including Cp.
    pub const fn transport(&self) -> FluidTransportState {
        self.transport
    }
    /// Isobaric-to-isochoric heat-capacity ratio, at least one.
    pub const fn heat_capacity_ratio(&self) -> f64 {
        self.heat_capacity_ratio
    }
    /// Isochoric specific heat [J/(kg K)].
    pub const fn specific_heat_cv_j_kg_k(&self) -> f64 {
        self.specific_heat_cv
    }
}

/// Sound and diffusion from one original fluid evaluation, with explicit Cv.
/// An acoustic/transport owner lacking Cv does not implicitly acquire this
/// capability: in particular, a lossless liquid is not silently assigned gamma.
pub trait ThermoviscousFluid: FluidAcoustics + FluidTransport {
    /// Checked shared state. Original source and phase errors remain refusals.
    fn thermoviscous_state(&self) -> Result<ThermoviscousFluidState, MaterialError>;
}

impl<T: ThermoviscousFluid + ?Sized> ThermoviscousFluid for &T {
    fn thermoviscous_state(&self) -> Result<ThermoviscousFluidState, MaterialError> {
        T::thermoviscous_state(*self)
    }
}
impl FluidAcoustics for ThermoviscousFluidState {
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        Ok(self.acoustic)
    }
}
impl FluidTransport for ThermoviscousFluidState {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        Ok(self.transport)
    }
}
impl ThermoviscousFluid for ThermoviscousFluidState {
    fn thermoviscous_state(&self) -> Result<ThermoviscousFluidState, MaterialError> {
        Ok(*self)
    }
}
impl ThermoviscousFluid for GasState {
    fn thermoviscous_state(&self) -> Result<ThermoviscousFluidState, MaterialError> {
        // Both views retain the gas owner's existing coherence admission. Stale
        // public gas fields must not be laundered into a raw numerical tuple.
        ThermoviscousFluidState::from_heat_capacity_ratio(
            self.acoustic_state()?,
            self.transport_state()?,
            self.gamma,
        )
    }
}
macro_rules! gas_owner {
    ($($ty:ty),+ $(,)?) => {$ (
        impl ThermoviscousFluid for $ty {
            fn thermoviscous_state(&self) -> Result<ThermoviscousFluidState, MaterialError> {
                self.state().thermoviscous_state()
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

impl ThermoviscousFluid for ResolvedThermochemicalRealGasState {
    fn thermoviscous_state(&self) -> Result<ThermoviscousFluidState, MaterialError> {
        // source().state() is an ideal reference, not this PR branch's response.
        // Keep its actual density, sound speed and caloric departure together.
        ThermoviscousFluidState::try_new(
            self.acoustic_state()?,
            self.transport_state()?,
            self.specific_heat_cv_j_kg_k(),
        )
    }
}

/// A source liquid's complete acoustic/transport bundle plus its selected Cv.
/// All properties share one query, original claims, receipts and refresh path.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedThermoviscousLiquidState {
    acoustic_liquid: ResolvedAcousticLiquidState,
    cv_requirement: ScalarPropertyRequirement,
    state: ThermoviscousFluidState,
    identity: ContentHash,
}

impl ResolvedThermoviscousLiquidState {
    /// Add Cv already present in the original bundle. This selects no new claim
    /// and retains the original liquid's optional density-expansion path.
    pub fn try_from_acoustic(
        acoustic_liquid: ResolvedAcousticLiquidState,
        cv_requirement: &ScalarPropertyRequirement,
    ) -> Result<Self, MaterialStatePointError> {
        require_cv(acoustic_liquid.liquid_requirements(), cv_requirement)?;
        let cv = acoustic_liquid.material().property(cv_requirement.name())
            .filter(|property| property.requirement() == cv_requirement)
            .ok_or_else(|| MaterialStatePointError::InvalidRequirement {
                property: cv_requirement.name().into(),
                reason: "thermoviscous liquid requires the exact selected Cv in its original bundle",
            })?.value_si();
        let state = ThermoviscousFluidState::try_new(
            *acoustic_liquid.acoustic().state(),
            acoustic_liquid.transport,
            cv,
        )
        .map_err(|_| MaterialStatePointError::InvalidDerived {
            quantity: "coherent thermoviscous liquid with finite Cp >= Cv > 0",
        })?;
        let mut identity = DomainHasher::new("org.frankensim.fs-material.thermoviscous-liquid.v1");
        identity.update(acoustic_liquid.acoustic().identity().as_bytes());
        for r in acoustic_liquid
            .liquid_requirements()
            .properties()
            .iter()
            .chain(std::iter::once(cv_requirement))
        {
            identity.update(&(r.name().len() as u64).to_le_bytes());
            identity.update(r.name().as_bytes());
        }
        identity.update(&[match acoustic_liquid.viscosity_basis() {
            LiquidViscosityBasis::Dynamic => 0,
            LiquidViscosityBasis::Kinematic => 1,
        }]);
        for v in [
            state.specific_heat_cv,
            state.heat_capacity_ratio,
            state.transport.specific_heat_j_kg_k(),
            state.transport.dynamic_viscosity_pa_s(),
            state.transport.thermal_conductivity_w_m_k(),
        ] {
            identity.update(&v.to_bits().to_le_bytes());
        }
        Ok(Self {
            acoustic_liquid,
            cv_requirement: cv_requirement.clone(),
            state,
            identity: identity.finalize(),
        })
    }
    /// Complete source bundle, including original query and every usage receipt.
    pub const fn material(&self) -> &ResolvedMaterialStatePoint {
        self.acoustic_liquid.material()
    }
    /// Original liquid and acoustic owners, including any density-expansion path.
    pub const fn acoustic_liquid(&self) -> &ResolvedAcousticLiquidState {
        &self.acoustic_liquid
    }
    /// Exact selected isochoric heat-capacity role and admissibility.
    pub const fn cv_requirement(&self) -> &ScalarPropertyRequirement {
        &self.cv_requirement
    }
    /// Checked acoustic/caloric/transport response at this source state.
    pub const fn state(&self) -> &ThermoviscousFluidState {
        &self.state
    }
    /// Identity of original source/query receipts, selected roles and response.
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    /// Refresh the complete original bundle atomically. Unsupported pressure,
    /// temperature, phase or composition does not change this admitted state.
    pub fn at_query_point(
        &self,
        card: &MaterialCard,
        point: &QueryPoint,
    ) -> Result<Self, MaterialStatePointError> {
        Self::try_from_acoustic(
            self.acoustic_liquid.at_query_point(card, point)?,
            &self.cv_requirement,
        )
    }
}
impl FluidAcoustics for ResolvedThermoviscousLiquidState {
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        Ok(self.state.acoustic)
    }
}
impl FluidTransport for ResolvedThermoviscousLiquidState {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        Ok(self.state.transport)
    }
}
impl ThermalFluidTransport for ResolvedThermoviscousLiquidState {
    fn thermal_transport_state(&self) -> Result<ThermalFluidTransportState, MaterialError> {
        self.acoustic_liquid.thermal_transport_state()
    }
}
impl ThermoviscousFluid for ResolvedThermoviscousLiquidState {
    fn thermoviscous_state(&self) -> Result<ThermoviscousFluidState, MaterialError> {
        Ok(self.state)
    }
}

/// Resolve canonical liquid transport, an explicit acoustic parameter and Cv at
/// one original query. An absent Cv claim is an error, never a gas-law fallback.
pub fn resolve_thermoviscous_liquid_state(
    card: &MaterialCard,
    point: &QueryPoint,
    parameter: &FluidAcousticParameter,
    viscosity_basis: LiquidViscosityBasis,
    selection: MaterialPropertySelection,
) -> Result<ResolvedThermoviscousLiquidState, MaterialStatePointError> {
    resolve_thermoviscous_liquid_state_with_requirements(
        card,
        point,
        parameter,
        &LiquidPropertyRequirements::canonical(viscosity_basis)?,
        &positive_requirement(
            SPECIFIC_HEAT_CAPACITY_CV_PROPERTY,
            SPECIFIC_HEAT_CAPACITY_DIMS,
        )?,
        selection,
    )
}

/// Select exact original source keys for all six physical roles together.
pub fn resolve_thermoviscous_liquid_state_with_requirements(
    card: &MaterialCard,
    point: &QueryPoint,
    parameter: &FluidAcousticParameter,
    liquid_requirements: &LiquidPropertyRequirements,
    cv_requirement: &ScalarPropertyRequirement,
    selection: MaterialPropertySelection,
) -> Result<ResolvedThermoviscousLiquidState, MaterialStatePointError> {
    resolve_bundle(
        card,
        point,
        parameter,
        liquid_requirements,
        cv_requirement,
        selection,
        false,
    )
}

/// Resolve transport, adiabatic acoustics, Cv and known signed volumetric
/// expansion in ONE original bundle. Full pressure/velocity/temperature
/// consumers need this extra derivative; Cp/Cv alone cannot determine its sign.
/// No expansion is inferred from the material name or from an ideal-gas law.
/// The original liquid owner admits its dimensional/quantity meaning and
/// preserves it through the same `at_query_point` refresh as the other roles.
pub fn resolve_thermoviscous_liquid_state_with_expansion(
    card: &MaterialCard,
    point: &QueryPoint,
    parameter: &FluidAcousticParameter,
    viscosity_basis: LiquidViscosityBasis,
    selection: MaterialPropertySelection,
) -> Result<ResolvedThermoviscousLiquidState, MaterialStatePointError> {
    resolve_bundle(
        card,
        point,
        parameter,
        &LiquidPropertyRequirements::canonical(viscosity_basis)?,
        &positive_requirement(
            SPECIFIC_HEAT_CAPACITY_CV_PROPERTY,
            SPECIFIC_HEAT_CAPACITY_DIMS,
        )?,
        selection,
        true,
    )
}

#[allow(clippy::too_many_arguments)]
fn resolve_bundle(
    card: &MaterialCard,
    point: &QueryPoint,
    parameter: &FluidAcousticParameter,
    liquid_requirements: &LiquidPropertyRequirements,
    cv_requirement: &ScalarPropertyRequirement,
    selection: MaterialPropertySelection,
    expansion: bool,
) -> Result<ResolvedThermoviscousLiquidState, MaterialStatePointError> {
    require_scalar(parameter.requirement(), parameter.dims())?;
    require_cv(liquid_requirements, cv_requirement)?;
    let mut requirements = liquid_requirements.properties().to_vec();
    requirements.push(parameter.requirement().clone());
    requirements.push(cv_requirement.clone());
    if expansion {
        requirements.push(ScalarPropertyRequirement::try_new(
            VOLUMETRIC_EXPANSION_PROPERTY,
            VOLUMETRIC_EXPANSION_DIMS,
            ScalarAdmissibility::Finite,
        )?);
    }
    let material = resolve_material_state_point(card, point, &requirements, selection)?;
    let acoustic = ResolvedFluidAcousticState::try_from_material(
        material,
        &liquid_requirements.properties()[0],
        parameter,
    )?;
    ResolvedThermoviscousLiquidState::try_from_acoustic(
        ResolvedAcousticLiquidState::admit(acoustic, liquid_requirements)?,
        cv_requirement,
    )
}

fn require_cv(
    liquid: &LiquidPropertyRequirements,
    cv: &ScalarPropertyRequirement,
) -> Result<(), MaterialStatePointError> {
    require_scalar(cv, SPECIFIC_HEAT_CAPACITY_DIMS)?;
    if liquid.properties().iter().any(|r| r.name() == cv.name()) {
        return Err(MaterialStatePointError::InvalidRequirement {
            property: cv.name().into(),
            reason: "Cv must be an independently selected source role, distinct from Cp and transport keys",
        });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::thermochemical_gas::nonideal::{
        PengRobinsonSpeciesBinding, ThermochemicalPengRobinsonModel,
    };
    use crate::thermochemical_gas::tests::{model, point};
    use fs_qty::{Pressure, Temperature};
    use fs_thermochem::SpeciesId;
    use fs_thermochem::eos::peng_robinson::PengRobinsonPhaseV1;

    #[test]
    fn real_gas_thermoviscous_response_keeps_its_actual_caloric_departure() {
        let bindings = vec![
            PengRobinsonSpeciesBinding {
                species: SpeciesId::new("B").unwrap(),
                critical_temperature: Temperature::new(150.687),
                critical_pressure: Pressure::new(4.863e6),
                acentric_factor: -0.00219,
            },
            PengRobinsonSpeciesBinding {
                species: SpeciesId::new("A").unwrap(),
                critical_temperature: Temperature::new(190.564),
                critical_pressure: Pressure::new(4.5992e6),
                acentric_factor: 0.01142,
            },
        ];
        let model =
            ThermochemicalPengRobinsonModel::new(model(vec![0.75, 0.25], 0.0), bindings, vec![])
                .unwrap();
        let dense = model
            .evaluate_real_gas(&point(300.0, 5.0e6), PengRobinsonPhaseV1::Stable)
            .unwrap();
        let state = dense.thermoviscous_state().unwrap();
        assert_eq!(state.acoustic().density_kg_m3(), dense.density_kg_m3());
        assert_eq!(state.acoustic().sound_speed_m_s(), dense.sound_speed_m_s());
        assert_eq!(
            state.transport().specific_heat_j_kg_k(),
            dense.specific_heat_cp_j_kg_k()
        );
        assert_eq!(
            state.specific_heat_cv_j_kg_k(),
            dense.specific_heat_cv_j_kg_k()
        );
        assert!((state.heat_capacity_ratio() / dense.gamma() - 1.0).abs() < 2e-14);
        assert_ne!(
            state.acoustic().density_kg_m3(),
            dense.source().state().density
        );
        assert_ne!(
            state.specific_heat_cv_j_kg_k(),
            dense.source().state().specific_heat_cp / dense.source().state().gamma
        );
    }
}
