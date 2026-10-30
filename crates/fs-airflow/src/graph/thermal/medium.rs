//! One material state supplies hydraulic loading, convection and advected heat.
//!
//! These are frozen-property comparisons, not a temperature/composition update
//! law: replacing the medium rebuilds the loss and channel coefficients. The
//! existing flow, correlation, Robin and exponential channel owners still own
//! their equations and admission. In particular, a material swap cannot turn
//! an out-of-domain correlation into a supported one. Fan density similarity
//! requires an explicit reference medium and applicability through [`fan`];
//! no pump correction or liquid phase validity is inferred.

pub mod convection;
pub mod fan;
pub mod series;
pub mod steady_gas;
pub use series::{FlowHeatSection, MaterialFlowPath, MaterialHeatExchange, MAX_MATERIAL_FLOW_SECTIONS};

use fs_convection::{
    CorrelationError, CorrelationId, CorrelationRobinBoundary, HeatTransferCoefficient,
    NusseltEvaluation, ThermalConductivity, ThermalDirection,
};
use fs_evidence::Evidence;
use fs_material::MaterialError;
use fs_material::fluid_transport::{
    FluidTransport, FluidTransportState, ThermalFluidTransport, ThermalFluidTransportState,
};
use fs_material::gas::{GasState, ResolvedGasState};
use fs_material::gas_mixture::ResolvedGasMixtureState;
use fs_material::liquid::ResolvedLiquidState;
use fs_material::moist_air::ResolvedMoistAirState;
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;
use fs_material::thermochemical_gas::nonideal::ResolvedThermochemicalRealGasState;
use fs_qty::{Area, Density, Length, Temperature, VolumetricFlowRate};

use crate::conjugate::{AirPath, AirSegment};
use crate::{AirflowError, LossElement, sharp_edged_orifice_loss};
use super::transport::TransportAir;
use convection::{ForcedConvectionSettings, MaterialConvectionError, evaluate_forced_from_material};

/// Borrow the actual material authority; no second property registry or presets.
///
/// Source-resolved variants retain their full bundles. `DeclaredGas` retains
/// only authored/evaluated coefficients and does not acquire source evidence.
/// The legacy `AirPath`/`TransportAir` names do not impose air constants: their
/// constant-density, constant-cp sensible-heat equations also apply to an
/// explicitly single-phase liquid. Phase validity remains the caller's duty.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum FrozenFlowMedium<'a> {
    /// A single-phase liquid, including either source viscosity basis.
    Liquid(&'a ResolvedLiquidState),
    /// A source-resolved gas and its transport-model choice.
    Gas(&'a ResolvedGasState),
    /// A source-resolved arbitrary frozen gas composition.
    GasMixture(&'a ResolvedGasMixtureState),
    /// The source-resolved humid-air specialization.
    MoistAir(&'a ResolvedMoistAirState),
    /// Shared species thermodynamics and transport at one frozen composition.
    /// The original species, caloric, EOS and transport receipts remain bound
    /// to the hydraulic and thermal consumers; no gas-spec reconstruction.
    ThermochemicalGas(&'a ResolvedThermochemicalGasState),
    /// A selected stable nonideal EOS branch with its original NASA caloric
    /// reference and transport materials. The frozen-property network consumes
    /// the same PR density/Cp as convection; it adds no phase-change or
    /// compressible-flow capability to the native channel approximation.
    ThermochemicalRealGas(&'a ResolvedThermochemicalRealGasState),
    /// Explicitly authored gas coefficients; not promoted to sourced data.
    DeclaredGas(&'a GasState),
}

impl FluidTransport for FrozenFlowMedium<'_> {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        match self {
            Self::Liquid(value) => value.transport_state(),
            Self::Gas(value) => value.transport_state(),
            Self::GasMixture(value) => value.transport_state(),
            Self::MoistAir(value) => value.transport_state(),
            Self::ThermochemicalGas(value) => value.transport_state(),
            Self::ThermochemicalRealGas(value) => value.transport_state(),
            Self::DeclaredGas(value) => value.transport_state(),
        }
    }
}

impl ThermalFluidTransport for FrozenFlowMedium<'_> {
    fn thermal_transport_state(&self) -> Result<ThermalFluidTransportState, MaterialError> {
        match self {
            Self::Liquid(value) => value.thermal_transport_state(),
            Self::Gas(value) => value.thermal_transport_state(),
            Self::GasMixture(value) => value.thermal_transport_state(),
            Self::MoistAir(value) => value.thermal_transport_state(),
            Self::ThermochemicalGas(value) => value.thermal_transport_state(),
            Self::ThermochemicalRealGas(value) => value.thermal_transport_state(),
            Self::DeclaredGas(value) => value.thermal_transport_state(),
        }
    }
}

/// Independently authored flow and heat-exchange geometry in coherent SI.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FlowHeatGeometry {
    /// Flow-normal area used by velocity, not the wetted area.
    pub flow_area: Area,
    /// Correlation length (hydraulic diameter for internal-flow correlations).
    pub characteristic_length: Length,
    /// Stream-wise heated length; supplies the correlation's L/D ratio.
    pub heated_length: Length,
    /// Actual wall area participating in the exchange.
    pub wetted_area: Area,
    /// Explicit rectangular aspect ratio when required by the chosen card.
    pub aspect_ratio: Option<f64>,
}

/// Material, geometry, correlation or downstream solver refusal.
#[derive(Debug)]
pub enum MediumFlowError {
    /// Invalid physical input or unrepresentable derived value.
    Invalid { field: &'static str, value_bits: u64 },
    /// Shared material-transport admission refused this state.
    Material(MaterialError),
    /// A uniform-heat-flux correlation is incompatible with this constant-wall-temperature channel.
    UnsupportedCorrelation(CorrelationId),
    /// The existing correlation refused this material/flow regime.
    Correlation(CorrelationError),
    /// The existing loss/channel owner refused admission.
    Airflow(AirflowError),
}
impl core::fmt::Display for MediumFlowError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material flow coupling: {self:?}")
    }
}
impl std::error::Error for MediumFlowError {}
impl From<MaterialError> for MediumFlowError {
    fn from(error: MaterialError) -> Self { Self::Material(error) }
}
impl From<MaterialConvectionError> for MediumFlowError {
    fn from(error: MaterialConvectionError) -> Self {
        match error {
            MaterialConvectionError::Material(error) => Self::Material(error),
            MaterialConvectionError::Correlation(error) => Self::Correlation(error),
        }
    }
}
impl From<CorrelationError> for MediumFlowError {
    fn from(error: CorrelationError) -> Self { Self::Correlation(error) }
}
impl From<AirflowError> for MediumFlowError {
    fn from(error: AirflowError) -> Self { Self::Airflow(error) }
}
fn positive(field: &'static str, value: f64) -> Result<f64, MediumFlowError> {
    if value.is_finite() && value > 0.0 { Ok(value) }
    else { Err(MediumFlowError::Invalid { field, value_bits: value.to_bits() }) }
}

impl<'a> FrozenFlowMedium<'a> {
    /// Supply the SAME density and cp to the existing heat-transport network.
    /// A detached legacy tuple carries no material authority; retain this view.
    pub fn transport_air(self) -> Result<TransportAir, MediumFlowError> {
        let transport = self.transport_state()?;
        Ok(TransportAir {
            density: Density::new(transport.density_kg_m3()),
            specific_heat_j_kg_k: transport.specific_heat_j_kg_k(),
        })
    }

    /// Bind the existing sharp-edged-orifice loss to this medium's density.
    /// Its Reynolds plateau still requires the ordinary post-solve regime audit.
    /// The original orifice evidence and the separate material evidence survive.
    pub fn orifice_loss(self, name: &str, area: Area) -> Result<MaterialFlowLoss<'a>, MediumFlowError> {
        let transport = self.transport_state()?;
        let element = sharp_edged_orifice_loss(name, area, Density::new(transport.density_kg_m3()))?;
        Ok(MaterialFlowLoss { medium: self, element })
    }

    /// Build an actual forced-convection channel, not just sampled properties.
    ///
    /// Re, Pr, h and mass-flow heat capacity all come from this SAME snapshot.
    /// Flow is positive in the declared inlet-to-outlet direction. The chosen
    /// correlation retains all its original domain gates and model evidence.
    /// Coefficients are frozen at the resolved material point even though the
    /// downstream sensible-heat model computes a temperature rise; that is an
    /// explicit approximation, not source validity throughout the heated path.
    pub fn channel(
        self, region: &str, flow: VolumetricFlowRate, inlet: Temperature,
        geometry: FlowHeatGeometry, correlation: CorrelationId, direction: ThermalDirection,
    ) -> Result<MaterialFlowChannel<'a>, MediumFlowError> {
        if matches!(correlation, CorrelationId::CircularDuctLaminarChf | CorrelationId::RectangularDuctLaminarChf) {
            return Err(MediumFlowError::UnsupportedCorrelation(correlation));
        }
        let q = positive("volume flow", flow.value())?;
        let area = positive("flow area", geometry.flow_area.value())?;
        let length = positive("characteristic length", geometry.characteristic_length.value())?;
        let heated = positive("heated length", geometry.heated_length.value())?;
        positive("wetted area", geometry.wetted_area.value())?;
        let velocity = positive("velocity", q / area)?;
        let result = evaluate_forced_from_material(
            correlation, &self, velocity, geometry.characteristic_length,
            ForcedConvectionSettings {
                length_over_hydraulic_diameter: Some(positive("heated length ratio", heated / length)?),
                aspect_ratio: geometry.aspect_ratio,
                direction,
            },
        )?;
        let (transport, evaluation, coefficient) = result.into_parts();
        let segment = AirSegment::new(region, geometry.wetted_area.value(), coefficient.value.value())?;
        let path = AirPath::new(
            inlet.value(), positive("mass flow", transport.density_kg_m3() * q)?,
            transport.specific_heat_j_kg_k(), vec![segment],
        )?;
        Ok(MaterialFlowChannel { medium: self, transport, geometry, evaluation, coefficient, path })
    }
}

/// Existing hydraulic loss kept beside the actual medium that produced it.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialFlowLoss<'a> {
    medium: FrozenFlowMedium<'a>,
    element: LossElement,
}
impl<'a> MaterialFlowLoss<'a> {
    /// Underlying loss for the existing series/parallel or graph solver.
    #[must_use]
    pub const fn element(&self) -> &LossElement { &self.element }
    /// Source authority; cloning a detached element does not copy this binding.
    #[must_use]
    pub const fn medium(&self) -> FrozenFlowMedium<'a> { self.medium }
}

/// Source-bound convection, conduction-boundary lowering and transported heat.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialFlowChannel<'a> {
    medium: FrozenFlowMedium<'a>,
    transport: FluidTransportState,
    geometry: FlowHeatGeometry,
    evaluation: NusseltEvaluation,
    coefficient: Evidence<HeatTransferCoefficient>,
    path: AirPath,
}
impl<'a> MaterialFlowChannel<'a> {
    /// Actual channel consumed by the existing conjugate coupling drivers.
    #[must_use]
    pub const fn path(&self) -> &AirPath { &self.path }
    /// Original source medium (including full receipts for resolved variants).
    #[must_use]
    pub const fn medium(&self) -> FrozenFlowMedium<'a> { self.medium }
    /// Exact snapshot shared by momentum and heat-transfer coefficients.
    #[must_use]
    pub const fn transport(&self) -> &FluidTransportState { &self.transport }
    /// Original correlation evidence; material uncertainty is not promoted into it.
    #[must_use]
    pub const fn coefficient(&self) -> &Evidence<HeatTransferCoefficient> { &self.coefficient }
    /// Original Nusselt evaluation, including admitted dimensionless groups.
    #[must_use]
    pub const fn correlation(&self) -> &NusseltEvaluation { &self.evaluation }
    /// Lower exactly this coefficient through the existing conduction Robin seam.
    /// A conjugate solve supplies its effective reference, not a second coolant.
    pub fn robin_boundary(&self, reference: Temperature) -> Result<CorrelationRobinBoundary, MediumFlowError> {
        Ok(self.evaluation.robin_boundary(
            ThermalConductivity::new(self.transport.thermal_conductivity_w_m_k()),
            self.geometry.characteristic_length, reference,
        )?)
    }
}
