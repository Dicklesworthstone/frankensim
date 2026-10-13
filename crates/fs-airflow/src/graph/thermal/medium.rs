//! One material state supplies hydraulic loading, convection and advected heat.
//!
//! These are frozen-property comparisons, not a temperature/composition update
//! law: replacing the medium rebuilds the loss and channel coefficients. The
//! existing flow, correlation, Robin and exponential channel owners still own
//! their equations and admission. In particular, a material swap cannot turn
//! an out-of-domain correlation into a supported one. Nothing here silently
//! rescales a measured fan/pump curve or assumes a liquid cannot boil.

use fs_convection::{
    CorrelationError, CorrelationId, CorrelationInputs, CorrelationRobinBoundary,
    HeatTransferCoefficient, NusseltEvaluation, ThermalConductivity, ThermalDirection,
};
use fs_evidence::Evidence;
use fs_material::gas::{GasState, ResolvedGasState};
use fs_material::gas_mixture::ResolvedGasMixtureState;
use fs_material::liquid::ResolvedLiquidState;
use fs_material::moist_air::ResolvedMoistAirState;
use fs_qty::{Area, Density, Length, Temperature, VolumetricFlowRate};

use crate::conjugate::{AirPath, AirSegment};
use crate::{AirflowError, LossElement, sharp_edged_orifice_loss};
use super::transport::TransportAir;

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
    /// Explicitly authored gas coefficients; not promoted to sourced data.
    DeclaredGas(&'a GasState),
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
    // Return one coherent source tuple. Derived Pr is never a caller override.
    fn properties(self) -> Result<(f64, f64, f64, f64), MediumFlowError> {
        let (rho, cp, mu, k) = match self {
            Self::Liquid(liquid) => (
                liquid.density_kg_m3(), liquid.specific_heat_j_kg_k(),
                liquid.dynamic_viscosity_pa_s(), liquid.thermal_conductivity_w_m_k(),
            ),
            gas => {
                let state = match gas {
                    Self::Gas(value) => value.state(),
                    Self::GasMixture(value) => value.state(),
                    Self::MoistAir(value) => value.state(),
                    Self::DeclaredGas(value) => value,
                    Self::Liquid(_) => unreachable!("liquid handled above"),
                };
                (state.density, state.specific_heat_cp, state.dynamic_viscosity, state.thermal_conductivity)
            }
        };
        Ok((positive("density", rho)?, positive("specific heat", cp)?,
            positive("dynamic viscosity", mu)?, positive("thermal conductivity", k)?))
    }

    /// Supply the SAME density and cp to the existing heat-transport network.
    /// A detached legacy tuple carries no material authority; retain this view.
    pub fn transport_air(self) -> Result<TransportAir, MediumFlowError> {
        let (rho, cp, _, _) = self.properties()?;
        Ok(TransportAir { density: Density::new(rho), specific_heat_j_kg_k: cp })
    }

    /// Bind the existing sharp-edged-orifice loss to this medium's density.
    /// Its Reynolds plateau still requires the ordinary post-solve regime audit.
    /// The original orifice evidence and the separate material evidence survive.
    pub fn orifice_loss(self, name: &str, area: Area) -> Result<MaterialFlowLoss<'a>, MediumFlowError> {
        let (rho, _, _, _) = self.properties()?;
        let element = sharp_edged_orifice_loss(name, area, Density::new(rho))?;
        Ok(MaterialFlowLoss { medium: self, element })
    }

    /// Build an actual forced-convection channel, not just sampled properties.
    ///
    /// Re, Pr, h and mass-flow heat capacity all come from this SAME medium.
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
        let (rho, cp, mu, k) = self.properties()?;
        let q = positive("volume flow", flow.value())?;
        let area = positive("flow area", geometry.flow_area.value())?;
        let length = positive("characteristic length", geometry.characteristic_length.value())?;
        let heated = positive("heated length", geometry.heated_length.value())?;
        positive("wetted area", geometry.wetted_area.value())?;
        let velocity = positive("velocity", q / area)?;
        let reynolds = positive("Reynolds number", (rho / mu) * velocity * length)?;
        let prandtl = positive("Prandtl number", (mu / k) * cp)?;
        let mut inputs = CorrelationInputs::forced(reynolds, prandtl)
            .with_length_ratio(positive("heated length ratio", heated / length)?)
            .with_direction(direction);
        if let Some(aspect) = geometry.aspect_ratio {
            inputs = inputs.with_aspect_ratio(positive("aspect ratio", aspect)?);
        }
        let evaluation = fs_convection::evaluate(correlation, inputs)?;
        let coefficient = evaluation.heat_transfer_coefficient(
            ThermalConductivity::new(k), geometry.characteristic_length,
        )?;
        let segment = AirSegment::new(region, geometry.wetted_area.value(), coefficient.value.value())?;
        let path = AirPath::new(inlet.value(), positive("mass flow", rho * q)?, cp, vec![segment])?;
        Ok(MaterialFlowChannel { medium: self, geometry, evaluation, coefficient, path })
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
    /// Original correlation evidence; material uncertainty is not promoted into it.
    #[must_use]
    pub const fn coefficient(&self) -> &Evidence<HeatTransferCoefficient> { &self.coefficient }
    /// Original Nusselt evaluation, including admitted dimensionless groups.
    #[must_use]
    pub const fn correlation(&self) -> &NusseltEvaluation { &self.evaluation }
    /// Lower exactly this coefficient through the existing conduction Robin seam.
    /// A conjugate solve supplies its effective reference, not a second coolant.
    pub fn robin_boundary(&self, reference: Temperature) -> Result<CorrelationRobinBoundary, MediumFlowError> {
        let (_, _, _, k) = self.medium.properties()?;
        Ok(self.evaluation.robin_boundary(ThermalConductivity::new(k), self.geometry.characteristic_length, reference)?)
    }
}
