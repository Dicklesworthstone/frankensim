//! Material-bound series exchange using the existing single-channel compiler.
//!
//! There is one inlet and one conserved mass/capacity rate. A local wall region
//! is not a fresh coolant supply. Material properties remain frozen; this does
//! not infer entrance development, phase validity or hydraulic redistribution.
//! The explicit [`enthalpy`] march instead consumes a source-backed liquid
//! caloric chart, while retaining the same frozen hydraulic and convection data.

pub mod enthalpy;
pub mod hydraulic;

use fs_convection::{
    CorrelationId, CorrelationRobinBoundary, HeatTransferCoefficient,
    NusseltEvaluation, ThermalConductivity, ThermalDirection,
};
use fs_evidence::Evidence;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_qty::{Temperature, VolumetricFlowRate};

use crate::conjugate::AirPath;
use super::{FlowHeatGeometry, FrozenFlowMedium, MaterialFlowChannel, MediumFlowError, positive};

/// Bound on synchronous compilation, checked before allocation.
pub const MAX_MATERIAL_FLOW_SECTIONS: usize = 1024;

/// One independently described solid/fluid exchange, in stream-wise order.
/// Each correlation's local geometric and entrance assumptions remain explicit;
/// placing sections in series does not reset a physical entrance length.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FlowHeatSection<'a> {
    /// Unique region name accepted by the existing solid/air coupling driver.
    pub region: &'a str,
    /// Local flow and wall geometry.
    pub geometry: FlowHeatGeometry,
    /// Explicit correlation; not selected from the material name.
    pub correlation: CorrelationId,
    /// Wall-to-fluid heating/cooling convention used by the selected law.
    pub direction: ThermalDirection,
}

impl<'a> FrozenFlowMedium<'a> {
    /// Compile several heated regions into ONE conserved-flow path.
    ///
    /// The inlet is applied only once. The existing `AirPath::march` passes each
    /// outlet to the next region, so downstream heating is never computed from
    /// a fresh cold inlet. Local area/diameter changes update velocity, Re and h;
    /// mass flow and heat capacity remain common to the whole series. Each
    /// region retains its own correlation and material binding for Robin rows.
    /// Feed `path()` directly into the existing conjugate solvers.
    ///
    /// # Errors
    /// Refuses empty/oversized inputs before allocation, duplicate region names,
    /// unsupported correlations, invalid geometry and all material/solver errors.
    /// A late region refusal returns no partial path and mutates no input.
    pub fn channel_series(
        self, flow: VolumetricFlowRate, inlet: Temperature, sections: &[FlowHeatSection<'_>],
    ) -> Result<MaterialFlowPath<'a>, MediumFlowError> {
        if sections.is_empty() || sections.len() > MAX_MATERIAL_FLOW_SECTIONS {
            return Err(MediumFlowError::Invalid {
                field: "section count must be in 1..=MAX_MATERIAL_FLOW_SECTIONS",
                value_bits: (sections.len() as f64).to_bits(),
            });
        }
        let state = self.transport_state()?;
        let q = positive("volume flow", flow.value())?;
        positive("absolute inlet temperature", inlet.value())?;
        let mass_flow = positive("mass flow", state.density_kg_m3() * q)?;
        let mut exchanges = Vec::with_capacity(sections.len());
        let mut segments = Vec::with_capacity(sections.len());
        for section in sections {
            // Reuse the actual production compiler and all of its source,
            // geometry, constant-wall-temperature and correlation admission.
            // Its isolated probe path is not exposed as another coolant inlet.
            let MaterialFlowChannel { medium, transport, geometry, evaluation, coefficient, path } =
                self.channel(section.region, flow, inlet, section.geometry, section.correlation, section.direction)?;
            segments.extend_from_slice(path.segments());
            exchanges.push(MaterialHeatExchange { medium, transport, geometry, evaluation, coefficient });
        }
        // The existing path owner checks duplicate/invalid region names and
        // total exchange area/capacity. Do not implement a second admission rule.
        let path = AirPath::new(inlet.value(), mass_flow, state.specific_heat_j_kg_k(), segments)?;
        Ok(MaterialFlowPath { medium: self, exchanges, path, mass_flow_kg_per_s: mass_flow })
    }
}

/// A local heat-transfer law retaining its original material and correlation.
/// This is not a separate inlet or mass stream.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialHeatExchange<'a> {
    medium: FrozenFlowMedium<'a>,
    transport: FluidTransportState,
    geometry: FlowHeatGeometry,
    evaluation: NusseltEvaluation,
    coefficient: Evidence<HeatTransferCoefficient>,
}
impl<'a> MaterialHeatExchange<'a> {
    /// Original source medium, with receipts for source-resolved variants.
    #[must_use]
    pub const fn medium(&self) -> FrozenFlowMedium<'a> { self.medium }
    /// The same checked snapshot used by all other flow and heat calculations.
    #[must_use]
    pub const fn transport_state(&self) -> &FluidTransportState { &self.transport }
    /// Local flow and wall geometry.
    #[must_use]
    pub const fn geometry(&self) -> &FlowHeatGeometry { &self.geometry }
    /// Original correlation evidence; material uncertainty is not promoted into it.
    #[must_use]
    pub const fn coefficient(&self) -> &Evidence<HeatTransferCoefficient> { &self.coefficient }
    /// Original Nusselt evaluation, including admitted dimensionless groups.
    #[must_use]
    pub const fn correlation(&self) -> &NusseltEvaluation { &self.evaluation }
    /// Lower this coefficient with the effective reference supplied by a
    /// conjugate solve. Conductivity and length cannot be independently replaced.
    pub fn robin_boundary(&self, reference: Temperature) -> Result<CorrelationRobinBoundary, MediumFlowError> {
        positive("absolute Robin reference temperature", reference.value())?;
        Ok(self.evaluation.robin_boundary(
            ThermalConductivity::new(self.transport.thermal_conductivity_w_m_k()),
            self.geometry.characteristic_length, reference,
        )?)
    }
}

/// One inlet and one capacity rate across multiple source-bound wall regions.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialFlowPath<'a> {
    medium: FrozenFlowMedium<'a>,
    exchanges: Vec<MaterialHeatExchange<'a>>,
    path: AirPath,
    mass_flow_kg_per_s: f64,
}
impl<'a> MaterialFlowPath<'a> {
    /// One conserved mass flow used by the complete path, kg/s.
    #[must_use]
    pub const fn mass_flow_kg_per_s(&self) -> f64 { self.mass_flow_kg_per_s }
    /// The complete path, directly accepted by existing conjugate solvers.
    #[must_use]
    pub const fn path(&self) -> &AirPath { &self.path }
    /// Original source medium; no detached row gains its authority.
    #[must_use]
    pub const fn medium(&self) -> FrozenFlowMedium<'a> { self.medium }
    /// Per-region laws in exactly `path().segments()` order.
    #[must_use]
    pub fn exchanges(&self) -> &[MaterialHeatExchange<'a>] { &self.exchanges }
}
