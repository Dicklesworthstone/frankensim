//! Work from solved pressure and heat-transfer inputs cannot diverge silently.
use super::SeriesPressureFlow;
use crate::graph::thermal::medium::convection::{
    ForcedConvectionSettings, MaterialConvection, MaterialConvectionError,
    evaluate_forced_from_material,
};
use fs_convection::CorrelationId;
use fs_material::fluid_transport::FluidTransportState;
use fs_qty::Length;

#[derive(Debug)]
pub enum PressureConvectionError {
    InvalidSection(usize),
    Convection(MaterialConvectionError),
}
impl core::fmt::Display for PressureConvectionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "pressure-derived convection: {self:?}")
    }
}
impl std::error::Error for PressureConvectionError {}

/// Original pressure solution and existing convection result evaluated with
/// its EXACT material snapshot. The snapshot alone gains no source authority;
/// the borrowed pressure owner retains the original typed medium and geometry.
#[derive(Debug)]
pub struct PressureConvection<'a, F> {
    flow: &'a SeriesPressureFlow<F>,
    section: usize,
    speed_m_s: f64,
    convection: MaterialConvection<'a, FluidTransportState>,
}
impl<'a, F> PressureConvection<'a, F> {
    pub const fn flow(&self) -> &'a SeriesPressureFlow<F> { self.flow }
    pub const fn section(&self) -> usize { self.section }
    pub fn medium(&self) -> &F { &self.flow.sections()[self.section].medium }
    pub const fn speed_m_s(&self) -> f64 { self.speed_m_s }
    /// Original correlation, evidence, h and Robin-boundary lowering. The
    /// material input referenced here is the retained hydraulic snapshot.
    pub const fn convection(&self) -> &MaterialConvection<'a, FluidTransportState> { &self.convection }
}

impl<F> SeriesPressureFlow<F> {
    /// Lower ONE solved section into the existing forced-convection evaluator.
    ///
    /// Speed is |local volume flow| / ORIGINAL flow area; callers cannot enter
    /// an unrelated velocity, density, viscosity or conductivity. No original
    /// source is queried a second time. Changing source data and solving again
    /// changes hydraulic flow, Reynolds, Prandtl and heat transfer together.
    /// The result retains a borrow of this complete original pressure solution.
    ///
    /// Correlation length and settings remain explicit HEAT geometry. They need
    /// not equal the hydraulic Reynolds length, but the caller must justify
    /// their physical applicability. The chosen correlation retains its own
    /// support, including any zero-flow refusal; no epsilon velocity, fallback
    /// correlation or heat-transfer coefficient is fabricated.
    ///
    /// This uses a single frozen state for both physics consumers. It does NOT
    /// query a new film temperature, integrate a thermal state, deposit pressure
    /// losses, solve a pump or refresh an evolving material history.
    pub fn forced_convection(
        &self,
        section: usize,
        correlation: CorrelationId,
        characteristic_length: Length,
        settings: ForcedConvectionSettings,
    ) -> Result<PressureConvection<'_, F>, PressureConvectionError> {
        let state = self.evaluations().get(section)
            .ok_or(PressureConvectionError::InvalidSection(section))?;
        let speed_m_s = state.volume_flow_m3_s.abs() / state.law.flow_area().value();
        let convection = evaluate_forced_from_material(correlation, &state.transport,
            speed_m_s, characteristic_length, settings).map_err(PressureConvectionError::Convection)?;
        Ok(PressureConvection { flow: self, section, speed_m_s, convection })
    }
}
