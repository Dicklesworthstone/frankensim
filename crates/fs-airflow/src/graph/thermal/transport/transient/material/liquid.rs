//! Original liquid calorics and transport on distributed finite graph storage.

use super::convection::{GraphConvectionError, GraphConvectiveLink, GraphSourceHeat};
use super::*;
use crate::graph::thermal::medium::series::enthalpy::admit_liquid_caloric_source;
use fs_material::liquid::ResolvedLiquidState;
use fs_material::liquid_temperature::LiquidTemperatureLaw;
use fs_material::phase::ResolvedHeatCapacityCurve;

/// Prepared original source and exact registered caloric object. The original
/// temperature law retains material-card selection, pressure and every
/// nonthermal coordinate; no numerical Cp(T)*T replacement or new claim occurs.
pub struct LiquidGraphSources<'a, 'card> {
    law: &'a LiquidTemperatureLaw<'card>,
    caloric: &'a ResolvedHeatCapacityCurve,
    core: BoundGraphCore<'a>,
}

impl<'a, 'card> LiquidGraphSources<'a, 'card> {
    pub fn try_new(
        cx: &Cx<'_>,
        law: &'a LiquidTemperatureLaw<'card>,
        caloric: &'a ResolvedHeatCapacityCurve,
        thermal: &DepositionSimulation<'_>,
        binding: SourceGraphBinding<'a>,
        limits: SourceGraphLimits,
    ) -> Result<Self, SourceGraphError> {
        poll(cx)?;
        admit_liquid_caloric_source(cx, law.reference_state(), caloric, law.temperature_axis())
            .map_err(SourceGraphError::LiquidCaloric)?;
        let reference = law
            .reference_state()
            .transport_state()
            .map_err(SourceGraphError::Material)?;
        let core =
            BoundGraphCore::try_new(cx, caloric.curve(), reference, thermal, binding, limits)?;
        Ok(Self { law, caloric, core })
    }

    pub const fn law(&self) -> &'a LiquidTemperatureLaw<'card> {
        self.law
    }
    pub const fn caloric(&self) -> &'a ResolvedHeatCapacityCurve {
        self.caloric
    }
    pub fn inventory(&self) -> &'a DistributedGraphInventory {
        self.core.inventory
    }

    /// Re-evaluate every native junction and pipe temperature through the
    /// original source. Unsupported viscosity/density/transport is a refusal,
    /// including at zero flow; no previous snapshot or speed floor is reused.
    pub fn resolve(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
        trial_volume_flows_m3_s: &[f64],
    ) -> Result<ResolvedSourceGraph<ResolvedLiquidState>, SourceGraphError> {
        self.core
            .resolve(cx, thermal, trial_volume_flows_m3_s, |t| {
                self.law.at_temperature(t).map_err(SourceGraphError::Liquid)
            })
    }

    /// Run this with the actual native PHS report inside the driver's exchange
    /// callback. Original support is checked at average, old, new and crossing
    /// rates. Native graph and original section admission are both necessary:
    /// aggregated gravity can hide cancelling local work.
    pub fn admit_actual(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
        sources: &GraphSourceState<ResolvedLiquidState>,
        report: &InertialGraphReport,
    ) -> Result<GraphSourceAdmission, SourceGraphError> {
        self.core.admit_actual(cx, thermal, sources, report)
    }

    /// Admit actual/envelope hydraulics, then evaluate reciprocal finite-cell
    /// convection from the same original liquid source at bulk/film temperature.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_links(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
        sources: &GraphSourceState<ResolvedLiquidState>,
        report: &InertialGraphReport,
        bindings: &[GraphConvectiveLink],
        maximum_links: usize,
    ) -> Result<GraphSourceHeat<ResolvedLiquidState>, GraphConvectionError> {
        let admission = self
            .admit_actual(cx, thermal, sources, report)
            .map_err(GraphConvectionError::Source)?;
        let convection = convection::evaluate_links(
            &self.core,
            cx,
            thermal,
            report,
            bindings,
            maximum_links,
            |t| self.law.at_temperature(t).map_err(SourceGraphError::Liquid),
        )?;
        Ok(GraphSourceHeat {
            admission,
            convection,
        })
    }
}
