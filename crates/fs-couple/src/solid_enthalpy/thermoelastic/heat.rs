//! Complete original material heat transfer and reciprocal native solid motion.
//!
//! Conduction, thermal interfaces, integrated layers and enclosure radiation
//! retain the prepared heat network. Equilibrium or Newmark inertia and local
//! reversible heat retain the original thermoelastic trajectory. Both see the
//! same finite caloric cells and source context; neither is advanced twice.

mod assembly;
pub use assembly::ThermoelasticHeatAssemblyError;

use super::dynamic::{DynamicThermoelasticCandidate, DynamicThermoelasticSolid};
use super::*;
use crate::circuit_enthalpy::contacts::heat::{
    HeatNetworkError, MaterialHeatObserver, PreparedHeatNetwork,
    solid::{SolidHeatAudit, closed_account},
};
use fs_thermal::observers::{CoupledObservers, ObserverCouplingError};

/// The original failing physical owner is retained, including observer failures.
#[derive(Debug)]
pub enum ThermoelasticHeatError<E = Infallible> {
    Heat(HeatNetworkError),
    Binding(SolidEnthalpyError),
    Coupled(ThermoelasticFeedbackError<ObserverCouplingError<HeatNetworkError, E>>),
    Invalid(&'static str),
    Cancelled,
}
impl<E: core::fmt::Debug> core::fmt::Display for ThermoelasticHeatError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "complete material heat/reciprocal solid: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for ThermoelasticHeatError<E> {}

// Caloric identity alone is insufficient: two valid queries of one card can
// refer to different pressure, composition, processing or reference conditions.
// Conductivity/source selection remains entirely with the prepared network.
fn admit<E>(
    source: &SolidEnthalpy<'_>,
    thermal: &DepositionSimulation<'_>,
    prepared: &PreparedHeatNetwork<'_, '_>,
    cx: &Cx<'_>,
) -> Result<(), ThermoelasticHeatError<E>> {
    let check = || {
        cx.checkpoint()
            .map_err(|_| ThermoelasticHeatError::Cancelled)
    };
    check()?;
    let network = prepared
        .thermal_material_network()
        .map_err(ThermoelasticHeatError::Heat)?;
    for expected in &source.cells {
        check()?;
        let native = &expected.source;
        let actual = network
            .cells
            .iter()
            .find(|cell| cell.cell == native.cell)
            .ok_or(ThermoelasticHeatError::Invalid(
                "reciprocal solid cell lacks its original heat-network binding",
            ))?;
        if expected.reference_cell.material_id != actual.material_id
            || !std::ptr::eq(native.curve, actual.curve)
            || native.point.material.content_hash() != actual.point.material.content_hash()
            || native.point.temperature_axis != actual.point.temperature_axis
            || native.point.liquid_fraction_axis != actual.point.liquid_fraction_axis
        {
            return Err(ThermoelasticHeatError::Invalid(
                "heat and reciprocal solid require the same material and caloric owner",
            ));
        }
        let cell = thermal
            .cells()
            .get(&native.cell)
            .ok_or(ThermoelasticHeatError::Invalid(
                "inactive reciprocal solid thermal cell",
            ))?;
        let phase = native
            .curve
            .state_at_specific_enthalpy(cell.specific_h_j_kg)
            .map_err(|error| ThermoelasticHeatError::Binding(error.into()))?;
        let mechanical = native
            .point
            .at(phase)
            .map_err(|error| ThermoelasticHeatError::Binding(error.into()))?;
        let transport = actual
            .point
            .at(phase)
            .map_err(|error| ThermoelasticHeatError::Binding(error.into()))?;
        if mechanical != transport {
            return Err(ThermoelasticHeatError::Invalid(
                "heat and reciprocal solid material conditions disagree",
            ));
        }
    }
    check()
}

impl<'s, 'm> ThermoelasticSolid<'s, 'm> {
    /// Advance compatible elastic equilibrium and the complete original heat
    /// network as one source-consistent candidate. Each fixed-point trial starts
    /// from the same accepted thermal/observer history; local reversible heat,
    /// source potentials and actual nodal work keep their native accounting.
    ///
    /// The additional closed audit removes enclosure roundoff from external
    /// heat and admits both corrected residuals under the SAME original joule
    /// tolerance. Internal radiation is not a second energy source. Electrical
    /// bindings refuse because this trajectory does not own circuit work.
    ///
    /// `inner` is the caller's independent observer, not a previously returned
    /// composed heat observer. Retain its `second` history on continuation.
    /// Existing fixed geometry/supports, sub-solidus small-strain/source ranges
    /// and explicit thermal-substep accuracy remain unchanged. No caller-owned
    /// state is changed on source refusal, cancellation or a failed energy gate.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn advance_heat_network_candidate<'t, 'p, 'n, 'a, 'scope, O: DepositionObserver>(
        &self,
        thermal: &DepositionSimulation<'t>,
        prepared: &'p PreparedHeatNetwork<'n, 'a>,
        inner: &O,
        target_s: f64,
        loads: &TetLoadCase3<'_, '_>,
        thermal_config: StepConfig,
        config: ThermoelasticSolve,
        cx: &'p Cx<'scope>,
    ) -> Result<
        (
            ThermoelasticCandidate<
                's,
                'm,
                't,
                CoupledObservers<MaterialHeatObserver<'p, 'n, 'a, 'scope>, O>,
            >,
            SolidHeatAudit,
        ),
        ThermoelasticHeatError<O::Error>,
    > {
        admit(self.source, thermal, prepared, cx)?;
        let observer = prepared
            .observe(cx, inner.clone(), config.energy_tolerance_j)
            .map_err(ThermoelasticHeatError::Heat)?;
        let candidate = self
            .advance(
                thermal,
                &observer,
                target_s,
                loads,
                thermal_config,
                config,
                cx,
            )
            .map_err(ThermoelasticHeatError::Coupled)?;
        let report = &candidate.report;
        let audit = closed_account(
            report.external_heat_j,
            report.thermal_residual_j,
            report.combined_residual_j,
            candidate.observer.first.closed_radiation_net_j(),
            config.energy_tolerance_j,
        )
        .map_err(ThermoelasticHeatError::Invalid)?;
        cx.checkpoint()
            .map_err(|_| ThermoelasticHeatError::Cancelled)?;
        Ok((candidate, audit))
    }
}

impl<'s, 'm> DynamicThermoelasticSolid<'s, 'm> {
    /// Continue actual Newmark inertia, full thermoelastic feedback and original
    /// material heat transport together. The loaded displacement, velocity,
    /// acceleration, source potentials and step counter are preserved; this does
    /// not replace a vibrating state with a new static equilibrium.
    ///
    /// The same source admission and closed enclosure balance as the static
    /// entry point apply. Native kinetic storage and exported Rayleigh loss stay
    /// in the original total balance; neither is silently deposited as heat.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn advance_heat_network_candidate<'t, 'p, 'n, 'a, 'scope, O: DepositionObserver>(
        &self,
        thermal: &DepositionSimulation<'t>,
        prepared: &'p PreparedHeatNetwork<'n, 'a>,
        inner: &O,
        target_s: f64,
        loads: &TetLoadCase3<'_, '_>,
        thermal_config: StepConfig,
        config: ThermoelasticSolve,
        cx: &'p Cx<'scope>,
    ) -> Result<
        (
            DynamicThermoelasticCandidate<
                's,
                'm,
                't,
                CoupledObservers<MaterialHeatObserver<'p, 'n, 'a, 'scope>, O>,
            >,
            SolidHeatAudit,
        ),
        ThermoelasticHeatError<O::Error>,
    > {
        admit(self.source, thermal, prepared, cx)?;
        let observer = prepared
            .observe(cx, inner.clone(), config.energy_tolerance_j)
            .map_err(ThermoelasticHeatError::Heat)?;
        let candidate = self
            .advance(
                thermal,
                &observer,
                target_s,
                loads,
                thermal_config,
                config,
                cx,
            )
            .map_err(ThermoelasticHeatError::Coupled)?;
        let report = &candidate.report;
        let audit = closed_account(
            report.external_heat_j,
            report.thermal_residual_j,
            report.combined_residual_j,
            candidate.observer.first.closed_radiation_net_j(),
            config.energy_tolerance_j,
        )
        .map_err(ThermoelasticHeatError::Invalid)?;
        cx.checkpoint()
            .map_err(|_| ThermoelasticHeatError::Cancelled)?;
        Ok((candidate, audit))
    }
}
