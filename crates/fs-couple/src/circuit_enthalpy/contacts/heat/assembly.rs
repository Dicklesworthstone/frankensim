//! One complete heat evaluation shared by native coupled physical owners.
use super::{
    CircuitEnthalpyError, ContactThermalCoefficients, DepositionSimulation, HeatNetworkError,
    HeatNetworkSample, PreparedHeatNetwork, ThermalContactPressure, as_circuit_error,
};
use fs_exec::Cx;

/// A prepared heat network together with its additional physical heat laws.
///
/// The original prepared network remains authoritative for material identity,
/// source selection, resistor destinations and electrical-contact heat routes.
/// An implementation can add prescribed reservoirs or disjoint native finite
/// links, retaining their original evaluations in its associated sample.
///
/// Every call resolves the SAME proposed thermal state. Returned source powers
/// already include each added heat path once, but must exclude Joule heat and
/// reversible field heat owned by the consuming circuit/field driver. Finite
/// links transfer heat internally through the conductance map. Closed-enclosure
/// radiation stays in the original sample so its numerical net can be excluded
/// from external supply by the existing physical energy account.
///
/// Keep this provider across accepted steps to preserve source pins. Evaluation
/// must be repeatable and must not advance external physical history. A refusal
/// returns no partial coefficients; cancellation uses the native Cancelled
/// variant. This interface adds neither a solver nor a temperature/history owner.
pub trait HeatAssembly {
    type Sample: AsRef<HeatNetworkSample>;
    type Error: From<HeatNetworkError>;

    /// The exact original prepared owner, including its admission and routing.
    fn prepared_heat_network(&self) -> &PreparedHeatNetwork<'_, '_>;

    /// Admit an outer prescribed-fluid convection layer on the named cells.
    /// Composed providers must reject overlap with their existing reservoir
    /// declarations even when sampled heat flux is zero. Native boundary and
    /// exposed-area checks still belong to the actual convection evaluator.
    fn admit_additional_reservoir_convection(
        &self,
        cx: &Cx<'_>,
        cells: &std::collections::BTreeSet<fs_thermal::enthalpy::CellId>,
    ) -> Result<(), CircuitEnthalpyError<Self::Error>> {
        cx.checkpoint().map_err(|_| CircuitEnthalpyError::Cancelled)?;
        if cells.is_empty() {
            Ok(())
        } else {
            Err(CircuitEnthalpyError::Source(HeatNetworkError::Invalid(
                "heat assembly has not declared prescribed-reservoir convection ownership",
            ).into()))
        }
    }

    /// Additional explicit-source rate bounds [W/K] at the exact resolved sample.
    ///
    /// Include the original closed-radiation bound and every prescribed
    /// reservoir source. Native finite-link conductances are already handled
    /// by the thermal owner and must not appear again here. The original
    /// frozen-coefficient interpretation applies; this is not a nonlinear
    /// temporal-error estimate.
    ///
    /// Implicit coupled drivers do not require this hook. An explicit observer
    /// refuses an undeclared bound instead of silently treating an extra source
    /// as unconditionally stable. Evaluation must not advance any state.
    fn explicit_source_conductance_bounds_w_k(
        &self,
        _sample: &Self::Sample,
        _checkpoint: &mut dyn FnMut() -> Result<(), fs_thermal::enthalpy::DepositionError>,
    ) -> Result<std::collections::BTreeMap<fs_thermal::enthalpy::CellId, f64>,
        CircuitEnthalpyError<Self::Error>> {
        Err(CircuitEnthalpyError::Source(HeatNetworkError::Invalid(
            "heat assembly has not declared its explicit source-rate bounds",
        ).into()))
    }

    /// Complete transport at one original candidate state, with typed failures.
    fn resolve_heat(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<ContactThermalCoefficients<Self::Sample>, CircuitEnthalpyError<Self::Error>>;

    /// Complete heat evaluation with explicitly owned mechanical pressures on
    /// original thermal interfaces. Composed providers must forward these to
    /// their original contact resolver while retaining all additional sources
    /// and finite links exactly once. A provider that has not implemented that
    /// composition refuses a nonempty request instead of freezing pressure.
    fn resolve_heat_with_contact_pressures(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
        pressures: &[ThermalContactPressure<'_>],
    ) -> Result<ContactThermalCoefficients<Self::Sample>, CircuitEnthalpyError<Self::Error>> {
        cx.checkpoint().map_err(|_| CircuitEnthalpyError::Cancelled)?;
        if pressures.is_empty() {
            self.resolve_heat(cx, thermal)
        } else {
            Err(CircuitEnthalpyError::Source(HeatNetworkError::Invalid(
                "heat assembly does not own mechanical contact-pressure continuation",
            ).into()))
        }
    }
}

impl AsRef<HeatNetworkSample> for HeatNetworkSample {
    fn as_ref(&self) -> &HeatNetworkSample {
        self
    }
}

impl HeatAssembly for PreparedHeatNetwork<'_, '_> {
    type Sample = HeatNetworkSample;
    type Error = HeatNetworkError;

    fn prepared_heat_network(&self) -> &PreparedHeatNetwork<'_, '_> {
        self
    }

    fn admit_additional_reservoir_convection(
        &self,
        cx: &Cx<'_>,
        _cells: &std::collections::BTreeSet<fs_thermal::enthalpy::CellId>,
    ) -> Result<(), CircuitEnthalpyError<Self::Error>> {
        // This owner adds bulk/contact transport and radiation, not a fluid
        // reservoir. Actual cell boundaries and occupied surface coverage are
        // admitted by the original convection owner at the candidate state.
        cx.checkpoint().map_err(|_| CircuitEnthalpyError::Cancelled)
    }

    fn explicit_source_conductance_bounds_w_k(
        &self,
        sample: &Self::Sample,
        checkpoint: &mut dyn FnMut() -> Result<(), fs_thermal::enthalpy::DepositionError>,
    ) -> Result<std::collections::BTreeMap<fs_thermal::enthalpy::CellId, f64>,
        CircuitEnthalpyError<Self::Error>> {
        sample.explicit_source_conductance_bounds_w_k(checkpoint).map_err(as_circuit_error)
    }

    fn resolve_heat(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<ContactThermalCoefficients<Self::Sample>, CircuitEnthalpyError<Self::Error>> {
        self.resolve(cx, thermal).map_err(as_circuit_error)
    }

    fn resolve_heat_with_contact_pressures(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
        pressures: &[ThermalContactPressure<'_>],
    ) -> Result<ContactThermalCoefficients<Self::Sample>, CircuitEnthalpyError<Self::Error>> {
        PreparedHeatNetwork::resolve_heat_with_contact_pressures(self, cx, thermal, pressures)
            .map_err(as_circuit_error)
    }
}
