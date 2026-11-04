//! Two-wall gas layers in the original complete material heat assembly.

use super::super::heat::ConvectionHeatError;
use super::{
    ResolvedEnclosedGasLayers, SourceEnclosedGasLayer, evaluate_source_enclosed_gas_layers,
};
use fs_couple::circuit_enthalpy::{
    CircuitEnthalpyError, CircuitThermalCoefficients,
    contacts::{
        ContactThermalCoefficients,
        heat::{
            HeatAssembly, HeatNetworkSample, PreparedHeatNetwork, ThermalContactPressure,
            ThermalLinkMetric,
        },
    },
};
use fs_exec::Cx;
use fs_material::fluid_transport::ThermalFluidTransport;
use fs_qty::Temperature;
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use std::collections::{BTreeMap, BTreeSet};

/// All original solid/electrical/radiation and gas-gap sources at one trial.
#[derive(Debug)]
pub struct EnclosedGasHeatSample<F> {
    pub heat: HeatNetworkSample,
    /// Total gas conductances already in the native reciprocal link map.
    /// Their powers are internal exchange, never a second heat-source map.
    pub gas_layers: ResolvedEnclosedGasLayers<F>,
}

impl<F> AsRef<HeatNetworkSample> for EnclosedGasHeatSample<F> {
    fn as_ref(&self) -> &HeatNetworkSample {
        &self.heat
    }
}

/// Complete heat assembly with quasi-steady horizontal or inclined gas gaps.
///
/// Both wall caloric states, material resistance, contact laws and enclosure
/// radiation are evaluated at the SAME proposed thermal endpoint. Gas storage
/// is negligible; the original link transfers its total gas heat reciprocally.
/// No gas temperature field, source pair or additional time integrator exists.
///
/// First evaluate `layers` at the accepted state with
/// [`evaluate_source_enclosed_gas_layers`], then prepare the original
/// `MaterialHeatNetwork` using `prepare_with_link_conductances` and that result's
/// conductance map. Fourier bindings must exclude every gas-gap link; an added
/// material contact resistance on the same link also refuses. Retain that
/// prepared owner here. Preparation pins material laws, not gas coefficients:
/// every subsequent trial re-evaluates both actual wall temperatures and the
/// same gas resolver. Missing or competing link owners remain native failures.
///
/// Pass this assembly to the existing `HeatAssembly` circuit/field drivers or
/// `observe_heat_assembly`. Wrap it in `SourceConvectionHeatLayer` to add outside
/// prescribed-fluid reservoirs on the remaining exposed wall area. Native gap
/// conductances enter the existing thermal timestep bound once; external-source
/// bounds retain their original radiation/reservoir owner.
///
/// Mechanical contact pressures and disjoint current contact metrics pass to
/// the prepared network. Metrics affecting either gap wall refuse because they
/// provide no current gap, plate spans or wall orientation. Rebind explicitly
/// after a geometry change. All original air, aspect-ratio, density-departure,
/// uniform-wall and negligible-gas-storage requirements remain in force.
///
/// Keep immutable gas laws, composition and the declared pressure model in
/// `resolve_gas` across trials and accepted steps. The callback must be repeatable
/// and preserve source support; it must not advance external physical history.
pub struct SourceEnclosedGasHeatAssembly<'p, 'n, 'a, 'b, R> {
    pub heat: &'p PreparedHeatNetwork<'n, 'a>,
    pub layers: &'b [SourceEnclosedGasLayer],
    pub maximum_layers: usize,
    pub resolve_gas: R,
}

impl<R, F: ThermalFluidTransport, E> HeatAssembly
    for SourceEnclosedGasHeatAssembly<'_, '_, '_, '_, R>
where
    R: Fn(&SourceEnclosedGasLayer, &DepositionSimulation<'_>, Temperature) -> Result<F, E>,
{
    type Sample = EnclosedGasHeatSample<F>;
    type Error = ConvectionHeatError<E>;

    fn prepared_heat_network(&self) -> &PreparedHeatNetwork<'_, '_> {
        self.heat
    }

    fn admit_additional_reservoir_convection(
        &self,
        cx: &Cx<'_>,
        cells: &BTreeSet<CellId>,
    ) -> Result<(), CircuitEnthalpyError<Self::Error>> {
        self.heat
            .admit_additional_reservoir_convection(cx, cells)
            .map_err(|error| error.map_source(ConvectionHeatError::from))
    }

    fn explicit_source_conductance_bounds_w_k(
        &self,
        sample: &Self::Sample,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, CircuitEnthalpyError<Self::Error>> {
        // The native link already owns the gas G and its CFL contribution.
        self.heat
            .explicit_source_conductance_bounds_w_k(&sample.heat, checkpoint)
            .map_err(|error| error.map_source(ConvectionHeatError::from))
    }

    fn resolve_heat(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<ContactThermalCoefficients<Self::Sample>, CircuitEnthalpyError<Self::Error>> {
        self.resolve_heat_with_contact_metrics(cx, thermal, &[], &[])
    }

    fn resolve_heat_with_contact_pressures(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
        pressures: &[ThermalContactPressure<'_>],
    ) -> Result<ContactThermalCoefficients<Self::Sample>, CircuitEnthalpyError<Self::Error>> {
        self.resolve_heat_with_contact_metrics(cx, thermal, pressures, &[])
    }

    fn resolve_heat_with_contact_metrics(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
        pressures: &[ThermalContactPressure<'_>],
        metrics: &[ThermalLinkMetric],
    ) -> Result<ContactThermalCoefficients<Self::Sample>, CircuitEnthalpyError<Self::Error>> {
        cx.checkpoint()
            .map_err(|_| CircuitEnthalpyError::Cancelled)?;
        let gas_layers = evaluate_source_enclosed_gas_layers(
            thermal,
            self.layers,
            self.maximum_layers,
            &self.resolve_gas,
            || cx.checkpoint().is_err(),
        )
        .map_err(ConvectionHeatError::from)
        .map_err(circuit_error)?;
        let coefficients = self
            .heat
            .resolve_with_link_conductances_and_contact_metrics(
                cx,
                thermal,
                gas_layers.conductances_w_k(),
                pressures,
                metrics,
            )
            .map_err(ConvectionHeatError::from)
            .map_err(circuit_error)?;
        cx.checkpoint()
            .map_err(|_| CircuitEnthalpyError::Cancelled)?;
        let bulk = coefficients.bulk;
        Ok(ContactThermalCoefficients {
            contacts: coefficients.contacts,
            bulk: CircuitThermalCoefficients {
                resistors: bulk.resistors,
                sources_w: bulk.sources_w,
                conductances_w_k: bulk.conductances_w_k,
                source: EnclosedGasHeatSample {
                    heat: bulk.source,
                    gas_layers,
                },
            },
        })
    }
}

fn circuit_error<E>(error: ConvectionHeatError<E>) -> CircuitEnthalpyError<ConvectionHeatError<E>> {
    match error {
        ConvectionHeatError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(other),
    }
}
