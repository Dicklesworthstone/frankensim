//! Complete source-resolved coolant assembly for original coupled drivers.
//!
//! This view retains the prepared material heat owner and the fluid resolvers;
//! it neither owns a physical state nor introduces another time integrator.

use super::{
    ConvectionHeatError, ConvectionHeatSample, SourceConvectionHeatNetwork, circuit_error,
};
use crate::graph::thermal::medium::convection::deposition::source::{
    SourceConvectiveCell, links::SourceConvectiveLink,
};
use fs_conduction::transient::deposition::DepositionSimulation;
use fs_couple::circuit_enthalpy::{
    CircuitEnthalpyError,
    contacts::{
        ContactThermalCoefficients,
        heat::{HeatAssembly, PreparedHeatNetwork},
    },
};
use fs_exec::Cx;
use fs_material::fluid_transport::ThermalFluidTransport;
use fs_qty::Temperature;

impl<R, F: ThermalFluidTransport, E> HeatAssembly
    for SourceConvectionHeatNetwork<'_, '_, '_, '_, R>
where
    R: Fn(&SourceConvectiveCell, &DepositionSimulation<'_>, Temperature) -> Result<F, E>,
{
    type Sample = ConvectionHeatSample<F>;
    type Error = ConvectionHeatError<E>;

    fn prepared_heat_network(&self) -> &PreparedHeatNetwork<'_, '_> {
        self.heat
    }

    fn explicit_source_conductance_bounds_w_k(
        &self,
        sample: &Self::Sample,
        checkpoint: &mut dyn FnMut() -> Result<(), fs_thermal::enthalpy::DepositionError>,
    ) -> Result<std::collections::BTreeMap<fs_thermal::enthalpy::CellId, f64>,
        CircuitEnthalpyError<Self::Error>> {
        sample.explicit_source_conductance_bounds_w_k(checkpoint)
            .map_err(ConvectionHeatError::<E>::from).map_err(circuit_error)
    }

    fn resolve_heat(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<ContactThermalCoefficients<Self::Sample>, CircuitEnthalpyError<Self::Error>> {
        self.resolve(cx, thermal).map_err(circuit_error)
    }
}

/// Borrowed complete heat assembly including the original finite-fluid links.
///
/// Link powers remain reciprocal native exchange. Prescribed reservoirs retain
/// their separate external heat account. The original heat sample, all fluid
/// evaluations and their source identities survive in the accepted driver report.
/// Keep the same view, source models and bindings for continued physical steps.
pub struct SourceConvectionHeatAssembly<'s, 'l, N, Q> {
    network: &'s N,
    links: &'l [SourceConvectiveLink],
    maximum_links: usize,
    resolve_links: Q,
}

impl<R> SourceConvectionHeatNetwork<'_, '_, '_, '_, R> {
    /// Bind finite-fluid convection for any original `HeatAssembly` consumer.
    ///
    /// This performs no source query or physical advance. Every consumer trial
    /// calls the existing complete resolver with that actual proposed state.
    /// Complete disjoint link ownership, source support, budget and cancellation
    /// admission therefore remain with the same original owners.
    ///
    /// Both fluid resolvers must be repeatable and retain their material/source
    /// identity across trials. Accepted physical histories must not be mutated
    /// by a callback; external callback effects cannot be rolled back.
    pub fn with_links<'s, 'l, Q, L: ThermalFluidTransport, E>(
        &'s self,
        links: &'l [SourceConvectiveLink],
        maximum_links: usize,
        resolve_links: Q,
    ) -> SourceConvectionHeatAssembly<'s, 'l, Self, Q>
    where
        Q: Fn(&SourceConvectiveLink, &DepositionSimulation<'_>, Temperature) -> Result<L, E>,
    {
        SourceConvectionHeatAssembly {
            network: self,
            links,
            maximum_links,
            resolve_links,
        }
    }
}

impl<R, Q, F: ThermalFluidTransport, L: ThermalFluidTransport, E> HeatAssembly
    for SourceConvectionHeatAssembly<'_, '_, SourceConvectionHeatNetwork<'_, '_, '_, '_, R>, Q>
where
    R: Fn(&SourceConvectiveCell, &DepositionSimulation<'_>, Temperature) -> Result<F, E>,
    Q: Fn(&SourceConvectiveLink, &DepositionSimulation<'_>, Temperature) -> Result<L, E>,
{
    type Sample = ConvectionHeatSample<F, L>;
    type Error = ConvectionHeatError<E>;

    fn prepared_heat_network(&self) -> &PreparedHeatNetwork<'_, '_> {
        self.network.heat
    }

    fn explicit_source_conductance_bounds_w_k(
        &self,
        sample: &Self::Sample,
        checkpoint: &mut dyn FnMut() -> Result<(), fs_thermal::enthalpy::DepositionError>,
    ) -> Result<std::collections::BTreeMap<fs_thermal::enthalpy::CellId, f64>,
        CircuitEnthalpyError<Self::Error>> {
        sample.explicit_source_conductance_bounds_w_k(checkpoint)
            .map_err(ConvectionHeatError::<E>::from).map_err(circuit_error)
    }

    fn resolve_heat(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<ContactThermalCoefficients<Self::Sample>, CircuitEnthalpyError<Self::Error>> {
        self.network
            .resolve_with_links(
                cx,
                thermal,
                self.links,
                self.maximum_links,
                &self.resolve_links,
            )
            .map_err(circuit_error)
    }
}
