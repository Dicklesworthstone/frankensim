//! Complete native heat assembly with source-resolved conducting-wall films.
//!
//! The original prepared heat network owns every other heat mechanism. Generic
//! HeatAssembly consumers retain their existing state, residual and time stepper.

use super::super::super::ConvectionSourceError;
use super::super::{
    SourceCellConvectionError,
    links::{
        SourceConvectiveLink,
        conjugate::{
            PreparedSourceConjugateLinks, ResolvedSourceConjugateHeat, SourceConjugateError,
        },
    },
};
use super::composition::admit_reference_convection_geometry;
use fs_couple::circuit_enthalpy::{
    CircuitEnthalpyError, CircuitThermalCoefficients,
    contacts::{
        ContactThermalCoefficients,
        heat::{
            HeatAssembly, HeatNetworkError, HeatNetworkSample, PreparedHeatNetwork,
            ThermalContactPressure, ThermalLinkMetric,
        },
    },
};
use fs_exec::Cx;
use fs_material::fluid_transport::ThermalFluidTransport;
use fs_material::thermal::conduction::{LayerConductionConfig, source::LayerSourceError};
use fs_qty::Temperature;
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug)]
pub enum SourceConjugateHeatError<E> {
    Heat(HeatNetworkError),
    Conjugate(SourceConjugateError<E>),
    Cancelled,
}
impl<E> From<HeatNetworkError> for SourceConjugateHeatError<E> {
    fn from(error: HeatNetworkError) -> Self {
        match error {
            HeatNetworkError::Cancelled => Self::Cancelled,
            other => Self::Heat(other),
        }
    }
}
impl<E> From<SourceConjugateError<E>> for SourceConjugateHeatError<E> {
    fn from(error: SourceConjugateError<E>) -> Self {
        match error {
            SourceConjugateError::Thermal(DepositionError::Cancelled)
            | SourceConjugateError::Layer(LayerSourceError::Thermal(DepositionError::Cancelled))
            | SourceConjugateError::Convection(SourceCellConvectionError::Convection(
                ConvectionSourceError::Cancelled,
            )) => Self::Cancelled,
            other => Self::Conjugate(other),
        }
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for SourceConjugateHeatError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source conjugate heat assembly: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for SourceConjugateHeatError<E> {}

/// Full original heat-network sample and each actual material/film balance.
/// The combined finite-link conductances are internal exchange, never an
/// external heat supply; no bare-film coefficient is separately applied.
#[derive(Debug)]
pub struct SourceConjugateHeatSample<F> {
    pub heat: HeatNetworkSample,
    pub conjugate: ResolvedSourceConjugateHeat<F>,
}
impl<F> AsRef<HeatNetworkSample> for SourceConjugateHeatSample<F> {
    fn as_ref(&self) -> &HeatNetworkSample {
        &self.heat
    }
}

/// A borrowed prepared heat network and its finite wall/fluid material paths.
///
/// Prepare `heat` with the conjugate owner's initial combined conductance map
/// through `prepare_with_link_conductances`; its other Fourier/contact owners
/// must exclude these same links. Every native trial re-evaluates surface T,
/// original fluid properties and the complete heat network on the SAME state.
/// The original reciprocal links own both energy debits and explicit rate
/// bounds. Native circuit/field/mechanical drivers keep their usual energy gates.
///
/// Prescribed reservoirs may compose around this provider through the existing
/// `SourceConvectionHeatLayer`; original contact coverage excludes these faces.
/// Disjoint mechanical contact pressures and current path metrics pass through
/// to the original heat owner. A metric touching a conjugate endpoint refuses
/// rather than silently retaining obsolete wall/film geometry.
///
/// Keep the same source-pinned material layers and repeatable fluid resolver
/// across trials. Callback evaluation must not advance physical history.
pub struct SourceConjugateHeatAssembly<'p, 'n, 'a, 's, 'l, 'm, R> {
    pub heat: &'p PreparedHeatNetwork<'n, 'a>,
    pub links: &'s PreparedSourceConjugateLinks<'l, 'm>,
    pub layer_config: LayerConductionConfig,
    pub resolve_fluid: R,
}

impl<R, F: ThermalFluidTransport, E> HeatAssembly
    for SourceConjugateHeatAssembly<'_, '_, '_, '_, '_, '_, R>
where
    R: Fn(&SourceConvectiveLink, &DepositionSimulation<'_>, Temperature) -> Result<F, E>,
{
    type Sample = SourceConjugateHeatSample<F>;
    type Error = SourceConjugateHeatError<E>;

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
            .map_err(|error| error.map_source(SourceConjugateHeatError::from))
    }

    fn explicit_source_conductance_bounds_w_k(
        &self,
        sample: &Self::Sample,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, CircuitEnthalpyError<Self::Error>> {
        // Finite conjugate G is already in the native conductance matrix.
        // Only the original source/radiation rate allowance belongs here.
        sample
            .heat
            .explicit_source_conductance_bounds_w_k(checkpoint)
            .map_err(SourceConjugateHeatError::from)
            .map_err(circuit_error)
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
        admit_reference_convection_geometry(
            cx,
            thermal,
            &[],
            self.links.convection_bindings(),
            metrics,
        )
        .map_err(SourceConjugateHeatError::from)
        .map_err(circuit_error)?;
        let conjugate = self
            .links
            .resolve(
                thermal,
                self.layer_config,
                |binding, state, temperature| (self.resolve_fluid)(binding, state, temperature),
                || cx.checkpoint().map_err(|_| DepositionError::Cancelled),
            )
            .map_err(SourceConjugateHeatError::from)
            .map_err(circuit_error)?;
        let coefficients = self
            .heat
            .resolve_with_link_conductances_and_contact_metrics(
                cx,
                thermal,
                conjugate.conductances_w_k(),
                pressures,
                metrics,
            )
            .map_err(SourceConjugateHeatError::from)
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
                source: SourceConjugateHeatSample {
                    heat: bulk.source,
                    conjugate,
                },
            },
        })
    }
}

fn circuit_error<E>(
    error: SourceConjugateHeatError<E>,
) -> CircuitEnthalpyError<SourceConjugateHeatError<E>> {
    match error {
        SourceConjugateHeatError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(other),
    }
}
