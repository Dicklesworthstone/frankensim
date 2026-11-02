//! Add prescribed-fluid convection without replacing a complete heat owner.
use super::{
    ConvectionHeatError, ConvectionHeatSample, ResolvedCellConvectionSources,
    SourceCellConvectionError, SourceConvectiveCell, add_convection_bounds, circuit_error,
    evaluate_source_cell_convection,
};
use super::super::links::SourceConvectiveLink;
use fs_conduction::transient::deposition::{DepositionError, DepositionSimulation};
use fs_couple::circuit_enthalpy::{
    CircuitEnthalpyError, CircuitThermalCoefficients,
    contacts::{
        ContactThermalCoefficients,
        heat::{HeatAssembly, HeatNetworkError, PreparedHeatNetwork, ThermalContactPressure,
            ThermalLinkMetric},
    },
};
use fs_exec::Cx;
use fs_material::fluid_transport::ThermalFluidTransport;
use fs_qty::Temperature;
use fs_thermal::enthalpy::CellId;
use std::collections::{BTreeMap, BTreeSet};

/// Source-resolved reservoirs on an already complete physical heat assembly.
///
/// The original assembly resolves exactly once per thermal trial. Its full
/// sample and error remain nested, together with every original coolant/contact/
/// radiation receipt. Finite links stay under that inner owner's conductance
/// map; this layer adds only prescribed-reservoir powers on exposed cell area.
/// It introduces no state, material evaluator, heat solver or time integrator.
///
/// The inner assembly must explicitly admit additional reservoir-convection
/// ownership. Existing convection on a requested cell refuses even at zero
/// flux; radiation and native finite/contact links have separate ownership.
/// Original contact coverage still masks the exposed reservoir area. A cell
/// can therefore exchange with finite coolant and a prescribed fluid on its
/// remaining surface, without inventing a second finite fluid inventory.
///
/// Retain the same prepared assembly and source-pinned fluid resolvers across
/// trials and accepted steps. The callback sees the actual candidate state and
/// must preserve material identity and support; it must not advance history.
/// Generic transport alone cannot establish that identity for the caller.
/// Rate bounds delegate to the inner owner and add this layer once; an inner
/// refusal to declare explicit bounds remains a refusal of explicit evolution.
/// Native contact pressures and internal link metrics pass to the same inner
/// owner. A metric on a reservoir cell refuses: it supplies neither transformed
/// exposed area nor the correlation's current characteristic length.
pub struct SourceConvectionHeatLayer<'p, 'b, A: ?Sized, R> {
    pub inner: &'p A,
    pub cells: &'b [SourceConvectiveCell],
    pub maximum_cells: usize,
    pub resolve_fluid: R,
}

// Shared with the original specialized convection adapters. Ownership follows
// declarations, never a sampled power or conductance that can vanish at equilibrium.
pub(super) fn admit_disjoint_reservoirs(
    cx: &Cx<'_>,
    owned: &[SourceConvectiveCell],
    additional: &BTreeSet<CellId>,
) -> Result<(), HeatNetworkError> {
    for binding in owned {
        cx.checkpoint().map_err(|_| HeatNetworkError::Cancelled)?;
        if additional.contains(&binding.cell) {
            return Err(HeatNetworkError::Invalid(
                "prescribed reservoir convection already owns this cell's exposed boundary",
            ));
        }
    }
    cx.checkpoint().map_err(|_| HeatNetworkError::Cancelled)
}

// A selected internal metric supplies no transformed exterior surface or
// characteristic length. Reference reservoirs and finite-fluid links must be
// disjoint from mapped cells, even when their present heat flux is zero.
pub(super) fn admit_reference_convection_geometry(
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'_>,
    cells: &[SourceConvectiveCell],
    links: &[SourceConvectiveLink],
    metrics: &[ThermalLinkMetric],
) -> Result<(), HeatNetworkError> {
    cx.checkpoint().map_err(|_| HeatNetworkError::Cancelled)?;
    if metrics.is_empty() {
        return Ok(());
    }
    if metrics.len() > thermal.links().len() {
        return Err(HeatNetworkError::Invalid(
            "contact metric count exceeds original thermal links",
        ));
    }
    let mut mapped = BTreeSet::new();
    for metric in metrics {
        cx.checkpoint().map_err(|_| HeatNetworkError::Cancelled)?;
        let link = thermal.links().get(&metric.link).ok_or(HeatNetworkError::Invalid(
            "contact metric names no original thermal link",
        ))?;
        mapped.insert(link.a);
        mapped.insert(link.b);
    }
    for binding in cells {
        cx.checkpoint().map_err(|_| HeatNetworkError::Cancelled)?;
        if mapped.contains(&binding.cell) {
            return Err(HeatNetworkError::Invalid(
                "contact metric cannot infer transformed reservoir boundary geometry",
            ));
        }
    }
    for binding in links {
        cx.checkpoint().map_err(|_| HeatNetworkError::Cancelled)?;
        let link = thermal.links().get(&binding.link).ok_or(HeatNetworkError::Invalid(
            "finite convection names no original thermal link",
        ))?;
        if mapped.contains(&link.a) || mapped.contains(&link.b) {
            return Err(HeatNetworkError::Invalid(
                "contact metric cannot infer transformed finite-coolant exchange geometry",
            ));
        }
    }
    cx.checkpoint().map_err(|_| HeatNetworkError::Cancelled)
}

impl<A, R, F, E> HeatAssembly for SourceConvectionHeatLayer<'_, '_, A, R>
where
    A: HeatAssembly + ?Sized,
    F: ThermalFluidTransport,
    R: Fn(&SourceConvectiveCell, &DepositionSimulation<'_>, Temperature) -> Result<F, E>,
{
    type Sample = ConvectionHeatSample<F, F, A::Sample>;
    type Error = ConvectionHeatError<E, A::Error>;

    fn prepared_heat_network(&self) -> &PreparedHeatNetwork<'_, '_> {
        self.inner.prepared_heat_network()
    }

    fn admit_additional_reservoir_convection(
        &self,
        cx: &Cx<'_>,
        cells: &BTreeSet<CellId>,
    ) -> Result<(), CircuitEnthalpyError<Self::Error>> {
        admit_disjoint_reservoirs(cx, self.cells, cells)
            .map_err(Self::Error::from)
            .map_err(circuit_error)?;
        self.inner
            .admit_additional_reservoir_convection(cx, cells)
            .map_err(|error| error.map_source(ConvectionHeatError::Heat))
    }

    fn explicit_source_conductance_bounds_w_k(
        &self,
        sample: &Self::Sample,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, CircuitEnthalpyError<Self::Error>> {
        let inner = self
            .inner
            .explicit_source_conductance_bounds_w_k(&sample.heat, checkpoint)
            .map_err(|error| error.map_source(ConvectionHeatError::Heat))?;
        add_convection_bounds(inner, &sample.convection, checkpoint)
            .map_err(Self::Error::from)
            .map_err(circuit_error)
    }

    fn resolve_heat(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<ContactThermalCoefficients<Self::Sample>, CircuitEnthalpyError<Self::Error>> {
        self.resolve_heat_with_contact_pressures(cx, thermal, &[])
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
        if self.cells.len() > self.maximum_cells || self.cells.len() > thermal.cells().len() {
            return Err(CircuitEnthalpyError::Source(
                HeatNetworkError::Invalid("convective cells exceed their explicit work bound")
                    .into(),
            ));
        }
        let mut ids = BTreeSet::new();
        for binding in self.cells {
            cx.checkpoint()
                .map_err(|_| CircuitEnthalpyError::Cancelled)?;
            if !ids.insert(binding.cell) {
                return Err(CircuitEnthalpyError::Source(
                    HeatNetworkError::Invalid("duplicate prescribed-reservoir convection cell")
                        .into(),
                ));
            }
        }
        admit_reference_convection_geometry(cx, thermal, self.cells, &[], metrics)
            .map_err(Self::Error::from)
            .map_err(circuit_error)?;
        self.inner
            .admit_additional_reservoir_convection(cx, &ids)
            .map_err(|error| error.map_source(ConvectionHeatError::Heat))?;
        let mut coefficients = self
            .inner
            .resolve_heat_with_contact_metrics(cx, thermal, pressures, metrics)
            .map_err(|error| error.map_source(ConvectionHeatError::Heat))?;
        let convection: ResolvedCellConvectionSources<F> = evaluate_source_cell_convection(
            thermal,
            self.cells,
            self.maximum_cells,
            |binding, temperature| (self.resolve_fluid)(binding, thermal, temperature),
            || cx.checkpoint().is_err(),
        )
        .map_err(Self::Error::from)
        .map_err(circuit_error)?;
        convection
            .convection()
            .add_to(&mut coefficients.bulk.sources_w, || {
                cx.checkpoint().is_err()
            })
            .map_err(SourceCellConvectionError::<E>::from)
            .map_err(Self::Error::from)
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
                source: ConvectionHeatSample {
                    heat: bulk.source,
                    convection,
                    finite_convection: None,
                },
            },
        })
    }
}
