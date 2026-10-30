//! Full material heat assembly and finite/prescribed convection on native steps.
use super::*;
use fs_thermal::enthalpy::{CellId, DepositionObserver, Event};
use fs_thermal::observers::CoupledObservers;
use std::collections::BTreeMap;
use std::sync::Arc;

// One immutable starting-state evaluation serves the native rate hooks.
// It is discarded at acceptance, events, or the start of a new advance.
struct SourceEpoch<F, L> {
    time_s: f64,
    sample: Arc<ConvectionHeatSample<F, L>>,
    conductances_w_k: BTreeMap<u64, f64>,
    additional_conductances_w_k: BTreeMap<CellId, f64>,
}

/// Complete thermal-only assembly on the original explicit enthalpy timeline.
///
/// The native thermal cells retain mass, caloric laws and temperature history.
/// Their links apply finite-fluid exchange and its CFL contribution exactly
/// once. Prescribed reservoirs and enclosure radiation share the same material
/// source epoch; their conductance bounds further tighten the native step.
///
/// Callback clones must deeply copy any mutable physical history. Immutable
/// prepared material laws may be shared. External effects or shared mutable
/// state cannot be rolled back by the thermal timeline.
pub struct SourceConvectionHeatObserver<'p, 'n, 'a, 'b, 'scope, R, Q, F, L> {
    network: SourceConvectionHeatNetwork<'p, 'n, 'a, 'b, R>,
    cx: &'p Cx<'scope>,
    links: &'p [SourceConvectiveLink],
    maximum_links: usize,
    resolve_links: Q,
    pending: Option<Arc<SourceEpoch<F, L>>>,
    applied: Option<Arc<ConvectionHeatSample<F, L>>>,
    endpoint: Option<Arc<ConvectionHeatSample<F, L>>>,
    external_convection_heat_j: f64,
    closed_radiation_net_j: f64,
    maximum_closed_radiation_defect_j: f64,
}
impl<R: Clone, Q: Clone, F, L> Clone
    for SourceConvectionHeatObserver<'_, '_, '_, '_, '_, R, Q, F, L>
{
    fn clone(&self) -> Self {
        Self {
            network: SourceConvectionHeatNetwork {
                heat: self.network.heat,
                cells: self.network.cells,
                maximum_cells: self.network.maximum_cells,
                resolve_fluid: self.network.resolve_fluid.clone(),
            },
            cx: self.cx,
            links: self.links,
            maximum_links: self.maximum_links,
            resolve_links: self.resolve_links.clone(),
            pending: self.pending.clone(),
            applied: self.applied.clone(),
            endpoint: self.endpoint.clone(),
            external_convection_heat_j: self.external_convection_heat_j,
            closed_radiation_net_j: self.closed_radiation_net_j,
            maximum_closed_radiation_defect_j: self.maximum_closed_radiation_defect_j,
        }
    }
}
impl<R, Q, F, L> SourceConvectionHeatObserver<'_, '_, '_, '_, '_, R, Q, F, L> {
    /// Original source evaluation actually applied at the START of the latest
    /// accepted substep. Events do not rewrite that historical observation.
    pub fn applied_sample(&self) -> Option<&ConvectionHeatSample<F, L>> {
        self.applied.as_deref()
    }
    /// Fresh material sources at the latest accepted endpoint. Events
    /// invalidate this until the next substep admits the completed active set.
    pub fn sample(&self) -> Option<&ConvectionHeatSample<F, L>> {
        self.endpoint.as_deref()
    }
    /// Signed prescribed-convection heat in the latest candidate advance.
    /// It is already included once in the native external-heat audit.
    pub const fn external_convection_heat_j(&self) -> f64 {
        self.external_convection_heat_j
    }
    /// Signed numerical net of closed enclosure exchange in the latest advance.
    /// Subtract this once from the thermal audit when interpreting external
    /// heat. The configured gate bounds it before any endpoint is accepted.
    pub const fn closed_radiation_net_j(&self) -> f64 {
        self.closed_radiation_net_j
    }
}

impl<'n, 'a, 'b, R> SourceConvectionHeatNetwork<'_, 'n, 'a, 'b, R> {
    /// Compose this whole heat assembly with an existing physical observer.
    ///
    /// The first returned owner applies material conduction/contact/layers,
    /// closed radiation and both convection kinds at the same native substeps;
    /// the second keeps the caller's history. Electrical resistor/contact/basis
    /// bindings refuse: use the coupled circuit owner when electrical work is
    /// present. Competing native link owners refuse through CoupledObservers.
    ///
    /// Prepare the heat network with the same finite-link conductance partition
    /// used by `resolve_with_links`. Empty links permit reservoir-only exchange.
    /// Caloric/source laws, geometry, area coverage and relative flow retain the
    /// original adapter's domain assumptions. No new solver or state is added.
    ///
    /// The explicit source coefficients and emissivities are frozen per step.
    /// Combined conductance bounds do not establish nonlinear stability or
    /// temporal accuracy; refine max_dt_s. Source failures, work exhaustion and
    /// cancellation discard the entire thermal/observer candidate together.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn observe_with_links<'p, 'scope, Q, F, L, E, O>(
        &'p self,
        cx: &'p Cx<'scope>,
        links: &'p [SourceConvectiveLink],
        maximum_links: usize,
        resolve_links: Q,
        inner: O,
        maximum_closed_radiation_defect_j: f64,
    ) -> Result<
        CoupledObservers<SourceConvectionHeatObserver<'p, 'n, 'a, 'b, 'scope, R, Q, F, L>, O>,
        ConvectionHeatError<E>,
    >
    where
        R: Clone
            + Fn(&SourceConvectiveCell, &DepositionSimulation<'_>, Temperature) -> Result<F, E>,
        Q: Clone
            + Fn(&SourceConvectiveLink, &DepositionSimulation<'_>, Temperature) -> Result<L, E>,
        F: ThermalFluidTransport,
        L: ThermalFluidTransport,
    {
        cx.checkpoint()
            .map_err(|_| ConvectionHeatError::Cancelled)?;
        self.heat.admit_thermal_only()?;
        if !maximum_closed_radiation_defect_j.is_finite() || maximum_closed_radiation_defect_j < 0.0
        {
            return Err(HeatNetworkError::Invalid(
                "finite nonnegative radiation energy allowance required",
            )
            .into());
        }
        Ok(CoupledObservers::new(
            SourceConvectionHeatObserver {
                network: SourceConvectionHeatNetwork {
                    heat: self.heat,
                    cells: self.cells,
                    maximum_cells: self.maximum_cells,
                    resolve_fluid: self.resolve_fluid.clone(),
                },
                cx,
                links,
                maximum_links,
                resolve_links,
                pending: None,
                applied: None,
                endpoint: None,
                external_convection_heat_j: 0.0,
                closed_radiation_net_j: 0.0,
                maximum_closed_radiation_defect_j,
            },
            inner,
        ))
    }
}

impl<R, Q, F, L, E> SourceConvectionHeatObserver<'_, '_, '_, '_, '_, R, Q, F, L>
where
    R: Fn(&SourceConvectiveCell, &DepositionSimulation<'_>, Temperature) -> Result<F, E>,
    Q: Fn(&SourceConvectiveLink, &DepositionSimulation<'_>, Temperature) -> Result<L, E>,
    F: ThermalFluidTransport,
    L: ThermalFluidTransport,
{
    fn resolve_checked(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ContactThermalCoefficients<ConvectionHeatSample<F, L>>, ConvectionHeatError<E>>
    {
        // The original evaluators poll Cx inside bounded work. Keep the
        // timeline's richer checkpoint failures before/after that evaluation.
        checkpoint()?;
        let result = self.network.resolve_with_links(
            self.cx,
            thermal,
            self.links,
            self.maximum_links,
            &self.resolve_links,
        )?;
        checkpoint()?;
        Ok(result)
    }

    fn pending(&self, time_s: f64) -> Result<&SourceEpoch<F, L>, ConvectionHeatError<E>> {
        self.pending
            .as_deref()
            .filter(|epoch| epoch.time_s == time_s)
            .ok_or_else(|| {
                HeatNetworkError::Invalid("complete heat transfer lacks its current source epoch")
                    .into()
            })
    }
}

impl<R, Q, F, L, E> DepositionObserver
    for SourceConvectionHeatObserver<'_, '_, '_, '_, '_, R, Q, F, L>
where
    R: Clone + Fn(&SourceConvectiveCell, &DepositionSimulation<'_>, Temperature) -> Result<F, E>,
    Q: Clone + Fn(&SourceConvectiveLink, &DepositionSimulation<'_>, Temperature) -> Result<L, E>,
    F: ThermalFluidTransport,
    L: ThermalFluidTransport,
{
    type Error = ConvectionHeatError<E>;

    fn begin(&mut self) -> Result<(), Self::Error> {
        self.cx
            .checkpoint()
            .map_err(|_| ConvectionHeatError::Cancelled)?;
        self.pending = None;
        self.external_convection_heat_j = 0.0;
        self.closed_radiation_net_j = 0.0;
        Ok(())
    }

    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        let coefficients = self.resolve_checked(thermal, checkpoint)?.bulk;
        let source = coefficients.source;
        let additional = source.explicit_source_conductance_bounds_w_k(&mut *checkpoint)?;
        if !source.heat.radiation_net_power_w().is_finite()
            || !source
                .convection
                .convection()
                .total_heat_into_cells_w()
                .is_finite()
        {
            return Err(HeatNetworkError::Invalid("nonfinite complete heat source balance").into());
        }
        checkpoint()?;
        self.pending = Some(Arc::new(SourceEpoch {
            time_s: thermal.time_s(),
            sample: Arc::new(source),
            conductances_w_k: coefficients.conductances_w_k,
            additional_conductances_w_k: additional,
        }));
        Ok(coefficients.sources_w)
    }

    fn conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<u64, f64>, Self::Error> {
        checkpoint()?;
        // Native sources/conductances hooks see the same unchanged start state.
        // Return the exact G belonging to the already retained source samples.
        Ok(self.pending(thermal.time_s())?.conductances_w_k.clone())
    }

    fn limit_timestep_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        proposed_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<f64, Self::Error> {
        checkpoint()?;
        let epoch = self.pending(thermal.time_s())?;
        Ok(thermal.limit_timestep_for_additional_conductances(
            proposed_s,
            &epoch.additional_conductances_w_k,
            checkpoint,
        )?)
    }

    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        self.cx
            .checkpoint()
            .map_err(|_| ConvectionHeatError::Cancelled)?;
        // Same-time births arrive individually; admit the completed active set
        // on its next native source epoch, not midway through a batch.
        self.pending = None;
        self.endpoint = None;
        Ok(())
    }

    fn step(
        &mut self,
        start_s: f64,
        duration_s: f64,
        _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        checkpoint()?;
        let epoch = self.pending(start_s)?;
        if !duration_s.is_finite() || duration_s <= 0.0 {
            return Err(HeatNetworkError::Invalid(
                "complete heat accounting needs the actual source interval",
            )
            .into());
        }
        let external = self.external_convection_heat_j
            + duration_s
                * epoch
                    .sample
                    .convection
                    .convection()
                    .total_heat_into_cells_w();
        let radiation =
            self.closed_radiation_net_j + duration_s * epoch.sample.heat.radiation_net_power_w();
        if !external.is_finite() {
            return Err(
                HeatNetworkError::Invalid("integrated convective heat is unrepresentable").into(),
            );
        }
        if !radiation.is_finite() || radiation.abs() > self.maximum_closed_radiation_defect_j {
            return Err(ConvectionHeatError::ClosedEnergy {
                residual_j: radiation,
                tolerance_j: self.maximum_closed_radiation_defect_j,
            });
        }
        let endpoint = self.resolve_checked(thermal, checkpoint)?.bulk.source;
        checkpoint()?;
        self.applied = self.pending.take().map(|epoch| epoch.sample.clone());
        self.endpoint = Some(Arc::new(endpoint));
        self.external_convection_heat_j = external;
        self.closed_radiation_net_j = radiation;
        Ok(())
    }
}
