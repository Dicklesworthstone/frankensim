//! Complete sourced heat transfer on the existing thermal observer timeline.
//! Prepared networks and complete assemblies share this one observer and state.
use std::collections::BTreeMap;
use std::sync::Arc;
use fs_thermal::enthalpy::{CellId, DepositionObserver, Event};
use fs_thermal::observers::CoupledObservers;
use super::*;

/// Original prepared-network observer; its source and error types are unchanged.
pub type MaterialHeatObserver<'p, 'n, 'a, 'scope> =
    HeatSourceObserver<'p, 'scope, PreparedHeatNetwork<'n, 'a>, HeatNetworkSample, HeatNetworkError>;

/// Complete-assembly observer, retaining the provider's original sample/error.
pub type AssemblyHeatObserver<'p, 'scope, A> =
    HeatSourceObserver<'p, 'scope, A, <A as HeatAssembly>::Sample,
        CircuitEnthalpyError<<A as HeatAssembly>::Error>>;

/// Shared thermal observer for the two source-preserving entry points above.
///
/// Construction stays with `PreparedHeatNetwork::observe` or
/// `observe_heat_assembly`. Immutable laws are re-evaluated at actual thermal
/// states; no independent temperature field or future coefficient cache exists.
/// Cloning retains only immutable source references and original observer
/// history. Neither the provider nor its sample/error needs to implement Clone.
pub struct HeatSourceObserver<'p, 'scope, N: ?Sized, S, E> {
    source: &'p N,
    resolve: fn(&N, &Cx<'_>, &DepositionSimulation<'_>)
        -> Result<ContactThermalCoefficients<S>, E>,
    bounds: fn(&N, &S, &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<BTreeMap<CellId, f64>, E>,
    error: fn(HeatNetworkError) -> E,
    cx: &'p Cx<'scope>,
    sample: Option<Arc<S>>,
    source_epoch: Option<(f64, f64)>,
    source_conductances_w_k: BTreeMap<CellId, f64>,
    closed_radiation_net_j: f64,
    maximum_closed_radiation_defect_j: f64,
}

impl<N: ?Sized, S, E> Clone for HeatSourceObserver<'_, '_, N, S, E> {
    fn clone(&self) -> Self {
        Self {
            source: self.source, resolve: self.resolve, bounds: self.bounds, error: self.error,
            cx: self.cx, sample: self.sample.clone(), source_epoch: self.source_epoch,
            source_conductances_w_k: self.source_conductances_w_k.clone(),
            closed_radiation_net_j: self.closed_radiation_net_j,
            maximum_closed_radiation_defect_j: self.maximum_closed_radiation_defect_j,
        }
    }
}

impl<N: ?Sized, S: AsRef<HeatNetworkSample>, E> HeatSourceObserver<'_, '_, N, S, E> {
    /// Complete source results at the last accepted endpoint/event. Rates used
    /// to advance it were evaluated at its substep START, as required by the
    /// original explicit thermal integrator.
    pub fn sample(&self) -> Option<&S> { self.sample.as_deref() }

    /// Integrated internal radiative roundoff in the most recent advance.
    /// Subtract this ONCE from the thermal audit's external source account.
    pub const fn closed_radiation_net_j(&self) -> f64 { self.closed_radiation_net_j }

    fn resolve_checked(
        &self, thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ContactThermalCoefficients<S>, E> {
        checkpoint().map_err(HeatNetworkError::from).map_err(self.error)?;
        let result = (self.resolve)(self.source, self.cx, thermal)?;
        checkpoint().map_err(HeatNetworkError::from).map_err(self.error)?;
        Ok(result)
    }
}

#[allow(clippy::too_many_arguments, clippy::type_complexity)]
fn observer<'p, 'scope, N: ?Sized, S, E>(
    source: &'p N,
    prepared: &PreparedHeatNetwork<'_, '_>,
    resolve: fn(&N, &Cx<'_>, &DepositionSimulation<'_>)
        -> Result<ContactThermalCoefficients<S>, E>,
    bounds: fn(&N, &S, &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<BTreeMap<CellId, f64>, E>,
    error: fn(HeatNetworkError) -> E,
    cx: &'p Cx<'scope>,
    maximum_closed_radiation_defect_j: f64,
) -> Result<HeatSourceObserver<'p, 'scope, N, S, E>, HeatNetworkError> {
    check(cx)?;
    if !maximum_closed_radiation_defect_j.is_finite()
        || maximum_closed_radiation_defect_j < 0.0 {
        return Err(HeatNetworkError::Invalid("finite nonnegative radiation energy allowance required"));
    }
    prepared.admit_thermal_only()?;
    Ok(HeatSourceObserver {
        source, resolve, bounds, error, cx, sample: None, source_epoch: None,
        source_conductances_w_k: BTreeMap::new(), closed_radiation_net_j: 0.0,
        maximum_closed_radiation_defect_j,
    })
}

fn assembly_error<E: From<HeatNetworkError>>(error: HeatNetworkError) -> CircuitEnthalpyError<E> {
    match error {
        HeatNetworkError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(other.into()),
    }
}

/// Compose one COMPLETE material heat assembly with the caller's independent
/// observer on the original finite thermal timeline.
///
/// Prescribed coolant reservoirs and disjoint finite coolant links retain their
/// original source evaluations in `sample()`. Their powers/conductances are
/// supplied once, together with the prepared material network. No extra thermal
/// step or material resampling at a different physical state is introduced.
///
/// Native thermal timestep bounds include all finite-link conductances,
/// prescribed-reservoir bounds and the existing frozen-emissivity radiation
/// bound. This is first-order explicit
/// transport, requiring physical/thermal refinement; source support alone is
/// not a temporal-error estimate. Electrical bindings refuse without a current
/// owner. Competing conductance owners refuse through CoupledObservers.
///
/// Retain the assembly to preserve pins. On continuation retain the returned
/// `second` observer for independent history; each new first observer starts
/// interval accounting afresh. Late source/checkpoint errors publish no state.
#[allow(clippy::type_complexity)]
pub fn observe_heat_assembly<'p, 'scope, A: HeatAssembly + ?Sized, O>(
    assembly: &'p A,
    cx: &'p Cx<'scope>,
    inner: O,
    maximum_closed_radiation_defect_j: f64,
) -> Result<CoupledObservers<AssemblyHeatObserver<'p, 'scope, A>, O>, HeatNetworkError> {
    Ok(CoupledObservers::new(observer(
        assembly, assembly.prepared_heat_network(),
        |source, cx, thermal| source.resolve_heat(cx, thermal),
        |source, sample, checkpoint| source.explicit_source_conductance_bounds_w_k(sample, checkpoint),
        assembly_error::<A::Error>, cx, maximum_closed_radiation_defect_j,
    )?, inner))
}

/// Preserve the original binding/admission failure or native coupled failure.
#[derive(Debug)]
pub enum HeatAssemblyAdvanceError<B, C> {
    Binding(B),
    Coupled(C),
}
impl<B: core::fmt::Debug, C: core::fmt::Debug> core::fmt::Display
    for HeatAssemblyAdvanceError<B, C> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "complete material heat assembly: {self:?}")
    }
}
impl<B: std::error::Error + 'static, C: std::error::Error + 'static> std::error::Error
    for HeatAssemblyAdvanceError<B, C> {}

impl<'n, 'a> PreparedHeatNetwork<'n, 'a> {
    /// Compose sourced conduction, thermal-interface resistance and radiation
    /// with another observer on the SAME accepted thermal substeps.
    ///
    /// The original prepared-network sample and typed failures are preserved.
    /// Complete `HeatAssembly` providers use `observe_heat_assembly`; both
    /// entry points share this exact timeline, radiation account and source
    /// epoch admission. Neither adds a circuit or a parallel time integrator.
    ///
    /// The first returned observer retains source results and the bounded net
    /// radiation defect; the second retains the caller's original history.
    /// Competing conductance owners refuse through CoupledObservers. Topology
    /// changes require explicit rebinding; refine the native thermal timestep
    /// for temperature-dependent transport and radiation.
    pub fn observe<'p, 'scope, O>(
        &'p self, cx: &'p Cx<'scope>, inner: O,
        maximum_closed_radiation_defect_j: f64,
    ) -> Result<CoupledObservers<MaterialHeatObserver<'p, 'n, 'a, 'scope>, O>, HeatNetworkError> {
        Ok(CoupledObservers::new(observer(
            self, self, |source, cx, thermal| source.resolve(cx, thermal),
            |_, sample, checkpoint| sample.explicit_source_conductance_bounds_w_k(checkpoint),
            core::convert::identity, cx, maximum_closed_radiation_defect_j,
        )?, inner))
    }
}

impl<N: ?Sized, S: AsRef<HeatNetworkSample>, E: From<DepositionError>> DepositionObserver
    for HeatSourceObserver<'_, '_, N, S, E> {
    type Error = E;

    fn begin(&mut self) -> Result<(), Self::Error> {
        check(self.cx).map_err(self.error)?;
        self.source_epoch = None;
        self.source_conductances_w_k.clear();
        self.closed_radiation_net_j = 0.0;
        Ok(())
    }

    fn sources(&mut self, thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<BTreeMap<CellId, f64>, Self::Error> {
        let evaluated = self.resolve_checked(thermal, checkpoint)?;
        let original = evaluated.bulk.source.as_ref();
        let net = original.radiation_net_power_w();
        if !net.is_finite() {
            return Err((self.error)(HeatNetworkError::Invalid("nonfinite internal radiation balance")));
        }
        let bounds = (self.bounds)(self.source, &evaluated.bulk.source, checkpoint)?;
        self.source_epoch = Some((thermal.time_s(), net));
        self.source_conductances_w_k = bounds;
        Ok(evaluated.bulk.sources_w)
    }

    fn conductances(&mut self, thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<BTreeMap<u64, f64>, Self::Error> {
        Ok(self.resolve_checked(thermal, checkpoint)?.bulk.conductances_w_k)
    }

    fn limit_timestep_s(
        &mut self, thermal: &DepositionSimulation<'_>, proposed_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<f64, Self::Error> {
        checkpoint().map_err(HeatNetworkError::from).map_err(self.error)?;
        if self.source_epoch.map(|(time, _)| time) != Some(thermal.time_s()) {
            return Err((self.error)(HeatNetworkError::Invalid(
                "radiation timestep lacks its current source epoch")));
        }
        thermal.limit_timestep_for_additional_conductances(
            proposed_s, &self.source_conductances_w_k, checkpoint,
        ).map_err(HeatNetworkError::from).map_err(self.error)
    }

    fn event(&mut self, _: &Event, thermal: &DepositionSimulation<'_>)
        -> Result<(), Self::Error> {
        let sample = (self.resolve)(self.source, self.cx, thermal)?.bulk.source;
        check(self.cx).map_err(self.error)?;
        self.sample = Some(Arc::new(sample));
        self.source_epoch = None;
        self.source_conductances_w_k.clear();
        Ok(())
    }

    fn step(&mut self, start_s: f64, duration_s: f64, _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<(), Self::Error> {
        checkpoint().map_err(HeatNetworkError::from).map_err(self.error)?;
        let (epoch, net_w) = self.source_epoch.ok_or_else(|| (self.error)(HeatNetworkError::Invalid(
            "accepted thermal substep lacks its radiation source epoch")))?;
        if epoch != start_s || !duration_s.is_finite() || duration_s <= 0.0 {
            return Err((self.error)(HeatNetworkError::Invalid(
                "radiation accounting needs the actual source substep")));
        }
        let defect = self.closed_radiation_net_j + duration_s * net_w;
        if !defect.is_finite() || defect.abs() > self.maximum_closed_radiation_defect_j {
            return Err((self.error)(HeatNetworkError::Invalid(
                "closed radiation cannot supply unbudgeted external energy")));
        }
        let sample = self.resolve_checked(thermal, checkpoint)?.bulk.source;
        self.sample = Some(Arc::new(sample));
        self.closed_radiation_net_j = defect;
        self.source_epoch = None;
        self.source_conductances_w_k.clear();
        Ok(())
    }
}
