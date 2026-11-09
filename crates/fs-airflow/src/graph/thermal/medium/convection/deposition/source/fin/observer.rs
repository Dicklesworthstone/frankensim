//! The original explicit timeline owns all fin, base and fluid heat storage.
use super::*;
use fs_thermal::enthalpy::{DepositionObserver, Event};
use std::sync::Arc;

/// Source fin coefficients on actual native enthalpy substeps.
///
/// Every start and accepted endpoint re-resolves the original solid and fluid
/// sources. The native timeline applies reciprocal link heat and its outgoing
/// conductance CFL. Refine its timestep and axial mesh separately; the frozen
/// coefficient CFL is not a nonlinear time-error certificate.
///
/// Clone must deeply copy any mutable callback state. Immutable shared source
/// laws are supported; external callback effects cannot be rolled back. Compose
/// other disjoint physical owners through the existing CoupledObservers.
pub struct SourceTransientFinObserver<'p, 't, 'a, 'm, R, F> {
    model: &'p PreparedSourceTransientFin<'t, 'a, 'm>,
    resolve_fluid: R,
    pending: Option<(f64, Arc<SourceTransientFinSample<F>>)>,
    applied: Option<Arc<SourceTransientFinSample<F>>>,
    endpoint: Option<Arc<SourceTransientFinSample<F>>>,
}
impl<R: Clone, F> Clone for SourceTransientFinObserver<'_, '_, '_, '_, R, F> {
    fn clone(&self) -> Self {
        Self {
            model: self.model,
            resolve_fluid: self.resolve_fluid.clone(),
            pending: self.pending.clone(),
            applied: self.applied.clone(),
            endpoint: self.endpoint.clone(),
        }
    }
}
impl<'p, 't, 'a, 'm, R, F> SourceTransientFinObserver<'p, 't, 'a, 'm, R, F> {
    /// Source admission occurs at actual native substeps. Call model.resolve
    /// directly to admit a zero-duration diagnostic without evolving any state.
    pub fn new(model: &'p PreparedSourceTransientFin<'t, 'a, 'm>, resolve_fluid: R) -> Self {
        Self {
            model,
            resolve_fluid,
            pending: None,
            applied: None,
            endpoint: None,
        }
    }
    /// Actual sources used for the most recent accepted explicit interval.
    pub fn applied_sample(&self) -> Option<&SourceTransientFinSample<F>> {
        self.applied.as_deref()
    }
    /// Fresh sources at the most recent accepted endpoint; events invalidate it.
    pub fn sample(&self) -> Option<&SourceTransientFinSample<F>> {
        self.endpoint.as_deref()
    }
}
impl<R, F, E> DepositionObserver for SourceTransientFinObserver<'_, '_, '_, '_, R, F>
where
    R: Clone + FnMut(&SourceConvectiveLink, &DepositionSimulation<'_>, Temperature) -> Result<F, E>,
    F: ThermalFluidTransport,
{
    type Error = SourceTransientFinError<E>;

    fn begin(&mut self) -> Result<(), Self::Error> {
        self.pending = None;
        Ok(())
    }

    fn conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<u64, f64>, Self::Error> {
        let source = self
            .model
            .resolve(thermal, &mut self.resolve_fluid, &mut *checkpoint)?;
        checkpoint()?;
        self.pending = Some((thermal.time_s(), Arc::new(source.payload)));
        Ok(source.conductances_w_k)
    }

    fn event(&mut self, _: &Event, thermal: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        self.pending = None;
        self.endpoint = None;
        // A target-time event still has to preserve the prepared geometry and
        // calorics, even when no later source substep will run. Native event
        // dispatch polls cancellation around this bounded admission.
        self.model.admit_materials(thermal, &mut || Ok(()))
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
        let epoch = self.pending.as_ref().map(|(epoch, _)| *epoch).ok_or(
            SourceTransientFinError::Invalid("fin interval lacks its original source epoch"),
        )?;
        if epoch != start_s || !duration_s.is_finite() || duration_s <= 0.0 {
            return Err(SourceTransientFinError::Invalid(
                "fin sources need the actual native interval",
            ));
        }
        let endpoint = self
            .model
            .resolve(thermal, &mut self.resolve_fluid, &mut *checkpoint)?;
        checkpoint()?;
        self.applied = self.pending.take().map(|(_, sample)| sample);
        self.endpoint = Some(Arc::new(endpoint.payload));
        Ok(())
    }
}
