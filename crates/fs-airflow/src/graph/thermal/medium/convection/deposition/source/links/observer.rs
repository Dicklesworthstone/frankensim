//! Finite fluid/wall convection on the original explicit thermal timeline.
use super::super::observer::with_checkpoint;
use super::*;
use fs_conduction::transient::deposition::{DepositionObserver, Event};
use std::collections::BTreeMap;
use std::sync::Arc;

/// Native conductance owner for material-dependent finite fluid interfaces.
///
/// Compose with independent heat, advection or mechanical observers through
/// the existing `CoupledObservers`. One original link carries both sides of
/// each exchange, its area coverage and its contribution to the thermal CFL.
/// No external heat rate or second fluid thermal history is introduced.
///
/// The resolver must obey the original observer's clone/rollback contract:
/// retain immutable material laws and deeply copy any mutable physical state.
/// Its material, pressure and composition must agree with the bound fluid's
/// caloric law. This adapter does not infer those physical relationships.
///
/// Coefficients are reevaluated at every actual explicit substep, and the
/// original native timestep bound includes their conductance. Nonlinear
/// coefficient-freezing error and spatial lumping still require refinement.
pub struct SourceLinkConvectionObserver<'b, R, F> {
    bindings: &'b [SourceConvectiveLink],
    maximum_links: usize,
    resolve: R,
    pending: Option<(f64, Arc<ResolvedLinkConvection<F>>)>,
    applied: Option<Arc<ResolvedLinkConvection<F>>>,
    endpoint: Option<Arc<ResolvedLinkConvection<F>>>,
}
impl<R: Clone, F> Clone for SourceLinkConvectionObserver<'_, R, F> {
    fn clone(&self) -> Self {
        Self {
            bindings: self.bindings,
            maximum_links: self.maximum_links,
            resolve: self.resolve.clone(),
            pending: self.pending.clone(),
            applied: self.applied.clone(),
            endpoint: self.endpoint.clone(),
        }
    }
}
impl<'b, R, F> SourceLinkConvectionObserver<'b, R, F> {
    /// Bind existing native links without changing time, area or thermal state.
    pub fn new(bindings: &'b [SourceConvectiveLink], maximum_links: usize, resolve: R) -> Self {
        Self {
            bindings,
            maximum_links,
            resolve,
            pending: None,
            applied: None,
            endpoint: None,
        }
    }
    /// Actual starting-state sources/conductances used on the last accepted
    /// substep. Later events do not rewrite this historical evaluation.
    pub fn applied_sample(&self) -> Option<&ResolvedLinkConvection<F>> {
        self.applied.as_deref()
    }
    /// Fresh material evaluation at the last accepted endpoint. Events
    /// invalidate this value until the next substep resolves the new state.
    pub fn sample(&self) -> Option<&ResolvedLinkConvection<F>> {
        self.endpoint.as_deref()
    }
}
impl<R, F, E> SourceLinkConvectionObserver<'_, R, F>
where
    R: FnMut(&SourceConvectiveLink, &DepositionSimulation<'_>, Temperature) -> Result<F, E>,
    F: ThermalFluidTransport,
{
    fn resolve_checked(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ResolvedLinkConvection<F>, SourceCellConvectionError<E>> {
        with_checkpoint(checkpoint, |cancelled| {
            evaluate_source_link_convection(
                thermal,
                self.bindings,
                self.maximum_links,
                &mut self.resolve,
                cancelled,
            )
        })
    }
}
impl<R, F, E> DepositionObserver for SourceLinkConvectionObserver<'_, R, F>
where
    R: Clone + FnMut(&SourceConvectiveLink, &DepositionSimulation<'_>, Temperature) -> Result<F, E>,
    F: ThermalFluidTransport,
{
    type Error = SourceCellConvectionError<E>;

    fn begin(&mut self) -> Result<(), Self::Error> {
        self.pending = None;
        Ok(())
    }
    fn conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<u64, f64>, Self::Error> {
        let source = self.resolve_checked(thermal, checkpoint)?;
        let values = source.conductances_w_k().clone();
        checkpoint()?;
        self.pending = Some((thermal.time_s(), Arc::new(source)));
        Ok(values)
    }
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        // Same-time births arrive individually. Source admission occurs once
        // the original timeline has completed their active state.
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
        let (epoch, _) = self.pending.as_ref().ok_or(ConvectionSourceError::Invalid(
            "accepted fluid-link step lacks its conductance epoch",
        ))?;
        if *epoch != start_s || !duration_s.is_finite() || duration_s <= 0.0 {
            return Err(ConvectionSourceError::Invalid(
                "fluid-link conductances need the actual thermal substep",
            )
            .into());
        }
        let endpoint = self.resolve_checked(thermal, checkpoint)?;
        checkpoint()?;
        self.applied = self.pending.take().map(|(_, sample)| sample);
        self.endpoint = Some(Arc::new(endpoint));
        Ok(())
    }
}
