//! Quasi-steady gas-gap conductance on the original finite-wall timeline.
use super::super::observer::with_checkpoint;
use super::*;
use fs_conduction::transient::deposition::{DepositionObserver, Event};
use std::sync::Arc;

/// Source-dependent gas gaps on the original explicit finite-wall integrator.
///
/// Each native link carries one reciprocal heat, excludes its area from ambient
/// exchange, and contributes the TOTAL gas conductance to the original CFL.
/// The gas source is resolved at both actual starting and accepted endpoint
/// wall states. An unsupported endpoint refuses the native step atomically.
/// Nonlinear coefficient-freezing and uniform-wall errors still need refinement.
///
/// Compose independent observers using the original `CoupledObservers`; do
/// not assign another conductance owner to a bound link. Resolver clones must
/// preserve immutable gas laws and deeply copy mutable physical state, as
/// required by the original observer's rollback contract.
pub struct SourceEnclosedGasLayerObserver<'b, R, F> {
    bindings: &'b [SourceEnclosedGasLayer],
    maximum_layers: usize,
    resolve: R,
    pending: Option<(f64, Arc<ResolvedEnclosedGasLayers<F>>)>,
    applied: Option<Arc<ResolvedEnclosedGasLayers<F>>>,
    endpoint: Option<Arc<ResolvedEnclosedGasLayers<F>>>,
}
impl<R: Clone, F> Clone for SourceEnclosedGasLayerObserver<'_, R, F> {
    fn clone(&self) -> Self {
        Self {
            bindings: self.bindings,
            maximum_layers: self.maximum_layers,
            resolve: self.resolve.clone(),
            pending: self.pending.clone(),
            applied: self.applied.clone(),
            endpoint: self.endpoint.clone(),
        }
    }
}
impl<'b, R, F> SourceEnclosedGasLayerObserver<'b, R, F> {
    pub fn new(bindings: &'b [SourceEnclosedGasLayer], maximum_layers: usize, resolve: R) -> Self {
        Self {
            bindings,
            maximum_layers,
            resolve,
            pending: None,
            applied: None,
            endpoint: None,
        }
    }
    /// Starting-state coefficient/source actually applied on the last substep.
    pub fn applied_sample(&self) -> Option<&ResolvedEnclosedGasLayers<F>> {
        self.applied.as_deref()
    }
    /// Fresh source evaluation at the last accepted endpoint; events clear it.
    pub fn sample(&self) -> Option<&ResolvedEnclosedGasLayers<F>> {
        self.endpoint.as_deref()
    }
}
impl<R, F, E> SourceEnclosedGasLayerObserver<'_, R, F>
where
    R: FnMut(&SourceEnclosedGasLayer, &DepositionSimulation<'_>, Temperature) -> Result<F, E>,
    F: ThermalFluidTransport,
{
    fn resolve_checked(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ResolvedEnclosedGasLayers<F>, SourceCellConvectionError<E>> {
        with_checkpoint(checkpoint, |cancelled| {
            evaluate_source_enclosed_gas_layers(
                thermal,
                self.bindings,
                self.maximum_layers,
                &mut self.resolve,
                cancelled,
            )
        })
    }
}
impl<R, F, E> DepositionObserver for SourceEnclosedGasLayerObserver<'_, R, F>
where
    R: Clone
        + FnMut(&SourceEnclosedGasLayer, &DepositionSimulation<'_>, Temperature) -> Result<F, E>,
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
            "accepted gas-layer step lacks its conductance epoch",
        ))?;
        if *epoch != start_s || !duration_s.is_finite() || duration_s <= 0.0 {
            return Err(ConvectionSourceError::Invalid(
                "gas-layer conductance needs the actual thermal substep",
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
