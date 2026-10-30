//! Material convection on the original explicit enthalpy substeps.
//! The timeline owns integration, events, caloric storage and acceptance.

use std::collections::BTreeMap;
use std::sync::Arc;

use fs_conduction::transient::deposition::{
    CellId, DepositionError, DepositionObserver, DepositionSimulation, Event,
};
use fs_material::fluid_transport::ThermalFluidTransport;
use fs_qty::Temperature;

use super::super::ConvectionSourceError;
use super::{
    ResolvedCellConvectionSources, SourceCellConvectionError, SourceConvectiveCell,
    evaluate_source_cell_convection,
};

/// Re-query fluid material properties on the SAME thermal timeline as other
/// deposition observers. Compose this with existing heat, solid or contact
/// observers through `fs_thermal::observers::CoupledObservers`.
///
/// Every substep evaluates its original start-state wall temperature and exposed
/// area, then the selected reservoir or film property temperature. The original
/// thermal integrator applies that power once. The accepted endpoint retains a
/// fresh source evaluation; its power is NOT retroactively used for the step.
///
/// The prescribed fluid reservoir is external, with no invented finite heat
/// capacity. Pressure, composition and source selection remain with `resolve`.
/// Retain a prepared material law or gas model in that resolver to preserve
/// source identity. Its Clone must deeply copy any mutable physical history;
/// shared immutable laws are valid, shared mutable state or external effects
/// cannot be rolled back by the native thermal owner.
///
/// This follows the existing explicit timeline, not an implicit thermal solve.
/// Its start-state hA tightens the native conductive timestep through the same
/// caloric slope/mass bound. This is only a frozen-coefficient bound, not a
/// general nonlinear-stability or time-error guarantee: still refine
/// `StepConfig::max_dt_s` for the convective source. Bindings must name
/// active cells whenever sources are evaluated. Events invalidate the retained
/// sample; the next substep reads the updated contact coverage and boundaries.
pub struct SourceConvectionObserver<'b, R, F> {
    bindings: &'b [SourceConvectiveCell],
    maximum_cells: usize,
    resolve: R,
    sample: Option<Arc<ResolvedCellConvectionSources<F>>>,
    source_epoch: Option<(f64, f64)>,
    source_conductances_w_k: BTreeMap<CellId, f64>,
    external_heat_j: f64,
}

// Source receipts are immutable, so cloning a candidate does not require the
// resolved material type itself to implement Clone.
impl<R: Clone, F> Clone for SourceConvectionObserver<'_, R, F> {
    fn clone(&self) -> Self {
        Self {
            bindings: self.bindings,
            maximum_cells: self.maximum_cells,
            resolve: self.resolve.clone(),
            sample: self.sample.clone(),
            source_epoch: self.source_epoch,
            source_conductances_w_k: self.source_conductances_w_k.clone(),
            external_heat_j: self.external_heat_j,
        }
    }
}

impl<'b, R, F> SourceConvectionObserver<'b, R, F> {
    /// Bind the existing source evaluator without advancing thermal state.
    /// Physical admission occurs at the actual native source/endpoint queries.
    pub fn new(bindings: &'b [SourceConvectiveCell], maximum_cells: usize, resolve: R) -> Self {
        Self {
            bindings,
            maximum_cells,
            resolve,
            sample: None,
            source_epoch: None,
            source_conductances_w_k: BTreeMap::new(),
            external_heat_j: 0.0,
        }
    }

    /// Material receipts and heat rates at the last accepted substep endpoint.
    /// None before the first substep or after an event invalidates that sample.
    pub fn sample(&self) -> Option<&ResolvedCellConvectionSources<F>> {
        self.sample.as_deref()
    }

    /// Signed heat into the bodies during the most recent candidate advance.
    /// This is already included ONCE in the native EnergyAudit::external_heat_j;
    /// it identifies convection's share, not an additional energy contribution.
    pub const fn external_heat_j(&self) -> f64 {
        self.external_heat_j
    }
}

fn invalid<E>(message: &'static str) -> SourceCellConvectionError<E> {
    SourceCellConvectionError::Convection(ConvectionSourceError::Invalid(message))
}

impl<R, F, E> SourceConvectionObserver<'_, R, F>
where
    R: FnMut(&SourceConvectiveCell, Temperature) -> Result<F, E>,
    F: ThermalFluidTransport,
{
    fn resolve_checked(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ResolvedCellConvectionSources<F>, SourceCellConvectionError<E>> {
        with_checkpoint(checkpoint, |cancelled| evaluate_source_cell_convection(
            thermal, self.bindings, self.maximum_cells, &mut self.resolve, cancelled))

    }
}

impl<R, F, E> DepositionObserver for SourceConvectionObserver<'_, R, F>
where
    R: Clone + FnMut(&SourceConvectiveCell, Temperature) -> Result<F, E>,
    F: ThermalFluidTransport,
{
    type Error = SourceCellConvectionError<E>;

    fn begin(&mut self) -> Result<(), Self::Error> {
        self.source_epoch = None;
        self.source_conductances_w_k.clear();
        self.external_heat_j = 0.0;
        Ok(())
    }

    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        let evaluated = self.resolve_checked(thermal, checkpoint)?;
        let heat = evaluated.convection();
        let power_w = heat.total_heat_into_cells_w();
        let sources = heat.sources_w().clone();
        let mut conductances = BTreeMap::new();
        for sample in heat.samples() {
            checkpoint()?;
            let g = sample.evaluation.as_ref().map_or(0.0, |evaluation| {
                evaluation.coefficient.value.value() * sample.exposed_area_m2
            });
            conductances.insert(sample.cell, g);
        }
        checkpoint()?;
        self.source_epoch = Some((thermal.time_s(), power_w));
        self.source_conductances_w_k = conductances;
        Ok(sources)
    }

    fn limit_timestep_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        proposed_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<f64, Self::Error> {
        checkpoint()?;
        if self.source_epoch.map(|(time, _)| time) != Some(thermal.time_s()) {
            return Err(invalid("convection timestep lacks its current source epoch"));
        }
        Ok(thermal.limit_timestep_for_additional_conductances(
            proposed_s, &self.source_conductances_w_k, checkpoint,
        )?)
    }

    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        // Events arrive individually, including consecutive same-time births.
        // Do not query an incomplete active set between those events.
        self.sample = None;
        self.source_epoch = None;
        self.source_conductances_w_k.clear();
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
        let (epoch, power_w) = self.source_epoch.ok_or_else(|| {
            invalid("accepted convection substep lacks its original source epoch")
        })?;
        if epoch != start_s || !duration_s.is_finite() || duration_s <= 0.0 {
            return Err(invalid(
                "convection heat accounting needs the actual source interval",
            ));
        }
        let external_heat_j = self.external_heat_j + duration_s * power_w;
        if !external_heat_j.is_finite() {
            return Err(invalid("integrated convective heat is unrepresentable"));
        }
        let sample = self.resolve_checked(thermal, checkpoint)?;
        checkpoint()?;
        self.sample = Some(Arc::new(sample));
        self.external_heat_j = external_heat_j;
        self.source_epoch = None;
        self.source_conductances_w_k.clear();
        Ok(())
    }
}

// Preserve native budget/domain errors through both source evaluators' Boolean
// cancellation seam. No source query or candidate state escapes an interruption.
pub(super) fn with_checkpoint<T, E>(
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    evaluate: impl FnOnce(&mut dyn FnMut() -> bool)
        -> Result<T, SourceCellConvectionError<E>>,
) -> Result<T, SourceCellConvectionError<E>> {
    checkpoint()?;
    let mut interrupted = None;
    let result = evaluate(&mut || {
        if interrupted.is_some() { return true; }
        if let Err(error) = checkpoint() {
            interrupted = Some(error);
            return true;
        }
        false
    });
    if let Some(error) = interrupted { return Err(error.into()); }
    checkpoint()?;
    result
}
