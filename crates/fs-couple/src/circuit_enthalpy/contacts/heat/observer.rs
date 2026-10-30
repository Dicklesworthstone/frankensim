//! Complete sourced heat transfer on the existing thermal observer timeline.
//! No dummy circuit, second thermal advance or alternate radiation/contact law.
use std::collections::BTreeMap;
use std::sync::Arc;
use fs_thermal::enthalpy::{CellId, DepositionObserver, Event};
use fs_thermal::observers::CoupledObservers;
use super::*;

/// Thermal-only view of one prepared material heat network. All mutable state
/// is cloned by the original timeline. The shared sample is immutable, not a
/// second temperature field or a cache of coefficients used in later steps.
#[derive(Clone)]
pub struct MaterialHeatObserver<'p, 'n, 'a, 'scope> {
    prepared: &'p PreparedHeatNetwork<'n, 'a>,
    cx: &'p Cx<'scope>,
    sample: Option<Arc<HeatNetworkSample>>,
    source_epoch: Option<(f64, f64)>,
    source_conductances_w_k: BTreeMap<CellId, f64>,
    closed_radiation_net_j: f64,
    maximum_closed_radiation_defect_j: f64,
}

impl MaterialHeatObserver<'_, '_, '_, '_> {
    /// Complete material, interface and radiation results at the last accepted
    /// endpoint/event. Source rates used to advance it were evaluated at its
    /// substep START, as required by the existing explicit thermal integrator.
    pub fn sample(&self) -> Option<&HeatNetworkSample> { self.sample.as_deref() }

    /// Integrated roundoff-net radiation in the most recent candidate advance.
    /// This is internal-exchange error, not physical external heat. The thermal
    /// audit includes it in generated sources; subtract it ONCE when interpreting
    /// that audit as external heat. It has already been bounded before acceptance.
    pub const fn closed_radiation_net_j(&self) -> f64 { self.closed_radiation_net_j }

    fn resolve_checked(
        &self, thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ContactThermalCoefficients<HeatNetworkSample>, HeatNetworkError> {
        checkpoint()?;
        // Original evaluators poll this Cx within their bounded work. Also keep
        // the timeline's richer checkpoint errors before and after evaluation.
        let result = self.prepared.resolve(self.cx, thermal)?;
        checkpoint()?;
        Ok(result)
    }
}

impl<'n, 'a> PreparedHeatNetwork<'n, 'a> {
    /// Compose sourced conduction, thermal-interface resistance and radiation
    /// with another observer on the SAME accepted thermal substeps.
    ///
    /// This enables existing thermoplastic solids, contact histories and other
    /// DepositionObserver consumers to use the complete heat network, without
    /// adding a circuit or a parallel stepping implementation. Electrical work
    /// is deliberately not supported: resistor/contact/basis bindings refuse.
    ///
    /// Sources use the original explicit timeline, NOT the circuit driver's
    /// implicit residual. The same native timestep owner combines conduction
    /// with a conservative frozen-emissivity radiation rate. This is not a
    /// temporal-error or temperature-dependent-emissivity derivative bound;
    /// callers must still refine max_dt_s. The existing solid's outer
    /// thermoplastic iteration remains unchanged.
    ///
    /// The first returned observer retains source results and the bounded net
    /// radiation defect; the second retains the caller's original history.
    /// Competing conductance owners refuse through CoupledObservers. Topology
    /// changes require complete rebinding at a boundary, not hidden reselection.
    pub fn observe<'p, 'scope, O>(
        &'p self, cx: &'p Cx<'scope>, inner: O,
        maximum_closed_radiation_defect_j: f64,
    ) -> Result<CoupledObservers<MaterialHeatObserver<'p, 'n, 'a, 'scope>, O>, HeatNetworkError> {
        check(cx)?;
        if !maximum_closed_radiation_defect_j.is_finite()
            || maximum_closed_radiation_defect_j < 0.0 {
            return Err(HeatNetworkError::Invalid("finite nonnegative radiation energy allowance required"));
        }
        self.admit_thermal_only()?;
        Ok(CoupledObservers::new(MaterialHeatObserver {
            prepared: self, cx, sample: None, source_epoch: None,
            source_conductances_w_k: BTreeMap::new(),
            closed_radiation_net_j: 0.0, maximum_closed_radiation_defect_j,
        }, inner))
    }
}

impl DepositionObserver for MaterialHeatObserver<'_, '_, '_, '_> {
    type Error = HeatNetworkError;

    fn begin(&mut self) -> Result<(), Self::Error> {
        check(self.cx)?;
        self.source_epoch = None;
        self.source_conductances_w_k.clear();
        self.closed_radiation_net_j = 0.0;
        Ok(())
    }

    fn sources(&mut self, thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<BTreeMap<CellId, f64>, Self::Error> {
        let evaluated = self.resolve_checked(thermal, checkpoint)?;
        let net = evaluated.bulk.source.radiation_net_power_w();
        if !net.is_finite() {
            return Err(HeatNetworkError::Invalid("nonfinite internal radiation balance"));
        }
        let bounds = evaluated.bulk.source.radiation.as_ref()
            .map(|radiation| radiation.conductance_bounds_w_k(&mut *checkpoint))
            .transpose()?.unwrap_or_default();
        self.source_epoch = Some((thermal.time_s(), net));
        self.source_conductances_w_k = bounds;
        Ok(evaluated.bulk.sources_w)
    }

    fn conductances(&mut self, thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<BTreeMap<u64, f64>, Self::Error> {
        // Re-query immutable laws at this actual state. Do not reuse coefficients
        // merely because a caller supplies a different state at the same time.
        Ok(self.resolve_checked(thermal, checkpoint)?.bulk.conductances_w_k)
    }

    fn limit_timestep_s(
        &mut self, thermal: &DepositionSimulation<'_>, proposed_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<f64, Self::Error> {
        checkpoint()?;
        if self.source_epoch.map(|(time, _)| time) != Some(thermal.time_s()) {
            return Err(HeatNetworkError::Invalid(
                "radiation timestep lacks its current source epoch"));
        }
        Ok(thermal.limit_timestep_for_additional_conductances(
            proposed_s, &self.source_conductances_w_k, checkpoint,
        )?)
    }

    fn event(&mut self, _: &Event, thermal: &DepositionSimulation<'_>)
        -> Result<(), Self::Error> {
        let sample = self.prepared.resolve(self.cx, thermal)?.bulk.source;
        check(self.cx)?;
        self.sample = Some(Arc::new(sample));
        self.source_epoch = None;
        self.source_conductances_w_k.clear();
        Ok(())
    }

    fn step(&mut self, start_s: f64, duration_s: f64, _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<(), Self::Error> {
        checkpoint()?;
        let (epoch, net_w) = self.source_epoch.ok_or(HeatNetworkError::Invalid(
            "accepted thermal substep lacks its radiation source epoch"))?;
        if epoch != start_s || !duration_s.is_finite() || duration_s <= 0.0 {
            return Err(HeatNetworkError::Invalid("radiation accounting needs the actual source substep"));
        }
        let defect = self.closed_radiation_net_j + duration_s * net_w;
        if !defect.is_finite() || defect.abs() > self.maximum_closed_radiation_defect_j {
            return Err(HeatNetworkError::Invalid(
                "closed radiation cannot supply unbudgeted external energy"));
        }
        let sample = self.resolve_checked(thermal, checkpoint)?.bulk.source;
        // All fallible work precedes publication into this candidate observer.
        self.sample = Some(Arc::new(sample));
        self.closed_radiation_net_j = defect;
        self.source_epoch = None;
        self.source_conductances_w_k.clear();
        Ok(())
    }
}
