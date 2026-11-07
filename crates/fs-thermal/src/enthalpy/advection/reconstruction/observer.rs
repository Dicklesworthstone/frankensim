//! Limited face transport on the original explicit thermal timeline.

use super::{Advected, LimitedEnthalpyReconstruction, Sum, entropy, finite, invalid};
use crate::enthalpy::advection::{AdvectionAudit, MassBoundary, MassFlow};
use crate::enthalpy::{
    CellId, DepositionError, DepositionObserver, DepositionSimulation, Event, HeatTransfer,
};
use std::collections::BTreeMap;
use std::sync::Arc;

/// Bounded axial face enthalpies on the existing explicit thermal steps.
///
/// The original observer owns sources, mass balance, liquid carriage, paired
/// heat transfers and accepted audits. This adapter replaces only selected
/// internal donor powers using the same reconstruction as the implicit lane.
/// It shares the immutable plan across transactional observer clones.
///
/// If `a` is the decay budget represented by the prior native/wrapped observer
/// bound and `b` the largest outgoing mass rate divided by cell mass, the bound is
/// `dt <= 0.8/(a + 2*b)`. The actual nonuniform half-cell distances bound BOTH
/// traces of each donor; outgoing mass Courant at most one half then expresses
/// pure transport as a convex combination of those traces and incoming traces.
/// The original timeline retains its 0.8 safety margin and Euler time accuracy.
/// Nonlinear source stiffness and temporal accuracy still require refinement.
///
/// Nonempty plans retain the original geometry/material/stencil admissions.
/// Stationary species, selective outlets and upwind entropy histories require
/// separate reconstructed closures. The original bulk thermal entropy FACE
/// limiter is supported; it supplies no finite-step entropy history.
#[derive(Clone)]
pub struct ReconstructedAdvected<O> {
    advection: Advected<O>,
    reconstruction: Option<Arc<LimitedEnthalpyReconstruction>>,
}

impl<O> Advected<O> {
    /// Select bounded axial reconstruction on the original explicit timeline.
    /// Configure the original boundaries/carrier before attaching this plan.
    pub fn with_reconstruction(
        self,
        reconstruction: impl Into<Arc<LimitedEnthalpyReconstruction>>,
    ) -> Result<ReconstructedAdvected<O>, DepositionError> {
        self.with_optional_reconstruction(Some(reconstruction.into()))
    }

    /// Select a shared plan, or preserve the exact original upwind observer.
    /// `None` and empty plans use the original flux, timestep and audit paths.
    pub fn with_optional_reconstruction(
        self,
        reconstruction: Option<Arc<LimitedEnthalpyReconstruction>>,
    ) -> Result<ReconstructedAdvected<O>, DepositionError> {
        let reconstruction = reconstruction.filter(|plan| !plan.is_empty());
        if let Some(plan) = &reconstruction {
            // Both the authored plan and this owner's flow graph are immutable.
            // Compile their compatibility once, without per-step topology maps.
            plan.admit_flow(&self, &mut || Ok(()))?;
        }
        Ok(ReconstructedAdvected {
            advection: self,
            reconstruction,
        })
    }
}

impl<O> ReconstructedAdvected<O> {
    /// Original wrapped physical observer and its accepted history.
    pub fn inner(&self) -> &O {
        self.advection.inner()
    }
    /// Recover the original wrapped physical observer.
    pub fn into_inner(self) -> O {
        self.advection.into_inner()
    }
    /// Accepted transport diagnostics, including actual reconstructed powers.
    pub fn audit(&self) -> AdvectionAudit {
        self.advection.audit()
    }
    /// Original signed mass-flow connections.
    pub fn flows(&self) -> &[MassFlow] {
        self.advection.flows()
    }
    /// Original external mass and enthalpy streams.
    pub fn boundaries(&self) -> &[MassBoundary] {
        self.advection.boundaries()
    }
    /// Original donor-cell transport owner, including its carrier bindings.
    pub fn advection(&self) -> &Advected<O> {
        &self.advection
    }
    /// Recover the original donor-cell owner without reconstruction.
    pub fn into_advection(self) -> Advected<O> {
        self.advection
    }
    /// Active plan; absent for both `None` and empty-plan selections.
    pub fn reconstruction(&self) -> Option<&LimitedEnthalpyReconstruction> {
        self.reconstruction.as_deref()
    }
}

impl<O: DepositionObserver> DepositionObserver for ReconstructedAdvected<O> {
    type Error = O::Error;

    fn begin(&mut self) -> Result<(), Self::Error> {
        self.advection.begin()
    }

    fn next_transition_time_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Option<f64>, Self::Error> {
        self.advection.next_transition_time_s(thermal, check)
    }

    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        if let Some(plan) = &self.reconstruction {
            plan.admit_cells(thermal, check)?;
        }
        self.advection.sources(thermal, check)
    }

    fn conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<u64, f64>, Self::Error> {
        self.advection.conductances(thermal, check)
    }

    fn ambient_conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.advection.ambient_conductances(thermal, check)
    }

    fn reservoir_conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.advection.reservoir_conductances(thermal, check)
    }

    fn heat_transfers(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<HeatTransfer>, Self::Error> {
        let Some(plan) = &self.reconstruction else {
            return self.advection.heat_transfers(thermal, check);
        };
        // The native call consumes this epoch and appends its own transfers
        // AFTER any inner observer transfers. Preserve that prefix unchanged.
        let boundary_power = self
            .advection
            .boundary_pending
            .as_ref()
            .map(|pending| pending.absolute_power_w);
        let mut transfers = self.advection.heat_transfers(thermal, check)?;
        let count = self
            .advection
            .flows
            .iter()
            .filter(|flow| flow.kg_s != 0.0)
            .count();
        let offset = transfers
            .len()
            .checked_sub(count)
            .ok_or(DepositionError::Arithmetic)?;
        let mut own_transfers = transfers[offset..].iter_mut();
        let mut absolute_power = Sum::default();
        absolute_power.add(boundary_power.ok_or(DepositionError::Arithmetic)?)?;
        for flow in &self.advection.flows {
            check()?;
            let Some((up, down, rate)) = flow.directed() else {
                continue;
            };
            let transfer = own_transfers.next().ok_or(DepositionError::Arithmetic)?;
            if transfer.a != up || transfer.b != down {
                return Err(invalid("native upwind transfer order changed").into());
            }
            let mut increment = plan.face_increment(thermal, &self.advection, up, down)?;
            if plan.thermal_entropy_limit && increment.value != 0.0 {
                increment = entropy::limit_face(thermal, up, down, increment, false, None)?;
            }
            let original_power = transfer.power_w;
            if increment.value != 0.0 {
                let correction = finite(rate * increment.value)?;
                transfer.power_w = if increment.guard_power_roundoff {
                    entropy::apply_correction(original_power, correction)?
                } else {
                    finite(original_power + correction)?
                };
            }
            // Count the actual represented correction in the gross diagnostic;
            // it is already in the paired internal watts, never external heat.
            let donor_h = self.advection.carried_h(&thermal.cells()[&up])?;
            let gross_power =
                finite(finite(rate * donor_h)? + finite(transfer.power_w - original_power)?)?;
            absolute_power.add(gross_power.abs())?;
        }
        self.advection
            .pending
            .as_mut()
            .ok_or(DepositionError::Arithmetic)?
            .absolute_power_w = absolute_power.total()?;
        check()?;
        Ok(transfers)
    }

    fn limit_timestep_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        proposed_s: f64,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<f64, Self::Error> {
        let dt = self
            .advection
            .limit_timestep_s(thermal, proposed_s, check)?;
        if self.reconstruction.is_none() {
            return Ok(dt);
        }
        let decay = self
            .advection
            .pending
            .as_ref()
            .ok_or(DepositionError::Arithmetic)?
            .decay_per_s;
        if decay == 0.0 {
            return Ok(dt);
        }
        // Native advection already adds b to the shared decay budget. Add one
        // more b, rather than taking an independent minimum and spending the
        // conductive allowance twice. Retain b itself for the actual audit CFL.
        let extra_limit = 0.8 / decay;
        let result = if dt <= extra_limit {
            dt / (1.0 + dt / extra_limit)
        } else {
            extra_limit / (1.0 + extra_limit / dt)
        };
        if !result.is_finite() || result <= 0.0 || result > dt {
            return Err(DepositionError::Arithmetic.into());
        }
        Ok(result)
    }

    fn requires_step_interval(&self) -> bool {
        self.advection.requires_step_interval()
    }

    fn step_sources(
        &mut self,
        thermal_start: &DepositionSimulation<'_>,
        start_s: f64,
        dt_s: f64,
        end_s: f64,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.advection
            .step_sources(thermal_start, start_s, dt_s, end_s, check)
    }

    fn event(
        &mut self,
        event: &Event,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<(), Self::Error> {
        self.advection.event(event, thermal)
    }

    fn step(
        &mut self,
        start_s: f64,
        dt_s: f64,
        old: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.advection.step(start_s, dt_s, old, thermal, check)
    }
}
