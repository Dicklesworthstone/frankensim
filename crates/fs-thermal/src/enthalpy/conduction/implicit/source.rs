//! Local heat-generation feedback on the original implicit conduction graph.

use super::{ImplicitConductionChain, ImplicitConductionError};
use crate::enthalpy::{CellId, DepositionError, DepositionSimulation, ThermalLaw};
use std::collections::BTreeMap;

/// Signed power INTO the original cell and its derivative with respect to that
/// cell's specific stored energy. Positive power heats, negative power cools.
/// This is total cell power [W], not a volumetric density [W/m³].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImplicitHeatSourceState {
    pub power_w: f64,
    pub power_enthalpy_derivative_w_kg_j: f64,
}

impl ImplicitHeatSourceState {
    pub(super) fn admit(self) -> Result<(), DepositionError> {
        if !self.power_w.is_finite() || !self.power_enthalpy_derivative_w_kg_j.is_finite() {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(())
    }
}

/// Immutable local source law evaluated at the SAME state as heat conduction.
///
/// A physical source owner supplies its complete signed power and analytic h
/// derivative. Temperature feedback includes its own dT/dh chain rule; a
/// phase-dependent source can change across a latent plateau even at constant
/// temperature. Every evaluation must enforce the source's validity domain.
/// Speculative calls may neither advance kinetics nor mutate physical history.
///
/// A source may also depend on prescribed physical time. The native endpoint
/// and SDIRK stages supply their actual time, with h derivatives holding time
/// fixed. Continuous slope changes must be declared so fixed and adaptive
/// comparisons end exactly at them. Discontinuous powers require native events
/// or interval-mean observers; an endpoint callback does not supply jump sides.
///
/// This local closure does not infer an electrical circuit, reaction inventory,
/// source entropy, or cross-cell feedback. A coupled source needs its
/// corresponding physical owner. Independent authored and observed source watts
/// remain additive; the caller must not supply the same power through both.
/// Cross-cell source feedback uses
/// [`super::coupled_source::ImplicitCoupledHeatSourceLaw`] with a complete
/// sparse derivative pattern on the original caloric objects.
pub trait ImplicitHeatSourceLaw {
    fn caloric_law(&self) -> &dyn ThermalLaw;
    fn state_at_specific_enthalpy(
        &self,
        specific_h_j_kg: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitHeatSourceState, DepositionError>;

    /// Complete source power and its h derivative at the physical stage time.
    /// Autonomous laws retain their existing query and behavior.
    fn state_at_specific_enthalpy_and_time(
        &self,
        specific_h_j_kg: f64,
        time_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitHeatSourceState, DepositionError> {
        checkpoint()?;
        if !time_s.is_finite() || time_s < 0.0 {
            return Err(DepositionError::InvalidInput(
                "invalid implicit source time",
            ));
        }
        self.state_at_specific_enthalpy(specific_h_j_kg, checkpoint)
    }

    /// The next finite, strictly future change of a continuous source's slope.
    /// The native owner clips every method interval and adaptive comparison to
    /// this time. None means no further transition, not a discontinuous jump.
    fn next_transition_time_s(&self, _after_s: f64) -> Option<f64> {
        None
    }
}

/// One sample of a continuous, prescribed source-power multiplier.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HeatSourceScheduleKnot {
    pub time_s: f64,
    pub power_multiplier: f64,
}

/// Bounded continuous power history shared by any number of original sources.
///
/// Linear interpolation preserves nonnegative multipliers within the declared
/// closed time range. Queries outside that range refuse; constant tails must
/// be authored explicitly. Duplicate times, negative or nonfinite samples
/// refuse; a sharp jump requires a native event instead.
/// Scaling power rather than temperature retains the original enthalpy storage
/// and any nonlinear or phase-dependent absorption/generation feedback.
#[derive(Clone, Debug, PartialEq)]
pub struct HeatSourceSchedule {
    knots: Vec<HeatSourceScheduleKnot>,
}

impl HeatSourceSchedule {
    pub fn new(
        knots: Vec<HeatSourceScheduleKnot>,
        maximum_knots: usize,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, DepositionError> {
        checkpoint()?;
        if knots.len() < 2 {
            return Err(DepositionError::InvalidInput(
                "heat-source schedule requires at least two knots",
            ));
        }
        if knots.len() > maximum_knots {
            return Err(DepositionError::Budget("heat-source schedule samples"));
        }
        for (i, knot) in knots.iter().enumerate() {
            checkpoint()?;
            if !knot.time_s.is_finite()
                || knot.time_s < 0.0
                || !knot.power_multiplier.is_finite()
                || knot.power_multiplier < 0.0
                || (i != 0 && knot.time_s <= knots[i - 1].time_s)
            {
                return Err(DepositionError::InvalidInput(
                    "heat-source schedule needs ordered times and nonnegative finite power",
                ));
            }
        }
        checkpoint()?;
        Ok(Self { knots })
    }

    pub fn knots(&self) -> &[HeatSourceScheduleKnot] {
        &self.knots
    }

    pub fn next_transition_time_s(&self, after_s: f64) -> Option<f64> {
        self.knots
            .get(self.knots.partition_point(|knot| knot.time_s <= after_s))
            .map(|knot| knot.time_s)
    }

    /// Deterministic logarithmic lookup within the explicit source time domain.
    /// No source state or history is changed, and neither endpoint extrapolates.
    pub fn multiplier_at(&self, time_s: f64) -> Result<f64, DepositionError> {
        if !time_s.is_finite()
            || time_s < self.knots[0].time_s
            || time_s > self.knots[self.knots.len() - 1].time_s
        {
            return Err(DepositionError::InvalidInput(
                "heat-source query is outside the prescribed time domain",
            ));
        }
        let upper = self.knots.partition_point(|knot| knot.time_s <= time_s);
        let left = self.knots[upper - 1];
        if upper == self.knots.len() || time_s == left.time_s {
            return Ok(left.power_multiplier);
        }
        let right = self.knots[upper];
        let fraction = (time_s - left.time_s) / (right.time_s - left.time_s);
        // A convex sum avoids overflowing a difference between endpoint levels.
        super::finite((1.0 - fraction) * left.power_multiplier + fraction * right.power_multiplier)
    }
}

/// Apply a prescribed power history to one immutable original source law.
///
/// The same multiplier scales total watts and the complete h tangent at every
/// Newton trial. Off intervals still query/admit the original source domain;
/// switching power off cannot conceal invalid material or caloric support.
/// The borrowed schedule may be shared across cells without copying its knots.
#[derive(Clone, Copy)]
pub struct ScheduledHeatSource<'a> {
    source: &'a dyn ImplicitHeatSourceLaw,
    schedule: &'a HeatSourceSchedule,
}

impl<'a> ScheduledHeatSource<'a> {
    pub fn new(source: &'a dyn ImplicitHeatSourceLaw, schedule: &'a HeatSourceSchedule) -> Self {
        Self { source, schedule }
    }
}

impl ImplicitHeatSourceLaw for ScheduledHeatSource<'_> {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self.source.caloric_law()
    }

    fn state_at_specific_enthalpy(
        &self,
        _specific_h_j_kg: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitHeatSourceState, DepositionError> {
        checkpoint()?;
        Err(DepositionError::InvalidInput(
            "scheduled implicit heat source requires physical time",
        ))
    }

    fn state_at_specific_enthalpy_and_time(
        &self,
        specific_h_j_kg: f64,
        time_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitHeatSourceState, DepositionError> {
        checkpoint()?;
        let multiplier = self.schedule.multiplier_at(time_s)?;
        let state =
            self.source
                .state_at_specific_enthalpy_and_time(specific_h_j_kg, time_s, checkpoint)?;
        state.admit()?;
        let scaled = ImplicitHeatSourceState {
            power_w: state.power_w * multiplier,
            power_enthalpy_derivative_w_kg_j: state.power_enthalpy_derivative_w_kg_j * multiplier,
        };
        scaled.admit()?;
        checkpoint()?;
        Ok(scaled)
    }

    fn next_transition_time_s(&self, after_s: f64) -> Option<f64> {
        earlier_transition(
            self.source.next_transition_time_s(after_s),
            self.schedule.next_transition_time_s(after_s),
            after_s,
        )
    }
}

// Preserve invalid owner transitions for the native admission check. f64::min
// alone could silently discard a NaN and step across an undeclared source edge.
fn earlier_transition(first: Option<f64>, second: Option<f64>, after_s: f64) -> Option<f64> {
    match (first, second) {
        (Some(time), _) if !time.is_finite() || time <= after_s => Some(time),
        (_, Some(time)) if !time.is_finite() || time <= after_s => Some(time),
        (Some(a), Some(b)) => Some(a.min(b)),
        (Some(time), None) | (None, Some(time)) => Some(time),
        (None, None) => None,
    }
}

impl ImplicitConductionChain<'_> {
    pub(super) fn next_forcing_transition_time_s(&self, after_s: f64) -> Option<f64> {
        let local = self.heat_sources.values().fold(
            self.next_advection_transition_time_s(after_s),
            |next, source| {
                earlier_transition(next, source.next_transition_time_s(after_s), after_s)
            },
        );
        let coupled = self.coupled_heat_source.as_ref()
            .and_then(|source| source.law.next_transition_time_s(after_s));
        earlier_transition(local, coupled, after_s)
    }
}

impl<'a> ImplicitConductionChain<'a> {
    /// Attach a bounded subset of local endpoint heat sources.
    ///
    /// Original caloric objects, cell identities and source domains are checked
    /// before returning. Empty, duplicate, unknown or repeated attachments
    /// refuse. Sources join the existing joule residual and diagonal Newton
    /// tangent in backward Euler, both SDIRK2 stages and adaptive advancement.
    /// They do not introduce another timestepper or an extra energy credit.
    /// Positive feedback can destroy equilibrium or make a chosen interval
    /// inadmissible; the existing nonlinear retry and work limits still apply.
    pub fn with_heat_sources(
        mut self,
        thermal: &DepositionSimulation<'_>,
        sources: &[(CellId, &'a dyn ImplicitHeatSourceLaw)],
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, ImplicitConductionError> {
        checkpoint()?;
        if !self.heat_sources.is_empty() || sources.is_empty() || sources.len() > self.ids.len() {
            return Err(ImplicitConductionError::Invalid(
                "one bounded nonempty implicit heat-source attachment is required",
            ));
        }
        self.admit(thermal, &mut checkpoint)?;
        let mut bindings = BTreeMap::new();
        for &(id, law) in sources {
            checkpoint()?;
            let cell = thermal
                .cells
                .get(&id)
                .ok_or(ImplicitConductionError::Invalid(
                    "implicit heat source requires an active original cell",
                ))?;
            if !thermal.uses_material_law(cell.material_id, law.caloric_law())
                || bindings.insert(id, law).is_some()
            {
                return Err(ImplicitConductionError::Invalid(
                    "implicit heat source requires a unique cell and its original caloric object",
                ));
            }
        }
        self.heat_sources = bindings;
        self.endpoint(thermal, &mut checkpoint)?;
        checkpoint()?;
        Ok(self)
    }
}
