//! Source-derived longitudinal stop geometry on the original performance clock.
use super::*;
use fs_exec::CancelGate;
use crate::render::schedule::force::coupled::ModalCouplingError;
use super::motion::{GeometryKnot, prepare_path};
pub use super::super::stops::stations::StopStationMotionLimits;

/// One endpoint of a prescribed path along each stop's original string span.
#[derive(Clone, Debug, PartialEq)]
pub struct StopStationKnot {
    /// Elapsed endpoint sample, positive and strictly increasing. Sample zero
    /// is the ORIGINAL source station at every contact, not zero position.
    pub sample: u64,
    /// Interior fractions of full supported spans, in original stop order.
    /// These are geometric positions, never speaking lengths or desired pitches.
    pub station_fractions: Vec<f64>,
}
impl GeometryKnot for StopStationKnot {
    fn sample(&self) -> u64 { self.sample }
    fn coordinates(&self) -> &[f64] { &self.station_fractions }
}

/// Bounded piecewise-linear longitudinal trajectories of existing compliant stops.
#[derive(Clone, Debug, PartialEq)]
pub struct StopStationProgram {
    /// Exact authored endpoints; the final position is held.
    pub knots: Vec<StopStationKnot>,
    /// Physical speed/distance and per-sample geometry-work admission.
    pub limits: StopStationMotionLimits,
    /// Maximum declared knots, checked before runtime allocation.
    pub max_knots: usize,
    /// Maximum scalar station values across the complete source path.
    pub max_coordinate_values: usize,
}
impl StopStationProgram {
    /// Complete path admission using one initial station and actual full span
    /// length per stop. Native loaders can call this before material pack I/O;
    /// binding repeats it against the realized source bank. Every step also
    /// validates its rounded endpoint and the original contact compliance.
    pub fn validate(&self, initial: &[f64], lengths_m: &[f64], samples: u64,
        rate: u32, gate: &CancelGate) -> Result<(), RenderError>
    {
        if gate.is_requested() { return Err(RenderError::Coupled(ModalCouplingError::Cancelled)); }
        self.limits.validate()?;
        let stops = initial.len();
        if stops == 0 || lengths_m.len() != stops || samples == 0 || rate == 0
            || samples > (1_u64 << 53) || self.knots.is_empty() || self.knots.len() > self.max_knots
            || stops.checked_mul(self.knots.len()).is_none_or(|n| n > self.max_coordinate_values)
            || initial.iter().any(|s| !s.is_finite() || !(0.0 < *s && *s < 1.0))
            || lengths_m.iter().any(|l| !l.is_finite() || *l <= 0.0) {
            return Err(render_invalid("stop station paths need actual interior geometry, a bounded clock and nonempty knots"));
        }
        let mut previous: Option<&StopStationKnot> = None;
        for knot in &self.knots {
            if gate.is_requested() { return Err(RenderError::Coupled(ModalCouplingError::Cancelled)); }
            let start = previous.map_or(0, |k| k.sample);
            if knot.sample <= start || knot.sample > samples || knot.station_fractions.len() != stops {
                return Err(render_invalid("station knots must be ordered complete endpoints inside the performance"));
            }
            let steps = (knot.sample-start) as f64;
            for (i, &end) in knot.station_fractions.iter().enumerate() {
                let begin = previous.map_or(initial[i], |k| k.station_fractions[i]);
                let step = (end-begin)*lengths_m[i]/steps;
                let speed = step*f64::from(rate);
                if !end.is_finite() || !(0.0 < end && end < 1.0) || !step.is_finite()
                    || step.abs() > self.limits.maximum_step_m || !speed.is_finite()
                    || speed.abs() > self.limits.maximum_speed_m_s {
                    return Err(render_invalid("station trajectory exceeds its interior span, distance or speed limit"));
                }
            }
            previous = Some(knot);
        }
        Ok(())
    }
}

pub(super) struct BoundStopStations {
    source: StopStationProgram,
    initial: Vec<f64>,
    pub(super) candidate: Vec<f64>,
    pub(super) next: usize,
}
impl BoundStopStations {
    pub(super) fn prepare(&mut self, endpoint: u64) -> usize {
        prepare_path(&self.source.knots, self.next, endpoint, Some(&self.initial), &mut self.candidate)
    }
}

impl MultiBowPerformance {
    /// Bind longitudinal stop paths to the SAME steps as bow, external-force
    /// and normal-support controls. The full source state and original pressure
    /// renderer survive. No frequency retuning, state projection or audio-layer
    /// pitch shift is used. Callback partitioning does not alter the path.
    ///
    /// A contact supplying a bow's normal load keeps its ORIGINAL station in
    /// this independent stop path. A separately bound bow-station program owns
    /// that whole patch and supplies its actual normal/tangent endpoint together.
    /// Independent finger/fret stops may move alongside it.
    pub fn with_stop_station_motion(mut self, program: StopStationProgram,
        gate: &CancelGate) -> Result<Self, RenderError>
    {
        if gate.is_requested() { return Err(RenderError::Coupled(ModalCouplingError::Cancelled)); }
        self.validate_sample_count(0)?;
        if self.bank.system.samples_rendered() != 0 || self.stop_station_motion.is_some()
            || self.bank.stop_station_motion_limits().is_some() {
            return Err(render_invalid("stop station paths bind once before the physical performance begins"));
        }
        let initial: Vec<_> = self.bank.source.stops.iter().map(|s| s.location.station_fraction).collect();
        let lengths: Vec<_> = self.bank.source.stops.iter().map(|s| self.bank.cards[s.location.string].length_m).collect();
        program.validate(&initial, &lengths, self.config.samples,
            self.bank.source.bank().options().sample_rate_hz, gate)?;
        for binding in &self.bank.source.bow_contact_loads {
            if program.knots.iter().any(|k| k.station_fractions[binding.stop].to_bits() != initial[binding.stop].to_bits()) {
                return Err(render_invalid("a contact-loaded normal follows its bow program; independent stop paths must retain its original station"));
            }
        }
        self.bank = self.bank.with_stop_station_motion(program.limits, gate)?;
        self.stop_station_motion = Some(BoundStopStations {
            source: program, candidate: initial.clone(), initial, next: 0,
        });
        Ok(self)
    }

    /// Original longitudinal path and limits, without a mutable geometry bypass.
    #[must_use]
    pub fn stop_station_program(&self) -> Option<&StopStationProgram> {
        self.stop_station_motion.as_ref().map(|m| &m.source)
    }
}
