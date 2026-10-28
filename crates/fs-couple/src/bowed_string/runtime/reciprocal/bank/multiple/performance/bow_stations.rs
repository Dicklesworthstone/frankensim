//! Source-bound bow-position trajectories on the original performance clock.
pub use super::super::stations::BowStationMotionLimits;
use super::motion::{GeometryKnot, prepare_path};
use super::*;
use crate::render::schedule::force::coupled::ModalCouplingError;
use fs_exec::CancelGate;

/// One physical endpoint of every original bow/contact location.
#[derive(Clone, Debug, PartialEq)]
pub struct BowStationKnot {
    /// Positive sample endpoint relative to performance start. The implicit
    /// zero endpoint is each contact's original source location.
    pub sample: u64,
    /// Interior fractions of original supported spans, in actual bow order.
    /// For a ribbon all cells follow one translation at its original width.
    pub station_fractions: Vec<f64>,
}
impl GeometryKnot for BowStationKnot {
    fn sample(&self) -> u64 {
        self.sample
    }
    fn coordinates(&self) -> &[f64] {
        &self.station_fractions
    }
}

/// A finite piecewise-linear bow-position path, independent of speed/load tracks.
#[derive(Clone, Debug, PartialEq)]
pub struct BowStationProgram {
    /// Strictly increasing positive endpoints; the final location is held.
    pub knots: Vec<BowStationKnot>,
    /// Physical migration and complete per-step geometry work admission.
    pub limits: BowStationMotionLimits,
    /// Maximum number of authored knots.
    pub max_knots: usize,
    /// Maximum scalar station values across the complete path.
    pub max_coordinate_values: usize,
}
impl BowStationProgram {
    /// Admit the complete path before material I/O or numerical realization.
    /// Actual source stations and span lengths must be supplied, one per bow.
    /// Binding repeats this against the real bank and checks ribbon geometry;
    /// each physical step also checks the rounded endpoint and source map.
    pub fn validate(
        &self,
        initial: &[f64],
        lengths_m: &[f64],
        samples: u64,
        rate: u32,
        gate: &CancelGate,
    ) -> Result<(), RenderError> {
        if gate.is_requested() {
            return Err(RenderError::Coupled(ModalCouplingError::Cancelled));
        }
        self.limits.validate().map_err(RenderError::Friction)?;
        if initial.is_empty()
            || lengths_m.len() != initial.len()
            || samples == 0
            || rate == 0
            || samples > (1_u64 << 53)
            || self.knots.is_empty()
            || self.knots.len() > self.max_knots
            || initial
                .len()
                .checked_mul(self.knots.len())
                .is_none_or(|n| n > self.max_coordinate_values)
            || initial
                .iter()
                .any(|s| !s.is_finite() || !(0.0 < *s && *s < 1.0))
            || lengths_m.iter().any(|l| !l.is_finite() || *l <= 0.0)
        {
            return Err(render_invalid(
                "bow station paths need actual interior geometry, a bounded clock and nonempty knots",
            ));
        }
        let mut previous: Option<&BowStationKnot> = None;
        for knot in &self.knots {
            if gate.is_requested() {
                return Err(RenderError::Coupled(ModalCouplingError::Cancelled));
            }
            let start = previous.map_or(0, |p| p.sample);
            if knot.sample <= start
                || knot.sample > samples
                || knot.station_fractions.len() != initial.len()
            {
                return Err(render_invalid(
                    "bow station knots must be ordered complete endpoints inside the performance",
                ));
            }
            let steps = (knot.sample - start) as f64;
            for (i, &end) in knot.station_fractions.iter().enumerate() {
                let begin = previous.map_or(initial[i], |p| p.station_fractions[i]);
                let distance = (end - begin) * lengths_m[i] / steps;
                let speed = distance * f64::from(rate);
                if !end.is_finite()
                    || !(0.0 < end && end < 1.0)
                    || !distance.is_finite()
                    || distance.abs() > self.limits.maximum_step_m
                    || !speed.is_finite()
                    || speed.abs() > self.limits.maximum_speed_m_s
                {
                    return Err(render_invalid(
                        "bow trajectory exceeds its interior span, distance or speed limit",
                    ));
                }
            }
            previous = Some(knot);
        }
        Ok(())
    }
}

pub(super) struct BoundBowStations {
    source: BowStationProgram,
    initial: Vec<f64>,
    pub(super) candidate: Vec<f64>,
    pub(super) next: usize,
}
impl BoundBowStations {
    pub(super) fn prepare(&mut self, endpoint: u64) -> usize {
        prepare_path(
            &self.source.knots,
            self.next,
            endpoint,
            Some(&self.initial),
            &mut self.candidate,
        )
    }
}

impl MultiBowPerformance {
    /// Bind bow-position motion to the SAME samples as bow loads/speeds, forces
    /// and independently driven normal stops. All initial vibration and source
    /// histories survive. The original pressure renderer and its optional
    /// exterior observer consume the resulting mechanical state.
    ///
    /// A finite ribbon's complete interval must stay inside its original string;
    /// every spatial cell translates by the same amount, to rounding precision.
    /// Hair contact coordinates, support frame and material free length remain
    /// those in the original description. A contact-loaded normal follows its
    /// bow automatically through the same product-rule contact-work owner.
    /// Independently authored stop paths keep mapped normal entries at their
    /// original stations; this program owns their actual shared trajectories.
    /// Retained shear/thermal layers still require patch-history transport and
    /// refuse migration.
    pub fn with_bow_station_motion(
        mut self,
        program: BowStationProgram,
        gate: &CancelGate,
    ) -> Result<Self, RenderError> {
        if gate.is_requested() {
            return Err(RenderError::Coupled(ModalCouplingError::Cancelled));
        }
        self.validate_sample_count(0)?;
        if self.bank.system.samples_rendered() != 0
            || self.bow_station_motion.is_some()
            || self.bank.bow_station_motion_limits().is_some()
        {
            return Err(render_invalid(
                "bow station paths bind once before the physical performance begins",
            ));
        }
        let initial: Vec<_> = self
            .bank
            .source
            .bows
            .iter()
            .map(|b| b.location.station_fraction)
            .collect();
        let lengths: Vec<_> = self
            .bank
            .source
            .bows
            .iter()
            .map(|b| self.bank.cards[b.location.string].length_m)
            .collect();
        program.validate(
            &initial,
            &lengths,
            self.config.samples,
            self.bank.source.bank().options().sample_rate_hz,
            gate,
        )?;
        // An independently admitted stop path may retain tighter normal-map
        // limits. Validate the bow-owned subset against that same trajectory
        // owner now, rather than discovering a conflict after audio begins.
        let contacts = &self.bank.source.bow_contact_loads;
        if let Some(limits) = self.bank.stop_station_motion_limits().filter(|_| !contacts.is_empty()) {
            let normal_program = super::stations::StopStationProgram {
                knots: program.knots.iter().map(|knot| super::stations::StopStationKnot {
                    sample: knot.sample,
                    station_fractions: contacts.iter().map(|binding| knot.station_fractions[binding.bow]).collect(),
                }).collect(),
                limits, max_knots: program.max_knots, max_coordinate_values: program.max_coordinate_values,
            };
            let initial: Vec<_> = contacts.iter().map(|binding| initial[binding.bow]).collect();
            let lengths: Vec<_> = contacts.iter().map(|binding| lengths[binding.bow]).collect();
            normal_program.validate(&initial, &lengths, self.config.samples,
                self.bank.source.bank().options().sample_rate_hz, gate)?;
        }
        if let Some(ribbon) = self.bank.source.ribbon() {
            for knot in &program.knots {
                super::super::stations::validate_ribbon_stations(
                    ribbon,
                    self.bank.cards[ribbon.centre().string].length_m,
                    &knot.station_fractions,
                )
                .map_err(RenderError::Friction)?;
            }
        }
        self.bank = self.bank.with_bow_station_motion(program.limits, gate)?;
        self.bow_station_motion = Some(BoundBowStations {
            source: program,
            candidate: initial.clone(),
            initial,
            next: 0,
        });
        Ok(self)
    }

    /// Original finite path and physical limits; no mutable geometry bypass.
    #[must_use]
    pub fn bow_station_program(&self) -> Option<&BowStationProgram> {
        self.bow_station_motion.as_ref().map(|p| &p.source)
    }
}
