//! Sample-timed bow performances over the existing persistent string solver.
//!
//! Compilation uses `GestureSchedule::sample_value`, including interrupted ramps;
//! it is not another interpolator. Control tick k applies at ceil(k * audio / control)
//! relative to the voice clock at binding. Each three-coordinate bow value is held
//! until the next tick. State, contact physics, substeps and acoustic observations
//! remain owned by `BowedStringState`. No note-to-frequency mapping is introduced.

use super::{BowedRenderOutcome, BowedRunError, BowedSample, BowedStringState};
use fs_exec::CancelGate;
use fs_scenario::gesture::{GestureError, GestureSchedule, GestureTarget, GestureValue};

/// One complete bow input, applied before its absolute voice-clock sample.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ScheduledBowControl {
    /// Absolute sample index, including samples completed before binding.
    pub sample: u64,
    /// Signed transverse bow speed [m/s]; reversals are physical inputs.
    pub velocity_m_s: f64,
    /// Nonnegative compressive load [N]; zero lifts the bow without silencing it.
    pub normal_force_n: f64,
    /// Station as the fraction of speaking length used by the existing bow model.
    pub station: f64,
}

/// Source, admission and runtime failures remain distinguishable.
#[derive(Debug)]
pub enum BowedScheduleError {
    /// The source sampler refused a track or value.
    Gesture(GestureError),
    /// The physical voice refused admission, a control, or a sample.
    Voice(BowedRunError),
    /// Clocks, sizes, horizons or the selected track are incompatible.
    Invalid { /// What must be corrected.
        what: &'static str },
    /// The stateless source sampler would exceed the explicit visit budget.
    WorkBudget { /// Conservative track/event visits required.
        required: u128, /// Caller-supplied limit.
        allowed: u64 },
    /// More distinct bow controls are needed than the storage budget allows.
    EventBudget { /// Caller-supplied maximum number of controls.
        allowed: usize },
}

impl core::fmt::Display for BowedScheduleError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Gesture(error) => write!(f, "bow schedule source: {error}"),
            Self::Voice(error) => write!(f, "bow schedule voice: {error:?}"),
            Self::Invalid { what } => write!(f, "bow schedule: {what}"),
            Self::WorkBudget { required, allowed } => write!(f,
                "bow schedule needs at most {required} sampling visits; budget is {allowed}"),
            Self::EventBudget { allowed } => write!(f,
                "bow schedule exceeds the {allowed}-control storage budget"),
        }
    }
}
impl std::error::Error for BowedScheduleError {}

fn invalid(what: &'static str) -> BowedScheduleError {
    BowedScheduleError::Invalid { what }
}

/// A finite, explicitly selected bow track bound to one retained physical voice.
///
/// Other source tracks are not executed: callers select the track for this voice
/// by ID. Exclusive ownership prevents stepping or changing inputs behind the
/// schedule. An already vibrating voice is supported; schedule time zero binds
/// to its current sample, without resetting its modes or its physical clock.
///
/// Compilation admits the complete horizon before any input changes. Interface
/// applicability depends on the future relative string velocity, so those checks
/// stay in the original solver at execution time. A physical refusal poisons the
/// voice; discard that callback, never retry partial output as a completed block.
/// The solver's one-way bridge/body approximation and allocation boundary remain.
/// This adapter does not claim general mixed-instrument rendering or real-time V&V.
pub struct ScheduledBowedRenderer {
    state: BowedStringState,
    controls: Vec<ScheduledBowControl>,
    next_control: usize,
    end_sample: u64,
}

impl ScheduledBowedRenderer {
    /// Compile one admitted bow track using the ACTUAL destination audio rate.
    ///
    /// `samples` is the additional, half-open performance horizon. `max_work`
    /// bounds tick count times worst-case source track/event visits; `max_events`
    /// bounds emitted changes, including the initial value for a nonempty run.
    /// Held values are coalesced bitwise. Compilation may allocate, playback does
    /// not allocate scheduler storage. The hosted solver retains its own behavior.
    ///
    /// # Errors
    /// Unknown/non-bow tracks, invalid clocks, clock overflow, exhausted budgets,
    /// allocation failure, source sampling refusal, or an already poisoned voice.
    pub fn new(
        state: BowedStringState,
        schedule: &GestureSchedule,
        track_id: &str,
        samples: u64,
        max_work: u64,
        max_events: usize,
    ) -> Result<Self, BowedScheduleError> {
        state.check_live().map_err(BowedScheduleError::Voice)?;
        let control_rate = schedule.control_rate_hz;
        let audio_rate = state.sample_rate_hz();
        if control_rate == 0 || audio_rate == 0 || control_rate > audio_rate {
            return Err(invalid("require 0 < control rate <= the voice audio rate"));
        }
        let track = schedule.tracks().iter().find(|track| track.id == track_id)
            .ok_or_else(|| BowedScheduleError::Gesture(GestureError::UnknownControlId {
                requested: track_id.to_string(),
            }))?;
        if !matches!(track.target, GestureTarget::BowStroke { .. }) {
            return Err(invalid("the selected track must target a bow stroke"));
        }
        let start = state.samples_rendered();
        let end_sample = start.checked_add(samples)
            .ok_or_else(|| invalid("performance horizon overflows the voice sample clock"))?;
        let audio = u128::from(audio_rate);
        let control = u128::from(control_rate);
        // ceil(k*a/c) < samples iff k*a <= (samples-1)*c.
        let ticks = if samples == 0 { 0 } else {
            u128::from(samples - 1) * control / audio + 1
        };
        let per_tick = schedule.tracks().len() as u128 + track.events.len() as u128 + 1;
        let required = ticks.checked_mul(per_tick)
            .ok_or_else(|| invalid("source sampling work bound overflows"))?;
        if required > u128::from(max_work) {
            return Err(BowedScheduleError::WorkBudget { required, allowed: max_work });
        }
        let mut controls = Vec::new();
        let mut previous = None;
        for tick in 0..ticks as u64 {
            let GestureValue::Bow { velocity_m_per_s, normal_force_n, station } =
                schedule.sample_value(track_id, tick).map_err(BowedScheduleError::Gesture)?
            else {
                return Err(invalid("bow track sampled a non-bow value"));
            };
            let bits = [velocity_m_per_s.to_bits(), normal_force_n.to_bits(), station.to_bits()];
            if previous == Some(bits) {
                continue;
            }
            if controls.len() == max_events {
                return Err(BowedScheduleError::EventBudget { allowed: max_events });
            }
            controls.try_reserve(1)
                .map_err(|_| invalid("cannot allocate compiled bow controls"))?;
            // Bounds above prove offset < samples and start+offset <= end_sample.
            let offset = (u128::from(tick) * audio).div_ceil(control) as u64;
            controls.push(ScheduledBowControl {
                sample: start + offset,
                velocity_m_s: velocity_m_per_s,
                normal_force_n,
                station,
            });
            previous = Some(bits);
        }
        Ok(Self { state, controls, next_control: 0, end_sample })
    }

    /// Read-only access to the current physical state and its lifetime clock.
    #[must_use]
    pub const fn state(&self) -> &BowedStringState { &self.state }

    /// Remaining samples in the admitted performance horizon.
    #[must_use]
    pub fn remaining_samples(&self) -> u64 { self.end_sample - self.state.samples_rendered() }

    /// Controls actually applied, with callback-independent sample timestamps.
    #[must_use]
    pub fn applied_controls(&self) -> &[ScheduledBowControl] {
        &self.controls[..self.next_control]
    }

    /// Admitted future controls, unchanged by cancellation.
    #[must_use]
    pub fn pending_controls(&self) -> &[ScheduledBowControl] {
        &self.controls[self.next_control..]
    }

    /// Recover the voice without resetting vibration. A poisoned voice stays
    /// poisoned even when the failure occurred while applying a scheduled input.
    #[must_use]
    pub fn into_state(self) -> BowedStringState { self.state }

    fn validate_window(&self, len: usize) -> Result<(), BowedScheduleError> {
        self.state.check_live().map_err(BowedScheduleError::Voice)?;
        let count = u64::try_from(len).map_err(|_| invalid("output exceeds the sample clock"))?;
        if count > self.remaining_samples() {
            return Err(invalid("output exceeds the compiled performance horizon"));
        }
        Ok(())
    }

    fn render_block<T>(
        &mut self,
        out: &mut [T],
        render: fn(&mut BowedStringState, &mut [T]) -> Result<(), BowedRunError>,
    ) -> Result<(), BowedScheduleError> {
        self.state.validate_request(out.len()).map_err(BowedScheduleError::Voice)?;
        self.validate_window(out.len())?;
        let end = self.state.samples_rendered() + out.len() as u64;
        let mut offset = 0;
        while offset < out.len() {
            let now = self.state.samples_rendered();
            if let Some(event) = self.controls.get(self.next_control).copied() {
                if event.sample == now {
                    if let Err(error) = self.state.set_bow(
                        event.velocity_m_s, event.normal_force_n, event.station,
                    ) {
                        self.state.poisoned = true;
                        return Err(BowedScheduleError::Voice(error));
                    }
                    self.next_control += 1;
                }
            }
            let stop = self.controls.get(self.next_control).map_or(end, |e| e.sample.min(end));
            // The compiler emits at most one complete bow input per sample.
            let count = (stop - now) as usize;
            if let Err(error) = render(&mut self.state, &mut out[offset..offset + count]) {
                self.state.poisoned = true;
                return Err(BowedScheduleError::Voice(error));
            }
            offset += count;
        }
        Ok(())
    }

    /// Render endpoint mechanics, splitting only at actual bow changes.
    /// Whole-block shape/horizon admission precedes any input or output change.
    pub fn block(&mut self, out: &mut [BowedSample]) -> Result<(), BowedScheduleError> {
        self.render_block(out, BowedStringState::block)
    }

    /// Render actual compact-body observer pressure [Pa]. A rigid termination
    /// refuses BEFORE consuming controls or changing output, never faking sound
    /// by relabeling a mechanical velocity channel.
    pub fn pressure_block(&mut self, out: &mut [f64]) -> Result<(), BowedScheduleError> {
        if !self.state.has_radiation() {
            return Err(invalid("pressure output requires an attached physical body"));
        }
        self.render_block(out, BowedStringState::pressure_block)
    }

    fn render_gated<T>(
        &mut self,
        gate: &CancelGate,
        out: &mut [T],
        block_len: usize,
        render: fn(&mut BowedStringState, &mut [T]) -> Result<(), BowedRunError>,
    ) -> Result<BowedRenderOutcome, BowedScheduleError> {
        self.state.validate_request(block_len).map_err(BowedScheduleError::Voice)?;
        self.validate_window(out.len())?;
        let mut samples = 0;
        for block in out.chunks_mut(block_len) {
            if gate.is_requested() {
                return Ok(BowedRenderOutcome::Cancelled { samples });
            }
            self.render_block(block, render)?;
            samples += block.len();
        }
        Ok(BowedRenderOutcome::Completed { samples })
    }

    /// Render mechanics with cancellation checked before every host callback.
    /// A cancelled suffix and all its controls remain untouched and resumable.
    pub fn render_under_gate(
        &mut self, gate: &CancelGate, out: &mut [BowedSample], block_len: usize,
    ) -> Result<BowedRenderOutcome, BowedScheduleError> {
        self.render_gated(gate, out, block_len, BowedStringState::block)
    }

    /// Render physical observer pressure with the same cancellation/resume rule.
    pub fn pressure_under_gate(
        &mut self, gate: &CancelGate, out: &mut [f64], block_len: usize,
    ) -> Result<BowedRenderOutcome, BowedScheduleError> {
        if !self.state.has_radiation() {
            return Err(invalid("pressure output requires an attached physical body"));
        }
        self.render_gated(gate, out, block_len, BowedStringState::pressure_block)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bowed_string::{BowGesture, BowedRunConfig, BowedStringCard, FrictionIsland, Termination};
    use crate::stribeck_friction::StribeckFriction;
    use fs_scenario::gesture::{GestureEvent, GestureTrack};

    fn config() -> BowedRunConfig {
        BowedRunConfig {
            card: BowedStringCard {
                length_m: 0.65, tension_n: 60.0, linear_density_kg_m: 6.0e-4,
                bending_stiffness_n_m2: 0.0, viscous_bending_n_m2_s: 0.0,
                mode_count: 16,
                zetas: (0..16).map(|k| 1.0e-3 * (1.0 + 0.55 * f64::from(k))).collect(),
                sample_rate_hz: 48_000,
            },
            island: FrictionIsland::Stribeck(StribeckFriction::try_new(0.8, 0.4, 0.04).unwrap()),
            gesture: BowGesture::admit(0.45, 3.9, 0.11).unwrap(),
            steps: 1024, subsamples: 16, termination: Termination::Rigid, listener_m: 1.0,
        }
    }

    fn bow(v: f64, f: f64, x: f64) -> GestureValue {
        GestureValue::Bow { velocity_m_per_s: v, normal_force_n: f, station: x }
    }

    fn performance() -> GestureSchedule {
        GestureSchedule::try_new(700, vec![GestureTrack {
            id: "bow".into(), target: GestureTarget::BowStroke { string: 0 },
            initial: bow(0.45, 3.9, 0.11),
            events: vec![
                GestureEvent { time_s: 0.0, transition_s: 5.0 / 700.0, value: bow(-0.3, 2.0, 0.2) },
                GestureEvent { time_s: 3.0 / 700.0, transition_s: 3.0 / 700.0, value: bow(0.25, 1.0, 0.15) },
                GestureEvent { time_s: 7.0 / 700.0, transition_s: 0.0, value: bow(-0.1, 0.0, 0.15) },
                GestureEvent { time_s: 9.0 / 700.0, transition_s: 1.0 / 700.0, value: bow(0.35, 2.0, 0.12) },
                GestureEvent { time_s: 12.0 / 700.0, transition_s: 0.0, value: bow(0.0, 0.0, 0.12) },
            ],
        }]).unwrap()
    }

    fn bits(s: BowedSample) -> [u64; 4] {
        [s.bow_point_velocity_m_s.to_bits(), s.relative_velocity_m_s.to_bits(),
            s.bridge_force_n.to_bits(), s.total_modal_energy_j.to_bits()]
    }

    #[test]
    fn nondivisor_bow_performance_matches_samplewise_physics_across_partitions() {
        let source = performance();
        let decoded = GestureSchedule::from_canonical_bytes(&source.to_canonical_bytes()).unwrap();
        let mut direct = BowedStringState::new(&config(), 1024).unwrap();
        for _ in 0..19 { direct.step().unwrap(); }
        let initial_energy = direct.total_modal_energy_j().to_bits();
        let mut expected = Vec::new();
        let mut previous = None;
        for sample in 0..1024_u64 {
            // Inverse sample clock, independent of the compiler's ceil mapping.
            let tick = sample * 700 / 48_000;
            let GestureValue::Bow { velocity_m_per_s, normal_force_n, station } =
                source.sample_value("bow", tick).unwrap() else { panic!("bow value") };
            let key = [velocity_m_per_s.to_bits(), normal_force_n.to_bits(), station.to_bits()];
            if previous != Some(key) {
                direct.set_bow(velocity_m_per_s, normal_force_n, station).unwrap();
                previous = Some(key);
            }
            expected.push(bits(direct.step().unwrap()));
        }
        assert!(expected.iter().any(|s| f64::from_bits(s[3]) > 1e-10));
        for schedule in [&source, &decoded] {
            for partition in [1, 37, 256, 1024] {
                let mut state = BowedStringState::new(&config(), 1024).unwrap();
                for _ in 0..19 { state.step().unwrap(); }
                let mut render = ScheduledBowedRenderer::new(state, schedule, "bow", 1024, 10000, 64).unwrap();
                assert_eq!(render.state().total_modal_energy_j().to_bits(), initial_energy);
                assert_eq!(render.state().samples_rendered(), 19);
                assert_eq!(render.pending_controls()[1].sample, 19 + 69);
                let mut out = vec![BowedSample::default(); 1024];
                for block in out.chunks_mut(partition) { render.block(block).unwrap(); }
                assert_eq!(out.into_iter().map(bits).collect::<Vec<_>>(), expected);
                assert_eq!(render.remaining_samples(), 0);
                assert_eq!(render.state().samples_rendered(), 1043);
                assert!(render.pending_controls().is_empty());
                assert_eq!(render.state().total_modal_energy_j().to_bits(), direct.total_modal_energy_j().to_bits());
            }
        }
    }

    #[test]
    fn physical_pressure_matches_direct_observer_and_survives_bow_release() {
        use crate::thin_plate::CompactBody;
        use fs_material::gas::{GasSpec, GasState};
        use fs_scenario::RadiatingPlate;

        let mut cfg = config();
        cfg.termination = Termination::PlateOnePort {
            body: Box::new(CompactBody::from_radiator(RadiatingPlate {
                area_m2: 3.0e-3, mass_kg: 0.15, frequency_hz: 280.0, damping_ratio: 0.02,
            }).unwrap()),
            ambient: GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap(),
        };
        let source = performance();
        let mut direct = BowedStringState::new(&cfg, 1024).unwrap();
        let mut expected = Vec::new();
        let mut previous = None;
        for sample in 0..1024_u64 {
            let GestureValue::Bow { velocity_m_per_s, normal_force_n, station } =
                source.sample_value("bow", sample * 700 / 48_000).unwrap()
                else { panic!("bow value") };
            let key = [velocity_m_per_s.to_bits(), normal_force_n.to_bits(), station.to_bits()];
            if previous != Some(key) {
                direct.set_bow(velocity_m_per_s, normal_force_n, station).unwrap();
                previous = Some(key);
            }
            expected.push(direct.step().unwrap().radiated_pressure_pa.unwrap().to_bits());
        }
        // Final release begins at ceil(12 * 48000 / 700) = 823, not a buffer edge.
        assert!(expected[823..].iter().any(|&p| f64::from_bits(p).abs() > 1e-12));
        for partition in [1, 37, 256, 1024] {
            let mut render = ScheduledBowedRenderer::new(
                BowedStringState::new(&cfg, 1024).unwrap(), &source, "bow", 1024, 10000, 64,
            ).unwrap();
            let mut out = vec![0.0; 1024];
            assert_eq!(render.pressure_under_gate(&CancelGate::new(), &mut out, partition).unwrap(),
                BowedRenderOutcome::Completed { samples: 1024 });
            assert_eq!(out.into_iter().map(f64::to_bits).collect::<Vec<_>>(), expected);
        }
    }

    #[test]
    fn boundary_events_and_cancelled_callbacks_remain_pending_until_resumed() {
        let make = || ScheduledBowedRenderer::new(
            BowedStringState::new(&config(), 1024).unwrap(), &performance(), "bow", 1024, 10000, 64,
        ).unwrap();
        let mut render = make();
        let mut prefix = vec![BowedSample::default(); 69];
        render.block(&mut prefix).unwrap();
        assert_eq!(render.pending_controls()[0].sample, 69);
        let pending = render.pending_controls().to_vec();
        let applied = render.applied_controls().to_vec();
        let energy = render.state().total_modal_energy_j().to_bits();
        let sentinel = BowedSample { bridge_force_n: 12345.0, ..BowedSample::default() };
        let mut suffix = vec![sentinel; 955];
        let gate = CancelGate::new();
        gate.request();
        assert_eq!(render.render_under_gate(&gate, &mut suffix, 37).unwrap(),
            BowedRenderOutcome::Cancelled { samples: 0 });
        assert_eq!(suffix, vec![sentinel; 955]);
        assert_eq!(render.pending_controls(), pending);
        assert_eq!(render.applied_controls(), applied);
        assert_eq!(render.state().total_modal_energy_j().to_bits(), energy);
        assert_eq!(render.state().samples_rendered(), 69);
        assert_eq!(render.render_under_gate(&CancelGate::new(), &mut suffix, 37).unwrap(),
            BowedRenderOutcome::Completed { samples: 955 });
        prefix.extend(suffix);
        let mut reference = vec![BowedSample::default(); 1024];
        make().block(&mut reference).unwrap();
        assert_eq!(prefix.into_iter().map(bits).collect::<Vec<_>>(), reference.into_iter().map(bits).collect::<Vec<_>>());
    }

    #[test]
    fn invalid_blocks_pressure_and_horizons_do_not_consume_controls() {
        let mut render = ScheduledBowedRenderer::new(
            BowedStringState::new(&config(), 8).unwrap(), &performance(), "bow", 10, 1000, 64,
        ).unwrap();
        let pending = render.pending_controls().to_vec();
        let sentinel = BowedSample { bridge_force_n: 123.0, ..BowedSample::default() };
        let mut out = [sentinel; 11];
        assert!(render.block(&mut []).is_err());
        assert!(render.block(&mut out[..9]).is_err());
        assert!(render.render_under_gate(&CancelGate::new(), &mut out, 8).is_err());
        assert_eq!(out, [sentinel; 11]);
        let mut pressure = [123.0; 8];
        assert!(render.pressure_block(&mut pressure).is_err());
        assert!(render.pressure_under_gate(&CancelGate::new(), &mut pressure, 8).is_err());
        assert_eq!(pressure, [123.0; 8]);
        assert_eq!(render.pending_controls(), pending);
        assert_eq!(render.state().samples_rendered(), 0);
        render.block(&mut out[..8]).unwrap();
        assert!(render.block(&mut out[..3]).is_err());
        assert_eq!(render.remaining_samples(), 2);
    }

    #[test]
    fn source_clocks_work_storage_and_poisoned_voices_refuse_at_binding() {
        let bind = |s: &GestureSchedule, id: &str, samples, work, events| {
            ScheduledBowedRenderer::new(BowedStringState::new(&config(), 8).unwrap(), s, id, samples, work, events)
        };
        let mut source = performance();
        assert!(matches!(bind(&source, "absent", 0, 0, 0), Err(BowedScheduleError::Gesture(_))));
        assert!(matches!(bind(&source, "bow", 1024, 104, 64),
            Err(BowedScheduleError::WorkBudget { required: 105, allowed: 104 })));
        assert!(matches!(bind(&source, "bow", 1024, 10000, 0), Err(BowedScheduleError::EventBudget { allowed: 0 })));
        assert!(bind(&source, "bow", u64::MAX, 10000, 64).is_err());
        assert_eq!(bind(&source, "bow", 0, 0, 0).unwrap().remaining_samples(), 0);
        for rate in [0, 48001] {
            source.control_rate_hz = rate;
            assert!(matches!(bind(&source, "bow", 0, 0, 0), Err(BowedScheduleError::Invalid { .. })));
        }
        let pressure = GestureSchedule::try_new(700, vec![GestureTrack {
            id: "pressure".into(), target: GestureTarget::BlowingPressure,
            initial: GestureValue::PressurePa(0.0), events: Vec::new(),
        }]).unwrap();
        assert!(matches!(bind(&pressure, "pressure", 0, 0, 0), Err(BowedScheduleError::Invalid { .. })));
        let mut state = BowedStringState::new(&config(), 8).unwrap();
        state.poisoned = true;
        assert!(matches!(ScheduledBowedRenderer::new(state, &performance(), "bow", 8, 1000, 64),
            Err(BowedScheduleError::Voice(BowedRunError::Poisoned))));
    }

    #[test]
    fn held_controls_coalesce_but_the_work_budget_counts_all_ticks() {
        let source = GestureSchedule::try_new(700, vec![GestureTrack {
            id: "held".into(), target: GestureTarget::BowStroke { string: 0 },
            initial: bow(0.45, 3.9, 0.11), events: Vec::new(),
        }]).unwrap();
        let make = || BowedStringState::new(&config(), 1024).unwrap();
        assert!(matches!(ScheduledBowedRenderer::new(make(), &source, "held", 1024, 29, 1),
            Err(BowedScheduleError::WorkBudget { required: 30, allowed: 29 })));
        let render = ScheduledBowedRenderer::new(make(), &source, "held", 1024, 30, 1).unwrap();
        assert_eq!(render.pending_controls().len(), 1);
        assert_eq!(render.pending_controls()[0].sample, 0);
    }

    #[test]
    fn failed_physical_callback_stays_poisoned_after_recovering_the_voice() {
        let mut cfg = config();
        cfg.island = FrictionIsland::ViscousOnly { viscous_n_s_per_m: f64::MAX };
        let mut render = ScheduledBowedRenderer::new(
            BowedStringState::new(&cfg, 8).unwrap(), &performance(), "bow", 8, 1000, 64,
        ).unwrap();
        assert!(render.block(&mut [BowedSample::default(); 8]).is_err());
        assert!(matches!(render.block(&mut [BowedSample::default(); 8]),
            Err(BowedScheduleError::Voice(BowedRunError::Poisoned))));
        assert!(matches!(render.into_state().step(), Err(BowedRunError::Poisoned)));
    }
}
