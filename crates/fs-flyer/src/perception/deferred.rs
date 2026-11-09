//! Acquisition after control is causal only for a strictly delayed channel.
//! The current output has already read an older slot. Completing the current
//! slot before publishing a successful tick adds neither latency nor history.

use super::{N_CUES, PerceivedCues, PerceptionModelSpec, PerceptionState, Refusal};

const HEAVE: usize = 2;

// A move-only token binds completion to one channel, tick and exact ring slot.
pub(crate) struct DeferredHeaveSample {
    tick: u64,
    channel: usize,
    slot: usize,
    ring_len: usize,
}

fn invalid() -> Refusal {
    Refusal {
        code: "perception-deferred-input-invalid",
        message: "heave acquisition needs a finite sample at its reserved delayed tick/slot".into(),
        ranked_repairs: vec![
            "use a positive heave delay and complete the matching accepted physics tick".into(),
        ],
    }
}

impl PerceptionModelSpec {
    fn reserve_heave_sample(
        &self,
        state: &PerceptionState,
        sample_tick: u64,
    ) -> Result<DeferredHeaveSample, Refusal> {
        let delay = self.cues[HEAVE].delay_ticks;
        if delay == 0
            || delay > super::MAX_DELAY_TICKS
            || state.tick != sample_tick
            || state.rings.get(HEAVE).map(Vec::len) != Some(delay + 1)
        {
            return Err(invalid());
        }
        Ok(DeferredHeaveSample {
            tick: sample_tick,
            channel: HEAVE,
            slot: sample_tick as usize % (delay + 1),
            ring_len: delay + 1,
        })
    }

    pub(crate) fn admit_deferred_heave(
        &self,
        state: &PerceptionState,
        sample_tick: u64,
    ) -> Result<(), Refusal> {
        self.reserve_heave_sample(state, sample_tick).map(|_| ())
    }

    pub(crate) fn step_with_deferred_heave(
        &self,
        state: &mut PerceptionState,
        sample_tick: u64,
        mut raw: [f64; N_CUES],
    ) -> Result<(PerceivedCues, DeferredHeaveSample), Refusal> {
        let pending = self.reserve_heave_sample(state, sample_tick)?;
        raw[HEAVE] = 0.0;
        Ok((self.step(state, raw)?, pending))
    }

    pub(crate) fn complete_heave_sample(
        &self,
        state: &mut PerceptionState,
        pending: DeferredHeaveSample,
        raw: f64,
    ) -> Result<(), Refusal> {
        let delay = self.cues[HEAVE].delay_ticks;
        if !raw.is_finite()
            || delay == 0
            || delay > super::MAX_DELAY_TICKS
            || pending.channel != HEAVE
            || pending.tick.checked_add(1) != Some(state.tick)
            || pending.ring_len != delay + 1
            || pending.slot != pending.tick as usize % (delay + 1)
            || state.rings.get(HEAVE).map(Vec::len) != Some(pending.ring_len)
        {
            return Err(invalid());
        }
        state.rings[HEAVE][pending.slot] = raw;
        Ok(())
    }

    pub(crate) fn latest_heave_input(&self, state: &PerceptionState) -> Option<(u64, f64)> {
        let tick = state.tick.checked_sub(1)?;
        let ring = state.rings.get(HEAVE)?;
        if ring.len() != self.cues[HEAVE].delay_ticks + 1 {
            return None;
        }
        Some((tick, *ring.get(tick as usize % ring.len())?))
    }

    /// Reconstruct the last emitted cues from the retained filters and the
    /// tick-addressed remnant. No clock, ring, filter or random stream advances.
    /// The state must belong to this model; fresh state has no output yet.
    ///
    /// # Errors
    /// The same remnant-stream refusal as [`Self::step`].
    pub fn last_output(&self, state: &PerceptionState) -> Result<Option<PerceivedCues>, Refusal> {
        let Some(tick) = state.tick.checked_sub(1) else {
            return Ok(None);
        };
        let mut stream = self.remnant_stream(tick)?;
        let values = core::array::from_fn(|i| {
            let remnant = self.cues[i].remnant_sigma * stream.next_normal();
            let gained = self.cues[i].gain * state.filters[i];
            gained + remnant
        });
        Ok(Some(PerceivedCues { values, tick }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::perception::{PERCEPTION_HZ, perception_v1};

    #[test]
    fn acquisition_after_control_has_exactly_six_ticks_of_delay() {
        let mut model = perception_v1(41);
        for cue in &mut model.cues {
            cue.remnant_sigma = 0.0;
        }
        let mut state = model.init().unwrap();
        let mut direct = state.clone();
        for tick in 0..14 {
            let sample = if tick == 0 { 3.0 } else { 0.0 };
            let (output, pending) = model
                .step_with_deferred_heave(&mut state, tick, [0.0; N_CUES])
                .unwrap();
            if tick < 6 {
                assert_eq!(output.values[HEAVE], 0.0);
            } else if tick == 6 {
                let expected = 3.0 * (1.0 - (-1.0 / (PERCEPTION_HZ * 0.08)).exp());
                assert!((output.values[HEAVE] - expected).abs() < 2e-14);
            }
            model
                .complete_heave_sample(&mut state, pending, sample)
                .unwrap();
            let mut raw = [0.0; N_CUES];
            raw[HEAVE] = sample;
            assert_eq!(output, model.step(&mut direct, raw).unwrap());
            assert_eq!(
                state, direct,
                "acquisition order must not alter the delay/filter"
            );
            let before = state.to_bytes();
            assert_eq!(model.last_output(&state).unwrap(), Some(output));
            assert_eq!(before, state.to_bytes());
        }
    }

    #[test]
    fn read_only_output_preserves_the_real_remnant_sequence() {
        let model = perception_v1(47);
        let mut state = model.init().unwrap();
        for _ in 0..12 {
            let output = model.step(&mut state, [0.3; N_CUES]).unwrap();
            let before = state.to_bytes();
            assert_eq!(model.last_output(&state).unwrap(), Some(output));
            assert_eq!(model.last_output(&state).unwrap(), Some(output));
            assert_eq!(before, state.to_bytes());
        }
    }

    #[test]
    fn completion_refuses_wrong_channel_tick_slot_shape_and_nonfinite_input() {
        for fault in 0..6 {
            let model = perception_v1(43);
            let mut state = model.init().unwrap();
            let (_, mut pending) = model
                .step_with_deferred_heave(&mut state, 0, [0.0; N_CUES])
                .unwrap();
            let raw = match fault {
                0 => {
                    pending.channel = 1;
                    2.0
                }
                1 => {
                    pending.tick += 1;
                    2.0
                }
                2 => {
                    pending.slot += 1;
                    2.0
                }
                3 => {
                    pending.ring_len += 1;
                    2.0
                }
                4 => f64::NAN,
                _ => f64::INFINITY,
            };
            let before = state.clone();
            assert!(
                model
                    .complete_heave_sample(&mut state, pending, raw)
                    .is_err()
            );
            assert_eq!(state, before);
        }
        let mut model = perception_v1(43);
        model.cues[HEAVE].delay_ticks = 0;
        let mut state = model.init().unwrap();
        let before = state.clone();
        assert!(
            model
                .step_with_deferred_heave(&mut state, 0, [0.0; N_CUES])
                .is_err()
        );
        assert_eq!(state, before, "zero-delay feedback cannot be deferred");
    }
}
