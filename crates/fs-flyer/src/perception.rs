//! Deterministic pilot-perception service (bead wf-root-guzez.5.16.1,
//! E4.6c-i). Plan §5.1.4 + Round-3 revision: perception is a
//! DETERMINISTIC fs-flyer service at the 120 Hz baseline with a frozen
//! causal order — and it is NEVER renderer-fed, which here is
//! STRUCTURAL: the step signature admits only the physics cue vector
//! and the tick; there is no render input to leak (V-16b's
//! cross-backend invariance follows from the type, and the browser
//! matrix execution of that claim lives in the E0.6 lane).
//!
//! Per cue: integer-tick delay ring → first-order filter (EXACT
//! exponential update) → gain, plus a deterministic remnant drawn from
//! a philox counter-addressed stream whose draw index is BOUND TO THE
//! TICK (fs-rand doctrine) — checkpoint/resume reproduces the identical
//! remnant sequence.

use crate::Refusal;
use fs_math::det;
use fs_rand::{Stream, StreamKey};

/// The perception baseline rate [Hz] (plan Round-5 decision).
pub const PERCEPTION_HZ: f64 = 120.0;

/// Delay-ring cap [ticks] (2 s at 120 Hz — an absurd-input guard).
pub const MAX_DELAY_TICKS: usize = 240;

/// Number of cues (frozen order: pitch attitude, pitch rate, heave
/// acceleration, roll attitude, roll rate, yaw rate).
pub const N_CUES: usize = 6;

/// Draws consumed per tick (one normal per cue — the tick→index map).
pub const DRAWS_PER_TICK: u64 = N_CUES as u64;

/// One cue channel's spec.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CueSpec {
    /// Perceptual delay [ticks at 120 Hz].
    pub delay_ticks: usize,
    /// First-order filter time constant [s] (> 0).
    pub filter_tau_s: f64,
    /// Perceptual gain.
    pub gain: f64,
    /// Remnant one-sigma on this cue (post-filter units).
    pub remnant_sigma: f64,
}

/// The perception model spec (enters PerceptionModelId).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PerceptionModelSpec {
    /// Per-cue channels (frozen order).
    pub cues: [CueSpec; N_CUES],
    /// Remnant stream identity.
    pub stream_key: StreamKey,
}

/// The registered v1 spec: vestibular/visual delay classes (~150 ms
/// attitude, ~80 ms rates), light smoothing, modest remnant.
#[must_use]
pub fn perception_v1(seed: u64) -> PerceptionModelSpec {
    let att = CueSpec {
        delay_ticks: 18, // 150 ms
        filter_tau_s: 0.06,
        gain: 1.0,
        remnant_sigma: 0.002,
    };
    let rate = CueSpec {
        delay_ticks: 10, // ~83 ms
        filter_tau_s: 0.04,
        gain: 1.0,
        remnant_sigma: 0.004,
    };
    let accel = CueSpec {
        delay_ticks: 6,
        filter_tau_s: 0.08,
        gain: 1.0,
        remnant_sigma: 0.02,
    };
    PerceptionModelSpec {
        cues: [att, rate, accel, att, rate, rate],
        stream_key: StreamKey {
            seed,
            kernel: 0x4643_5031, // "FCP1" — perception kernel id
            tile: 0,
        },
    }
}

/// Perception state (rings + filters + tick).
#[derive(Clone, Debug, PartialEq)]
pub struct PerceptionState {
    rings: Vec<Vec<f64>>,
    filters: [f64; N_CUES],
    tick: u64,
}

/// One tick's perceived cues.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PerceivedCues {
    /// Filtered, delayed, gained cue values + remnant (frozen order).
    pub values: [f64; N_CUES],
    /// The tick these belong to.
    pub tick: u64,
}

impl PerceptionModelSpec {
    /// Admit the spec.
    ///
    /// # Errors
    /// `perception-spec-invalid` (delay above [`MAX_DELAY_TICKS`] —
    /// tested at cap AND cap+1; non-positive tau; non-finite gain;
    /// negative remnant).
    pub fn admit(&self) -> Result<(), Refusal> {
        for (i, c) in self.cues.iter().enumerate() {
            let ok = c.delay_ticks <= MAX_DELAY_TICKS
                && c.filter_tau_s.is_finite()
                && c.filter_tau_s > 0.0
                && c.gain.is_finite()
                && c.remnant_sigma.is_finite()
                && c.remnant_sigma >= 0.0;
            if !ok {
                return Err(Refusal {
                    code: "perception-spec-invalid",
                    message: format!("cue {i}: {c:?}"),
                    ranked_repairs: vec![format!(
                        "delay <= {MAX_DELAY_TICKS}; tau > 0; finite gain; remnant >= 0"
                    )],
                });
            }
        }
        Ok(())
    }

    /// Fresh state at tick 0 (rings zero-filled: trim-relative cues).
    ///
    /// # Errors
    /// Admission refusals.
    pub fn init(&self) -> Result<PerceptionState, Refusal> {
        self.admit()?;
        Ok(PerceptionState {
            rings: self
                .cues
                .iter()
                .map(|c| vec![0.0; c.delay_ticks + 1])
                .collect(),
            filters: [0.0; N_CUES],
            tick: 0,
        })
    }

    /// One 120 Hz step: push raw cues, read the delayed samples, filter
    /// exactly, add the tick-addressed remnant. The ONLY inputs are the
    /// physics cue vector and the state — renderer feeding is
    /// impossible by construction.
    ///
    /// # Errors
    /// `perception-input-invalid` (non-finite cue).
    pub fn step(
        &self,
        state: &mut PerceptionState,
        raw: [f64; N_CUES],
    ) -> Result<PerceivedCues, Refusal> {
        if raw.iter().any(|v| !v.is_finite()) {
            return Err(Refusal {
                code: "perception-input-invalid",
                message: format!("{raw:?}"),
                ranked_repairs: vec!["check the physics cue assembly".into()],
            });
        }
        let dt = 1.0 / PERCEPTION_HZ;
        // Tick-addressed remnant: draw index = tick * DRAWS_PER_TICK,
        // one normal per cue, in cue order (frozen causal order).
        let mut stream = Stream::resume(fs_rand::StreamCheckpoint {
            checkpoint_version: fs_rand::STREAM_CHECKPOINT_VERSION,
            stream_semantics_version: fs_rand::STREAM_SEMANTICS_VERSION,
            key: self.stream_key,
            index: state.tick * DRAWS_PER_TICK,
        })
        .map_err(|e| Refusal {
            code: "perception-stream-invalid",
            message: format!("{e:?}"),
            ranked_repairs: vec![
                "stream checkpoint arithmetic overflow — tick out of range".into(),
            ],
        })?;
        let mut values = [0.0f64; N_CUES];
        for i in 0..N_CUES {
            let c = &self.cues[i];
            let d = c.delay_ticks;
            let t = state.tick as usize;
            // Exactly-D-tick delay: the sample written at tick t−D lives
            // in slot (t−D) mod (D+1) = (t+1) mod (D+1); read it BEFORE
            // writing this tick's sample into slot t mod (D+1). D = 0 is
            // the passthrough case (read what we are about to write).
            let delayed = if d == 0 {
                raw[i]
            } else {
                state.rings[i][(t + 1) % (d + 1)]
            };
            state.rings[i][t % (d + 1)] = raw[i];
            // Exact first-order filter toward the delayed sample.
            let a = det::exp(-dt / c.filter_tau_s);
            state.filters[i] = delayed + (state.filters[i] - delayed) * a;
            let remnant = c.remnant_sigma * stream.next_normal();
            values[i] = c.gain * state.filters[i] + remnant;
        }
        let out = PerceivedCues {
            values,
            tick: state.tick,
        };
        state.tick += 1;
        Ok(out)
    }

    /// Checkpoint the state (rings + filters + tick are the complete
    /// causal state; the remnant stream is reconstructed from the tick).
    #[must_use]
    pub fn checkpoint(&self, state: &PerceptionState) -> PerceptionState {
        state.clone()
    }
}
