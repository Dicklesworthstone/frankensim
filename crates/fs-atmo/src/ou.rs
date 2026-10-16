//! Exact-discrete OU amplitude evolution + StationaryOuPathV1 (bead
//! wf-root-guzez.4.6.1, E3.3b-i). Plan §5.4:
//!
//!   a_{k,n+1} = ρ_k·a_{k,n} + σ_k·ξ_{k,n}
//!
//! with ρ = e^(−Δt/τ) and σ_innov² = σ_st²·(1 − ρ²) — the EXACT discrete
//! transition of the continuous OU process, so the marginal at every tick
//! is exactly the stationary law and step-size refinement is not an
//! approximation axis.
//!
//! Innovations are philox COUNTER-ADDRESSED: ξ(mode, tick) is a pure
//! function of (seed, OU_KERNEL, mode, tick) via `StreamCheckpoint::current`
//! at index 2·(tick − anchor). Consequences, each battery-proven:
//! - the state is SEQUENTIAL and checkpointable, and resuming from a
//!   checkpoint reproduces the continuous path BITWISE (the innovations a
//!   resumed path consumes are the same pure function);
//! - `StationaryOuPathV1`: one exact stationary draw at the FIXED anchor
//!   tick −3840 (−32 s at 120 Hz, plan Round-4); every pre-roll window is
//!   a suffix of ONE physical path, so CRN comparisons survive changes in
//!   when sampling starts.
//!
//! `sample(x, tick)` NEVER hashes a fresh amplitude per query: queries read
//! the evolved sequential state; only `advance_to` moves it.

use crate::{Refusal, TICK_HZ, refuse};
use fs_math::det;
use fs_rand::{Stream, StreamCheckpoint, StreamKey};

/// Registered fs-rand kernel id for OU innovations ("OU__").
pub const OU_KERNEL: u32 = 0x4F555F5F;
/// The fixed StationaryOuPathV1 anchor tick (−32 s at 120 Hz).
pub const STATIONARY_ANCHOR_TICK: i64 = -3840;
/// Mode cap (matches the field cap; refusals at cap AND cap+1).
pub const MAX_OU_MODES: usize = 256;
/// Cap on a single `advance_to` span [ticks] (~1 hour at 120 Hz).
pub const MAX_ADVANCE_TICKS: i64 = 432_000;

/// Per-mode OU parameters (derived from a correlation time).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OuMode {
    /// One-tick autoregression factor ρ = e^(−Δt/τ).
    pub rho: f64,
    /// Innovation std σ_innov = σ_st·√(1 − ρ²).
    pub sigma_innov: f64,
    /// Stationary std σ_st.
    pub sigma_stationary: f64,
}

impl OuMode {
    /// Exact-discrete parameters for stationary std `sigma_st` and
    /// correlation time `tau_s` at the fixed 120 Hz tick.
    ///
    /// # Errors
    /// `ou-params-invalid` (non-finite, σ < 0, or τ ≤ 0).
    pub fn from_correlation_time(sigma_st: f64, tau_s: f64) -> Result<OuMode, Refusal> {
        if !(sigma_st.is_finite() && tau_s.is_finite()) || sigma_st < 0.0 || tau_s <= 0.0 {
            return Err(refuse(
                "ou-params-invalid",
                format!("sigma_st {sigma_st}, tau {tau_s}"),
                "sigma_st >= 0 and tau > 0 (integral-scale / mean-speed gives tau)",
            ));
        }
        let rho = det::exp(-1.0 / (TICK_HZ * tau_s));
        let sigma_innov = sigma_st * det::sqrt(1.0 - rho * rho);
        Ok(OuMode {
            rho,
            sigma_innov,
            sigma_stationary: sigma_st,
        })
    }
}

/// The sequential OU path state: per-mode amplitudes at `tick`.
#[derive(Clone, Debug, PartialEq)]
pub struct StationaryOuPath {
    seed: u64,
    modes: Vec<OuMode>,
    /// Current tick (state = amplitudes AT this tick).
    tick: i64,
    /// Per-mode amplitudes at `tick`.
    amplitudes: Vec<f64>,
}

/// A checkpoint of the OU path (transport data; part of CheckpointStateV1's
/// atmosphere block downstream).
#[derive(Clone, Debug, PartialEq)]
pub struct OuCheckpoint {
    /// Seed the path was built with.
    pub seed: u64,
    /// Tick the amplitudes belong to.
    pub tick: i64,
    /// Per-mode amplitudes at `tick`.
    pub amplitudes: Vec<f64>,
}

/// The counter-addressed innovation ξ(mode, tick): a pure function.
/// Index 2·(tick − anchor) — each normal consumes exactly two raw draws
/// (Box–Muller), so consecutive ticks use disjoint counter pairs.
fn innovation(seed: u64, mode: u32, tick: i64) -> f64 {
    let offset = u64::try_from(tick - STATIONARY_ANCHOR_TICK)
        .expect("tick >= anchor is enforced by construction");
    let key = StreamKey {
        seed,
        kernel: OU_KERNEL,
        tile: mode,
    };
    let mut s = Stream::resume(StreamCheckpoint::current(key, 2 * offset))
        .expect("current-version checkpoint always resumes");
    s.next_normal()
}

impl StationaryOuPath {
    /// Open the path at the FIXED anchor with the exact stationary draw:
    /// a(anchor) = σ_st·ξ(mode, anchor).
    ///
    /// # Errors
    /// `mode-count-invalid` (0 or above [`MAX_OU_MODES`], tested at cap
    /// AND cap+1).
    pub fn stationary_at_anchor(
        seed: u64,
        modes: Vec<OuMode>,
    ) -> Result<StationaryOuPath, Refusal> {
        if modes.is_empty() || modes.len() > MAX_OU_MODES {
            return Err(refuse(
                "mode-count-invalid",
                format!("{} modes outside [1, {MAX_OU_MODES}]", modes.len()),
                "match the turbulence field's mode count",
            ));
        }
        let amplitudes = modes
            .iter()
            .enumerate()
            .map(|(i, m)| m.sigma_stationary * innovation(seed, i as u32, STATIONARY_ANCHOR_TICK))
            .collect();
        Ok(StationaryOuPath {
            seed,
            modes,
            tick: STATIONARY_ANCHOR_TICK,
            amplitudes,
        })
    }

    /// Current tick.
    #[must_use]
    pub fn tick(&self) -> i64 {
        self.tick
    }

    /// Amplitudes at the current tick (one per mode).
    #[must_use]
    pub fn amplitudes(&self) -> &[f64] {
        &self.amplitudes
    }

    /// Advance the sequential state to `target_tick` (exact transition per
    /// tick; queries between ticks do not exist — the tick IS the clock).
    ///
    /// # Errors
    /// `ou-advance-backwards` (reduced time never runs backwards; resume
    /// from a checkpoint instead); `ou-advance-span-exceeded` above
    /// [`MAX_ADVANCE_TICKS`] (tested at cap AND cap+1).
    pub fn advance_to(&mut self, target_tick: i64) -> Result<(), Refusal> {
        if target_tick < self.tick {
            return Err(refuse(
                "ou-advance-backwards",
                format!("target {target_tick} behind current {}", self.tick),
                "resume an earlier checkpoint; the path never rewinds",
            ));
        }
        if target_tick - self.tick > MAX_ADVANCE_TICKS {
            return Err(refuse(
                "ou-advance-span-exceeded",
                format!(
                    "span {} exceeds {MAX_ADVANCE_TICKS}",
                    target_tick - self.tick
                ),
                "advance in bounded chunks (checkpointing between)",
            ));
        }
        while self.tick < target_tick {
            let next = self.tick + 1;
            for (i, (a, m)) in self.amplitudes.iter_mut().zip(&self.modes).enumerate() {
                let xi = innovation(self.seed, i as u32, next);
                *a = m.rho * *a + m.sigma_innov * xi;
            }
            self.tick = next;
        }
        Ok(())
    }

    /// Snapshot the sequential state.
    #[must_use]
    pub fn checkpoint(&self) -> OuCheckpoint {
        OuCheckpoint {
            seed: self.seed,
            tick: self.tick,
            amplitudes: self.amplitudes.clone(),
        }
    }

    /// Resume from a checkpoint (modes are re-supplied — they are model
    /// parameters, not path state).
    ///
    /// # Errors
    /// `ou-checkpoint-invalid` (mode-count mismatch, non-finite amplitude,
    /// or a tick before the anchor).
    pub fn resume(
        checkpoint: OuCheckpoint,
        modes: Vec<OuMode>,
    ) -> Result<StationaryOuPath, Refusal> {
        if checkpoint.amplitudes.len() != modes.len()
            || checkpoint.tick < STATIONARY_ANCHOR_TICK
            || !checkpoint.amplitudes.iter().all(|a| a.is_finite())
        {
            return Err(refuse(
                "ou-checkpoint-invalid",
                format!(
                    "{} amplitudes vs {} modes at tick {}",
                    checkpoint.amplitudes.len(),
                    modes.len(),
                    checkpoint.tick
                ),
                "a checkpoint pairs with the mode set it was written under",
            ));
        }
        Ok(StationaryOuPath {
            seed: checkpoint.seed,
            modes,
            tick: checkpoint.tick,
            amplitudes: checkpoint.amplitudes,
        })
    }
}
