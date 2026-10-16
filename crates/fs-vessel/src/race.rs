//! Stage 4b: the E-RACED CANDIDATE SCREEN — and, first-class, the
//! vessel's own DECLARED RACING CONVENTION.
//!
//! Candidate lips are scored by the (cheap, deterministic) growth
//! objective at the nominal fluid and then screened through the shared
//! `fs-race` core under a noisy validator proxy, so dominated lips die
//! early with anytime validity and their kills are ledgered.
//!
//! The convention this module owns is the part that used to live only
//! in `tests/battery.rs`, which meant no auditor outside this crate
//! could drive it (bead `frankensim-extreal-program-f85xj.2.31`): the
//! vessel's screening losses are scaled by [`SCREEN_SCALE`] and the
//! declared paired-loss support is DATA-DERIVED — the fixture-wide
//! base spread plus the full jitter width, scaled the same way (see
//! [`declared_span`]). That is materially different from the
//! ornithoid's convention (`fs_ornith::screen`, which normalizes onto
//! a fixed ceiling and declares a constant span), and a cross-consumer
//! audit can only compare the two if BOTH are public. This module is
//! the vessel's half of that comparison.
//!
//! The 200× scale is a measured design decision, not decoration:
//! unscaled growth gaps (~1e-3) starve the PairwiseRace betting
//! e-process (measured: 0 eliminations in 400 rounds).

use crate::stability::{VesselProfile, growth_objective};
pub use fs_race::{LossSpan, RaceError};

/// Orr–Sommerfeld stations used by the screening loss.
pub const SCREEN_STATIONS: usize = 3;
/// Orr–Sommerfeld modes used by the screening loss.
pub const SCREEN_MODES: usize = 3;
/// Nominal pour rate the screen scores candidates at.
pub const SCREEN_RATE: f64 = 1.0;
/// Nominal viscosity the screen scores candidates at.
pub const SCREEN_VISCOSITY: f64 = 1.0;
/// The vessel's declared loss NORMALIZATION multiplier. Losses handed
/// to the race are `SCREEN_SCALE × (base + jitter)`.
pub const SCREEN_SCALE: f64 = 200.0;
/// Total width of the deterministic per-observation validator jitter,
/// in UNSCALED objective units (so a paired difference of jitters
/// spans at most this much).
pub const SCREEN_JITTER_WIDTH: f64 = 1e-4;

/// The e-raced lip screen's record — the vessel's P7 evidence row.
#[derive(Debug, Clone, PartialEq)]
pub struct LipScreenReport {
    /// Winner index into the candidate slice.
    pub winner: usize,
    /// Candidates eliminated early.
    pub eliminated: usize,
    /// Race evaluations actually spent.
    pub evaluations_used: u64,
    /// What a fixed-N tournament would have spent.
    pub fixed_n_equivalent: u64,
    /// Deterministic base screening losses (growth objective; lower is
    /// better), BEFORE the vessel's normalization.
    pub losses: Vec<f64>,
    /// The declared paired-loss support this race actually ran under.
    pub declared_span: f64,
}

/// A screen that cannot issue a statistically valid verdict.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ScreenError {
    /// A race needs at least two candidates.
    TooFewCandidates {
        /// Number of supplied candidates.
        count: usize,
    },
    /// A base loss is non-finite, so no declared support exists for it.
    NonFiniteLoss {
        /// Candidate that carries the invalid base loss.
        candidate: usize,
        /// IEEE-754 bits of the invalid value.
        value_bits: u64,
    },
    /// The data-derived span is not a valid [`LossSpan`] (overflow to
    /// non-finite). Fail closed: no span, no race, no verdict.
    InvalidSpan {
        /// IEEE-754 bits of the rejected span.
        span_bits: u64,
    },
    /// The shared race core refused.
    Race(
        /// The structured race refusal.
        RaceError,
    ),
}

impl core::fmt::Display for ScreenError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            ScreenError::TooFewCandidates { count } => {
                write!(f, "a lip screen needs at least two candidates; got {count}")
            }
            ScreenError::NonFiniteLoss {
                candidate,
                value_bits,
            } => write!(
                f,
                "candidate {candidate} has a non-finite base loss ({}), so the vessel convention declares no support",
                f64::from_bits(*value_bits)
            ),
            ScreenError::InvalidSpan { span_bits } => write!(
                f,
                "the data-derived vessel loss span is invalid ({})",
                f64::from_bits(*span_bits)
            ),
            ScreenError::Race(err) => write!(f, "e-raced lip screen: {err}"),
        }
    }
}

impl std::error::Error for ScreenError {}

impl From<RaceError> for ScreenError {
    fn from(err: RaceError) -> Self {
        ScreenError::Race(err)
    }
}

/// The vessel's deterministic base screening losses: the min-max
/// certified modal growth of each candidate lip at the nominal fluid
/// (lower is better).
#[must_use]
pub fn screening_losses(lips: &[f64]) -> Vec<f64> {
    lips.iter()
        .map(|&lip| {
            growth_objective(
                &VesselProfile::carafe(lip),
                SCREEN_RATE,
                SCREEN_VISCOSITY,
                SCREEN_STATIONS,
                SCREEN_MODES,
            )
        })
        .collect()
}

/// The vessel's DECLARED paired-loss support for a base-loss table:
/// the fixture-wide base spread plus the full jitter width, scaled by
/// [`SCREEN_SCALE`]. Data-derived on purpose — the growth gaps are
/// tiny and fixture-dependent, so a constant span would either be a
/// loose (power-destroying) over-declaration or an unsound
/// under-declaration.
///
/// # Errors
/// [`ScreenError::TooFewCandidates`] below two candidates,
/// [`ScreenError::NonFiniteLoss`] for a non-finite base loss, and
/// [`ScreenError::InvalidSpan`] if the derived span is not finite and
/// strictly positive.
pub fn declared_span(base: &[f64]) -> Result<LossSpan, ScreenError> {
    if base.len() < 2 {
        return Err(ScreenError::TooFewCandidates { count: base.len() });
    }
    for (candidate, &value) in base.iter().enumerate() {
        if !value.is_finite() {
            return Err(ScreenError::NonFiniteLoss {
                candidate,
                value_bits: value.to_bits(),
            });
        }
    }
    let base_span = base.iter().copied().fold(f64::NEG_INFINITY, f64::max)
        - base.iter().copied().fold(f64::INFINITY, f64::min);
    let span = SCREEN_SCALE * (base_span + SCREEN_JITTER_WIDTH);
    LossSpan::new(span).map_err(|_| ScreenError::InvalidSpan {
        span_bits: span.to_bits(),
    })
}

/// The vessel's deterministic per-observation validator jitter: a
/// hashed counter (no RNG state), total width [`SCREEN_JITTER_WIDTH`]
/// in unscaled objective units.
fn jitter(candidate: usize, round: u64, seed: u64) -> f64 {
    let mut h = (candidate as u64) << 32 ^ round ^ seed;
    h ^= h << 13;
    h ^= h >> 7;
    h ^= h << 17;
    #[allow(clippy::cast_precision_loss)]
    {
        ((h >> 11) as f64 / (1u64 << 53) as f64 - 0.5) * SCREEN_JITTER_WIDTH
    }
}

/// Race an EXPLICIT base-loss table under the vessel's convention.
///
/// This is the vessel's public racing wrapper, split from
/// [`screen_lips`] so an auditor can drive the vessel's convention
/// over a loss table it did not generate — e.g. another flagship's
/// table, for a cross-consumer convention comparison.
///
/// # Errors
/// [`ScreenError`] as documented on [`declared_span`], plus
/// [`ScreenError::Race`] if the shared race core refuses (an
/// observation outside the declared support, a non-finite loss, or an
/// unregistered candidate gate).
pub fn race_base_losses(base: &[f64], seed: u64) -> Result<LipScreenReport, ScreenError> {
    let span = declared_span(base)?;
    let kills = fs_exec::KillRegistry::new();
    for candidate in 0..base.len() {
        let _ = kills.register(candidate as u64);
    }
    let mut loss = |i: usize, t: u64| (base[i] + jitter(i, t, seed)) * SCREEN_SCALE;
    let out = fs_race::race_field(
        &mut loss,
        base.len(),
        fs_race::RaceSettings::new(span),
        &kills,
    )?;
    Ok(LipScreenReport {
        winner: out.winner,
        eliminated: out.eliminated.len(),
        evaluations_used: out.evaluations_used,
        fixed_n_equivalent: out.fixed_n_equivalent,
        losses: base.to_vec(),
        declared_span: span.get(),
    })
}

/// Screen candidate lips: score them with the growth objective, then
/// e-race them under the vessel's declared convention.
///
/// # Errors
/// [`ScreenError`] exactly as [`race_base_losses`].
pub fn screen_lips(lips: &[f64], seed: u64) -> Result<LipScreenReport, ScreenError> {
    race_base_losses(&screening_losses(lips), seed)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn declared_span_refuses_degenerate_tables() {
        assert_eq!(
            declared_span(&[1.0]),
            Err(ScreenError::TooFewCandidates { count: 1 })
        );
        assert!(matches!(
            declared_span(&[1.0, f64::NAN]),
            Err(ScreenError::NonFiniteLoss { candidate: 1, .. })
        ));
        assert!(matches!(
            declared_span(&[-f64::MAX, f64::MAX]),
            Err(ScreenError::InvalidSpan { .. })
        ));
    }

    #[test]
    fn declared_span_is_the_scaled_spread_plus_jitter() {
        let span = declared_span(&[0.25, 0.75, 0.5]).expect("finite table");
        assert_eq!(
            span.get().to_bits(),
            (SCREEN_SCALE * (0.5 + SCREEN_JITTER_WIDTH)).to_bits()
        );
    }
}
