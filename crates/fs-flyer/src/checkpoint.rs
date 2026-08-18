//! Structured determinism checkpoints (bead wf-root-guzez.4.10, E3.5).
//! Plan Round-1 (moved early, before physics churn): per-tick PER-SUBSYSTEM
//! digests so a divergence localizes to a subsystem in the tick's causal
//! order, not to a whole-run hash. Two runs that differ answer THREE
//! questions: which tick, which subsystem FIRST in causal order, and what
//! the expected/observed digests were.
//!
//! The channel set is the plan's §5 list (atmosphere, section loads,
//! circulation, propulsion, generalized loads, integrator state). Channels
//! not yet produced by real subsystems (E3.3/E4 land them) are recorded as
//! `absent` — absence is DATA and digests distinctly from every value, so
//! a subsystem silently dropping out of a run is itself a localized
//! divergence, never a quiet match.

use crate::Refusal;
use fs_blake3::hash_domain;

/// Identity domain for per-subsystem tick digests.
pub const SUBSYSTEM_DIGEST_DOMAIN: &str = "org.frankensim.fs-flyer.subsystem-tick.v1";

/// The declared per-tick causal order (plan §5): earlier entries feed
/// later ones. Localization reports the FIRST divergent subsystem in this
/// order — the causal root of the tick's divergence.
pub const SUBSYSTEM_ORDER: [&str; 6] = [
    "atmosphere",
    "section-loads",
    "circulation",
    "propulsion",
    "generalized-loads",
    "integrator-state",
];

/// One tick's structured checkpoint: per-subsystem digests in causal
/// order; `None` = the subsystem was absent this run (data, not a wildcard).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TickCheckpoint {
    /// Tick index.
    pub tick: u32,
    /// Digest (or absence) per SUBSYSTEM_ORDER slot.
    pub digests: [Option<String>; 6],
}

/// A localized divergence between two checkpoint streams.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Divergence {
    /// First divergent tick.
    pub tick: u32,
    /// First divergent subsystem in causal order at that tick.
    pub subsystem: &'static str,
    /// Expected digest ("absent" when the subsystem was absent).
    pub expected: String,
    /// Observed digest ("absent" likewise).
    pub observed: String,
}

/// Index of a subsystem name in the declared order.
///
/// # Errors
/// `subsystem-unknown` naming the registered set.
pub fn subsystem_index(name: &str) -> Result<usize, Refusal> {
    SUBSYSTEM_ORDER
        .iter()
        .position(|s| *s == name)
        .ok_or_else(|| Refusal {
            code: "subsystem-unknown",
            message: format!("subsystem {name:?} is not in the declared causal order"),
            ranked_repairs: vec![format!("registered subsystems: {SUBSYSTEM_ORDER:?}")],
        })
}

/// Digest one subsystem's payload bytes for one tick. The preimage binds
/// the tick, the subsystem name, and the exact bytes — the same bytes on a
/// different channel or tick digest differently.
#[must_use]
pub fn subsystem_digest(tick: u32, subsystem: &str, payload: &[u8]) -> String {
    let mut preimage = Vec::with_capacity(payload.len() + subsystem.len() + 8);
    preimage.extend_from_slice(&tick.to_le_bytes());
    preimage.extend_from_slice(&(subsystem.len() as u32).to_le_bytes());
    preimage.extend_from_slice(subsystem.as_bytes());
    preimage.extend_from_slice(payload);
    hash_domain(SUBSYSTEM_DIGEST_DOMAIN, &preimage).to_hex()
}

/// Builder for one tick's checkpoint: subsystems record their payloads in
/// any order; absent channels stay `None`.
#[derive(Clone, Debug)]
pub struct CheckpointBuilder {
    tick: u32,
    digests: [Option<String>; 6],
}

impl CheckpointBuilder {
    /// Start a checkpoint for `tick`.
    #[must_use]
    pub fn new(tick: u32) -> Self {
        CheckpointBuilder {
            tick,
            digests: Default::default(),
        }
    }

    /// Record a subsystem's payload bytes (f64 payloads enter as exact
    /// little-endian bit patterns upstream).
    ///
    /// # Errors
    /// `subsystem-unknown`; `subsystem-duplicate` (each channel records at
    /// most once per tick — a double write is a wiring bug, not a merge).
    pub fn record(&mut self, subsystem: &str, payload: &[u8]) -> Result<(), Refusal> {
        let idx = subsystem_index(subsystem)?;
        if self.digests[idx].is_some() {
            return Err(Refusal {
                code: "subsystem-duplicate",
                message: format!(
                    "subsystem {subsystem:?} already recorded for tick {}",
                    self.tick
                ),
                ranked_repairs: vec!["one record per subsystem per tick; fix the wiring".into()],
            });
        }
        self.digests[idx] = Some(subsystem_digest(self.tick, subsystem, payload));
        Ok(())
    }

    /// Finish the tick's checkpoint.
    #[must_use]
    pub fn finish(self) -> TickCheckpoint {
        TickCheckpoint {
            tick: self.tick,
            digests: self.digests,
        }
    }
}

/// Compare two checkpoint streams; `Ok(())` iff identical. A mismatch is
/// LOCALIZED: first divergent tick, then the first divergent subsystem in
/// the declared causal order at that tick.
///
/// # Errors
/// `checkpoint-stream-length-mismatch`; `checkpoint-diverged` carrying the
/// [`Divergence`] payload in its message (machine fields also returned via
/// the structured variant [`first_divergence`]).
pub fn verify_streams(
    expected: &[TickCheckpoint],
    observed: &[TickCheckpoint],
) -> Result<(), Refusal> {
    if expected.len() != observed.len() {
        return Err(Refusal {
            code: "checkpoint-stream-length-mismatch",
            message: format!(
                "expected {} ticks, observed {}",
                expected.len(),
                observed.len()
            ),
            ranked_repairs: vec!["compare equal-extent runs (same end_tick_exclusive)".into()],
        });
    }
    match first_divergence(expected, observed) {
        None => Ok(()),
        Some(d) => Err(Refusal {
            code: "checkpoint-diverged",
            message: format!(
                "first divergence at tick {} in subsystem {:?}: expected {}, observed {}",
                d.tick, d.subsystem, d.expected, d.observed
            ),
            ranked_repairs: vec![
                format!(
                    "the causal root is {:?} — inspect that subsystem's inputs at tick {}",
                    d.subsystem, d.tick
                ),
                "downstream subsystem mismatches at the same tick are consequences, not causes"
                    .into(),
            ],
        }),
    }
}

/// The first divergence between two equal-length streams, in (tick, causal
/// order) — `None` iff identical.
#[must_use]
pub fn first_divergence(
    expected: &[TickCheckpoint],
    observed: &[TickCheckpoint],
) -> Option<Divergence> {
    for (e, o) in expected.iter().zip(observed) {
        for (idx, name) in SUBSYSTEM_ORDER.iter().enumerate() {
            if e.digests[idx] != o.digests[idx] {
                let show = |d: &Option<String>| d.clone().unwrap_or_else(|| "absent".to_string());
                return Some(Divergence {
                    tick: e.tick,
                    subsystem: name,
                    expected: show(&e.digests[idx]),
                    observed: show(&o.digests[idx]),
                });
            }
        }
    }
    None
}

/// Pack an f64 slice as exact little-endian bits (the standard payload
/// encoding for physics channels).
#[must_use]
pub fn f64_payload(values: &[f64]) -> Vec<u8> {
    let mut out = Vec::with_capacity(values.len() * 8);
    for v in values {
        out.extend_from_slice(&v.to_bits().to_le_bytes());
    }
    out
}
