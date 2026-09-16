//! Standard-normal samples from a named fs-rand stream.
//!
//! `start_index` counts **draws**, not normals. Each normal consumes exactly
//! two draws via strict [`Stream::next_normal`] (Box–Muller). The ziggurat
//! path is never used. A bare empty `Vec<f64>` is not a refusal.

use super::{Stream, StreamCheckpoint, StreamKey};

/// Declared output budget for one `philox_normals` call (normals, not draws).
pub const PHILOX_NORMALS_MAX_COUNT: usize = 1_048_576;

/// Kernel identity baked into WASM `ok` envelopes.
pub const KERNEL_VERSION: &str = concat!("fs-rand philox_normals ", env!("CARGO_PKG_VERSION"));

/// Whether a [`Refusal`] maps to a registered refusal code or an execution outcome.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TargetKind {
    RefusalCode,
    ExecutionOutcome,
}

/// Typed refusal / budget miss. Never represented as `[]` or `NaN`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Refusal {
    pub code: &'static str,
    pub message: String,
    pub ranked_repairs: Vec<&'static str>,
    pub details: String,
    pub target_kind: TargetKind,
}

/// Opaque admitted request. Forged values are re-checked by [`philox_normals_admitted`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PhiloxNormalsSpec {
    seed: u64,
    stream_kernel: u32,
    tile: u32,
    start_index: u64,
    count: usize,
}

impl PhiloxNormalsSpec {
    #[must_use]
    pub const fn seed(&self) -> u64 {
        self.seed
    }
    #[must_use]
    pub const fn stream_kernel(&self) -> u32 {
        self.stream_kernel
    }
    #[must_use]
    pub const fn tile(&self) -> u32 {
        self.tile
    }
    #[must_use]
    pub const fn start_index(&self) -> u64 {
        self.start_index
    }
    #[must_use]
    pub const fn count(&self) -> usize {
        self.count
    }
}

fn overflow_refusal(start_index: u64, count: usize) -> Refusal {
    let draws = (count as u64).saturating_mul(2);
    Refusal {
        code: "stream-index-overflow",
        message: "This request would exceed the random stream's 64-bit draw counter.".to_string(),
        ranked_repairs: vec![
            "reduce count",
            "lower start_index",
            "start_index=18446744073709551613 count=1 is the last accepted pair",
        ],
        details: format!(
            "{{\"startIndex\":\"{start_index}\",\"draws\":\"{draws}\",\"maxIndex\":\"18446744073709551615\"}}"
        ),
        target_kind: TargetKind::RefusalCode,
    }
}

/// Admit a request. Draws nothing.
pub fn admit_philox_normals(
    seed: u64,
    stream_kernel: u32,
    tile: u32,
    start_index: u64,
    count: usize,
) -> Result<PhiloxNormalsSpec, Refusal> {
    if count == 0 {
        return Err(Refusal {
            code: "invalid-parameter",
            message: "These inputs do not meet the calculation's stated requirements.".to_string(),
            ranked_repairs: vec!["request count >= 1"],
            details: "{\"name\":\"count\",\"value\":0}".to_string(),
            target_kind: TargetKind::RefusalCode,
        });
    }
    if count > PHILOX_NORMALS_MAX_COUNT {
        return Err(Refusal {
            code: "budget-exhausted",
            message: format!(
                "Requested {count} normals exceeds the declared budget of {PHILOX_NORMALS_MAX_COUNT}."
            ),
            ranked_repairs: vec![
                "request fewer normals",
                "issue several calls with advancing start_index",
            ],
            details: format!(
                "{{\"requested\":{count},\"allowed\":{PHILOX_NORMALS_MAX_COUNT},\"unit\":\"normals\"}}"
            ),
            target_kind: TargetKind::ExecutionOutcome,
        });
    }
    let Some(draws) = (count as u64).checked_mul(2) else {
        return Err(overflow_refusal(start_index, count));
    };
    if start_index.checked_add(draws).is_none() {
        return Err(overflow_refusal(start_index, count));
    }
    Ok(PhiloxNormalsSpec {
        seed,
        stream_kernel,
        tile,
        start_index,
        count,
    })
}

/// Generate normals for an admitted spec. Re-runs admission so a forged spec still refuses.
pub fn philox_normals_admitted(spec: &PhiloxNormalsSpec) -> Result<Vec<f64>, Refusal> {
    let spec = admit_philox_normals(
        spec.seed,
        spec.stream_kernel,
        spec.tile,
        spec.start_index,
        spec.count,
    )?;
    let key = StreamKey {
        seed: spec.seed,
        kernel: spec.stream_kernel,
        tile: spec.tile,
    };
    let mut stream = Stream::resume(StreamCheckpoint::current(key, spec.start_index))
        .expect("current() checkpoints use this build's versions");
    let mut out = Vec::with_capacity(spec.count);
    for _ in 0..spec.count {
        let z = stream.next_normal();
        if !z.is_finite() {
            return Err(Refusal {
                code: "invalid-parameter",
                message: "A generated normal was not finite.".to_string(),
                ranked_repairs: vec!["retry with a different start_index"],
                details: format!("{{\"index\":\"{}\"}}", stream.index()),
                target_kind: TargetKind::RefusalCode,
            });
        }
        out.push(z);
    }
    debug_assert_eq!(out.len(), spec.count);
    Ok(out)
}

/// Admit then generate. Never returns an empty buffer as a stand-in for a refusal.
pub fn philox_normals(
    seed: u64,
    stream_kernel: u32,
    tile: u32,
    start_index: u64,
    count: usize,
) -> Result<Vec<f64>, Refusal> {
    let spec = admit_philox_normals(seed, stream_kernel, tile, start_index, count)?;
    philox_normals_admitted(&spec)
}
