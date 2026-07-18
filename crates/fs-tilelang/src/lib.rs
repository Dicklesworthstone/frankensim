//! fs-tilelang — the safe tile-kernel DSL runtime (plan patch Rev C).
//! Layer: L0 SUBSTRATE.
//!
//! The `kernel!` macro (re-exported from fs-tilelang-macros) lowers
//! ONE restricted kernel body into: a scalar reference variant, a
//! lane-shaped variant (chunked loops the autovectorizer maps onto the
//! resolved SIMD tier — per-element arithmetic is IDENTICAL, so the
//! variants are bitwise-equal by construction), kernel METADATA
//! (arithmetic intensity for the roofline harness and autotuner — P6:
//! every kernel ships its intensity analysis), and generated
//! G0 tier-equivalence + G5 determinism twin tests.
//!
//! This crate is the runtime the generated code targets: metadata
//! types, lane-width resolution (once, via fs-substrate dispatch —
//! never in hot loops), and deterministic/fast reduction combiners.

use core::fmt::{self, Write as _};
pub use fs_tilelang_macros::kernel;

/// Maximum admitted UTF-8 bytes in a public kernel name.
pub const MAX_KERNEL_NAME_BYTES: usize = 256;
/// Maximum admitted UTF-8 bytes in a case or verdict log label.
pub const MAX_LOG_LABEL_BYTES: usize = 128;
/// Maximum bytes emitted by one admitted metadata JSON object.
pub const MAX_METADATA_JSON_BYTES: usize = 2048;
/// Maximum bytes emitted by one admitted structured log record.
pub const MAX_LOG_RECORD_BYTES: usize = 4096;

const METADATA_JSON_FIXED_CAPACITY: usize = 256;
const LOG_RECORD_FIXED_CAPACITY: usize = 128;

fn write_json_string(out: &mut String, value: &str) {
    for ch in value.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\u{0008}' => out.push_str("\\b"),
            '\u{000c}' => out.push_str("\\f"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            '\u{0000}'..='\u{001f}' | '\u{007f}'..='\u{009f}' | '\u{2028}' | '\u{2029}' => {
                write!(out, "\\u{:04x}", u32::from(ch)).expect("writing to a String cannot fail");
            }
            _ => out.push(ch),
        }
    }
}

fn escaped_json_len(value: &str) -> Result<usize, MetadataRenderError> {
    let mut bytes = 0_usize;
    for ch in value.chars() {
        let encoded = match ch {
            '"' | '\\' | '\u{0008}' | '\u{000c}' | '\n' | '\r' | '\t' => 2,
            '\u{0000}'..='\u{001f}' | '\u{007f}'..='\u{009f}' | '\u{2028}' | '\u{2029}' => 6,
            _ => ch.len_utf8(),
        };
        bytes = bytes
            .checked_add(encoded)
            .ok_or(MetadataRenderError::ProjectedSizeOverflow)?;
    }
    Ok(bytes)
}

fn validate_text(
    field: &'static str,
    value: &str,
    limit: usize,
) -> Result<usize, MetadataRenderError> {
    if value.is_empty() {
        return Err(MetadataRenderError::EmptyText { field });
    }
    if value.len() > limit {
        return Err(MetadataRenderError::TextTooLong {
            field,
            actual: value.len(),
            limit,
        });
    }
    escaped_json_len(value)
}

/// Typed refusal from bounded kernel-metadata or structured-log rendering.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum MetadataRenderError {
    /// A semantic text field was empty.
    EmptyText {
        /// Stable field identity.
        field: &'static str,
    },
    /// A semantic text field exceeded its byte budget.
    TextTooLong {
        /// Stable field identity.
        field: &'static str,
        /// Observed UTF-8 byte count.
        actual: usize,
        /// Maximum admitted UTF-8 byte count.
        limit: usize,
    },
    /// Checked projected-size arithmetic overflowed.
    ProjectedSizeOverflow,
}

impl MetadataRenderError {
    fn refusal_json(self) -> String {
        match self {
            Self::EmptyText { field } => format!(
                "{{\"kernel_metadata\":\"refused\",\"rule\":\"empty-text\",\"field\":\"{field}\"}}"
            ),
            Self::TextTooLong {
                field,
                actual,
                limit,
            } => format!(
                "{{\"kernel_metadata\":\"refused\",\"rule\":\"text-too-long\",\"field\":\"{field}\",\"actual_bytes\":{actual},\"limit_bytes\":{limit}}}"
            ),
            Self::ProjectedSizeOverflow => {
                "{\"kernel_metadata\":\"refused\",\"rule\":\"projected-size-overflow\"}".to_owned()
            }
        }
    }
}

impl fmt::Display for MetadataRenderError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyText { field } => write!(f, "{field} must not be empty"),
            Self::TextTooLong {
                field,
                actual,
                limit,
            } => write!(
                f,
                "{field} is {actual} UTF-8 bytes, exceeding the {limit}-byte limit"
            ),
            Self::ProjectedSizeOverflow => f.write_str("projected JSON size overflowed"),
        }
    }
}

impl std::error::Error for MetadataRenderError {}

/// Determinism class a kernel declares (part of its metadata).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DeterminismClass {
    /// Bitwise identical across tiers, batch shapes, and runs.
    BitwiseAllTiers,
    /// Deterministic per tier; reductions envelope-bounded across
    /// tiers (the fs-simd reduction class).
    PerTier,
}

/// Reduction flavor: the DETERMINISTIC variant combines fixed-width
/// chunk partials in index order (a fixed-shape tree keyed by logical
/// position — never by worker); the FAST variant may reassociate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReductionKind {
    /// No reduction output.
    None,
    /// Fixed-shape deterministic sum.
    DeterministicSum,
    /// Reassociation-permitted sum (bit-pattern NOT part of any
    /// contract; must agree with the deterministic variant within an
    /// envelope).
    FastSum,
}

/// Static per-kernel metadata emitted by the macro.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct KernelMeta {
    /// Kernel name (the `kernel!` declaration name). Public literals are
    /// admitted by serializers as nonempty and at most
    /// [`MAX_KERNEL_NAME_BYTES`] UTF-8 bytes.
    pub name: &'static str,
    /// Floating-point operations per processed element (macro-time
    /// count of arithmetic operators and mul_add calls in the body;
    /// mul_add counts as 2).
    pub flops_per_elem: u32,
    /// Bytes moved per processed element (8 per f64 read/write
    /// buffer, 4 per u32 index buffer).
    pub bytes_per_elem: u32,
    /// Declared halo (elements skipped at each end for stencils).
    pub halo: u32,
    /// Reduction flavor.
    pub reduction: ReductionKind,
    /// Determinism class.
    pub determinism: DeterminismClass,
}

impl KernelMeta {
    /// Arithmetic intensity in FLOP/byte (the roofline x-axis).
    #[must_use]
    pub fn intensity(&self) -> f64 {
        f64::from(self.flops_per_elem) / f64::from(self.bytes_per_elem.max(1))
    }

    fn validated_name_escaped_len(&self) -> Result<usize, MetadataRenderError> {
        validate_text("kernel", self.name, MAX_KERNEL_NAME_BYTES)
    }

    fn write_descr(&self, out: &mut String) {
        out.push_str("{\"kernel\":\"");
        write_json_string(out, self.name);
        write!(
            out,
            "\",\"flops_per_elem\":{},\"bytes_per_elem\":{},\"intensity\":{:.4},\"halo\":{},\"reduction\":\"{:?}\",\"determinism\":\"{:?}\"}}",
            self.flops_per_elem,
            self.bytes_per_elem,
            self.intensity(),
            self.halo,
            self.reduction,
            self.determinism,
        )
        .expect("writing to a String cannot fail");
    }

    /// Fallibly render one bounded JSON metadata object
    /// (roofline/autotuner food, ledger-ready).
    #[must_use]
    pub fn try_descr(&self) -> Result<String, MetadataRenderError> {
        let escaped_name_len = self.validated_name_escaped_len()?;
        let capacity = METADATA_JSON_FIXED_CAPACITY
            .checked_add(escaped_name_len)
            .ok_or(MetadataRenderError::ProjectedSizeOverflow)?;
        if capacity > MAX_METADATA_JSON_BYTES {
            return Err(MetadataRenderError::ProjectedSizeOverflow);
        }
        let mut out = String::with_capacity(capacity);
        self.write_descr(&mut out);
        if out.len() > MAX_METADATA_JSON_BYTES {
            return Err(MetadataRenderError::ProjectedSizeOverflow);
        }
        Ok(out)
    }

    /// One bounded JSON metadata object. Ordinary macro-generated names and
    /// legacy names whose escaping was already canonical retain their
    /// historical exact bytes. Invalid public struct literals produce a
    /// bounded, structured refusal object that never repeats attacker text;
    /// callers needing typed refusal should use [`Self::try_descr`].
    #[must_use]
    pub fn descr(&self) -> String {
        self.try_descr()
            .unwrap_or_else(MetadataRenderError::refusal_json)
    }

    /// Render the authoritative bounded outer log record with metadata as a
    /// nested JSON OBJECT, never as an interpolated quoted JSON string.
    #[must_use]
    pub fn render_log_record(
        &self,
        case: &str,
        verdict: &str,
    ) -> Result<String, MetadataRenderError> {
        let escaped_name_len = self.validated_name_escaped_len()?;
        let escaped_case_len = validate_text("case", case, MAX_LOG_LABEL_BYTES)?;
        let escaped_verdict_len = validate_text("verdict", verdict, MAX_LOG_LABEL_BYTES)?;
        let capacity = LOG_RECORD_FIXED_CAPACITY
            .checked_add(METADATA_JSON_FIXED_CAPACITY)
            .and_then(|bytes| bytes.checked_add(escaped_name_len))
            .and_then(|bytes| bytes.checked_add(escaped_case_len))
            .and_then(|bytes| bytes.checked_add(escaped_verdict_len))
            .ok_or(MetadataRenderError::ProjectedSizeOverflow)?;
        if capacity > MAX_LOG_RECORD_BYTES {
            return Err(MetadataRenderError::ProjectedSizeOverflow);
        }

        let mut out = String::with_capacity(capacity);
        out.push_str("{\"suite\":\"fs-tilelang\",\"case\":\"");
        write_json_string(&mut out, case);
        out.push_str("\",\"verdict\":\"");
        write_json_string(&mut out, verdict);
        out.push_str("\",\"detail\":");
        self.write_descr(&mut out);
        out.push('}');
        if out.len() > MAX_LOG_RECORD_BYTES {
            return Err(MetadataRenderError::ProjectedSizeOverflow);
        }
        Ok(out)
    }
}

/// Lane width for the RESOLVED SIMD tier (elements of f64 per lane
/// group): Scalar = 1, NEON = 2, AVX2 = 4, AVX-512 = 8. Resolved once
/// per call site through fs-substrate's cached dispatch — generated
/// kernels hoist this out of their loops.
#[must_use]
pub fn lane_width() -> usize {
    match fs_substrate::dispatch_tier() {
        fs_substrate::SimdTier::Scalar => 1,
        fs_substrate::SimdTier::Neon => 2,
        fs_substrate::SimdTier::Avx2 => 4,
        fs_substrate::SimdTier::Avx512 => 8,
    }
}

/// Chunk quantum for deterministic reductions: partials are formed
/// over fixed 64-element chunks and combined in index order. The
/// shape is a function of LENGTH ONLY — never of tier or thread — so
/// deterministic-sum results are bitwise identical everywhere.
pub const REDUCTION_CHUNK: usize = 64;

/// Fixed-shape deterministic sum: per-chunk sequential partials
/// combined in chunk-index order.
#[must_use]
pub fn deterministic_sum(values: &[f64]) -> f64 {
    let mut total = 0.0f64;
    for chunk in values.chunks(REDUCTION_CHUNK) {
        let mut partial = 0.0f64;
        for &v in chunk {
            partial += v;
        }
        total += partial;
    }
    total
}

/// Reassociation-permitted sum (currently a straight fold; the
/// contract permits future lane-parallel reassociation, which is why
/// its bit pattern is NOT part of any golden).
#[must_use]
pub fn fast_sum(values: &[f64]) -> f64 {
    values.iter().sum()
}

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
