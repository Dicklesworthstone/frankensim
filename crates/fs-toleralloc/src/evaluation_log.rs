//! Canonical, content-addressed logging for one correlated-stack evaluation.
//!
//! Version one is deliberately a ledger-ready value, not a ledger write. The
//! identity remains an unratified candidate until its schema is promoted into
//! the workspace identity-authority registry.

use core::mem::size_of;

use fs_blake3::{ContentHash, hash_domain};
use fs_evidence::COLOR_ALGEBRA_VERSION;

use super::{
    AdmittedCorrelationModel, ColorRank, CorrelatedStackError, CorrelatedStackReceipt,
    CorrelatedStackTerm, MAX_CORRELATED_STACK_TERM_NAME_BYTES_V1, MAX_CORRELATED_STACK_TERMS_V1,
    propagate_correlated_stack,
};

/// Version of the canonical evaluation-log byte schema.
pub const CORRELATED_STACK_EVALUATION_LOG_SCHEMA_V1: u32 = 1;

/// Version of the correlated-stack numeric algorithm bound into the log.
///
/// This is distinct from the byte-schema version so a future arithmetic change
/// cannot silently reuse a version-one evaluation identity.
pub const CORRELATED_STACK_EVALUATION_ALGORITHM_V1: u32 = 1;

/// Domain-separated BLAKE3 context for the unratified version-one log identity.
pub const CORRELATED_STACK_EVALUATION_IDENTITY_DOMAIN_V1: &str =
    "org.frankensim.fs-toleralloc.correlated-stack-evaluation-log.v1";

const CORRELATED_STACK_EVALUATION_PREIMAGE_MAGIC_V1: &[u8; 8] = b"FSTLOGV1";
const U32_BYTES: usize = size_of::<u32>();
const U64_BYTES: usize = size_of::<u64>();
const F64_BYTES: usize = size_of::<f64>();
const HASH_BYTES: usize = 32;
const COLOR_TAG_BYTES: usize = 1;
const PUBLISHED_QUANTITY_COUNT: usize = 5;
const TERM_FIXED_BYTES: usize = U64_BYTES + U64_BYTES + F64_BYTES + COLOR_TAG_BYTES + F64_BYTES;

/// Largest canonical preimage admitted by the version-one wrapper.
///
/// The bound covers a maximum-width namespace, a dense `128 x 128` factor,
/// 128 maximum-width term names, and every fixed-width field. It is a byte
/// admission bound, not a persistence or transport guarantee.
pub const MAX_CORRELATED_STACK_EVALUATION_LOG_BYTES_V1: usize =
    CORRELATED_STACK_EVALUATION_PREIMAGE_MAGIC_V1.len()
        + 3 * U32_BYTES
        + U64_BYTES
        + CORRELATED_STACK_EVALUATION_IDENTITY_DOMAIN_V1.len()
        + U64_BYTES
        + super::MAX_CORRELATION_MODEL_NAMESPACE_BYTES_V1
        + U64_BYTES
        + HASH_BYTES
        + U64_BYTES
        + U64_BYTES
        + MAX_CORRELATED_STACK_TERMS_V1 * MAX_CORRELATED_STACK_TERMS_V1 * F64_BYTES
        + F64_BYTES
        + U64_BYTES
        + MAX_CORRELATED_STACK_TERMS_V1
            * (TERM_FIXED_BYTES + MAX_CORRELATED_STACK_TERM_NAME_BYTES_V1)
        + PUBLISHED_QUANTITY_COUNT * F64_BYTES;

/// Nominal identity of one version-one correlated-stack evaluation log.
///
/// This type prevents accidental interchange with an arbitrary content hash.
/// It does not confer authority, authenticity, execution proof, or durability.
#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct CorrelatedStackEvaluationIdV1(ContentHash);

impl CorrelatedStackEvaluationIdV1 {
    /// Exact BLAKE3 digest bytes.
    #[must_use]
    pub fn as_bytes(&self) -> &[u8; 32] {
        self.0.as_bytes()
    }

    /// Lowercase hexadecimal rendering.
    #[must_use]
    pub fn to_hex(&self) -> String {
        self.0.to_hex()
    }
}

impl core::fmt::Debug for CorrelatedStackEvaluationIdV1 {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        formatter
            .debug_tuple("CorrelatedStackEvaluationIdV1")
            .field(&self.0)
            .finish()
    }
}

impl core::fmt::Display for CorrelatedStackEvaluationIdV1 {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        core::fmt::Display::fmt(&self.0, formatter)
    }
}

/// Failure before an atomic evaluation-log value can be published.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CorrelatedStackEvaluationLogErrorV1 {
    /// The underlying correlated-stack evaluator refused the model or terms.
    Stack(CorrelatedStackError),
    /// Checked canonical-preimage size arithmetic overflowed.
    CanonicalSizeOverflow,
    /// The bounded canonical-preimage allocation failed.
    CanonicalAllocation {
        /// Exact byte capacity requested from the allocator.
        required_bytes: usize,
    },
    /// The encoder and its checked size model disagreed.
    CanonicalLengthMismatch {
        /// Size predicted before allocation.
        expected: usize,
        /// Size actually emitted by the encoder.
        actual: usize,
    },
}

impl From<CorrelatedStackError> for CorrelatedStackEvaluationLogErrorV1 {
    fn from(error: CorrelatedStackError) -> Self {
        Self::Stack(error)
    }
}

/// Self-contained result of one successfully logged correlated-stack call.
///
/// The retained receipt contains the exact admitted model, ordered terms, and
/// all five published quantities. The canonical preimage is exposed so a
/// caller can independently verify the candidate identity before deciding how
/// or whether to persist it.
#[derive(Debug, Clone, PartialEq)]
pub struct CorrelatedStackEvaluationLogV1 {
    receipt: CorrelatedStackReceipt,
    canonical_preimage: Vec<u8>,
    identity: CorrelatedStackEvaluationIdV1,
}

impl CorrelatedStackEvaluationLogV1 {
    /// Complete non-forgeable evaluation receipt.
    #[must_use]
    pub const fn receipt(&self) -> &CorrelatedStackReceipt {
        &self.receipt
    }

    /// Exact admitted correlation model used by the evaluation.
    #[must_use]
    pub const fn model(&self) -> &AdmittedCorrelationModel {
        self.receipt.model()
    }

    /// Exact caller-supplied terms in positional factor order.
    #[must_use]
    pub fn terms(&self) -> &[CorrelatedStackTerm] {
        self.receipt.terms()
    }

    /// Versioned canonical identity preimage.
    #[must_use]
    pub fn canonical_preimage(&self) -> &[u8] {
        &self.canonical_preimage
    }

    /// Nominal candidate identity of the complete canonical preimage.
    #[must_use]
    pub const fn identity(&self) -> CorrelatedStackEvaluationIdV1 {
        self.identity
    }

    /// First-order standard deviation under an independence assumption.
    #[must_use]
    pub const fn independent_standard_deviation(&self) -> f64 {
        self.receipt.independent_standard_deviation()
    }

    /// First-order variance under an independence assumption.
    #[must_use]
    pub const fn independent_variance(&self) -> f64 {
        self.receipt.independent_variance()
    }

    /// First-order standard deviation under the admitted factor.
    #[must_use]
    pub const fn correlated_standard_deviation(&self) -> f64 {
        self.receipt.correlated_standard_deviation()
    }

    /// First-order variance under the admitted factor.
    #[must_use]
    pub const fn correlated_variance(&self) -> f64 {
        self.receipt.correlated_variance()
    }

    /// Signed correlated-minus-independent binary64 variance delta.
    #[must_use]
    pub const fn correlation_variance_delta(&self) -> f64 {
        self.receipt.correlation_variance_delta()
    }
}

fn color_rank_tag_v1(rank: ColorRank) -> u8 {
    match rank {
        ColorRank::Estimated => 1,
        ColorRank::Validated => 2,
        ColorRank::Verified => 3,
    }
}

fn checked_add(length: usize, amount: usize) -> Result<usize, CorrelatedStackEvaluationLogErrorV1> {
    length
        .checked_add(amount)
        .ok_or(CorrelatedStackEvaluationLogErrorV1::CanonicalSizeOverflow)
}

fn canonical_preimage_len(
    receipt: &CorrelatedStackReceipt,
) -> Result<usize, CorrelatedStackEvaluationLogErrorV1> {
    let model = receipt.model();
    let mut length = CORRELATED_STACK_EVALUATION_PREIMAGE_MAGIC_V1.len();
    length = checked_add(length, 3 * U32_BYTES)?;
    length = checked_add(length, U64_BYTES)?;
    length = checked_add(length, CORRELATED_STACK_EVALUATION_IDENTITY_DOMAIN_V1.len())?;
    length = checked_add(length, U64_BYTES)?;
    length = checked_add(length, model.namespace().len())?;
    length = checked_add(length, U64_BYTES + HASH_BYTES + U64_BYTES + U64_BYTES)?;
    let factor_bytes = model
        .lower_factor()
        .len()
        .checked_mul(F64_BYTES)
        .ok_or(CorrelatedStackEvaluationLogErrorV1::CanonicalSizeOverflow)?;
    length = checked_add(length, factor_bytes)?;
    length = checked_add(length, F64_BYTES + U64_BYTES)?;
    for term in receipt.terms() {
        length = checked_add(length, TERM_FIXED_BYTES)?;
        length = checked_add(length, term.name.len())?;
    }
    checked_add(length, PUBLISHED_QUANTITY_COUNT * F64_BYTES)
}

fn push_u32(output: &mut Vec<u8>, value: u32) {
    output.extend_from_slice(&value.to_le_bytes());
}

fn push_u64(output: &mut Vec<u8>, value: u64) {
    output.extend_from_slice(&value.to_le_bytes());
}

fn push_usize(
    output: &mut Vec<u8>,
    value: usize,
) -> Result<(), CorrelatedStackEvaluationLogErrorV1> {
    let value = u64::try_from(value)
        .map_err(|_| CorrelatedStackEvaluationLogErrorV1::CanonicalSizeOverflow)?;
    push_u64(output, value);
    Ok(())
}

fn push_f64(output: &mut Vec<u8>, value: f64) {
    push_u64(output, value.to_bits());
}

fn push_len_bytes(
    output: &mut Vec<u8>,
    value: &[u8],
) -> Result<(), CorrelatedStackEvaluationLogErrorV1> {
    push_usize(output, value.len())?;
    output.extend_from_slice(value);
    Ok(())
}

fn encode_canonical_preimage(
    receipt: &CorrelatedStackReceipt,
) -> Result<Vec<u8>, CorrelatedStackEvaluationLogErrorV1> {
    let required_bytes = canonical_preimage_len(receipt)?;
    if required_bytes > MAX_CORRELATED_STACK_EVALUATION_LOG_BYTES_V1 {
        return Err(CorrelatedStackEvaluationLogErrorV1::CanonicalSizeOverflow);
    }
    let mut output = Vec::new();
    output
        .try_reserve_exact(required_bytes)
        .map_err(|_| CorrelatedStackEvaluationLogErrorV1::CanonicalAllocation { required_bytes })?;

    output.extend_from_slice(CORRELATED_STACK_EVALUATION_PREIMAGE_MAGIC_V1);
    push_u32(&mut output, CORRELATED_STACK_EVALUATION_LOG_SCHEMA_V1);
    push_u32(&mut output, CORRELATED_STACK_EVALUATION_ALGORITHM_V1);
    push_u32(&mut output, COLOR_ALGEBRA_VERSION);
    push_len_bytes(
        &mut output,
        CORRELATED_STACK_EVALUATION_IDENTITY_DOMAIN_V1.as_bytes(),
    )?;

    let model = receipt.model();
    push_len_bytes(&mut output, model.namespace().as_bytes())?;
    push_u64(&mut output, model.schema_version().get());
    output.extend_from_slice(&model.semantic_digest());
    push_usize(&mut output, model.dimension())?;
    push_usize(&mut output, model.lower_factor().len())?;
    for &factor in model.lower_factor() {
        push_f64(&mut output, factor);
    }
    push_f64(&mut output, model.max_row_norm_defect());

    push_usize(&mut output, receipt.terms().len())?;
    for (ordinal, term) in receipt.terms().iter().enumerate() {
        push_usize(&mut output, ordinal)?;
        push_len_bytes(&mut output, term.name.as_bytes())?;
        push_f64(&mut output, term.signed_sensitivity);
        output.push(color_rank_tag_v1(term.sensitivity_color));
        push_f64(&mut output, term.standard_deviation);
    }
    push_f64(&mut output, receipt.independent_standard_deviation());
    push_f64(&mut output, receipt.independent_variance());
    push_f64(&mut output, receipt.correlated_standard_deviation());
    push_f64(&mut output, receipt.correlated_variance());
    push_f64(&mut output, receipt.correlation_variance_delta());

    if output.len() != required_bytes {
        return Err(
            CorrelatedStackEvaluationLogErrorV1::CanonicalLengthMismatch {
                expected: required_bytes,
                actual: output.len(),
            },
        );
    }
    Ok(output)
}

/// Evaluate a correlated stack and atomically construct its canonical log.
///
/// The existing raw [`propagate_correlated_stack`] API remains available and
/// is not retroactively logged. This wrapper publishes a value only after the
/// stack evaluation, complete preimage construction, and candidate identity
/// derivation all succeed.
///
/// # Errors
///
/// Returns the underlying stack refusal, checked preimage-size failure,
/// bounded allocation failure, or an internal encoder-size disagreement. No
/// partial log or identity is returned on any error.
pub fn propagate_correlated_stack_logged(
    model: &AdmittedCorrelationModel,
    terms: &[CorrelatedStackTerm],
) -> Result<CorrelatedStackEvaluationLogV1, CorrelatedStackEvaluationLogErrorV1> {
    let receipt = propagate_correlated_stack(model, terms)?;
    let canonical_preimage = encode_canonical_preimage(&receipt)?;
    let identity = CorrelatedStackEvaluationIdV1(hash_domain(
        CORRELATED_STACK_EVALUATION_IDENTITY_DOMAIN_V1,
        &canonical_preimage,
    ));
    Ok(CorrelatedStackEvaluationLogV1 {
        receipt,
        canonical_preimage,
        identity,
    })
}

#[cfg(test)]
mod tests {
    use std::num::NonZeroU64;

    use super::*;

    fn receipt() -> CorrelatedStackReceipt {
        let model = AdmittedCorrelationModel::try_new(
            "gear/log-private-seam",
            NonZeroU64::new(1).expect("one is nonzero"),
            [0x71; 32],
            2,
            vec![1.0, 0.0, 0.8, 0.6],
        )
        .expect("fixture model is admissible");
        propagate_correlated_stack(
            &model,
            &[
                CorrelatedStackTerm {
                    name: "pitch".to_string(),
                    signed_sensitivity: 1.0,
                    sensitivity_color: ColorRank::Verified,
                    standard_deviation: 1.0,
                },
                CorrelatedStackTerm {
                    name: "runout".to_string(),
                    signed_sensitivity: -0.5,
                    sensitivity_color: ColorRank::Validated,
                    standard_deviation: 2.0,
                },
            ],
        )
        .expect("fixture evaluates")
    }

    #[test]
    fn private_derived_and_defect_fields_move_the_canonical_preimage() {
        let baseline = receipt();
        let baseline_bytes = encode_canonical_preimage(&baseline).expect("baseline encodes");

        let mut mutations = Vec::new();
        let mut moved = baseline.clone();
        moved.model.max_row_norm_defect =
            f64::from_bits(moved.model.max_row_norm_defect.to_bits().wrapping_add(1));
        mutations.push(moved);

        let mut moved = baseline.clone();
        moved.independent_standard_deviation = moved.independent_standard_deviation.next_up();
        mutations.push(moved);
        let mut moved = baseline.clone();
        moved.independent_variance = moved.independent_variance.next_up();
        mutations.push(moved);
        let mut moved = baseline.clone();
        moved.correlated_standard_deviation = moved.correlated_standard_deviation.next_up();
        mutations.push(moved);
        let mut moved = baseline.clone();
        moved.correlated_variance = moved.correlated_variance.next_up();
        mutations.push(moved);
        let mut moved = baseline;
        moved.correlation_variance_delta = moved.correlation_variance_delta.next_up();
        mutations.push(moved);

        for mutation in mutations {
            assert_ne!(
                encode_canonical_preimage(&mutation).expect("mutation encodes"),
                baseline_bytes,
            );
        }
    }

    #[test]
    fn color_rank_tags_are_explicit_and_version_bound() {
        assert_eq!(color_rank_tag_v1(ColorRank::Estimated), 1);
        assert_eq!(color_rank_tag_v1(ColorRank::Validated), 2);
        assert_eq!(color_rank_tag_v1(ColorRank::Verified), 3);
        assert_eq!(COLOR_ALGEBRA_VERSION, 2);
    }
}
