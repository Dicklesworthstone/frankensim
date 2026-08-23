//! fs-conform — the restriction-map plugin conformance SDK (plan addendum,
//! Proposal 7). Layer: L2.
//!
//! The restriction-map layer is where the hard engineering hides: the sheaf
//! organizes bookkeeping GIVEN the trace/conversion operators, and it will
//! faithfully propagate GARBAGE WITH CERTIFICATES ATTACHED if those operators
//! are bad (risk R6). This crate turns that weakest point into an ecosystem
//! play: third parties ship [`Converter`]s (chart-to-chart operators — the Rep
//! Router edges), and a CONFORMANCE SUITE auto-generated from the sheaf axioms
//! certifies each converter into a [`Tier`]. A converter reaches a tier ONLY by
//! passing all three axioms:
//!
//! 1. **Functoriality** — composition agrees (`f∘g == direct`) and identities
//!    act as identities.
//! 2. **Adjoint consistency** — the declared transpose really is the adjoint
//!    the ledger uses: `⟨A x, y⟩ == ⟨x, Aᵀ y⟩`.
//! 3. **Tolerance honesty** — against MANUFACTURED solutions with known
//!    interface traces, the exact error must not exceed the converter's
//!    DECLARED error plus the suite's explicit numerical tolerance. A converter
//!    that understates its error beyond that admitted tolerance FAILS.
//!
//! R6 mitigation: certification is applied to FIRST-PARTY converters with the
//! same severity as third-party ones. The certified tier is meant to be stamped
//! on every ledger entry the converter touches. The SDK control flow is
//! deterministic for a fixed callback transcript; robust evidence arithmetic
//! uses the shared `fs-math` double-double rung.
//!
//! TWO execution lanes exist, with different trust boundaries:
//!
//! * [`certify`] — the LEGACY trusted in-process path. It contains no faults
//!   and proves nothing about later behavior of the same object; its tiers
//!   carry the uncontained no-claim boundary in `CONTRACT.md` and are for
//!   first-party use only.
//! * [`certify_contained`] — the production lane (bead
//!   frankensim-contain-fs-conform-callbacks-6bc6g). Every callback runs
//!   under an admitted work envelope with panic containment and typed failure
//!   classes; transcripts bind inputs, outputs, seeds, budgets, and an
//!   immutable implementation identity; a permuted replay must reproduce the
//!   transcript bitwise before any tier above `Rejected` is minted; and all
//!   axiom arithmetic then runs against the frozen transcript table rather
//!   than live user code. Containment limits remain honest: unwinding panics
//!   are contained, but abort, OOM, native UB, and nontermination are NOT,
//!   so a bounded-callback tier never claims hard isolation.

use fs_math::{dd::Dd, eft::two_sum};

/// A chart-to-chart converter (a Rep Router edge / restriction map). Kept
/// object-safe so heterogeneous third-party converters share one SDK surface.
pub trait Converter {
    /// A stable id (stamped alongside the tier on ledger entries).
    fn id(&self) -> &str;
    /// The source chart dimension.
    fn source_dim(&self) -> usize;
    /// The target chart dimension.
    fn target_dim(&self) -> usize;
    /// Apply the conversion (source → target).
    fn apply(&self, x: &[f64]) -> Vec<f64>;
    /// The DECLARED adjoint/transpose (target → source).
    fn adjoint(&self, y: &[f64]) -> Vec<f64>;
    /// The DECLARED error bound of the converter's error model.
    fn declared_error(&self) -> f64;
}

/// A manufactured solution: an input with its KNOWN exact converted output.
#[derive(Debug, Clone, PartialEq)]
pub struct ManufacturedCase {
    /// The source-chart input.
    pub input: Vec<f64>,
    /// The known-exact target-chart output.
    pub exact_output: Vec<f64>,
}

/// A functoriality witness: `after ∘ self` must equal `direct` on `probes`.
pub struct Composition<'a> {
    /// The converter applied AFTER `self` (target → C).
    pub after: &'a dyn Converter,
    /// The claimed direct converter (source → C).
    pub direct: &'a dyn Converter,
    /// Source-chart probe vectors.
    pub probes: Vec<Vec<f64>>,
}

/// The conformance suite for one converter.
pub struct ConformanceSuite<'a> {
    /// `(x, y)` pairs (source, target) for the adjoint identity.
    pub adjoint_pairs: Vec<(Vec<f64>, Vec<f64>)>,
    /// Manufactured tolerance-honesty cases.
    pub manufactured: Vec<ManufacturedCase>,
    /// An optional functoriality witness (composition).
    pub composition: Option<Composition<'a>>,
    /// An optional identity witness: probes on which a converter CLAIMED to be
    /// the identity map must act as one (`source_dim == target_dim`,
    /// `apply(x) == x`). `None` for converters that are not identities.
    pub identity: Option<Vec<Vec<f64>>>,
    /// Numerical tolerance for the axiom checks.
    pub tolerance: f64,
}

impl ConformanceSuite<'_> {
    /// An incomplete empty suite with the given numerical tolerance. Populate
    /// adjoint and manufactured evidence before calling [`certify`].
    #[must_use]
    pub fn new(tolerance: f64) -> ConformanceSuite<'static> {
        ConformanceSuite {
            adjoint_pairs: Vec::new(),
            manufactured: Vec::new(),
            composition: None,
            identity: None,
            tolerance,
        }
    }
}

/// The certified conformance tier.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Tier {
    /// Failed a hard axiom — NOT certified (do not trust its certificates).
    Rejected,
    /// Certified, coarse admitted error (`declared + suite tolerance`).
    Bronze,
    /// Certified, tight admitted error (`declared + suite tolerance`).
    Silver,
    /// Certified, very tight admitted error (`declared + suite tolerance`).
    Gold,
}

/// The conformance report for a converter.
#[derive(Debug, Clone, PartialEq)]
pub struct ConformanceReport {
    /// The converter id.
    pub converter: String,
    /// Did composition/identity hold? (`true` if no witness supplied.)
    pub functoriality: bool,
    /// Did the adjoint identity hold?
    pub adjoint_consistent: bool,
    /// Did the declared error plus suite tolerance contain the exact
    /// manufactured error?
    pub tolerance_honest: bool,
    /// An outward-rounded upper bound on the worst measured
    /// manufactured-solution error.
    pub measured_error: f64,
    /// The awarded tier.
    pub tier: Tier,
    /// Schema version + rung of the exact arithmetic that decided every
    /// axiom verdict in this report (bead frankensim-i8iva).
    pub arithmetic: ComparisonEvidence,
    /// Human-readable findings (reasons for any failure).
    pub findings: Vec<String>,
}

impl ConformanceReport {
    /// Was the converter certified (any tier above `Rejected`)?
    #[must_use]
    pub fn certified(&self) -> bool {
        self.tier != Tier::Rejected
    }
}

/// Stable schema version of the typed arithmetic evidence carried beside
/// every exact certificate verdict (bead frankensim-i8iva).
pub const CONFORM_ARITHMETIC_SCHEMA_VERSION: u32 = 1;

/// Explicit admitted per-dot dimension budget. The superaccumulator itself
/// allocates nothing and spans the full binary64 product exponent range;
/// this budget bounds the synchronous cold-path work envelope per dot.
pub const CONFORM_DOT_DIMENSION_BUDGET: usize = 1 << 20;

/// The arithmetic rung that decided a certificate comparison.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ArithmeticRung {
    /// Fixed-bin signed superaccumulator over the full binary64 product
    /// exponent range (`SUPERACC_BASE_EXPONENT ..= +2047`). Decides in exact
    /// real arithmetic with fixed stack storage and no allocation.
    ExactSuperaccumulator,
}

impl ArithmeticRung {
    /// Stable wire code under [`CONFORM_ARITHMETIC_SCHEMA_VERSION`].
    #[must_use]
    pub const fn code(self) -> u32 {
        match self {
            Self::ExactSuperaccumulator => 1,
        }
    }
}

/// Why an exact-arithmetic certificate refused to decide. A refusal is NOT a
/// measured failure: nothing is claimed about the converter either way.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ArithmeticRefusal {
    /// The two vectors of one dot product have different lengths.
    DimensionMismatch,
    /// A sampled value was non-finite.
    NonFiniteInput,
    /// A vector exceeded [`CONFORM_DOT_DIMENSION_BUDGET`].
    DimensionBudgetExceeded,
    /// A finite input had no exact integer lattice representation.
    LatticeRefusal,
    /// An exact product or carry left the fixed bin range.
    BinRangeExceeded,
}

/// Status of one exact certificate decision over a witness stream.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExactStatus {
    /// Every witness decided exactly within its tolerance.
    Holds,
    /// At least one witness decided exactly outside its tolerance.
    Violated,
    /// Exact arithmetic could not decide; nothing is claimed.
    Refused,
}

impl ExactStatus {
    /// Stable wire code under [`CONFORM_ARITHMETIC_SCHEMA_VERSION`]
    /// (`0` = refused/no-claim, `1` = holds, `2` = violated).
    #[must_use]
    pub const fn code(self) -> u32 {
        match self {
            Self::Refused => 0,
            Self::Holds => 1,
            Self::Violated => 2,
        }
    }
}

/// Typed arithmetic evidence for one exact certificate decision.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ComparisonEvidence {
    /// [`CONFORM_ARITHMETIC_SCHEMA_VERSION`] of this record's semantics.
    pub schema_version: u32,
    /// Rung that decided (or would have decided) every comparison.
    pub rung: ArithmeticRung,
    /// Witness count (adjoint pairs, probes, or manufactured cases).
    pub terms: usize,
    /// Component dimension each witness vector was held to.
    pub dimension: usize,
    /// Deterministic first refusal site — a flat witness/coordinate index
    /// into the stream — with the refusal reason, if any.
    pub first_refusal: Option<(usize, ArithmeticRefusal)>,
}

/// Outcome of one exact certificate check.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ExactOutcome {
    /// Typed verdict; `Refused` claims nothing either way.
    pub status: ExactStatus,
    /// Typed arithmetic evidence beside the verdict.
    pub evidence: ComparisonEvidence,
}

fn valid_tolerance(tol: f64) -> bool {
    tol.is_finite() && tol >= 0.0
}

fn finite_vector(values: &[f64]) -> bool {
    values.iter().all(|value| value.is_finite())
}

fn finite_dd(value: Dd) -> bool {
    value.hi.is_finite() && value.lo.is_finite()
}

fn zero_dd(value: Dd) -> bool {
    value.hi == 0.0 && value.lo == 0.0
}

fn nonnegative_dd(value: Dd) -> bool {
    finite_dd(value) && !value.lt(Dd::ZERO)
}

fn admitted_bound(declared: f64, tolerance: f64) -> Option<Dd> {
    if !valid_tolerance(declared) || !valid_tolerance(tolerance) {
        return None;
    }
    let declared = Dd::from_f64(declared);
    let tolerance = Dd::from_f64(tolerance);
    let bound = declared + tolerance;
    (nonnegative_dd(bound) && exact_dd_add_represented(declared, tolerance, bound)).then_some(bound)
}

/// Smallest finite `f64` that does not understate a non-negative normalized DD
/// value produced by the `fs-math` arithmetic rung.
fn measured_error_upper(value: Dd) -> Option<f64> {
    if !nonnegative_dd(value) {
        return None;
    }
    let rounded = value.to_f64();
    let upper = if value.lo > 0.0 {
        rounded.next_up()
    } else {
        rounded
    };
    upper.is_finite().then_some(upper)
}

/// Largest power of two no greater than a positive finite `f64`.
///
/// Scaling by a power of two preserves both DD components exactly whenever
/// the scaled value remains representable, which lets the caller detect the
/// underflow boundary by reconstruction instead of silently rounding it away.
fn scale_power(value: f64) -> Option<f64> {
    if !value.is_finite() || value <= 0.0 {
        return None;
    }
    let bits = value.to_bits();
    let exponent = bits & 0x7ff0_0000_0000_0000;
    if exponent != 0 {
        return Some(f64::from_bits(exponent));
    }
    let significand = bits & 0x000f_ffff_ffff_ffff;
    let highest_bit = 63_u32.checked_sub(significand.leading_zeros())?;
    Some(f64::from_bits(1_u64 << highest_bit))
}

/// Exact integer significand and binary exponent for one finite nonzero f64:
/// `abs(value) = significand * 2^exponent`.
fn float_lattice(value: f64) -> Option<(bool, u64, i32)> {
    if !value.is_finite() || value == 0.0 {
        return None;
    }
    let bits = value.to_bits();
    let negative = bits >> 63 != 0;
    let raw_exponent = ((bits >> 52) & 0x7ff) as i32;
    let fraction = bits & 0x000f_ffff_ffff_ffff;
    let (significand, exponent) = if raw_exponent == 0 {
        (fraction, -1074)
    } else {
        ((1_u64 << 52) | fraction, raw_exponent - 1023 - 52)
    };
    let trailing_zeros = significand.trailing_zeros();
    Some((
        negative,
        significand >> trailing_zeros,
        exponent + trailing_zeros as i32,
    ))
}

fn lattice_component_at(value: f64, base_exponent: i32) -> Option<i128> {
    if value == 0.0 {
        return Some(0);
    }
    let (negative, significand, exponent) = float_lattice(value)?;
    let shift = u32::try_from(exponent.checked_sub(base_exponent)?).ok()?;
    let magnitude = u128::from(significand).checked_shl(shift)?;
    if (magnitude >> shift) != u128::from(significand) {
        return None;
    }
    if magnitude > i128::MAX as u128 {
        return None;
    }
    let signed = magnitude as i128;
    Some(if negative { -signed } else { signed })
}

/// Does `represented.hi + represented.lo` equal the exact real product of
/// the two f64 operands? This closes the FMA residual-underflow hole at the
/// normal/subnormal boundary, where a zero residual does not prove exactness.
fn exact_product_represented(left: f64, right: f64, represented: Dd) -> bool {
    let Some((left_negative, left_significand, left_exponent)) = float_lattice(left) else {
        return left == 0.0 && zero_dd(represented);
    };
    let Some((right_negative, right_significand, right_exponent)) = float_lattice(right) else {
        return right == 0.0 && zero_dd(represented);
    };
    let Some(base_exponent) = left_exponent.checked_add(right_exponent) else {
        return false;
    };
    let magnitude = u128::from(left_significand) * u128::from(right_significand);
    if magnitude > i128::MAX as u128 {
        return false;
    }
    let exact = if left_negative == right_negative {
        magnitude as i128
    } else {
        -(magnitude as i128)
    };
    let Some(hi) = lattice_component_at(represented.hi, base_exponent) else {
        return false;
    };
    let Some(lo) = lattice_component_at(represented.lo, base_exponent) else {
        return false;
    };
    hi.checked_add(lo) == Some(exact)
}

/// Exact zero test for `left + right - represented` using a fixed-capacity
/// Shewchuk grow-expansion. Six f64 inputs can produce at most six nonzero
/// components, so this cold certification check allocates nothing.
fn exact_dd_add_represented(left: Dd, right: Dd, represented: Dd) -> bool {
    fn grow(expansion: &mut [f64; 6], len: &mut usize, value: f64) -> bool {
        if value == 0.0 {
            return true;
        }
        let previous = *expansion;
        let previous_len = *len;
        *len = 0;
        let mut q = value;
        for component in previous.into_iter().take(previous_len) {
            let (sum, residual) = two_sum(q, component);
            if !sum.is_finite() || !residual.is_finite() {
                return false;
            }
            if residual != 0.0 {
                if *len == expansion.len() {
                    return false;
                }
                expansion[*len] = residual;
                *len += 1;
            }
            q = sum;
        }
        if q != 0.0 {
            if *len == expansion.len() {
                return false;
            }
            expansion[*len] = q;
            *len += 1;
        }
        true
    }

    if !finite_dd(left) || !finite_dd(right) || !finite_dd(represented) {
        return false;
    }
    let mut expansion = [0.0; 6];
    let mut len = 0usize;
    for component in [
        left.lo,
        left.hi,
        right.lo,
        right.hi,
        -represented.lo,
        -represented.hi,
    ] {
        if !grow(&mut expansion, &mut len, component) {
            return false;
        }
    }
    len == 0
}

// Binary64 products span exponents -2148 through 2047. Seventy limbs also
// leave more than 64 carry bits above that range, enough to sum the maximum
// number of coordinates addressable by usize without wrapping the witness.
const SUPERACC_BASE_EXPONENT: i32 = -2148;
const SUPERACC_LIMBS: usize = 70;

#[derive(Clone, Copy)]
struct PositiveSuperacc {
    limbs: [u64; SUPERACC_LIMBS],
}

impl PositiveSuperacc {
    const ZERO: Self = Self {
        limbs: [0; SUPERACC_LIMBS],
    };

    fn add_word(&mut self, mut index: usize, mut word: u64) -> bool {
        while word != 0 {
            let Some(slot) = self.limbs.get_mut(index) else {
                return false;
            };
            let (sum, carry) = slot.overflowing_add(word);
            *slot = sum;
            word = u64::from(carry);
            index += 1;
        }
        true
    }

    fn add_shifted_u128(&mut self, value: u128, exponent: i32) -> bool {
        if value == 0 {
            return true;
        }
        let Some(bit_offset) = exponent.checked_sub(SUPERACC_BASE_EXPONENT) else {
            return false;
        };
        let Ok(bit_offset) = usize::try_from(bit_offset) else {
            return false;
        };
        let word_index = bit_offset / 64;
        let shift = (bit_offset % 64) as u32;
        let low = value as u64;
        let high = (value >> 64) as u64;
        let first = low << shift;
        let second = if shift == 0 {
            high
        } else {
            (high << shift) | (low >> (64 - shift))
        };
        let third = if shift == 0 { 0 } else { high >> (64 - shift) };
        self.add_word(word_index, first)
            && self.add_word(word_index + 1, second)
            && self.add_word(word_index + 2, third)
    }

    fn add_accumulator(&mut self, other: &Self) -> bool {
        for (index, &word) in other.limbs.iter().enumerate() {
            if !self.add_word(index, word) {
                return false;
            }
        }
        true
    }

    fn le(&self, other: &Self) -> bool {
        for (&left, &right) in self.limbs.iter().zip(&other.limbs).rev() {
            if left != right {
                return left < right;
            }
        }
        true
    }
}

fn add_exact_product(
    left: f64,
    right: f64,
    doubled: bool,
    positive: &mut PositiveSuperacc,
    negative: &mut PositiveSuperacc,
) -> bool {
    if left == 0.0 || right == 0.0 {
        return true;
    }
    let Some((left_negative, left_significand, left_exponent)) = float_lattice(left) else {
        return false;
    };
    let Some((right_negative, right_significand, right_exponent)) = float_lattice(right) else {
        return false;
    };
    let Some(exponent) = left_exponent
        .checked_add(right_exponent)
        .and_then(|value| value.checked_add(i32::from(doubled)))
    else {
        return false;
    };
    let magnitude = u128::from(left_significand) * u128::from(right_significand);
    if left_negative == right_negative {
        positive.add_shifted_u128(magnitude, exponent)
    } else {
        negative.add_shifted_u128(magnitude, exponent)
    }
}

/// Add one finite `f64` exactly into the signed bin pair.
fn add_exact_f64(
    value: f64,
    positive: &mut PositiveSuperacc,
    negative: &mut PositiveSuperacc,
) -> bool {
    if value == 0.0 {
        return true;
    }
    let Some((negative_sign, significand, exponent)) = float_lattice(value) else {
        return false;
    };
    let magnitude = u128::from(significand);
    if negative_sign {
        negative.add_shifted_u128(magnitude, exponent)
    } else {
        positive.add_shifted_u128(magnitude, exponent)
    }
}

/// Exact real-arithmetic signed sum with fixed stack storage.
///
/// Products and single values accumulate into separate positive/negative
/// unsigned bins, so no cancellation ever borrows through a bin and the
/// result is independent of accumulation order.
#[derive(Clone, Copy)]
struct ExactSignedSum {
    positive: PositiveSuperacc,
    negative: PositiveSuperacc,
}

impl ExactSignedSum {
    const ZERO: Self = Self {
        positive: PositiveSuperacc::ZERO,
        negative: PositiveSuperacc::ZERO,
    };

    fn add_product(&mut self, left: f64, right: f64) -> bool {
        add_exact_product(left, right, false, &mut self.positive, &mut self.negative)
    }

    fn add_exact_value(&mut self, value: f64) -> bool {
        add_exact_f64(value, &mut self.positive, &mut self.negative)
    }

    /// Exact real comparison `self <= other + slack`. Cross-adds the
    /// opposite-sign bins so no subtraction ever borrows:
    /// `self - other <= slack` iff `self.pos + other.neg + slack.neg`
    /// `<= other.pos + self.neg + slack.pos`.
    fn le_within_slack(&self, other: &Self, slack: &Self) -> bool {
        let mut left = self.positive;
        let mut right = other.positive;
        left.add_accumulator(&other.negative);
        left.add_accumulator(&slack.negative);
        right.add_accumulator(&self.negative);
        right.add_accumulator(&slack.positive);
        left.le(&right)
    }
}

fn add_exact_dd_square(
    value: Dd,
    positive: &mut PositiveSuperacc,
    negative: &mut PositiveSuperacc,
) -> bool {
    finite_dd(value)
        && add_exact_product(value.hi, value.hi, false, positive, negative)
        && add_exact_product(value.hi, value.lo, true, positive, negative)
        && add_exact_product(value.lo, value.lo, false, positive, negative)
}

/// Exact real-arithmetic comparison `sum_i (a_i-b_i)^2 <= bound^2` over the
/// full finite binary64 exponent range. Positive and negative component terms
/// are accumulated separately, so DD tails of either sign never borrow through
/// an unsigned bin or depend on vector order.
fn squared_norm_le_bound(
    a: &[f64],
    b: &[f64],
    bound: Dd,
) -> Result<bool, (usize, ArithmeticRefusal)> {
    let refuse = |site, reason| Err((site, reason));
    if a.len() != b.len() {
        return refuse(0, ArithmeticRefusal::DimensionMismatch);
    }
    if !finite_vector(a) || !finite_vector(b) {
        return refuse(0, ArithmeticRefusal::NonFiniteInput);
    }
    if !nonnegative_dd(bound) {
        return refuse(0, ArithmeticRefusal::LatticeRefusal);
    }
    let mut norm_positive = PositiveSuperacc::ZERO;
    let mut norm_negative = PositiveSuperacc::ZERO;
    for (index, (&left, &right)) in a.iter().zip(b).enumerate() {
        let delta = Dd::from_f64(left) - Dd::from_f64(right);
        if !exact_dd_add_represented(Dd::from_f64(left), -Dd::from_f64(right), delta)
            || !add_exact_dd_square(delta, &mut norm_positive, &mut norm_negative)
        {
            return refuse(index, ArithmeticRefusal::LatticeRefusal);
        }
    }
    let mut bound_positive = PositiveSuperacc::ZERO;
    let mut bound_negative = PositiveSuperacc::ZERO;
    if !add_exact_dd_square(bound, &mut bound_positive, &mut bound_negative) {
        return refuse(a.len(), ArithmeticRefusal::BinRangeExceeded);
    }

    // norm_pos - norm_neg <= bound_pos - bound_neg
    // iff norm_pos + bound_neg <= bound_pos + norm_neg.
    let mut left = norm_positive;
    let mut right = bound_positive;
    if !left.add_accumulator(&bound_negative) || !right.add_accumulator(&norm_negative) {
        return refuse(a.len(), ArithmeticRefusal::BinRangeExceeded);
    }
    Ok(left.le(&right))
}

fn dist_upper(a: &[f64], b: &[f64]) -> Option<f64> {
    if a.len() != b.len() || !finite_vector(a) || !finite_vector(b) {
        return None;
    }
    let mut largest_delta = Dd::ZERO;
    let mut nonzero_deltas = 0usize;
    for (&left, &right) in a.iter().zip(b) {
        let signed_delta = Dd::from_f64(left) - Dd::from_f64(right);
        if !exact_dd_add_represented(Dd::from_f64(left), -Dd::from_f64(right), signed_delta) {
            return None;
        }
        let delta = signed_delta.abs();
        if !nonnegative_dd(delta) {
            return None;
        }
        if !zero_dd(delta) {
            nonzero_deltas = nonzero_deltas.checked_add(1)?;
        }
        if largest_delta.lt(delta) {
            largest_delta = delta;
        }
    }
    if zero_dd(largest_delta) {
        return Some(0.0);
    }
    if nonzero_deltas == 1 {
        // The Euclidean norm of one nonzero coordinate is its exact absolute
        // DD difference; do not manufacture square/sqrt rounding here.
        return measured_error_upper(largest_delta);
    }

    // LAPACK xLASSQ scaling avoids overflow/underflow in the norm, while the
    // power-of-two divisor preserves DD difference residuals exactly. If even
    // the normalized DD representation loses a nonzero term, this gate cannot
    // certify the evidence and fails closed.
    let scale = scale_power(if largest_delta.hi > 0.0 {
        largest_delta.hi
    } else {
        largest_delta.lo
    })?;
    let scale_dd = Dd::from_f64(scale);
    let mut normalized_squared = Dd::ZERO;
    for (&left, &right) in a.iter().zip(b) {
        let delta = Dd::from_f64(left) - Dd::from_f64(right);
        if !exact_dd_add_represented(Dd::from_f64(left), -Dd::from_f64(right), delta) {
            return None;
        }
        if zero_dd(delta) {
            continue;
        }
        let ratio = delta / scale_dd;
        if !finite_dd(ratio) || zero_dd(ratio) || ratio * scale_dd != delta {
            return None;
        }
        // This rung proves a square only when the exact power-of-two-scaled
        // coordinate is one f64. A two-component ratio needs the exact or
        // outward-enclosed successor path; accepting it here would silently
        // omit the lo^2 term in fs-math's approximate DD multiplication.
        if ratio.lo != 0.0 {
            return None;
        }
        let term = ratio * ratio;
        if !finite_dd(term) || zero_dd(term) || !exact_product_represented(ratio.hi, ratio.hi, term)
        {
            return None;
        }
        let next = normalized_squared + term;
        if !finite_dd(next)
            || !exact_dd_add_represented(normalized_squared, term, next)
            || (next == normalized_squared)
            || (!zero_dd(normalized_squared) && next == term)
        {
            return None;
        }
        normalized_squared = next;
    }
    // Convert the exact two-component squared sum to an outward f64 bound,
    // then use correctly-rounded sqrt plus one successor step. The final
    // power-of-two rescale is exact in DD or refuses. This is the public
    // measured-error projection; the superaccumulator comparison above is the
    // independent boolean authority.
    let normalized_squared_upper = measured_error_upper(normalized_squared)?;
    let rounded_root = normalized_squared_upper.sqrt();
    let rounded_root_square = Dd::from_f64(rounded_root) * Dd::from_f64(rounded_root);
    let root_is_exact = exact_product_represented(rounded_root, rounded_root, rounded_root_square)
        && exact_dd_add_represented(rounded_root_square, -normalized_squared, Dd::ZERO);
    let normalized_distance_upper = if root_is_exact {
        rounded_root
    } else {
        rounded_root.next_up()
    };
    if !normalized_distance_upper.is_finite() {
        return None;
    }
    let scaled = Dd::from_f64(normalized_distance_upper) * scale_dd;
    if !finite_dd(scaled) || !exact_product_represented(normalized_distance_upper, scale, scaled) {
        return None;
    }
    measured_error_upper(scaled)
}

/// Check adjoint consistency `⟨A x, y⟩ == ⟨x, Aᵀ y⟩` over the pairs,
/// deciding every comparison in EXACT real arithmetic (bead frankensim-i8iva).
///
/// The boolean verdict is `true` only for [`ExactStatus::Holds`]; an
/// [`ExactStatus::Refused`] outcome claims nothing and also returns `false`,
/// but the typed evidence names the deterministic first refusal site so an
/// arithmetic uncertainty is never reported as an ordinary measured failure.
#[must_use]
pub fn check_adjoint_with_evidence(
    c: &dyn Converter,
    pairs: &[(Vec<f64>, Vec<f64>)],
    tol: f64,
) -> ExactOutcome {
    let evidence_base = |first_refusal: Option<(usize, ArithmeticRefusal)>| ComparisonEvidence {
        schema_version: CONFORM_ARITHMETIC_SCHEMA_VERSION,
        rung: ArithmeticRung::ExactSuperaccumulator,
        terms: pairs.len(),
        dimension: c.source_dim(),
        first_refusal,
    };
    if pairs.is_empty() || !valid_tolerance(tol) {
        return ExactOutcome {
            status: ExactStatus::Refused,
            evidence: evidence_base(None),
        };
    }
    let (source_dim, target_dim) = (c.source_dim(), c.target_dim());
    // The non-negative tolerance is itself an exact f64 lattice value.
    let mut slack = ExactSignedSum::ZERO;
    if !slack.add_exact_value(tol) {
        return ExactOutcome {
            status: ExactStatus::Refused,
            evidence: evidence_base(Some((0, ArithmeticRefusal::LatticeRefusal))),
        };
    }
    let mut all_consistent = true;
    for (pair_index, (x, y)) in pairs.iter().enumerate() {
        let refuse = |reason| ExactOutcome {
            status: ExactStatus::Refused,
            evidence: evidence_base(Some((pair_index, reason))),
        };
        if x.len() != source_dim || y.len() != target_dim {
            return refuse(ArithmeticRefusal::DimensionMismatch);
        }
        if x.len() > CONFORM_DOT_DIMENSION_BUDGET || y.len() > CONFORM_DOT_DIMENSION_BUDGET {
            return refuse(ArithmeticRefusal::DimensionBudgetExceeded);
        }
        if !finite_vector(x) || !finite_vector(y) {
            return refuse(ArithmeticRefusal::NonFiniteInput);
        }
        let applied = c.apply(x);
        let adjoint = c.adjoint(y);
        if applied.len() != target_dim || adjoint.len() != source_dim {
            return refuse(ArithmeticRefusal::DimensionMismatch);
        }
        if !finite_vector(&applied) || !finite_vector(&adjoint) {
            return refuse(ArithmeticRefusal::NonFiniteInput);
        }
        let mut lhs = ExactSignedSum::ZERO;
        for (&a, &b) in applied.iter().zip(y) {
            if !lhs.add_product(a, b) {
                return refuse(ArithmeticRefusal::BinRangeExceeded);
            }
        }
        let mut rhs = ExactSignedSum::ZERO;
        for (&a, &b) in x.iter().zip(&adjoint) {
            if !rhs.add_product(a, b) {
                return refuse(ArithmeticRefusal::BinRangeExceeded);
            }
        }
        // |lhs - rhs| <= tol  iff  lhs <= rhs + tol AND rhs <= lhs + tol.
        let within = lhs.le_within_slack(&rhs, &slack) && rhs.le_within_slack(&lhs, &slack);
        if !within {
            all_consistent = false;
            break;
        }
    }
    ExactOutcome {
        status: if all_consistent {
            ExactStatus::Holds
        } else {
            ExactStatus::Violated
        },
        evidence: evidence_base(None),
    }
}

/// Boolean adjoint-consistency gate with identical admission semantics to the
/// historical DD path: any structural refusal or exact inconsistency fails.
#[must_use]
pub fn check_adjoint(c: &dyn Converter, pairs: &[(Vec<f64>, Vec<f64>)], tol: f64) -> bool {
    check_adjoint_with_evidence(c, pairs, tol).status == ExactStatus::Holds
}

fn check_tolerance_honesty_with_declared(
    c: &dyn Converter,
    cases: &[ManufacturedCase],
    tol: f64,
    declared: f64,
) -> (bool, f64) {
    let (outcome, measured_upper) = check_tolerance_honesty_evidence_inner(c, cases, tol, declared);
    (outcome.status == ExactStatus::Holds, measured_upper)
}
/// Internal honesty evaluation carrying typed arithmetic evidence beside the
/// boolean verdict and the outward measured-error projection.
fn check_tolerance_honesty_evidence_inner(
    c: &dyn Converter,
    cases: &[ManufacturedCase],
    tol: f64,
    declared: f64,
) -> (ExactOutcome, f64) {
    let evidence_base = |first_refusal: Option<(usize, ArithmeticRefusal)>| ComparisonEvidence {
        schema_version: CONFORM_ARITHMETIC_SCHEMA_VERSION,
        rung: ArithmeticRung::ExactSuperaccumulator,
        terms: cases.len(),
        dimension: c.source_dim(),
        first_refusal,
    };
    if cases.is_empty() || !valid_tolerance(tol) || !declared.is_finite() || declared < 0.0 {
        return (
            ExactOutcome {
                status: ExactStatus::Refused,
                evidence: evidence_base(None),
            },
            f64::INFINITY,
        );
    }
    let Some(admitted_bound) = admitted_bound(declared, tol) else {
        return (
            ExactOutcome {
                status: ExactStatus::Refused,
                evidence: evidence_base(Some((0, ArithmeticRefusal::LatticeRefusal))),
            },
            f64::INFINITY,
        );
    };
    let (source_dim, target_dim) = (c.source_dim(), c.target_dim());
    let mut any_violation = false;
    let mut measured_upper = 0.0_f64;
    for (case_index, case) in cases.iter().enumerate() {
        if case.input.len() != source_dim
            || case.exact_output.len() != target_dim
            || !finite_vector(&case.input)
            || !finite_vector(&case.exact_output)
        {
            return (
                ExactOutcome {
                    status: ExactStatus::Refused,
                    evidence: evidence_base(Some((
                        case_index,
                        ArithmeticRefusal::DimensionMismatch,
                    ))),
                },
                f64::INFINITY,
            );
        }
        let applied = c.apply(&case.input);
        if applied.len() != target_dim {
            return (
                ExactOutcome {
                    status: ExactStatus::Refused,
                    evidence: evidence_base(Some((
                        case_index,
                        ArithmeticRefusal::DimensionMismatch,
                    ))),
                },
                f64::INFINITY,
            );
        }
        let within_bound = match squared_norm_le_bound(&applied, &case.exact_output, admitted_bound)
        {
            Ok(within) => within,
            Err((coordinate, reason)) => {
                return (
                    ExactOutcome {
                        status: ExactStatus::Refused,
                        evidence: evidence_base(Some((case_index + coordinate, reason))),
                    },
                    f64::INFINITY,
                );
            }
        };
        let Some(error_upper) = dist_upper(&applied, &case.exact_output) else {
            return (
                ExactOutcome {
                    status: ExactStatus::Refused,
                    evidence: evidence_base(Some((case_index, ArithmeticRefusal::LatticeRefusal))),
                },
                f64::INFINITY,
            );
        };
        if measured_upper < error_upper {
            measured_upper = error_upper;
        }
        if !within_bound {
            // Keep scanning: `measured_upper` must remain the honest
            // worst-case outward projection across ALL cases even when the
            // verdict is already Violated.
            any_violation = true;
        }
    }
    (
        ExactOutcome {
            status: if any_violation {
                ExactStatus::Violated
            } else {
                ExactStatus::Holds
            },
            evidence: evidence_base(None),
        },
        measured_upper,
    )
}

/// Check tolerance honesty; returns `(honest, outward_worst_measured_error)`.
#[must_use]
pub fn check_tolerance_honesty(
    c: &dyn Converter,
    cases: &[ManufacturedCase],
    tol: f64,
) -> (bool, f64) {
    check_tolerance_honesty_with_declared(c, cases, tol, c.declared_error())
}

/// Tolerance honesty with typed arithmetic evidence: `ExactStatus::Refused`
/// means the exact comparison could not decide and claims nothing; the
/// measured error remains an outward projection over every case.
#[must_use]
pub fn check_tolerance_honesty_with_evidence(
    c: &dyn Converter,
    cases: &[ManufacturedCase],
    tol: f64,
) -> (ExactOutcome, f64) {
    check_tolerance_honesty_evidence_inner(c, cases, tol, c.declared_error())
}

/// Functoriality witness (`after ∘ self == direct` on every probe) decided in
/// exact real arithmetic, with typed arithmetic evidence beside the verdict.
pub fn check_functoriality_with_evidence(
    c: &dyn Converter,
    comp: &Composition,
    tol: f64,
) -> ExactOutcome {
    let evidence_base = |first_refusal: Option<(usize, ArithmeticRefusal)>| ComparisonEvidence {
        schema_version: CONFORM_ARITHMETIC_SCHEMA_VERSION,
        rung: ArithmeticRung::ExactSuperaccumulator,
        terms: comp.probes.len(),
        dimension: c.source_dim(),
        first_refusal,
    };
    let refused = |first_refusal| ExactOutcome {
        status: ExactStatus::Refused,
        evidence: evidence_base(first_refusal),
    };
    if comp.probes.is_empty()
        || !valid_tolerance(tol)
        || c.target_dim() != comp.after.source_dim()
        || c.source_dim() != comp.direct.source_dim()
        || comp.after.target_dim() != comp.direct.target_dim()
    {
        return refused(None);
    }
    let (source_dim, middle_dim, target_dim) =
        (c.source_dim(), c.target_dim(), comp.after.target_dim());
    let mut all_hold = true;
    for (probe_index, x) in comp.probes.iter().enumerate() {
        if x.len() != source_dim {
            return refused(Some((probe_index, ArithmeticRefusal::DimensionMismatch)));
        }
        if !finite_vector(x) {
            return refused(Some((probe_index, ArithmeticRefusal::NonFiniteInput)));
        }
        let middle = c.apply(x);
        if middle.len() != middle_dim || !finite_vector(&middle) {
            return refused(Some((probe_index, ArithmeticRefusal::DimensionMismatch)));
        }
        let composed = comp.after.apply(&middle);
        let direct = comp.direct.apply(x);
        if composed.len() != target_dim || direct.len() != target_dim {
            return refused(Some((probe_index, ArithmeticRefusal::DimensionMismatch)));
        }
        match squared_norm_le_bound(&composed, &direct, Dd::from_f64(tol)) {
            Ok(true) => {}
            Ok(false) => all_hold = false,
            Err((coordinate, reason)) => {
                return refused(Some((probe_index * source_dim + coordinate, reason)));
            }
        }
    }
    ExactOutcome {
        status: if all_hold {
            ExactStatus::Holds
        } else {
            ExactStatus::Violated
        },
        evidence: evidence_base(None),
    }
}

/// Boolean functoriality gate over the exact superaccumulator rung.
#[must_use]
pub fn check_functoriality(c: &dyn Converter, comp: &Composition, tol: f64) -> bool {
    check_functoriality_with_evidence(c, comp, tol).status == ExactStatus::Holds
}

/// Check that a converter claiming to be an identity acts as one, with typed
/// arithmetic evidence beside the verdict.
#[must_use]
pub fn check_identity_with_evidence(
    c: &dyn Converter,
    probes: &[Vec<f64>],
    tol: f64,
) -> ExactOutcome {
    let evidence_base = |first_refusal: Option<(usize, ArithmeticRefusal)>| ComparisonEvidence {
        schema_version: CONFORM_ARITHMETIC_SCHEMA_VERSION,
        rung: ArithmeticRung::ExactSuperaccumulator,
        terms: probes.len(),
        dimension: c.source_dim(),
        first_refusal,
    };
    if probes.is_empty() || !valid_tolerance(tol) || c.source_dim() != c.target_dim() {
        return ExactOutcome {
            status: ExactStatus::Refused,
            evidence: evidence_base(None),
        };
    }
    let dim = c.source_dim();
    let mut all_hold = true;
    for (probe_index, x) in probes.iter().enumerate() {
        if x.len() != dim {
            return ExactOutcome {
                status: ExactStatus::Refused,
                evidence: evidence_base(Some((probe_index, ArithmeticRefusal::DimensionMismatch))),
            };
        }
        if !finite_vector(x) {
            return ExactOutcome {
                status: ExactStatus::Refused,
                evidence: evidence_base(Some((probe_index, ArithmeticRefusal::NonFiniteInput))),
            };
        }
        let applied = c.apply(x);
        if applied.len() != dim {
            return ExactOutcome {
                status: ExactStatus::Refused,
                evidence: evidence_base(Some((probe_index, ArithmeticRefusal::DimensionMismatch))),
            };
        }
        match squared_norm_le_bound(&applied, x, Dd::from_f64(tol)) {
            Ok(true) => {}
            Ok(false) => all_hold = false,
            Err((coordinate, reason)) => {
                return ExactOutcome {
                    status: ExactStatus::Refused,
                    evidence: evidence_base(Some((probe_index * dim + coordinate, reason))),
                };
            }
        }
    }
    ExactOutcome {
        status: if all_hold {
            ExactStatus::Holds
        } else {
            ExactStatus::Violated
        },
        evidence: evidence_base(None),
    }
}

/// Boolean identity gate over the exact superaccumulator rung.
#[must_use]
pub fn check_identity(c: &dyn Converter, probes: &[Vec<f64>], tol: f64) -> bool {
    check_identity_with_evidence(c, probes, tol).status == ExactStatus::Holds
}

/// Push the honest adjoint-consistency finding, distinguishing an exact
/// arithmetic refusal (which claims nothing) from a genuine violation.
fn push_adjoint_finding(findings: &mut Vec<String>, outcome: Option<&ExactOutcome>) {
    findings.push(match outcome {
        None => "adjoint consistency: no witness pairs supplied".to_string(),
        Some(outcome) => match (outcome.status, outcome.evidence.first_refusal) {
            (ExactStatus::Refused, Some((site, reason))) => format!(
                "adjoint consistency: exact arithmetic refused at pair {site} \
                 ({reason:?}); nothing is claimed about the converter"
            ),
            (ExactStatus::Refused, None) => {
                "adjoint consistency: exact arithmetic refused before any pair; \
                 nothing is claimed about the converter"
                    .to_string()
            }
            _ => "adjoint consistency: <Ax,y> != <x,Aᵀy> (declared transpose is not \
                  the adjoint, decided in exact arithmetic)"
                .to_string(),
        },
    });
}

fn evaluate_functoriality(c: &dyn Converter, suite: &ConformanceSuite) -> (bool, Vec<String>) {
    let mut findings = Vec::new();
    let composition_ok = match &suite.composition {
        Some(comp) if comp.probes.is_empty() => {
            findings.push("functoriality: supplied composition has no probes".to_string());
            false
        }
        Some(comp) => {
            let ok = check_functoriality(c, comp, suite.tolerance);
            if !ok {
                findings.push(format!(
                    "functoriality: {} ∘ {} != direct",
                    comp.after.id(),
                    c.id()
                ));
            }
            ok
        }
        None => true,
    };
    let identity_ok = match &suite.identity {
        Some(probes) if probes.is_empty() => {
            findings.push("identity: supplied identity witness has no probes".to_string());
            false
        }
        Some(probes) => {
            let ok = check_identity(c, probes, suite.tolerance);
            if !ok {
                findings.push(format!(
                    "identity: {} claims to be an identity but apply(x) != x",
                    c.id()
                ));
            }
            ok
        }
        None => true,
    };
    (composition_ok && identity_ok, findings)
}

/// Certify a converter against its suite through the legacy trusted
/// in-process path. It reaches a tier ABOVE `Rejected` only by passing every
/// supplied axiom; the tier level reflects how tight an (honestly met within
/// the suite tolerance) error model it declares. Adjoint and manufactured
/// evidence must be non-empty; any supplied composition or identity witness
/// must carry at least one probe.
///
/// TRUST BOUNDARY: this path executes arbitrary user code directly and
/// contains nothing — no fault containment, no budget admission, and no
/// replay evidence bind the observed transcript to later invocations. It is
/// reserved for FIRST-PARTY converters under the operator's own trust
/// decision, and any tier it mints carries the uncontained no-claim boundary
/// documented in `CONTRACT.md`. Production third-party certification MUST go
/// through [`certify_contained`], which executes every callback under an
/// admitted envelope, binds a canonical transcript to an immutable
/// implementation identity, and proves deterministic replay before awarding
/// any tier.
#[must_use]
pub fn certify(c: &dyn Converter, suite: &ConformanceSuite) -> ConformanceReport {
    let declared_error = c.declared_error();
    if !valid_tolerance(suite.tolerance) || !declared_error.is_finite() || declared_error < 0.0 {
        let findings = vec![
            "admission: tolerance and declared error must be finite and non-negative".to_string(),
        ];
        return ConformanceReport {
            converter: c.id().to_string(),
            functoriality: false,
            adjoint_consistent: false,
            tolerance_honest: false,
            measured_error: f64::INFINITY,
            tier: Tier::Rejected,
            arithmetic: ComparisonEvidence {
                schema_version: CONFORM_ARITHMETIC_SCHEMA_VERSION,
                rung: ArithmeticRung::ExactSuperaccumulator,
                terms: 0,
                dimension: c.source_dim(),
                first_refusal: None,
            },
            findings,
        };
    }
    assemble_axiom_report(c, suite, declared_error, Vec::new())
}

/// Shared tail of both certification lanes: run the three axiom checks over
/// the supplied converter view, assemble typed arithmetic evidence, findings,
/// and the tier. `pre_findings` carries lane-specific admissions ahead of the
/// axiom findings. No user code runs here for the contained lane — its
/// converter view is a frozen transcript table.
fn assemble_axiom_report(
    c: &dyn Converter,
    suite: &ConformanceSuite,
    declared_error: f64,
    mut findings: Vec<String>,
) -> ConformanceReport {
    let (functoriality, functoriality_findings) = evaluate_functoriality(c, suite);
    findings.extend(functoriality_findings);
    let adjoint_outcome = if suite.adjoint_pairs.is_empty() {
        None
    } else {
        Some(check_adjoint_with_evidence(
            c,
            &suite.adjoint_pairs,
            suite.tolerance,
        ))
    };
    let adjoint_consistent =
        matches!(&adjoint_outcome, Some(outcome) if outcome.status == ExactStatus::Holds);
    if !adjoint_consistent {
        push_adjoint_finding(&mut findings, adjoint_outcome.as_ref());
    }

    let (tolerance_honest, measured_error) = if suite.manufactured.is_empty() {
        (false, f64::INFINITY)
    } else {
        check_tolerance_honesty_with_declared(
            c,
            &suite.manufactured,
            suite.tolerance,
            declared_error,
        )
    };
    if !tolerance_honest {
        findings.push(if suite.manufactured.is_empty() {
            "tolerance honesty: no manufactured cases supplied".to_string()
        } else {
            format!(
                "tolerance honesty: evidence exceeds or cannot be enclosed within declared \
                 {declared_error:.3e} + suite tolerance {:.3e} (outward measured error \
                 {measured_error:.3e})",
                suite.tolerance
            )
        });
    }

    let tier = if functoriality && adjoint_consistent && tolerance_honest {
        tier_for_admitted_error(declared_error, suite.tolerance)
    } else {
        Tier::Rejected
    };

    let arithmetic = adjoint_outcome.map_or(
        ComparisonEvidence {
            schema_version: CONFORM_ARITHMETIC_SCHEMA_VERSION,
            rung: ArithmeticRung::ExactSuperaccumulator,
            terms: 0,
            dimension: c.source_dim(),
            first_refusal: None,
        },
        |outcome| outcome.evidence,
    );
    ConformanceReport {
        converter: c.id().to_string(),
        functoriality,
        adjoint_consistent,
        tolerance_honest,
        measured_error,
        tier,
        arithmetic,
        findings,
    }
}

/// The tier awarded to a converter that passed every axiom, by its exact
/// admitted bound (`declared + suite tolerance`). Charging the tolerance to the
/// tier prevents a loose verification policy from laundering a weak guarantee
/// into Gold.
fn tier_for_admitted_error(declared: f64, tolerance: f64) -> Tier {
    let Some(admitted) = admitted_bound(declared, tolerance) else {
        return Tier::Rejected;
    };
    if !Dd::from_f64(1e-6).lt(admitted) {
        Tier::Gold
    } else if !Dd::from_f64(1e-3).lt(admitted) {
        Tier::Silver
    } else {
        Tier::Bronze
    }
}

/// FNV-1a fold of one u64 word (matches the fs-la golden-hash mechanism).
fn receipt_fold(acc: &mut u64, word: u64) {
    for byte in word.to_le_bytes() {
        *acc ^= u64::from(byte);
        *acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
    }
}

/// Deterministic G5 receipt over a canonical arithmetic transcript
/// (bead frankensim-i8iva).
///
/// Feeds every `ExactStatus` code, evidence field, refusal site, and the RAW
/// superaccumulator limb words for a fixed witness set covering the G0
/// boundary classes (integer cancellation, subnormal underflow products,
/// mixed ~600-exponent scales, f64::MAX extrema, exact equality at tol=0,
/// one-ULP violation). Integer-only folding: the value is identical across
/// debug/release and both reference ISAs. A changed value means certificate
/// bit semantics moved — re-freeze only under docs/GOLDEN_POLICY.md with a
/// plausible root cause.
#[must_use]
pub fn arithmetic_receipt_hash() -> u64 {
    const TOL: f64 = 1e-9;
    type Witness = (usize, usize, Vec<f64>, Vec<f64>, f64);

    struct Scale {
        source_dim: usize,
        target_dim: usize,
    }
    impl Converter for Scale {
        fn id(&self) -> &'static str {
            "receipt-scale"
        }
        fn source_dim(&self) -> usize {
            self.source_dim
        }
        fn target_dim(&self) -> usize {
            self.target_dim
        }
        fn apply(&self, x: &[f64]) -> Vec<f64> {
            x.iter()
                .enumerate()
                .map(|(i, &v)| v * f64::from(i as u32 + 2))
                .collect()
        }
        fn adjoint(&self, y: &[f64]) -> Vec<f64> {
            y.iter()
                .enumerate()
                .map(|(i, &v)| v * f64::from(i as u32 + 2))
                .collect()
        }
        fn declared_error(&self) -> f64 {
            0.0
        }
    }

    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;

    // (source_dim, target_dim, x, y, tol) — distinct decision classes.
    let witnesses: [Witness; 6] = [
        // Exact cancellation to zero on both dot sides.
        (2, 2, vec![9.0, -9.0], vec![5.0, 5.0], TOL),
        // Underflowing products (2^-1078): f64 flushes to zero, bins do not.
        (
            2,
            2,
            vec![2f64.powi(-420), 2f64.powi(-420)],
            vec![2f64.powi(-420), 2f64.powi(-420)],
            TOL,
        ),
        // Mixed scales spanning ~600 binary exponents.
        (
            2,
            2,
            vec![2f64.powi(300), -2f64.powi(-300)],
            vec![2f64.powi(300), 2f64.powi(300)],
            TOL,
        ),
        // Opposite-sign extremum product.
        (1, 1, vec![f64::MAX], vec![-1.0], TOL),
        // One-ULP dishonest transpose: exactly Violated at tol = 0.
        (1, 1, vec![1024.0], vec![1024.0], 0.0),
        // Non-finite witness: typed refusal at site 0.
        (1, 1, vec![1.0], vec![f64::NAN], TOL),
    ];

    for (index, (source_dim, target_dim, x, y, tol)) in witnesses.iter().enumerate() {
        let converter = Scale {
            source_dim: *source_dim,
            target_dim: *target_dim,
        };
        let outcome = check_adjoint_with_evidence(&converter, &[(x.clone(), y.clone())], *tol);
        receipt_fold(&mut acc, u64::from(outcome.status.code()));
        receipt_fold(
            &mut acc,
            u64::from(outcome.evidence.schema_version)
                | (u64::from(outcome.evidence.rung.code()) << 32),
        );
        receipt_fold(&mut acc, outcome.evidence.terms as u64);
        receipt_fold(&mut acc, outcome.evidence.dimension as u64);
        match outcome.evidence.first_refusal {
            Some((site, reason)) => {
                receipt_fold(&mut acc, 1);
                receipt_fold(&mut acc, site as u64);
                receipt_fold(&mut acc, reason_code(reason));
            }
            None => receipt_fold(&mut acc, 0),
        }
        let _ = index;
    }

    // Raw canonical accumulator limbs from a fixed mixed-scale dot:
    // integer words, order-fixed, identical on every ISA.
    let mut lhs = ExactSignedSum::ZERO;
    let mut rhs = ExactSignedSum::ZERO;
    let applied = [2.0, -6.0];
    let adjoint = [1.0, 12.0];
    for (&a, &b) in applied.iter().zip(&cy_fixture()) {
        assert!(lhs.add_product(a, b));
    }
    for (&a, &b) in cx_fixture().iter().zip(&adjoint) {
        assert!(rhs.add_product(a, b));
    }
    for i in 0..SUPERACC_LIMBS {
        receipt_fold(&mut acc, lhs.positive.limbs[i]);
        receipt_fold(&mut acc, lhs.negative.limbs[i]);
        receipt_fold(&mut acc, rhs.positive.limbs[i]);
        receipt_fold(&mut acc, rhs.negative.limbs[i]);
    }
    acc
}

fn cx_fixture() -> [f64; 2] {
    [2f64.powi(300), -2f64.powi(-300)]
}

fn cy_fixture() -> [f64; 2] {
    [2f64.powi(300), 2f64.powi(300)]
}

/// Golden sentinel for [`arithmetic_receipt_hash`] (bead frankensim-i8iva).
///
/// Reproduced identically on arm64-macOS and x86_64-Linux, in both debug and
/// release profiles (four runs, rch remote evidence retained in the bead).
pub const ARITHMETIC_RECEIPT_GOLDEN_HASH: u64 = 0x82a1_6112_5e2e_dee7;

/// Stable integer code for [`ArithmeticRefusal`] (declaration order).
#[must_use]
const fn reason_code(reason: ArithmeticRefusal) -> u64 {
    match reason {
        ArithmeticRefusal::DimensionMismatch => 1,
        ArithmeticRefusal::NonFiniteInput => 2,
        ArithmeticRefusal::DimensionBudgetExceeded => 3,
        ArithmeticRefusal::LatticeRefusal => 4,
        ArithmeticRefusal::BinRangeExceeded => 5,
    }
}

// ---------------------------------------------------------------------------
// Contained execution boundary
// (bead frankensim-contain-fs-conform-callbacks-6bc6g)
//
// The legacy [`Converter`] trait executes arbitrary caller code with no fault
// containment, no resource admission, and no binding between the transcript a
// certification run observed and whatever the same object does later. This
// section adds the certifiable lane: every callback invocation runs under an
// admitted work envelope with panic containment and typed failure classes,
// produces a canonical transcript record, and must replay byte-identically
// under witness permutation before any tier above `Rejected` can be minted.
//
// HONEST CONTAINMENT LIMITS: `catch_unwind` contains Rust panics only. It
// CANNOT contain `abort`, allocation failure (OOM), native undefined
// behavior, or nontermination. A bounded-callback tier therefore claims
// fault-typed, budget-bounded, replay-verified TRANSCRIPT evidence — never
// hard isolation. A real process/guest isolation boundary remains explicitly
// out of scope and feature-gated future work.
// ---------------------------------------------------------------------------

/// Stable schema version of the contained-execution protocol.
pub const CONTAINED_PROTOCOL_SCHEMA_VERSION: u32 = 1;

/// Refusal reasons for sealing an [`ImplementationIdentity`]. Admission is
/// fail-closed: an identity that cannot be validated never reaches execution.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IdentityRefusal {
    /// The converter id string was empty.
    EmptyConverterId,
    /// Declared error was not finite.
    NonFiniteDeclaredError,
    /// Declared error was negative.
    NegativeDeclaredError,
    /// Source dimension was zero.
    ZeroSourceDimension,
    /// Target dimension was zero.
    ZeroTargetDimension,
    /// `max_calls` was zero: at least one call must be admitted.
    ZeroMaxCalls,
    /// The output budget cannot hold one declared output vector.
    InsufficientOutputBudget {
        /// Largest output length one call can produce.
        required: usize,
    },
    /// The work budget cannot fund even the smallest legal call.
    WorkEnvelopeTooSmall {
        /// Minimum work units one minimal call needs.
        required: usize,
    },
}

/// Declared per-pass execution envelope for one converter. Both certification
/// passes (original and permuted replay) receive a fresh budget of this size;
/// a pass that exceeds it fails closed with a typed fault.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct WorkEnvelope {
    /// Maximum number of callback invocations admitted per pass.
    pub max_calls: usize,
    /// Maximum total f64 elements read plus written per pass.
    pub max_work_units: usize,
    /// Maximum length of any single output vector.
    pub max_output_len: usize,
}

/// Deterministic seed policy. A converter that derives randomness MUST derive
/// it from this fixed seed; ambient entropy would be caught by the replay
/// audit as nondeterminism, but declaring the policy up front makes the
/// identity self-describing and folds the seed into every receipt.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SeedPolicy {
    /// One fixed u64 seed for the whole certification run.
    Fixed(u64),
}

/// Immutable, content-bound identity of one converter implementation. The
/// digest folds every field with distinct tags via FNV-1a — the same
/// deterministic-correlation class as the arithmetic receipts (NOT
/// cryptographic; see the crate contract). Any change to id, dimensions,
/// declared error, seed policy, or envelope changes the digest and therefore
/// every transcript receipt derived from it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ImplementationIdentity {
    /// Stable converter id stamped on the tier.
    pub converter_id: String,
    /// Source chart dimension.
    pub source_dim: usize,
    /// Target chart dimension.
    pub target_dim: usize,
    /// Raw bit pattern of the declared error bound.
    pub declared_error_bits: u64,
    /// Declared deterministic seed policy.
    pub seed_policy: SeedPolicy,
    /// Declared per-pass work envelope.
    pub envelope: WorkEnvelope,
    /// FNV-1a digest over all fields above plus the protocol version.
    pub digest: u64,
}

impl ImplementationIdentity {
    /// Validate every field and seal the immutable identity. Fails closed on
    /// any field that could make later admission or execution ill-defined.
    pub fn seal(
        converter_id: &str,
        source_dim: usize,
        target_dim: usize,
        declared_error: f64,
        seed_policy: SeedPolicy,
        envelope: WorkEnvelope,
    ) -> Result<ImplementationIdentity, IdentityRefusal> {
        if converter_id.is_empty() {
            return Err(IdentityRefusal::EmptyConverterId);
        }
        if !declared_error.is_finite() {
            return Err(IdentityRefusal::NonFiniteDeclaredError);
        }
        if declared_error < 0.0 {
            return Err(IdentityRefusal::NegativeDeclaredError);
        }
        if source_dim == 0 {
            return Err(IdentityRefusal::ZeroSourceDimension);
        }
        if target_dim == 0 {
            return Err(IdentityRefusal::ZeroTargetDimension);
        }
        if envelope.max_calls == 0 {
            return Err(IdentityRefusal::ZeroMaxCalls);
        }
        let widest_output = source_dim.max(target_dim);
        if envelope.max_output_len < widest_output {
            return Err(IdentityRefusal::InsufficientOutputBudget {
                required: widest_output,
            });
        }
        let min_work = source_dim.saturating_add(target_dim);
        if envelope.max_work_units < min_work {
            return Err(IdentityRefusal::WorkEnvelopeTooSmall { required: min_work });
        }
        let mut acc = u64::from(CONTAINED_PROTOCOL_SCHEMA_VERSION)
            .wrapping_mul(0x0000_0100_0000_01b3)
            .wrapping_add(0x2d5c_1a4b);
        receipt_fold(&mut acc, 0x6964_656e_7469_7479); // tag "identity"
        fold_str(&mut acc, converter_id);
        receipt_fold(&mut acc, source_dim as u64);
        receipt_fold(&mut acc, target_dim as u64);
        receipt_fold(&mut acc, declared_error.to_bits());
        receipt_fold(
            &mut acc,
            match seed_policy {
                SeedPolicy::Fixed(seed) => seed,
            },
        );
        receipt_fold(&mut acc, envelope.max_calls as u64);
        receipt_fold(&mut acc, envelope.max_work_units as u64);
        receipt_fold(&mut acc, envelope.max_output_len as u64);
        Ok(ImplementationIdentity {
            converter_id: converter_id.to_string(),
            source_dim,
            target_dim,
            declared_error_bits: declared_error.to_bits(),
            seed_policy,
            envelope,
            digest: acc,
        })
    }

    /// The declared error bound reconstructed from its sealed bits.
    #[must_use]
    pub fn declared_error(&self) -> f64 {
        f64::from_bits(self.declared_error_bits)
    }

    /// Re-seal from fields and compare digests. Runtime use verifies an
    /// incoming identity this way before accepting any tier, closing the
    /// time-of-check versus time-of-use hole on hand-built identities.
    #[must_use]
    pub fn verify(&self) -> bool {
        Self::seal(
            &self.converter_id,
            self.source_dim,
            self.target_dim,
            self.declared_error(),
            self.seed_policy,
            self.envelope,
        )
        .is_ok_and(|resealed| resealed.digest == self.digest)
    }
}

fn fold_str(acc: &mut u64, value: &str) {
    receipt_fold(acc, value.len() as u64);
    for byte in value.as_bytes() {
        receipt_fold(acc, u64::from(*byte));
    }
}

/// Which role a converter plays in a certification run. The main converter is
/// the candidate being certified; composition witnesses execute two further
/// converters whose transcripts are bound into the same receipt.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum TranscriptRole {
    Main = 0,
    After = 1,
    Direct = 2,
}

/// Which protocol method produced one transcript record.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum TranscriptKind {
    Apply = 0,
    Adjoint = 1,
}

/// One canonical transcript record: role, kind, input bits, output bits.
/// Bit patterns (not float equality) are compared so `-0.0` and `+0.0`
/// inputs/outputs can never silently alias.
#[derive(Debug, Clone, PartialEq, Eq)]
struct TranscriptRecord {
    role: TranscriptRole,
    kind: TranscriptKind,
    input_bits: Vec<u64>,
    output_bits: Vec<u64>,
}

impl TranscriptRecord {
    fn key(&self) -> (u8, u8, &[u64], &[u64]) {
        (
            self.role as u8,
            self.kind as u8,
            &self.input_bits[..],
            &self.output_bits[..],
        )
    }
}

fn sort_records(records: &mut [TranscriptRecord]) {
    records.sort_by(|left, right| left.key().cmp(&right.key()));
}

/// Typed failure of one bounded callback invocation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CallFault {
    /// Input vector had the wrong dimension.
    InputDimension {
        /// Expected dimension.
        expected: usize,
        /// Observed dimension.
        got: usize,
    },
    /// First non-finite input element site.
    NonFiniteInput {
        /// Element index.
        index: usize,
    },
    /// Output vector had the wrong dimension.
    OutputDimension {
        /// Expected dimension.
        expected: usize,
        /// Observed dimension.
        got: usize,
    },
    /// First non-finite output element site.
    NonFiniteOutput {
        /// Element index.
        index: usize,
    },
    /// The call budget of this pass ran out.
    CallBudgetExhausted,
    /// The work budget of this pass ran out.
    WorkBudgetExhausted,
    /// An output vector exceeded the declared output budget.
    OutputBudgetExceeded {
        /// Declared maximum output length.
        max: usize,
    },
    /// The callback panicked. The panic payload is deliberately dropped: its
    /// contents are arbitrary host state, not deterministic diagnostics.
    /// This contains unwinding only — never abort, OOM, UB, or loops.
    Panicked,
}

/// Budget meter charged by bounded calls. Implementors of
/// [`ContainedConverter`] drive it through [`CallBudget::begin_call`] and
/// [`CallBudget::complete_call`]; the certification executor creates one per
/// pass from the sealed envelope.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CallBudget {
    calls_left: usize,
    work_left: usize,
}

impl CallBudget {
    /// A fresh budget for one pass over the declared envelope.
    #[must_use]
    pub fn new(envelope: WorkEnvelope) -> CallBudget {
        CallBudget {
            calls_left: envelope.max_calls,
            work_left: envelope.max_work_units,
        }
    }

    /// Admit one call reading `input_len` input elements. Charges the call
    /// count and the input read against the budget.
    ///
    /// # Errors
    /// [`CallFault::CallBudgetExhausted`] when no calls remain,
    /// [`CallFault::WorkBudgetExhausted`] when the input cannot be funded.
    pub fn begin_call(&mut self, input_len: usize) -> Result<(), CallFault> {
        if self.calls_left == 0 {
            return Err(CallFault::CallBudgetExhausted);
        }
        if self.work_left < input_len {
            return Err(CallFault::WorkBudgetExhausted);
        }
        self.calls_left -= 1;
        self.work_left -= input_len;
        Ok(())
    }

    /// Charge the observed output length after the callback returned.
    pub fn complete_call(
        &mut self,
        envelope: WorkEnvelope,
        output_len: usize,
    ) -> Result<(), CallFault> {
        if output_len > envelope.max_output_len {
            return Err(CallFault::OutputBudgetExceeded {
                max: envelope.max_output_len,
            });
        }
        if self.work_left < output_len {
            return Err(CallFault::WorkBudgetExhausted);
        }
        self.work_left -= output_len;
        Ok(())
    }
}

/// The certifiable converter protocol. Every method is fallible and
/// budget-metered; outputs land in caller-provided storage instead of fresh
/// allocations chosen by untrusted code.
pub trait ContainedConverter {
    /// The immutable implementation identity this object is bound to.
    fn identity(&self) -> &ImplementationIdentity;
    /// Bounded forward application (source → target).
    ///
    /// # Errors
    /// Returns a typed [`CallFault`] for dimension/finiteness violations,
    /// budget exhaustion, or a contained panic. `out` is cleared first and
    /// left empty unless the call completed successfully.
    fn apply_bounded(
        &self,
        x: &[f64],
        out: &mut Vec<f64>,
        budget: &mut CallBudget,
    ) -> Result<(), CallFault>;
    /// Bounded declared adjoint (target → source), same contract as
    /// [`ContainedConverter::apply_bounded`].
    ///
    /// # Errors
    /// Returns a typed [`CallFault`] under the same conditions.
    fn adjoint_bounded(
        &self,
        y: &[f64],
        out: &mut Vec<f64>,
        budget: &mut CallBudget,
    ) -> Result<(), CallFault>;
}

fn validate_input(input: &[f64], expected: usize) -> Result<(), CallFault> {
    if input.len() != expected {
        return Err(CallFault::InputDimension {
            expected,
            got: input.len(),
        });
    }
    for (index, value) in input.iter().enumerate() {
        if !value.is_finite() {
            return Err(CallFault::NonFiniteInput { index });
        }
    }
    Ok(())
}

fn validate_output(output: &[f64], expected: usize) -> Result<(), CallFault> {
    if output.len() != expected {
        return Err(CallFault::OutputDimension {
            expected,
            got: output.len(),
        });
    }
    for (index, value) in output.iter().enumerate() {
        if !value.is_finite() {
            return Err(CallFault::NonFiniteOutput { index });
        }
    }
    Ok(())
}

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|value| value.to_bits()).collect()
}

fn bits_eq(left: &[f64], right: &[u64]) -> bool {
    left.len() == right.len()
        && left
            .iter()
            .zip(right.iter())
            .all(|(a, b)| a.to_bits() == *b)
}

/// Adapter binding a legacy [`Converter`] trait object into the contained
/// protocol. Dimension and finiteness validation, budget charging, output
/// storage, and panic containment happen HERE, not inside trusted code.
pub struct BoundedCallback<'a> {
    inner: &'a dyn Converter,
    identity: ImplementationIdentity,
}

impl<'a> BoundedCallback<'a> {
    /// Bind a legacy converter under a sealed identity. Fails closed when the
    /// declaration cannot be validated.
    ///
    /// # Errors
    /// Returns [`IdentityRefusal`] when any identity field is invalid.
    pub fn bind(
        inner: &'a dyn Converter,
        seed_policy: SeedPolicy,
        envelope: WorkEnvelope,
    ) -> Result<BoundedCallback<'a>, IdentityRefusal> {
        let identity = ImplementationIdentity::seal(
            inner.id(),
            inner.source_dim(),
            inner.target_dim(),
            inner.declared_error(),
            seed_policy,
            envelope,
        )?;
        Ok(BoundedCallback { inner, identity })
    }

    fn invoke(
        &self,
        kind: TranscriptKind,
        input: &[f64],
        out: &mut Vec<f64>,
    ) -> Result<(), CallFault> {
        let expected_in = match kind {
            TranscriptKind::Apply => self.identity.source_dim,
            TranscriptKind::Adjoint => self.identity.target_dim,
        };
        let expected_out = match kind {
            TranscriptKind::Apply => self.identity.target_dim,
            TranscriptKind::Adjoint => self.identity.source_dim,
        };
        validate_input(input, expected_in)?;
        out.clear();
        let outcome = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| match kind {
            TranscriptKind::Apply => self.inner.apply(input),
            TranscriptKind::Adjoint => self.inner.adjoint(input),
        }));
        match outcome {
            Ok(produced) => {
                *out = produced;
            }
            Err(payload) => {
                drop(payload);
                out.clear();
                return Err(CallFault::Panicked);
            }
        }
        validate_output(out, expected_out)
    }
}

impl ContainedConverter for BoundedCallback<'_> {
    fn identity(&self) -> &ImplementationIdentity {
        &self.identity
    }

    fn apply_bounded(
        &self,
        x: &[f64],
        out: &mut Vec<f64>,
        _budget: &mut CallBudget,
    ) -> Result<(), CallFault> {
        self.invoke(TranscriptKind::Apply, x, out)
    }

    fn adjoint_bounded(
        &self,
        y: &[f64],
        out: &mut Vec<f64>,
        _budget: &mut CallBudget,
    ) -> Result<(), CallFault> {
        self.invoke(TranscriptKind::Adjoint, y, out)
    }
}

/// Containment class declared by the caller for this certification run.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ContainmentClass {
    /// First-party operator code running through the same bounded path.
    /// Severity is uniform (R6), but the isolation claim stays "bounded
    /// callbacks in-process", never "hard isolation".
    FirstPartyBoundedCallbacks,
    /// The production third-party lane.
    ThirdPartyBoundedCallbacks,
    /// Exploratory adapters. Never executes and never mints a tier.
    ExploratoryUncontained,
}

/// Caller-supplied policy for one contained certification run.
#[derive(Clone, Copy)]
pub struct ContainmentPolicy<'a> {
    /// Declared containment class for this run.
    pub class: ContainmentClass,
    /// Cooperative cancellation probe checked between witness steps. When it
    /// returns true the run fails closed with
    /// [`ExecutionFault::Cancelled`]; cancellation drains cleanly and never
    /// leaves a partial tier behind.
    pub cancelled: Option<&'a dyn Fn() -> bool>,
}

/// Typed failure of a whole contained execution pass.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ExecutionFault {
    /// The cancellation probe fired; the run drained without publishing.
    Cancelled,
    /// The suite supplied a composition witness but no contained binding.
    MissingWitnessBinding,
    /// A bounded call failed with the typed fault at this step index.
    Call {
        /// Index of the failing step within the pass plan.
        step_index: usize,
        /// The typed per-call fault.
        fault: CallFault,
    },
    /// Original and permuted-replay transcripts diverged at this position of
    /// the canonically sorted record sets — typed nondeterminism evidence.
    NondeterministicReplay {
        /// Position of the first differing sorted record.
        first_divergent_record: usize,
    },
}

/// Content-addressed evidence about HOW the transcripts were produced. Every
/// field is deterministic; identical runs across threads, ISAs, and processes
/// produce equal receipts.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExecutionReceipt {
    /// Protocol schema version.
    pub protocol_schema_version: u32,
    /// Declared containment class of the run.
    pub class: ContainmentClass,
    /// Digest of the certified implementation's sealed identity.
    pub identity_digest: u64,
    /// Sorted digests of composition-witness identities (after/direct).
    pub auxiliary_identity_digests: Vec<u64>,
    /// FNV-1a fold over the full sorted transcript set plus identity context.
    pub transcript_receipt: u64,
    /// Did the permuted replay reproduce every record bitwise?
    pub replay_verified: bool,
    /// Did both passes complete with no fault and no cancellation?
    pub drained_cleanly: bool,
    /// First typed failure, when the run did not drain cleanly.
    pub first_fault: Option<ExecutionFault>,
}

/// Result of [`certify_contained`]: the ordinary axiom report plus the
/// execution receipt that says whether the report may be believed.
#[derive(Debug, Clone, PartialEq)]
pub struct ContainedCertification {
    /// The conformance verdict. Its tier is `Rejected` whenever the execution
    /// receipt shows any fault, cancellation, or nondeterminism.
    pub report: ConformanceReport,
    /// The typed execution evidence bound to this run.
    pub execution: ExecutionReceipt,
}

/// Why one planned step did not complete.
enum StepAbort {
    /// Cancellation probe fired before the step ran.
    Cancelled,
    /// The bounded call failed with this typed fault.
    Call(CallFault),
}

#[allow(clippy::needless_pass_by_value)]
fn step_fault(step_index: usize, abort: StepAbort) -> ExecutionFault {
    match abort {
        StepAbort::Cancelled => ExecutionFault::Cancelled,
        StepAbort::Call(fault) => ExecutionFault::Call { step_index, fault },
    }
}

#[allow(clippy::too_many_arguments)]
fn run_one(
    converter: &dyn ContainedConverter,
    role: TranscriptRole,
    kind: TranscriptKind,
    input: &[f64],
    records: &mut Vec<TranscriptRecord>,
    scratch: &mut Vec<f64>,
    budget: &mut CallBudget,
    cancelled: &dyn Fn() -> bool,
) -> Result<Vec<f64>, StepAbort> {
    if cancelled() {
        return Err(StepAbort::Cancelled);
    }
    // Executor-side admission and metering: enforced HERE for every
    // implementor of the trait, independent of adapter cooperation.
    let identity = converter.identity();
    let expected_in = match kind {
        TranscriptKind::Apply => identity.source_dim,
        TranscriptKind::Adjoint => identity.target_dim,
    };
    if let Err(fault) = validate_input(input, expected_in) {
        return Err(StepAbort::Call(fault));
    }
    if let Err(fault) = budget.begin_call(input.len()) {
        return Err(StepAbort::Call(fault));
    }
    let outcome = match kind {
        TranscriptKind::Apply => converter.apply_bounded(input, scratch, budget),
        TranscriptKind::Adjoint => converter.adjoint_bounded(input, scratch, budget),
    };
    outcome.map_err(StepAbort::Call)?;
    let output = std::mem::take(scratch);
    let expected_out = match kind {
        TranscriptKind::Apply => identity.target_dim,
        TranscriptKind::Adjoint => identity.source_dim,
    };
    if let Err(fault) = validate_output(&output, expected_out) {
        *scratch = output;
        return Err(StepAbort::Call(fault));
    }
    if let Err(fault) = budget.complete_call(identity.envelope, output.len()) {
        *scratch = output;
        return Err(StepAbort::Call(fault));
    }
    records.push(TranscriptRecord {
        role,
        kind,
        input_bits: bits(input),
        output_bits: bits(&output),
    });
    Ok(output)
}

fn witness_iter<'a, T>(items: &'a [T], reverse: bool) -> Box<dyn Iterator<Item = &'a T> + 'a> {
    if reverse {
        Box::new(items.iter().rev())
    } else {
        Box::new(items.iter())
    }
}

#[allow(clippy::too_many_lines)]
fn execute_pass(
    main: &dyn ContainedConverter,
    after: Option<&dyn ContainedConverter>,
    direct: Option<&dyn ContainedConverter>,
    suite: &ConformanceSuite,
    policy: ContainmentPolicy<'_>,
    reverse_witnesses: bool,
) -> Result<Vec<TranscriptRecord>, ExecutionFault> {
    let noop = || false;
    let cancelled: &dyn Fn() -> bool = policy.cancelled.unwrap_or(&noop);
    let mut budget = CallBudget::new(main.identity().envelope);
    let mut scratch: Vec<f64> = Vec::new();
    let mut records: Vec<TranscriptRecord> = Vec::new();
    let mut step_index: usize = 0;

    // Adjoint pairs: Apply(x_i), then Adjoint(y_i).
    for pair in witness_iter(&suite.adjoint_pairs, reverse_witnesses) {
        run_one(
            main,
            TranscriptRole::Main,
            TranscriptKind::Apply,
            &pair.0,
            &mut records,
            &mut scratch,
            &mut budget,
            cancelled,
        )
        .map_err(|abort| step_fault(step_index, abort))?;
        step_index += 1;
        run_one(
            main,
            TranscriptRole::Main,
            TranscriptKind::Adjoint,
            &pair.1,
            &mut records,
            &mut scratch,
            &mut budget,
            cancelled,
        )
        .map_err(|abort| step_fault(step_index, abort))?;
        step_index += 1;
    }

    // Manufactured tolerance-honesty cases.
    for case in witness_iter(&suite.manufactured, reverse_witnesses) {
        run_one(
            main,
            TranscriptRole::Main,
            TranscriptKind::Apply,
            &case.input,
            &mut records,
            &mut scratch,
            &mut budget,
            cancelled,
        )
        .map_err(|abort| step_fault(step_index, abort))?;
        step_index += 1;
    }

    // Composition probes: Main(p) -> After(main_out) -> Direct(p). The chain
    // is data-dependent by definition, so it is built during execution; both
    // passes build the same structure because replay compares canonicalized
    // record multisets, not raw orders.
    let probes: &[Vec<f64>] = suite
        .composition
        .as_ref()
        .map_or(&[], |composition| composition.probes.as_slice());
    for probe in witness_iter(probes, reverse_witnesses) {
        let mid = run_one(
            main,
            TranscriptRole::Main,
            TranscriptKind::Apply,
            probe,
            &mut records,
            &mut scratch,
            &mut budget,
            cancelled,
        )
        .map_err(|abort| step_fault(step_index, abort))?;
        step_index += 1;
        let Some(after_converter) = after else {
            return Err(ExecutionFault::MissingWitnessBinding);
        };
        run_one(
            after_converter,
            TranscriptRole::After,
            TranscriptKind::Apply,
            &mid,
            &mut records,
            &mut scratch,
            &mut budget,
            cancelled,
        )
        .map_err(|abort| step_fault(step_index, abort))?;
        step_index += 1;
        let Some(direct_converter) = direct else {
            return Err(ExecutionFault::MissingWitnessBinding);
        };
        run_one(
            direct_converter,
            TranscriptRole::Direct,
            TranscriptKind::Apply,
            probe,
            &mut records,
            &mut scratch,
            &mut budget,
            cancelled,
        )
        .map_err(|abort| step_fault(step_index, abort))?;
        step_index += 1;
    }

    // Identity witness probes.
    if let Some(identity_probes) = &suite.identity {
        for probe in witness_iter(identity_probes, reverse_witnesses) {
            run_one(
                main,
                TranscriptRole::Main,
                TranscriptKind::Apply,
                probe,
                &mut records,
                &mut scratch,
                &mut budget,
                cancelled,
            )
            .map_err(|abort| step_fault(step_index, abort))?;
            step_index += 1;
        }
    }

    if cancelled() {
        return Err(ExecutionFault::Cancelled);
    }
    Ok(records)
}

/// Read-only view over one frozen transcript set. Implements [`Converter`]
/// purely by table lookup so every downstream axiom check runs against the
/// FROZEN transcript instead of live user code — zero callbacks execute
/// during certification arithmetic.
struct ReplayTable {
    identity: ImplementationIdentity,
    role: TranscriptRole,
    records: Vec<TranscriptRecord>,
}

impl ReplayTable {
    fn lookup(&self, kind: TranscriptKind, input: &[f64]) -> Vec<f64> {
        for record in &self.records {
            if record.role == self.role && record.kind == kind && bits_eq(input, &record.input_bits)
            {
                return record
                    .output_bits
                    .iter()
                    .map(|word| f64::from_bits(*word))
                    .collect();
            }
        }
        Vec::new()
    }
}

impl Converter for ReplayTable {
    fn id(&self) -> &str {
        &self.identity.converter_id
    }

    fn source_dim(&self) -> usize {
        self.identity.source_dim
    }

    fn target_dim(&self) -> usize {
        self.identity.target_dim
    }

    fn apply(&self, x: &[f64]) -> Vec<f64> {
        self.lookup(TranscriptKind::Apply, x)
    }

    fn adjoint(&self, y: &[f64]) -> Vec<f64> {
        self.lookup(TranscriptKind::Adjoint, y)
    }

    fn declared_error(&self) -> f64 {
        self.identity.declared_error()
    }
}

/// Contained bindings for both auxiliary converters of a composition
/// witness, sealed under one seed/envelope policy.
pub struct CompositionBindings<'a> {
    /// The converter applied AFTER the candidate (target chart onward).
    pub after: BoundedCallback<'a>,
    /// The claimed direct converter (source chart to the composed target).
    pub direct: BoundedCallback<'a>,
}

/// Bind a composition witness's auxiliary converters for the contained lane.
///
/// # Errors
/// Returns [`IdentityRefusal`] when either auxiliary identity is invalid.
pub fn bind_composition<'a>(
    composition: &'a Composition<'a>,
    seed_policy: SeedPolicy,
    envelope: WorkEnvelope,
) -> Result<CompositionBindings<'a>, IdentityRefusal> {
    Ok(CompositionBindings {
        after: BoundedCallback::bind(composition.after, seed_policy, envelope)?,
        direct: BoundedCallback::bind(composition.direct, seed_policy, envelope)?,
    })
}

fn fold_record(acc: &mut u64, record: &TranscriptRecord) {
    receipt_fold(acc, u64::from(record.role as u8));
    receipt_fold(acc, u64::from(record.kind as u8));
    receipt_fold(acc, record.input_bits.len() as u64);
    for word in &record.input_bits {
        receipt_fold(acc, *word);
    }
    receipt_fold(acc, record.output_bits.len() as u64);
    for word in &record.output_bits {
        receipt_fold(acc, *word);
    }
}

fn transcript_receipt(
    class: ContainmentClass,
    identity_digest: u64,
    auxiliary: &[u64],
    seed_policy: SeedPolicy,
    sorted_records: &[TranscriptRecord],
) -> u64 {
    let mut acc = 0x243f_6a88_85a3_08d3;
    receipt_fold(&mut acc, u64::from(CONTAINED_PROTOCOL_SCHEMA_VERSION));
    receipt_fold(
        &mut acc,
        match class {
            ContainmentClass::FirstPartyBoundedCallbacks => 1,
            ContainmentClass::ThirdPartyBoundedCallbacks => 2,
            ContainmentClass::ExploratoryUncontained => 3,
        },
    );
    receipt_fold(&mut acc, identity_digest);
    receipt_fold(&mut acc, auxiliary.len() as u64);
    for digest in auxiliary {
        receipt_fold(&mut acc, *digest);
    }
    receipt_fold(
        &mut acc,
        match seed_policy {
            SeedPolicy::Fixed(seed) => seed,
        },
    );
    receipt_fold(&mut acc, sorted_records.len() as u64);
    for record in sorted_records {
        fold_record(&mut acc, record);
    }
    acc
}

fn first_divergence(a: &[TranscriptRecord], b: &[TranscriptRecord]) -> Option<usize> {
    if a.len() != b.len() {
        return Some(a.len().min(b.len()));
    }
    for (index, (left, right)) in a.iter().zip(b.iter()).enumerate() {
        if left.key() != right.key() {
            return Some(index);
        }
    }
    None
}

fn refused_report(identity: &ImplementationIdentity, finding: String) -> ConformanceReport {
    ConformanceReport {
        converter: identity.converter_id.clone(),
        functoriality: false,
        adjoint_consistent: false,
        tolerance_honest: false,
        measured_error: f64::INFINITY,
        tier: Tier::Rejected,
        arithmetic: ComparisonEvidence {
            schema_version: CONFORM_ARITHMETIC_SCHEMA_VERSION,
            rung: ArithmeticRung::ExactSuperaccumulator,
            terms: 0,
            dimension: identity.source_dim,
            first_refusal: None,
        },
        findings: vec![finding],
    }
}

fn describe_fault(fault: &ExecutionFault) -> String {
    match fault {
        ExecutionFault::Cancelled => {
            "containment: cancelled before a clean drain; no tier published".to_string()
        }
        ExecutionFault::Call { step_index, fault } => {
            format!("containment: bounded call at step {step_index} failed: {fault:?}")
        }
        ExecutionFault::NondeterministicReplay {
            first_divergent_record,
        } => format!(
            "containment: permuted replay diverged at sorted record \
             {first_divergent_record}; typed nondeterminism evidence"
        ),
        ExecutionFault::MissingWitnessBinding => {
            "containment: composition witness lacks a contained binding".to_string()
        }
    }
}

fn fail_closed(
    identity: &ImplementationIdentity,
    class: ContainmentClass,
    auxiliary_identity_digests: Vec<u64>,
    fault: ExecutionFault,
) -> ContainedCertification {
    let finding = describe_fault(&fault);
    let execution = ExecutionReceipt {
        protocol_schema_version: CONTAINED_PROTOCOL_SCHEMA_VERSION,
        class,
        identity_digest: identity.digest,
        auxiliary_identity_digests,
        transcript_receipt: 0,
        replay_verified: false,
        drained_cleanly: false,
        first_fault: Some(fault),
    };
    ContainedCertification {
        report: refused_report(identity, finding),
        execution,
    }
}

/// Certify a converter through the CONTAINED lane: every callback runs under
/// an admitted envelope with panic containment and typed failure classes,
/// every invocation is transcribed, the permuted replay must reproduce the
/// transcript bitwise, and the axiom checks then run against the frozen
/// transcript table — never against live user code. A tier above `Rejected`
/// is minted ONLY when execution drained cleanly and replay verified; any
/// fault, cancellation, or nondeterminism fails closed with typed evidence
/// and no partial positive tier.
///
/// The declared error bound is taken from the sealed IDENTITY, not from the
/// live object, so a mutable implementation cannot change its declaration
/// between certification and use.
#[allow(clippy::too_many_lines)]
#[must_use]
pub fn certify_contained(
    candidate: &dyn ContainedConverter,
    witnesses: Option<&CompositionBindings<'_>>,
    suite: &ConformanceSuite,
    policy: ContainmentPolicy<'_>,
) -> ContainedCertification {
    let identity = candidate.identity();
    let auxiliary_identity_digests: Vec<u64> = witnesses.map_or_else(Vec::new, |bindings| {
        let mut digests = vec![
            bindings.after.identity.digest,
            bindings.direct.identity.digest,
        ];
        digests.sort_unstable();
        digests
    });

    // Admission: exploratory adapters never execute and never mint tiers.
    if policy.class == ContainmentClass::ExploratoryUncontained {
        return ContainedCertification {
            report: refused_report(
                identity,
                "admission: exploratory uncontained adapters cannot enter the certified lane"
                    .to_string(),
            ),
            execution: refused_receipt(policy.class, identity, auxiliary_identity_digests),
        };
    }
    // Admission: re-seal verification closes TOCTOU on hand-built identities.
    if !identity.verify() {
        return ContainedCertification {
            report: refused_report(
                identity,
                "admission: implementation identity failed re-seal verification".to_string(),
            ),
            execution: refused_receipt(policy.class, identity, auxiliary_identity_digests),
        };
    }
    // Admission: a composition witness in the suite needs contained bindings.
    if suite.composition.is_some() && witnesses.is_none() {
        return ContainedCertification {
            report: refused_report(
                identity,
                "admission: composition witness present without contained bindings".to_string(),
            ),
            execution: refused_receipt(policy.class, identity, auxiliary_identity_digests),
        };
    }

    let (after_conv, direct_conv): (
        Option<&dyn ContainedConverter>,
        Option<&dyn ContainedConverter>,
    ) = match witnesses {
        Some(bindings) => (Some(&bindings.after), Some(&bindings.direct)),
        None => (None, None),
    };

    // Pass A: canonical witness order. Pass B: permuted witnesses, fresh
    // budgets. Either pass aborting fails closed with typed evidence.
    let sorted_a = match run_both_passes(candidate, after_conv, direct_conv, suite, policy) {
        Ok(records) => records,
        Err(fault) => {
            return fail_closed(identity, policy.class, auxiliary_identity_digests, fault);
        }
    };

    let receipt_value = transcript_receipt(
        policy.class,
        identity.digest,
        &auxiliary_identity_digests,
        identity.seed_policy,
        &sorted_a,
    );

    // Frozen tables: certification arithmetic reads transcripts only.
    let main_table = ReplayTable {
        identity: identity.clone(),
        role: TranscriptRole::Main,
        records: sorted_a.clone(),
    };
    let after_table = ReplayTable {
        identity: witnesses.map_or_else(
            || identity.clone(),
            |bindings| bindings.after.identity.clone(),
        ),
        role: TranscriptRole::After,
        records: sorted_a.clone(),
    };
    let direct_table = ReplayTable {
        identity: witnesses.map_or_else(
            || identity.clone(),
            |bindings| bindings.direct.identity.clone(),
        ),
        role: TranscriptRole::Direct,
        records: sorted_a,
    };
    let assembled = ConformanceSuite {
        adjoint_pairs: suite.adjoint_pairs.clone(),
        manufactured: suite.manufactured.clone(),
        composition: suite.composition.as_ref().map(|_| Composition {
            after: &after_table,
            direct: &direct_table,
            probes: suite
                .composition
                .as_ref()
                .map_or(Vec::new(), |composition| composition.probes.clone()),
        }),
        identity: suite.identity.clone(),
        tolerance: suite.tolerance,
    };
    let report = assemble_axiom_report(
        &main_table,
        &assembled,
        identity.declared_error(),
        Vec::new(),
    );
    let execution = ExecutionReceipt {
        protocol_schema_version: CONTAINED_PROTOCOL_SCHEMA_VERSION,
        class: policy.class,
        identity_digest: identity.digest,
        auxiliary_identity_digests,
        transcript_receipt: receipt_value,
        replay_verified: true,
        drained_cleanly: true,
        first_fault: None,
    };
    ContainedCertification { report, execution }
}

fn refused_receipt(
    class: ContainmentClass,
    identity: &ImplementationIdentity,
    auxiliary_identity_digests: Vec<u64>,
) -> ExecutionReceipt {
    ExecutionReceipt {
        protocol_schema_version: CONTAINED_PROTOCOL_SCHEMA_VERSION,
        class,
        identity_digest: identity.digest,
        auxiliary_identity_digests,
        transcript_receipt: 0,
        replay_verified: false,
        drained_cleanly: false,
        first_fault: None,
    }
}

fn run_both_passes(
    candidate: &dyn ContainedConverter,
    after: Option<&dyn ContainedConverter>,
    direct: Option<&dyn ContainedConverter>,
    suite: &ConformanceSuite,
    policy: ContainmentPolicy<'_>,
) -> Result<Vec<TranscriptRecord>, ExecutionFault> {
    let pass_a = execute_pass(candidate, after, direct, suite, policy, false)?;
    let pass_b = execute_pass(candidate, after, direct, suite, policy, true)?;
    let mut sorted_a = pass_a;
    sort_records(&mut sorted_a);
    let mut sorted_b = pass_b;
    sort_records(&mut sorted_b);
    if let Some(first_divergent_record) = first_divergence(&sorted_a, &sorted_b) {
        return Err(ExecutionFault::NondeterministicReplay {
            first_divergent_record,
        });
    }
    Ok(sorted_a)
}

/// Convenience: bind a legacy trait object under a sealed identity and run
/// the contained lane in one call. Suites carrying a composition witness
/// additionally need [`bind_composition`] and the plain [`certify_contained`]
/// entry point.
///
/// # Errors
/// Returns [`IdentityRefusal`] when the candidate's identity is invalid.
pub fn certify_contained_legacy(
    candidate: &dyn Converter,
    suite: &ConformanceSuite,
    seed_policy: SeedPolicy,
    envelope: WorkEnvelope,
    policy: ContainmentPolicy<'_>,
) -> Result<ContainedCertification, IdentityRefusal> {
    let bound = BoundedCallback::bind(candidate, seed_policy, envelope)?;
    Ok(certify_contained(&bound, None, suite, policy))
}
#[cfg(test)]
mod arithmetic_tests {
    use super::*;

    /// G5 sentinel over exact-arithmetic certificate bits (bead
    /// frankensim-i8iva). Frozen after four-way reproduction: arm64-macOS and
    /// x86_64-Linux (rch worker vmi1227854), each in debug AND release,
    /// all agreeing on 0x82a1_6112_5e2e_dee7 at committed tree 8e1134e7+
    /// (receipt fn added in this commit; integer-only folding, no float
    /// control flow). Re-pin only per docs/GOLDEN_POLICY.md.
    #[test]
    fn arithmetic_receipt_is_bit_stable_across_isas_and_profiles() {
        let receipt = arithmetic_receipt_hash();
        assert_eq!(
            receipt, ARITHMETIC_RECEIPT_GOLDEN_HASH,
            "certificate bit semantics moved: {receipt:#018x} vs              {ARITHMETIC_RECEIPT_GOLDEN_HASH:#018x} — re-freeze only under \
             docs/GOLDEN_POLICY.md with a plausible root cause"
        );
    }

    #[test]
    fn superacc_shift_and_long_carry_boundaries_are_exact() {
        let mut cross_word = PositiveSuperacc::ZERO;
        assert!(cross_word.add_shifted_u128((1_u128 << 64) | 1, SUPERACC_BASE_EXPONENT + 63));
        assert_eq!(cross_word.limbs[0], 1_u64 << 63);
        assert_eq!(cross_word.limbs[1], 1_u64 << 63);
        assert!(cross_word.limbs[2..].iter().all(|&limb| limb == 0));

        let mut three_words = PositiveSuperacc::ZERO;
        assert!(three_words.add_shifted_u128((1_u128 << 127) | 1, SUPERACC_BASE_EXPONENT + 1));
        assert_eq!(three_words.limbs[0], 2);
        assert_eq!(three_words.limbs[1], 0);
        assert_eq!(three_words.limbs[2], 1);
        assert!(three_words.limbs[3..].iter().all(|&limb| limb == 0));

        let mut carry_chain = PositiveSuperacc::ZERO;
        carry_chain.limbs[..SUPERACC_LIMBS - 1].fill(u64::MAX);
        assert!(carry_chain.add_word(0, 1));
        assert!(
            carry_chain.limbs[..SUPERACC_LIMBS - 1]
                .iter()
                .all(|&limb| limb == 0)
        );
        assert_eq!(carry_chain.limbs[SUPERACC_LIMBS - 1], 1);

        let mut full = PositiveSuperacc {
            limbs: [u64::MAX; SUPERACC_LIMBS],
        };
        assert!(
            !full.add_word(0, 1),
            "capacity overflow must refuse instead of wrapping into a certificate"
        );
    }

    #[test]
    fn exact_squared_comparison_handles_signed_tails_and_full_range() {
        let min_subnormal = f64::from_bits(1);
        let below_one = Dd {
            hi: 1.0,
            lo: -min_subnormal,
        };
        let above_one = Dd {
            hi: 1.0,
            lo: min_subnormal,
        };
        assert_eq!(squared_norm_le_bound(&[1.0], &[0.0], below_one), Ok(false));
        assert_eq!(squared_norm_le_bound(&[1.0], &[0.0], above_one), Ok(true));

        assert_eq!(
            squared_norm_le_bound(&[f64::MAX], &[0.0], Dd::from_f64(f64::MAX)),
            Ok(true)
        );
        assert_eq!(
            squared_norm_le_bound(&[f64::MAX], &[0.0], Dd::from_f64(f64::MAX.next_down())),
            Ok(false)
        );

        let mut maximum_square = PositiveSuperacc::ZERO;
        let mut negative = PositiveSuperacc::ZERO;
        assert!(add_exact_product(
            f64::MAX,
            f64::MAX,
            false,
            &mut maximum_square,
            &mut negative
        ));
        assert!(negative.limbs.iter().all(|&limb| limb == 0));
        let top_index = usize::try_from(2047 - SUPERACC_BASE_EXPONENT).unwrap() / 64;
        let top_bit = u32::try_from((2047 - SUPERACC_BASE_EXPONENT) % 64).unwrap();
        assert_ne!(maximum_square.limbs[top_index] & (1_u64 << top_bit), 0);
        assert!(
            maximum_square.limbs[top_index + 1..]
                .iter()
                .all(|&limb| limb == 0)
        );
    }

    /// Independent i128 oracle: for integer-valued f64 inputs the exact real
    /// dot product equals this integer, so `<= slack` decisions must agree
    /// bit-for-bit with the superaccumulator across cancellation, mixed
    /// scales, and every permutation.
    #[test]
    fn exact_signed_sum_matches_integer_oracle() {
        // Every entry must be an integer with |v| <= 2^31 so `as i64` and
        // the i128 product are EXACT; the i128 sum is the independent
        // oracle. Subnormal and extreme-scale cases are covered by the
        // integration battery against their hand-computed lattice values.
        let cases: [(&[f64], &[f64]); 7] = [
            (&[3.0, -5.0, 7.0], &[2.0, 4.0, -6.0]),
            // Exact internal cancellation to a nonzero remainder.
            (&[1.0, 1.0], &[1024.0, -1024.0]),
            // Cancellation to exactly zero.
            (&[9.0, -9.0], &[5.0, 5.0]),
            // Scale mix within the exact-integer domain.
            (
                &[2f64.powi(30), -(2f64.powi(20))],
                &[2f64.powi(20), 2f64.powi(10)],
            ),
            (&[-3.0, 0.0, 11.0], &[0.0, 7.0, -2.0]),
            (&[-1.0, -1.0, -1.0], &[1.0, 2.0, 3.0]),
            (&[12345.0, -12345.0, 1.0], &[6789.0, 6789.0, 0.0]),
        ];
        for (a, b) in cases {
            let oracle: i128 = a
                .iter()
                .zip(b)
                .map(|(&x, &y)| i128::from(x as i64) * i128::from(y as i64))
                .sum();
            let mut exact = ExactSignedSum::ZERO;
            for (&x, &y) in a.iter().zip(b) {
                assert!(exact.add_product(x, y));
            }
            // slack = k flips the verdict exactly at k: verdict == (oracle <= k).
            for slack_units in 0..=3_i128 {
                let mut slack = ExactSignedSum::ZERO;
                assert!(slack.add_exact_value(slack_units as f64));
                let mut zero = ExactSignedSum::ZERO;
                assert!(zero.add_exact_value(0.0));
                let expected = oracle <= slack_units;
                assert_eq!(
                    exact.le_within_slack(&zero, &slack),
                    expected,
                    "oracle {oracle} vs slack {slack_units}"
                );
            }
        }
    }
    #[test]
    fn exact_sum_verdict_is_permutation_invariant() {
        let a = [2f64.powi(300), -2f64.powi(80), 5.0, -2f64.powi(-200)];
        let b = [2f64.powi(40), 3.0, -2f64.powi(150), 7.0];
        let mut base = ExactSignedSum::ZERO;
        for (&x, &y) in a.iter().zip(&b) {
            assert!(base.add_product(x, y));
        }
        let mut slack = ExactSignedSum::ZERO;
        assert!(slack.add_exact_value(4.0));
        let permutations: [[usize; 4]; 4] =
            [[0, 1, 2, 3], [3, 2, 1, 0], [1, 3, 0, 2], [2, 0, 3, 1]];
        for order in permutations {
            let mut shuffled = ExactSignedSum::ZERO;
            for &index in &order {
                assert!(shuffled.add_product(a[index], b[index]));
            }
            assert_eq!(
                shuffled.le_within_slack(&ExactSignedSum::ZERO, &slack),
                base.le_within_slack(&ExactSignedSum::ZERO, &slack),
            );
        }
    }
}
