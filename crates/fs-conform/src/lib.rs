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
//! R6 mitigation: [`certify`] is applied to FIRST-PARTY converters with the
//! same severity as third-party ones. The certified tier is meant to be stamped
//! on every ledger entry the converter touches. The SDK control flow is
//! deterministic for a fixed callback transcript; robust evidence arithmetic
//! uses the shared `fs-math` double-double rung. The current object-safe trait
//! does not itself contain callback faults or prove that a third-party
//! implementation is pure or deterministic; see the crate contract no-claim
//! boundary.

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

fn nonnegative_dd_le_f64(value: Dd, bound: f64) -> bool {
    nonnegative_dd(value) && valid_tolerance(bound) && !Dd::from_f64(bound).lt(value)
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
        .and_then(|value| value.checked_add(if doubled { 1 } else { 0 }))
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
fn squared_norm_le_bound(a: &[f64], b: &[f64], bound: Dd) -> Option<bool> {
    if a.len() != b.len() || !finite_vector(a) || !finite_vector(b) || !nonnegative_dd(bound) {
        return None;
    }
    let mut norm_positive = PositiveSuperacc::ZERO;
    let mut norm_negative = PositiveSuperacc::ZERO;
    for (&left, &right) in a.iter().zip(b) {
        let delta = Dd::from_f64(left) - Dd::from_f64(right);
        if !exact_dd_add_represented(Dd::from_f64(left), -Dd::from_f64(right), delta)
            || !add_exact_dd_square(delta, &mut norm_positive, &mut norm_negative)
        {
            return None;
        }
    }
    let mut bound_positive = PositiveSuperacc::ZERO;
    let mut bound_negative = PositiveSuperacc::ZERO;
    if !add_exact_dd_square(bound, &mut bound_positive, &mut bound_negative) {
        return None;
    }

    // norm_pos - norm_neg <= bound_pos - bound_neg
    // iff norm_pos + bound_neg <= bound_pos + norm_neg.
    let mut left = norm_positive;
    let mut right = bound_positive;
    if !left.add_accumulator(&bound_negative) || !right.add_accumulator(&norm_negative) {
        return None;
    }
    Some(left.le(&right))
}

fn dot(a: &[f64], b: &[f64]) -> Option<Dd> {
    if a.len() != b.len() || !finite_vector(a) || !finite_vector(b) {
        return None;
    }
    let mut total = Dd::ZERO;
    for (&x, &y) in a.iter().zip(b) {
        let leading_product = x * y;
        if !leading_product.is_finite() || (x != 0.0 && y != 0.0 && leading_product == 0.0) {
            return None;
        }
        let product = Dd::from_f64(x) * Dd::from_f64(y);
        let next = total + product;
        if !finite_dd(product)
            || !exact_product_represented(x, y, product)
            || !finite_dd(next)
            || !exact_dd_add_represented(total, product, next)
            || (!zero_dd(product) && next == total)
            || (!zero_dd(total) && next == product)
        {
            return None;
        }
        total = next;
    }
    Some(total)
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

/// Check adjoint consistency `⟨A x, y⟩ == ⟨x, Aᵀ y⟩` over the pairs.
#[must_use]
pub fn check_adjoint(c: &dyn Converter, pairs: &[(Vec<f64>, Vec<f64>)], tol: f64) -> bool {
    if pairs.is_empty() || !valid_tolerance(tol) {
        return false;
    }
    let (source_dim, target_dim) = (c.source_dim(), c.target_dim());
    pairs.iter().all(|(x, y)| {
        if x.len() != source_dim || y.len() != target_dim || !finite_vector(x) || !finite_vector(y)
        {
            return false;
        }
        let applied = c.apply(x);
        let adjoint = c.adjoint(y);
        if applied.len() != target_dim || adjoint.len() != source_dim {
            return false;
        }
        let (Some(lhs), Some(rhs)) = (dot(&applied, y), dot(x, &adjoint)) else {
            return false;
        };
        let delta = lhs - rhs;
        exact_dd_add_represented(lhs, -rhs, delta) && nonnegative_dd_le_f64(delta.abs(), tol)
    })
}

fn check_tolerance_honesty_with_declared(
    c: &dyn Converter,
    cases: &[ManufacturedCase],
    tol: f64,
    declared: f64,
) -> (bool, f64) {
    if cases.is_empty() || !valid_tolerance(tol) || !declared.is_finite() || declared < 0.0 {
        return (false, f64::INFINITY);
    }
    let Some(admitted_bound) = admitted_bound(declared, tol) else {
        return (false, f64::INFINITY);
    };
    let (source_dim, target_dim) = (c.source_dim(), c.target_dim());
    let mut measured_upper = 0.0_f64;
    let mut all_within_bound = true;
    for case in cases {
        if case.input.len() != source_dim
            || case.exact_output.len() != target_dim
            || !finite_vector(&case.input)
            || !finite_vector(&case.exact_output)
        {
            return (false, f64::INFINITY);
        }
        let applied = c.apply(&case.input);
        if applied.len() != target_dim {
            return (false, f64::INFINITY);
        }
        let Some(within_bound) =
            squared_norm_le_bound(&applied, &case.exact_output, admitted_bound)
        else {
            return (false, f64::INFINITY);
        };
        all_within_bound &= within_bound;
        let Some(error_upper) = dist_upper(&applied, &case.exact_output) else {
            return (false, f64::INFINITY);
        };
        if measured_upper < error_upper {
            measured_upper = error_upper;
        }
    }
    (all_within_bound, measured_upper)
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

/// Check functoriality: `after(self(x)) == direct(x)` on the probes.
#[must_use]
pub fn check_functoriality(c: &dyn Converter, comp: &Composition, tol: f64) -> bool {
    if comp.probes.is_empty()
        || !valid_tolerance(tol)
        || c.target_dim() != comp.after.source_dim()
        || c.source_dim() != comp.direct.source_dim()
        || comp.after.target_dim() != comp.direct.target_dim()
    {
        return false;
    }
    let (source_dim, middle_dim, target_dim) =
        (c.source_dim(), c.target_dim(), comp.after.target_dim());
    comp.probes.iter().all(|x| {
        if x.len() != source_dim || !finite_vector(x) {
            return false;
        }
        let middle = c.apply(x);
        if middle.len() != middle_dim || !finite_vector(&middle) {
            return false;
        }
        let composed = comp.after.apply(&middle);
        let direct = comp.direct.apply(x);
        if composed.len() != target_dim || direct.len() != target_dim {
            return false;
        }
        squared_norm_le_bound(&composed, &direct, Dd::from_f64(tol)) == Some(true)
    })
}

/// Check that a converter claiming to be an identity acts as one.
#[must_use]
pub fn check_identity(c: &dyn Converter, probes: &[Vec<f64>], tol: f64) -> bool {
    if probes.is_empty() || !valid_tolerance(tol) || c.source_dim() != c.target_dim() {
        return false;
    }
    let dim = c.source_dim();
    probes.iter().all(|x| {
        if x.len() != dim || !finite_vector(x) {
            return false;
        }
        let applied = c.apply(x);
        applied.len() == dim && squared_norm_le_bound(&applied, x, Dd::from_f64(tol)) == Some(true)
    })
}

/// Certify a converter against its suite. It reaches a tier ABOVE `Rejected`
/// only by passing every supplied axiom; the tier level reflects how tight an
/// (honestly met within the suite tolerance) error model it declares. Adjoint
/// and manufactured evidence must be non-empty; any supplied composition or
/// identity witness must carry at least one probe.
#[must_use]
pub fn certify(c: &dyn Converter, suite: &ConformanceSuite) -> ConformanceReport {
    let mut findings = Vec::new();
    let declared_error = c.declared_error();
    if !valid_tolerance(suite.tolerance) || !declared_error.is_finite() || declared_error < 0.0 {
        findings.push(
            "admission: tolerance and declared error must be finite and non-negative".to_string(),
        );
        return ConformanceReport {
            converter: c.id().to_string(),
            functoriality: false,
            adjoint_consistent: false,
            tolerance_honest: false,
            measured_error: f64::INFINITY,
            tier: Tier::Rejected,
            findings,
        };
    }

    // Functoriality: composition agrees AND (if the converter claims to be an
    // identity) it acts as the identity.
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
    let functoriality = composition_ok && identity_ok;

    let adjoint_consistent =
        !suite.adjoint_pairs.is_empty() && check_adjoint(c, &suite.adjoint_pairs, suite.tolerance);
    if !adjoint_consistent {
        findings.push(if suite.adjoint_pairs.is_empty() {
            "adjoint consistency: no witness pairs supplied".to_string()
        } else {
            "adjoint consistency: <Ax,y> != <x,Aᵀy> (declared transpose is not the adjoint)"
                .to_string()
        });
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

    ConformanceReport {
        converter: c.id().to_string(),
        functoriality,
        adjoint_consistent,
        tolerance_honest,
        measured_error,
        tier,
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

#[cfg(test)]
mod arithmetic_tests {
    use super::*;

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
        assert_eq!(
            squared_norm_le_bound(&[1.0], &[0.0], below_one),
            Some(false)
        );
        assert_eq!(squared_norm_le_bound(&[1.0], &[0.0], above_one), Some(true));

        assert_eq!(
            squared_norm_le_bound(&[f64::MAX], &[0.0], Dd::from_f64(f64::MAX)),
            Some(true)
        );
        assert_eq!(
            squared_norm_le_bound(&[f64::MAX], &[0.0], Dd::from_f64(f64::MAX.next_down())),
            Some(false)
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
}
