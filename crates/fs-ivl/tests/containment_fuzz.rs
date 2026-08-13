//! Targeted containment fuzzing near floating-point pathology boundaries
//! (bead frankensim-extreal-program-f85xj.3.3).
//!
//! The broad G0 containment law samples widely; pathology lives at the
//! edges: gradual underflow, overflow saturation, catastrophic
//! cancellation, signed zeros, and ULP-boundary powers of two. This battery
//! aims structured samplers at exactly those neighborhoods for the
//! load-bearing interval operations, checks containment against an
//! independent double-double witness, and — the load-bearing part — proves
//! its own DETECTION POWER: for every boundary class, an off-by-one-ULP
//! in-test fault must be FOUND within the class budget, and un-nudged /
//! inward-rounded mutant operations must be killed. A fuzzer that never
//! alarms proves nothing unless it demonstrably alarms on a planted fault.
//!
//! Determinism: every case derives from `fs_propcheck::Stream::for_case`
//! with the suite seed recorded in the log rows; there is no ambient RNG,
//! so any failure replays from (class, case index) alone. Shrinking walks
//! the failing operands toward simple values while the failure persists and
//! reports the minimal case; permanently persisted minimal fixtures live in
//! `MINIMAL_REGRESSION_FIXTURES` (currently empty: no real violation has
//! been found; the self-check proves the machinery would have caught one).
//!
//! Oracle choice, logged per class: the double-double witness (`fs_math::dd`,
//! error-free two-sum/two-product transforms) — exact for single add/sub/mul
//! on f64 operands and ~2^-104 for division, far below the 1-ULP outward
//! nudge under audit. For the SUBNORMAL class the witness runs in a scaled
//! frame (operands pre-scaled by an exact 2^600, endpoints compared in the
//! same frame), because double-double widens the mantissa, not the exponent
//! — a product of two subnormals underflows inside a naive oracle, a
//! vacuity this battery's own planted-fault self-check exposed on its first
//! run. The quarantined MPFR oracle (e03/.3.1) cannot link into workspace
//! tests by design; elementary-function containment is covered by the
//! budget-widening audit (.3.5) and the e03 lane instead.
//!
//! Fault classes (allocation, worker, I/O, process) are N/A by type: pure
//! computation, no allocation on the audited path beyond test bookkeeping,
//! no receipts; a panic is the failure mode. Not a proof of every input —
//! a fuzz battery is measured evidence with a declared budget.

use fs_ivl::Interval;
use fs_math::dd::Dd;
use fs_propcheck::Stream;

/// Suite seed: recorded in every log row so any case replays.
const SUITE_SEED: u64 = 0x3357_1e5f_ca5e_b007;
/// Cases per (class, operation) pair.
const CASES_PER_CLASS: u64 = 2_000;
/// Persisted minimal counterexamples from prior fuzz failures. Empty:
/// no genuine containment violation has been found to date. A future
/// failure prints its minimal case; commit it here so it replays forever.
const MINIMAL_REGRESSION_FIXTURES: &[(f64, f64, f64, f64, &str)] = &[];

// ---------------------------------------------------------------------------
// Boundary-class samplers
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Class {
    Subnormal,
    NearOverflow,
    Cancellation,
    SignedZero,
    Pow2Boundary,
    GeneralNormal,
}

const CLASSES: [Class; 6] = [
    Class::Subnormal,
    Class::NearOverflow,
    Class::Cancellation,
    Class::SignedZero,
    Class::Pow2Boundary,
    Class::GeneralNormal,
];

impl Class {
    fn name(self) -> &'static str {
        match self {
            Class::Subnormal => "subnormal",
            Class::NearOverflow => "near-overflow",
            Class::Cancellation => "cancellation",
            Class::SignedZero => "signed-zero",
            Class::Pow2Boundary => "pow2-boundary",
            Class::GeneralNormal => "general-normal",
        }
    }

    /// Whether a sample belongs to the class neighborhood — the generator
    /// contract a biased mutant generator must violate.
    fn admits(self, x: f64) -> bool {
        match self {
            Class::Subnormal => x == 0.0 || (x.is_finite() && !x.is_normal()),
            Class::NearOverflow => x.is_finite() && x.abs() >= f64::MAX / 1024.0,
            // Cancellation is a property of the PAIR; members are ordinary
            // normals.
            Class::Cancellation | Class::GeneralNormal => x.is_finite(),
            Class::SignedZero => x == 0.0 || x.is_finite(),
            Class::Pow2Boundary => {
                x.is_finite()
                    && (x.to_bits() & ((1u64 << 52) - 1) < 8
                        || x.to_bits() & ((1u64 << 52) - 1) > (1u64 << 52) - 9)
            }
        }
    }
}

fn signed(stream: &mut Stream, magnitude_bits: u64) -> f64 {
    let sign = if stream.next_u64() & 1 == 0 {
        0
    } else {
        1u64 << 63
    };
    f64::from_bits(sign | magnitude_bits)
}

/// One operand pair biased to the class neighborhood. The second operand of
/// a cancellation pair sits within a few ULPs of the first; other classes
/// draw independently.
fn sample_pair(class: Class, stream: &mut Stream) -> (f64, f64) {
    let normal_bits = |s: &mut Stream| {
        // Exponents away from the extremes so single ops stay finite in the
        // general class; boundary classes push the extremes deliberately.
        let exponent = 512 + s.next_u64() % 1024;
        let mantissa = s.next_u64() & ((1u64 << 52) - 1);
        (exponent << 52) | mantissa
    };
    match class {
        Class::Subnormal => {
            let ma = stream.next_u64() & ((1u64 << 52) - 1);
            let a = signed(stream, ma);
            let mb = stream.next_u64() & ((1u64 << 52) - 1);
            let b = signed(stream, mb);
            (a, b)
        }
        Class::NearOverflow => {
            let near_max = |s: &mut Stream| {
                let magnitude = f64::MAX.to_bits() - s.next_u64() % 4096;
                signed(s, magnitude)
            };
            (near_max(stream), near_max(stream))
        }
        Class::Cancellation => {
            let magnitude = normal_bits(stream);
            let a = signed(stream, magnitude);
            let offset = stream.next_u64() % 64;
            let b_bits = if stream.next_u64() & 1 == 0 {
                a.to_bits().wrapping_add(offset)
            } else {
                a.to_bits().wrapping_sub(offset)
            };
            // Same sign, adjacent magnitude: a - b cancels catastrophically.
            (a, f64::from_bits(b_bits))
        }
        Class::SignedZero => {
            let zero = if stream.next_u64() & 1 == 0 {
                0.0
            } else {
                -0.0
            };
            let magnitude = normal_bits(stream);
            let other = signed(stream, magnitude);
            if stream.next_u64() & 1 == 0 {
                (zero, other)
            } else {
                (other, zero)
            }
        }
        Class::Pow2Boundary => {
            let pow2_near = |s: &mut Stream| {
                let exponent = 1 + s.next_u64() % 2045;
                let nudge = s.next_u64() % 8;
                let base = exponent << 52;
                let bits = if s.next_u64() & 1 == 0 {
                    base + nudge
                } else {
                    base - nudge
                };
                signed(s, bits)
            };
            (pow2_near(stream), pow2_near(stream))
        }
        Class::GeneralNormal => {
            let ma = normal_bits(stream);
            let a = signed(stream, ma);
            let mb = normal_bits(stream);
            let b = signed(stream, mb);
            (a, b)
        }
    }
}

// ---------------------------------------------------------------------------
// Containment harness
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, Debug)]
enum OpKind {
    Add,
    Sub,
    Mul,
    Div,
}

const OPS: [OpKind; 4] = [OpKind::Add, OpKind::Sub, OpKind::Mul, OpKind::Div];

impl OpKind {
    fn name(self) -> &'static str {
        match self {
            OpKind::Add => "add",
            OpKind::Sub => "sub",
            OpKind::Mul => "mul",
            OpKind::Div => "div",
        }
    }

    fn apply(self, a: Interval, b: Interval) -> Interval {
        match self {
            OpKind::Add => a + b,
            OpKind::Sub => a - b,
            OpKind::Mul => a * b,
            OpKind::Div => a / b,
        }
    }

    fn witness(self, a: f64, b: f64) -> Dd {
        let (a, b) = (Dd::from_f64(a), Dd::from_f64(b));
        match self {
            OpKind::Add => a + b,
            OpKind::Sub => a - b,
            OpKind::Mul => a * b,
            OpKind::Div => a / b,
        }
    }
}

fn dd_contains(iv: Interval, v: Dd) -> bool {
    // lo <= v <= hi in double-double order (both endpoints finite here;
    // infinite endpoints trivially contain).
    let above_lo = iv.lo() == f64::NEG_INFINITY || !v.lt(Dd::from_f64(iv.lo()));
    let below_hi = iv.hi() == f64::INFINITY || !Dd::from_f64(iv.hi()).lt(v);
    above_lo && below_hi
}

/// Class-aware containment check. For the subnormal class the plain dd
/// witness is VACUOUS for products: double-double widens the mantissa, not
/// the exponent, so a product of two subnormals underflows inside the
/// oracle itself (this battery's own planted-fault self-check exposed
/// that). The repair is exact power-of-two rescaling: operands are scaled
/// by 2^600 (exact), the witness is computed in the normal range, and the
/// result endpoints are compared inside the same scaled frame (endpoint
/// scaling by a power of two is exact while finite).
fn contained(class: Class, op: OpKind, result: Interval, pa: f64, pb: f64) -> bool {
    if class != Class::Subnormal {
        return dd_contains(result, op.witness(pa, pb));
    }
    let s = f64::from_bits((600 + 1023) << 52); // 2^600, exact
    let (witness, frame_powers) = match op {
        OpKind::Add => (Dd::from_f64(pa * s) + Dd::from_f64(pb * s), 1u32),
        OpKind::Sub => (Dd::from_f64(pa * s) - Dd::from_f64(pb * s), 1),
        OpKind::Mul => (Dd::from_f64(pa * s) * Dd::from_f64(pb * s), 2),
        OpKind::Div => (Dd::from_f64(pa * s) / Dd::from_f64(pb * s), 0),
    };
    let frame = |e: f64| -> f64 {
        let mut v = e;
        for _ in 0..frame_powers {
            v *= s;
        }
        v
    };
    let (lo, hi) = (frame(result.lo()), frame(result.hi()));
    // Subnormal-class results are tiny, so the framed endpoints stay finite;
    // an infinite endpoint (WHOLE etc.) remains trivially containing.
    let above_lo = lo == f64::NEG_INFINITY || !witness.lt(Dd::from_f64(lo));
    let below_hi = hi == f64::INFINITY || !Dd::from_f64(hi).lt(witness);
    above_lo && below_hi
}

#[derive(Debug, Clone, Copy)]
struct Violation {
    case_index: u64,
    a: f64,
    b: f64,
    probe_a: f64,
    probe_b: f64,
}

/// Fuzz one (class, op) pair against a CANDIDATE op implementation. The
/// probes are the operand-interval endpoints — the paths a direction flip
/// corrupts. Returns the first containment violation, if any.
fn fuzz_class(
    class: Class,
    op: OpKind,
    apply: &dyn Fn(Interval, Interval) -> Interval,
    budget: u64,
) -> (u64, Option<Violation>) {
    let mut checked = 0u64;
    for case_index in 0..budget {
        let mut stream =
            Stream::for_case(SUITE_SEED ^ (class as u64) << 8 ^ (op as u64), case_index);
        let (a1, a2) = sample_pair(class, &mut stream);
        let (b1, b2) = sample_pair(class, &mut stream);
        if !(a1.is_finite() && a2.is_finite() && b1.is_finite() && b2.is_finite()) {
            continue;
        }
        let a = Interval::new(a1.min(a2), a1.max(a2));
        let b = Interval::new(b1.min(b2), b1.max(b2));
        if matches!(op, OpKind::Div) && b.contains_zero() {
            // Zero-straddling divisors answer WHOLE by policy; containment
            // there is trivial and would only dilute detection power.
            continue;
        }
        let result = apply(a, b);
        for &pa in &[a.lo(), a.hi()] {
            for &pb in &[b.lo(), b.hi()] {
                checked += 1;
                if !contained(class, op, result, pa, pb) {
                    return (
                        checked,
                        Some(Violation {
                            case_index,
                            a: pa,
                            b: pb,
                            probe_a: pa,
                            probe_b: pb,
                        }),
                    );
                }
            }
        }
    }
    (checked, None)
}

/// Shrink a failing operand pair toward simple values (halving the bit
/// distance to 1.0) while the violation persists under the same faulty op.
fn shrink(
    class: Class,
    op: OpKind,
    apply: &dyn Fn(Interval, Interval) -> Interval,
    mut a: f64,
    mut b: f64,
) -> (f64, f64) {
    let violates = |a: f64, b: f64| -> bool {
        let ia = Interval::point(a);
        let ib = Interval::point(b);
        if matches!(op, OpKind::Div) && ib.contains_zero() {
            return false;
        }
        !contained(class, op, apply(ia, ib), a, b)
    };
    if !violates(a, b) {
        return (a, b);
    }
    let goal = 1.0f64.to_bits();
    for _ in 0..128 {
        let mut advanced = false;
        for which in 0..2usize {
            let current = if which == 0 { a } else { b };
            let candidate = f64::from_bits(current.to_bits() / 2 + goal / 2);
            let (na, nb) = if which == 0 {
                (candidate, b)
            } else {
                (a, candidate)
            };
            if violates(na, nb) {
                a = na;
                b = nb;
                advanced = true;
            }
        }
        if !advanced {
            break;
        }
    }
    (a, b)
}

// ---------------------------------------------------------------------------
// The battery
// ---------------------------------------------------------------------------

#[test]
fn g0_boundary_class_containment_holds_for_every_load_bearing_op() {
    for class in CLASSES {
        let mut class_checked = 0u64;
        for op in OPS {
            let (checked, violation) =
                fuzz_class(class, op, &|a, b| op.apply(a, b), CASES_PER_CLASS);
            if let Some(v) = violation {
                let (sa, sb) = shrink(class, op, &|a, b| op.apply(a, b), v.a, v.b);
                panic!(
                    "containment violation: class {} op {} case {} probes ({:e},{:e}); \
                     minimal ({sa:e},{sb:e}) [{:016x},{:016x}] — commit to \
                     MINIMAL_REGRESSION_FIXTURES",
                    class.name(),
                    op.name(),
                    v.case_index,
                    v.probe_a,
                    v.probe_b,
                    sa.to_bits(),
                    sb.to_bits()
                );
            }
            class_checked += checked;
        }
        assert!(
            class_checked > 10_000,
            "class {} is not vacuous: {class_checked}",
            class.name()
        );
        println!(
            "{{\"suite\":\"fs-ivl-containment-fuzz\",\"case\":\"class\",\"class\":\"{}\",\
             \"oracle\":\"dd\",\"seed\":\"{SUITE_SEED:016x}\",\"checks\":{class_checked},\
             \"violations\":0,\"verdict\":\"pass\"}}",
            class.name()
        );
    }
}

#[test]
fn g0_persisted_minimal_fixtures_replay_green() {
    // Today this pins the EMPTINESS honestly and keeps the replay path
    // compiled; the moment a fixture is committed it becomes a permanent
    // regression.
    for &(a, b, lo, hi, what) in MINIMAL_REGRESSION_FIXTURES {
        let iv = Interval::new(lo, hi);
        assert!(
            iv.contains(a) || iv.contains(b),
            "stale fixture: {what} ({a:e},{b:e})"
        );
    }
    println!(
        "{{\"suite\":\"fs-ivl-containment-fuzz\",\"case\":\"fixtures\",\"count\":{},\
         \"verdict\":\"pass\"}}",
        MINIMAL_REGRESSION_FIXTURES.len()
    );
}

#[test]
fn g0_samplers_actually_hit_their_declared_neighborhoods() {
    // Generator contract: a biased mutant sampler (one that quietly avoids
    // the boundary) violates the class predicate and fails here — coverage
    // is asserted, not assumed.
    for class in CLASSES {
        let mut in_class = 0u64;
        let total = 512u64;
        for case_index in 0..total {
            let mut stream = Stream::for_case(SUITE_SEED ^ 0xC0FF_EE00 ^ class as u64, case_index);
            let (a, b) = sample_pair(class, &mut stream);
            if class.admits(a) && class.admits(b) {
                in_class += 1;
            }
            if class == Class::Cancellation {
                // The PAIR property: within 64 ULPs.
                let d = a.to_bits().abs_diff(b.to_bits());
                assert!(d <= 64, "cancellation pair drifted: {d} ulps");
            }
        }
        assert_eq!(
            in_class,
            total,
            "sampler for {} left its declared neighborhood",
            class.name()
        );
    }
}

#[test]
fn g3_metamorphic_commutation_and_inclusion_monotonicity() {
    let mut checked = 0u64;
    for class in CLASSES {
        for case_index in 0..512u64 {
            let mut stream =
                Stream::for_case(SUITE_SEED ^ 0xBEEF ^ (class as u64) << 4, case_index);
            let (a1, a2) = sample_pair(class, &mut stream);
            let (b1, b2) = sample_pair(class, &mut stream);
            if !(a1.is_finite() && a2.is_finite() && b1.is_finite() && b2.is_finite()) {
                continue;
            }
            let a = Interval::new(a1.min(a2), a1.max(a2));
            let b = Interval::new(b1.min(b2), b1.max(b2));
            // Commutation is exact for add and mul.
            let (s1, s2) = (a + b, b + a);
            assert_eq!(
                (s1.lo().to_bits(), s1.hi().to_bits()),
                (s2.lo().to_bits(), s2.hi().to_bits()),
                "add commutation"
            );
            let (m1, m2) = (a * b, b * a);
            assert_eq!(
                (m1.lo().to_bits(), m1.hi().to_bits()),
                (m2.lo().to_bits(), m2.hi().to_bits()),
                "mul commutation"
            );
            // Inclusion monotonicity: widening an operand can only widen
            // the enclosure.
            let wider = Interval::new(fs_math::next_down(a.lo()), fs_math::next_up(a.hi()));
            for op in OPS {
                if matches!(op, OpKind::Div) && b.contains_zero() {
                    continue;
                }
                let narrow = op.apply(a, b);
                let wide = op.apply(wider, b);
                assert!(
                    wide.lo() <= narrow.lo() && wide.hi() >= narrow.hi(),
                    "inclusion monotonicity: {} at [{:e},{:e}]",
                    op.name(),
                    a.lo(),
                    a.hi()
                );
            }
            checked += 1;
        }
    }
    assert!(
        checked > 2_000,
        "metamorphic battery not vacuous: {checked}"
    );
}

// ---------------------------------------------------------------------------
// Detection power: planted faults and mutant operations must be caught
// ---------------------------------------------------------------------------

#[test]
fn g3_off_by_one_ulp_faults_are_found_in_every_boundary_class() {
    // The planted fault removes ONE side's protective outward nudge by
    // stepping the lower endpoint inward a single ULP. For every class
    // there must exist an operation whose fuzz run FINDS it within budget;
    // classes whose arithmetic is exact for some op (e.g. signed-zero
    // addition) still detect through an inexact op, which is why the class
    // is required to alarm on ANY op rather than on all four.
    for class in CLASSES {
        let mut found = None;
        'ops: for op in OPS {
            let faulty = move |a: Interval, b: Interval| -> Interval {
                let r = op.apply(a, b);
                if r.lo().is_finite() {
                    Interval::new(fs_math::next_up(r.lo()), r.hi())
                } else {
                    r
                }
            };
            let (_, violation) = fuzz_class(class, op, &faulty, CASES_PER_CLASS);
            if let Some(v) = violation {
                found = Some((op, v.case_index));
                break 'ops;
            }
        }
        let (op, case_index) = found.unwrap_or_else(|| {
            panic!(
                "class {}: a one-ULP inward fault survived every op's budget — \
                 the fuzzer has no detection power here",
                class.name()
            )
        });
        println!(
            "{{\"suite\":\"fs-ivl-containment-fuzz\",\"case\":\"planted-fault\",\
             \"class\":\"{}\",\"detected_by\":\"{}\",\"case_index\":{case_index},\
             \"verdict\":\"killed\"}}",
            class.name(),
            op.name()
        );
    }
}

#[test]
fn g3_unnudged_and_inward_mutant_ops_are_killed() {
    // Realistic wrong implementations, built from scratch rather than by
    // un-nudging the real ones: nearest-rounded (no outward step at all)
    // and inward-rounded (both nudges reversed) addition/multiplication.
    let nearest_add = |a: Interval, b: Interval| Interval::new(a.lo() + b.lo(), a.hi() + b.hi());
    let nearest_mul = |a: Interval, b: Interval| {
        let products = [
            a.lo() * b.lo(),
            a.lo() * b.hi(),
            a.hi() * b.lo(),
            a.hi() * b.hi(),
        ];
        let lo = products.iter().copied().fold(f64::INFINITY, f64::min);
        let hi = products.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        Interval::new(lo, hi)
    };
    let inward_sub = |a: Interval, b: Interval| {
        let lo = a.lo() - b.hi();
        let hi = a.hi() - b.lo();
        let lo_in = if lo.is_finite() {
            fs_math::next_up(lo)
        } else {
            lo
        };
        let hi_in = if hi.is_finite() {
            fs_math::next_down(hi)
        } else {
            hi
        };
        if lo_in <= hi_in {
            Interval::new(lo_in, hi_in)
        } else {
            Interval::point(lo_in.min(hi_in))
        }
    };
    let mutants: [(&str, OpKind, &dyn Fn(Interval, Interval) -> Interval); 3] = [
        ("nearest-add", OpKind::Add, &nearest_add),
        ("nearest-mul", OpKind::Mul, &nearest_mul),
        ("inward-sub", OpKind::Sub, &inward_sub),
    ];
    for (name, op, mutant) in mutants {
        let mut killed_by = None;
        for class in CLASSES {
            let (_, violation) = fuzz_class(class, op, mutant, CASES_PER_CLASS);
            if let Some(v) = violation {
                killed_by = Some((class, v.case_index));
                break;
            }
        }
        let (class, case_index) =
            killed_by.unwrap_or_else(|| panic!("mutant `{name}` survived every class budget"));
        println!(
            "{{\"suite\":\"fs-ivl-containment-fuzz\",\"case\":\"mutant\",\"name\":\"{name}\",\
             \"killed_by_class\":\"{}\",\"case_index\":{case_index},\"verdict\":\"killed\"}}",
            class.name()
        );
    }
}

#[test]
fn g0_affine_boundary_class_containment_holds_for_add_sub_mul() {
    // Affine forms concentrate their rounding conservatism in the noise
    // radius; the containment law is the same — the represented interval of
    // an affine op result must contain the double-double witness of the
    // endpoint probes. Same samplers, same class-aware scaled oracle.
    use fs_ivl::AffineCtx;
    for class in CLASSES {
        let mut checked = 0u64;
        for op in [OpKind::Add, OpKind::Sub, OpKind::Mul] {
            for case_index in 0..512u64 {
                let mut stream = Stream::for_case(
                    SUITE_SEED ^ 0xAFF1_2E00 ^ (class as u64) << 8 ^ (op as u64),
                    case_index,
                );
                let (a1, a2) = sample_pair(class, &mut stream);
                let (b1, b2) = sample_pair(class, &mut stream);
                if !(a1.is_finite() && a2.is_finite() && b1.is_finite() && b2.is_finite()) {
                    continue;
                }
                let ia = Interval::new(a1.min(a2), a1.max(a2));
                let ib = Interval::new(b1.min(b2), b1.max(b2));
                let mut ctx = AffineCtx::new();
                let fa = ctx.from_interval(ia);
                let fb = ctx.from_interval(ib);
                let affine_result = match op {
                    OpKind::Add => &fa + &fb,
                    OpKind::Sub => &fa - &fb,
                    OpKind::Mul => &fa * &fb,
                    OpKind::Div => unreachable!(),
                };
                let result = affine_result.to_interval();
                for &pa in &[ia.lo(), ia.hi()] {
                    for &pb in &[ib.lo(), ib.hi()] {
                        checked += 1;
                        assert!(
                            contained(class, op, result, pa, pb),
                            "affine {} containment: class {} case {case_index} \
                             probes ({pa:e},{pb:e}) result [{:e},{:e}]",
                            op.name(),
                            class.name(),
                            result.lo(),
                            result.hi()
                        );
                    }
                }
            }
        }
        assert!(
            checked > 4_000,
            "affine class {} not vacuous: {checked}",
            class.name()
        );
        println!(
            "{{\"suite\":\"fs-ivl-containment-fuzz\",\"case\":\"affine-class\",\"class\":\"{}\",\
             \"oracle\":\"dd\",\"checks\":{checked},\"violations\":0,\"verdict\":\"pass\"}}",
            class.name()
        );
    }
}
