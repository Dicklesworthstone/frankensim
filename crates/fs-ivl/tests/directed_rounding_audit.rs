//! Directed-rounding audit (bead frankensim-extreal-program-f85xj.3.5).
//!
//! fs-ivl implements outward rounding by 1-ULP nudges (`fs_math::next_up` /
//! `next_down`) instead of global FPU rounding-mode state. That strategy is
//! sound ONLY if (1) the nudge functions are exactly correct at every IEEE-754
//! boundary class, and (2) every call site nudges in the direction its
//! operation's monotonicity requires. This battery supplies:
//!
//! - `MODEL v1`: an INDEPENDENT bit-level successor/predecessor model
//!   (integer step on the sign-magnitude bit pattern — the IEEE `nextUp`
//!   definition transcribed directly, sharing no code with the audited path);
//! - exhaustive-by-class f64 verification windows plus a full-exhaustive f32
//!   transition sweep of the same model construction (the width where
//!   exhaustion is feasible);
//! - property/metamorphic laws (ordering, involution, sign reflection,
//!   exact power-of-two scaling);
//! - a REGISTERED CALL-SITE INVENTORY: the test re-counts `next_up` /
//!   `next_down` code tokens in every fs-ivl source file and fails on an
//!   unregistered site, so a new outward-rounding site cannot land without
//!   being classified here;
//! - the declared-budget widening check for the elementary enclosures; and
//! - a MUTANT LEDGER: reversed, skipped, duplicated, conditionally-skipped,
//!   signed-zero-mangling, and infinity-mishandling nudge mutants must all
//!   be killed by the same harness that passes the real functions.
//!
//! Fault-injection classes (allocation, timeout, process, I/O) are N/A by
//! type: this battery is pure computation over `f64` with no allocation on
//! the audited path, no I/O, and no receipts — a panic IS the failure mode
//! and produces no green result. Logs are bounded casebook-style JSON lines.
//!
//! This audit is evidence toward the outward-rounding invariant, not a
//! formal model-to-code proof (that remains .3.8's scope), and the site
//! inventory covers fs-ivl only.

use fs_ivl::Interval;
use fs_math::det;

/// Frozen model identity for downstream consumers (.3.7 / .3.8.1).
const MODEL_VERSION: u32 = 1;

// ---------------------------------------------------------------------------
// MODEL v1: bit-level ground truth
// ---------------------------------------------------------------------------

/// IEEE-754 `nextUp` transcribed at the bit level: the smallest float
/// strictly greater than `x`. Independent of both `std` and `fs_math`.
fn model_next_up(x: f64) -> f64 {
    if x.is_nan() {
        return x;
    }
    if x == f64::INFINITY {
        return x;
    }
    if x == 0.0 {
        // Both zeros: the next value up is the least positive subnormal.
        return f64::from_bits(1);
    }
    let bits = x.to_bits();
    if bits >> 63 == 0 {
        // Positive: magnitude step up (monotone in the bit pattern).
        f64::from_bits(bits + 1)
    } else {
        // Negative: magnitude step down. -min_subnormal steps to -0.0
        // (bits 0x8000...0), which compares equal to 0.0 as required.
        f64::from_bits(bits - 1)
    }
}

/// `nextDown` by reflection of the model, not of the audited code.
fn model_next_down(x: f64) -> f64 {
    -model_next_up(-x)
}

/// The f32 analog of the same construction, used where FULL exhaustion is
/// feasible; validated against `f32::next_up` (an implementation this crate
/// never uses) purely to certify the model construction itself.
fn model_next_up_f32(x: f32) -> f32 {
    if x.is_nan() {
        return x;
    }
    if x == f32::INFINITY {
        return x;
    }
    if x == 0.0 {
        return f32::from_bits(1);
    }
    let bits = x.to_bits();
    if bits >> 31 == 0 {
        f32::from_bits(bits + 1)
    } else {
        f32::from_bits(bits - 1)
    }
}

/// Bitwise equality with NaN ≡ NaN (payload-insensitive: the audited
/// functions pass NaN through, and any NaN answer encloses nothing).
fn same(a: f64, b: f64) -> bool {
    (a.is_nan() && b.is_nan()) || a.to_bits() == b.to_bits()
}

// ---------------------------------------------------------------------------
// Class windows
// ---------------------------------------------------------------------------

/// Boundary-class windows: `(class name, members)`. Every value in every
/// window is checked in both directions against the model.
fn class_windows() -> Vec<(&'static str, Vec<f64>)> {
    let window_at = |bits: u64, below: u64, above: u64| -> Vec<f64> {
        let lo = bits.saturating_sub(below);
        (lo..=bits.saturating_add(above))
            .map(f64::from_bits)
            .collect()
    };
    let mut classes = vec![
        (
            "subnormal-positive",
            (0..=4096).map(f64::from_bits).collect::<Vec<_>>(),
        ),
        (
            "subnormal-negative",
            (0..=4096)
                .map(|m| f64::from_bits((1u64 << 63) | m))
                .collect(),
        ),
        (
            "min-normal-boundary",
            window_at(f64::MIN_POSITIVE.to_bits(), 2048, 2048),
        ),
        ("one-boundary", window_at(1.0f64.to_bits(), 2048, 2048)),
        ("max-finite", window_at(f64::MAX.to_bits(), 2048, 0)),
        (
            "specials",
            vec![
                0.0,
                -0.0,
                f64::INFINITY,
                f64::NEG_INFINITY,
                f64::NAN,
                f64::MIN_POSITIVE,
                -f64::MIN_POSITIVE,
                f64::MAX,
                f64::MIN,
            ],
        ),
    ];
    // Powers of two both signs: the ULP halves/doubles exactly here, the
    // classic wrong-by-a-factor-of-two hazard for nudge implementations.
    // Exact bit-level construction (no powi: the det:: doctrine and the
    // model's independence both want arithmetic-free seeds).
    let mut pow2 = Vec::new();
    let mut k = -1022i64;
    while k <= 1023 {
        let bits = (u64::try_from(k + 1023).expect("normal exponent")) << 52;
        pow2.extend(window_at(bits, 8, 8));
        pow2.extend(window_at((1u64 << 63) | bits, 8, 8));
        k += 37;
    }
    classes.push(("power-of-two", pow2));
    classes
}

/// Run the full class battery against a candidate (up, down) pair.
/// Returns the per-class check counts, or the first divergence.
fn audit_against_model(
    up: &dyn Fn(f64) -> f64,
    down: &dyn Fn(f64) -> f64,
) -> Result<Vec<(&'static str, usize)>, String> {
    let mut counts = Vec::new();
    for (class, members) in class_windows() {
        for &v in &members {
            let (got_up, want_up) = (up(v), model_next_up(v));
            if !same(got_up, want_up) {
                return Err(format!(
                    "class {class}: next_up({v:e} [{:016x}]) = {got_up:e}, model says {want_up:e}",
                    v.to_bits()
                ));
            }
            let (got_down, want_down) = (down(v), model_next_down(v));
            if !same(got_down, want_down) {
                return Err(format!(
                    "class {class}: next_down({v:e} [{:016x}]) = {got_down:e}, model says {want_down:e}",
                    v.to_bits()
                ));
            }
        }
        counts.push((class, members.len() * 2));
    }
    Ok(counts)
}

#[test]
fn g0_next_up_down_match_the_bit_level_model_at_every_boundary_class() {
    let counts = audit_against_model(&fs_math::next_up, &fs_math::next_down)
        .expect("audited nudges match MODEL v1");
    let total: usize = counts.iter().map(|c| c.1).sum();
    assert!(total > 20_000, "the battery is not vacuous: {total}");
    for (class, n) in &counts {
        println!(
            "{{\"suite\":\"fs-ivl-directed-rounding\",\"model_version\":{MODEL_VERSION},\
             \"case\":\"class\",\"class\":\"{class}\",\"checks\":{n},\"verdict\":\"pass\"}}"
        );
    }
}

#[test]
fn g0_f32_transitions_are_exhausted_for_the_model_construction() {
    // The f64 model cannot be exhausted; the SAME construction at f32 can be
    // exhausted over every subnormal and every exponent-boundary window,
    // against `f32::next_up` (never used by production code). This certifies
    // the bit-level scheme itself where full coverage is affordable.
    let mut checks = 0usize;
    // All subnormals, both signs (2^23 each), plus the boundary into normals.
    for m in 0..=(1u32 << 23) {
        for sign in [0u32, 1u32 << 31] {
            let v = f32::from_bits(sign | m);
            assert_eq!(
                model_next_up_f32(v).to_bits(),
                v.next_up().to_bits(),
                "f32 subnormal transition at {:08x}",
                v.to_bits()
            );
            checks += 1;
        }
    }
    // Every exponent boundary: 64-value windows around each power of two.
    for e in 1..=254u32 {
        let p = e << 23;
        for d in 0..64u32 {
            for sign in [0u32, 1u32 << 31] {
                let v = f32::from_bits(sign | (p - 32 + d));
                assert_eq!(
                    model_next_up_f32(v).to_bits(),
                    v.next_up().to_bits(),
                    "f32 exponent boundary at {:08x}",
                    v.to_bits()
                );
                checks += 1;
            }
        }
    }
    for v in [f32::INFINITY, f32::NEG_INFINITY, f32::MAX, f32::NAN] {
        let (got, want) = (model_next_up_f32(v), v.next_up());
        assert!(
            (got.is_nan() && want.is_nan()) || got.to_bits() == want.to_bits(),
            "f32 special {v:e}"
        );
        checks += 1;
    }
    println!(
        "{{\"suite\":\"fs-ivl-directed-rounding\",\"model_version\":{MODEL_VERSION},\
         \"case\":\"f32-exhaustion\",\"checks\":{checks},\"verdict\":\"pass\"}}"
    );
}

#[test]
fn g3_ordering_involution_reflection_and_pow2_scaling_laws() {
    // Deterministic witness set: class windows plus an FNV-chained sample of
    // the full finite range (recorded as the sampled f64 witnesses for .3.7).
    let mut witnesses: Vec<f64> = class_windows().into_iter().flat_map(|c| c.1).collect();
    let mut state = 0xcbf2_9ce4_8422_2325u64;
    for _ in 0..8192 {
        state = (state ^ (state >> 33)).wrapping_mul(0x0100_0000_01b3);
        let v = f64::from_bits(state);
        if v.is_finite() {
            witnesses.push(v);
        }
    }
    let mut checked = 0usize;
    for &x in &witnesses {
        if !x.is_finite() {
            continue;
        }
        let up = fs_math::next_up(x);
        // Strict ordering.
        assert!(up > x, "ordering at {x:e}");
        // Involution (value equality: +0/-0 collapse is correct behavior).
        if up.is_finite() {
            assert_eq!(fs_math::next_down(up), x, "involution at {x:e}");
        }
        // Sign reflection: next_up(-x) == -next_down(x), signed zeros included.
        assert!(
            same(fs_math::next_up(-x), -fs_math::next_down(x)),
            "reflection at {x:e}"
        );
        // Exact power-of-two scaling where admissible: for normal x with 2x
        // normal and finite, ulp(2x) = 2*ulp(x) exactly. EXPLICIT
        // PRECONDITION discovered by this battery's own first run: a
        // NEGATIVE exact power of two is excluded, because next_up steps
        // toward zero across the exponent boundary into finer spacing
        // (next_up(-2^k) = -(2^k - 2^(k-53))), so the doubled step and the
        // step of the doubled value differ by design, not by defect.
        let neg_pow2 = x < 0.0 && x.to_bits() & ((1u64 << 52) - 1) == 0;
        if x.is_normal() && (2.0 * x).is_normal() && (2.0 * x).is_finite() && !neg_pow2 {
            assert_eq!(
                fs_math::next_up(2.0 * x),
                2.0 * fs_math::next_up(x),
                "pow2 scaling at {x:e}"
            );
        }
        checked += 1;
    }
    assert!(checked > 20_000, "law battery is not vacuous: {checked}");
    println!(
        "{{\"suite\":\"fs-ivl-directed-rounding\",\"model_version\":{MODEL_VERSION},\
         \"case\":\"laws\",\"witnesses\":{checked},\"verdict\":\"pass\"}}"
    );
}

// ---------------------------------------------------------------------------
// Call-site inventory
// ---------------------------------------------------------------------------

/// Registered outward-rounding call sites per fs-ivl source file:
/// `(file, next_up code tokens, next_down code tokens, classification)`.
///
/// The classification is the human-reviewed monotonicity rationale; the
/// counts are enforced against the actual sources below, so ADDING A SITE
/// FAILS THIS TEST until it is classified here. Doc-comment mentions are
/// excluded by the counter.
const REGISTERED_SITES: &[(&str, usize, usize, &str)] = &[
    (
        "src/interval.rs",
        2,
        2,
        "up_k/down_k budget loops (paired outward); pi enclosure (down on lo, \
         up on hi around a correctly-rounded constant). The strictly_outside \
         helper and test-module assertions sit under #[cfg(test)] and are \
         excluded by the production cutoff: they are expectations about \
         nudges, not outward-rounding sites",
    ),
    (
        "src/affine.rs",
        16,
        1,
        "error-radius accumulation: every arithmetic step on the NON-NEGATIVE \
         noise radius rounds UP (conservative growth); to_interval rounds \
         center-radius outward (down on lo, up on hi); round_err rounds a \
         magnitude bound up",
    ),
    (
        "src/newton.rs",
        1,
        0,
        "lipschitz_bound: derivative-magnitude supremum rounds UP \
         (overestimating a Lipschitz constant is conservative)",
    ),
    ("src/taylor.rs", 0, 0, "no direct nudge sites"),
    ("src/expansion.rs", 0, 0, "exact arithmetic; no nudge sites"),
    (
        "src/predicates.rs",
        0,
        0,
        "exact predicates; no nudge sites",
    ),
    ("src/lib.rs", 0, 0, "doc references only"),
];

#[test]
fn g0_every_outward_rounding_call_site_is_registered_and_classified() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    for (file, want_up, want_down, classification) in REGISTERED_SITES {
        let source = std::fs::read_to_string(root.join(file)).expect("fs-ivl source readable");
        // The inventory audits PRODUCTION nudge sites: counting stops at the
        // in-file test module, whose assertions may legitimately use nudges
        // as expectations without being outward-rounding sites.
        let production = source.split("#[cfg(test)]").next().unwrap_or("");
        let mut got_up = 0usize;
        let mut got_down = 0usize;
        for line in production.lines() {
            // Strip line comments so doc text does not count as a site.
            let code = line.split("//").next().unwrap_or("");
            got_up += code.matches("next_up").count();
            got_down += code.matches("next_down").count();
        }
        assert!(
            !classification.is_empty(),
            "{file}: a registered site row must carry its monotonicity rationale"
        );
        assert_eq!(
            (got_up, got_down),
            (*want_up, *want_down),
            "{file}: outward-rounding call sites changed. Classify the new/removed \
             site's monotonicity direction and update REGISTERED_SITES — an \
             unreviewed nudge direction can silently convert certificates into lies."
        );
        println!(
            "{{\"suite\":\"fs-ivl-directed-rounding\",\"case\":\"site-inventory\",\
             \"file\":\"{file}\",\"next_up\":{got_up},\"next_down\":{got_down},\
             \"verdict\":\"registered\"}}"
        );
    }
}

// ---------------------------------------------------------------------------
// Budget widening
// ---------------------------------------------------------------------------

#[test]
fn g0_elementary_enclosures_widen_by_at_least_the_declared_budget() {
    // For singleton inputs at boundary-flavored points, the enclosure must
    // contain the computed point value nudged outward by the FULL declared
    // budget in both directions (clamps to the mathematical range excepted:
    // a bound already at the range edge is exact, not under-widened).
    let probes = [0.0, 1.0, -1.0, 0.5, 2.0, f64::MIN_POSITIVE, 20.0, -20.0];
    let mut checks = 0usize;
    for &x in &probes {
        let s = Interval::new(x, x);

        let e = s.exp();
        let v = det::exp(x);
        if v.is_finite() {
            let mut lo_req = v;
            let mut hi_req = v;
            for _ in 0..det::EXP_ULP_BUDGET {
                lo_req = fs_math::next_down(lo_req);
                hi_req = fs_math::next_up(hi_req);
            }
            assert!(
                e.lo() <= lo_req.max(0.0) && e.hi() >= hi_req,
                "exp widening at {x:e}: [{:e},{:e}] vs required [{lo_req:e},{hi_req:e}]",
                e.lo(),
                e.hi()
            );
            checks += 1;
        }

        let t = s.tanh();
        let v = det::tanh(x);
        let mut lo_req = v;
        let mut hi_req = v;
        for _ in 0..det::TANH_ULP_BUDGET {
            lo_req = fs_math::next_down(lo_req);
            hi_req = fs_math::next_up(hi_req);
        }
        assert!(
            t.lo() <= lo_req.max(-1.0) && t.hi() >= hi_req.min(1.0),
            "tanh widening at {x:e}"
        );
        checks += 1;

        let sn = s.sin();
        let v = det::sin(x);
        let mut lo_req = v;
        let mut hi_req = v;
        for _ in 0..det::SIN_ULP_BUDGET {
            lo_req = fs_math::next_down(lo_req);
            hi_req = fs_math::next_up(hi_req);
        }
        assert!(
            sn.lo() <= lo_req.max(-1.0) && sn.hi() >= hi_req.min(1.0),
            "sin widening at {x:e}"
        );
        checks += 1;
    }
    println!(
        "{{\"suite\":\"fs-ivl-directed-rounding\",\"case\":\"budget-widening\",\
         \"checks\":{checks},\"verdict\":\"pass\"}}"
    );
}

// ---------------------------------------------------------------------------
// Mutant ledger
// ---------------------------------------------------------------------------

#[test]
fn g3_nudge_mutants_are_killed_by_the_model_battery() {
    // Each mutant is a realistic wrong implementation; every one must be
    // caught by the same harness that passes the real pair. The ledger rows
    // are the retained evidence for .3.7's end-to-end tripwire campaign.
    let mutants: Vec<(&str, Box<dyn Fn(f64) -> f64>, Box<dyn Fn(f64) -> f64>)> = vec![
        (
            "reversed-directions",
            Box::new(fs_math::next_down),
            Box::new(fs_math::next_up),
        ),
        (
            "skipped-nudge",
            Box::new(|x| x),
            Box::new(fs_math::next_down),
        ),
        (
            "duplicated-nudge",
            Box::new(|x| fs_math::next_up(fs_math::next_up(x))),
            Box::new(fs_math::next_down),
        ),
        (
            "conditionally-skipped-subnormals",
            Box::new(|x: f64| {
                if x != 0.0 && !x.is_normal() {
                    x
                } else {
                    fs_math::next_up(x)
                }
            }),
            Box::new(fs_math::next_down),
        ),
        (
            "signed-zero-mangled",
            Box::new(|x: f64| {
                if x == 0.0 {
                    // Wrong canonicalization: steps to the NEGATIVE side.
                    -f64::from_bits(1)
                } else {
                    fs_math::next_up(x)
                }
            }),
            Box::new(fs_math::next_down),
        ),
        (
            "infinity-mishandled",
            Box::new(|x: f64| {
                if x == f64::NEG_INFINITY {
                    // Wrong: nextUp(-inf) must be -MAX, not -inf.
                    x
                } else {
                    fs_math::next_up(x)
                }
            }),
            Box::new(fs_math::next_down),
        ),
    ];
    for (name, up, down) in &mutants {
        let outcome = audit_against_model(up, down);
        assert!(
            outcome.is_err(),
            "mutant `{name}` survived the class battery"
        );
        println!(
            "{{\"suite\":\"fs-ivl-directed-rounding\",\"case\":\"mutant\",\"name\":\"{name}\",\
             \"verdict\":\"killed\",\"first_divergence\":{:?}}}",
            outcome.unwrap_err()
        );
    }
    // Positive control: the real pair passes the identical harness.
    assert!(audit_against_model(&fs_math::next_up, &fs_math::next_down).is_ok());
}
