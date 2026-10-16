//! G0/G3/G5 full-trajectory replay for the public Strong-Wolfe line search
//! (7tv.21.50).
//!
//! Five closed-form fixtures with dyadic callback tables retain every
//! callback's alpha, value, and derivative bits together with every public
//! `WolfeOutcome` field. An
//! independent expanded-form oracle derives the complete expected callback
//! tables, evaluation accounting, Wolfe predicates, accepted-callback
//! correspondence, and failure projection. A `StreamKey`-selected one-bit
//! corruption of a successful returned alpha is refused while stale, after
//! resealing against the retained reference, and by the semantic gate.
//!
//! This is fixed-input, same-process replay evidence for these five curves.
//! It makes no claim about arbitrary-function convergence or existence,
//! optimizer quality, Riemannian semantics, public hard-budget or `Cx`
//! behavior, non-finite containment, cross-process or cross-ISA authority,
//! persistence or authentication, or performance.

use fs_ascent::{WolfeOutcome, strong_wolfe};
use fs_obs::ident::{IDENT_SCHEMA_VERSION, IdentityBuilder, ReplayIdentity, check_version};
use fs_obs::{Emitter, Event, EventKind, Severity};
use fs_rand::StreamKey;
use std::panic::catch_unwind;

const SUITE: &str = "fs-ascent/strong-wolfe-study-replay";
const CASE: &str = "five-closed-form-full-trajectories";
const RED_CASE: &str = "seeded-success-alpha-one-bit-corruption";

const FIXTURE_IDENTITY_KIND: &str = "fs-ascent-strong-wolfe-study-fixture-v1";
const RESULT_IDENTITY_KIND: &str = "fs-ascent-strong-wolfe-study-result-v1";
const FIXTURE_SCHEMA_VERSION: u32 = 1;
const RESULT_SCHEMA_VERSION: u32 = 1;

const REPLAY_SEED: u64 = 0x57_01_FE_21_50;
const REPLAY_KERNEL: u32 = 0x5750;
const REPLAY_TILE: u32 = 0;

const MUTATION_SEED: u64 = 0xBA_DA_55_21_50;
const MUTATION_KERNEL: u32 = 0xD150;
const MUTATION_TILE: u32 = 1;
const MUTATION_MANTISSA_BITS: [u32; 8] = [0, 1, 2, 3, 4, 5, 6, 7];

const EXPANSION_CAP: usize = 20;
const ZOOM_CAP: usize = 40;
const TINY_CENTER: f64 = 2.842_170_943_040_401e-14; // 2^-45.
const TINY_CENTER_SQUARED: f64 = 8.077_935_669_463_161e-28; // 2^-90.

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Curve {
    UnitWell,
    ThirteenFourthsWell,
    DescendingLine,
    TinyCenterWell,
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct CaseSpec {
    id: &'static str,
    curve: Curve,
    f0: f64,
    dphi0: f64,
    alpha_init: f64,
    c1: f64,
    c2: f64,
}

const CASES: [CaseSpec; 5] = [
    CaseSpec {
        id: "immediate-acceptance",
        curve: Curve::UnitWell,
        f0: 1.0,
        dphi0: -2.0,
        alpha_init: 1.0,
        c1: 0.25,
        c2: 0.5,
    },
    CaseSpec {
        id: "armijo-triggered-zoom",
        curve: Curve::UnitWell,
        f0: 1.0,
        dphi0: -2.0,
        alpha_init: 4.0,
        c1: 0.25,
        c2: 0.5,
    },
    CaseSpec {
        id: "expansion-then-derivative-sign-zoom",
        curve: Curve::ThirteenFourthsWell,
        f0: 10.5625, // 169/16.
        dphi0: -6.5,
        alpha_init: 0.5,
        c1: 0.01,
        c2: 0.1,
    },
    CaseSpec {
        id: "deterministic-20-step-expansion-cap",
        curve: Curve::DescendingLine,
        f0: 0.0,
        dphi0: -1.0,
        alpha_init: 1.0,
        c1: 0.25,
        c2: 0.5,
    },
    CaseSpec {
        id: "deterministic-40-step-zoom-cap",
        curve: Curve::TinyCenterWell,
        f0: TINY_CENTER_SQUARED,
        dphi0: -5.684_341_886_080_802e-14, // -2^-44.
        alpha_init: 1.0,
        c1: 0.25,
        c2: 0.5,
    },
];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct CallbackBits {
    alpha: u64,
    value: u64,
    derivative: u64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct OutcomeBits {
    alpha: u64,
    f_new: u64,
    evals: usize,
    success: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct CaseRecord {
    id: &'static str,
    callbacks: Vec<CallbackBits>,
    outcome: OutcomeBits,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct StudyRecord {
    cases: Vec<CaseRecord>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct StudyRun {
    replay_key: StreamKey,
    fixture: ReplayIdentity,
    record: StudyRecord,
    result: ReplayIdentity,
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum AdmissionError {
    PayloadIdentityMismatch { declared: u64, computed: u64 },
    ReferenceIdentityMismatch { expected: u64, found: u64 },
    SemanticInconsistency(Vec<String>),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Mutation {
    seed: u64,
    kernel: u32,
    tile: u32,
    selector_draws: u64,
    successful_case_slot: usize,
    case_index: usize,
    mantissa_bit_slot: usize,
    mantissa_bit: u32,
    before: u64,
    after: u64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct SeededCorruption {
    run: StudyRun,
    mutation: Mutation,
    stale_error: AdmissionError,
    reference_error: AdmissionError,
    semantic_error: AdmissionError,
}

fn usize_u64(value: usize) -> u64 {
    u64::try_from(value).expect("fixed fixture cardinality fits u64")
}

fn curve_name(curve: Curve) -> &'static str {
    match curve {
        Curve::UnitWell => "unit-well",
        Curve::ThirteenFourthsWell => "thirteen-fourths-well",
        Curve::DescendingLine => "descending-line",
        Curve::TinyCenterWell => "two-to-minus-forty-five-center-well",
    }
}

fn producer_definition(curve: Curve) -> &'static str {
    match curve {
        Curve::UnitWell => "r=a-1;phi=r*r;dphi=2*r",
        Curve::ThirteenFourthsWell => "r=a-13/4;phi=r*r;dphi=2*r",
        Curve::DescendingLine => "phi=-a;dphi=-1",
        Curve::TinyCenterWell => "r=a-2^-45;phi=r*r;dphi=2*r",
    }
}

fn oracle_definition(curve: Curve) -> &'static str {
    match curve {
        Curve::UnitWell => "phi=a*a-2*a+1;dphi=2*a-2",
        Curve::ThirteenFourthsWell => "phi=a*a-(13/2)*a+169/16;dphi=2*a-13/2",
        Curve::DescendingLine => "phi=0-a;dphi=0-1",
        Curve::TinyCenterWell => "phi=a*a-2*(2^-45)*a+2^-90;dphi=2*a-2^-44",
    }
}

/// Formula passed to the production line search.
fn evaluate_producer(curve: Curve, alpha: f64) -> (f64, f64) {
    match curve {
        Curve::UnitWell => {
            let residual = alpha - 1.0;
            (residual * residual, 2.0 * residual)
        }
        Curve::ThirteenFourthsWell => {
            let residual = alpha - 3.25;
            (residual * residual, 2.0 * residual)
        }
        Curve::DescendingLine => (-alpha, -1.0),
        Curve::TinyCenterWell => {
            let residual = alpha - TINY_CENTER;
            (residual * residual, 2.0 * residual)
        }
    }
}

/// Algebraically independent expanded/factored forms used by admission only.
fn evaluate_oracle(curve: Curve, alpha: f64) -> (f64, f64) {
    match curve {
        Curve::UnitWell => (alpha * alpha - 2.0 * alpha + 1.0, 2.0 * alpha - 2.0),
        Curve::ThirteenFourthsWell => (alpha * alpha - 6.5 * alpha + 10.5625, 2.0 * alpha - 6.5),
        Curve::DescendingLine => (0.0 - alpha, 0.0 - 1.0),
        Curve::TinyCenterWell => (
            alpha * alpha - (TINY_CENTER + TINY_CENTER) * alpha + TINY_CENTER_SQUARED,
            2.0 * alpha - (TINY_CENTER + TINY_CENTER),
        ),
    }
}

/// Cancellation-safe analytic evaluation for a returned alpha that need not
/// belong to the exact callback table. This is intentionally separate from
/// both the callback producer and the expanded-form table oracle: near a
/// quadratic's root, `a*a - 2*a + 1` can round to zero even when `(a-1)^2`
/// is nonzero, which would make a low-mantissa corruption invisible.
fn evaluate_returned_alpha(curve: Curve, alpha: f64) -> (f64, f64) {
    match curve {
        Curve::UnitWell => {
            let displacement = alpha - 1.0;
            (displacement * displacement, displacement + displacement)
        }
        Curve::ThirteenFourthsWell => {
            let displacement = alpha - 3.25;
            (displacement * displacement, displacement + displacement)
        }
        Curve::DescendingLine => (-alpha, -1.0),
        Curve::TinyCenterWell => {
            let displacement = alpha - TINY_CENTER;
            (displacement * displacement, displacement + displacement)
        }
    }
}

fn callback_bits(alpha: f64, value: f64, derivative: f64) -> CallbackBits {
    CallbackBits {
        alpha: alpha.to_bits(),
        value: value.to_bits(),
        derivative: derivative.to_bits(),
    }
}

fn outcome_bits(outcome: &WolfeOutcome) -> OutcomeBits {
    OutcomeBits {
        alpha: outcome.alpha.to_bits(),
        f_new: outcome.f_new.to_bits(),
        evals: outcome.evals,
        success: outcome.success,
    }
}

fn exact_pow2(exponent: i32) -> f64 {
    assert!((-1_022..=1_023).contains(&exponent));
    let biased = u64::try_from(exponent + 1_023).expect("bounded exponent is nonnegative");
    f64::from_bits(biased << 52)
}

/// Exact alpha tables derived from the public state-machine branches and caps.
fn expected_alphas(case_index: usize) -> Vec<f64> {
    match case_index {
        0 => vec![1.0],
        1 => vec![4.0, 2.0, 1.0],
        2 => vec![0.5, 1.0, 2.0, 4.0, 3.0],
        3 => (0..EXPANSION_CAP)
            .map(|exponent| {
                exact_pow2(i32::try_from(exponent).expect("small expansion exponent fits i32"))
            })
            .collect(),
        4 => {
            let mut alphas = Vec::with_capacity(ZOOM_CAP + 1);
            alphas.push(1.0);
            alphas.extend((1..=ZOOM_CAP).map(|depth| {
                exact_pow2(-i32::try_from(depth).expect("small zoom depth fits i32"))
            }));
            alphas
        }
        _ => panic!("unknown fixed case index {case_index}"),
    }
}

fn expected_callbacks(case_index: usize, spec: CaseSpec) -> Vec<CallbackBits> {
    expected_alphas(case_index)
        .into_iter()
        .map(|alpha| {
            let (value, derivative) = evaluate_oracle(spec.curve, alpha);
            callback_bits(alpha, value, derivative)
        })
        .collect()
}

fn expected_outcome(case_index: usize, spec: CaseSpec) -> OutcomeBits {
    match case_index {
        0 => OutcomeBits {
            alpha: 1.0f64.to_bits(),
            f_new: 0.0f64.to_bits(),
            evals: 1,
            success: true,
        },
        1 => OutcomeBits {
            alpha: 1.0f64.to_bits(),
            f_new: 0.0f64.to_bits(),
            evals: 3,
            success: true,
        },
        2 => OutcomeBits {
            alpha: 3.0f64.to_bits(),
            f_new: 0.0625f64.to_bits(),
            evals: 5,
            success: true,
        },
        3 => OutcomeBits {
            alpha: 0.0f64.to_bits(),
            f_new: spec.f0.to_bits(),
            evals: EXPANSION_CAP,
            success: false,
        },
        4 => OutcomeBits {
            alpha: 0.0f64.to_bits(),
            f_new: spec.f0.to_bits(),
            evals: ZOOM_CAP + 1,
            success: false,
        },
        _ => panic!("unknown fixed case index {case_index}"),
    }
}

fn replay_key() -> StreamKey {
    StreamKey {
        seed: REPLAY_SEED,
        kernel: REPLAY_KERNEL,
        tile: REPLAY_TILE,
    }
}

fn fixture_identity(key: StreamKey) -> ReplayIdentity {
    check_version(IDENT_SCHEMA_VERSION).expect("declared identity schema is current");
    let mut builder = IdentityBuilder::new(FIXTURE_IDENTITY_KIND)
        .u64("fixture-schema-version", u64::from(FIXTURE_SCHEMA_VERSION))
        .u64("result-schema-version", u64::from(RESULT_SCHEMA_VERSION))
        .u64("identity-schema-version", u64::from(IDENT_SCHEMA_VERSION))
        .str("algorithm", "fs_ascent::strong_wolfe")
        .str(
            "state-machine",
            "expand-factor=2;expand-cap=20;zoom=ordered-endpoint-midpoint;zoom-cap=40;width-break=1e-16",
        )
        .str(
            "retained-fields",
            "ordered-callback-alpha-value-derivative-bits;outcome-alpha-f_new-evals-success",
        )
        .str(
            "oracle",
            "independent-exact-dyadic-callback-tables-and-expanded-curve-forms-v1",
        )
        .str("units", "alpha-value-derivative-all-dimensionless")
        .str("algorithm-randomness", "none")
        .str(
            "replay-key-role",
            "identity-bound-same-seed-study-provenance;not-consumed-by-strong-wolfe",
        )
        .u64("replay-seed", key.seed)
        .u64("replay-kernel", u64::from(key.kernel))
        .u64("replay-tile", u64::from(key.tile))
        .str("fs-ascent-version", fs_ascent::VERSION)
        .str("fs-obs-version", fs_obs::VERSION)
        .str("fs-rand-version", fs_rand::VERSION)
        .u64(
            "fs-rand-stream-semantics-version",
            u64::from(fs_rand::STREAM_SEMANTICS_VERSION),
        )
        .u64("mutation-seed", MUTATION_SEED)
        .u64("mutation-kernel", u64::from(MUTATION_KERNEL))
        .u64("mutation-tile", u64::from(MUTATION_TILE))
        .str(
            "mutation-target-policy",
            "StreamKey-select-case-from-ordered-successful-indices-[0,1,2]-then-bit-from-ordered-low-mantissa-bits-[0,1,2,3,4,5,6,7]",
        )
        .str(
            "no-claims",
            "arbitrary-function-convergence-or-existence;optimizer-quality;Riemannian-semantics;public-hard-budget-or-Cx-behavior;non-finite-containment;cross-process-or-ISA-authority;persistence-or-authentication;performance",
        )
        .u64("ordered-case-count", usize_u64(CASES.len()));

    for (case_index, spec) in CASES.into_iter().enumerate() {
        builder = builder
            .u64("case-index", usize_u64(case_index))
            .str("case-id", spec.id)
            .str("curve", curve_name(spec.curve))
            .str("producer-definition", producer_definition(spec.curve))
            .str("oracle-definition", oracle_definition(spec.curve))
            .f64_bits("f0", spec.f0)
            .f64_bits("dphi0", spec.dphi0)
            .f64_bits("alpha-init", spec.alpha_init)
            .f64_bits("c1", spec.c1)
            .f64_bits("c2", spec.c2);
    }
    builder.finish()
}

fn result_identity(fixture: &ReplayIdentity, record: &StudyRecord) -> ReplayIdentity {
    let mut builder = IdentityBuilder::new(RESULT_IDENTITY_KIND)
        .u64("result-schema-version", u64::from(RESULT_SCHEMA_VERSION))
        .child("fixture-root", fixture)
        .bytes("fixture-canonical-bytes", fixture.canonical_bytes())
        .u64("case-count", usize_u64(record.cases.len()));
    for (case_index, case) in record.cases.iter().enumerate() {
        builder = builder
            .u64("case-index", usize_u64(case_index))
            .str("case-id", case.id)
            .u64("callback-count", usize_u64(case.callbacks.len()));
        for (callback_index, callback) in case.callbacks.iter().enumerate() {
            builder = builder
                .u64("callback-index", usize_u64(callback_index))
                .f64_bits("callback-alpha", f64::from_bits(callback.alpha))
                .f64_bits("callback-value", f64::from_bits(callback.value))
                .f64_bits("callback-derivative", f64::from_bits(callback.derivative));
        }
        builder = builder
            .f64_bits("outcome-alpha", f64::from_bits(case.outcome.alpha))
            .f64_bits("outcome-f-new", f64::from_bits(case.outcome.f_new))
            .u64("outcome-evals", usize_u64(case.outcome.evals))
            .flag("outcome-success", case.outcome.success);
    }
    builder.finish()
}

fn run_study(key: StreamKey) -> StudyRun {
    let mut cases = Vec::with_capacity(CASES.len());
    for spec in CASES {
        let mut callbacks = Vec::new();
        let mut callback = |alpha: f64| {
            let (value, derivative) = evaluate_producer(spec.curve, alpha);
            callbacks.push(callback_bits(alpha, value, derivative));
            (value, derivative)
        };
        let outcome = strong_wolfe(
            &mut callback,
            spec.f0,
            spec.dphi0,
            spec.alpha_init,
            spec.c1,
            spec.c2,
        );
        drop(callback);
        cases.push(CaseRecord {
            id: spec.id,
            callbacks,
            outcome: outcome_bits(&outcome),
        });
    }

    let record = StudyRecord { cases };
    let fixture = fixture_identity(key);
    let result = result_identity(&fixture, &record);
    StudyRun {
        replay_key: key,
        fixture,
        record,
        result,
    }
}

fn validate_payload(run: &StudyRun) -> Result<(), AdmissionError> {
    let expected_fixture = fixture_identity(run.replay_key);
    if run.fixture != expected_fixture {
        return Err(AdmissionError::PayloadIdentityMismatch {
            declared: run.fixture.root(),
            computed: expected_fixture.root(),
        });
    }
    let computed_result = result_identity(&run.fixture, &run.record);
    if run.result != computed_result {
        return Err(AdmissionError::PayloadIdentityMismatch {
            declared: run.result.root(),
            computed: computed_result.root(),
        });
    }
    Ok(())
}

fn armijo_holds(spec: CaseSpec, alpha: f64, value: f64) -> bool {
    value <= spec.f0 + spec.c1 * alpha * spec.dphi0
}

fn curvature_holds(spec: CaseSpec, derivative: f64) -> bool {
    derivative.abs() <= spec.c2 * spec.dphi0.abs()
}

fn push_path_specific_violations(
    case_index: usize,
    spec: CaseSpec,
    case: &CaseRecord,
    violations: &mut Vec<String>,
) {
    let callback = |index: usize| {
        let retained = case.callbacks[index];
        (
            f64::from_bits(retained.alpha),
            f64::from_bits(retained.value),
            f64::from_bits(retained.derivative),
        )
    };
    match case_index {
        0 => {
            let (alpha, value, derivative) = callback(0);
            if !armijo_holds(spec, alpha, value) || !curvature_holds(spec, derivative) {
                violations.push("case[0]-immediate-callback-not-strong-wolfe".to_string());
            }
        }
        1 => {
            let (first_alpha, first_value, _) = callback(0);
            let (second_alpha, second_value, _) = callback(1);
            if armijo_holds(spec, first_alpha, first_value)
                || armijo_holds(spec, second_alpha, second_value)
            {
                violations.push("case[1]-callbacks-did-not-retain-armijo-zoom".to_string());
            }
        }
        2 => {
            let expansion_prefix_is_valid = (0..3).all(|index| {
                let (alpha, value, derivative) = callback(index);
                armijo_holds(spec, alpha, value)
                    && derivative < 0.0
                    && !curvature_holds(spec, derivative)
            });
            let (turn_alpha, turn_value, turn_derivative) = callback(3);
            let (_, previous_value, _) = callback(2);
            let (accepted_alpha, accepted_value, accepted_derivative) = callback(4);
            if !expansion_prefix_is_valid
                || !armijo_holds(spec, turn_alpha, turn_value)
                || turn_value >= previous_value
                || turn_derivative <= 0.0
                || curvature_holds(spec, turn_derivative)
                || !armijo_holds(spec, accepted_alpha, accepted_value)
                || !curvature_holds(spec, accepted_derivative)
            {
                violations.push("case[2]-expansion-did-not-enter-derivative-sign-zoom".to_string());
            }
        }
        3 => {
            if case.callbacks.iter().any(|retained| {
                let alpha = f64::from_bits(retained.alpha);
                let value = f64::from_bits(retained.value);
                let derivative = f64::from_bits(retained.derivative);
                !armijo_holds(spec, alpha, value)
                    || derivative >= 0.0
                    || curvature_holds(spec, derivative)
            }) {
                violations.push("case[3]-expansion-cap-precondition-broke".to_string());
            }
        }
        4 => {
            if case.callbacks.iter().any(|retained| {
                armijo_holds(
                    spec,
                    f64::from_bits(retained.alpha),
                    f64::from_bits(retained.value),
                )
            }) {
                violations.push("case[4]-zoom-cap-armijo-failure-broke".to_string());
            }
            let final_width = exact_pow2(-i32::try_from(ZOOM_CAP).expect("cap fits i32"));
            if final_width < 1.0e-16 {
                violations.push("case[4]-zoom-width-break-preempted-cap".to_string());
            }
        }
        _ => violations.push(format!("case[{case_index}]-unexpected-index")),
    }
}

#[allow(clippy::too_many_lines)] // One gate audits every retained callback and outcome field.
fn semantic_violations(record: &StudyRecord) -> Vec<String> {
    let mut violations = Vec::new();
    if record.cases.len() != CASES.len() {
        violations.push(format!(
            "case-count:{}!=expected-{}",
            record.cases.len(),
            CASES.len()
        ));
        return violations;
    }

    for (case_index, (spec, case)) in CASES.into_iter().zip(&record.cases).enumerate() {
        if case.id != spec.id {
            violations.push(format!(
                "case[{case_index}]-id:{}!=expected-{}",
                case.id, spec.id
            ));
        }

        let expected_callbacks = expected_callbacks(case_index, spec);
        if case.callbacks != expected_callbacks {
            violations.push(format!("case[{case_index}]-exact-callback-table-mismatch"));
        }
        if case.outcome != expected_outcome(case_index, spec) {
            violations.push(format!("case[{case_index}]-exact-outcome-table-mismatch"));
        }
        if case.outcome.evals != case.callbacks.len() {
            violations.push(format!(
                "case[{case_index}]-evaluation-accounting:{}!={}",
                case.outcome.evals,
                case.callbacks.len()
            ));
        }

        for (callback_index, callback) in case.callbacks.iter().enumerate() {
            let alpha = f64::from_bits(callback.alpha);
            let (oracle_value, oracle_derivative) = evaluate_oracle(spec.curve, alpha);
            if callback.value != oracle_value.to_bits() {
                violations.push(format!(
                    "case[{case_index}]-callback[{callback_index}]-value-oracle-mismatch"
                ));
            }
            if callback.derivative != oracle_derivative.to_bits() {
                violations.push(format!(
                    "case[{case_index}]-callback[{callback_index}]-derivative-oracle-mismatch"
                ));
            }
        }

        if case.outcome.success {
            let returned_alpha = f64::from_bits(case.outcome.alpha);
            let accepted_callback = case
                .callbacks
                .last()
                .filter(|callback| callback.alpha == case.outcome.alpha);
            match accepted_callback {
                Some(callback) => {
                    if callback.value != case.outcome.f_new {
                        violations.push(format!(
                            "case[{case_index}]-accepted-callback-value!=outcome-f-new"
                        ));
                    }
                }
                None => violations.push(format!(
                    "case[{case_index}]-success-alpha-not-retained-callback"
                )),
            }

            let (oracle_value, oracle_derivative) =
                evaluate_returned_alpha(spec.curve, returned_alpha);
            if oracle_value.to_bits() != case.outcome.f_new {
                violations.push(format!(
                    "case[{case_index}]-returned-value-does-not-match-curve-at-alpha"
                ));
            }
            if !armijo_holds(spec, returned_alpha, oracle_value) {
                violations.push(format!("case[{case_index}]-strong-wolfe-armijo-failed"));
            }
            if !curvature_holds(spec, oracle_derivative) {
                violations.push(format!("case[{case_index}]-strong-wolfe-curvature-failed"));
            }
        } else {
            if case.outcome.alpha != 0.0f64.to_bits() || case.outcome.f_new != spec.f0.to_bits() {
                violations.push(format!(
                    "case[{case_index}]-failure-projection-is-not-alpha-zero-f0"
                ));
            }
            if case.callbacks.iter().any(|callback| {
                armijo_holds(
                    spec,
                    f64::from_bits(callback.alpha),
                    f64::from_bits(callback.value),
                ) && curvature_holds(spec, f64::from_bits(callback.derivative))
            }) {
                violations.push(format!(
                    "case[{case_index}]-failure-skipped-a-retained-strong-wolfe-callback"
                ));
            }
        }
        push_path_specific_violations(case_index, spec, case, &mut violations);
    }
    violations
}

fn validate_semantics(run: &StudyRun) -> Result<(), AdmissionError> {
    let violations = semantic_violations(&run.record);
    if violations.is_empty() {
        Ok(())
    } else {
        Err(AdmissionError::SemanticInconsistency(violations))
    }
}

fn admit_reference(run: &StudyRun, reference: &StudyRun) -> Result<(), AdmissionError> {
    validate_payload(run)?;
    if run.result == reference.result {
        Ok(())
    } else {
        Err(AdmissionError::ReferenceIdentityMismatch {
            expected: reference.result.root(),
            found: run.result.root(),
        })
    }
}

fn reseal(run: &mut StudyRun) {
    run.result = result_identity(&run.fixture, &run.record);
}

fn exact_one_bit_outcome_delta(
    reference: &StudyRun,
    mutant: &StudyRun,
    mutation: Mutation,
) -> bool {
    let Some(mask) = 1u64.checked_shl(mutation.mantissa_bit) else {
        return false;
    };
    if reference.replay_key != mutant.replay_key
        || reference.fixture != mutant.fixture
        || mutation.before ^ mutation.after != mask
        || reference.record.cases.len() != mutant.record.cases.len()
    {
        return false;
    }
    let Some(reference_case) = reference.record.cases.get(mutation.case_index) else {
        return false;
    };
    let Some(mutant_case) = mutant.record.cases.get(mutation.case_index) else {
        return false;
    };
    if reference_case.callbacks != mutant_case.callbacks
        || reference_case.outcome.alpha != mutation.before
        || mutant_case.outcome.alpha != mutation.after
    {
        return false;
    }
    let mut expected = reference.record.clone();
    expected.cases[mutation.case_index].outcome.alpha = mutation.after;
    expected == mutant.record
}

fn seeded_corruption(reference: &StudyRun) -> SeededCorruption {
    let mut selector = StreamKey {
        seed: MUTATION_SEED,
        kernel: MUTATION_KERNEL,
        tile: MUTATION_TILE,
    }
    .stream();
    let successful_cases = [0usize, 1, 2];
    let successful_case_slot =
        usize::try_from(selector.next_below(usize_u64(successful_cases.len())))
            .expect("selected successful-case slot fits usize");
    let case_index = successful_cases[successful_case_slot];
    let mantissa_bit_slot =
        usize::try_from(selector.next_below(usize_u64(MUTATION_MANTISSA_BITS.len())))
            .expect("selected mantissa-bit slot fits usize");
    let mantissa_bit = MUTATION_MANTISSA_BITS[mantissa_bit_slot];
    let selector_draws = selector.index();

    let mut run = reference.clone();
    let before = run.record.cases[case_index].outcome.alpha;
    let after = before ^ (1u64 << mantissa_bit);
    run.record.cases[case_index].outcome.alpha = after;
    assert!(run.record.cases[case_index].outcome.success);
    assert!(f64::from_bits(after).is_finite());
    assert!(f64::from_bits(after) > 0.0);

    let stale_error = validate_payload(&run).expect_err("unsealed mutation must refuse");
    reseal(&mut run);
    let reference_error = admit_reference(&run, reference)
        .expect_err("resealed mutation must not match retained reference");
    let semantic_error = validate_semantics(&run)
        .expect_err("resealed returned-alpha mutation must remain semantically invalid");
    SeededCorruption {
        run,
        mutation: Mutation {
            seed: MUTATION_SEED,
            kernel: MUTATION_KERNEL,
            tile: MUTATION_TILE,
            selector_draws,
            successful_case_slot,
            case_index,
            mantissa_bit_slot,
            mantissa_bit,
            before,
            after,
        },
        stale_error,
        reference_error,
        semantic_error,
    }
}

fn green_receipt(run: &StudyRun) -> Event {
    let callback_count: usize = run
        .record
        .cases
        .iter()
        .map(|case| case.callbacks.len())
        .sum();
    let mut emitter = Emitter::new(SUITE, CASE);
    emitter.emit(
        Severity::Info,
        EventKind::Custom {
            name: "strong-wolfe-full-trajectory-replay-receipt".to_string(),
            json: format!(
                concat!(
                    "{{\"fixture_identity\":\"{}\",\"result_identity\":\"{}\",",
                    "\"replay_seed\":\"0x{:016x}\",\"cases\":{},\"callbacks\":{},",
                    "\"case_evals\":[1,3,5,20,41],",
                    "\"versions\":{{\"fs_ascent\":\"{}\",\"fs_obs\":\"{}\",\"fs_rand\":\"{}\"}},",
                    "\"no_claims\":[\"arbitrary-function-convergence-or-existence\",",
                    "\"optimizer-quality\",\"Riemannian-semantics\",",
                    "\"public-hard-budget-or-Cx-behavior\",\"non-finite-containment\",",
                    "\"cross-process-or-ISA-authority\",\"persistence-or-authentication\",",
                    "\"performance\"]}}"
                ),
                run.fixture.hex(),
                run.result.hex(),
                run.replay_key.seed,
                run.record.cases.len(),
                callback_count,
                fs_ascent::VERSION,
                fs_obs::VERSION,
                fs_rand::VERSION,
            ),
        },
        None,
    )
}

fn green_verdict(run: &StudyRun) -> Event {
    let mut emitter = Emitter::new(SUITE, format!("{CASE}/verdict"));
    emitter.emit(
        Severity::Info,
        EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: CASE.to_string(),
            pass: true,
            detail: format!(
                "fixture={}; result={}; replay_seed=0x{:016x}; five exact callback tables and every WolfeOutcome field admitted",
                run.fixture.hex(),
                run.result.hex(),
                run.replay_key.seed,
            ),
            seed: run.replay_key.seed,
        },
        None,
    )
}

fn corruption_event(reference: &StudyRun, corruption: &SeededCorruption) -> Event {
    let detail = format!(
        "reference={}; mutant={}; seed=0x{:016x}; kernel=0x{:04x}; tile={}; selector_draws={}; successful_case_slot={}; case_index={}; case={}; mantissa_bit_slot={}; mantissa_bit={}; before=0x{:016x}; after=0x{:016x}; callbacks_unchanged=true; stale={:?}; reference_gate={:?}; semantic_gate={:?}",
        reference.result.hex(),
        corruption.run.result.hex(),
        corruption.mutation.seed,
        corruption.mutation.kernel,
        corruption.mutation.tile,
        corruption.mutation.selector_draws,
        corruption.mutation.successful_case_slot,
        corruption.mutation.case_index,
        reference.record.cases[corruption.mutation.case_index].id,
        corruption.mutation.mantissa_bit_slot,
        corruption.mutation.mantissa_bit,
        corruption.mutation.before,
        corruption.mutation.after,
        corruption.stale_error,
        corruption.reference_error,
        corruption.semantic_error,
    );
    let mut emitter = Emitter::new(SUITE, RED_CASE);
    emitter.emit(
        Severity::Error,
        EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: RED_CASE.to_string(),
            pass: false,
            detail,
            seed: MUTATION_SEED,
        },
        None,
    )
}

fn assert_mergeable(event: &Event) {
    let EventKind::ConformanceCase {
        case, pass, detail, ..
    } = &event.kind
    else {
        panic!("merge gate accepts only ConformanceCase evidence");
    };
    assert!(*pass, "merge gate refused {case}: {detail}");
}

fn assert_event_replays(first: &Event, second: &Event, label: &str) {
    assert_eq!(
        first.content_identity().canonical_bytes(),
        second.content_identity().canonical_bytes(),
        "{label} content identity must replay byte-for-byte"
    );
    for event in [first, second] {
        fs_obs::lint_failure_record(event).expect("evidence retains replay inputs");
        fs_obs::validate_line(&event.to_jsonl()).expect("evidence is wire-valid");
        let receipt = event.content_identity_receipt();
        event
            .admit_content_identity(&receipt)
            .expect("fresh evidence content identity admits exactly");
    }
}

#[test]
#[allow(clippy::too_many_lines)] // One causal test spans replay and every refusal gate.
fn strong_wolfe_full_trajectories_replay_and_seeded_alpha_corruption_is_refused() {
    let key = replay_key();
    let original = run_study(key);
    let replay = run_study(key);
    for run in [&original, &replay] {
        assert_eq!(validate_payload(run), Ok(()));
        assert_eq!(validate_semantics(run), Ok(()));
        assert_eq!(admit_reference(run, &original), Ok(()));
    }
    assert_eq!(original, replay, "same-key study must replay exactly");
    assert_eq!(original.record, replay.record);
    assert_eq!(
        original.fixture.canonical_bytes(),
        replay.fixture.canonical_bytes(),
        "same-seed fixture identity must replay byte-for-byte"
    );
    assert_eq!(
        original.result.canonical_bytes(),
        replay.result.canonical_bytes(),
        "same-seed complete result identity must replay byte-for-byte"
    );

    let first_receipt = green_receipt(&original);
    let second_receipt = green_receipt(&replay);
    assert_event_replays(&first_receipt, &second_receipt, "green receipt");
    println!("{}", first_receipt.to_jsonl());

    let first_green = green_verdict(&original);
    let second_green = green_verdict(&replay);
    assert_event_replays(&first_green, &second_green, "green verdict");
    assert_mergeable(&first_green);
    assert_mergeable(&second_green);
    println!("{}", first_green.to_jsonl());

    let first = seeded_corruption(&original);
    let second = seeded_corruption(&replay);
    assert_eq!(first, second, "seeded corruption must replay exactly");
    assert!(
        exact_one_bit_outcome_delta(&original, &first.run, first.mutation),
        "mutation must change exactly one successful returned-alpha bit and no callback"
    );
    assert_eq!(
        original.record.cases[first.mutation.case_index].callbacks,
        first.run.record.cases[first.mutation.case_index].callbacks,
        "the complete callback trace must remain unchanged"
    );
    assert_eq!(
        validate_payload(&first.run),
        Ok(()),
        "resealed mutation must be internally self-consistent"
    );
    assert!(matches!(
        &first.stale_error,
        AdmissionError::PayloadIdentityMismatch { declared, computed }
            if *declared == original.result.root()
                && *computed == first.run.result.root()
    ));
    assert!(matches!(
        &first.reference_error,
        AdmissionError::ReferenceIdentityMismatch { expected, found }
            if *expected == original.result.root()
                && *found == first.run.result.root()
    ));
    let AdmissionError::SemanticInconsistency(violations) = &first.semantic_error else {
        panic!("seeded alpha corruption must carry semantic violations");
    };
    assert!(violations.iter().any(|violation| {
        violation
            == &format!(
                "case[{}]-success-alpha-not-retained-callback",
                first.mutation.case_index
            )
    }));
    assert!(violations.iter().any(|violation| {
        violation
            == &format!(
                "case[{}]-returned-value-does-not-match-curve-at-alpha",
                first.mutation.case_index
            )
    }));

    let first_red = corruption_event(&original, &first);
    let second_red = corruption_event(&replay, &second);
    assert_event_replays(&first_red, &second_red, "stable-red evidence");
    println!("{}", first_red.to_jsonl());

    let panic = catch_unwind(|| assert_mergeable(&first_red))
        .expect_err("merge gate must refuse seeded Strong-Wolfe corruption");
    let message = panic
        .downcast_ref::<String>()
        .map(String::as_str)
        .or_else(|| panic.downcast_ref::<&str>().copied())
        .expect("merge-gate panic carries text");
    assert!(message.contains(RED_CASE));
    assert!(message.contains(&format!("0x{MUTATION_SEED:016x}")));
    assert!(message.contains("callbacks_unchanged=true"));
    assert!(message.contains("PayloadIdentityMismatch"));
    assert!(message.contains("ReferenceIdentityMismatch"));
    assert!(message.contains("SemanticInconsistency"));
}
