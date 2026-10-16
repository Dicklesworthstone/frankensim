//! FrankenScipy special-function oracle evidence for BEDROCK `erf`/`erfc`.
//!
//! Exact IEEE fixtures cover signed zero and infinities, while disclosed
//! finite fixtures compare `fs-math` with the test-only `fsci-special` oracle.
//! The `erf` fixture excludes the known cancellation-sensitive band
//! `1.5 < |x| < 3.5`; its budget-grade
//! evidence remains the independent Taylor/continued-fraction cross-check in
//! `extensions_battery.rs`. The positive `erfc` points through 26 are dyadic
//! or integral, so `x*x` is exact in both implementations. Points at and above
//! `sqrt(CEPHES_MAXLOG) ~= 26.615` are excluded because the oracle zeros there
//! before `fs-math`'s subnormal-preserving 27.5 cutoff.
//!
//! This is portable G0 oracle evidence. It makes no Python-SciPy execution,
//! correctly-rounded, performance-parity, arbitrary-domain, fresh dual-ISA,
//! or full-G5 claim. NaN propagation is structural only: payload equality is
//! intentionally not claimed.

use core::fmt::Write as _;
use std::panic::catch_unwind;

use fs_casebook::{CASEBOOK_RECORD_VERSION, CaseOutcome, Suite, ToleranceSpec, fnv1a64};
use fs_math::{det, ulp_distance};
use fsci_special::{erf_scalar, erfc_scalar};

const SUITE: &str = "bedrock/fs-math-frankenscipy-special-oracle-v1";
const ORACLE_VERSION: &str = "fsci-special/0.1.0";
const ORACLE_API: &str = "fs_math::det::{erf,erfc}+fsci_special::{erf_scalar,erfc_scalar}:v1";
const ERF_ULP_BOUND: u64 = 12;
const ERFC_ULP_BOUND: u64 = 16;
const CANONICAL_NAN_BITS: u64 = 0x7ff8_0000_0000_0000;
const ERF_POINTS: [f64; 14] = [
    -6.0, -4.0, -3.5, -1.5, -1.0, -0.5, -0.125, 0.125, 0.5, 1.0, 1.5, 3.5, 4.0, 6.0,
];
const ERFC_POINTS: [f64; 14] = [
    0.125, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 6.0, 8.0, 12.0, 20.0, 26.0,
];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct ExpectedBits {
    erf: u64,
    erfc: u64,
}

const SPECIAL_POINTS_BITS: [u64; 4] = [
    0x0000_0000_0000_0000,
    0x8000_0000_0000_0000,
    0x7ff0_0000_0000_0000,
    0xfff0_0000_0000_0000,
];

const SPECIAL_EXPECTED: [ExpectedBits; 4] = [
    ExpectedBits {
        erf: 0x0000_0000_0000_0000,
        erfc: 0x3ff0_0000_0000_0000,
    },
    ExpectedBits {
        erf: 0x8000_0000_0000_0000,
        erfc: 0x3ff0_0000_0000_0000,
    },
    ExpectedBits {
        erf: 0x3ff0_0000_0000_0000,
        erfc: 0x0000_0000_0000_0000,
    },
    ExpectedBits {
        erf: 0xbff0_0000_0000_0000,
        erfc: 0x4000_0000_0000_0000,
    },
];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct SampleBits {
    frankensim_erf: u64,
    frankenscipy_erf: u64,
    frankensim_erfc: u64,
    frankenscipy_erfc: u64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct OracleMeasurement {
    frankensim: Vec<u64>,
    frankenscipy: Vec<u64>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Function {
    Erf,
    Erfc,
}

impl Function {
    fn label(self) -> &'static str {
        match self {
            Self::Erf => "erf",
            Self::Erfc => "erfc",
        }
    }

    fn points(self) -> &'static [f64] {
        match self {
            Self::Erf => &ERF_POINTS,
            Self::Erfc => &ERFC_POINTS,
        }
    }

    fn domain(self) -> &'static str {
        match self {
            Self::Erf => "fixed-dyadic-or-integral;abs(x)<=1.5-or-3.5<=abs(x)<=6",
            Self::Erfc => "fixed-positive-dyadic-or-integral;0.125<=x<=26;exact-x-square",
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct Corruption {
    seed: u64,
    point: usize,
    function: &'static str,
    bit: u32,
}

fn push_u64(bytes: &mut Vec<u8>, value: u64) {
    bytes.extend_from_slice(&value.to_le_bytes());
}

fn push_len(bytes: &mut Vec<u8>, value: usize) {
    push_u64(
        bytes,
        u64::try_from(value).expect("special-oracle fixture lengths fit u64"),
    );
}

fn push_text(bytes: &mut Vec<u8>, value: &str) {
    push_len(bytes, value.len());
    bytes.extend_from_slice(value.as_bytes());
}

fn push_u64s(bytes: &mut Vec<u8>, values: &[u64]) {
    push_len(bytes, values.len());
    for &value in values {
        push_u64(bytes, value);
    }
}

fn push_f64s(bytes: &mut Vec<u8>, values: &[f64]) {
    push_len(bytes, values.len());
    for value in values {
        push_u64(bytes, value.to_bits());
    }
}

fn push_expected(bytes: &mut Vec<u8>, expected: &[ExpectedBits]) {
    push_len(bytes, expected.len());
    for value in expected {
        push_u64(bytes, value.erf);
        push_u64(bytes, value.erfc);
    }
}

fn push_nested(bytes: &mut Vec<u8>, label: &str, frame: &[u8]) {
    push_text(bytes, label);
    push_len(bytes, frame.len());
    bytes.extend_from_slice(frame);
}

fn hex_bytes(bytes: &[u8]) -> String {
    let mut encoded = String::with_capacity(bytes.len() * 2);
    for byte in bytes {
        write!(&mut encoded, "{byte:02x}").expect("writing to String cannot fail");
    }
    encoded
}

fn digest_bits(bits: &[u64]) -> u64 {
    let mut bytes = Vec::with_capacity(bits.len() * 8);
    for &value in bits {
        bytes.extend_from_slice(&value.to_le_bytes());
    }
    fnv1a64(&bytes)
}

fn common_frame_prefix(domain: &[u8]) -> Vec<u8> {
    let mut bytes = domain.to_vec();
    push_text(&mut bytes, "encoding");
    push_text(
        &mut bytes,
        "length-prefixed-little-endian-u64-and-f64-bits:v1",
    );
    push_text(&mut bytes, "casebook-record-version");
    push_u64(&mut bytes, u64::from(CASEBOOK_RECORD_VERSION));
    push_text(&mut bytes, "fs-math-version");
    push_text(&mut bytes, env!("CARGO_PKG_VERSION"));
    push_text(&mut bytes, "oracle-version");
    push_text(&mut bytes, ORACLE_VERSION);
    push_text(&mut bytes, "bound-apis");
    push_text(&mut bytes, ORACLE_API);
    bytes
}

fn special_inputs(expected: &[ExpectedBits]) -> Vec<u8> {
    let mut bytes = common_frame_prefix(b"bedrock:fs-math:frankenscipy-special-values:v1");
    push_text(&mut bytes, "point-bits:+0,-0,+inf,-inf");
    push_u64s(&mut bytes, &SPECIAL_POINTS_BITS);
    push_text(&mut bytes, "expected-erf-then-erfc-bits");
    push_expected(&mut bytes, expected);
    push_text(&mut bytes, "outputs-per-backend");
    push_len(&mut bytes, SPECIAL_POINTS_BITS.len() * 2);
    push_text(&mut bytes, "exact-complement-identity-per-point:v1");
    bytes
}

fn nan_inputs() -> Vec<u8> {
    let mut bytes = common_frame_prefix(b"bedrock:fs-math:frankenscipy-nan-policy:v1");
    push_text(&mut bytes, "canonical-qnan-input-bits");
    push_u64(&mut bytes, CANONICAL_NAN_BITS);
    push_text(&mut bytes, "policy");
    push_text(
        &mut bytes,
        "all-four-results-is-nan;payload-and-sign-not-compared",
    );
    push_text(&mut bytes, "input-point-count");
    push_len(&mut bytes, 1);
    push_text(&mut bytes, "backend-function-output-count");
    push_len(&mut bytes, 4);
    bytes
}

fn oracle_inputs(function: Function, points: &[f64], ulp_bound: u64) -> Vec<u8> {
    let mut bytes = common_frame_prefix(b"bedrock:fs-math:frankenscipy-finite-oracle:v1");
    push_text(&mut bytes, "function");
    push_text(&mut bytes, function.label());
    push_text(&mut bytes, "fixture");
    push_text(&mut bytes, "fixed-dyadic-or-integral-v1");
    push_text(&mut bytes, "admitted-domain");
    push_text(&mut bytes, function.domain());
    push_text(&mut bytes, "point-bits");
    push_f64s(&mut bytes, points);
    push_text(&mut bytes, "ulp-bound");
    push_u64(&mut bytes, ulp_bound);
    push_text(
        &mut bytes,
        "no-claims:erf-band-1.5<abs(x)<3.5;erfc-x>=sqrt(CEPHES_MAXLOG)",
    );
    bytes
}

fn corruption_inputs(
    seed: u64,
    point: usize,
    bit: u32,
    canonical: &[ExpectedBits],
    corrupted: &[ExpectedBits],
) -> Vec<u8> {
    let special = special_inputs(canonical);
    let mut bytes = common_frame_prefix(b"bedrock:fs-math:seeded-special-reference-corruption:v1");
    push_u64(&mut bytes, seed);
    push_len(&mut bytes, point);
    push_text(&mut bytes, "function");
    push_text(&mut bytes, "erf");
    push_u64(&mut bytes, u64::from(bit));
    push_nested(&mut bytes, "nested-canonical-special-frame", &special);
    push_text(&mut bytes, "canonical-reference-bits");
    push_expected(&mut bytes, canonical);
    push_text(&mut bytes, "corrupted-reference-bits");
    push_expected(&mut bytes, corrupted);
    bytes
}

fn evaluate(x: f64) -> SampleBits {
    SampleBits {
        frankensim_erf: det::erf(x).to_bits(),
        frankenscipy_erf: erf_scalar(x).to_bits(),
        frankensim_erfc: det::erfc(x).to_bits(),
        frankenscipy_erfc: erfc_scalar(x).to_bits(),
    }
}

fn measure_samples(points: &[f64]) -> Vec<SampleBits> {
    points.iter().copied().map(evaluate).collect()
}

fn backend_special_bits(samples: &[SampleBits], frankensim: bool) -> Vec<u64> {
    let mut bits = Vec::with_capacity(samples.len() * 2);
    for sample in samples {
        if frankensim {
            bits.extend_from_slice(&[sample.frankensim_erf, sample.frankensim_erfc]);
        } else {
            bits.extend_from_slice(&[sample.frankenscipy_erf, sample.frankenscipy_erfc]);
        }
    }
    bits
}

fn corruption_context(corruption: Option<Corruption>) -> String {
    corruption.map_or_else(
        || "mode=canonical".to_owned(),
        |value| {
            format!(
                "seed=0x{:016x}; point={}; function={}; bit={}",
                value.seed, value.point, value.function, value.bit,
            )
        },
    )
}

fn special_outcome(
    reference: [ExpectedBits; SPECIAL_POINTS_BITS.len()],
    corruption: Option<Corruption>,
    input_frame: &[u8],
) -> CaseOutcome {
    let points = SPECIAL_POINTS_BITS.map(f64::from_bits);
    let first = measure_samples(&points);
    let replay = measure_samples(&points);
    let inputs_hex = hex_bytes(input_frame);
    let context = corruption_context(corruption);
    if first != replay {
        return CaseOutcome::fail(format!(
            "{context}; stage=same-run-replay; first={first:016x?}; replay={replay:016x?}; inputs_hex={inputs_hex}"
        ))
        .with_evidence("crates/fs-math/CONTRACT.md#determinism-class");
    }
    if first.len() != SPECIAL_POINTS_BITS.len() {
        return CaseOutcome::fail(format!(
            "{context}; stage=special-output-length; expected={}; actual={}; outputs={first:016x?}; inputs_hex={inputs_hex}",
            SPECIAL_POINTS_BITS.len(),
            first.len(),
        ))
        .with_evidence("crates/fs-math/CONTRACT.md#invariants");
    }
    for (index, (sample, expected)) in first.iter().zip(reference).enumerate() {
        if sample.frankensim_erf != expected.erf
            || sample.frankenscipy_erf != expected.erf
            || sample.frankensim_erfc != expected.erfc
            || sample.frankenscipy_erfc != expected.erfc
        {
            return CaseOutcome::fail(format!(
                "{context}; stage=ieee-special-known-answer; point={index}; input_bits=0x{:016x}; frankensim_erf=0x{:016x}; frankenscipy_erf=0x{:016x}; expected_erf=0x{:016x}; frankensim_erfc=0x{:016x}; frankenscipy_erfc=0x{:016x}; expected_erfc=0x{:016x}; inputs_hex={inputs_hex}",
                SPECIAL_POINTS_BITS[index],
                sample.frankensim_erf,
                sample.frankenscipy_erf,
                expected.erf,
                sample.frankensim_erfc,
                sample.frankenscipy_erfc,
                expected.erfc,
            ))
            .with_evidence("crates/fs-math/CONTRACT.md#invariants")
            .with_evidence("constellation.lock:frankenscipy-0.1.0");
        }
        let frankensim_sum =
            f64::from_bits(sample.frankensim_erf) + f64::from_bits(sample.frankensim_erfc);
        let frankenscipy_sum =
            f64::from_bits(sample.frankenscipy_erf) + f64::from_bits(sample.frankenscipy_erfc);
        if frankensim_sum.to_bits() != 1.0_f64.to_bits()
            || frankenscipy_sum.to_bits() != 1.0_f64.to_bits()
        {
            return CaseOutcome::fail(format!(
                "{context}; stage=exact-special-complement; point={index}; input_bits=0x{:016x}; frankensim_sum_bits=0x{:016x}; frankenscipy_sum_bits=0x{:016x}; expected_sum_bits=0x{:016x}; inputs_hex={inputs_hex}",
                SPECIAL_POINTS_BITS[index],
                frankensim_sum.to_bits(),
                frankenscipy_sum.to_bits(),
                1.0_f64.to_bits(),
            ))
            .with_evidence("crates/fs-math/CONTRACT.md#invariants");
        }
    }
    let frankensim_digest = digest_bits(&backend_special_bits(&first, true));
    let frankenscipy_digest = digest_bits(&backend_special_bits(&first, false));
    CaseOutcome::pass(format!(
        "points=4; outputs_per_backend=8; exact_complements=8; frankensim_output_digest={frankensim_digest:016x}; frankenscipy_output_digest={frankenscipy_digest:016x}; same_run=identical"
    ))
    .with_evidence("crates/fs-math/CONTRACT.md#invariants")
    .with_evidence("constellation.lock:frankenscipy-0.1.0")
}

fn nan_outcome(input_frame: &[u8]) -> CaseOutcome {
    let point = f64::from_bits(CANONICAL_NAN_BITS);
    let first = measure_samples(&[point]);
    let replay = measure_samples(&[point]);
    let inputs_hex = hex_bytes(input_frame);
    if first != replay {
        return CaseOutcome::fail(format!(
            "stage=same-run-nan-replay; first={first:016x?}; replay={replay:016x?}; inputs_hex={inputs_hex}"
        ))
        .with_evidence("crates/fs-math/CONTRACT.md#determinism-class");
    }
    let [sample] = first.as_slice() else {
        return CaseOutcome::fail(format!(
            "stage=nan-output-length; expected=1; actual={}; outputs={first:016x?}; inputs_hex={inputs_hex}",
            first.len(),
        ))
        .with_evidence("crates/fs-math/CONTRACT.md#error-model");
    };
    let outputs = [
        sample.frankensim_erf,
        sample.frankenscipy_erf,
        sample.frankensim_erfc,
        sample.frankenscipy_erfc,
    ];
    if outputs.iter().any(|bits| !f64::from_bits(*bits).is_nan()) {
        return CaseOutcome::fail(format!(
            "stage=nan-propagation; input_bits=0x{CANONICAL_NAN_BITS:016x}; frankensim_erf=0x{:016x}; frankenscipy_erf=0x{:016x}; frankensim_erfc=0x{:016x}; frankenscipy_erfc=0x{:016x}; inputs_hex={inputs_hex}",
            outputs[0], outputs[1], outputs[2], outputs[3],
        ))
        .with_evidence("crates/fs-math/CONTRACT.md#error-model");
    }
    CaseOutcome::pass(format!(
        "input=canonical_qnan; results=4_nan; payload_claim=none; output_digest={:016x}; same_run=identical",
        digest_bits(&outputs),
    ))
    .with_evidence("crates/fs-math/CONTRACT.md#error-model")
}

fn measure_oracle(function: Function, points: &[f64]) -> OracleMeasurement {
    let mut frankensim = Vec::with_capacity(points.len());
    let mut frankenscipy = Vec::with_capacity(points.len());
    for &point in points {
        match function {
            Function::Erf => {
                frankensim.push(det::erf(point).to_bits());
                frankenscipy.push(erf_scalar(point).to_bits());
            }
            Function::Erfc => {
                frankensim.push(det::erfc(point).to_bits());
                frankenscipy.push(erfc_scalar(point).to_bits());
            }
        }
    }
    OracleMeasurement {
        frankensim,
        frankenscipy,
    }
}

fn oracle_outcome(
    function: Function,
    points: &[f64],
    ulp_bound: u64,
    input_frame: &[u8],
) -> CaseOutcome {
    let first = measure_oracle(function, points);
    let replay = measure_oracle(function, points);
    let inputs_hex = hex_bytes(input_frame);
    if first != replay {
        return CaseOutcome::fail(format!(
            "stage=same-run-{}-replay; first_frankensim={:016x?}; first_frankenscipy={:016x?}; replay_frankensim={:016x?}; replay_frankenscipy={:016x?}; inputs_hex={inputs_hex}",
            function.label(),
            first.frankensim,
            first.frankenscipy,
            replay.frankensim,
            replay.frankenscipy,
        ))
        .with_evidence("crates/fs-math/CONTRACT.md#determinism-class");
    }
    if first.frankensim.len() != points.len() || first.frankenscipy.len() != points.len() {
        return CaseOutcome::fail(format!(
            "stage={}-output-length; expected={}; frankensim_len={}; frankenscipy_len={}; frankensim={:016x?}; frankenscipy={:016x?}; inputs_hex={inputs_hex}",
            function.label(),
            points.len(),
            first.frankensim.len(),
            first.frankenscipy.len(),
            first.frankensim,
            first.frankenscipy,
        ))
        .with_evidence("crates/fs-math/CONTRACT.md#invariants");
    }

    let mut max_ulps = 0_u64;
    let mut max_abs = 0.0_f64;
    for (index, ((&frankensim_bits, &frankenscipy_bits), &point)) in first
        .frankensim
        .iter()
        .zip(&first.frankenscipy)
        .zip(points)
        .enumerate()
    {
        let frankensim = f64::from_bits(frankensim_bits);
        let frankenscipy = f64::from_bits(frankenscipy_bits);
        let ulps = ulp_distance(frankensim, frankenscipy);
        let abs = (frankensim - frankenscipy).abs();
        max_ulps = max_ulps.max(ulps);
        max_abs = max_abs.max(abs);
        if !frankensim.is_finite()
            || !frankenscipy.is_finite()
            || !abs.is_finite()
            || ulps > ulp_bound
        {
            return CaseOutcome::fail(format!(
                "stage={}-oracle-agreement; fixture=fixed; point={index}; input_bits=0x{:016x}; frankensim_bits=0x{frankensim_bits:016x}; frankenscipy_bits=0x{frankenscipy_bits:016x}; ulps={ulps}; ulp_bound={ulp_bound}; abs={abs:.17e}; domain={}; inputs_hex={inputs_hex}",
                function.label(),
                point.to_bits(),
                function.domain(),
            ))
            .with_evidence("crates/fs-math/CONTRACT.md#public-types-and-semantics")
            .with_evidence("constellation.lock:frankenscipy-0.1.0");
        }
    }

    CaseOutcome::pass(format!(
        "function={}; fixture=fixed; points={}; domain={}; ulp_bound={ulp_bound}; max_ulps={max_ulps}; max_abs={max_abs:.17e}; frankensim_output_digest={:016x}; frankenscipy_output_digest={:016x}; same_run=identical",
        function.label(),
        points.len(),
        function.domain(),
        digest_bits(&first.frankensim),
        digest_bits(&first.frankenscipy),
    ))
    .with_evidence("crates/fs-math/CONTRACT.md#public-types-and-semantics")
    .with_evidence("constellation.lock:frankenscipy-0.1.0")
}

#[test]
#[allow(clippy::too_many_lines)] // Keep every canonical frame beside its registered case.
fn frankenscipy_special_oracle_casebook_emits_replay_complete_green_records() {
    assert_eq!(CASEBOOK_RECORD_VERSION, 1);
    let special_frame = special_inputs(&SPECIAL_EXPECTED);
    let nan_frame = nan_inputs();
    let erf_points = Function::Erf.points();
    let erfc_points = Function::Erfc.points();
    let erf_frame = oracle_inputs(Function::Erf, erf_points, ERF_ULP_BOUND);
    let erfc_frame = oracle_inputs(Function::Erfc, erfc_points, ERFC_ULP_BOUND);
    let special_digest = fnv1a64(&special_frame);
    let nan_digest = fnv1a64(&nan_frame);
    let erf_digest = fnv1a64(&erf_frame);
    let erfc_digest = fnv1a64(&erfc_frame);

    // Canonical frame pins are reconstructed independently from the disclosed
    // LE encoding and generator; drift here is an input/schema change.
    assert_eq!(special_digest, 0x9e82_a07f_ef5e_77cd);
    assert_eq!(nan_digest, 0xea07_b7ba_4630_b7ba);
    assert_eq!(erf_digest, 0x6357_985e_ff78_4c7e);
    assert_eq!(erfc_digest, 0x3735_a50a_7ae0_0f90);

    let report = Suite::new(SUITE)
        .case(
            "ieee-special-values-known-answer",
            special_digest,
            ToleranceSpec::Exact,
            move || special_outcome(SPECIAL_EXPECTED, None, &special_frame),
        )
        .case(
            "nan-propagation-policy",
            nan_digest,
            ToleranceSpec::Structural,
            move || nan_outcome(&nan_frame),
        )
        .case(
            "fixed-erf-strong-domain-oracle-agreement",
            erf_digest,
            ToleranceSpec::Ulps(u32::try_from(ERF_ULP_BOUND).expect("erf ULP bound fits u32")),
            move || oracle_outcome(Function::Erf, erf_points, ERF_ULP_BOUND, &erf_frame),
        )
        .case(
            "fixed-erfc-strong-domain-oracle-agreement",
            erfc_digest,
            ToleranceSpec::Ulps(u32::try_from(ERFC_ULP_BOUND).expect("erfc ULP bound fits u32")),
            move || oracle_outcome(Function::Erfc, erfc_points, ERFC_ULP_BOUND, &erfc_frame),
        )
        .run();

    report.assert_green();
    assert_eq!(
        report
            .records
            .iter()
            .map(|record| record.case.as_str())
            .collect::<Vec<_>>(),
        [
            "ieee-special-values-known-answer",
            "nan-propagation-policy",
            "fixed-erf-strong-domain-oracle-agreement",
            "fixed-erfc-strong-domain-oracle-agreement",
        ]
    );
    assert!(
        report.records[0]
            .json_line()
            .contains("frankensim_output_digest=")
    );
    assert!(
        report.records[0]
            .json_line()
            .contains("frankenscipy_output_digest=")
    );
    assert!(report.records[0].json_line().contains("0701f60ae5d0a5e5"));
}

#[test]
fn disclosed_seeded_special_reference_corruption_turns_suite_red() {
    const CORRUPTION_SEED: u64 = 0x5AEC_0000;
    let point = (CORRUPTION_SEED & 0x3) as usize;
    let bit = ((CORRUPTION_SEED >> 2) & 0x3f) as u32;
    assert_eq!((point, bit), (0, 0));
    let corruption = Corruption {
        seed: CORRUPTION_SEED,
        point,
        function: "erf",
        bit,
    };
    let mut corrupted = SPECIAL_EXPECTED;
    corrupted[point].erf ^= 1_u64 << bit;
    assert_eq!(corrupted[point].erf, SPECIAL_EXPECTED[point].erf ^ 1);

    let inputs = corruption_inputs(CORRUPTION_SEED, point, bit, &SPECIAL_EXPECTED, &corrupted);
    let inputs_digest = fnv1a64(&inputs);
    assert_eq!(inputs_digest, 0xb85f_44ff_0d7d_3109);

    let make_report = || {
        let input_frame = inputs.clone();
        Suite::new(SUITE)
            .case(
                "seeded-special-expected-reference-corruption",
                inputs_digest,
                ToleranceSpec::Exact,
                move || special_outcome(corrupted, Some(corruption), &input_frame),
            )
            .run()
    };
    let first = make_report();
    let replay = make_report();
    let first_failures = first.failures();
    let replay_failures = replay.failures();
    let [first_failure] = first_failures.as_slice() else {
        panic!("the disclosed special-reference corruption must produce exactly one failure");
    };
    let [replay_failure] = replay_failures.as_slice() else {
        panic!("the replayed special-reference corruption must produce exactly one failure");
    };
    assert_eq!(first_failure.json_line(), replay_failure.json_line());
    assert_eq!(
        first_failure.case,
        "seeded-special-expected-reference-corruption"
    );
    assert_eq!(first_failure.inputs_digest, "b85f44ff0d7d3109");
    assert!(
        first_failure
            .details
            .contains("stage=ieee-special-known-answer")
    );
    assert!(
        first_failure
            .details
            .contains(&format!("seed=0x{CORRUPTION_SEED:016x}"))
    );
    assert!(
        first_failure
            .details
            .contains("point=0; function=erf; bit=0")
    );
    assert!(first_failure.details.contains("inputs_hex="));
    assert!(
        first_failure
            .json_line()
            .contains("\"tolerance\":\"exact\",\"pass\":false")
    );

    let panic = catch_unwind(|| first.assert_green())
        .expect_err("the merge gate must reject the disclosed special-reference corruption");
    let message = panic
        .downcast_ref::<String>()
        .map(String::as_str)
        .or_else(|| panic.downcast_ref::<&str>().copied())
        .expect("casebook panic carries text");
    assert!(message.contains("seeded-special-expected-reference-corruption"));
    assert!(message.contains(&format!("seed=0x{CORRUPTION_SEED:016x}")));
}
