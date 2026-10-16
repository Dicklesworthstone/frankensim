//! FrankenScipy FFT oracle evidence for the BEDROCK transform spine.
//!
//! A unit-impulse KAT is exact, while small power-of-two complex and real
//! fixtures compare `fs-fft` with FrankenScipy's deliberately independent
//! O(n^2) `NaiveDft` backend. The bounded records use the same scaled-complex
//! agreement margins as the existing in-crate oracle batteries.
//!
//! This is portable G0 agreement evidence, not a forward-error certificate.
//! It does not run Python SciPy, cover inverse/c2r, general-n, N-D, or DCT
//! transforms, claim performance parity, or establish fresh dual-ISA G5.
//! FrankenScipy's naive oracle uses platform trigonometry, so finite output
//! digests are observed same-run evidence and are intentionally not pinned.

use core::fmt::Write as _;
use std::panic::{AssertUnwindSafe, catch_unwind};

use fs_casebook::{CASEBOOK_RECORD_VERSION, CaseOutcome, Suite, ToleranceSpec, fnv1a64};
use fs_fft::{C64, Fft, RealFft, TRANSFORM_BIT_SEMANTICS_VERSION, VERSION as FS_FFT_VERSION};
use fs_math::VERSION as FS_MATH_VERSION;
use fsci_fft::{BackendKind, FftOptions, Normalization, WorkerPolicy, fft};

const SUITE: &str = "bedrock/fs-fft-frankenscipy-oracle-v1";
const ORACLE_VERSION: &str = "fsci-fft/0.1.0";
const ORACLE_PIN: &str = "9e271fd734465e2b2ff755aa73ea66a7217d619b";
const ORACLE_API: &str = "fsci_fft::{fft,FftOptions,BackendKind::NaiveDft,Normalization::Backward,WorkerPolicy::Exact}:v1";
const OPTIONS: &str = "mode=Strict(default);normalization=Backward;workers=Exact(1);backend=NaiveDft;check_finite=true;overwrite_input=false";
const SIZES: [usize; 4] = [4, 8, 16, 32];
const FIXTURES_PER_SIZE: usize = 2;
const COMPLEX_SEED: u64 = 0xF17F_C04D_0000_0001;
const REAL_SEED: u64 = 0xF17F_2EA1_0000_0001;
const LCG_MULTIPLIER: u64 = 6_364_136_223_846_793_005;
const LCG_INCREMENT: u64 = 1_442_695_040_888_963_407;
const COMPLEX_BOUND: f64 = 1.0e-12;
const REAL_BOUND: f64 = 1.0e-11;
const KAT_N: usize = 8;
const KAT_EXPECTED: [[u64; 2]; KAT_N] = [[0x3ff0_0000_0000_0000, 0x0000_0000_0000_0000]; KAT_N];

type ComplexBits = [u64; 2];

#[derive(Debug, Clone)]
struct ComplexFixture {
    n: usize,
    input: Vec<ComplexBits>,
}

#[derive(Debug, Clone)]
struct RealFixture {
    n: usize,
    input: Vec<u64>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct Measurement {
    frankensim: Result<Vec<ComplexBits>, String>,
    frankenscipy: Result<Vec<ComplexBits>, String>,
}

#[derive(Debug, Clone, Copy)]
struct Corruption {
    seed: u64,
    output: usize,
    component: usize,
    bit: u32,
}

fn push_u64(bytes: &mut Vec<u8>, value: u64) {
    bytes.extend_from_slice(&value.to_le_bytes());
}

fn push_len(bytes: &mut Vec<u8>, value: usize) {
    push_u64(
        bytes,
        u64::try_from(value).expect("FFT oracle fixture lengths fit u64"),
    );
}

fn push_text(bytes: &mut Vec<u8>, value: &str) {
    push_len(bytes, value.len());
    bytes.extend_from_slice(value.as_bytes());
}

fn push_complex_bits(bytes: &mut Vec<u8>, values: &[ComplexBits]) {
    push_len(bytes, values.len());
    for [re, im] in values {
        push_u64(bytes, *re);
        push_u64(bytes, *im);
    }
}

fn push_real_bits(bytes: &mut Vec<u8>, values: &[u64]) {
    push_len(bytes, values.len());
    for &value in values {
        push_u64(bytes, value);
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

fn digest_complex_bits(values: &[ComplexBits]) -> u64 {
    let mut bytes = Vec::with_capacity(values.len() * 16);
    for [re, im] in values {
        bytes.extend_from_slice(&re.to_le_bytes());
        bytes.extend_from_slice(&im.to_le_bytes());
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
    push_text(&mut bytes, "fs-fft-version");
    push_text(&mut bytes, FS_FFT_VERSION);
    push_text(&mut bytes, "fs-fft-transform-bit-semantics-version");
    push_u64(&mut bytes, u64::from(TRANSFORM_BIT_SEMANTICS_VERSION));
    push_text(&mut bytes, "fs-math-version");
    push_text(&mut bytes, FS_MATH_VERSION);
    push_text(&mut bytes, "oracle-version");
    push_text(&mut bytes, ORACLE_VERSION);
    push_text(&mut bytes, "oracle-pin");
    push_text(&mut bytes, ORACLE_PIN);
    push_text(&mut bytes, "oracle-api");
    push_text(&mut bytes, ORACLE_API);
    push_text(&mut bytes, "oracle-options");
    push_text(&mut bytes, OPTIONS);
    bytes
}

fn kat_fixture() -> ComplexFixture {
    let mut input = vec![[0, 0]; KAT_N];
    input[0] = [1.0_f64.to_bits(), 0];
    ComplexFixture { n: KAT_N, input }
}

fn kat_inputs(expected: &[ComplexBits; KAT_N]) -> Vec<u8> {
    let mut bytes = common_frame_prefix(b"bedrock:fs-fft:frankenscipy-unit-impulse-kat:v1");
    push_text(&mut bytes, "analytic-unit-impulse");
    push_len(&mut bytes, KAT_N);
    push_text(&mut bytes, "input-bits");
    push_complex_bits(&mut bytes, &kat_fixture().input);
    push_text(&mut bytes, "expected-output-bits");
    push_complex_bits(&mut bytes, expected);
    push_text(&mut bytes, "expected-property:fft(delta_0)=all-(1,+0):v1");
    bytes
}

fn complex_inputs(fixtures: &[ComplexFixture]) -> Vec<u8> {
    let mut bytes = common_frame_prefix(b"bedrock:fs-fft:frankenscipy-complex-naive:v1");
    push_text(&mut bytes, "lcg64-top53-centered-v1");
    for value in [COMPLEX_SEED, LCG_MULTIPLIER, LCG_INCREMENT] {
        push_u64(&mut bytes, value);
    }
    push_text(&mut bytes, "sizes");
    push_len(&mut bytes, SIZES.len());
    for size in SIZES {
        push_len(&mut bytes, size);
    }
    push_text(&mut bytes, "fixtures-per-size");
    push_len(&mut bytes, FIXTURES_PER_SIZE);
    push_text(&mut bytes, "scaled-complex-bound");
    push_u64(&mut bytes, COMPLEX_BOUND.to_bits());
    push_text(
        &mut bytes,
        "metric=abs(delta)/max(1,sum(abs(complex-input)))",
    );
    push_len(&mut bytes, fixtures.len());
    for fixture in fixtures {
        push_len(&mut bytes, fixture.n);
        push_complex_bits(&mut bytes, &fixture.input);
    }
    bytes
}

fn real_inputs(fixtures: &[RealFixture]) -> Vec<u8> {
    let mut bytes = common_frame_prefix(b"bedrock:fs-fft:frankenscipy-real-naive:v1");
    push_text(&mut bytes, "lcg64-top53-centered-v1");
    for value in [REAL_SEED, LCG_MULTIPLIER, LCG_INCREMENT] {
        push_u64(&mut bytes, value);
    }
    push_text(&mut bytes, "sizes");
    push_len(&mut bytes, SIZES.len());
    for size in SIZES {
        push_len(&mut bytes, size);
    }
    push_text(&mut bytes, "fixtures-per-size");
    push_len(&mut bytes, FIXTURES_PER_SIZE);
    push_text(&mut bytes, "scaled-complex-bound");
    push_u64(&mut bytes, REAL_BOUND.to_bits());
    push_text(&mut bytes, "metric=abs(delta)/max(1,sum(abs(real-input)))");
    push_len(&mut bytes, fixtures.len());
    for fixture in fixtures {
        push_len(&mut bytes, fixture.n);
        push_real_bits(&mut bytes, &fixture.input);
    }
    bytes
}

fn corruption_inputs(
    corruption: Corruption,
    canonical: &[ComplexBits; KAT_N],
    corrupted: &[ComplexBits; KAT_N],
) -> Vec<u8> {
    let kat = kat_inputs(canonical);
    let mut bytes = common_frame_prefix(b"bedrock:fs-fft:seeded-impulse-reference-corruption:v1");
    push_u64(&mut bytes, corruption.seed);
    push_len(&mut bytes, corruption.output);
    push_len(&mut bytes, corruption.component);
    push_u64(&mut bytes, u64::from(corruption.bit));
    push_nested(&mut bytes, "nested-canonical-kat", &kat);
    push_text(&mut bytes, "canonical-reference-bits");
    push_complex_bits(&mut bytes, canonical);
    push_text(&mut bytes, "corrupted-reference-bits");
    push_complex_bits(&mut bytes, corrupted);
    bytes
}

fn next_centered(state: &mut u64) -> f64 {
    *state = state
        .wrapping_mul(LCG_MULTIPLIER)
        .wrapping_add(LCG_INCREMENT);
    ((*state >> 11) as f64) / ((1_u64 << 53) as f64) - 0.5
}

fn complex_fixtures() -> Vec<ComplexFixture> {
    let mut state = COMPLEX_SEED;
    let mut fixtures = Vec::with_capacity(SIZES.len() * FIXTURES_PER_SIZE);
    for n in SIZES {
        for _ in 0..FIXTURES_PER_SIZE {
            let input = (0..n)
                .map(|_| {
                    [
                        next_centered(&mut state).to_bits(),
                        next_centered(&mut state).to_bits(),
                    ]
                })
                .collect();
            fixtures.push(ComplexFixture { n, input });
        }
    }
    fixtures
}

fn real_fixtures() -> Vec<RealFixture> {
    let mut state = REAL_SEED;
    let mut fixtures = Vec::with_capacity(SIZES.len() * FIXTURES_PER_SIZE);
    for n in SIZES {
        for _ in 0..FIXTURES_PER_SIZE {
            let input = (0..n)
                .map(|_| next_centered(&mut state).to_bits())
                .collect();
            fixtures.push(RealFixture { n, input });
        }
    }
    fixtures
}

fn oracle_options() -> FftOptions {
    FftOptions::default()
        .with_normalization(Normalization::Backward)
        .with_workers(WorkerPolicy::Exact(1))
        .with_backend(BackendKind::NaiveDft)
        .with_check_finite(true)
}

fn panic_message(payload: &(dyn core::any::Any + Send)) -> String {
    payload
        .downcast_ref::<String>()
        .cloned()
        .or_else(|| {
            payload
                .downcast_ref::<&str>()
                .map(|text| (*text).to_owned())
        })
        .unwrap_or_else(|| "non-text panic payload".to_owned())
}

fn c64_bits(value: C64) -> ComplexBits {
    [value.re.to_bits(), value.im.to_bits()]
}

fn measure_complex(fixture: &ComplexFixture) -> Measurement {
    let frankensim = catch_unwind(AssertUnwindSafe(|| {
        let plan = Fft::new(fixture.n);
        let mut data: Vec<C64> = fixture
            .input
            .iter()
            .map(|[re, im]| C64::new(f64::from_bits(*re), f64::from_bits(*im)))
            .collect();
        let mut scratch = vec![C64::default(); fixture.n];
        plan.forward(&mut data, &mut scratch);
        data.into_iter().map(c64_bits).collect()
    }))
    .map_err(|payload| format!("fs-fft panicked: {}", panic_message(&*payload)));
    let oracle_input: Vec<(f64, f64)> = fixture
        .input
        .iter()
        .map(|[re, im]| (f64::from_bits(*re), f64::from_bits(*im)))
        .collect();
    let frankenscipy = fft(&oracle_input, &oracle_options())
        .map(|output| {
            output
                .into_iter()
                .map(|(re, im)| [re.to_bits(), im.to_bits()])
                .collect()
        })
        .map_err(|error| format!("{error:?}"));
    Measurement {
        frankensim,
        frankenscipy,
    }
}

fn measure_real(fixture: &RealFixture) -> Measurement {
    let input: Vec<f64> = fixture
        .input
        .iter()
        .map(|bits| f64::from_bits(*bits))
        .collect();
    let frankensim = catch_unwind(AssertUnwindSafe(|| {
        RealFft::new(fixture.n)
            .forward(&input)
            .into_iter()
            .map(c64_bits)
            .collect()
    }))
    .map_err(|payload| format!("fs-fft panicked: {}", panic_message(&*payload)));
    // Deliberately use the full complex O(n^2) DFT and truncate its Hermitian
    // half. Calling the oracle's rfft would share fs-fft's pack-and-untangle
    // structure and could let the same algebraic bug agree with itself.
    let oracle_input: Vec<(f64, f64)> = input.iter().map(|&value| (value, 0.0)).collect();
    let frankenscipy = fft(&oracle_input, &oracle_options())
        .map(|output| {
            output
                .into_iter()
                .take(fixture.n / 2 + 1)
                .map(|(re, im)| [re.to_bits(), im.to_bits()])
                .collect()
        })
        .map_err(|error| format!("{error:?}"));
    Measurement {
        frankensim,
        frankenscipy,
    }
}

fn corruption_context(corruption: Option<Corruption>) -> String {
    corruption.map_or_else(
        || "mode=canonical".to_owned(),
        |value| {
            format!(
                "seed=0x{:016x}; output={}; component={}; bit={}",
                value.seed, value.output, value.component, value.bit,
            )
        },
    )
}

fn kat_outcome(
    reference: [ComplexBits; KAT_N],
    corruption: Option<Corruption>,
    input_frame: &[u8],
) -> CaseOutcome {
    let fixture = kat_fixture();
    let first = measure_complex(&fixture);
    let replay = measure_complex(&fixture);
    let inputs_hex = hex_bytes(input_frame);
    let context = corruption_context(corruption);
    if first != replay {
        return CaseOutcome::fail(format!(
            "{context}; stage=same-run-replay; first={first:016x?}; replay={replay:016x?}; inputs_hex={inputs_hex}"
        ))
        .with_evidence("crates/fs-fft/CONTRACT.md#determinism-class");
    }
    let (Ok(frankensim), Ok(frankenscipy)) = (&first.frankensim, &first.frankenscipy) else {
        return CaseOutcome::fail(format!(
            "{context}; stage=impulse-execution; measurement={first:016x?}; inputs_hex={inputs_hex}"
        ))
        .with_evidence("crates/fs-fft/CONTRACT.md#error-model");
    };
    if frankensim.len() != KAT_N || frankenscipy.len() != KAT_N {
        return CaseOutcome::fail(format!(
            "{context}; stage=impulse-output-length; expected={KAT_N}; frankensim_len={}; frankenscipy_len={}; frankensim={frankensim:016x?}; frankenscipy={frankenscipy:016x?}; inputs_hex={inputs_hex}",
            frankensim.len(), frankenscipy.len(),
        ))
        .with_evidence("crates/fs-fft/CONTRACT.md#invariants");
    }
    if frankensim.as_slice() != reference.as_slice()
        || frankenscipy.as_slice() != reference.as_slice()
    {
        return CaseOutcome::fail(format!(
            "{context}; stage=unit-impulse-known-answer; frankensim={frankensim:016x?}; frankenscipy={frankenscipy:016x?}; reference={reference:016x?}; inputs_hex={inputs_hex}"
        ))
        .with_evidence("crates/fs-fft/CONTRACT.md#invariants")
        .with_evidence("constellation.lock:frankenscipy-0.1.0");
    }
    CaseOutcome::pass(format!(
        "n={KAT_N}; impulse=delta_0; outputs_per_backend={KAT_N}; frankensim_output_digest={:016x}; frankenscipy_output_digest={:016x}; same_run=identical",
        digest_complex_bits(frankensim),
        digest_complex_bits(frankenscipy),
    ))
    .with_evidence("crates/fs-fft/CONTRACT.md#invariants")
    .with_evidence("constellation.lock:frankenscipy-0.1.0")
}

fn complex_input_scale(fixture: &ComplexFixture) -> f64 {
    fixture
        .input
        .iter()
        .map(|[re, im]| f64::from_bits(*re).hypot(f64::from_bits(*im)))
        .sum::<f64>()
        .max(1.0)
}

fn real_input_scale(fixture: &RealFixture) -> f64 {
    fixture
        .input
        .iter()
        .map(|bits| f64::from_bits(*bits).abs())
        .sum::<f64>()
        .max(1.0)
}

fn scaled_complex_error(actual: ComplexBits, reference: ComplexBits, input_scale: f64) -> f64 {
    let actual = (f64::from_bits(actual[0]), f64::from_bits(actual[1]));
    let reference = (f64::from_bits(reference[0]), f64::from_bits(reference[1]));
    let delta = (actual.0 - reference.0).hypot(actual.1 - reference.1);
    delta / input_scale
}

fn bounded_outcome(
    kind: &str,
    expected_len: impl Fn(usize) -> usize,
    measurements: impl Fn() -> Vec<Measurement>,
    fixture_sizes: &[usize],
    input_scales: &[f64],
    bound: f64,
    input_frame: &[u8],
) -> CaseOutcome {
    let first = measurements();
    let replay = measurements();
    let inputs_hex = hex_bytes(input_frame);
    if first != replay {
        return CaseOutcome::fail(format!(
            "stage=same-run-{kind}-replay; first={first:016x?}; replay={replay:016x?}; inputs_hex={inputs_hex}"
        ))
        .with_evidence("crates/fs-fft/CONTRACT.md#determinism-class");
    }
    if first.len() != fixture_sizes.len() || first.len() != input_scales.len() {
        return CaseOutcome::fail(format!(
            "stage={kind}-fixture-count; expected={}; scales={}; actual={}; measurements={first:016x?}; inputs_hex={inputs_hex}",
            fixture_sizes.len(), input_scales.len(), first.len(),
        ))
        .with_evidence("crates/fs-fft/CONTRACT.md#invariants");
    }

    let mut max_scaled = 0.0_f64;
    let mut frankensim_outputs = Vec::new();
    let mut frankenscipy_outputs = Vec::new();
    for (fixture, ((measurement, &n), &input_scale)) in first
        .iter()
        .zip(fixture_sizes)
        .zip(input_scales)
        .enumerate()
    {
        let (Ok(frankensim), Ok(frankenscipy)) =
            (&measurement.frankensim, &measurement.frankenscipy)
        else {
            return CaseOutcome::fail(format!(
                "stage={kind}-execution; fixture={fixture}; n={n}; measurement={measurement:016x?}; inputs_hex={inputs_hex}"
            ))
            .with_evidence("crates/fs-fft/CONTRACT.md#error-model");
        };
        let expected = expected_len(n);
        if frankensim.len() != expected || frankenscipy.len() != expected {
            return CaseOutcome::fail(format!(
                "stage={kind}-output-length; fixture={fixture}; n={n}; expected={expected}; frankensim_len={}; frankenscipy_len={}; frankensim={frankensim:016x?}; frankenscipy={frankenscipy:016x?}; inputs_hex={inputs_hex}",
                frankensim.len(), frankenscipy.len(),
            ))
            .with_evidence("crates/fs-fft/CONTRACT.md#invariants");
        }
        for (bin, (&actual, &reference)) in frankensim.iter().zip(frankenscipy).enumerate() {
            let scaled = scaled_complex_error(actual, reference, input_scale);
            max_scaled = max_scaled.max(scaled);
            if !scaled.is_finite() || scaled > bound {
                return CaseOutcome::fail(format!(
                    "stage={kind}-oracle-agreement; fixture={fixture}; n={n}; bin={bin}; frankensim_bits={actual:016x?}; frankenscipy_bits={reference:016x?}; input_scale={input_scale:.17e}; scaled={scaled:.17e}; bound={bound:.17e}; inputs_hex={inputs_hex}"
                ))
                .with_evidence("crates/fs-fft/CONTRACT.md#conformance-tests")
                .with_evidence("constellation.lock:frankenscipy-0.1.0");
            }
        }
        frankensim_outputs.extend_from_slice(frankensim);
        frankenscipy_outputs.extend_from_slice(frankenscipy);
    }
    CaseOutcome::pass(format!(
        "kind={kind}; fixtures={}; sizes=4,8,16,32; scaled_bound={bound:.17e}; max_scaled={max_scaled:.17e}; frankensim_output_digest={:016x}; frankenscipy_output_digest={:016x}; same_run=identical",
        first.len(),
        digest_complex_bits(&frankensim_outputs),
        digest_complex_bits(&frankenscipy_outputs),
    ))
    .with_evidence("crates/fs-fft/CONTRACT.md#conformance-tests")
    .with_evidence("constellation.lock:frankenscipy-0.1.0")
}

#[test]
fn frankenscipy_fft_oracle_casebook_emits_replay_complete_green_records() {
    assert_eq!(CASEBOOK_RECORD_VERSION, 1);
    assert_eq!(TRANSFORM_BIT_SEMANTICS_VERSION, 1);
    let kat_frame = kat_inputs(&KAT_EXPECTED);
    let complex_fixtures = complex_fixtures();
    let real_fixtures = real_fixtures();
    let complex_frame = complex_inputs(&complex_fixtures);
    let real_frame = real_inputs(&real_fixtures);
    let kat_digest = fnv1a64(&kat_frame);
    let complex_digest = fnv1a64(&complex_frame);
    let real_digest = fnv1a64(&real_frame);
    assert_eq!(
        (kat_frame.len(), kat_digest),
        (1_069, 0x4bf7_d7d8_797e_8598)
    );
    assert_eq!(
        (complex_frame.len(), complex_digest),
        (2_948, 0xc6d0_f04b_076b_ffb1)
    );
    assert_eq!(
        (real_frame.len(), real_digest),
        (1_982, 0x0208_2957_e971_d951)
    );

    let complex_sizes: Vec<usize> = complex_fixtures.iter().map(|fixture| fixture.n).collect();
    let real_sizes: Vec<usize> = real_fixtures.iter().map(|fixture| fixture.n).collect();
    let complex_scales: Vec<f64> = complex_fixtures.iter().map(complex_input_scale).collect();
    let real_scales: Vec<f64> = real_fixtures.iter().map(real_input_scale).collect();
    let report = Suite::new(SUITE)
        .case(
            "unit-impulse-complex-known-answer",
            kat_digest,
            ToleranceSpec::Exact,
            move || kat_outcome(KAT_EXPECTED, None, &kat_frame),
        )
        .case(
            "bounded-complex-naive-oracle-agreement",
            complex_digest,
            ToleranceSpec::RelativeLe(COMPLEX_BOUND),
            move || {
                bounded_outcome(
                    "complex-fft",
                    |n| n,
                    || complex_fixtures.iter().map(measure_complex).collect(),
                    &complex_sizes,
                    &complex_scales,
                    COMPLEX_BOUND,
                    &complex_frame,
                )
            },
        )
        .case(
            "bounded-real-naive-oracle-agreement",
            real_digest,
            ToleranceSpec::RelativeLe(REAL_BOUND),
            move || {
                bounded_outcome(
                    "real-rfft",
                    |n| n / 2 + 1,
                    || real_fixtures.iter().map(measure_real).collect(),
                    &real_sizes,
                    &real_scales,
                    REAL_BOUND,
                    &real_frame,
                )
            },
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
            "unit-impulse-complex-known-answer",
            "bounded-complex-naive-oracle-agreement",
            "bounded-real-naive-oracle-agreement",
        ]
    );
    assert!(
        report.records[0]
            .details
            .contains("frankensim_output_digest=")
    );
    assert!(
        report.records[0]
            .details
            .contains("frankenscipy_output_digest=")
    );
}

#[test]
fn disclosed_seeded_impulse_reference_corruption_turns_suite_red() {
    const CORRUPTION_SEED: u64 = 0xF17F_0000;
    let output = (CORRUPTION_SEED & 0x7) as usize;
    let component = ((CORRUPTION_SEED >> 3) & 0x1) as usize;
    let bit = ((CORRUPTION_SEED >> 4) & 0x3f) as u32;
    assert_eq!((output, component, bit), (0, 0, 0));
    let corruption = Corruption {
        seed: CORRUPTION_SEED,
        output,
        component,
        bit,
    };
    let mut corrupted = KAT_EXPECTED;
    corrupted[output][component] ^= 1_u64 << bit;
    let inputs = corruption_inputs(corruption, &KAT_EXPECTED, &corrupted);
    let inputs_digest = fnv1a64(&inputs);
    assert_eq!(
        (inputs.len(), inputs_digest),
        (2_141, 0x5c17_6ee7_569e_7b16)
    );

    let make_report = || {
        let input_frame = inputs.clone();
        Suite::new(SUITE)
            .case(
                "seeded-impulse-expected-reference-corruption",
                inputs_digest,
                ToleranceSpec::Exact,
                move || kat_outcome(corrupted, Some(corruption), &input_frame),
            )
            .run()
    };
    let first = make_report();
    let replay = make_report();
    let first_failures = first.failures();
    let replay_failures = replay.failures();
    let [first_failure] = first_failures.as_slice() else {
        panic!("the disclosed impulse corruption must produce exactly one failure");
    };
    let [replay_failure] = replay_failures.as_slice() else {
        panic!("the replayed impulse corruption must produce exactly one failure");
    };
    assert_eq!(first_failure.json_line(), replay_failure.json_line());
    assert!(
        first_failure
            .details
            .contains("stage=unit-impulse-known-answer")
    );
    assert!(
        first_failure
            .details
            .contains(&format!("seed=0x{CORRUPTION_SEED:016x}"))
    );
    assert!(
        first_failure
            .details
            .contains("output=0; component=0; bit=0")
    );
    assert!(first_failure.details.contains("inputs_hex="));
    assert!(
        first_failure
            .json_line()
            .contains("\"tolerance\":\"exact\",\"pass\":false")
    );
    let panic = catch_unwind(|| first.assert_green())
        .expect_err("the merge gate must reject the disclosed impulse corruption");
    let message = panic
        .downcast_ref::<String>()
        .map(String::as_str)
        .or_else(|| panic.downcast_ref::<&str>().copied())
        .expect("casebook panic carries text");
    assert!(message.contains("seeded-impulse-expected-reference-corruption"));
    assert!(message.contains(&format!("seed=0x{CORRUPTION_SEED:016x}")));
}
