//! FrankenScipy sparse-baseline oracle evidence for BEDROCK.
//!
//! Identical canonical CSR fixtures execute through FrankenSim's fused,
//! fixed-order SpMV and FrankenScipy's SciPy-compatible reference matvec.
//! The dyadic fixture is bit-exact; the seeded finite battery uses a declared
//! absolute oracle-agreement margin because the two implementations deliberately
//! differ in fused versus separately rounded multiply-add arithmetic.
//!
//! This is portable G0 oracle evidence. It makes no throughput, full-G5,
//! fresh dual-ISA, or production dependency claim: `fsci-sparse` is test-only.

use core::fmt::Write as _;
use std::panic::{AssertUnwindSafe, catch_unwind};

use fs_casebook::{CASEBOOK_RECORD_VERSION, CaseOutcome, Suite, ToleranceSpec, fnv1a64};
use fs_sparse::Csr;
use fsci_sparse::{CsrMatrix, Shape2D, SparseError};

const SUITE: &str = "bedrock/fs-sparse-frankenscipy-oracle-v1";
const ORACLE_VERSION: &str = "fsci-sparse/0.1.0";
const ORACLE_API: &str = "CsrMatrix::from_components+matvec:v1";
const KAT_NROWS: usize = 4;
const KAT_NCOLS: usize = 5;
const KAT_ROW_PTR: [usize; 5] = [0, 3, 3, 5, 8];
const KAT_COLUMNS: [usize; 8] = [0, 2, 4, 1, 3, 0, 2, 4];
const KAT_VALUES: [f64; 8] = [2.0, -0.5, 4.0, -1.5, 0.25, 0.125, -2.0, 3.0];
const KAT_X: [f64; 5] = [1.0, -2.0, 0.5, 4.0, -1.0];
const KAT_EXPECTED_BITS: [u64; 4] = [
    0xc002_0000_0000_0000,
    0,
    0x4010_0000_0000_0000,
    0xc00f_0000_0000_0000,
];
const SEEDED_ROOT: u64 = 0x5A25_5A25_D15C_0001;
const SEEDED_CASES: usize = 12;
const SEEDED_NROWS: usize = 7;
const SEEDED_NCOLS: usize = 11;
const SEEDED_NNZ_PER_ROW: usize = 5;
const SEEDED_ABS_BOUND: f64 = 32.0 * f64::EPSILON;
const POISON_BITS: u64 = 0x7ff8_0000_0000_0000;

#[derive(Debug, Clone)]
struct Fixture {
    nrows: usize,
    ncols: usize,
    row_ptr: Vec<usize>,
    columns: Vec<usize>,
    values: Vec<f64>,
    x: Vec<f64>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct Measurement {
    frankensim: Vec<u64>,
    frankenscipy: Vec<u64>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct RefusalMeasurement {
    frankensim_panicked: bool,
    frankensim_message: String,
    frankensim_output_bits: Vec<u64>,
    frankenscipy_typed_shape_error: bool,
    frankenscipy_message: String,
}

#[derive(Debug, Clone, Copy)]
struct Corruption {
    seed: u64,
    output: usize,
    bit: u32,
}

fn push_u64(bytes: &mut Vec<u8>, value: u64) {
    bytes.extend_from_slice(&value.to_le_bytes());
}

fn push_len(bytes: &mut Vec<u8>, value: usize) {
    push_u64(
        bytes,
        u64::try_from(value).expect("oracle fixture lengths fit u64"),
    );
}

fn push_text(bytes: &mut Vec<u8>, value: &str) {
    push_len(bytes, value.len());
    bytes.extend_from_slice(value.as_bytes());
}

fn push_usizes(bytes: &mut Vec<u8>, values: &[usize]) {
    push_len(bytes, values.len());
    for &value in values {
        push_len(bytes, value);
    }
}

fn push_f64s(bytes: &mut Vec<u8>, values: &[f64]) {
    push_len(bytes, values.len());
    for value in values {
        push_u64(bytes, value.to_bits());
    }
}

fn push_u64s(bytes: &mut Vec<u8>, values: &[u64]) {
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

fn push_fixture(bytes: &mut Vec<u8>, fixture: &Fixture) {
    push_len(bytes, fixture.nrows);
    push_len(bytes, fixture.ncols);
    push_text(bytes, "row-pointers");
    push_usizes(bytes, &fixture.row_ptr);
    push_text(bytes, "column-indices");
    push_usizes(bytes, &fixture.columns);
    push_text(bytes, "value-bits");
    push_f64s(bytes, &fixture.values);
    push_text(bytes, "x-bits");
    push_f64s(bytes, &fixture.x);
}

fn hex_bytes(bytes: &[u8]) -> String {
    let mut encoded = String::with_capacity(bytes.len() * 2);
    for byte in bytes {
        write!(&mut encoded, "{byte:02x}").expect("writing to String cannot fail");
    }
    encoded
}

fn kat_fixture() -> Fixture {
    Fixture {
        nrows: KAT_NROWS,
        ncols: KAT_NCOLS,
        row_ptr: KAT_ROW_PTR.to_vec(),
        columns: KAT_COLUMNS.to_vec(),
        values: KAT_VALUES.to_vec(),
        x: KAT_X.to_vec(),
    }
}

fn next_u64(state: &mut u64) -> u64 {
    *state = state
        .wrapping_mul(6_364_136_223_846_793_005)
        .wrapping_add(1_442_695_040_888_963_407);
    *state
}

fn next_centered(state: &mut u64) -> f64 {
    ((next_u64(state) >> 11) as f64) / ((1_u64 << 53) as f64) - 0.5
}

fn seeded_fixtures() -> Vec<Fixture> {
    let mut state = SEEDED_ROOT;
    let mut fixtures = Vec::with_capacity(SEEDED_CASES);
    for _ in 0..SEEDED_CASES {
        let mut row_ptr = Vec::with_capacity(SEEDED_NROWS + 1);
        let mut columns = Vec::with_capacity(SEEDED_NROWS * SEEDED_NNZ_PER_ROW);
        let mut values = Vec::with_capacity(columns.capacity());
        row_ptr.push(0);
        for _ in 0..SEEDED_NROWS {
            let base = (next_u64(&mut state) % SEEDED_NCOLS as u64) as usize;
            let stride = 1 + (next_u64(&mut state) % (SEEDED_NCOLS as u64 - 1)) as usize;
            let mut row: Vec<(usize, f64)> = (0..SEEDED_NNZ_PER_ROW)
                .map(|offset| {
                    (
                        (base + offset * stride) % SEEDED_NCOLS,
                        next_centered(&mut state),
                    )
                })
                .collect();
            row.sort_unstable_by_key(|(column, _)| *column);
            for (column, value) in row {
                columns.push(column);
                values.push(value);
            }
            row_ptr.push(columns.len());
        }
        let x = (0..SEEDED_NCOLS)
            .map(|_| next_centered(&mut state))
            .collect();
        fixtures.push(Fixture {
            nrows: SEEDED_NROWS,
            ncols: SEEDED_NCOLS,
            row_ptr,
            columns,
            values,
            x,
        });
    }
    fixtures
}

fn panic_message(payload: &(dyn core::any::Any + Send)) -> String {
    payload
        .downcast_ref::<String>()
        .cloned()
        .or_else(|| {
            payload
                .downcast_ref::<&str>()
                .map(|value| (*value).to_owned())
        })
        .unwrap_or_else(|| "non-text panic payload".to_owned())
}

fn measure(fixture: &Fixture) -> Result<Measurement, String> {
    let frankensim = catch_unwind(AssertUnwindSafe(|| {
        let matrix = Csr::from_parts(
            fixture.nrows,
            fixture.ncols,
            fixture.row_ptr.clone(),
            fixture.columns.clone(),
            fixture.values.clone(),
        );
        let mut output = vec![f64::from_bits(POISON_BITS); fixture.nrows];
        matrix.spmv(&fixture.x, &mut output);
        output
    }))
    .map_err(|payload| format!("fs-sparse panicked: {}", panic_message(&*payload)))?;

    let oracle = CsrMatrix::from_components(
        Shape2D::new(fixture.nrows, fixture.ncols),
        fixture.values.clone(),
        fixture.columns.clone(),
        fixture.row_ptr.clone(),
        false,
    )
    .and_then(|matrix| matrix.matvec(&fixture.x))
    .map_err(|error| format!("{ORACLE_VERSION} refused canonical fixture: {error:?}"))?;

    Ok(Measurement {
        frankensim: frankensim.into_iter().map(f64::to_bits).collect(),
        frankenscipy: oracle.into_iter().map(f64::to_bits).collect(),
    })
}

fn common_frame_prefix(domain: &[u8]) -> Vec<u8> {
    let mut bytes = domain.to_vec();
    push_text(&mut bytes, "fs-sparse-version");
    push_text(&mut bytes, env!("CARGO_PKG_VERSION"));
    push_text(&mut bytes, "oracle-version");
    push_text(&mut bytes, ORACLE_VERSION);
    push_text(&mut bytes, "oracle-api");
    push_text(&mut bytes, ORACLE_API);
    bytes
}

fn kat_inputs() -> Vec<u8> {
    let mut bytes = common_frame_prefix(b"bedrock:fs-sparse:frankenscipy-spmv-kat:v1");
    push_text(
        &mut bytes,
        "canonical-csr;ascending-columns;fs-sparse=fused;fsci=mul-then-add",
    );
    push_fixture(&mut bytes, &kat_fixture());
    push_text(&mut bytes, "expected-y-bits");
    push_u64s(&mut bytes, &KAT_EXPECTED_BITS);
    push_text(&mut bytes, "dyadic-products-and-sums-exact:v1");
    bytes
}

fn seeded_inputs() -> Vec<u8> {
    let mut bytes = common_frame_prefix(b"bedrock:fs-sparse:frankenscipy-seeded-spmv:v1");
    push_text(&mut bytes, "lcg64-v1");
    push_u64(&mut bytes, SEEDED_ROOT);
    push_u64(&mut bytes, 6_364_136_223_846_793_005);
    push_u64(&mut bytes, 1_442_695_040_888_963_407);
    push_len(&mut bytes, SEEDED_CASES);
    push_len(&mut bytes, SEEDED_NROWS);
    push_len(&mut bytes, SEEDED_NCOLS);
    push_len(&mut bytes, SEEDED_NNZ_PER_ROW);
    push_text(&mut bytes, "absolute-oracle-agreement-bound");
    push_u64(&mut bytes, SEEDED_ABS_BOUND.to_bits());
    let fixtures = seeded_fixtures();
    push_len(&mut bytes, fixtures.len());
    for fixture in &fixtures {
        push_fixture(&mut bytes, fixture);
    }
    bytes
}

fn refusal_inputs() -> Vec<u8> {
    let mut bytes = common_frame_prefix(b"bedrock:fs-sparse:oracle-shape-refusal:v1");
    push_fixture(&mut bytes, &kat_fixture());
    push_text(&mut bytes, "malformed-x-bits");
    push_f64s(&mut bytes, &[1.0, 2.0, 3.0, 4.0]);
    push_text(&mut bytes, "frankensim-contract");
    push_text(&mut bytes, "panic:spmv-x-length-must-equal-ncols-5");
    push_text(&mut bytes, "frankenscipy-contract");
    push_text(
        &mut bytes,
        "SparseError::IncompatibleShape{x length 4 != cols 5}",
    );
    push_text(&mut bytes, "frankensim-output-before-refusal");
    push_u64s(&mut bytes, &[POISON_BITS; KAT_NROWS]);
    bytes
}

fn exact_outcome(
    reference: [u64; KAT_NROWS],
    corruption: Option<Corruption>,
    input_frame: &[u8],
) -> CaseOutcome {
    let first = measure(&kat_fixture());
    let replay = measure(&kat_fixture());
    let inputs_hex = hex_bytes(input_frame);
    let context = corruption.map_or_else(
        || "mode=canonical".to_owned(),
        |corruption| {
            format!(
                "seed=0x{:016x}; output={}; bit={}",
                corruption.seed, corruption.output, corruption.bit
            )
        },
    );
    if first != replay {
        return CaseOutcome::fail(format!(
            "{context}; stage=same-run-replay; first={first:?}; replay={replay:?}; inputs_hex={inputs_hex}"
        ))
        .with_evidence("crates/fs-sparse/CONTRACT.md#determinism-class");
    }
    let Ok(measurement) = first else {
        return CaseOutcome::fail(format!(
            "{context}; stage=canonical-construction-or-execution; result={first:?}; inputs_hex={inputs_hex}"
        ))
        .with_evidence("crates/fs-sparse/CONTRACT.md#error-model");
    };
    if measurement.frankensim.as_slice() != reference.as_slice()
        || measurement.frankenscipy.as_slice() != reference.as_slice()
    {
        return CaseOutcome::fail(format!(
            "{context}; stage=dyadic-known-answer; frankensim={:016x?}; frankenscipy={:016x?}; reference={reference:016x?}; inputs_hex={inputs_hex}",
            measurement.frankensim, measurement.frankenscipy,
        ))
        .with_evidence("crates/fs-sparse/CONTRACT.md#invariants")
        .with_evidence("constellation.lock:frankenscipy-0.1.0");
    }
    CaseOutcome::pass(
        "shape=4x5; nnz=8; y=[-2.25,+0,4,-3.875]; fs-sparse=exact; fsci-sparse=exact; same_run=identical",
    )
    .with_evidence("crates/fs-sparse/CONTRACT.md#invariants")
    .with_evidence("constellation.lock:frankenscipy-0.1.0")
}

fn ordered_bits(value: f64) -> u64 {
    let bits = value.to_bits();
    if bits & (1_u64 << 63) == 0 {
        bits | (1_u64 << 63)
    } else {
        !bits
    }
}

fn ulp_distance(left: f64, right: f64) -> u64 {
    ordered_bits(left).abs_diff(ordered_bits(right))
}

fn seeded_outcome(input_frame: &[u8]) -> CaseOutcome {
    let fixtures = seeded_fixtures();
    let first: Result<Vec<_>, _> = fixtures.iter().map(measure).collect();
    let replay: Result<Vec<_>, _> = fixtures.iter().map(measure).collect();
    let inputs_hex = hex_bytes(input_frame);
    if first != replay {
        return CaseOutcome::fail(format!(
            "stage=same-run-replay; first={first:?}; replay={replay:?}; inputs_hex={inputs_hex}"
        ))
        .with_evidence("crates/fs-sparse/CONTRACT.md#determinism-class");
    }
    let Ok(measurements) = first else {
        return CaseOutcome::fail(format!(
            "stage=seeded-construction-or-execution; result={first:?}; inputs_hex={inputs_hex}"
        ))
        .with_evidence("crates/fs-sparse/CONTRACT.md#error-model");
    };

    let mut max_abs = 0.0_f64;
    let mut max_ulps = 0_u64;
    let mut frankensim_output_bytes = Vec::with_capacity(SEEDED_CASES * SEEDED_NROWS * 8);
    let mut frankenscipy_output_bytes = Vec::with_capacity(SEEDED_CASES * SEEDED_NROWS * 8);
    for (case_index, measurement) in measurements.iter().enumerate() {
        let expected_rows = fixtures[case_index].nrows;
        if measurement.frankensim.len() != expected_rows
            || measurement.frankenscipy.len() != expected_rows
        {
            return CaseOutcome::fail(format!(
                "stage=seeded-output-shape; seed=0x{SEEDED_ROOT:016x}; case={case_index}; expected_rows={expected_rows}; frankensim_rows={}; frankenscipy_rows={}; fixture={:?}; inputs_hex={inputs_hex}",
                measurement.frankensim.len(),
                measurement.frankenscipy.len(),
                fixtures[case_index],
            ))
            .with_evidence("crates/fs-sparse/CONTRACT.md#invariants")
            .with_evidence("constellation.lock:frankenscipy-0.1.0");
        }
        for (row, (&frankensim_bits, &frankenscipy_bits)) in measurement
            .frankensim
            .iter()
            .zip(&measurement.frankenscipy)
            .enumerate()
        {
            let frankensim = f64::from_bits(frankensim_bits);
            let frankenscipy = f64::from_bits(frankenscipy_bits);
            let abs = (frankensim - frankenscipy).abs();
            let ulps = ulp_distance(frankensim, frankenscipy);
            frankensim_output_bytes.extend_from_slice(&frankensim_bits.to_le_bytes());
            frankenscipy_output_bytes.extend_from_slice(&frankenscipy_bits.to_le_bytes());
            max_abs = max_abs.max(abs);
            max_ulps = max_ulps.max(ulps);
            if !abs.is_finite() || abs > SEEDED_ABS_BOUND {
                return CaseOutcome::fail(format!(
                    "stage=seeded-oracle-margin; seed=0x{SEEDED_ROOT:016x}; case={case_index}; row={row}; frankensim_bits=0x{frankensim_bits:016x}; frankenscipy_bits=0x{frankenscipy_bits:016x}; abs={abs:.17e}; abs_bound={SEEDED_ABS_BOUND:.17e}; ulps={ulps}; fixture={:?}; inputs_hex={inputs_hex}",
                    fixtures[case_index],
                ))
                .with_evidence("crates/fs-sparse/CONTRACT.md#invariants")
                .with_evidence("constellation.lock:frankenscipy-0.1.0");
            }
        }
    }
    let frankensim_digest = fnv1a64(&frankensim_output_bytes);
    let frankenscipy_digest = fnv1a64(&frankenscipy_output_bytes);

    CaseOutcome::pass(format!(
        "seed=0x{SEEDED_ROOT:016x}; cases={SEEDED_CASES}; outputs={}; frankensim_output_digest={frankensim_digest:016x}; frankenscipy_output_digest={frankenscipy_digest:016x}; max_abs={max_abs:.17e}; abs_bound={SEEDED_ABS_BOUND:.17e}; max_ulps={max_ulps}; same_run=identical",
        SEEDED_CASES * SEEDED_NROWS,
    ))
    .with_evidence("crates/fs-sparse/CONTRACT.md#invariants")
    .with_evidence("constellation.lock:frankenscipy-0.1.0")
}

fn refusal_measurement() -> Result<RefusalMeasurement, String> {
    let fixture = kat_fixture();
    let matrix = catch_unwind(AssertUnwindSafe(|| {
        Csr::from_parts(
            fixture.nrows,
            fixture.ncols,
            fixture.row_ptr.clone(),
            fixture.columns.clone(),
            fixture.values.clone(),
        )
    }))
    .map_err(|payload| {
        format!(
            "fs-sparse canonical construction panicked: {}",
            panic_message(&*payload)
        )
    })?;
    let malformed_x = [1.0, 2.0, 3.0, 4.0];
    let mut output = [f64::from_bits(POISON_BITS); KAT_NROWS];
    let primary_refusal = catch_unwind(AssertUnwindSafe(|| {
        matrix.spmv(&malformed_x, &mut output);
    }));
    let (frankensim_panicked, frankensim_message) = match primary_refusal {
        Ok(()) => (false, "accepted malformed vector".to_owned()),
        Err(payload) => (true, panic_message(&*payload)),
    };

    let oracle = CsrMatrix::from_components(
        Shape2D::new(fixture.nrows, fixture.ncols),
        fixture.values,
        fixture.columns,
        fixture.row_ptr,
        false,
    )
    .map_err(|error| format!("{ORACLE_VERSION} refused canonical fixture: {error:?}"))?;
    let (frankenscipy_typed_shape_error, frankenscipy_message) = match oracle.matvec(&malformed_x) {
        Err(SparseError::IncompatibleShape { message }) => {
            let exact = message == "x length 4 != cols 5";
            (exact, message)
        }
        Err(error) => (false, format!("wrong error: {error:?}")),
        Ok(output) => (false, format!("accepted malformed vector: {output:?}")),
    };

    Ok(RefusalMeasurement {
        frankensim_panicked,
        frankensim_message,
        frankensim_output_bits: output.map(f64::to_bits).to_vec(),
        frankenscipy_typed_shape_error,
        frankenscipy_message,
    })
}

fn refusal_outcome(input_frame: &[u8]) -> CaseOutcome {
    let first = refusal_measurement();
    let replay = refusal_measurement();
    let inputs_hex = hex_bytes(input_frame);
    if first != replay {
        return CaseOutcome::fail(format!(
            "stage=same-run-refusal-replay; first={first:?}; replay={replay:?}; inputs_hex={inputs_hex}"
        ))
        .with_evidence("crates/fs-sparse/CONTRACT.md#error-model");
    }
    let Ok(measurement) = first else {
        return CaseOutcome::fail(format!(
            "stage=refusal-setup; result={first:?}; inputs_hex={inputs_hex}"
        ))
        .with_evidence("crates/fs-sparse/CONTRACT.md#error-model");
    };
    let output_unchanged = measurement
        .frankensim_output_bits
        .iter()
        .all(|bits| *bits == POISON_BITS);
    let primary_message_is_actionable = measurement
        .frankensim_message
        .contains("spmv: x length must equal ncols 5");
    if !measurement.frankensim_panicked
        || !primary_message_is_actionable
        || !output_unchanged
        || !measurement.frankenscipy_typed_shape_error
    {
        return CaseOutcome::fail(format!(
            "stage=malformed-vector-refusal; measurement={measurement:?}; primary_message_actionable={primary_message_is_actionable}; output_unchanged={output_unchanged}; inputs_hex={inputs_hex}"
        ))
        .with_evidence("crates/fs-sparse/CONTRACT.md#error-model")
        .with_evidence("constellation.lock:frankenscipy-0.1.0");
    }

    CaseOutcome::pass(
        "x_len=4; ncols=5; fs-sparse=actionable-panic-before-output; output=unchanged; fsci-sparse=typed-IncompatibleShape; same_run=identical",
    )
    .with_evidence("crates/fs-sparse/CONTRACT.md#error-model")
    .with_evidence("constellation.lock:frankenscipy-0.1.0")
}

#[test]
fn frankenscipy_sparse_oracle_suite_emits_replay_complete_green_records() {
    assert_eq!(CASEBOOK_RECORD_VERSION, 1);
    let kat_frame = kat_inputs();
    let seeded_frame = seeded_inputs();
    let refusal_frame = refusal_inputs();
    let kat_digest = fnv1a64(&kat_frame);
    let seeded_digest = fnv1a64(&seeded_frame);
    let refusal_digest = fnv1a64(&refusal_frame);
    assert_eq!(kat_digest, 0xd7c1_b8ae_3077_7f71);
    assert_eq!(seeded_digest, 0x735e_a2d3_24e9_67bc);
    assert_eq!(refusal_digest, 0x3d11_85b7_e82a_131f);

    let report = Suite::new(SUITE)
        .case(
            "dyadic-csr-spmv-known-answer",
            kat_digest,
            ToleranceSpec::Exact,
            move || exact_outcome(KAT_EXPECTED_BITS, None, &kat_frame),
        )
        .case(
            "seeded-csr-spmv-oracle-margin",
            seeded_digest,
            ToleranceSpec::AbsoluteLe(SEEDED_ABS_BOUND),
            move || seeded_outcome(&seeded_frame),
        )
        .case(
            "malformed-vector-refusal-policy",
            refusal_digest,
            ToleranceSpec::Structural,
            move || refusal_outcome(&refusal_frame),
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
            "dyadic-csr-spmv-known-answer",
            "seeded-csr-spmv-oracle-margin",
            "malformed-vector-refusal-policy",
        ]
    );
    assert_eq!(
        report.records[0].json_line(),
        format!(
            concat!(
                "{{\"casebook\":{},\"suite\":\"bedrock/fs-sparse-frankenscipy-oracle-v1\",",
                "\"case\":\"dyadic-csr-spmv-known-answer\",\"inputs_digest\":\"d7c1b8ae30777f71\",",
                "\"tolerance\":\"exact\",\"pass\":true,",
                "\"details\":\"shape=4x5; nnz=8; y=[-2.25,+0,4,-3.875]; fs-sparse=exact; fsci-sparse=exact; same_run=identical\",",
                "\"evidence\":[\"crates/fs-sparse/CONTRACT.md#invariants\",",
                "\"constellation.lock:frankenscipy-0.1.0\"]}}"
            ),
            CASEBOOK_RECORD_VERSION,
        ),
        "the sparse-oracle record schema and field order are contract"
    );
}

#[test]
fn disclosed_seeded_oracle_reference_corruption_turns_the_suite_red() {
    const CORRUPTION_SEED: u64 = 0x5A25_0000;
    let output = (CORRUPTION_SEED & 0x3) as usize;
    let bit = ((CORRUPTION_SEED >> 2) & 0x3f) as u32;
    assert_eq!((output, bit), (0, 0));
    let corruption = Corruption {
        seed: CORRUPTION_SEED,
        output,
        bit,
    };
    let mut corrupted = KAT_EXPECTED_BITS;
    corrupted[output] ^= 1_u64 << bit;
    assert_eq!(corrupted[0], KAT_EXPECTED_BITS[0] ^ 1);

    let kat_frame = kat_inputs();
    let mut inputs = common_frame_prefix(b"bedrock:fs-sparse:seeded-oracle-corruption:v1");
    push_u64(&mut inputs, CORRUPTION_SEED);
    push_len(&mut inputs, output);
    push_u64(&mut inputs, u64::from(bit));
    push_nested(&mut inputs, "nested-dyadic-kat", &kat_frame);
    push_text(&mut inputs, "canonical-reference-bits");
    push_u64s(&mut inputs, &KAT_EXPECTED_BITS);
    push_text(&mut inputs, "corrupted-reference-bits");
    push_u64s(&mut inputs, &corrupted);
    let inputs_digest = fnv1a64(&inputs);
    assert_eq!(inputs_digest, 0x323f_3591_ddcb_0d6e);

    let make_report = || {
        let input_frame = inputs.clone();
        Suite::new(SUITE)
            .case(
                "seeded-dyadic-oracle-reference-corruption",
                inputs_digest,
                ToleranceSpec::Exact,
                move || exact_outcome(corrupted, Some(corruption), &input_frame),
            )
            .run()
    };
    let first = make_report();
    let replay = make_report();
    let first_failures = first.failures();
    let replay_failures = replay.failures();
    let [first_failure] = first_failures.as_slice() else {
        panic!("the disclosed oracle corruption must produce exactly one failure");
    };
    let [replay_failure] = replay_failures.as_slice() else {
        panic!("the replayed oracle corruption must produce exactly one failure");
    };
    assert_eq!(first_failure.json_line(), replay_failure.json_line());
    assert_eq!(first_failure.inputs_digest, "323f3591ddcb0d6e");
    assert!(first_failure.details.contains("stage=dyadic-known-answer"));
    assert!(
        first_failure
            .details
            .contains(&format!("seed=0x{CORRUPTION_SEED:016x}"))
    );
    assert!(first_failure.details.contains("inputs_hex="));
    assert!(
        first_failure
            .json_line()
            .contains("\"tolerance\":\"exact\",\"pass\":false")
    );

    let panic = catch_unwind(|| first.assert_green())
        .expect_err("the merge gate must reject the disclosed oracle corruption");
    let message = panic_message(&*panic);
    assert!(message.contains("seeded-dyadic-oracle-reference-corruption"));
    assert!(message.contains(&format!("seed=0x{CORRUPTION_SEED:016x}")));
}
