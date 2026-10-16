//! Structured G5 evidence for clone-checkpoint eigensolver replay.
//!
//! Production Lanczos and LOBPCG states are advanced to disclosed cut points,
//! cloned, resumed, and compared bit-for-bit with uninterrupted runs. The
//! complete public `EigenPair` surface (value, residual, and vector) is bound
//! into the output digest.
//!
//! This is same-build, same-ISA evidence for one finite matrix-free stencil.
//! It makes no convergence-order, accuracy, performance, cross-ISA,
//! persistence, cancellation, or caller-operator claim.

use core::fmt::Write as _;
use std::panic::{AssertUnwindSafe, catch_unwind};

use fs_casebook::{CASEBOOK_RECORD_VERSION, CaseOutcome, Suite, ToleranceSpec, fnv1a64};
use fs_la::VERSION as FS_LA_VERSION;
use fs_la::eigen::{EigenPair, LanczosState, LobpcgState, lanczos_run, lobpcg_run};
use fs_math::VERSION as FS_MATH_VERSION;

const SUITE: &str = "bedrock/fs-la-eigensolver-replay-v1";
const N: usize = 96;
const PAIRS: usize = 3;
const LANCZOS_PREFIX_PAIRS: usize = 1;
const LANCZOS_FIRST: usize = 11;
const LANCZOS_SECOND: usize = 19;
const LOBPCG_FIRST: usize = 9;
const LOBPCG_SECOND: usize = 13;
const RED_SEED: u64 = 0x6A5E_C0DE_E16E_0001;
const GREEN_FRAME_LEN: usize = 633;
const GREEN_FRAME_DIGEST: u64 = 0xb366_ba87_7b69_ccf8;
const RED_FRAME_LEN: usize = 1_037;
const RED_FRAME_DIGEST: u64 = 0xb5d1_f2ca_f03c_802f;

#[derive(Debug, Clone, PartialEq, Eq)]
struct PairBits {
    value: u64,
    residual: u64,
    vector: Vec<u64>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct LanczosSnapshot {
    dimension: usize,
    exhausted: bool,
    pairs: Vec<PairBits>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct LobpcgSnapshot {
    iterations: usize,
    pairs: Vec<PairBits>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct Evaluation {
    lanczos: LanczosSnapshot,
    lobpcg: LobpcgSnapshot,
}

fn push_u64(bytes: &mut Vec<u8>, value: u64) {
    bytes.extend_from_slice(&value.to_le_bytes());
}

fn push_len(bytes: &mut Vec<u8>, value: usize) {
    push_u64(
        bytes,
        u64::try_from(value).expect("Casebook fixture lengths fit u64"),
    );
}

fn push_text(bytes: &mut Vec<u8>, value: &str) {
    push_len(bytes, value.len());
    bytes.extend_from_slice(value.as_bytes());
}

fn push_nested(bytes: &mut Vec<u8>, label: &str, nested: &[u8]) {
    push_text(bytes, label);
    push_len(bytes, nested.len());
    bytes.extend_from_slice(nested);
}

fn hex_bytes(bytes: &[u8]) -> String {
    let mut encoded = String::with_capacity(bytes.len() * 2);
    for byte in bytes {
        write!(&mut encoded, "{byte:02x}").expect("writing to String cannot fail");
    }
    encoded
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

fn laplacian_1d(input: &[f64], output: &mut [f64]) {
    assert_eq!(input.len(), N, "fixture operator input length");
    assert_eq!(output.len(), N, "fixture operator output length");
    for index in 0..N {
        let mut value = 2.0 * input[index];
        if index > 0 {
            value -= input[index - 1];
        }
        if index + 1 < N {
            value -= input[index + 1];
        }
        output[index] = value;
    }
}

fn identity_preconditioner(residual: &[f64], output: &mut [f64]) {
    output.copy_from_slice(residual);
}

fn pair_bits(pair: &EigenPair) -> PairBits {
    PairBits {
        value: pair.value.to_bits(),
        residual: pair.residual.to_bits(),
        vector: pair.vector.iter().map(|value| value.to_bits()).collect(),
    }
}

fn pairs_bits(pairs: &[EigenPair]) -> Vec<PairBits> {
    pairs.iter().map(pair_bits).collect()
}

fn first_pair_mismatch(left: &[PairBits], right: &[PairBits]) -> Option<String> {
    if left.len() != right.len() {
        return Some(format!(
            "field=pair-count; left={}; right={}",
            left.len(),
            right.len()
        ));
    }
    for (pair, (left_pair, right_pair)) in left.iter().zip(right).enumerate() {
        if left_pair.value != right_pair.value {
            return Some(format!(
                "field=value; pair={pair}; left_bits=0x{:016x}; right_bits=0x{:016x}",
                left_pair.value, right_pair.value
            ));
        }
        if left_pair.residual != right_pair.residual {
            return Some(format!(
                "field=residual; pair={pair}; left_bits=0x{:016x}; right_bits=0x{:016x}",
                left_pair.residual, right_pair.residual
            ));
        }
        if left_pair.vector.len() != right_pair.vector.len() {
            return Some(format!(
                "field=vector-length; pair={pair}; left={}; right={}",
                left_pair.vector.len(),
                right_pair.vector.len()
            ));
        }
        if let Some((index, (&left_bits, &right_bits))) = left_pair
            .vector
            .iter()
            .zip(&right_pair.vector)
            .enumerate()
            .find(|(_, (left_bits, right_bits))| left_bits != right_bits)
        {
            return Some(format!(
                "field=vector; pair={pair}; index={index}; left_bits=0x{left_bits:016x}; right_bits=0x{right_bits:016x}"
            ));
        }
    }
    None
}

fn validate_pairs(solver: &str, pairs: &[PairBits]) -> Result<(), String> {
    if pairs.len() != PAIRS {
        return Err(format!(
            "stage=structural-evidence; solver={solver}; field=pair-count; computed={}; expected={PAIRS}",
            pairs.len()
        ));
    }
    for (pair, evidence) in pairs.iter().enumerate() {
        if evidence.vector.len() != N {
            return Err(format!(
                "stage=structural-evidence; solver={solver}; pair={pair}; field=vector-length; computed={}; expected={N}",
                evidence.vector.len()
            ));
        }
        for (field, bits) in [("value", evidence.value), ("residual", evidence.residual)] {
            if !f64::from_bits(bits).is_finite() {
                return Err(format!(
                    "stage=finite-evidence; solver={solver}; pair={pair}; field={field}; bits=0x{bits:016x}"
                ));
            }
        }
        if let Some((index, &bits)) = evidence
            .vector
            .iter()
            .enumerate()
            .find(|(_, bits)| !f64::from_bits(**bits).is_finite())
        {
            return Err(format!(
                "stage=finite-evidence; solver={solver}; pair={pair}; field=vector; index={index}; bits=0x{bits:016x}"
            ));
        }
    }
    Ok(())
}

fn run_lanczos() -> Result<LanczosSnapshot, String> {
    let mut prefix = LanczosState::new(N);
    let _ = lanczos_run(
        &laplacian_1d,
        &mut prefix,
        LANCZOS_FIRST,
        LANCZOS_PREFIX_PAIRS,
        false,
    );
    if prefix.dim() != LANCZOS_FIRST || prefix.exhausted() {
        return Err(format!(
            "stage=genuine-split; solver=lanczos; dimension={}; expected={LANCZOS_FIRST}; exhausted={}",
            prefix.dim(),
            prefix.exhausted()
        ));
    }

    let mut resumed = prefix.clone();
    let resumed_pairs = pairs_bits(&lanczos_run(
        &laplacian_1d,
        &mut resumed,
        LANCZOS_SECOND,
        PAIRS,
        false,
    ));
    let mut straight = LanczosState::new(N);
    let straight_pairs = pairs_bits(&lanczos_run(
        &laplacian_1d,
        &mut straight,
        LANCZOS_FIRST + LANCZOS_SECOND,
        PAIRS,
        false,
    ));
    let expected_dimension = LANCZOS_FIRST + LANCZOS_SECOND;
    if resumed.dim() != expected_dimension || straight.dim() != expected_dimension {
        return Err(format!(
            "stage=genuine-second-segment; solver=lanczos; requested_first={LANCZOS_FIRST}; requested_second={LANCZOS_SECOND}; expected_dimension={expected_dimension}; split_dimension={}; straight_dimension={}; split_exhausted={}; straight_exhausted={}",
            resumed.dim(),
            straight.dim(),
            resumed.exhausted(),
            straight.exhausted()
        ));
    }
    if resumed.dim() != straight.dim() || resumed.exhausted() != straight.exhausted() {
        return Err(format!(
            "stage=split-vs-straight; solver=lanczos; field=state; split_dimension={}; straight_dimension={}; split_exhausted={}; straight_exhausted={}",
            resumed.dim(),
            straight.dim(),
            resumed.exhausted(),
            straight.exhausted()
        ));
    }
    if let Some(mismatch) = first_pair_mismatch(&resumed_pairs, &straight_pairs) {
        return Err(format!(
            "stage=split-vs-straight; solver=lanczos; {mismatch}"
        ));
    }
    validate_pairs("lanczos", &resumed_pairs)?;
    Ok(LanczosSnapshot {
        dimension: resumed.dim(),
        exhausted: resumed.exhausted(),
        pairs: resumed_pairs,
    })
}

fn run_lobpcg() -> Result<LobpcgSnapshot, String> {
    let mut prefix = LobpcgState::new(N, PAIRS);
    let _ = lobpcg_run(
        &laplacian_1d,
        &mut prefix,
        LOBPCG_FIRST,
        false,
        &identity_preconditioner,
    );
    if prefix.iters != LOBPCG_FIRST {
        return Err(format!(
            "stage=genuine-split; solver=lobpcg; iterations={}; expected={LOBPCG_FIRST}",
            prefix.iters
        ));
    }

    let mut resumed = prefix.clone();
    let resumed_pairs = pairs_bits(&lobpcg_run(
        &laplacian_1d,
        &mut resumed,
        LOBPCG_SECOND,
        false,
        &identity_preconditioner,
    ));
    let mut straight = LobpcgState::new(N, PAIRS);
    let straight_pairs = pairs_bits(&lobpcg_run(
        &laplacian_1d,
        &mut straight,
        LOBPCG_FIRST + LOBPCG_SECOND,
        false,
        &identity_preconditioner,
    ));
    let expected_iterations = LOBPCG_FIRST + LOBPCG_SECOND;
    if resumed.iters != expected_iterations || straight.iters != expected_iterations {
        return Err(format!(
            "stage=genuine-second-segment; solver=lobpcg; requested_first={LOBPCG_FIRST}; requested_second={LOBPCG_SECOND}; expected_iterations={expected_iterations}; split_iterations={}; straight_iterations={}",
            resumed.iters, straight.iters
        ));
    }
    if resumed.iters != straight.iters {
        return Err(format!(
            "stage=split-vs-straight; solver=lobpcg; field=iterations; split={}; straight={}",
            resumed.iters, straight.iters
        ));
    }
    if let Some(mismatch) = first_pair_mismatch(&resumed_pairs, &straight_pairs) {
        return Err(format!(
            "stage=split-vs-straight; solver=lobpcg; {mismatch}"
        ));
    }
    validate_pairs("lobpcg", &resumed_pairs)?;
    Ok(LobpcgSnapshot {
        iterations: resumed.iters,
        pairs: resumed_pairs,
    })
}

fn evaluate() -> Result<Evaluation, String> {
    catch_unwind(AssertUnwindSafe(|| {
        Ok(Evaluation {
            lanczos: run_lanczos()?,
            lobpcg: run_lobpcg()?,
        })
    }))
    .map_err(|payload| {
        format!(
            "stage=eigensolver-execution; panic={}",
            panic_message(&*payload)
        )
    })?
}

fn push_pairs(bytes: &mut Vec<u8>, pairs: &[PairBits]) {
    push_len(bytes, pairs.len());
    for pair in pairs {
        push_u64(bytes, pair.value);
        push_u64(bytes, pair.residual);
        push_len(bytes, pair.vector.len());
        for &bits in &pair.vector {
            push_u64(bytes, bits);
        }
    }
}

fn evaluation_frame(evaluation: &Evaluation) -> Vec<u8> {
    let mut bytes = b"bedrock:fs-la-eigensolver-output:v1".to_vec();
    push_text(&mut bytes, "lanczos-dimension-exhausted-pairs");
    push_len(&mut bytes, evaluation.lanczos.dimension);
    push_u64(
        &mut bytes,
        u64::from(u8::from(evaluation.lanczos.exhausted)),
    );
    push_pairs(&mut bytes, &evaluation.lanczos.pairs);
    push_text(&mut bytes, "lobpcg-iterations-pairs");
    push_len(&mut bytes, evaluation.lobpcg.iterations);
    push_pairs(&mut bytes, &evaluation.lobpcg.pairs);
    bytes
}

fn evaluation_digest(evaluation: &Evaluation) -> u64 {
    fnv1a64(&evaluation_frame(evaluation))
}

fn first_evaluation_mismatch(left: &Evaluation, right: &Evaluation) -> Option<String> {
    if left.lanczos.dimension != right.lanczos.dimension {
        return Some(format!(
            "solver=lanczos; field=dimension; first={}; replay={}",
            left.lanczos.dimension, right.lanczos.dimension
        ));
    }
    if left.lanczos.exhausted != right.lanczos.exhausted {
        return Some(format!(
            "solver=lanczos; field=exhausted; first={}; replay={}",
            left.lanczos.exhausted, right.lanczos.exhausted
        ));
    }
    if let Some(mismatch) = first_pair_mismatch(&left.lanczos.pairs, &right.lanczos.pairs) {
        return Some(format!("solver=lanczos; {mismatch}"));
    }
    if left.lobpcg.iterations != right.lobpcg.iterations {
        return Some(format!(
            "solver=lobpcg; field=iterations; first={}; replay={}",
            left.lobpcg.iterations, right.lobpcg.iterations
        ));
    }
    first_pair_mismatch(&left.lobpcg.pairs, &right.lobpcg.pairs)
        .map(|mismatch| format!("solver=lobpcg; {mismatch}"))
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
    push_text(&mut bytes, "fs-la-version");
    push_text(&mut bytes, FS_LA_VERSION);
    push_text(&mut bytes, "fs-math-version");
    push_text(&mut bytes, FS_MATH_VERSION);
    bytes
}

fn green_inputs() -> Vec<u8> {
    let mut bytes = common_frame_prefix(b"bedrock:fs-la-eigensolver-replay:v1");
    push_text(&mut bytes, "operator");
    push_text(
        &mut bytes,
        "matrix-free-dirichlet-laplacian-1d:y[i]=2*x[i]-x[i-1]-x[i+1]:v1",
    );
    push_text(&mut bytes, "dimension-pairs-largest");
    push_len(&mut bytes, N);
    push_len(&mut bytes, PAIRS);
    push_u64(&mut bytes, 0);
    push_text(&mut bytes, "lanczos-prefix-pairs-and-split-steps");
    push_len(&mut bytes, LANCZOS_PREFIX_PAIRS);
    push_len(&mut bytes, LANCZOS_FIRST);
    push_len(&mut bytes, LANCZOS_SECOND);
    push_text(&mut bytes, "lobpcg-split-iterations");
    push_len(&mut bytes, LOBPCG_FIRST);
    push_len(&mut bytes, LOBPCG_SECOND);
    push_text(&mut bytes, "lobpcg-preconditioner");
    push_text(&mut bytes, "identity-copy:v1");
    push_text(&mut bytes, "checkpoint-policy");
    push_text(
        &mut bytes,
        "clone-state-at-first-cut;resume-second-cut;compare-uninterrupted-total:v1",
    );
    bytes
}

fn corruption_coordinates() -> (usize, u32) {
    let output_count = u64::try_from(PAIRS).expect("fixture pair count fits u64");
    let output = usize::try_from(RED_SEED % output_count).expect("derived output fits usize");
    let bit = u32::try_from((RED_SEED >> 16) % 52).expect("derived mantissa bit fits u32");
    (output, bit)
}

fn red_inputs(green: &[u8]) -> Vec<u8> {
    let (output, bit) = corruption_coordinates();
    let mut bytes = common_frame_prefix(b"bedrock:fs-la-eigensolver-red:v1");
    push_nested(&mut bytes, "nested-green-input-frame", green);
    push_text(&mut bytes, "corruption-seed-lanczos-value-mantissa-bit");
    push_u64(&mut bytes, RED_SEED);
    push_len(&mut bytes, output);
    push_u64(&mut bytes, u64::from(bit));
    push_text(
        &mut bytes,
        "policy=flip-one-derived-lanczos-eigenvalue-reference-mantissa-bit:v1",
    );
    bytes
}

fn green_outcome(input_frame: &[u8]) -> CaseOutcome {
    let inputs_hex = hex_bytes(input_frame);
    let first = match evaluate() {
        Ok(evaluation) => evaluation,
        Err(error) => {
            return CaseOutcome::fail(format!("{error}; inputs_hex={inputs_hex}"))
                .with_evidence("crates/fs-la/CONTRACT.md#determinism-class");
        }
    };
    let replay = match evaluate() {
        Ok(evaluation) => evaluation,
        Err(error) => {
            return CaseOutcome::fail(format!(
                "stage=same-run-replay; replay_error={error}; inputs_hex={inputs_hex}"
            ))
            .with_evidence("crates/fs-la/CONTRACT.md#determinism-class");
        }
    };
    if let Some(mismatch) = first_evaluation_mismatch(&first, &replay) {
        return CaseOutcome::fail(format!(
            "stage=same-run-replay; {mismatch}; first_digest={:016x}; replay_digest={:016x}; inputs_hex={inputs_hex}",
            evaluation_digest(&first),
            evaluation_digest(&replay)
        ))
        .with_evidence("crates/fs-la/CONTRACT.md#determinism-class");
    }
    CaseOutcome::pass(format!(
        "operator=laplacian-1d; n={N}; pairs={PAIRS}; lanczos_split={LANCZOS_FIRST}+{LANCZOS_SECOND}; lanczos_dimension={}; lanczos_exhausted={}; lobpcg_split={LOBPCG_FIRST}+{LOBPCG_SECOND}; lobpcg_iterations={}; output_digest={:016x}; same_run=identical",
        first.lanczos.dimension,
        first.lanczos.exhausted,
        first.lobpcg.iterations,
        evaluation_digest(&first)
    ))
    .with_evidence("crates/fs-la/CONTRACT.md#determinism-class")
}

fn red_outcome(input_frame: &[u8]) -> CaseOutcome {
    let inputs_hex = hex_bytes(input_frame);
    let first = match evaluate() {
        Ok(evaluation) => evaluation,
        Err(error) => {
            return CaseOutcome::fail(format!(
                "stage=red-prerequisite; error={error}; inputs_hex={inputs_hex}"
            ));
        }
    };
    let (output, bit) = corruption_coordinates();
    let actual = first.lanczos.pairs[output].value;
    let corrupted = actual ^ (1_u64 << bit);
    if actual == corrupted {
        return CaseOutcome::fail(format!(
            "stage=seeded-lanczos-reference-corruption; seed=0x{RED_SEED:016x}; output={output}; bit={bit}; error=derived-corruption-did-not-move-reference; inputs_hex={inputs_hex}"
        ));
    }
    CaseOutcome::fail(format!(
        "stage=seeded-lanczos-reference-corruption; seed=0x{RED_SEED:016x}; output={output}; bit={bit}; actual_bits=0x{actual:016x}; canonical_bits=0x{actual:016x}; corrupted_bits=0x{corrupted:016x}; inputs_hex={inputs_hex}"
    ))
    .with_evidence("crates/fs-la/tests/eigen_replay_casebook.rs#seeded-corruption")
}

#[test]
fn eigensolver_casebook_emits_replay_complete_green_record() {
    assert_eq!(CASEBOOK_RECORD_VERSION, 1);
    let inputs = green_inputs();
    let inputs_digest = fnv1a64(&inputs);
    assert_eq!(
        (inputs.len(), inputs_digest),
        (GREEN_FRAME_LEN, GREEN_FRAME_DIGEST)
    );
    let report = Suite::new(SUITE)
        .case(
            "lanczos-and-lobpcg-clone-checkpoint-bit-replay",
            inputs_digest,
            ToleranceSpec::Exact,
            move || green_outcome(&inputs),
        )
        .run();
    report.assert_green();
    let [record] = report.records.as_slice() else {
        panic!("the eigensolver replay suite must emit exactly one record");
    };
    assert_eq!(
        record.case,
        "lanczos-and-lobpcg-clone-checkpoint-bit-replay"
    );
    assert!(record.details.contains("lanczos_split=11+19"));
    assert!(record.details.contains("lobpcg_split=9+13"));
    assert!(record.details.contains("same_run=identical"));
}

#[test]
fn disclosed_seeded_lanczos_reference_corruption_turns_suite_red() {
    let green = green_inputs();
    let inputs = red_inputs(&green);
    let inputs_digest = fnv1a64(&inputs);
    let (output, bit) = corruption_coordinates();
    assert_eq!(
        (inputs.len(), inputs_digest),
        (RED_FRAME_LEN, RED_FRAME_DIGEST)
    );
    assert_eq!((output, bit), (2, 10));
    let make_report = || {
        let input_frame = inputs.clone();
        Suite::new(SUITE)
            .case(
                "seeded-lanczos-eigenvalue-reference-corruption",
                inputs_digest,
                ToleranceSpec::Exact,
                move || red_outcome(&input_frame),
            )
            .run()
    };
    let first = make_report();
    let replay = make_report();
    let first_failures = first.failures();
    let replay_failures = replay.failures();
    let [first_failure] = first_failures.as_slice() else {
        panic!("the disclosed corruption must produce exactly one failure");
    };
    let [replay_failure] = replay_failures.as_slice() else {
        panic!("the replayed corruption must produce exactly one failure");
    };
    assert_eq!(first_failure.json_line(), replay_failure.json_line());
    assert!(
        first_failure
            .details
            .contains("stage=seeded-lanczos-reference-corruption")
    );
    assert!(
        first_failure
            .details
            .contains(&format!("seed=0x{RED_SEED:016x}"))
    );
    assert!(first_failure.details.contains(&format!("output={output}")));
    assert!(first_failure.details.contains(&format!("bit={bit}")));
    assert!(first_failure.details.contains("inputs_hex="));
    assert!(
        first_failure
            .json_line()
            .contains("\"tolerance\":\"exact\",\"pass\":false")
    );
    let panic = catch_unwind(|| first.assert_green())
        .expect_err("the Casebook merge gate must reject the disclosed corruption");
    let message = panic
        .downcast_ref::<String>()
        .map(String::as_str)
        .or_else(|| panic.downcast_ref::<&str>().copied())
        .expect("Casebook panic carries text");
    assert!(message.contains("seeded-lanczos-eigenvalue-reference-corruption"));
    assert!(message.contains(&format!("seed=0x{RED_SEED:016x}")));
}
