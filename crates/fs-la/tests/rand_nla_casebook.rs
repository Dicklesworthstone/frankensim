//! Structured replay evidence for fs-la's keyed randomized-NLA surface.
//!
//! A dense positive-definite fixture is generated from logical fs-rand
//! coordinates under several simulated partitions and traversal orders. The
//! production range finder, randomized SVD, Hutchinson, and Hutch++ paths must
//! then replay every returned bit from the same seed, while a distinct seed
//! must independently move each API's stochastic snapshot.
//!
//! This is same-build, same-ISA G5 evidence for one disclosed finite fixture.
//! It is not fresh cross-ISA execution proof and makes no probability-coverage,
//! approximation-quality, performance, scheduling, or cancellation claim.

use core::fmt::Write as _;
use std::panic::{AssertUnwindSafe, catch_unwind};

use fs_casebook::{CASEBOOK_RECORD_VERSION, CaseOutcome, Suite, ToleranceSpec, fnv1a64};
use fs_la::VERSION as FS_LA_VERSION;
use fs_la::rand_nla::{RangeReport, TraceReport, hutch_pp, hutchinson, range_finder, rsvd};
use fs_rand::{
    STREAM_POSITION_IDENTITY_DOMAIN, STREAM_SEMANTICS_VERSION, Stream, StreamKey,
    VERSION as FS_RAND_VERSION,
};

const SUITE: &str = "bedrock/fs-la-randomized-nla-replay-v1";
const N: usize = 24;
const LATENT: usize = 6;
const RANK: usize = 5;
const OVERSAMPLE: usize = 3;
const Q_POWER: usize = 1;
const EFFECTIVE_RANK: usize = RANK + OVERSAMPLE;
const TRACE_PROBES: usize = 18;
const ROOT_SEED: u64 = 0x6A5E_5EED_4E4C_4101;
const INPUT_KERNEL: u32 = 0x4E4C_4149;
const ALGORITHM_SEED: u64 = 0x6A5E_A160_0000_0001;
const ALTERNATE_SEED: u64 = 0x6A5E_A160_0000_0002;
const RED_SEED: u64 = 0x6A5E_C0DE_4E4C_4101;
const LATENT_INPUT_DIGEST: u64 = 0x15b8_939e_267e_0b02;
const GREEN_FRAME_LEN: usize = 2_330;
const GREEN_FRAME_DIGEST: u64 = 0xe5e1_af26_678f_7f60;
const RED_FRAME_LEN: usize = 2_856;
const RED_FRAME_DIGEST: u64 = 0x12fb_5850_ef03_2d9f;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct GenerationPlan {
    partitions: usize,
    reverse_partitions: bool,
}

const GENERATION_PLANS: [GenerationPlan; 9] = [
    GenerationPlan {
        partitions: 1,
        reverse_partitions: false,
    },
    GenerationPlan {
        partitions: 2,
        reverse_partitions: false,
    },
    GenerationPlan {
        partitions: 2,
        reverse_partitions: true,
    },
    GenerationPlan {
        partitions: 3,
        reverse_partitions: false,
    },
    GenerationPlan {
        partitions: 3,
        reverse_partitions: true,
    },
    GenerationPlan {
        partitions: 5,
        reverse_partitions: false,
    },
    GenerationPlan {
        partitions: 5,
        reverse_partitions: true,
    },
    GenerationPlan {
        partitions: 8,
        reverse_partitions: false,
    },
    GenerationPlan {
        partitions: 8,
        reverse_partitions: true,
    },
];

#[derive(Debug, Clone, PartialEq, Eq)]
struct FixtureBits {
    latent: Vec<u64>,
    matrix: Vec<u64>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct RangeReportBits {
    rank: usize,
    estimate: u64,
    probes: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct RangeBits {
    basis: Vec<u64>,
    report: RangeReportBits,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct TraceBits {
    estimate: u64,
    probes: usize,
    variance: u64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct Snapshot {
    range: RangeBits,
    rsvd_u: Vec<u64>,
    rsvd_sigma: Vec<u64>,
    rsvd_v: Vec<u64>,
    rsvd_report: RangeReportBits,
    hutchinson: TraceBits,
    hutch_pp: TraceBits,
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

fn push_bits(bytes: &mut Vec<u8>, values: &[u64]) {
    push_len(bytes, values.len());
    for &value in values {
        push_u64(bytes, value);
    }
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

fn output_bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|value| value.to_bits()).collect()
}

fn digest_bits(values: &[u64]) -> u64 {
    let mut bytes = Vec::with_capacity(values.len().saturating_mul(8));
    for &value in values {
        bytes.extend_from_slice(&value.to_le_bytes());
    }
    fnv1a64(&bytes)
}

fn f64_values(bits: &[u64]) -> Vec<f64> {
    bits.iter().map(|&value| f64::from_bits(value)).collect()
}

fn generated_value_bits(row: usize, column: usize) -> u64 {
    let key = StreamKey {
        seed: ROOT_SEED,
        kernel: INPUT_KERNEL,
        tile: u32::try_from(row).expect("fixture row fits the logical tile field"),
    };
    let block = Stream::at(
        key,
        u64::try_from(column).expect("fixture column fits the stream index"),
    );
    let word = (u64::from(block[1]) << 32) | u64::from(block[0]);
    let top53 = word >> 11;
    // Exact dyadic ladder in [-1, 1): conversion and centering are exact,
    // followed by a power-of-two scale.
    ((top53 as f64 - 4_503_599_627_370_496.0) * 2.220_446_049_250_313e-16).to_bits()
}

fn materialize_latent(plan: GenerationPlan) -> Vec<u64> {
    assert!(
        plan.partitions > 0,
        "generation partitions must be positive"
    );
    let mut values = vec![0; N * LATENT];
    let mut fill_partition = |partition: usize| {
        for flat in (partition..values.len()).step_by(plan.partitions) {
            values[flat] = generated_value_bits(flat / LATENT, flat % LATENT);
        }
    };
    if plan.reverse_partitions {
        for partition in (0..plan.partitions).rev() {
            fill_partition(partition);
        }
    } else {
        for partition in 0..plan.partitions {
            fill_partition(partition);
        }
    }
    values
}

fn positive_definite_matrix(latent_bits: &[u64]) -> Vec<u64> {
    let latent = f64_values(latent_bits);
    let mut matrix = vec![0_u64; N * N];
    for row in 0..N {
        for column in row..N {
            let mut value = 0.0_f64;
            for component in 0..LATENT {
                value = latent[row * LATENT + component]
                    .mul_add(latent[column * LATENT + component], value);
            }
            if row == column {
                value += 0.25 + (row + 1) as f64 * 0.000_244_140_625;
            }
            let bits = value.to_bits();
            matrix[row * N + column] = bits;
            matrix[column * N + row] = bits;
        }
    }
    matrix
}

fn materialize_fixture(plan: GenerationPlan) -> FixtureBits {
    let latent = materialize_latent(plan);
    let matrix = positive_definite_matrix(&latent);
    FixtureBits { latent, matrix }
}

fn first_mismatch(left: &[u64], right: &[u64]) -> Option<(usize, u64, u64)> {
    left.iter()
        .zip(right)
        .enumerate()
        .find_map(|(index, (&lhs, &rhs))| (lhs != rhs).then_some((index, lhs, rhs)))
        .or_else(|| (left.len() != right.len()).then_some((left.len().min(right.len()), 0, 0)))
}

fn validate_generation() -> Result<FixtureBits, String> {
    let canonical = materialize_fixture(GENERATION_PLANS[0]);
    for plan in GENERATION_PLANS {
        let candidate = materialize_fixture(plan);
        for (stage, expected, computed) in [
            (
                "logical-latent-generation",
                canonical.latent.as_slice(),
                candidate.latent.as_slice(),
            ),
            (
                "positive-definite-materialization",
                canonical.matrix.as_slice(),
                candidate.matrix.as_slice(),
            ),
        ] {
            if let Some((index, reference, actual)) = first_mismatch(expected, computed) {
                return Err(format!(
                    "stage={stage}; partitions={}; reverse={}; index={index}; canonical_bits=0x{reference:016x}; computed_bits=0x{actual:016x}",
                    plan.partitions, plan.reverse_partitions,
                ));
            }
        }
    }
    Ok(canonical)
}

fn range_bits(basis: &[f64], report: &RangeReport) -> RangeBits {
    RangeBits {
        basis: output_bits(basis),
        report: range_report_bits(report),
    }
}

fn range_report_bits(report: &RangeReport) -> RangeReportBits {
    RangeReportBits {
        rank: report.rank,
        estimate: report.est_error.to_bits(),
        probes: report.probes,
    }
}

fn trace_bits(report: &TraceReport) -> TraceBits {
    TraceBits {
        estimate: report.estimate.to_bits(),
        probes: report.probes,
        variance: report.variance_est.to_bits(),
    }
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

fn evaluate(matrix_bits: &[u64], seed: u64) -> Result<Snapshot, String> {
    let matrix = f64_values(matrix_bits);
    catch_unwind(AssertUnwindSafe(|| {
        let (range_basis, range_report) =
            range_finder(&matrix, N, N, RANK, OVERSAMPLE, Q_POWER, seed);
        let (u, sigma, v, rsvd_report) = rsvd(&matrix, N, N, RANK, OVERSAMPLE, Q_POWER, seed);
        let hutchinson_report = hutchinson(&matrix, N, TRACE_PROBES, seed);
        let hutch_pp_report = hutch_pp(&matrix, N, TRACE_PROBES, seed);
        Snapshot {
            range: range_bits(&range_basis, &range_report),
            rsvd_u: output_bits(&u),
            rsvd_sigma: output_bits(&sigma),
            rsvd_v: output_bits(&v),
            rsvd_report: range_report_bits(&rsvd_report),
            hutchinson: trace_bits(&hutchinson_report),
            hutch_pp: trace_bits(&hutch_pp_report),
        }
    }))
    .map_err(|payload| {
        format!(
            "stage=randomized-nla-execution; seed=0x{seed:016x}; panic={}",
            panic_message(&*payload)
        )
    })
}

fn finite_bits(bits: u64) -> bool {
    f64::from_bits(bits).is_finite()
}

fn validate_snapshot(snapshot: &Snapshot) -> Result<(), String> {
    for (field, computed, expected) in [
        (
            "range_basis_len",
            snapshot.range.basis.len(),
            N * EFFECTIVE_RANK,
        ),
        ("range_rank", snapshot.range.report.rank, EFFECTIVE_RANK),
        ("range_probes", snapshot.range.report.probes, 8),
        ("rsvd_u_len", snapshot.rsvd_u.len(), N * RANK),
        ("rsvd_sigma_len", snapshot.rsvd_sigma.len(), RANK),
        ("rsvd_v_len", snapshot.rsvd_v.len(), N * RANK),
        (
            "rsvd_report_rank",
            snapshot.rsvd_report.rank,
            EFFECTIVE_RANK,
        ),
        ("rsvd_report_probes", snapshot.rsvd_report.probes, 8),
        (
            "hutchinson_probes",
            snapshot.hutchinson.probes,
            TRACE_PROBES,
        ),
        ("hutch_pp_probes", snapshot.hutch_pp.probes, TRACE_PROBES),
    ] {
        if computed != expected {
            return Err(format!(
                "stage=structural-evidence; field={field}; computed={computed}; expected={expected}"
            ));
        }
    }

    for (field, values) in [
        ("range_basis", snapshot.range.basis.as_slice()),
        ("rsvd_u", snapshot.rsvd_u.as_slice()),
        ("rsvd_sigma", snapshot.rsvd_sigma.as_slice()),
        ("rsvd_v", snapshot.rsvd_v.as_slice()),
    ] {
        if let Some((index, &bits)) = values
            .iter()
            .enumerate()
            .find(|(_, bits)| !finite_bits(**bits))
        {
            return Err(format!(
                "stage=finite-evidence; field={field}; index={index}; bits=0x{bits:016x}"
            ));
        }
    }
    for (field, bits) in [
        ("range_estimate", snapshot.range.report.estimate),
        ("rsvd_estimate", snapshot.rsvd_report.estimate),
        ("hutchinson_estimate", snapshot.hutchinson.estimate),
        ("hutchinson_variance", snapshot.hutchinson.variance),
        ("hutch_pp_estimate", snapshot.hutch_pp.estimate),
        ("hutch_pp_variance", snapshot.hutch_pp.variance),
    ] {
        if !finite_bits(bits) {
            return Err(format!(
                "stage=finite-evidence; field={field}; bits=0x{bits:016x}"
            ));
        }
    }
    Ok(())
}

fn snapshot_frame(snapshot: &Snapshot) -> Vec<u8> {
    let mut bytes = b"bedrock:fs-la-randomized-nla-snapshot:v1".to_vec();
    push_text(&mut bytes, "range-basis");
    push_bits(&mut bytes, &snapshot.range.basis);
    push_text(&mut bytes, "range-report-rank-estimate-probes");
    push_len(&mut bytes, snapshot.range.report.rank);
    push_u64(&mut bytes, snapshot.range.report.estimate);
    push_len(&mut bytes, snapshot.range.report.probes);
    push_text(&mut bytes, "rsvd-u-sigma-v");
    push_bits(&mut bytes, &snapshot.rsvd_u);
    push_bits(&mut bytes, &snapshot.rsvd_sigma);
    push_bits(&mut bytes, &snapshot.rsvd_v);
    push_text(&mut bytes, "rsvd-report-rank-estimate-probes");
    push_len(&mut bytes, snapshot.rsvd_report.rank);
    push_u64(&mut bytes, snapshot.rsvd_report.estimate);
    push_len(&mut bytes, snapshot.rsvd_report.probes);
    for (label, report) in [
        ("hutchinson", &snapshot.hutchinson),
        ("hutch-pp", &snapshot.hutch_pp),
    ] {
        push_text(&mut bytes, label);
        push_u64(&mut bytes, report.estimate);
        push_len(&mut bytes, report.probes);
        push_u64(&mut bytes, report.variance);
    }
    bytes
}

fn snapshot_digest(snapshot: &Snapshot) -> u64 {
    fnv1a64(&snapshot_frame(snapshot))
}

fn rsvd_equal(left: &Snapshot, right: &Snapshot) -> bool {
    left.rsvd_u == right.rsvd_u
        && left.rsvd_sigma == right.rsvd_sigma
        && left.rsvd_v == right.rsvd_v
        && left.rsvd_report == right.rsvd_report
}

fn first_api_mismatch(left: &Snapshot, right: &Snapshot) -> Option<&'static str> {
    if left.range != right.range {
        Some("range_finder")
    } else if !rsvd_equal(left, right) {
        Some("rsvd")
    } else if left.hutchinson != right.hutchinson {
        Some("hutchinson")
    } else if left.hutch_pp != right.hutch_pp {
        Some("hutch_pp")
    } else {
        None
    }
}

fn first_seed_insensitive_api(left: &Snapshot, right: &Snapshot) -> Option<&'static str> {
    if left.range == right.range {
        Some("range_finder")
    } else if rsvd_equal(left, right) {
        Some("rsvd")
    } else if left.hutchinson == right.hutchinson {
        Some("hutchinson")
    } else if left.hutch_pp == right.hutch_pp {
        Some("hutch_pp")
    } else {
        None
    }
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
    push_text(&mut bytes, "fs-rand-version");
    push_text(&mut bytes, FS_RAND_VERSION);
    push_text(&mut bytes, "stream-semantics-version");
    push_u64(&mut bytes, u64::from(STREAM_SEMANTICS_VERSION));
    push_text(&mut bytes, "stream-position-identity-domain");
    push_text(&mut bytes, STREAM_POSITION_IDENTITY_DOMAIN);
    bytes
}

fn green_inputs(fixture: &FixtureBits) -> Vec<u8> {
    let mut bytes = common_frame_prefix(b"bedrock:fs-la-randomized-nla-replay:v1");
    push_text(&mut bytes, "logical-coordinate-policy");
    push_text(
        &mut bytes,
        "kernel=fixed;tile=row-u32;index=column-u64;Stream::at;word=(block[1]<<32)|block[0]",
    );
    push_text(&mut bytes, "finite-f64-policy");
    push_text(
        &mut bytes,
        "top53=word>>11;value=(f64(top53)-2^52)*2^-52;exact-dyadic:[-1,1):v1",
    );
    push_text(&mut bytes, "positive-definite-policy");
    push_text(
        &mut bytes,
        "A=B*B^T+diag(0.25+(row+1)*2^-12);upper-triangle-fixed-component-order-mul_add;mirror-bits:v1",
    );
    push_text(&mut bytes, "root-seed-input-kernel");
    push_u64(&mut bytes, ROOT_SEED);
    push_u64(&mut bytes, u64::from(INPUT_KERNEL));
    push_text(&mut bytes, "algorithm-seeds");
    push_u64(&mut bytes, ALGORITHM_SEED);
    push_u64(&mut bytes, ALTERNATE_SEED);
    push_text(&mut bytes, "dimensions-and-algorithm-parameters");
    for value in [N, LATENT, RANK, OVERSAMPLE, Q_POWER, TRACE_PROBES] {
        push_len(&mut bytes, value);
    }
    push_text(&mut bytes, "algorithms");
    push_text(
        &mut bytes,
        "range_finder,rsvd,hutchinson,hutch_pp;full-returned-bit-snapshot",
    );
    push_text(&mut bytes, "generation-plans");
    push_len(&mut bytes, GENERATION_PLANS.len());
    for plan in GENERATION_PLANS {
        push_len(&mut bytes, plan.partitions);
        push_u64(&mut bytes, u64::from(u8::from(plan.reverse_partitions)));
    }
    push_text(&mut bytes, "latent-input-bits");
    push_bits(&mut bytes, &fixture.latent);
    bytes
}

fn corruption_coordinates() -> (usize, u32) {
    let output_count = u64::try_from(RANK).expect("fixture rank fits u64");
    let output = usize::try_from(RED_SEED % output_count).expect("derived output fits usize");
    let bit = u32::try_from((RED_SEED >> 16) % 52).expect("derived mantissa bit fits u32");
    (output, bit)
}

fn red_inputs(green: &[u8]) -> Vec<u8> {
    let (output, bit) = corruption_coordinates();
    let mut bytes = common_frame_prefix(b"bedrock:fs-la-randomized-nla-red:v1");
    push_nested(&mut bytes, "nested-green-input-frame", green);
    push_text(&mut bytes, "corruption-seed-sigma-output-mantissa-bit");
    push_u64(&mut bytes, RED_SEED);
    push_len(&mut bytes, output);
    push_u64(&mut bytes, u64::from(bit));
    push_text(
        &mut bytes,
        "policy=flip-one-derived-rsvd-sigma-reference-mantissa-bit:v1",
    );
    bytes
}

fn green_outcome(input_frame: &[u8]) -> CaseOutcome {
    let inputs_hex = hex_bytes(input_frame);
    let fixture = match validate_generation() {
        Ok(fixture) => fixture,
        Err(error) => {
            return CaseOutcome::fail(format!("{error}; inputs_hex={inputs_hex}"))
                .with_evidence("crates/fs-rand/CONTRACT.md#determinism-class");
        }
    };
    let first = match evaluate(&fixture.matrix, ALGORITHM_SEED) {
        Ok(snapshot) => snapshot,
        Err(error) => {
            return CaseOutcome::fail(format!("{error}; inputs_hex={inputs_hex}"))
                .with_evidence("crates/fs-la/CONTRACT.md#determinism-class");
        }
    };
    if let Err(error) = validate_snapshot(&first) {
        return CaseOutcome::fail(format!("{error}; inputs_hex={inputs_hex}"))
            .with_evidence("crates/fs-la/CONTRACT.md#public-types-and-semantics");
    }
    let replay = match evaluate(&fixture.matrix, ALGORITHM_SEED) {
        Ok(snapshot) => snapshot,
        Err(error) => {
            return CaseOutcome::fail(format!(
                "stage=same-seed-replay; replay_error={error}; inputs_hex={inputs_hex}"
            ))
            .with_evidence("crates/fs-la/CONTRACT.md#determinism-class");
        }
    };
    if let Some(api) = first_api_mismatch(&first, &replay) {
        return CaseOutcome::fail(format!(
            "stage=same-seed-replay; api={api}; first_digest={:016x}; replay_digest={:016x}; inputs_hex={inputs_hex}",
            snapshot_digest(&first),
            snapshot_digest(&replay),
        ))
        .with_evidence("crates/fs-la/CONTRACT.md#determinism-class");
    }
    let alternate = match evaluate(&fixture.matrix, ALTERNATE_SEED) {
        Ok(snapshot) => snapshot,
        Err(error) => {
            return CaseOutcome::fail(format!(
                "stage=alternate-seed; error={error}; inputs_hex={inputs_hex}"
            ))
            .with_evidence("crates/fs-la/CONTRACT.md#determinism-class");
        }
    };
    if let Err(error) = validate_snapshot(&alternate) {
        return CaseOutcome::fail(format!(
            "stage=alternate-seed; {error}; inputs_hex={inputs_hex}"
        ))
        .with_evidence("crates/fs-la/CONTRACT.md#public-types-and-semantics");
    }
    if let Some(api) = first_seed_insensitive_api(&first, &alternate) {
        return CaseOutcome::fail(format!(
            "stage=seed-separation; api={api}; canonical_seed=0x{ALGORITHM_SEED:016x}; alternate_seed=0x{ALTERNATE_SEED:016x}; canonical_digest={:016x}; alternate_digest={:016x}; inputs_hex={inputs_hex}",
            snapshot_digest(&first),
            snapshot_digest(&alternate),
        ))
        .with_evidence("crates/fs-la/src/rand_nla.rs");
    }

    CaseOutcome::pass(format!(
        "shape={N}x{N}; latent={LATENT}; rank={RANK}; oversample={OVERSAMPLE}; q_power={Q_POWER}; trace_probes={TRACE_PROBES}; generation_plans={}; matrix_digest={:016x}; snapshot_digest={:016x}; alternate_digest={:016x}; same_seed=identical; alternate_seed=distinct",
        GENERATION_PLANS.len(),
        digest_bits(&fixture.matrix),
        snapshot_digest(&first),
        snapshot_digest(&alternate),
    ))
    .with_evidence("crates/fs-la/CONTRACT.md#determinism-class")
    .with_evidence("crates/fs-rand/CONTRACT.md#determinism-class")
}

fn red_outcome(input_frame: &[u8]) -> CaseOutcome {
    let inputs_hex = hex_bytes(input_frame);
    let fixture = match validate_generation() {
        Ok(fixture) => fixture,
        Err(error) => {
            return CaseOutcome::fail(format!(
                "stage=red-generation-prerequisite; error={error}; inputs_hex={inputs_hex}"
            ));
        }
    };
    let first = match evaluate(&fixture.matrix, ALGORITHM_SEED) {
        Ok(snapshot) => snapshot,
        Err(error) => {
            return CaseOutcome::fail(format!(
                "stage=red-execution-prerequisite; error={error}; inputs_hex={inputs_hex}"
            ));
        }
    };
    let replay = match evaluate(&fixture.matrix, ALGORITHM_SEED) {
        Ok(snapshot) => snapshot,
        Err(error) => {
            return CaseOutcome::fail(format!(
                "stage=red-replay-prerequisite; error={error}; inputs_hex={inputs_hex}"
            ));
        }
    };
    if let Some(api) = first_api_mismatch(&first, &replay) {
        return CaseOutcome::fail(format!(
            "stage=red-same-seed-replay; api={api}; first_digest={:016x}; replay_digest={:016x}; inputs_hex={inputs_hex}",
            snapshot_digest(&first),
            snapshot_digest(&replay),
        ));
    }

    let (output, bit) = corruption_coordinates();
    let actual = first.rsvd_sigma[output];
    let corrupted = actual ^ (1_u64 << bit);
    if actual == corrupted {
        return CaseOutcome::fail(format!(
            "stage=seeded-rsvd-sigma-reference-corruption; seed=0x{RED_SEED:016x}; output={output}; bit={bit}; error=derived-corruption-did-not-move-reference; inputs_hex={inputs_hex}"
        ));
    }
    CaseOutcome::fail(format!(
        "stage=seeded-rsvd-sigma-reference-corruption; seed=0x{RED_SEED:016x}; output={output}; bit={bit}; actual_bits=0x{actual:016x}; canonical_bits=0x{actual:016x}; corrupted_bits=0x{corrupted:016x}; inputs_hex={inputs_hex}"
    ))
    .with_evidence("crates/fs-la/tests/rand_nla_casebook.rs#seeded-corruption")
}

#[test]
fn randomized_nla_casebook_emits_replay_complete_green_record() {
    assert_eq!(CASEBOOK_RECORD_VERSION, 1);
    assert_eq!(STREAM_SEMANTICS_VERSION, 1);
    let fixture = materialize_fixture(GENERATION_PLANS[0]);
    assert_eq!(digest_bits(&fixture.latent), LATENT_INPUT_DIGEST);
    let inputs = green_inputs(&fixture);
    let inputs_digest = fnv1a64(&inputs);
    assert_eq!(
        (inputs.len(), inputs_digest),
        (GREEN_FRAME_LEN, GREEN_FRAME_DIGEST)
    );

    let report = Suite::new(SUITE)
        .case(
            "keyed-randomized-nla-full-snapshot-replay",
            inputs_digest,
            ToleranceSpec::Exact,
            move || green_outcome(&inputs),
        )
        .run();
    report.assert_green();
    let [record] = report.records.as_slice() else {
        panic!("the randomized-NLA suite must emit exactly one record");
    };
    assert_eq!(record.case, "keyed-randomized-nla-full-snapshot-replay");
    assert!(record.details.contains("generation_plans=9"));
    assert!(record.details.contains("same_seed=identical"));
    assert!(record.details.contains("alternate_seed=distinct"));
}

#[test]
fn disclosed_seeded_rsvd_reference_corruption_turns_suite_red() {
    let fixture = materialize_fixture(GENERATION_PLANS[0]);
    let green = green_inputs(&fixture);
    let inputs = red_inputs(&green);
    let inputs_digest = fnv1a64(&inputs);
    let (output, bit) = corruption_coordinates();
    assert_eq!(
        (inputs.len(), inputs_digest),
        (RED_FRAME_LEN, RED_FRAME_DIGEST)
    );
    assert_eq!((output, bit), (4, 44));

    let make_report = || {
        let input_frame = inputs.clone();
        Suite::new(SUITE)
            .case(
                "seeded-rsvd-sigma-reference-corruption",
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
            .contains("stage=seeded-rsvd-sigma-reference-corruption")
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
    assert!(message.contains("seeded-rsvd-sigma-reference-corruption"));
    assert!(message.contains(&format!("seed=0x{RED_SEED:016x}")));
}
