//! Structured QMC replay evidence for the BEDROCK conformance parent
//! (6ys.18.18).
//!
//! This complements the Philox/checkpoint Casebook with production Sobol,
//! Owen-scrambling, and CBC-lattice paths.  The finite fixtures bind exact
//! inputs and returned bits; they do not certify all dimensions, seeds,
//! integrands, lattice sizes, ISAs, or performance regimes.

use fs_casebook::{
    CASEBOOK_RECORD_VERSION, CaseOutcome, Suite, SuiteReport, ToleranceSpec, fnv1a64,
};
use fs_rand::qmc::{Lattice, Sobol, baker};
use fs_rand::{STREAM_POSITION_IDENTITY_DOMAIN, STREAM_SEMANTICS_VERSION};
use std::fmt::Write as _;
use std::panic::catch_unwind;

const SUITE: &str = "fs-rand/qmc-replay-v1";
const INPUT_FRAME_VERSION: u32 = 1;
const INPUT_DOMAIN: &[u8] = b"fs-rand:qmc-casebook-input:v1";

const SOBOL_DIM: usize = 5;
const SOBOL_POINTS: u32 = 64;
const SOBOL_KAT: [f64; 8] = [0.0, 0.5, 0.75, 0.25, 0.375, 0.875, 0.625, 0.125];
const SOBOL_INPUTS_DIGEST: u64 = 0xE19F_4E88_D968_70EA;

const OWEN_DIM: usize = 4;
const OWEN_POINTS: u32 = 128;
const OWEN_SEED: u64 = 0x0A11_CE00_0000_0001;
const OWEN_ALT_SEED: u64 = 0x0A11_CE00_0000_0002;
const OWEN_STREAM_KERNEL: u32 = 0x0E11;
const OWEN_INPUTS_DIGEST: u64 = 0xD52F_0C16_0B18_9289;

const LATTICE_N: u32 = 257;
const LATTICE_DIM: usize = 5;
const BAKER_INPUTS: [f64; 5] = [0.0, 0.125, 0.5, 0.875, 1.0];
const BAKER_OUTPUTS: [f64; 5] = [0.0, 0.25, 1.0, 0.25, 0.0];
const LATTICE_INPUTS_DIGEST: u64 = 0x08CF_F502_6B21_3C41;

const CORRUPTION_SEED: u64 = 0x514D_4300_0000_0A05;
const CORRUPTION_INPUTS_DIGEST: u64 = 0x80D9_75B5_2D57_9C57;

fn usize_u64(value: usize) -> u64 {
    u64::try_from(value).expect("fixture cardinality fits u64")
}

fn push_text(bytes: &mut Vec<u8>, value: &str) {
    bytes.extend_from_slice(&usize_u64(value.len()).to_le_bytes());
    bytes.extend_from_slice(value.as_bytes());
}

fn push_f64(bytes: &mut Vec<u8>, value: f64) {
    bytes.extend_from_slice(&value.to_bits().to_le_bytes());
}

fn input_header(case: &str) -> Vec<u8> {
    let mut bytes = INPUT_DOMAIN.to_vec();
    push_text(&mut bytes, case);
    bytes.extend_from_slice(&INPUT_FRAME_VERSION.to_le_bytes());
    push_text(&mut bytes, fs_rand::VERSION);
    bytes.extend_from_slice(&STREAM_SEMANTICS_VERSION.to_le_bytes());
    push_text(&mut bytes, STREAM_POSITION_IDENTITY_DOMAIN);
    bytes
}

fn sobol_inputs() -> Vec<u8> {
    let mut bytes = input_header("sobol-gray-random-access");
    bytes.extend_from_slice(&usize_u64(SOBOL_DIM).to_le_bytes());
    bytes.extend_from_slice(&SOBOL_POINTS.to_le_bytes());
    bytes.extend_from_slice(&usize_u64(SOBOL_KAT.len()).to_le_bytes());
    for value in SOBOL_KAT {
        push_f64(&mut bytes, value);
    }
    bytes
}

fn owen_inputs() -> Vec<u8> {
    let mut bytes = input_header("owen-seeded-replay");
    bytes.extend_from_slice(&usize_u64(OWEN_DIM).to_le_bytes());
    bytes.extend_from_slice(&OWEN_POINTS.to_le_bytes());
    bytes.extend_from_slice(&OWEN_SEED.to_le_bytes());
    bytes.extend_from_slice(&OWEN_ALT_SEED.to_le_bytes());
    bytes.extend_from_slice(&OWEN_STREAM_KERNEL.to_le_bytes());
    push_text(&mut bytes, "nested-uniform-prefix-philox");
    push_text(&mut bytes, "exact-per-dimension-m7-stratification");
    bytes
}

fn lattice_inputs() -> Vec<u8> {
    let mut bytes = input_header("cbc-lattice-replay");
    bytes.extend_from_slice(&LATTICE_N.to_le_bytes());
    bytes.extend_from_slice(&usize_u64(LATTICE_DIM).to_le_bytes());
    push_text(&mut bytes, "cbc-b2-gamma1");
    bytes.extend_from_slice(&LATTICE_N.to_le_bytes());
    bytes.extend_from_slice(&usize_u64(BAKER_INPUTS.len()).to_le_bytes());
    for (&input, &output) in BAKER_INPUTS.iter().zip(&BAKER_OUTPUTS) {
        push_f64(&mut bytes, input);
        push_f64(&mut bytes, output);
    }
    bytes
}

fn output_digest(domain: &[u8], values: &[f64]) -> u64 {
    let mut bytes = domain.to_vec();
    bytes.extend_from_slice(&usize_u64(values.len()).to_le_bytes());
    for &value in values {
        push_f64(&mut bytes, value);
    }
    fnv1a64(&bytes)
}

fn owen_output_digest(primary: &[f64], alternate: &[f64]) -> u64 {
    let mut bytes = b"fs-rand:qmc:owen-output:v1".to_vec();
    bytes.extend_from_slice(&usize_u64(primary.len()).to_le_bytes());
    for &value in primary {
        push_f64(&mut bytes, value);
    }
    bytes.extend_from_slice(&usize_u64(alternate.len()).to_le_bytes());
    for &value in alternate {
        push_f64(&mut bytes, value);
    }
    fnv1a64(&bytes)
}

fn first_bits_mismatch(left: &[f64], right: &[f64], dim: usize) -> Option<String> {
    if left.len() != right.len() {
        return Some(format!("length:{}!={}", left.len(), right.len()));
    }
    for (index, (&a, &b)) in left.iter().zip(right).enumerate() {
        if a.to_bits() != b.to_bits() {
            return Some(format!(
                "point={}; coordinate={}; left=0x{:016x}; right=0x{:016x}",
                index / dim,
                index % dim,
                a.to_bits(),
                b.to_bits()
            ));
        }
    }
    None
}

fn qmc_failure(detail: String) -> CaseOutcome {
    CaseOutcome::fail(detail).with_evidence("crates/fs-rand/CONTRACT.md#qmc-qmc-module")
}

fn sobol_outcome() -> CaseOutcome {
    let kat = Sobol::new(1);
    let mut scalar = [0.0];
    for (point, expected) in SOBOL_KAT.into_iter().enumerate() {
        kat.point(u32::try_from(point).expect("small KAT index"), &mut scalar);
        if scalar[0].to_bits() != expected.to_bits() {
            return qmc_failure(format!(
                "stage=gray-code-kat; point={point}; computed=0x{:016x}; reference=0x{:016x}",
                scalar[0].to_bits(),
                expected.to_bits()
            ));
        }
    }

    let sobol = Sobol::new(SOBOL_DIM);
    let materialized = sobol.points(SOBOL_POINTS);
    let replay = Sobol::new(SOBOL_DIM).points(SOBOL_POINTS);
    let expected_len = usize::try_from(SOBOL_POINTS).expect("point count fits usize") * SOBOL_DIM;
    if materialized.len() != expected_len || replay.len() != expected_len {
        return qmc_failure(format!(
            "stage=frame-shape; materialized_len={}; replay_len={}; expected_len={expected_len}",
            materialized.len(),
            replay.len()
        ));
    }
    if let Some(mismatch) = first_bits_mismatch(&materialized, &replay, SOBOL_DIM) {
        return qmc_failure(format!("stage=independent-replay; {mismatch}"));
    }
    let mut point = vec![0.0; SOBOL_DIM];
    for index in 0..SOBOL_POINTS {
        sobol.point(index, &mut point);
        let start = usize::try_from(index).expect("point index fits usize") * SOBOL_DIM;
        if let Some(mismatch) =
            first_bits_mismatch(&materialized[start..start + SOBOL_DIM], &point, SOBOL_DIM)
        {
            return qmc_failure(format!(
                "stage=random-access-materialization; outer-point={index}; {mismatch}"
            ));
        }
    }

    let count = usize::try_from(SOBOL_POINTS).expect("point count fits usize");
    let mut bins = vec![vec![0_u32; count]; SOBOL_DIM];
    for row in materialized.chunks_exact(SOBOL_DIM) {
        for (dimension, &value) in row.iter().enumerate() {
            if !(0.0..1.0).contains(&value) {
                return qmc_failure(format!(
                    "stage=bounds; dimension={dimension}; value=0x{:016x}",
                    value.to_bits()
                ));
            }
            let bin = (value * f64::from(SOBOL_POINTS)) as usize;
            bins[dimension][bin] += 1;
        }
    }
    if let Some((dimension, bin, count)) = bins.iter().enumerate().find_map(|(dimension, row)| {
        row.iter()
            .enumerate()
            .find(|(_, count)| **count != 1)
            .map(|(bin, &count)| (dimension, bin, count))
    }) {
        return qmc_failure(format!(
            "stage=stratification; dimension={dimension}; bin={bin}; count={count}; expected=1"
        ));
    }

    CaseOutcome::pass(format!(
        "frame_version={INPUT_FRAME_VERSION}; fs_rand={}; stream_semantics={STREAM_SEMANTICS_VERSION}; dim={SOBOL_DIM}; points={SOBOL_POINTS}; gray_kat=8/8; random_access=exact; independent_replay=exact; stratification=m6; output_digest={:016x}",
        fs_rand::VERSION,
        output_digest(b"fs-rand:qmc:sobol-output:v1", &materialized)
    ))
    .with_evidence("crates/fs-rand/CONTRACT.md#qmc-qmc-module")
}

fn owen_outcome() -> CaseOutcome {
    let first = Sobol::scrambled(OWEN_DIM, OWEN_SEED).points(OWEN_POINTS);
    let replay = Sobol::scrambled(OWEN_DIM, OWEN_SEED).points(OWEN_POINTS);
    let alternate = Sobol::scrambled(OWEN_DIM, OWEN_ALT_SEED).points(OWEN_POINTS);
    let alternate_replay = Sobol::scrambled(OWEN_DIM, OWEN_ALT_SEED).points(OWEN_POINTS);
    let expected_len = usize::try_from(OWEN_POINTS).expect("point count fits usize") * OWEN_DIM;
    if first.len() != expected_len
        || replay.len() != expected_len
        || alternate.len() != expected_len
        || alternate_replay.len() != expected_len
    {
        return qmc_failure(format!(
            "stage=frame-shape; primary_len={}; primary_replay_len={}; alternate_len={}; alternate_replay_len={}; expected_len={expected_len}",
            first.len(),
            replay.len(),
            alternate.len(),
            alternate_replay.len()
        ));
    }
    if let Some(mismatch) = first_bits_mismatch(&first, &replay, OWEN_DIM) {
        return qmc_failure(format!("stage=same-seed-replay; {mismatch}"));
    }
    if let Some(mismatch) = first_bits_mismatch(&alternate, &alternate_replay, OWEN_DIM) {
        return qmc_failure(format!("stage=alternate-seed-replay; {mismatch}"));
    }
    if first_bits_mismatch(&first, &alternate, OWEN_DIM).is_none() {
        return qmc_failure(format!(
            "stage=seed-separation; input_seed=0x{OWEN_SEED:016x}; alternate_seed=0x{OWEN_ALT_SEED:016x}; outputs=identical"
        ));
    }

    let count = usize::try_from(OWEN_POINTS).expect("point count fits usize");
    let mut bins = vec![vec![0_u32; count]; OWEN_DIM];
    for row in first.chunks_exact(OWEN_DIM) {
        for (dimension, &value) in row.iter().enumerate() {
            if !(0.0..1.0).contains(&value) {
                return qmc_failure(format!(
                    "stage=bounds; dimension={dimension}; value=0x{:016x}",
                    value.to_bits()
                ));
            }
            let bin = (value * f64::from(OWEN_POINTS)) as usize;
            bins[dimension][bin] += 1;
        }
    }
    if let Some((dimension, bin, count)) = bins.iter().enumerate().find_map(|(dimension, row)| {
        row.iter()
            .enumerate()
            .find(|(_, count)| **count != 1)
            .map(|(bin, &count)| (dimension, bin, count))
    }) {
        return qmc_failure(format!(
            "stage=owen-stratification; dimension={dimension}; bin={bin}; count={count}; expected=1"
        ));
    }

    CaseOutcome::pass(format!(
        "frame_version={INPUT_FRAME_VERSION}; fs_rand={}; stream_semantics={STREAM_SEMANTICS_VERSION}; input_seed=0x{OWEN_SEED:016x}; alternate_seed=0x{OWEN_ALT_SEED:016x}; kernel=0x{OWEN_STREAM_KERNEL:04x}; dim={OWEN_DIM}; points={OWEN_POINTS}; both_seed_replays=exact; alternate_seed=distinct; stratification=m7; output_digest={:016x}",
        fs_rand::VERSION,
        owen_output_digest(&first, &alternate)
    ))
    .with_evidence("crates/fs-rand/CONTRACT.md#qmc-qmc-module")
}

fn gcd(mut a: u32, mut b: u32) -> u32 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

fn lattice_output_digest(lattice: &Lattice, error: f64, points: &[f64]) -> u64 {
    let mut bytes = b"fs-rand:qmc:cbc-output:v1".to_vec();
    bytes.extend_from_slice(&lattice.n.to_le_bytes());
    bytes.extend_from_slice(&usize_u64(lattice.z.len()).to_le_bytes());
    for &generator in &lattice.z {
        bytes.extend_from_slice(&generator.to_le_bytes());
    }
    push_f64(&mut bytes, error);
    bytes.extend_from_slice(&usize_u64(points.len()).to_le_bytes());
    for &value in points {
        push_f64(&mut bytes, value);
    }
    fnv1a64(&bytes)
}

#[allow(clippy::too_many_lines)] // Exhaustive public lattice/state replay audit.
fn lattice_outcome() -> CaseOutcome {
    let first = Lattice::cbc(LATTICE_N, LATTICE_DIM);
    let replay = Lattice::cbc(LATTICE_N, LATTICE_DIM);
    if first.n != replay.n || first.z != replay.z {
        return qmc_failure(format!(
            "stage=cbc-replay; first_n={}; replay_n={}; first_z={:?}; replay_z={:?}",
            first.n, replay.n, first.z, replay.z
        ));
    }
    if first.n != LATTICE_N || first.z.len() != LATTICE_DIM {
        return qmc_failure(format!(
            "stage=public-shape; n={}; z_len={}; expected_n={LATTICE_N}; expected_dim={LATTICE_DIM}",
            first.n,
            first.z.len()
        ));
    }
    for (dimension, &generator) in first.z.iter().enumerate() {
        if !(1..LATTICE_N).contains(&generator) || gcd(generator, LATTICE_N) != 1 {
            return qmc_failure(format!(
                "stage=generator-domain; dimension={dimension}; generator={generator}; n={LATTICE_N}"
            ));
        }
    }

    let error = first.korobov_error_sq();
    let replay_error = replay.korobov_error_sq();
    if error.to_bits() != replay_error.to_bits() || !error.is_finite() || error <= 0.0 {
        return qmc_failure(format!(
            "stage=error-replay; error=0x{:016x}; replay=0x{:016x}",
            error.to_bits(),
            replay_error.to_bits()
        ));
    }

    let point_count = usize::try_from(LATTICE_N).expect("point count fits usize");
    let mut points = Vec::with_capacity(point_count * LATTICE_DIM);
    let mut first_point = vec![0.0; LATTICE_DIM];
    let mut replay_point = vec![0.0; LATTICE_DIM];
    let mut residues = vec![vec![false; point_count]; LATTICE_DIM];
    for point in 0..LATTICE_N {
        first.point(point, &mut first_point);
        replay.point(point, &mut replay_point);
        for (dimension, ((&value, &replayed), &generator)) in first_point
            .iter()
            .zip(&replay_point)
            .zip(&first.z)
            .enumerate()
        {
            let residue =
                u32::try_from(u64::from(point) * u64::from(generator) % u64::from(LATTICE_N))
                    .expect("modular residue fits u32");
            let expected = f64::from(residue) / f64::from(LATTICE_N);
            if value.to_bits() != replayed.to_bits() || value.to_bits() != expected.to_bits() {
                return qmc_failure(format!(
                    "stage=lattice-point; point={point}; dimension={dimension}; generator={generator}; residue={residue}; computed=0x{:016x}; replay=0x{:016x}; expected=0x{:016x}",
                    value.to_bits(),
                    replayed.to_bits(),
                    expected.to_bits()
                ));
            }
            residues[dimension][usize::try_from(residue).expect("residue fits usize")] = true;
            points.push(value);
        }
    }
    if let Some((dimension, residue)) = residues.iter().enumerate().find_map(|(dimension, row)| {
        row.iter()
            .position(|present| !present)
            .map(|residue| (dimension, residue))
    }) {
        return qmc_failure(format!(
            "stage=residue-permutation; dimension={dimension}; missing_residue={residue}"
        ));
    }
    for (&input, &expected) in BAKER_INPUTS.iter().zip(&BAKER_OUTPUTS) {
        let computed = baker(input);
        if computed.to_bits() != expected.to_bits() {
            return qmc_failure(format!(
                "stage=baker-kat; input=0x{:016x}; computed=0x{:016x}; expected=0x{:016x}",
                input.to_bits(),
                computed.to_bits(),
                expected.to_bits()
            ));
        }
    }

    CaseOutcome::pass(format!(
        "frame_version={INPUT_FRAME_VERSION}; fs_rand={}; n={LATTICE_N}; dim={LATTICE_DIM}; z={:?}; error_bits=0x{:016x}; public_points=exact; residue_permutations={LATTICE_DIM}/{LATTICE_DIM}; baker_kat=5/5; output_digest={:016x}",
        fs_rand::VERSION,
        first.z,
        error.to_bits(),
        lattice_output_digest(&first, error, &points)
    ))
    .with_evidence("crates/fs-rand/CONTRACT.md#qmc-qmc-module")
}

fn corrupted_sobol_reference(seed: u64) -> ([f64; 8], usize, u32) {
    let point =
        usize::try_from(seed % usize_u64(SOBOL_KAT.len())).expect("corruption point fits usize");
    let bit = u32::try_from((seed >> 8) & 0x0f).expect("corruption bit fits u32");
    let mut reference = SOBOL_KAT;
    reference[point] = f64::from_bits(reference[point].to_bits() ^ (1_u64 << bit));
    (reference, point, bit)
}

fn corruption_inputs(reference: &[f64; 8], point: usize, bit: u32) -> Vec<u8> {
    let mut bytes = input_header("seeded-sobol-gray-kat-corruption");
    bytes.extend_from_slice(&CORRUPTION_SEED.to_le_bytes());
    bytes.extend_from_slice(&usize_u64(point).to_le_bytes());
    bytes.extend_from_slice(&bit.to_le_bytes());
    bytes.extend_from_slice(&usize_u64(reference.len()).to_le_bytes());
    for &value in reference {
        push_f64(&mut bytes, value);
    }
    bytes
}

fn bits_frame(values: &[f64]) -> String {
    let mut output = String::from("[");
    for (index, value) in values.iter().enumerate() {
        if index > 0 {
            output.push(',');
        }
        write!(&mut output, "0x{:016x}", value.to_bits()).expect("String writes are infallible");
    }
    output.push(']');
    output
}

fn corruption_outcome(reference: [f64; 8], point: usize, bit: u32) -> CaseOutcome {
    let sobol = Sobol::new(1);
    let mut computed = [0.0];
    for (index, &corrupted) in reference.iter().enumerate() {
        sobol.point(
            u32::try_from(index).expect("small KAT index"),
            &mut computed,
        );
        if computed[0].to_bits() != corrupted.to_bits() {
            return CaseOutcome::fail(format!(
                "mode=unscrambled; corruption_seed=0x{CORRUPTION_SEED:016x}; point={point}; bit={bit}; first_mismatch={index}; computed=0x{:016x}; canonical_reference=0x{:016x}; corrupted_reference=0x{:016x}; corrupted_frame={}",
                computed[0].to_bits(),
                SOBOL_KAT[index].to_bits(),
                corrupted.to_bits(),
                bits_frame(&reference)
            ))
            .with_evidence("crates/fs-rand/tests/qmc_replay_casebook.rs#seeded-corruption");
        }
    }
    CaseOutcome::pass("seeded Sobol reference corruption was not detected")
        .with_evidence("crates/fs-rand/tests/qmc_replay_casebook.rs#seeded-corruption")
}

fn seeded_corruption_report() -> SuiteReport {
    let (reference, point, bit) = corrupted_sobol_reference(CORRUPTION_SEED);
    assert_eq!((point, bit), (5, 10));
    assert!(reference.iter().all(|value| value.is_finite()));
    let inputs_digest = fnv1a64(&corruption_inputs(&reference, point, bit));
    assert_eq!(inputs_digest, CORRUPTION_INPUTS_DIGEST);
    Suite::new(SUITE)
        .case(
            "seeded-sobol-gray-reference-corruption",
            inputs_digest,
            ToleranceSpec::Exact,
            move || corruption_outcome(reference, point, bit),
        )
        .run()
}

#[test]
fn qmc_casebook_emits_replay_complete_green_records() {
    let sobol_digest = fnv1a64(&sobol_inputs());
    let owen_digest = fnv1a64(&owen_inputs());
    let lattice_digest = fnv1a64(&lattice_inputs());
    assert_eq!(sobol_digest, SOBOL_INPUTS_DIGEST);
    assert_eq!(owen_digest, OWEN_INPUTS_DIGEST);
    assert_eq!(lattice_digest, LATTICE_INPUTS_DIGEST);

    let report = Suite::new(SUITE)
        .case(
            "sobol-gray-random-access",
            sobol_digest,
            ToleranceSpec::Exact,
            sobol_outcome,
        )
        .case(
            "owen-seeded-replay",
            owen_digest,
            ToleranceSpec::Exact,
            owen_outcome,
        )
        .case(
            "cbc-lattice-replay",
            lattice_digest,
            ToleranceSpec::Exact,
            lattice_outcome,
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
            "sobol-gray-random-access",
            "owen-seeded-replay",
            "cbc-lattice-replay",
        ]
    );
    assert!(
        report
            .records
            .iter()
            .all(|record| record.version == CASEBOOK_RECORD_VERSION && record.pass)
    );
    assert_eq!(report.records[0].inputs_digest, "e19f4e88d96870ea");
    assert_eq!(report.records[1].inputs_digest, "d52f0c160b189289");
    assert_eq!(report.records[2].inputs_digest, "08cff5026b213c41");
    for record in &report.records {
        let line = record.json_line();
        assert!(line.contains("\"tolerance\":\"exact\",\"pass\":true"));
        assert!(line.contains("\"evidence\":[\"crates/fs-rand/CONTRACT.md#qmc-qmc-module\"]"));
    }
}

#[test]
fn disclosed_seeded_qmc_corruption_replays_and_is_refused() {
    let (first_reference, first_point, first_bit) = corrupted_sobol_reference(CORRUPTION_SEED);
    let (replay_reference, replay_point, replay_bit) = corrupted_sobol_reference(CORRUPTION_SEED);
    let first_frame = corruption_inputs(&first_reference, first_point, first_bit);
    let replay_frame = corruption_inputs(&replay_reference, replay_point, replay_bit);
    assert_eq!(
        first_frame, replay_frame,
        "the corruption seed must independently reconstruct identical canonical input bytes"
    );
    assert_eq!(fnv1a64(&first_frame), CORRUPTION_INPUTS_DIGEST);

    let first = seeded_corruption_report();
    let replay = seeded_corruption_report();
    let first_failures = first.failures();
    let replay_failures = replay.failures();
    let [first_failure] = first_failures.as_slice() else {
        panic!("the disclosed QMC corruption must produce exactly one failure");
    };
    let [replay_failure] = replay_failures.as_slice() else {
        panic!("the replayed QMC corruption must produce exactly one failure");
    };
    assert_eq!(first_failure.json_line(), replay_failure.json_line());
    assert_eq!(first_failure.inputs_digest, "80d975b52d579c57");
    assert!(first_failure.details.contains("mode=unscrambled"));
    assert!(
        first_failure
            .details
            .contains(&format!("corruption_seed=0x{CORRUPTION_SEED:016x}"))
    );
    assert!(first_failure.details.contains("point=5; bit=10"));
    assert!(first_failure.details.contains("first_mismatch=5"));
    assert!(first_failure.details.contains("corrupted_frame=["));

    let panic = catch_unwind(|| first.assert_green())
        .expect_err("the merge gate must reject the disclosed QMC corruption");
    let message = panic
        .downcast_ref::<String>()
        .map(String::as_str)
        .or_else(|| panic.downcast_ref::<&str>().copied())
        .expect("Casebook panic carries text");
    assert!(message.contains("seeded-sobol-gray-reference-corruption"));
    assert!(message.contains(&format!("0x{CORRUPTION_SEED:016x}")));
}
