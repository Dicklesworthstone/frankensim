//! Structured ASCENT conformance records for quality-diversity archives.
//!
//! The larger archive battery keeps broad behavioral coverage. These cases
//! expose the load-bearing deterministic laws through the shared fs-casebook
//! schema with canonical inputs and replay-complete failure diagnostics.

use core::fmt::Write as _;
use fs_archive::{CvtArchive, MapElites, novelty};
use fs_casebook::{CASEBOOK_RECORD_VERSION, CaseOutcome, Suite, ToleranceSpec, fnv1a64};

const SUITE: &str = "fs-archive/ascent-conformance-v1";
const MAP_EXPECTED: [bool; 5] = [true, false, true, false, true];
const CVT_ADD_EXPECTED: [bool; 4] = [true, false, true, true];

#[derive(Clone, Copy)]
struct MapOp {
    solution: f64,
    descriptor: [f64; 2],
    fitness: f64,
}

const MAP_OPS: [MapOp; 5] = [
    MapOp {
        solution: 1.0,
        descriptor: [0.10, 0.10],
        fitness: 5.0,
    },
    MapOp {
        solution: 2.0,
        descriptor: [0.20, 0.20],
        fitness: 3.0,
    },
    MapOp {
        solution: 3.0,
        descriptor: [0.24, 0.24],
        fitness: 9.0,
    },
    MapOp {
        solution: 4.0,
        descriptor: [0.23, 0.23],
        fitness: 9.0,
    },
    MapOp {
        solution: 5.0,
        descriptor: [0.90, 0.90],
        fitness: 1.0,
    },
];

const CENTROIDS: [[f64; 2]; 4] = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]];
const CVT_PROBES: [([f64; 2], usize); 5] = [
    ([0.10, 0.10], 0),
    ([0.90, 0.10], 1),
    ([0.10, 0.90], 2),
    ([0.50, 0.00], 0),
    ([0.50, 1.00], 2),
];
const CVT_OPS: [MapOp; 4] = [
    MapOp {
        solution: 1.0,
        descriptor: [0.10, 0.10],
        fitness: 4.0,
    },
    MapOp {
        solution: 2.0,
        descriptor: [0.20, 0.10],
        fitness: 4.0,
    },
    MapOp {
        solution: 3.0,
        descriptor: [0.20, 0.10],
        fitness: 6.0,
    },
    MapOp {
        solution: 4.0,
        descriptor: [0.90, 0.90],
        fitness: 2.0,
    },
];
const NOVELTY_OTHERS: [[f64; 2]; 3] = [[0.0, 0.0], [3.0, 4.0], [6.0, 8.0]];
const NOVELTY_CASES: [([f64; 2], usize, f64); 3] = [
    ([0.0, 0.0], 2, 2.5),
    ([3.0, 4.0], 2, 2.5),
    ([9.0, 12.0], 2, 7.5),
];

fn push_u64(bytes: &mut Vec<u8>, value: u64) {
    bytes.extend_from_slice(&value.to_le_bytes());
}

fn push_f64(bytes: &mut Vec<u8>, value: f64) {
    push_u64(bytes, value.to_bits());
}

fn push_len(bytes: &mut Vec<u8>, value: usize) {
    push_u64(
        bytes,
        u64::try_from(value).expect("conformance fixture lengths fit u64"),
    );
}

fn hex_bytes(bytes: &[u8]) -> String {
    let mut encoded = String::with_capacity(bytes.len() * 2);
    for byte in bytes {
        write!(&mut encoded, "{byte:02x}").expect("writing to String cannot fail");
    }
    encoded
}

fn append_map_fixture(bytes: &mut Vec<u8>, expected: [bool; 5]) {
    push_len(bytes, 2);
    for value in [0.0, 0.0, 1.0, 1.0] {
        push_f64(bytes, value);
    }
    push_u64(bytes, 4);
    push_u64(bytes, 4);
    push_len(bytes, MAP_OPS.len());
    for (op, accepted) in MAP_OPS.into_iter().zip(expected) {
        push_f64(bytes, op.solution);
        for value in op.descriptor {
            push_f64(bytes, value);
        }
        push_f64(bytes, op.fitness);
        bytes.push(u8::from(accepted));
    }
}

fn map_inputs() -> Vec<u8> {
    let mut bytes = b"fs-archive:map-elites-conformance:v1".to_vec();
    append_map_fixture(&mut bytes, MAP_EXPECTED);
    bytes
}

fn cvt_inputs() -> Vec<u8> {
    let mut bytes = b"fs-archive:cvt-conformance:v1".to_vec();
    push_len(&mut bytes, CENTROIDS.len());
    push_len(&mut bytes, 2);
    for centroid in CENTROIDS {
        for value in centroid {
            push_f64(&mut bytes, value);
        }
    }
    push_len(&mut bytes, CVT_PROBES.len());
    for (descriptor, expected) in CVT_PROBES {
        for value in descriptor {
            push_f64(&mut bytes, value);
        }
        push_len(&mut bytes, expected);
    }
    push_len(&mut bytes, CVT_OPS.len());
    for (op, accepted) in CVT_OPS.into_iter().zip(CVT_ADD_EXPECTED) {
        push_f64(&mut bytes, op.solution);
        for value in op.descriptor {
            push_f64(&mut bytes, value);
        }
        push_f64(&mut bytes, op.fitness);
        bytes.push(u8::from(accepted));
    }
    bytes
}

fn novelty_inputs() -> Vec<u8> {
    let mut bytes = b"fs-archive:novelty-known-answers:v1".to_vec();
    push_len(&mut bytes, NOVELTY_OTHERS.len());
    push_len(&mut bytes, 2);
    for other in NOVELTY_OTHERS {
        for value in other {
            push_f64(&mut bytes, value);
        }
    }
    push_len(&mut bytes, NOVELTY_CASES.len());
    for (query, k, expected) in NOVELTY_CASES {
        for value in query {
            push_f64(&mut bytes, value);
        }
        push_len(&mut bytes, k);
        push_f64(&mut bytes, expected);
    }
    // Bind both +infinity edge claims as explicit framed cases.
    push_len(&mut bytes, 2);
    bytes.push(0); // empty-neighbour case
    push_len(&mut bytes, 1);
    push_f64(&mut bytes, 0.0);
    push_len(&mut bytes, 0);
    push_len(&mut bytes, 2);
    push_f64(&mut bytes, f64::INFINITY);
    bytes.push(1); // k=0 case
    push_len(&mut bytes, 1);
    push_f64(&mut bytes, 0.0);
    push_len(&mut bytes, 1);
    push_f64(&mut bytes, 0.0);
    push_len(&mut bytes, 0);
    push_f64(&mut bytes, f64::INFINITY);
    bytes
}

fn grid() -> MapElites {
    MapElites::new(vec![0.0, 0.0], vec![1.0, 1.0], vec![4, 4])
}

fn map_elites_outcome() -> CaseOutcome {
    let inputs_hex = hex_bytes(&map_inputs());
    let mut archive = grid();
    let (mut prior_coverage, mut prior_qd) = (0.0, 0.0);
    for (index, (op, expected)) in MAP_OPS.into_iter().zip(MAP_EXPECTED).enumerate() {
        let accepted = archive.add(vec![op.solution], op.descriptor.to_vec(), op.fitness);
        if accepted != expected {
            return CaseOutcome::fail(format!(
                "operation={index}; descriptor_bits={:016x?}; fitness_bits=0x{:016x}; accepted={accepted}; expected={expected}; inputs_hex={inputs_hex}",
                op.descriptor.map(f64::to_bits),
                op.fitness.to_bits(),
            ))
            .with_evidence("crates/fs-archive/CONTRACT.md#invariants");
        }
        let (coverage, qd) = (archive.coverage(), archive.qd_score());
        if coverage < prior_coverage || qd < prior_qd {
            return CaseOutcome::fail(format!(
                "operation={index}; descriptor_bits={:016x?}; fitness_bits=0x{:016x}; prior_coverage_bits=0x{:016x}; coverage_bits=0x{:016x}; prior_qd_bits=0x{:016x}; qd_bits=0x{:016x}; inputs_hex={inputs_hex}",
                op.descriptor.map(f64::to_bits),
                op.fitness.to_bits(),
                prior_coverage.to_bits(),
                coverage.to_bits(),
                prior_qd.to_bits(),
                qd.to_bits(),
            ))
            .with_evidence("crates/fs-archive/CONTRACT.md#invariants");
        }
        (prior_coverage, prior_qd) = (coverage, qd);
    }

    let Some(elite) = archive.elite_at(&[0.10, 0.10]) else {
        return CaseOutcome::fail(format!(
            "replacement_niche_query_bits=[3fb999999999999a,3fb999999999999a]; elite=missing; inputs_hex={inputs_hex}"
        ))
            .with_evidence("crates/fs-archive/CONTRACT.md#invariants");
    };
    let final_facts = (
        archive.num_elites(),
        archive.coverage().to_bits(),
        archive.qd_score().to_bits(),
        archive.best().map(|best| best.fitness.to_bits()),
        elite.solution.first().map(|value| value.to_bits()),
    );
    let expected = (
        2,
        0.125_f64.to_bits(),
        10.0_f64.to_bits(),
        Some(9.0_f64.to_bits()),
        Some(3.0_f64.to_bits()),
    );
    if final_facts != expected {
        return CaseOutcome::fail(format!(
            "final_facts={final_facts:016x?}; expected={expected:016x?}; inputs_hex={inputs_hex}"
        ))
        .with_evidence("crates/fs-archive/CONTRACT.md#invariants");
    }
    CaseOutcome::pass(
        "operations=5; accepted=[true,false,true,false,true]; elites=2; coverage_bits=3fc0000000000000; qd_score_bits=4024000000000000; best_fitness_bits=4022000000000000",
    )
    .with_evidence("crates/fs-archive/CONTRACT.md#invariants")
}

fn cvt_outcome() -> CaseOutcome {
    let inputs_hex = hex_bytes(&cvt_inputs());
    let mut archive = CvtArchive::new(CENTROIDS.into_iter().map(|row| row.to_vec()).collect());
    for (probe, expected) in CVT_PROBES {
        let actual = archive.nearest_centroid(&probe);
        if actual != expected {
            return CaseOutcome::fail(format!(
                "probe_bits={:016x?}; nearest={actual}; expected={expected}; centroids_bits={:016x?}; inputs_hex={inputs_hex}",
                probe.map(f64::to_bits),
                CENTROIDS.map(|centroid| centroid.map(f64::to_bits)),
            ))
            .with_evidence("crates/fs-archive/CONTRACT.md#invariants");
        }
    }
    for (index, (op, expected)) in CVT_OPS.into_iter().zip(CVT_ADD_EXPECTED).enumerate() {
        let accepted = archive.add(vec![op.solution], op.descriptor.to_vec(), op.fitness);
        if accepted != expected {
            return CaseOutcome::fail(format!(
                "operation={index}; descriptor_bits={:016x?}; fitness_bits=0x{:016x}; accepted={accepted}; expected={expected}; inputs_hex={inputs_hex}",
                op.descriptor.map(f64::to_bits),
                op.fitness.to_bits(),
            ))
            .with_evidence("crates/fs-archive/CONTRACT.md#invariants");
        }
    }
    let final_facts = (
        archive.num_elites(),
        archive.coverage().to_bits(),
        archive.qd_score().to_bits(),
        archive.best().map(|best| best.fitness.to_bits()),
    );
    let expected = (
        2,
        0.5_f64.to_bits(),
        8.0_f64.to_bits(),
        Some(6.0_f64.to_bits()),
    );
    if final_facts != expected {
        return CaseOutcome::fail(format!(
            "final_facts={final_facts:016x?}; expected={expected:016x?}; inputs_hex={inputs_hex}"
        ))
        .with_evidence("crates/fs-archive/CONTRACT.md#invariants");
    }
    CaseOutcome::pass(
        "probes=5; ties=[0,2]; accepted=[true,false,true,true]; elites=2; coverage_bits=3fe0000000000000; qd_score_bits=4020000000000000; best_fitness_bits=4018000000000000",
    )
    .with_evidence("crates/fs-archive/CONTRACT.md#invariants")
}

fn novelty_outcome() -> CaseOutcome {
    let inputs_hex = hex_bytes(&novelty_inputs());
    let others = NOVELTY_OTHERS.map(|row| row.to_vec());
    for (index, (query, k, expected)) in NOVELTY_CASES.into_iter().enumerate() {
        let actual = novelty(&query, &others, k);
        if actual.to_bits() != expected.to_bits() {
            return CaseOutcome::fail(format!(
                "case={index}; query_bits={:016x?}; neighbours_bits={:016x?}; k={k}; actual_bits=0x{:016x}; expected_bits=0x{:016x}; inputs_hex={inputs_hex}",
                query.map(f64::to_bits),
                NOVELTY_OTHERS.map(|row| row.map(f64::to_bits)),
                actual.to_bits(),
                expected.to_bits(),
            ))
            .with_evidence("crates/fs-archive/CONTRACT.md#invariants");
        }
    }
    let empty = novelty(&[0.0], &[], 2);
    let k_zero = novelty(&[0.0], &[vec![0.0]], 0);
    if empty.to_bits() != f64::INFINITY.to_bits() || k_zero.to_bits() != f64::INFINITY.to_bits() {
        return CaseOutcome::fail(format!(
            "empty_case={{query_bits:[0000000000000000],neighbours:[],k:2,actual_bits:{:016x},expected_bits:7ff0000000000000}}; k0_case={{query_bits:[0000000000000000],neighbours:[[0000000000000000]],k:0,actual_bits:{:016x},expected_bits:7ff0000000000000}}; inputs_hex={inputs_hex}",
            empty.to_bits(),
            k_zero.to_bits(),
        ))
        .with_evidence("crates/fs-archive/CONTRACT.md#public-types-and-semantics");
    }
    CaseOutcome::pass(
        "known_answers_bits=[4004000000000000,4004000000000000,401e000000000000]; empty=+inf; k0=+inf",
    )
    .with_evidence("crates/fs-archive/CONTRACT.md#invariants")
}

#[test]
fn ascent_casebook_suite_emits_replay_complete_green_records() {
    let map_digest = fnv1a64(&map_inputs());
    let cvt_digest = fnv1a64(&cvt_inputs());
    let novelty_digest = fnv1a64(&novelty_inputs());
    assert_eq!(map_digest, 0x2962_c00f_dd1e_dc9b);
    assert_eq!(cvt_digest, 0x818f_e6bc_c470_bb12);
    assert_eq!(novelty_digest, 0x5a40_84a0_5b59_acc3);
    assert_eq!(CASEBOOK_RECORD_VERSION, 1);
    let report = Suite::new(SUITE)
        .case(
            "map-elites-strict-replacement",
            map_digest,
            ToleranceSpec::Exact,
            map_elites_outcome,
        )
        .case(
            "cvt-assignment-and-replacement",
            cvt_digest,
            ToleranceSpec::Exact,
            cvt_outcome,
        )
        .case(
            "novelty-known-answers",
            novelty_digest,
            ToleranceSpec::Exact,
            novelty_outcome,
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
            "map-elites-strict-replacement",
            "cvt-assignment-and-replacement",
            "novelty-known-answers"
        ]
    );
    assert_eq!(
        report.records[0].json_line(),
        concat!(
            "{\"casebook\":1,\"suite\":\"fs-archive/ascent-conformance-v1\",",
            "\"case\":\"map-elites-strict-replacement\",\"inputs_digest\":\"2962c00fdd1edc9b\",",
            "\"tolerance\":\"exact\",\"pass\":true,",
            "\"details\":\"operations=5; accepted=[true,false,true,false,true]; elites=2; coverage_bits=3fc0000000000000; qd_score_bits=4024000000000000; best_fitness_bits=4022000000000000\",",
            "\"evidence\":[\"crates/fs-archive/CONTRACT.md#invariants\"]}"
        )
    );
}

#[test]
fn disclosed_seeded_corruption_turns_archive_casebook_red() {
    const CORRUPTION_SEED: u64 = 0xA7C4_0001;
    let operation_count = u64::try_from(MAP_EXPECTED.len()).expect("the fixture length fits u64");
    let operation = usize::try_from(CORRUPTION_SEED % operation_count)
        .expect("the derived fixture index fits usize");
    assert_eq!(operation, 4);
    let mut corrupted_expected = MAP_EXPECTED;
    corrupted_expected[operation] = !corrupted_expected[operation];
    let mut inputs = b"fs-archive:seeded-map-acceptance-corruption:v1".to_vec();
    push_u64(&mut inputs, CORRUPTION_SEED);
    append_map_fixture(&mut inputs, corrupted_expected);
    let inputs_digest = fnv1a64(&inputs);
    assert_eq!(inputs_digest, 0x6772_924a_4426_af40);
    let inputs_hex = hex_bytes(&inputs);

    let report = Suite::new(SUITE)
        .case(
            "seeded-map-acceptance-corruption",
            inputs_digest,
            ToleranceSpec::Exact,
            move || {
                let mut archive = grid();
                for (index, (op, corrupted)) in MAP_OPS
                    .into_iter()
                    .zip(corrupted_expected)
                    .enumerate()
                {
                    let actual = archive.add(vec![op.solution], op.descriptor.to_vec(), op.fitness);
                    if actual != corrupted {
                        return CaseOutcome::fail(format!(
                            "seed=0x{CORRUPTION_SEED:016x}; operation={index}; descriptor_bits={:016x?}; fitness_bits=0x{:016x}; actual={actual}; corrupted_expected={corrupted}; canonical_expected={}; inputs_hex={inputs_hex}",
                            op.descriptor.map(f64::to_bits),
                            op.fitness.to_bits(),
                            MAP_EXPECTED[index],
                        ))
                        .with_evidence("crates/fs-archive/tests/conformance.rs::disclosed_seeded_corruption_turns_archive_casebook_red");
                    }
                }
                CaseOutcome::pass("seeded acceptance corruption escaped detection")
            },
        )
        .run();

    assert!(!report.all_passed(), "the corrupted oracle must turn red");
    let failures = report.failures();
    let [failure] = failures.as_slice() else {
        panic!("the seeded corruption must produce exactly one failure");
    };
    assert_eq!(failure.case, "seeded-map-acceptance-corruption");
    assert_eq!(failure.inputs_digest, "6772924a4426af40");
    assert!(
        failure
            .details
            .contains(&format!("seed=0x{CORRUPTION_SEED:016x}"))
    );
    assert!(failure.details.contains(&format!("operation={operation}")));
    let selected = MAP_OPS[operation];
    assert!(failure.details.contains(&format!(
        "descriptor_bits={:016x?}",
        selected.descriptor.map(f64::to_bits)
    )));
    assert!(failure.details.contains(&format!(
        "fitness_bits=0x{:016x}",
        selected.fitness.to_bits()
    )));
    assert!(failure.details.contains("actual=") && failure.details.contains("corrupted_expected="));
    assert!(failure.details.contains("canonical_expected=true"));
    assert!(failure.details.contains("inputs_hex="));
    assert!(
        failure
            .json_line()
            .contains("\"tolerance\":\"exact\",\"pass\":false")
    );

    let panic = std::panic::catch_unwind(|| report.assert_green())
        .expect_err("the merge-gate assertion must reject the seeded failure");
    let message = panic
        .downcast_ref::<String>()
        .map(String::as_str)
        .or_else(|| panic.downcast_ref::<&str>().copied())
        .expect("casebook panic carries text");
    assert!(message.contains("seeded-map-acceptance-corruption"));
    assert!(message.contains(&format!("seed=0x{CORRUPTION_SEED:016x}")));
}
