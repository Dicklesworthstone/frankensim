//! fs-casebook self-conformance (bead huq.5): the demo suite every agent
//! copies to add cases, plus the fail-closed self-tests the bead's
//! acceptance demands (an intentionally failing case must yield a
//! structured failure record and a non-green report).

use fs_casebook::{
    CASEBOOK_DISAGREEMENT_RECORD_VERSION, CASEBOOK_REPLAY_RECORD_VERSION, CaseOutcome, CaseRecord,
    DisagreementRecord, ReplayError, ReplaySpec, Suite, ToleranceSpec, fnv1a64,
};

/// The demo suite: exactly how a crate registers conformance cases.
/// Each case binds (stable id, inputs digest, tolerance claim, runner).
fn demo_suite() -> Suite {
    let fixture = b"casebook-demo-fixture-v1";
    Suite::new("casebook-demo")
        .case(
            "cb-demo-001-exact-roundtrip",
            fnv1a64(fixture),
            ToleranceSpec::Exact,
            || {
                let digest = fnv1a64(b"payload");
                let replay = fnv1a64(b"payload");
                if digest == replay {
                    CaseOutcome::pass(format!("digest {digest:#018x} replay-stable"))
                        .with_evidence("fixture:casebook-demo-fixture-v1")
                } else {
                    CaseOutcome::fail("digest not replay-stable")
                }
            },
        )
        .case(
            "cb-demo-002-tolerance-claim",
            fnv1a64(fixture),
            ToleranceSpec::AbsoluteLe(1e-12),
            || {
                let measured = (0.1f64 + 0.2) - 0.3;
                if measured.abs() <= 1e-12 {
                    CaseOutcome::pass(format!("measured {measured:e}"))
                } else {
                    CaseOutcome::fail(format!("measured {measured:e} above bound"))
                }
            },
        )
        .case(
            "cb-demo-003-structural-refusal",
            fnv1a64(b"refusal-probe"),
            ToleranceSpec::Structural,
            || {
                let refused = "".parse::<u32>().is_err();
                if refused {
                    CaseOutcome::pass("malformed input refused as typed error")
                } else {
                    CaseOutcome::fail("malformed input silently accepted")
                }
            },
        )
}

#[test]
fn demo_suite_runs_green_with_structured_records() {
    let report = demo_suite().run();
    report.assert_green();
    assert_eq!(report.records.len(), 3);
    assert!(
        report.replay_records.is_empty() && report.disagreements.is_empty(),
        "legacy Suite::case emits exactly its legacy case records"
    );
    for record in &report.records {
        assert_eq!(record.version, 1);
        assert_eq!(record.suite, "casebook-demo");
        assert_eq!(record.inputs_digest.len(), 16, "digest is 16 hex chars");
        let line = record.json_line();
        assert!(line.starts_with("{\"casebook\":1,\"suite\":\"casebook-demo\""));
        assert!(line.contains(&format!("\"case\":\"{}\"", record.case)));
        assert!(line.contains("\"pass\":true"));
        assert!(line.ends_with('}'));
    }
    assert!(
        report.records[0].evidence == vec!["fixture:casebook-demo-fixture-v1".to_owned()],
        "evidence pointers ride the record"
    );
}

#[test]
fn intentionally_failing_case_yields_structured_failure_and_non_green() {
    let report = Suite::new("casebook-selftest")
        .case("cb-self-001-pass", 1, ToleranceSpec::Structural, || {
            CaseOutcome::pass("ok")
        })
        .case(
            "cb-self-002-intentional-failure",
            2,
            ToleranceSpec::RelativeLe(0.01),
            || CaseOutcome::fail("measured 0.5 above declared 0.01").with_evidence("log:selftest"),
        )
        .run();
    assert!(!report.all_passed(), "failing case must not read green");
    let failures = report.failures();
    assert_eq!(failures.len(), 1);
    let record = failures[0];
    assert_eq!(record.case, "cb-self-002-intentional-failure");
    assert!(!record.pass);
    assert_eq!(record.tolerance, "rel<=1e-2");
    assert_eq!(record.evidence, vec!["log:selftest".to_owned()]);
    let line = record.json_line();
    assert!(line.contains("\"pass\":false"));
    assert!(line.contains("measured 0.5 above declared 0.01"));

    let panic = std::panic::catch_unwind(|| report.assert_green());
    let message = match panic {
        Ok(()) => panic!("assert_green accepted a failing suite"),
        Err(payload) => payload
            .downcast_ref::<String>()
            .cloned()
            .unwrap_or_default(),
    };
    assert!(
        message.contains("cb-self-002-intentional-failure"),
        "panic message must carry the failing record: {message}"
    );
}

#[test]
fn empty_suites_and_duplicate_ids_fail_closed() {
    let empty = Suite::new("casebook-empty").run();
    assert!(
        !empty.all_passed(),
        "a suite that ran nothing proved nothing"
    );

    let duplicated = Suite::new("casebook-dup")
        .case("cb-dup-001", 1, ToleranceSpec::Exact, || {
            CaseOutcome::pass("first")
        })
        .case("cb-dup-001", 2, ToleranceSpec::Exact, || {
            CaseOutcome::pass("second registration must not run as green")
        })
        .run();
    assert!(!duplicated.all_passed());
    let failures = duplicated.failures();
    assert_eq!(failures.len(), 1);
    assert!(failures[0].details.contains("duplicate case id"));
}

#[test]
fn duplicate_ids_cannot_cross_associate_replays_or_disagreements() {
    let first_disagreement = DisagreementRecord::first(
        "foreign-suite",
        "foreign-case",
        "implementation",
        "reference",
        b"first-implementation",
        b"first-reference",
    )
    .expect("first registration has a disagreement");
    let report = Suite::new("casebook-duplicate-owner")
        .case_replayable(
            "cb-duplicate-owner",
            ReplaySpec::new("run-first-registration", b"first-input".to_vec()),
            ToleranceSpec::Exact,
            move || {
                CaseOutcome::pass("first registration is red").with_disagreement(first_disagreement)
            },
        )
        .case_replayable(
            "cb-duplicate-owner",
            ReplaySpec::new("run-duplicate-registration", b"duplicate-input".to_vec()),
            ToleranceSpec::Exact,
            || panic!("a duplicate registration must not execute its case closure"),
        )
        .run();

    assert_eq!(report.failures().len(), 2);
    assert_eq!(report.replay_records.len(), 2);
    assert_eq!(report.disagreements.len(), 1);
    assert!(
        report.replay_for("cb-duplicate-owner").is_none(),
        "a string-only lookup must refuse an ambiguous duplicate id"
    );
    assert!(
        report.disagreements_for("cb-duplicate-owner").is_empty(),
        "a string-only lookup must not merge distinct registrations"
    );

    let panic = std::panic::catch_unwind(|| report.assert_green())
        .expect_err("both owning registrations are red");
    let message = panic
        .downcast_ref::<String>()
        .map(String::as_str)
        .or_else(|| panic.downcast_ref::<&str>().copied())
        .expect("Casebook refusal carries text");
    assert_eq!(message.matches("run-first-registration").count(), 1);
    assert_eq!(message.matches("run-duplicate-registration").count(), 1);
    assert_eq!(message.matches("\"casebook_disagreement\":1").count(), 1);
}

#[test]
fn json_escaping_keeps_records_one_line_and_parseable_shape() {
    let record = CaseRecord {
        version: 1,
        suite: "quote\"suite".to_owned(),
        case: "line\nbreak\tand\\slash".to_owned(),
        inputs_digest: "00000000deadbeef".to_owned(),
        tolerance: ToleranceSpec::Ulps(3).to_string(),
        pass: false,
        details: "control\u{1}char".to_owned(),
        evidence: vec!["path\\with\"quotes".to_owned()],
    };
    let line = record.json_line();
    assert!(!line.contains('\n'), "records are strictly one line");
    assert!(line.contains("quote\\\"suite"));
    assert!(line.contains("line\\nbreak\\tand\\\\slash"));
    assert!(line.contains("\\u0001"));
    assert!(line.contains("ulps<=3"));
    assert_eq!(
        line,
        "{\"casebook\":1,\"suite\":\"quote\\\"suite\",\"case\":\"line\\nbreak\\tand\\\\slash\",\"inputs_digest\":\"00000000deadbeef\",\"tolerance\":\"ulps<=3\",\"pass\":false,\"details\":\"control\\u0001char\",\"evidence\":[\"path\\\\with\\\"quotes\"]}"
    );
}

#[test]
fn digest_helper_is_the_canonical_fnv1a64() {
    assert_eq!(fnv1a64(b""), 0xcbf2_9ce4_8422_2325);
    assert_ne!(fnv1a64(b"a"), fnv1a64(b"b"));
    assert_eq!(fnv1a64(b"casebook"), fnv1a64(b"casebook"));
}

#[test]
fn replay_spec_round_trips_and_records_verify_frame_integrity() {
    let frame = b"casebook-replay-canonical-frame-v1".to_vec();
    let command = "cargo test -p fs-casebook --test casebook replay_spec_round_trips -- --exact";
    let spec = ReplaySpec::new(command, frame.clone());
    assert_eq!(spec.command(), command);
    assert_eq!(spec.canonical_inputs(), frame.as_slice());
    assert_eq!(spec.inputs_digest(), fnv1a64(&frame));

    let encoded = spec.canonical_inputs_hex();
    let reconstructed = ReplaySpec::from_hex(command, &encoded).expect("canonical hex decodes");
    assert_eq!(reconstructed, spec);

    let first = Suite::new("casebook-replay-selftest")
        .case_replayable("cb-replay-001", spec.clone(), ToleranceSpec::Exact, || {
            CaseOutcome::pass("replay companion retained")
        })
        .run();
    let replay = Suite::new("casebook-replay-selftest")
        .case_replayable("cb-replay-001", spec, ToleranceSpec::Exact, || {
            CaseOutcome::pass("replay companion retained")
        })
        .run();
    assert!(first.all_passed() && replay.all_passed());
    assert_eq!(
        first.records[0].inputs_digest,
        format!("{:016x}", fnv1a64(&frame))
    );
    let first_record = first
        .replay_for("cb-replay-001")
        .expect("replayable case emits a companion");
    let replay_record = replay
        .replay_for("cb-replay-001")
        .expect("replay emits the same companion");
    assert_eq!(first_record.version, CASEBOOK_REPLAY_RECORD_VERSION);
    assert_eq!(
        first_record.verify_and_decode().expect("record verifies"),
        frame
    );
    assert_eq!(first_record.json_line(), replay_record.json_line());
    assert_eq!(
        first_record.json_line(),
        "{\"casebook_replay\":1,\"suite\":\"casebook-replay-selftest\",\"case\":\"cb-replay-001\",\"command\":\"cargo test -p fs-casebook --test casebook replay_spec_round_trips -- --exact\",\"inputs_digest\":\"01ff42e80c555e77\",\"inputs_len\":34,\"canonical_inputs_hex\":\"63617365626f6f6b2d7265706c61792d63616e6f6e6963616c2d6672616d652d7631\"}"
    );

    let mut wrong_length = first_record.clone();
    wrong_length.inputs_len += 1;
    assert!(matches!(
        wrong_length.verify_and_decode(),
        Err(ReplayError::LengthMismatch { .. })
    ));

    let mut wrong_digest = first_record.clone();
    wrong_digest.inputs_digest = "0000000000000000".to_owned();
    assert!(matches!(
        wrong_digest.verify_and_decode(),
        Err(ReplayError::DigestMismatch { .. })
    ));
    let mut wrong_version = first_record.clone();
    wrong_version.version += 1;
    assert!(matches!(
        wrong_version.verify_and_decode(),
        Err(ReplayError::UnsupportedVersion { .. })
    ));
    let mut wrong_case = first_record.clone();
    wrong_case.case = "cb-replay-foreign".to_owned();
    assert!(matches!(
        wrong_case.verify_and_decode_for("casebook-replay-selftest", "cb-replay-001"),
        Err(ReplayError::IdentityMismatch { field: "case", .. })
    ));
    let mut empty_command = first_record.clone();
    empty_command.command.clear();
    assert!(matches!(
        empty_command.verify_and_decode(),
        Err(ReplayError::EmptyField { field: "command" })
    ));
    let mut uppercase_frame = first_record.clone();
    uppercase_frame.canonical_inputs_hex = uppercase_frame.canonical_inputs_hex.to_uppercase();
    assert!(matches!(
        uppercase_frame.verify_and_decode(),
        Err(ReplayError::NonCanonicalHex { .. })
    ));
    assert!(matches!(
        ReplaySpec::from_hex(command, "0"),
        Err(ReplayError::OddHexLength { length: 1 })
    ));
    assert!(matches!(
        ReplaySpec::from_hex(command, "0g"),
        Err(ReplayError::InvalidHex {
            offset: 1,
            byte: b'g'
        })
    ));
}

#[test]
fn a_tampered_replay_companion_cannot_leave_its_report_green() {
    let mut report = Suite::new("casebook-report-integrity")
        .case_replayable(
            "cb-report-integrity",
            ReplaySpec::new("run-report-integrity", b"bound-input".to_vec()),
            ToleranceSpec::Exact,
            || CaseOutcome::pass("case execution passed"),
        )
        .run();
    assert!(report.all_passed());

    report.replay_records[0].case = "foreign-case".to_owned();
    assert!(
        !report.all_passed(),
        "a valid case row with a foreign replay identity is not green"
    );
    let panic = std::panic::catch_unwind(|| report.assert_green())
        .expect_err("report integrity failure must trip the merge gate");
    let message = panic
        .downcast_ref::<String>()
        .map(String::as_str)
        .or_else(|| panic.downcast_ref::<&str>().copied())
        .expect("Casebook refusal carries text");
    assert!(message.contains("report-integrity failure"));
    assert!(message.contains("replay record case mismatch"));
}

#[test]
fn identical_frames_have_no_disagreement() {
    let frame = b"identical-frame";
    assert!(
        DisagreementRecord::first(
            "casebook-disagreement-selftest",
            "cb-disagreement-identical",
            "implementation",
            "reference",
            frame,
            frame,
        )
        .is_none()
    );
}

#[test]
fn seeded_byte_corruption_emits_stable_bug_report_and_fails_closed() {
    const SEED: u64 = 0xF5CB_0023_0000_0103;
    let reference = b"casebook-disagreement-reference-v1".to_vec();
    let offset = usize::try_from(SEED % reference.len() as u64).expect("offset fits usize");
    let bit = u32::try_from((SEED >> 8) % 8).expect("bit fits u32");
    let mut implementation = reference.clone();
    implementation[offset] ^= 1_u8 << bit;

    let first = DisagreementRecord::first(
        "casebook-disagreement-selftest",
        "cb-disagreement-seeded-byte",
        "frankensim",
        "frankenscipy",
        &implementation,
        &reference,
    )
    .expect("seeded corruption disagrees");
    let replay = DisagreementRecord::first(
        "casebook-disagreement-selftest",
        "cb-disagreement-seeded-byte",
        "frankensim",
        "frankenscipy",
        &implementation,
        &reference,
    )
    .expect("replayed corruption disagrees");
    assert_eq!(first, replay);
    assert_eq!(first.version(), CASEBOOK_DISAGREEMENT_RECORD_VERSION);
    assert_eq!(first.mismatch_kind(), "byte");
    assert_eq!(first.first_offset(), offset);
    assert_eq!(first.implementation_byte(), Some(implementation[offset]));
    assert_eq!(first.reference_byte(), Some(reference[offset]));
    assert_ne!(first.implementation_digest(), first.reference_digest());
    let line = first.json_line();
    assert_eq!(line, replay.json_line());
    assert!(!line.contains('\n'));
    assert_eq!(
        line,
        "{\"casebook_disagreement\":1,\"suite\":\"casebook-disagreement-selftest\",\"case\":\"cb-disagreement-seeded-byte\",\"implementation\":\"frankensim\",\"reference\":\"frankenscipy\",\"mismatch_kind\":\"byte\",\"implementation_len\":34,\"reference_len\":34,\"implementation_digest\":\"6a4cbec3a30ceded\",\"reference_digest\":\"6b3fe07ef2feb8af\",\"first_offset\":11,\"implementation_byte\":\"71\",\"reference_byte\":\"73\"}"
    );

    let replay_spec = ReplaySpec::new(
        "cargo test -p fs-casebook --test casebook seeded_byte_corruption -- --exact",
        reference,
    );
    let report = Suite::new("casebook-disagreement-selftest")
        .case_replayable(
            "cb-disagreement-seeded-byte",
            replay_spec,
            ToleranceSpec::Exact,
            move || {
                CaseOutcome::pass(format!("seed=0x{SEED:016x}; disclosed corruption"))
                    .with_disagreement(first)
            },
        )
        .run();
    assert!(
        !report.all_passed(),
        "attaching a disagreement cannot leave a case green"
    );
    assert_eq!(report.replay_records.len(), 1);
    assert_eq!(report.disagreements.len(), 1);
    assert_eq!(
        report
            .disagreements_for("cb-disagreement-seeded-byte")
            .len(),
        1
    );

    let panic = std::panic::catch_unwind(|| report.assert_green())
        .expect_err("assert_green must carry replay and disagreement records");
    let message = panic
        .downcast_ref::<String>()
        .map(String::as_str)
        .or_else(|| panic.downcast_ref::<&str>().copied())
        .expect("Casebook refusal carries text");
    assert!(message.contains("\"casebook_replay\":1"));
    assert!(message.contains("\"casebook_disagreement\":1"));
    assert!(message.contains(&format!("seed=0x{SEED:016x}")));
}

#[test]
fn length_boundary_disagreement_localizes_the_first_absent_byte() {
    let implementation = [0x10, 0x20, 0x30];
    let reference = [0x10, 0x20, 0x30, 0x40];
    let disagreement = DisagreementRecord::first(
        "casebook-disagreement-selftest",
        "cb-disagreement-length",
        "short-implementation",
        "long-reference",
        &implementation,
        &reference,
    )
    .expect("different lengths disagree");
    assert_eq!(disagreement.mismatch_kind(), "length-boundary");
    assert_eq!(disagreement.first_offset(), 3);
    assert_eq!(disagreement.implementation_byte(), None);
    assert_eq!(disagreement.reference_byte(), Some(0x40));
    let line = disagreement.json_line();
    assert!(line.contains("\"mismatch_kind\":\"length-boundary\""));
    assert!(line.contains("\"implementation_byte\":null"));
    assert!(line.contains("\"reference_byte\":\"40\""));
}

#[test]
fn disagreement_identity_is_bound_to_its_owning_case() {
    let disagreement = DisagreementRecord::first(
        "foreign-suite",
        "foreign-case",
        "implementation",
        "reference",
        b"different",
        b"reference",
    )
    .expect("frames disagree");
    let report = Suite::new("owning-suite")
        .case("owning-case", 7, ToleranceSpec::Exact, move || {
            CaseOutcome::pass("identity must be rebound").with_disagreement(disagreement)
        })
        .run();

    assert!(!report.all_passed());
    let retained = report
        .disagreements_for("owning-case")
        .into_iter()
        .next()
        .expect("bound disagreement remains discoverable");
    assert_eq!(retained.suite(), "owning-suite");
    assert_eq!(retained.case(), "owning-case");
    assert!(report.disagreements_for("foreign-case").is_empty());

    let panic = std::panic::catch_unwind(|| report.assert_green())
        .expect_err("bound disagreement keeps its owner red");
    let message = panic
        .downcast_ref::<String>()
        .map(String::as_str)
        .or_else(|| panic.downcast_ref::<&str>().copied())
        .expect("Casebook refusal carries text");
    assert!(message.contains("\"suite\":\"owning-suite\""));
    assert!(message.contains("\"case\":\"owning-case\""));
    assert!(!message.contains("foreign-case"));
}

#[test]
fn replay_and_disagreement_json_escape_all_free_text() {
    let suite = "quote\"suite\nline";
    let case = "case\twith\\slash";
    let disagreement = DisagreementRecord::first(
        suite,
        case,
        "impl\"name\nline",
        "ref\\name\tvalue",
        b"a",
        b"b",
    )
    .expect("frames disagree");
    let report = Suite::new(suite)
        .case_replayable(
            case,
            ReplaySpec::new("runner --label \"quoted\"\nnext\targ\\path", b"a".to_vec()),
            ToleranceSpec::Exact,
            move || CaseOutcome::fail("escaped companion rows").with_disagreement(disagreement),
        )
        .run();

    let replay_line = report.replay_records[0].json_line();
    let disagreement_line = report.disagreements[0].json_line();
    assert!(!replay_line.contains('\n'));
    assert!(!disagreement_line.contains('\n'));
    assert!(replay_line.contains("quote\\\"suite\\nline"));
    assert!(replay_line.contains("case\\twith\\\\slash"));
    assert!(replay_line.contains("runner --label \\\"quoted\\\"\\nnext\\targ\\\\path"));
    assert!(disagreement_line.contains("impl\\\"name\\nline"));
    assert!(disagreement_line.contains("ref\\\\name\\tvalue"));
}
