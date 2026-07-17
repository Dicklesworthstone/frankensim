//! fs-casebook self-conformance (bead huq.5): the demo suite every agent
//! copies to add cases, plus the fail-closed self-tests the bead's
//! acceptance demands (an intentionally failing case must yield a
//! structured failure record and a non-green report).

use fs_casebook::{CaseOutcome, CaseRecord, Suite, ToleranceSpec, fnv1a64};

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
}

#[test]
fn digest_helper_is_the_canonical_fnv1a64() {
    assert_eq!(fnv1a64(b""), 0xcbf2_9ce4_8422_2325);
    assert_ne!(fnv1a64(b"a"), fnv1a64(b"b"));
    assert_eq!(fnv1a64(b"casebook"), fnv1a64(b"casebook"));
}
