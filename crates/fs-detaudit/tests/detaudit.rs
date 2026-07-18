//! fs-detaudit engine self-conformance (bead 6nb.6 testing spec): seeded
//! violations per category — an arrival-order reduction is CAUGHT with
//! correct localization, cross-ISA fixtures classify with no silent rows,
//! envelope violations refuse, and the mode-delta tradeoff is observed.

use fs_detaudit::{
    AuditConfig, Classification, DivergenceClass, DivergencePolicy, IsaLedger, LedgerRow,
    StageHash, StagedTrace, Subject, WorkerMatrix, audit, classify_cross_isa, fnv1a64,
    measure_mode_delta,
};
use std::collections::BTreeMap;

fn stage(label: &str, hash: u64) -> StageHash {
    StageHash {
        label: label.to_owned(),
        hash,
    }
}

/// A deterministic three-stage subject: identical trace at any worker
/// count.
fn deterministic_subject() -> Subject {
    Subject {
        name: "selftest-deterministic",
        run: Box::new(|_workers| {
            let inputs: Vec<f64> = (1..=64).map(|i| 1.0 / f64::from(i)).collect();
            let input_bytes: Vec<u8> = inputs
                .iter()
                .flat_map(|x| x.to_bits().to_le_bytes())
                .collect();
            // Fixed pairwise order regardless of workers.
            let sum: f64 = inputs.iter().sum();
            StagedTrace {
                stages: vec![
                    stage("load", fnv1a64(&input_bytes)),
                    stage("reduce", fnv1a64(&sum.to_bits().to_le_bytes())),
                    stage("publish", fnv1a64(b"published")),
                ],
            }
        }),
    }
}

/// The bead's seeded violation: an arrival-order-dependent float
/// reduction — summation order depends on the worker count, so the
/// "reduce" stage hash drifts while "load" stays identical.
fn arrival_order_subject() -> Subject {
    Subject {
        name: "selftest-arrival-order",
        run: Box::new(|workers| {
            let inputs: Vec<f64> = (1..=64).map(|i| 1.0 / f64::from(i)).collect();
            let input_bytes: Vec<u8> = inputs
                .iter()
                .flat_map(|x| x.to_bits().to_le_bytes())
                .collect();
            // Simulated arrival order: stride the inputs by the worker
            // count (what a non-deterministic reduction tree does).
            let mut sum = 0.0f64;
            let stride = workers.max(1);
            for lane in 0..stride {
                let mut i = lane;
                while i < inputs.len() {
                    sum += inputs[i];
                    i += stride;
                }
            }
            StagedTrace {
                stages: vec![
                    stage("load", fnv1a64(&input_bytes)),
                    stage("reduce", fnv1a64(&sum.to_bits().to_le_bytes())),
                    stage("publish", fnv1a64(b"published")),
                ],
            }
        }),
    }
}

#[test]
fn deterministic_subject_audits_identical_across_the_matrix() {
    let config = AuditConfig {
        matrix: WorkerMatrix::explicit(vec![1, 2, 4, 6, 16]),
        repeats: 3,
    };
    let report = audit(&deterministic_subject(), &config);
    assert!(report.identical(), "deterministic subject must audit clean");
    assert_eq!(report.rows.len(), 5 * 3);
    assert!(report.rows.iter().all(|r| r.identical));
    for line in report.json_lines() {
        assert!(line.starts_with("{\"detaudit\":\"row\""));
        assert!(line.contains("\"identical\":true"));
    }
}

#[test]
fn audit_json_lines_escape_subjects_and_stage_labels() {
    const SUBJECT: &str = "probe\"\\name\nsecond\u{0001}";
    const STAGE: &str = "stage\"\\label\r\t\u{0008}\u{000c}\u{0002}";
    let subject = Subject {
        name: SUBJECT,
        run: Box::new(|workers| StagedTrace {
            stages: vec![stage(
                STAGE,
                u64::try_from(workers).expect("worker count fits u64"),
            )],
        }),
    };
    let report = audit(
        &subject,
        &AuditConfig {
            matrix: WorkerMatrix::explicit(vec![1, 2]),
            repeats: 1,
        },
    );
    assert_eq!(report.divergences.len(), 1);
    let lines = report.json_lines();
    for line in &lines {
        assert_eq!(line.lines().count(), 1, "{line}");
        assert!(!line.chars().any(|ch| ch < ' '), "{line}");
        assert!(
            line.contains("\"subject\":\"probe\\\"\\\\name\\nsecond\\u0001\""),
            "{line}"
        );
    }
    let divergence = lines
        .iter()
        .find(|line| line.contains("\"detaudit\":\"divergence\""))
        .expect("one divergence record");
    assert!(
        divergence.contains("\"stage_label\":\"stage\\\"\\\\label\\r\\t\\b\\f\\u0002\""),
        "{divergence}"
    );
}

#[test]
fn audit_json_lines_exhaustively_escape_the_json_c0_range() {
    const RAW: &str = "\u{0000}\u{0001}\u{0002}\u{0003}\u{0004}\u{0005}\u{0006}\u{0007}\u{0008}\u{0009}\u{000a}\u{000b}\u{000c}\u{000d}\u{000e}\u{000f}\u{0010}\u{0011}\u{0012}\u{0013}\u{0014}\u{0015}\u{0016}\u{0017}\u{0018}\u{0019}\u{001a}\u{001b}\u{001c}\u{001d}\u{001e}\u{001f}\"\\";
    const ESCAPED: &str = "\\u0000\\u0001\\u0002\\u0003\\u0004\\u0005\\u0006\\u0007\\b\\t\\n\\u000b\\f\\r\\u000e\\u000f\\u0010\\u0011\\u0012\\u0013\\u0014\\u0015\\u0016\\u0017\\u0018\\u0019\\u001a\\u001b\\u001c\\u001d\\u001e\\u001f\\\"\\\\";

    let subject = Subject {
        name: RAW,
        run: Box::new(|workers| StagedTrace {
            stages: vec![stage(
                RAW,
                u64::try_from(workers).expect("worker count fits u64"),
            )],
        }),
    };
    let report = audit(
        &subject,
        &AuditConfig {
            matrix: WorkerMatrix::explicit(vec![1, 2]),
            repeats: 1,
        },
    );

    assert_eq!(
        report.json_lines(),
        vec![
            format!(
                "{{\"detaudit\":\"row\",\"subject\":\"{ESCAPED}\",\"workers\":1,\"repeat\":0,\"final_hash\":\"0000000000000001\",\"identical\":true}}"
            ),
            format!(
                "{{\"detaudit\":\"row\",\"subject\":\"{ESCAPED}\",\"workers\":2,\"repeat\":0,\"final_hash\":\"0000000000000002\",\"identical\":false}}"
            ),
            format!(
                "{{\"detaudit\":\"divergence\",\"subject\":\"{ESCAPED}\",\"workers\":2,\"repeat\":0,\"first_stage\":0,\"stage_label\":\"{ESCAPED}\",\"baseline\":\"0000000000000001\",\"observed\":\"0000000000000002\"}}"
            ),
        ]
    );
}

#[test]
fn seeded_arrival_order_violation_is_caught_and_localized_to_the_reduce_stage() {
    let config = AuditConfig {
        matrix: WorkerMatrix::explicit(vec![1, 2, 4]),
        repeats: 2,
    };
    let report = audit(&arrival_order_subject(), &config);
    assert!(!report.identical(), "arrival-order drift must be caught");
    assert!(!report.divergences.is_empty());
    for divergence in &report.divergences {
        assert_eq!(
            divergence.first_stage, 1,
            "the locator must name the reduce stage, not load or publish"
        );
        assert_eq!(divergence.stage_label, "reduce");
        assert_ne!(divergence.baseline_hash, divergence.observed_hash);
        assert!(
            divergence.workers > 1,
            "the 1-worker baseline agrees with itself"
        );
    }
    // Repeats at the SAME worker count stay self-consistent (the drift is
    // worker-count-dependent, not run-to-run noise): every diverging
    // worker count diverges on both repeats, with identical bits.
    // (Not every reordering rounds differently — stride 2 over this
    // fixture happens to round identically — so assert over whichever
    // counts actually diverged.)
    let mut by_workers: std::collections::BTreeMap<usize, Vec<u64>> =
        std::collections::BTreeMap::new();
    for d in &report.divergences {
        by_workers
            .entry(d.workers)
            .or_default()
            .push(d.observed_hash);
    }
    assert!(!by_workers.is_empty());
    for (workers, hashes) in by_workers {
        assert_eq!(hashes.len(), 2, "both repeats at workers={workers} diverge");
        assert_eq!(
            hashes[0], hashes[1],
            "repeat-consistent at workers={workers}"
        );
    }
}

fn row(hash: u64, value: Option<f64>) -> LedgerRow {
    LedgerRow {
        hash,
        value_bits: value.map(f64::to_bits),
    }
}

fn ledger(isa: &str, rows: &[(&str, LedgerRow)]) -> IsaLedger {
    IsaLedger {
        isa: isa.to_owned(),
        rows: rows
            .iter()
            .map(|(name, r)| ((*name).to_owned(), r.clone()))
            .collect(),
    }
}

#[test]
fn cross_isa_report_classifies_every_declared_category() {
    let base = 1.000_000_000_000_1_f64;
    let one_ulp_up = f64::from_bits(base.to_bits() + 1);
    let a = ledger(
        "aarch64",
        &[
            ("kernel/identical", row(0xAAAA, None)),
            ("kernel/fma", row(0xB001, None)),
            (
                "libm/atan2",
                row(fnv1a64(&base.to_bits().to_le_bytes()), Some(base)),
            ),
        ],
    );
    let b = ledger(
        "x86-64",
        &[
            ("kernel/identical", row(0xAAAA, None)),
            ("kernel/fma", row(0xB002, None)),
            (
                "libm/atan2",
                row(
                    fnv1a64(&one_ulp_up.to_bits().to_le_bytes()),
                    Some(one_ulp_up),
                ),
            ),
        ],
    );
    let mut policy = DivergencePolicy::default();
    policy
        .declared
        .insert("kernel/fma".to_owned(), DivergenceClass::FmaContraction);
    policy.declared.insert(
        "libm/atan2".to_owned(),
        DivergenceClass::LibmUlp { max_ulps: 4 },
    );

    let report = classify_cross_isa(&a, &b, &policy);
    assert!(
        report.clean(),
        "every divergence is declared and in-envelope"
    );
    assert_eq!(report.rows.len(), 3);
    assert!(matches!(
        report.rows[1].classification,
        Classification::Identical
    ));
    let rendered = report.render_markdown();
    assert!(rendered.contains("fma-contraction"));
    assert!(rendered.contains("libm-ulp<=4: measured 1 ulps within envelope"));
    assert!(!rendered.contains("UNCLASSIFIED"));
    assert!(rendered.contains("clean (every divergence classified)"));
}

#[test]
fn undeclared_divergence_and_violated_envelopes_are_never_clean() {
    let a = ledger("aarch64", &[("kernel/reduction", row(0xC001, None))]);
    let b = ledger("x86-64", &[("kernel/reduction", row(0xC002, None))]);
    let report = classify_cross_isa(&a, &b, &DivergencePolicy::default());
    assert!(
        !report.clean(),
        "undeclared divergence must fail the report"
    );
    assert!(report.render_markdown().contains("UNCLASSIFIED"));

    let far = 2.0f64;
    let near = f64::from_bits(1.0f64.to_bits() + 3);
    let a = ledger("aarch64", &[("libm/exp", row(0x1, Some(1.0)))]);
    let b_violates = ledger("x86-64", &[("libm/exp", row(0x2, Some(far)))]);
    let b_within = ledger("x86-64", &[("libm/exp", row(0x2, Some(near)))]);
    let mut policy = DivergencePolicy::default();
    policy.declared.insert(
        "libm/exp".to_owned(),
        DivergenceClass::LibmUlp { max_ulps: 3 },
    );
    assert!(
        !classify_cross_isa(&a, &b_violates, &policy).clean(),
        "an envelope violation is a build failure, not a classification"
    );
    assert!(classify_cross_isa(&a, &b_within, &policy).clean());

    let missing = IsaLedger {
        isa: "x86-64".to_owned(),
        rows: BTreeMap::new(),
    };
    assert!(
        !classify_cross_isa(&a, &missing, &policy).clean(),
        "a missing artifact row is never silently fine"
    );
}

#[test]
fn mode_delta_records_the_reproducibility_loss_honestly() {
    use std::cell::Cell;
    let counter = Cell::new(0u64);
    let report = measure_mode_delta(
        "selftest-mode-delta",
        4,
        || fnv1a64(b"deterministic artifact"),
        move || {
            // Fast mode trades reproducibility away: each run hashes
            // differently (arrival-order surrogate).
            counter.set(counter.get() + 1);
            fnv1a64(&counter.get().to_le_bytes())
        },
    );
    assert!(report.deterministic_reproducible);
    assert!(
        !report.fast_reproducible,
        "the reproducibility loss must be OBSERVED, not assumed"
    );
    assert!(report.gain_ratio.is_finite() && report.gain_ratio > 0.0);
    assert!(report.deterministic_ns > 0 && report.fast_ns > 0);
}
