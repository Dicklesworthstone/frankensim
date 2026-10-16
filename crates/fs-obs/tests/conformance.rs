//! Conformance suite for fs-obs.
//!
//! Cases registered here define the crate's cross-implementation contract
//! (plan §13.3). The shared conformance harness (contract-conformance-infra
//! bead) will supersede this hand-rolled runner; the case shape is designed
//! so that migration is additive.

#[test]
fn obs_001_conformance_verdict_self_hosts() {
    let version_matches = fs_obs::VERSION == env!("CARGO_PKG_VERSION");
    let mut probe_emitter = fs_obs::Emitter::new("fs-obs/conformance", "obs-001/probe");
    let probe = probe_emitter.emit(
        fs_obs::Severity::Info,
        fs_obs::EventKind::ConformanceCase {
            suite: "fs-obs/conformance".to_string(),
            case: "obs-001/probe".to_string(),
            pass: true,
            detail: "self-hosting probe".to_string(),
            seed: 0,
        },
        None,
    );
    let kind_matches = probe.kind.kind_name() == "conformance_case";
    let failure_lint_accepts = fs_obs::lint_failure_record(&probe).is_ok();
    let probe_line = probe.to_jsonl();
    let wire_validates = fs_obs::validate_line(&probe_line).is_ok();
    let receipt = probe.content_identity_receipt();
    let identity_readmits = probe.admit_content_identity(&receipt).is_ok();
    let pass = version_matches
        && kind_matches
        && failure_lint_accepts
        && wire_validates
        && identity_readmits;
    let detail = format!(
        "fs-obs version {} self-host check: version_matches={version_matches}, \
         kind_matches={kind_matches}, failure_lint_accepts={failure_lint_accepts}, \
         wire_validates={wire_validates}, identity_readmits={identity_readmits}",
        fs_obs::VERSION,
    );
    let mut emitter = fs_obs::Emitter::new("fs-obs/conformance", "obs-001/self-host");
    let event = emitter.emit(
        if pass {
            fs_obs::Severity::Info
        } else {
            fs_obs::Severity::Error
        },
        fs_obs::EventKind::ConformanceCase {
            suite: "fs-obs/conformance".to_string(),
            case: "obs-001".to_string(),
            pass,
            detail: detail.clone(),
            seed: 0,
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("self-hosted verdict must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("self-hosted verdict must use the canonical wire schema");
    println!("{line}");
    assert!(pass, "obs-001: {detail}");
}

#[test]
fn obs_002_gated_emission_preserves_public_contract() {
    use core::cell::Cell;
    use core::num::NonZeroU64;

    let gate = fs_obs::EmissionGate::new(
        fs_obs::SamplingCadence::every(
            NonZeroU64::new(2).expect("two is a non-zero trace cadence"),
        ),
        fs_obs::SamplingCadence::never(),
    );
    let builds = Cell::new(0_u64);
    let mut emitter = fs_obs::Emitter::new("fs-obs/conformance", "obs-002/gate");

    let rejected = emitter.emit_gated(gate, fs_obs::Severity::Trace, 1, || {
        panic!("a rejected public opportunity must not evaluate its builder")
    });
    assert!(rejected.is_none());
    assert_eq!(builds.get(), 0);

    let trace = emitter
        .emit_gated(gate, fs_obs::Severity::Trace, 2, || {
            builds.set(builds.get() + 1);
            (
                fs_obs::EventKind::TileComplete {
                    tile: 2,
                    kernel: "obs-002".to_string(),
                },
                Some(20),
            )
        })
        .expect("even trace opportunity must be emitted");
    let warn = emitter
        .emit_gated(gate, fs_obs::Severity::Warn, u64::MAX, || {
            builds.set(builds.get() + 1);
            (
                fs_obs::EventKind::Cancellation {
                    reason: "budget".to_string(),
                },
                Some(30),
            )
        })
        .expect("warn events must bypass sampling");
    let direct = emitter.emit(
        fs_obs::Severity::Info,
        fs_obs::EventKind::Custom {
            name: "direct".to_string(),
            json: "{}".to_string(),
        },
        None,
    );

    assert_eq!(builds.get(), 2);
    assert_eq!((trace.seq, warn.seq, direct.seq), (0, 1, 2));
    for event in [&trace, &warn, &direct] {
        fs_obs::validate_line(&event.to_jsonl())
            .expect("gated and direct events must share the canonical wire schema");
        event
            .admit_content_identity(&event.content_identity_receipt())
            .expect("gated and direct events must share exact identity semantics");
    }
}
