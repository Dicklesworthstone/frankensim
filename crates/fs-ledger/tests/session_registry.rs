//! Cross-crate conformance tests for durable session-mutation registry semantics.

use fs_ledger::hash::hash_bytes;
use fs_ledger::session_registry::{
    MAX_SESSION_CLAIM_PAYLOAD_BYTES, MAX_SESSION_FLUSH_ENCODED_BYTES, MAX_SESSION_FLUSH_EVENTS,
    MAX_SESSION_FLUSH_TERMINALS, MAX_SESSION_TERMINAL_RECEIPT_BYTES, SessionMutationClaim,
    SessionMutationClaimResult, SessionTerminalBatch, SessionTerminalBatchResult,
    SessionTerminalGroup, SessionTerminalRow,
};
use fs_ledger::{EventRow, Ledger, LedgerError};

fn authority(seed: u64) -> fs_ledger::ContentHash {
    hash_bytes(&seed.to_le_bytes())
}

fn claim<'a>(
    ledger: &Ledger,
    authority: fs_ledger::ContentHash,
    payload: &'a [u8],
) -> SessionMutationClaim<'a> {
    SessionMutationClaim {
        authority,
        ledger_instance_id: ledger.instance_id(),
        governor_hash: hash_bytes(b"governor"),
        session_open_hash: hash_bytes(b"open"),
        kind: "test-atomic",
        session: 7,
        ledger_scope: "registry-test",
        generation: 3,
        causal_ordinal: None,
        payload,
    }
}

fn submission_claim<'a>(
    ledger: &Ledger,
    authority: fs_ledger::ContentHash,
    admission_ordinal: u64,
    payload: &'a [u8],
) -> SessionMutationClaim<'a> {
    SessionMutationClaim {
        kind: "submission",
        causal_ordinal: Some(admission_ordinal),
        ..claim(ledger, authority, payload)
    }
}

#[test]
fn claim_pending_terminal_and_exact_batch_replay_are_append_once() {
    let ledger = Ledger::open(":memory:").expect("fixture ledger");
    let authority = authority(1);
    let claim = submission_claim(&ledger, authority, 11, b"program-v1");
    let permit = match ledger
        .claim_session_mutation(&claim)
        .expect("fresh durable claim")
    {
        SessionMutationClaimResult::Claimed { permit } => permit,
        other => panic!("fresh claim returned {other:?}"),
    };
    assert!(matches!(
        ledger
            .claim_session_mutation(&claim)
            .expect("exact pending replay"),
        SessionMutationClaimResult::Pending { .. }
    ));
    assert_eq!(
        ledger
            .pending_session_mutation(
                claim.governor_hash,
                claim.session_open_hash,
                claim.kind,
                claim.session,
                claim.ledger_scope,
                claim.generation,
            )
            .expect("bounded Pending probe")
            .expect("Pending claim is visible")
            .authority,
        authority
    );

    let session_bytes = 7_u64.to_be_bytes();
    let event = EventRow {
        session: Some(&session_bytes),
        t: 11,
        kind: "session.idempotent-execution",
        payload: Some(r#"{"schema":"registry-test-v1"}"#),
    };
    let events = [event];
    let group = SessionTerminalGroup {
        terminal: SessionTerminalRow {
            claim,
            permit: Some(permit),
            receipt: b"typed-terminal-v1",
        },
        events: &events,
    };
    let groups = [group];
    let batch = SessionTerminalBatch { groups: &groups };
    let batch_id = match ledger
        .append_session_terminal_batch(&batch)
        .expect("terminal commit")
    {
        SessionTerminalBatchResult::Committed { batch_id, .. } => batch_id,
        other => panic!("fresh terminal batch returned {other:?}"),
    };
    assert_eq!(ledger.table_count("session_claims").unwrap(), 1);
    assert_eq!(ledger.table_count("session_terminals").unwrap(), 1);
    assert_eq!(ledger.table_count("session_terminal_events").unwrap(), 1);
    assert_eq!(ledger.table_count("events").unwrap(), 1);
    assert_eq!(ledger.table_count("session_flush_batches").unwrap(), 1);
    assert_eq!(
        ledger
            .session_terminal(&authority)
            .expect("verified terminal")
            .expect("terminal present")
            .receipt,
        b"typed-terminal-v1"
    );

    assert_eq!(
        ledger
            .append_session_terminal_batch(&batch)
            .expect("exact batch replay"),
        SessionTerminalBatchResult::Replayed { batch_id }
    );
    assert_eq!(ledger.table_count("events").unwrap(), 1);
    assert!(matches!(
        ledger
            .claim_session_mutation(&claim)
            .expect("terminal claim replay"),
        SessionMutationClaimResult::Terminal { .. }
    ));
    assert!(
        ledger
            .pending_session_mutation(
                claim.governor_hash,
                claim.session_open_hash,
                claim.kind,
                claim.session,
                claim.ledger_scope,
                claim.generation,
            )
            .expect("terminalized claim probe")
            .is_none()
    );
}

#[test]
fn immutable_claim_authority_rejects_foreign_ledger_governor_session_and_scope() {
    let ledger = Ledger::open(":memory:").expect("fixture ledger");
    let foreign = Ledger::open(":memory:").expect("foreign ledger");
    let original = claim(&ledger, authority(90), b"bound-payload");
    assert!(matches!(
        ledger
            .claim_session_mutation(&original)
            .expect("seed immutable claim"),
        SessionMutationClaimResult::Claimed { .. }
    ));

    let foreign_ledger = SessionMutationClaim {
        ledger_instance_id: foreign.instance_id(),
        ..original
    };
    assert!(matches!(
        ledger.claim_session_mutation(&foreign_ledger),
        Err(LedgerError::Invalid { field, .. })
            if field == "session_claim.ledger_instance_id"
    ));

    for (field, altered) in [
        (
            "governor",
            SessionMutationClaim {
                governor_hash: hash_bytes(b"foreign-governor"),
                ..original
            },
        ),
        (
            "session-open",
            SessionMutationClaim {
                session_open_hash: hash_bytes(b"foreign-open"),
                ..original
            },
        ),
        (
            "session",
            SessionMutationClaim {
                session: original.session + 1,
                ..original
            },
        ),
        (
            "scope",
            SessionMutationClaim {
                ledger_scope: "foreign-scope",
                ..original
            },
        ),
        (
            "payload",
            SessionMutationClaim {
                payload: b"altered-payload",
                ..original
            },
        ),
    ] {
        let error = ledger
            .claim_session_mutation(&altered)
            .expect_err("immutable authority cannot be rebound");
        assert!(
            error.to_string().contains("different claim identity"),
            "{field} rebinding returned {error}"
        );
    }

    assert!(matches!(
        ledger
            .claim_session_mutation(&original)
            .expect("exact identity remains replayable"),
        SessionMutationClaimResult::Pending { .. }
    ));
    assert_eq!(ledger.table_count("session_claims").unwrap(), 1);
    assert_eq!(ledger.table_count("session_terminals").unwrap(), 0);
}

#[test]
fn mixed_retry_batch_adds_a_bounded_verified_witness_for_existing_terminal() {
    let ledger = Ledger::open(":memory:").expect("fixture ledger");
    let first = claim(&ledger, authority(101), b"first-payload");
    let first_groups = [SessionTerminalGroup {
        terminal: SessionTerminalRow {
            claim: first,
            permit: None,
            receipt: b"first-receipt",
        },
        events: &[],
    }];
    ledger
        .append_session_terminal_batch(&SessionTerminalBatch {
            groups: &first_groups,
        })
        .expect("first terminal batch");

    let second = claim(&ledger, authority(102), b"second-payload");
    let mixed_groups = [
        SessionTerminalGroup {
            terminal: SessionTerminalRow {
                claim: first,
                permit: None,
                receipt: b"first-receipt",
            },
            events: &[],
        },
        SessionTerminalGroup {
            terminal: SessionTerminalRow {
                claim: second,
                permit: None,
                receipt: b"second-receipt",
            },
            events: &[],
        },
    ];
    let mixed = SessionTerminalBatch {
        groups: &mixed_groups,
    };
    assert!(matches!(
        ledger
            .append_session_terminal_batch(&mixed)
            .expect("mixed retry commits only the new terminal"),
        SessionTerminalBatchResult::Committed {
            terminals_inserted: 1,
            events_appended: 0,
            ..
        }
    ));
    assert_eq!(ledger.table_count("session_claims").unwrap(), 2);
    assert_eq!(ledger.table_count("session_terminals").unwrap(), 2);
    assert_eq!(ledger.table_count("session_flush_batches").unwrap(), 2);
    assert_eq!(
        ledger.table_count("session_flush_batch_members").unwrap(),
        3
    );
    assert_eq!(
        ledger
            .session_terminal(&first.authority)
            .expect("all first-terminal batch witnesses verify")
            .expect("first terminal remains readable")
            .receipt,
        b"first-receipt"
    );
    assert!(matches!(
        ledger
            .append_session_terminal_batch(&mixed)
            .expect("exact mixed-batch replay"),
        SessionTerminalBatchResult::Replayed { .. }
    ));
}

#[test]
fn mixed_batch_conflict_rolls_back_earlier_auto_claim_and_events() {
    let ledger = Ledger::open(":memory:").expect("fixture ledger");
    let conflicting_authority = authority(2);
    let original = claim(&ledger, conflicting_authority, b"original");
    assert!(matches!(
        ledger
            .claim_session_mutation(&original)
            .expect("seed pending claim"),
        SessionMutationClaimResult::Claimed { .. }
    ));

    let first_claim = claim(&ledger, authority(3), b"first");
    let altered_claim = claim(&ledger, conflicting_authority, b"altered");
    let session_bytes = 7_u64.to_be_bytes();
    let first_event = EventRow {
        session: Some(&session_bytes),
        t: 1,
        kind: "session.open",
        payload: Some(r#"{"schema":"first"}"#),
    };
    let first_events = [first_event];
    let groups = [
        SessionTerminalGroup {
            terminal: SessionTerminalRow {
                claim: first_claim,
                permit: None,
                receipt: b"first-terminal",
            },
            events: &first_events,
        },
        SessionTerminalGroup {
            terminal: SessionTerminalRow {
                claim: altered_claim,
                permit: None,
                receipt: b"conflicting-terminal",
            },
            events: &[],
        },
    ];
    let error = ledger
        .append_session_terminal_batch(&SessionTerminalBatch { groups: &groups })
        .expect_err("later authority conflict rolls back the whole batch");
    assert!(error.to_string().contains("different claim identity"));
    assert!(
        ledger
            .session_mutation_claim(&first_claim.authority)
            .unwrap()
            .is_none()
    );
    assert!(
        ledger
            .session_terminal(&first_claim.authority)
            .unwrap()
            .is_none()
    );
    assert_eq!(ledger.table_count("events").unwrap(), 0);
    assert_eq!(ledger.table_count("session_terminals").unwrap(), 0);
    assert_eq!(ledger.table_count("session_claims").unwrap(), 1);
}

#[test]
fn pending_claim_without_positive_permit_cannot_be_terminalized() {
    let ledger = Ledger::open(":memory:").expect("fixture ledger");
    let authority = authority(4);
    let claim = submission_claim(&ledger, authority, 12, b"may-have-run");
    let _permit = match ledger.claim_session_mutation(&claim).expect("fresh claim") {
        SessionMutationClaimResult::Claimed { permit } => permit,
        other => panic!("fresh claim returned {other:?}"),
    };
    let group = SessionTerminalGroup {
        terminal: SessionTerminalRow {
            claim,
            permit: None,
            receipt: b"invented-terminal",
        },
        events: &[],
    };
    let groups = [group];
    let error = ledger
        .append_session_terminal_batch(&SessionTerminalBatch { groups: &groups })
        .expect_err("pending claim without permit is indeterminate");
    assert!(error.to_string().contains("Indeterminate"));
    assert!(ledger.session_terminal(&authority).unwrap().is_none());
}

#[test]
fn submission_requires_preclaim_and_admission_ordinals_are_unique() {
    let ledger = Ledger::open(":memory:").expect("fixture ledger");
    let fabricated = submission_claim(&ledger, authority(5), 20, b"fabricated");
    let fabricated_groups = [SessionTerminalGroup {
        terminal: SessionTerminalRow {
            claim: fabricated,
            permit: None,
            receipt: b"fabricated-terminal",
        },
        events: &[],
    }];
    let error = ledger
        .append_session_terminal_batch(&SessionTerminalBatch {
            groups: &fabricated_groups,
        })
        .expect_err("submission terminal cannot bypass preclaim");
    assert!(error.to_string().contains("pre-execution claim"));
    assert!(
        ledger
            .session_mutation_claim(&fabricated.authority)
            .unwrap()
            .is_none()
    );

    let first = submission_claim(&ledger, authority(6), 21, b"first");
    assert!(matches!(
        ledger.claim_session_mutation(&first).expect("first claim"),
        SessionMutationClaimResult::Claimed { .. }
    ));
    let collision = submission_claim(&ledger, authority(7), 21, b"collision");
    let error = ledger
        .claim_session_mutation(&collision)
        .expect_err("one governor/kind ordinal has one owner");
    assert!(error.to_string().contains("already owned"));
    assert!(
        ledger
            .session_mutation_claim(&collision.authority)
            .unwrap()
            .is_none()
    );
    assert_eq!(
        ledger
            .session_mutation_claim_count(first.governor_hash)
            .expect("bounded governor count"),
        1
    );
}

#[test]
fn pause_and_submission_claims_fence_each_other_transactionally() {
    let ledger = Ledger::open(":memory:").expect("joint terminal ledger");
    let submission = SessionMutationClaim {
        generation: 3,
        ..submission_claim(&ledger, authority(8), 22, b"work")
    };
    let permit = match ledger
        .claim_session_mutation(&submission)
        .expect("Pending submission")
    {
        SessionMutationClaimResult::Claimed { permit } => permit,
        other => panic!("fresh claim returned {other:?}"),
    };
    let pause = SessionMutationClaim {
        authority: authority(9),
        kind: "pause-acknowledgement",
        generation: 4,
        causal_ordinal: Some(23),
        payload: b"checkpoint",
        ..claim(&ledger, authority(9), b"checkpoint")
    };
    let pause_only = [SessionTerminalGroup {
        terminal: SessionTerminalRow {
            claim: pause,
            permit: None,
            receipt: b"pause-terminal",
        },
        events: &[],
    }];
    let error = ledger
        .append_session_terminal_batch(&SessionTerminalBatch {
            groups: &pause_only,
        })
        .expect_err("Pending draining-generation work fences pause terminalization");
    assert!(error.to_string().contains("durably Pending"));

    let joint = [
        SessionTerminalGroup {
            terminal: SessionTerminalRow {
                claim: pause,
                permit: None,
                receipt: b"pause-terminal",
            },
            events: &[],
        },
        SessionTerminalGroup {
            terminal: SessionTerminalRow {
                claim: submission,
                permit: Some(permit),
                receipt: b"submission-terminal",
            },
            events: &[],
        },
    ];
    ledger
        .append_session_terminal_batch(&SessionTerminalBatch { groups: &joint })
        .expect("one transaction terminalizes work and pause together");

    let fenced = Ledger::open(":memory:").expect("successor fence ledger");
    let pause = SessionMutationClaim {
        authority: authority(10),
        kind: "pause-acknowledgement",
        generation: 4,
        causal_ordinal: Some(24),
        payload: b"checkpoint",
        ..claim(&fenced, authority(10), b"checkpoint")
    };
    let pause_groups = [SessionTerminalGroup {
        terminal: SessionTerminalRow {
            claim: pause,
            permit: None,
            receipt: b"pause-terminal",
        },
        events: &[],
    }];
    fenced
        .append_session_terminal_batch(&SessionTerminalBatch {
            groups: &pause_groups,
        })
        .expect("terminal pause fence");
    let stale = SessionMutationClaim {
        generation: 3,
        ..submission_claim(&fenced, authority(11), 25, b"late-work")
    };
    let error = fenced
        .claim_session_mutation(&stale)
        .expect_err("terminal successor pause rejects late old-generation work");
    assert!(error.to_string().contains("already fenced"));
}

#[test]
fn terminal_and_owned_events_verify_after_real_file_reopen() {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .expect("wall clock after epoch")
        .as_nanos();
    let path = std::env::temp_dir().join(format!(
        "fs-ledger-session-registry-{}-{unique}.db",
        std::process::id()
    ));
    let path = path.to_str().expect("UTF-8 temporary path").to_string();
    let authority = authority(5);
    let instance = {
        let ledger = Ledger::open(&path).expect("file ledger");
        let claim = claim(&ledger, authority, b"reopen-payload");
        let session_bytes = 7_u64.to_be_bytes();
        let event = EventRow {
            session: Some(&session_bytes),
            t: 12,
            kind: "session.meter-report",
            payload: Some(r#"{"schema":"reopen-v1"}"#),
        };
        let events = [event];
        let group = SessionTerminalGroup {
            terminal: SessionTerminalRow {
                claim,
                permit: None,
                receipt: b"reopen-receipt",
            },
            events: &events,
        };
        let groups = [group];
        ledger
            .append_session_terminal_batch(&SessionTerminalBatch { groups: &groups })
            .expect("atomic file commit");
        ledger.instance_id()
    };

    let reopened = Ledger::open(&path).expect("reopened file ledger");
    assert_eq!(reopened.instance_id(), instance);
    let terminal = reopened
        .session_terminal(&authority)
        .expect("bounded verified lookup")
        .expect("terminal survived reopen");
    assert_eq!(terminal.receipt, b"reopen-receipt");
    assert_eq!(terminal.event_count, 1);
    assert_eq!(reopened.table_count("events").unwrap(), 1);
}

#[test]
#[allow(clippy::too_many_lines)] // One exact-cap/limit+1 matrix shares fixtures and atomicity checks.
fn exact_terminal_count_and_field_caps_pass_while_limit_plus_one_is_atomic() {
    let ledger = Ledger::open(":memory:").expect("fixture ledger");
    let claims: Vec<_> = (0..MAX_SESSION_FLUSH_TERMINALS)
        .map(|index| claim(&ledger, authority(u64::try_from(index).unwrap() + 100), b""))
        .collect();
    let groups: Vec<_> = claims
        .iter()
        .map(|claim| SessionTerminalGroup {
            terminal: SessionTerminalRow {
                claim: *claim,
                permit: None,
                receipt: b"x",
            },
            events: &[],
        })
        .collect();
    ledger
        .append_session_terminal_batch(&SessionTerminalBatch { groups: &groups })
        .expect("exact terminal-count cap");
    assert_eq!(
        ledger.table_count("session_terminals").unwrap(),
        u64::try_from(MAX_SESSION_FLUSH_TERMINALS).unwrap()
    );

    let overflow_ledger = Ledger::open(":memory:").expect("overflow ledger");
    let overflow_claims: Vec<_> = (0..=MAX_SESSION_FLUSH_TERMINALS)
        .map(|index| {
            claim(
                &overflow_ledger,
                authority(u64::try_from(index).unwrap() + 10_000),
                b"",
            )
        })
        .collect();
    let overflow_groups: Vec<_> = overflow_claims
        .iter()
        .map(|claim| SessionTerminalGroup {
            terminal: SessionTerminalRow {
                claim: *claim,
                permit: None,
                receipt: b"x",
            },
            events: &[],
        })
        .collect();
    overflow_ledger
        .append_session_terminal_batch(&SessionTerminalBatch {
            groups: &overflow_groups,
        })
        .expect_err("terminal-count cap plus one");
    assert_eq!(overflow_ledger.table_count("session_claims").unwrap(), 0);

    let exact_payload = vec![0_u8; MAX_SESSION_CLAIM_PAYLOAD_BYTES];
    let payload_claim = claim(&overflow_ledger, authority(30_000), &exact_payload);
    assert!(matches!(
        overflow_ledger
            .claim_session_mutation(&payload_claim)
            .expect("exact payload cap"),
        SessionMutationClaimResult::Claimed { .. }
    ));
    let oversized_payload = vec![0_u8; MAX_SESSION_CLAIM_PAYLOAD_BYTES + 1];
    let oversized_claim = claim(&overflow_ledger, authority(30_001), &oversized_payload);
    overflow_ledger
        .claim_session_mutation(&oversized_claim)
        .expect_err("payload cap plus one");
    assert!(
        overflow_ledger
            .session_mutation_claim(&oversized_claim.authority)
            .unwrap()
            .is_none()
    );

    let exact_receipt = vec![0_u8; MAX_SESSION_TERMINAL_RECEIPT_BYTES];
    let receipt_claim = claim(&overflow_ledger, authority(30_002), b"");
    let receipt_group = SessionTerminalGroup {
        terminal: SessionTerminalRow {
            claim: receipt_claim,
            permit: None,
            receipt: &exact_receipt,
        },
        events: &[],
    };
    let receipt_groups = [receipt_group];
    overflow_ledger
        .append_session_terminal_batch(&SessionTerminalBatch {
            groups: &receipt_groups,
        })
        .expect("exact receipt cap");
    let oversized_receipt = vec![0_u8; MAX_SESSION_TERMINAL_RECEIPT_BYTES + 1];
    let oversized_receipt_claim = claim(&overflow_ledger, authority(30_003), b"");
    let oversized_receipt_group = SessionTerminalGroup {
        terminal: SessionTerminalRow {
            claim: oversized_receipt_claim,
            permit: None,
            receipt: &oversized_receipt,
        },
        events: &[],
    };
    let oversized_receipt_groups = [oversized_receipt_group];
    overflow_ledger
        .append_session_terminal_batch(&SessionTerminalBatch {
            groups: &oversized_receipt_groups,
        })
        .expect_err("receipt cap plus one");
    assert!(
        overflow_ledger
            .session_mutation_claim(&oversized_receipt_claim.authority)
            .unwrap()
            .is_none()
    );
}

fn exact_json_string(bytes: usize) -> String {
    assert!(bytes >= 2);
    format!("\"{}\"", "x".repeat(bytes - 2))
}

#[test]
#[allow(clippy::too_many_lines)] // One exact/limit+1 matrix shares the large bounded fixtures.
fn exact_event_count_and_aggregate_byte_caps_pass_and_limit_plus_one_is_atomic() {
    let count_ledger = Ledger::open(":memory:").expect("event-count ledger");
    let count_authority = authority(40_000);
    let count_claim = claim(&count_ledger, count_authority, b"");
    let session = 7_u64.to_be_bytes();
    let events: Vec<_> = (0..MAX_SESSION_FLUSH_EVENTS)
        .map(|index| EventRow {
            session: Some(session.as_slice()),
            t: i64::try_from(index).expect("bounded event ordinal"),
            kind: "e",
            payload: Some("{}"),
        })
        .collect();
    let groups = [SessionTerminalGroup {
        terminal: SessionTerminalRow {
            claim: count_claim,
            permit: None,
            receipt: b"x",
        },
        events: &events,
    }];
    count_ledger
        .append_session_terminal_batch(&SessionTerminalBatch { groups: &groups })
        .expect("exact event-count cap");
    assert_eq!(
        count_ledger
            .session_terminal(&count_authority)
            .unwrap()
            .unwrap()
            .event_count,
        MAX_SESSION_FLUSH_EVENTS
    );

    let overflow_count = Ledger::open(":memory:").expect("event-count overflow ledger");
    let overflow_authority = authority(40_001);
    let overflow_claim = claim(&overflow_count, overflow_authority, b"");
    let overflow_events: Vec<_> = (0..=MAX_SESSION_FLUSH_EVENTS)
        .map(|index| EventRow {
            session: Some(session.as_slice()),
            t: i64::try_from(index).expect("bounded event ordinal"),
            kind: "e",
            payload: Some("{}"),
        })
        .collect();
    let overflow_groups = [SessionTerminalGroup {
        terminal: SessionTerminalRow {
            claim: overflow_claim,
            permit: None,
            receipt: b"x",
        },
        events: &overflow_events,
    }];
    overflow_count
        .append_session_terminal_batch(&SessionTerminalBatch {
            groups: &overflow_groups,
        })
        .expect_err("event-count limit plus one");
    assert_eq!(overflow_count.table_count("session_claims").unwrap(), 0);
    assert_eq!(overflow_count.table_count("events").unwrap(), 0);

    // Conservative exact preimage: claim framing + actual kind/scope bytes +
    // payload, terminal framing + receipt, then two event framings of
    // 64 + 8 session bytes + one kind byte + JSON payload.
    let exact_ledger = Ledger::open(":memory:").expect("exact-byte ledger");
    let exact_payload = vec![0_u8; MAX_SESSION_CLAIM_PAYLOAD_BYTES];
    let exact_receipt = vec![0_u8; MAX_SESSION_TERMINAL_RECEIPT_BYTES];
    let exact_claim = claim(&exact_ledger, authority(40_002), &exact_payload);
    let fixed_bytes = 256
        + exact_claim.kind.len()
        + exact_claim.ledger_scope.len()
        + exact_payload.len()
        + 96
        + exact_receipt.len()
        + 2 * (64 + session.len() + "e".len());
    let remaining = MAX_SESSION_FLUSH_ENCODED_BYTES - fixed_bytes;
    let first_json_len = remaining.min(1024 * 1024);
    let second_json_len = remaining - first_json_len;
    let first_json = exact_json_string(first_json_len);
    let second_json = exact_json_string(second_json_len);
    let exact_events = [
        EventRow {
            session: Some(session.as_slice()),
            t: 1,
            kind: "e",
            payload: Some(first_json.as_str()),
        },
        EventRow {
            session: Some(session.as_slice()),
            t: 2,
            kind: "e",
            payload: Some(second_json.as_str()),
        },
    ];
    let exact_groups = [SessionTerminalGroup {
        terminal: SessionTerminalRow {
            claim: exact_claim,
            permit: None,
            receipt: &exact_receipt,
        },
        events: &exact_events,
    }];
    let exact_result = exact_ledger
        .append_session_terminal_batch(&SessionTerminalBatch {
            groups: &exact_groups,
        })
        .expect("exact aggregate-byte cap");
    assert!(matches!(
        exact_result,
        SessionTerminalBatchResult::Committed {
            events_appended: 2,
            ..
        }
    ));

    let overflow_bytes = Ledger::open(":memory:").expect("byte-overflow ledger");
    let overflow_payload = exact_payload;
    let overflow_receipt = exact_receipt;
    let overflow_second = exact_json_string(second_json_len + 1);
    let overflow_events = [
        EventRow {
            session: Some(session.as_slice()),
            t: 1,
            kind: "e",
            payload: Some(first_json.as_str()),
        },
        EventRow {
            session: Some(session.as_slice()),
            t: 2,
            kind: "e",
            payload: Some(overflow_second.as_str()),
        },
    ];
    let overflow_claim = claim(&overflow_bytes, authority(40_003), &overflow_payload);
    let overflow_groups = [SessionTerminalGroup {
        terminal: SessionTerminalRow {
            claim: overflow_claim,
            permit: None,
            receipt: &overflow_receipt,
        },
        events: &overflow_events,
    }];
    overflow_bytes
        .append_session_terminal_batch(&SessionTerminalBatch {
            groups: &overflow_groups,
        })
        .expect_err("aggregate-byte cap plus one");
    assert_eq!(overflow_bytes.table_count("session_claims").unwrap(), 0);
    assert_eq!(overflow_bytes.table_count("events").unwrap(), 0);
}
