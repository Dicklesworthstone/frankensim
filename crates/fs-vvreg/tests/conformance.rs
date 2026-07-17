//! Conformance battery for the G1/G2 benchmark & V&V registry (bead
//! frankensim-ext-benchmark-vv-registry-f1gq).
//!
//! Gauntlet coverage: G0 admission/framing/composition laws (fail-closed
//! citation gates, sealed receipts, mutation-sensitive identity,
//! deterministic serialization), plus the two fixtures the bead names
//! explicitly: the registry lint (an entry missing any load-bearing field
//! cannot be cited, with a typed refusal) and the CI proof that an
//! unpinned family-name citation fails.

use fs_vvreg::{
    AcceptanceEnvelope, CitationRefusal, ColorRank, ConsumptionRecord, ConsumptionRefusal,
    ConsumptionStatus, DeckPin, ENVELOPE_VERDICT_SCHEMA_VERSION, Edition, EnvelopeEvaluation,
    EnvelopeGateError, EnvelopeObservation, IntegrityFinding, LicenseState, MAX_BEAD_ID_LEN,
    MAX_LOOKUP_ID_LEN, MAX_QOIS_PER_ENTRY, OracleBinding, Qoi, Registry, RegistryEntry,
    RegistryTier, VVREG_VERSION, canonical_row, entry_digest, registry, validate_entry,
};

/// Distinct QoIs for the count/order mutation locks.
const QOI_A: Qoi = Qoi {
    name: "qoi_a",
    unit: "1",
    envelope: AcceptanceEnvelope::Tolerance {
        atol: 0.0,
        rtol: 1e-12,
    },
};

/// Second distinct QoI for the count/order mutation locks.
const QOI_B: Qoi = Qoi {
    name: "qoi_b",
    unit: "1",
    envelope: AcceptanceEnvelope::Tolerance {
        atol: 0.0,
        rtol: 1e-12,
    },
};

/// A fully pinned synthetic entry for gate-order probing.
const fn pinned_probe() -> RegistryEntry {
    RegistryEntry {
        id: "probe-pinned",
        tier: RegistryTier::G1Analytic,
        family: "probe",
        title: "fully pinned probe entry",
        edition: Edition::Exact {
            version: "probe v1",
        },
        source: "synthetic",
        license: LicenseState::Spdx { id: "MIT" },
        deck: DeckPin::AuthoredSpec {
            spec: "PROBE: canonical probe spec bytes.",
        },
        oracle: OracleBinding::SelfContained,
        qois: &[Qoi {
            name: "probe_qoi",
            unit: "1",
            envelope: AcceptanceEnvelope::Tolerance {
                atol: 0.0,
                rtol: 1e-12,
            },
        }],
        notes: "synthetic probe",
    }
}

#[test]
fn g1_seeds_split_into_citable_and_derivation_blocked() {
    let reg = registry();
    let mut citable_g1 = 0_usize;
    let mut derivation_blocked = Vec::new();
    for entry in reg.entries() {
        if entry.tier != RegistryTier::G1Analytic {
            continue;
        }
        match reg.cite(entry.id) {
            Ok(receipt) => {
                citable_g1 += 1;
                assert_eq!(receipt.entry_id(), entry.id);
                assert_eq!(receipt.registry_version(), VVREG_VERSION);
                let deck_digest = entry
                    .deck
                    .digest()
                    .expect("admitted G1 entries carry a pinned deck");
                assert_eq!(receipt.deck_digest(), deck_digest, "{}", entry.id);
                assert_eq!(receipt.entry_digest(), entry_digest(entry), "{}", entry.id);
                assert_eq!(receipt.registry_digest(), reg.digest(), "{}", entry.id);
            }
            Err(CitationRefusal::UnboundOracle { id, .. }) => derivation_blocked.push(id),
            Err(other) => panic!("G1 seed '{}' unexpected refusal: {other}", entry.id),
        }
    }
    assert_eq!(citable_g1, 12, "self-contained G1 analytic seeds");
    assert_eq!(
        derivation_blocked,
        [
            "g1-atkinson-cycle",
            "g1-bennett-linkage-mobility",
            "g1-geneva-closed-form",
            "g1-isentropic-nozzle",
            "g1-riemann-lax",
            "g1-riemann-sod"
        ],
        "decks that delegate load-bearing content are targets, not citable oracles"
    );
}

#[test]
fn unpinned_family_name_citation_fails_with_a_typed_refusal() {
    // THE fixture the bead demands: a family name (TEAM, NAFEMS, CFR,
    // IFToMM, ECN) is not an executable benchmark.
    let reg = registry();
    match reg.cite("g2-team-10") {
        Err(CitationRefusal::UnpinnedEdition { id }) => assert_eq!(id, "g2-team-10"),
        other => panic!("unpinned TEAM 10 must refuse on the edition gate, got {other:?}"),
    }
    // Every seeded G2 row is currently unpinned and must refuse.
    let mut g2_count = 0_usize;
    for entry in reg.entries() {
        if entry.tier != RegistryTier::G2Benchmark {
            continue;
        }
        g2_count += 1;
        assert_eq!(
            entry.oracle,
            OracleBinding::Unpinned,
            "an unpinned G2 target must not claim a self-contained oracle"
        );
        let refusal = reg
            .cite(entry.id)
            .expect_err("no seeded G2 row is citable until its deck is pinned");
        assert!(
            matches!(refusal, CitationRefusal::UnpinnedEdition { id } if id == entry.id),
            "G2 refusals fire on the documented first gate (edition): {entry:?}"
        );
    }
    assert_eq!(g2_count, 15, "the G2 benchmark seed families from the bead");
}

#[test]
fn unknown_and_bare_family_ids_are_refused() {
    let reg = registry();
    for bogus in ["team-10", "nafems", "g9-not-seeded"] {
        match reg.cite(bogus) {
            Err(CitationRefusal::UnknownEntry { id }) => assert_eq!(id, bogus),
            other => panic!("'{bogus}' must be unknown, got {other:?}"),
        }
    }
    for malformed in ["", "TEAM", "g2 team", "-g2-team"] {
        assert_eq!(
            reg.cite(malformed),
            Err(CitationRefusal::InvalidLookupId),
            "malformed lookup ids refuse without copying their input"
        );
    }
}

#[test]
fn oversized_lookup_ids_are_refused_before_the_copy() {
    let oversized = "a".repeat(MAX_LOOKUP_ID_LEN + 1);
    assert_eq!(
        registry().cite(&oversized),
        Err(CitationRefusal::OversizedLookupId {
            len: MAX_LOOKUP_ID_LEN + 1
        })
    );
    let at_cap = format!("a{}", "1".repeat(MAX_LOOKUP_ID_LEN - 1));
    assert!(matches!(
        registry().cite(&at_cap),
        Err(CitationRefusal::UnknownEntry { id }) if id == at_cap
    ));
}

#[test]
fn entry_validation_uses_the_same_bounded_slug_rule_as_lookup() {
    let mut malformed = pinned_probe();
    malformed.id = "TEAM";
    assert_eq!(
        validate_entry(&malformed),
        Err(CitationRefusal::InvalidLookupId)
    );

    let oversized = Box::leak("a".repeat(MAX_LOOKUP_ID_LEN + 1).into_boxed_str());
    malformed.id = oversized;
    assert_eq!(
        validate_entry(&malformed),
        Err(CitationRefusal::OversizedLookupId {
            len: MAX_LOOKUP_ID_LEN + 1
        })
    );
}

#[test]
fn registry_lint_partitions_citable_from_refused_and_finds_no_integrity_defects() {
    let reg = registry();
    let lint = reg.lint();
    assert_eq!(lint.citable.len(), 12);
    assert_eq!(lint.refused.len(), 21);
    assert!(
        lint.integrity.is_empty(),
        "seed integrity defects: {:?}",
        lint.integrity
    );
    assert_eq!(reg.entries().len(), 33);
    assert_eq!(reg.references().len(), 30);
    for id in &lint.citable {
        assert!(id.starts_with("g1-"), "citable id {id} must be a G1 seed");
    }
    let mut unpinned_editions = 0_usize;
    let mut unbound_oracles = 0_usize;
    for refusal in &lint.refused {
        match refusal {
            CitationRefusal::UnpinnedEdition { id } => {
                unpinned_editions += 1;
                assert!(id.starts_with("g2-"), "edition-gate refusals are G2 seeds");
            }
            CitationRefusal::UnboundOracle { id, .. } => {
                unbound_oracles += 1;
                assert!(id.starts_with("g1-"), "oracle-gate refusals are G1 seeds");
            }
            other => panic!("unexpected seed refusal kind: {other:?}"),
        }
    }
    assert_eq!(unpinned_editions, 15);
    assert_eq!(unbound_oracles, 6);
}

#[test]
fn synthetic_rows_and_registries_can_never_mint_receipts() {
    // THE forge regression: a fully pinned synthetic row passes the
    // validation gates...
    let probe = pinned_probe();
    assert_eq!(validate_entry(&probe), Ok(()));
    // ...but validation returns (), not a receipt, and a caller-built
    // registry holding it refuses citation outright: no path from
    // synthetic data to a `CitationReceipt` (and hence to the Verified
    // numerical cap) exists.
    let forged = Registry::build(vec![probe], Vec::new());
    assert_eq!(
        forged.cite("probe-pinned"),
        Err(CitationRefusal::UnauthoritativeRegistry)
    );
    // Receipts from the seeded registry bind its exact content digest.
    let reg = registry();
    let receipt = reg.cite("g1-otto-cycle").expect("citable G1 seed");
    assert_eq!(receipt.registry_digest(), reg.digest());
}

#[test]
fn duplicate_entry_ids_fail_closed_everywhere() {
    let a = pinned_probe();
    let mut b = pinned_probe();
    b.title = "a conflicting second definition of the same id";
    let reg = Registry::build(vec![a, b], Vec::new());
    // lint: nothing citable, every carrying row refused, integrity finding.
    let lint = reg.lint();
    assert!(lint.citable.is_empty(), "a duplicated id is never citable");
    assert_eq!(lint.refused.len(), 2);
    for refusal in &lint.refused {
        assert!(
            matches!(refusal, CitationRefusal::DuplicateEntry { id } if id == "probe-pinned"),
            "duplicated rows are refused as duplicates, got {refusal:?}"
        );
    }
    assert!(
        lint.integrity
            .contains(&IntegrityFinding::DuplicateEntryId { id: "probe-pinned" })
    );
}

/// Helper: identical reference key at a chosen index.
fn dup_key_reference(index: u32) -> fs_vvreg::PrimaryReference {
    fs_vvreg::PrimaryReference {
        index,
        key: "dup-key",
        citation: "c",
        locator: "l",
        anchors: "a",
        boundary: "b",
    }
}

#[test]
fn duplicate_reference_keys_are_found_even_at_different_indices() {
    let reg = Registry::build(Vec::new(), vec![dup_key_reference(1), dup_key_reference(7)]);
    assert!(
        reg.lint()
            .integrity
            .iter()
            .any(|f| matches!(f, IntegrityFinding::DuplicateReference { key: "dup-key" })),
        "same key at non-adjacent indices must still be a collision"
    );
}

#[test]
fn admission_gates_fire_in_the_documented_order() {
    let pinned = pinned_probe();
    assert!(validate_entry(&pinned).is_ok());

    let mut blank_title = pinned;
    blank_title.title = "  ";
    blank_title.edition = Edition::Unpinned; // text gate precedes edition gate
    assert!(matches!(
        validate_entry(&blank_title),
        Err(CitationRefusal::EmptyField { field: "title", .. })
    ));

    let mut blank_notes = pinned;
    blank_notes.notes = ""; // notes are load-bearing assumption boundaries
    assert!(matches!(
        validate_entry(&blank_notes),
        Err(CitationRefusal::EmptyField { field: "notes", .. })
    ));

    let mut unpinned_edition = pinned;
    unpinned_edition.edition = Edition::Unpinned;
    unpinned_edition.license = LicenseState::Unpinned; // edition precedes license
    assert!(matches!(
        validate_entry(&unpinned_edition),
        Err(CitationRefusal::UnpinnedEdition { .. })
    ));

    let mut unpinned_license = pinned;
    unpinned_license.license = LicenseState::Unpinned;
    unpinned_license.deck = DeckPin::Unpinned; // license precedes deck
    assert!(matches!(
        validate_entry(&unpinned_license),
        Err(CitationRefusal::UnpinnedLicense { .. })
    ));

    let mut unpinned_deck = pinned;
    unpinned_deck.deck = DeckPin::Unpinned;
    assert!(matches!(
        validate_entry(&unpinned_deck),
        Err(CitationRefusal::UnpinnedDeck { .. })
    ));

    let mut malformed_deck = pinned;
    malformed_deck.deck = DeckPin::External {
        digest_hex: "not-hex",
    };
    assert!(matches!(
        validate_entry(&malformed_deck),
        Err(CitationRefusal::MalformedDeckDigest { .. })
    ));

    let mut well_formed_external = pinned;
    well_formed_external.deck = DeckPin::External {
        digest_hex: "0000000000000000000000000000000000000000000000000000000000000000",
    };
    assert!(validate_entry(&well_formed_external).is_ok());

    let mut unbound_oracle = pinned;
    unbound_oracle.oracle = OracleBinding::DerivationRequired {
        obligation: "derive it",
    };
    unbound_oracle.qois = &[]; // oracle gate precedes QoI presence gate
    assert!(matches!(
        validate_entry(&unbound_oracle),
        Err(CitationRefusal::UnboundOracle {
            obligation: "derive it",
            ..
        })
    ));

    let mut unpinned_oracle = pinned;
    unpinned_oracle.oracle = OracleBinding::Unpinned;
    assert!(matches!(
        validate_entry(&unpinned_oracle),
        Err(CitationRefusal::UnpinnedOracle { .. })
    ));

    let mut no_qois = pinned;
    no_qois.qois = &[];
    assert!(matches!(
        validate_entry(&no_qois),
        Err(CitationRefusal::MissingQois { .. })
    ));

    let mut duplicate_qoi = pinned;
    duplicate_qoi.qois = &[
        Qoi {
            name: "q",
            unit: "1",
            envelope: AcceptanceEnvelope::Tolerance {
                atol: 0.0,
                rtol: 1e-9,
            },
        },
        Qoi {
            name: "q",
            unit: "m",
            envelope: AcceptanceEnvelope::Unpinned, // dedup precedes envelope gate
        },
    ];
    assert!(matches!(
        validate_entry(&duplicate_qoi),
        Err(CitationRefusal::DuplicateQoi { qoi: "q", .. })
    ));

    let mut unpinned_envelope = pinned;
    unpinned_envelope.qois = &[Qoi {
        name: "q",
        unit: "1",
        envelope: AcceptanceEnvelope::Unpinned,
    }];
    assert!(matches!(
        validate_entry(&unpinned_envelope),
        Err(CitationRefusal::UnpinnedEnvelope { qoi: "q", .. })
    ));

    // The QoI-count cap fires before every QoI traversal, bounding the row
    // count work even when each hostile QoI would fail the text gate.
    let hostile_qoi = Qoi {
        name: "",
        unit: "",
        envelope: AcceptanceEnvelope::Unpinned,
    };
    let mut too_many = pinned;
    too_many.qois = Box::leak(vec![hostile_qoi; MAX_QOIS_PER_ENTRY + 1].into_boxed_slice());
    assert!(matches!(
        validate_entry(&too_many),
        Err(CitationRefusal::TooManyQois { count, .. }) if count == MAX_QOIS_PER_ENTRY + 1
    ));
}

#[test]
fn invalid_envelopes_are_refused_with_reasons() {
    let cases: [(AcceptanceEnvelope, &str); 5] = [
        (
            AcceptanceEnvelope::Tolerance {
                atol: f64::INFINITY,
                rtol: 0.0,
            },
            "non-finite tolerance",
        ),
        (
            AcceptanceEnvelope::Tolerance {
                atol: -1.0,
                rtol: 0.0,
            },
            "negative tolerance",
        ),
        (
            AcceptanceEnvelope::Tolerance {
                atol: 0.0,
                rtol: 0.0,
            },
            "zero-width tolerance (declare an Interval for exact claims)",
        ),
        (
            AcceptanceEnvelope::Interval {
                lo: f64::NAN,
                hi: 1.0,
            },
            "non-finite bound",
        ),
        (
            AcceptanceEnvelope::Interval { lo: 2.0, hi: 1.0 },
            "inverted interval",
        ),
    ];
    for (envelope, expected_reason) in cases {
        let mut entry = pinned_probe();
        let qois = [Qoi {
            name: "q",
            unit: "1",
            envelope,
        }];
        // Leak is fine in a test: entries want 'static slices.
        entry.qois = Box::leak(Box::new(qois));
        match validate_entry(&entry) {
            Err(CitationRefusal::InvalidEnvelope { reason, .. }) => {
                assert_eq!(reason, expected_reason);
            }
            other => panic!("envelope {envelope:?} must refuse, got {other:?}"),
        }
    }
}

#[test]
fn executable_tolerance_gate_pins_boundary_and_seeded_violation_verdicts() {
    let reg = registry();
    let entry = reg
        .entry("g1-block-incline-stick-slip")
        .expect("seeded entry");
    let boundary = 1e-12_f64;
    let admitted = reg
        .check_acceptance_envelope(
            "g1-block-incline-stick-slip",
            "critical_angle",
            EnvelopeObservation::AgainstReference {
                reference: 0.0,
                computed: boundary,
            },
        )
        .expect("the inclusive tolerance boundary must pass");
    assert_eq!(admitted.attempt().entry_id(), entry.id);
    assert_eq!(admitted.attempt().qoi(), "critical_angle");
    assert_eq!(admitted.attempt().unit(), "rad");
    assert_eq!(admitted.attempt().entry_digest(), entry_digest(entry));
    assert_eq!(admitted.attempt().registry_digest(), reg.digest());
    assert_eq!(admitted.attempt().registry_version(), VVREG_VERSION);
    assert_eq!(admitted.margin(), 0.0);
    assert!(admitted.passed());
    assert_eq!(
        admitted.evaluation(),
        EnvelopeEvaluation::Tolerance {
            reference: 0.0,
            atol: boundary,
            rtol: 0.0,
            allowed: boundary,
            deviation: boundary,
        }
    );
    let boundary_bits = boundary.to_bits();
    let expected = format!(
        concat!(
            "{{\"vvreg\":\"acceptance-envelope-verdict\",\"schema\":{},",
            "\"registry_version\":{},\"registry_digest\":\"{}\",",
            "\"entry\":\"g1-block-incline-stick-slip\",\"entry_digest\":\"{}\",",
            "\"qoi\":\"critical_angle\",\"unit\":\"rad\",\"mode\":\"tolerance\",",
            "\"reference\":\"0x0000000000000000\",\"computed\":\"0x{:016x}\",",
            "\"atol\":\"0x{:016x}\",\"rtol\":\"0x0000000000000000\",",
            "\"allowed\":\"0x{:016x}\",\"deviation\":\"0x{:016x}\",",
            "\"margin\":\"0x0000000000000000\",\"pass\":true}}"
        ),
        ENVELOPE_VERDICT_SCHEMA_VERSION,
        VVREG_VERSION,
        reg.digest().to_hex(),
        entry_digest(entry).to_hex(),
        boundary_bits,
        boundary_bits,
        boundary_bits,
        boundary_bits,
    );
    assert_eq!(
        admitted.json_line(),
        expected,
        "fixed-field exact-bit golden"
    );

    let lower = reg
        .check_acceptance_envelope(
            "g1-block-incline-stick-slip",
            "critical_angle",
            EnvelopeObservation::AgainstReference {
                reference: 0.0,
                computed: -boundary,
            },
        )
        .expect("the negative inclusive tolerance boundary must pass");
    assert_eq!(lower.margin(), 0.0);

    // A disclosed deterministic seed chooses a positive 1..=8 ULP
    // perturbation beyond the exact boundary. This is the G2 meta-test: a
    // seeded QoI corruption must turn the executable envelope red.
    let seed = 0x6E_B3_2026_0717_u64;
    let ulps = (seed.rotate_left(17) ^ seed.rotate_right(9)) % 8 + 1;
    assert_eq!(ulps, 4, "the disclosed seed pins the corruption distance");
    let corrupted = f64::from_bits(boundary.to_bits() + ulps);
    let refusal = reg
        .check_acceptance_envelope(
            "g1-block-incline-stick-slip",
            "critical_angle",
            EnvelopeObservation::AgainstReference {
                reference: 0.0,
                computed: corrupted,
            },
        )
        .expect_err("a seeded boundary+ULP perturbation must fail the gate");
    let EnvelopeGateError::Violation { verdict } = refusal else {
        panic!("finite outside-envelope candidate must retain a failing verdict");
    };
    assert!(!verdict.passed());
    assert!(verdict.margin() < 0.0);
    let line = verdict.json_line();
    assert!(line.contains("\"entry\":\"g1-block-incline-stick-slip\""));
    assert!(line.contains("\"reference\":\"0x0000000000000000\""));
    assert!(line.contains(&format!("\"computed\":\"0x{:016x}\"", corrupted.to_bits())));
    assert!(line.contains("\"allowed\":"));
    assert!(line.contains("\"deviation\":"));
    assert!(line.contains("\"margin\":"));
    assert!(line.ends_with("\"pass\":false}"));
}

#[test]
fn executable_interval_gate_is_inclusive_and_mode_explicit() {
    let reg = registry();
    let boundary = 1.0_f64;
    let verdict = reg
        .check_acceptance_envelope(
            "g1-bennett-linkage-mobility",
            "mobility_dof",
            EnvelopeObservation::AgainstInterval { computed: boundary },
        )
        .expect("the singleton interval endpoint is inclusive");
    assert_eq!(verdict.margin(), 0.0);
    assert!(verdict.passed());
    assert!(verdict.json_line().contains("\"reference\":null"));
    assert_eq!(verdict.attempt().qoi(), "mobility_dof");

    let outside = f64::from_bits(1.0_f64.to_bits() + 1);
    let refusal = reg
        .check_acceptance_envelope(
            "g1-bennett-linkage-mobility",
            "mobility_dof",
            EnvelopeObservation::AgainstInterval { computed: outside },
        )
        .expect_err("one ULP beyond the singleton interval must refuse");
    assert!(matches!(
        refusal,
        EnvelopeGateError::Violation { verdict }
            if !verdict.passed() && verdict.margin() < 0.0
    ));

    let mismatch = reg
        .check_acceptance_envelope(
            "g1-bennett-linkage-mobility",
            "mobility_dof",
            EnvelopeObservation::AgainstReference {
                reference: 0.0,
                computed: 0.0,
            },
        )
        .expect_err("an interval cannot consume a reference observation");
    assert!(matches!(
        &mismatch,
        EnvelopeGateError::ModeMismatch {
            attempt,
            expected: "interval",
            got: "tolerance",
        } if attempt.entry_id() == "g1-bennett-linkage-mobility"
    ));
    let mismatch_line = mismatch.json_line();
    assert!(mismatch_line.contains("\"envelope_mode\":\"interval\""));
    assert!(mismatch_line.contains("\"observation_mode\":\"reference\""));
    assert!(mismatch_line.contains("\"outcome\":\"mode-mismatch\""));
}

#[test]
fn executable_envelope_gate_retains_nonfinite_and_overflow_attempts() {
    let reg = registry();
    for (reference, computed, field) in [
        (f64::NAN, 1.0, "reference"),
        (1.0, f64::INFINITY, "computed"),
    ] {
        let refusal = reg
            .check_acceptance_envelope(
                "g1-block-incline-stick-slip",
                "critical_angle",
                EnvelopeObservation::AgainstReference {
                    reference,
                    computed,
                },
            )
            .expect_err("non-finite observations must refuse");
        let EnvelopeGateError::NonFiniteInput {
            attempt,
            field: got,
        } = &refusal
        else {
            panic!("non-finite input must retain its registry-bound attempt");
        };
        assert_eq!(*got, field);
        let EnvelopeObservation::AgainstReference {
            reference: retained_reference,
            computed: retained_computed,
        } = attempt.observation()
        else {
            panic!("tolerance attempt must retain reference observation mode");
        };
        assert_eq!(retained_reference.to_bits(), reference.to_bits());
        assert_eq!(retained_computed.to_bits(), computed.to_bits());
        let line = refusal.json_line();
        assert!(line.contains("\"outcome\":\"non-finite-input\""));
        assert!(line.contains(&format!("\"reference\":\"0x{:016x}\"", reference.to_bits())));
        assert!(line.contains(&format!("\"computed\":\"0x{:016x}\"", computed.to_bits())));
    }

    let overflow = reg
        .check_acceptance_envelope(
            "g1-block-incline-stick-slip",
            "critical_angle",
            EnvelopeObservation::AgainstReference {
                reference: -f64::MAX,
                computed: f64::MAX,
            },
        )
        .expect_err("finite subtraction overflow must refuse");
    assert!(matches!(
        &overflow,
        EnvelopeGateError::ArithmeticOverflow {
            attempt,
            operation: "computed - reference",
        } if attempt.entry_digest()
            == entry_digest(reg.entry("g1-block-incline-stick-slip").expect("seeded"))
    ));
    let overflow_line = overflow.json_line();
    assert!(overflow_line.contains("\"reference\":\"0xffefffffffffffff\""));
    assert!(overflow_line.contains("\"computed\":\"0x7fefffffffffffff\""));
    assert!(overflow_line.contains("\"outcome\":\"arithmetic-overflow\""));
}

#[test]
fn executable_envelope_gate_refuses_unpinned_and_untrusted_bindings() {
    let reg = registry();
    assert!(matches!(
        reg.check_acceptance_envelope(
            "g2-team-10",
            "average_flux_density_probe",
            EnvelopeObservation::AgainstInterval { computed: 0.0 },
        ),
        Err(EnvelopeGateError::Registry(
            CitationRefusal::UnpinnedEnvelope {
                id: "g2-team-10",
                qoi: "average_flux_density_probe",
            }
        ))
    ));

    assert!(matches!(
        reg.check_acceptance_envelope(
            "g1-block-incline-stick-slip",
            "forged_looser_qoi",
            EnvelopeObservation::AgainstReference {
                reference: 0.0,
                computed: 0.0,
            },
        ),
        Err(EnvelopeGateError::UnknownQoi { id, qoi })
            if id == "g1-block-incline-stick-slip" && qoi == "forged_looser_qoi"
    ));
    assert!(matches!(
        reg.check_acceptance_envelope(
            "g1-does-not-exist",
            "critical_angle",
            EnvelopeObservation::AgainstReference {
                reference: 0.0,
                computed: 0.0,
            },
        ),
        Err(EnvelopeGateError::Registry(CitationRefusal::UnknownEntry { id }))
            if id == "g1-does-not-exist"
    ));

    let untrusted = Registry::build(vec![pinned_probe()], Vec::new());
    assert!(matches!(
        untrusted.check_acceptance_envelope(
            "probe-pinned",
            "probe_qoi",
            EnvelopeObservation::AgainstReference {
                reference: 1.0,
                computed: 1.0,
            },
        ),
        Err(EnvelopeGateError::Registry(
            CitationRefusal::UnauthoritativeRegistry
        ))
    ));
}

#[test]
fn canonical_serialization_golden_for_the_unpinned_team_10_row() {
    let reg = registry();
    let entry = reg.entry("g2-team-10").expect("seeded");
    let row = canonical_row(entry);
    let expected = concat!(
        "{\"id\":\"g2-team-10\",\"tier\":\"G2\",\"family\":\"TEAM\",",
        "\"title\":\"TEAM problem 10: steel plates around a coil, nonlinear transient eddy current\",",
        "\"edition\":null,",
        "\"source\":\"COMPUMAG TEAM benchmark suite, official problem definition 10\",",
        "\"license\":null,\"deck\":null,\"oracle\":null,",
        "\"qois\":[{\"name\":\"average_flux_density_probe\",\"unit\":\"T\",\"envelope\":null}],",
        "\"notes\":\"Family name only: exact geometry revision, excitation, material law, ",
        "circuit, QoI set, and acceptance data must be pinned before citation.\"}",
    );
    assert_eq!(row, expected);
}

#[test]
fn canonical_rows_are_deterministic_sorted_and_bit_exact_for_floats() {
    let reg = registry();
    let rows_a = reg.canonical_rows();
    let rows_b = reg.canonical_rows();
    assert_eq!(rows_a, rows_b, "serialization is deterministic");
    let mut sorted = rows_a.clone();
    sorted.sort();
    assert_eq!(rows_a, sorted, "rows are emitted in sorted-id order");
    let hertz = canonical_row(reg.entry("g1-hertz-sphere-plane").expect("seeded"));
    assert!(hertz.contains("\"tier\":\"G1\""));
    assert!(hertz.contains("{\"spdx\":\"MIT OR Apache-2.0\"}"));
    assert!(hertz.contains("\"oracle\":\"self-contained\""));
    assert!(hertz.contains("\"kind\":\"tolerance\""));
    // Floats are IEEE-754 bit tokens, never decimal formatting.
    assert!(hertz.contains("\"atol\":\"0x"));
    // The authored deck digest is present, variant-tagged, as 64 hex chars.
    let marker = "\"deck\":{\"authored\":\"";
    let deck_pos = hertz.find(marker).expect("authored deck digest");
    let digest_start = deck_pos + marker.len();
    let digest = &hertz[digest_start..digest_start + 64];
    assert!(digest.chars().all(|c| c.is_ascii_hexdigit()));
}

#[test]
fn deck_rows_and_identities_agree_on_variant_state_and_spelling() {
    // Valid external hex is normalized: case cannot fork identity or row.
    let mut upper = pinned_probe();
    upper.deck = DeckPin::External {
        digest_hex: "ABCDEFABCDEFABCDEFABCDEFABCDEFABCDEFABCDEFABCDEFABCDEFABCDEF1234",
    };
    let mut lower = pinned_probe();
    lower.deck = DeckPin::External {
        digest_hex: "abcdefabcdefabcdefabcdefabcdefabcdefabcdefabcdefabcdefabcdef1234",
    };
    assert_eq!(entry_digest(&upper), entry_digest(&lower));
    assert_eq!(canonical_row(&upper), canonical_row(&lower));
    assert!(canonical_row(&upper).contains("{\"external\":\"abcdef"));

    // Malformed external is a distinct state from unpinned, in both the
    // identity and the row.
    let mut malformed = pinned_probe();
    malformed.deck = DeckPin::External { digest_hex: "zz" };
    let mut unpinned = pinned_probe();
    unpinned.deck = DeckPin::Unpinned;
    assert_ne!(entry_digest(&malformed), entry_digest(&unpinned));
    assert!(canonical_row(&malformed).contains("{\"malformed\":\"zz\"}"));
    assert!(canonical_row(&unpinned).contains("\"deck\":null"));

    // An authored deck and an external deck carrying the SAME hex are
    // different identities and different row spellings (variant preserved).
    let authored = pinned_probe();
    let authored_hex = authored.deck.digest().expect("authored digest").to_hex();
    let mut external_twin = pinned_probe();
    external_twin.deck = DeckPin::External {
        digest_hex: Box::leak(authored_hex.into_boxed_str()),
    };
    assert_ne!(entry_digest(&authored), entry_digest(&external_twin));
    assert!(canonical_row(&authored).contains("{\"authored\":\""));
    assert!(canonical_row(&external_twin).contains("{\"external\":\""));
}

#[test]
fn registry_digest_is_deterministic_and_order_canonicalized() {
    let reg = registry();
    // Two INDEPENDENTLY built equivalent registries agree (not a self-
    // comparison), and input order cannot move the identity.
    let same_order = Registry::build(reg.entries().to_vec(), reg.references().to_vec());
    let mut reversed_entries: Vec<RegistryEntry> = reg.entries().to_vec();
    reversed_entries.reverse();
    let reversed = Registry::build(reversed_entries, reg.references().to_vec());
    let original_digest = reg.digest();
    assert_eq!(same_order.digest(), original_digest);
    assert_eq!(reversed.digest(), original_digest);
    assert_ne!(original_digest.as_bytes(), &[0_u8; 32]);

    // Even CONFLICTING duplicate-id rows land in one canonical order:
    // reversing the input cannot move the digest (content tie-break).
    let a = pinned_probe();
    let mut b = pinned_probe();
    b.title = "a conflicting second definition of the same id";
    let ab = Registry::build(vec![a, b], Vec::new());
    let ba = Registry::build(vec![b, a], Vec::new());
    assert_eq!(
        ab.digest(),
        ba.digest(),
        "duplicate-id tie-break is canonical"
    );
    // Same for duplicate reference keys with otherwise-different fields.
    let mut r1 = dup_key_reference(3);
    let mut r2 = dup_key_reference(3);
    r1.citation = "citation one";
    r2.citation = "citation two";
    let r12 = Registry::build(Vec::new(), vec![r1, r2]);
    let r21 = Registry::build(Vec::new(), vec![r2, r1]);
    assert_eq!(
        r12.digest(),
        r21.digest(),
        "reference tie-break is canonical"
    );
}

#[test]
fn entry_identity_is_mutation_sensitive_in_every_semantic_field() {
    let base = pinned_probe();
    let base_digest = entry_digest(&base);
    assert_eq!(base_digest, entry_digest(&base), "stable identity");

    let mut id = base;
    id.id = "probe-pinned-2";
    assert_ne!(entry_digest(&id), base_digest);

    let mut family = base;
    family.family = "different family";
    assert_ne!(entry_digest(&family), base_digest);

    let mut title = base;
    title.title = "different title";
    assert_ne!(entry_digest(&title), base_digest);

    let mut source = base;
    source.source = "different source";
    assert_ne!(entry_digest(&source), base_digest);

    let mut notes = base;
    notes.notes = "different notes";
    assert_ne!(entry_digest(&notes), base_digest);

    let mut tier = base;
    tier.tier = RegistryTier::G2Benchmark;
    assert_ne!(entry_digest(&tier), base_digest);

    let mut edition = base;
    edition.edition = Edition::Exact {
        version: "probe v2",
    };
    assert_ne!(entry_digest(&edition), base_digest);

    let mut license = base;
    license.license = LicenseState::Spdx { id: "Apache-2.0" };
    assert_ne!(entry_digest(&license), base_digest);

    let mut deck = base;
    deck.deck = DeckPin::AuthoredSpec {
        spec: "PROBE: different spec bytes.",
    };
    assert_ne!(entry_digest(&deck), base_digest);

    let mut oracle = base;
    oracle.oracle = OracleBinding::DerivationRequired {
        obligation: "derive it",
    };
    assert_ne!(entry_digest(&oracle), base_digest);

    let mut unpinned_oracle = base;
    unpinned_oracle.oracle = OracleBinding::Unpinned;
    assert_ne!(entry_digest(&unpinned_oracle), base_digest);
    assert_ne!(entry_digest(&unpinned_oracle), entry_digest(&oracle));

    let mut envelope = base;
    envelope.qois = &[Qoi {
        name: "probe_qoi",
        unit: "1",
        envelope: AcceptanceEnvelope::Tolerance {
            atol: 0.0,
            rtol: 1e-11,
        },
    }];
    assert_ne!(entry_digest(&envelope), base_digest);

    let mut qoi_name = base;
    qoi_name.qois = &[Qoi {
        name: "renamed_qoi",
        unit: "1",
        envelope: AcceptanceEnvelope::Tolerance {
            atol: 0.0,
            rtol: 1e-12,
        },
    }];
    assert_ne!(entry_digest(&qoi_name), base_digest);

    let mut qoi_unit = base;
    qoi_unit.qois = &[Qoi {
        name: "probe_qoi",
        unit: "m",
        envelope: AcceptanceEnvelope::Tolerance {
            atol: 0.0,
            rtol: 1e-12,
        },
    }];
    assert_ne!(entry_digest(&qoi_unit), base_digest);

    // QoI count and QoI order both move the identity.
    let mut one = base;
    one.qois = &[QOI_A];
    let mut two_ab = base;
    two_ab.qois = &[QOI_A, QOI_B];
    let mut two_ba = base;
    two_ba.qois = &[QOI_B, QOI_A];
    assert_ne!(entry_digest(&one), entry_digest(&two_ab), "count moves it");
    assert_ne!(
        entry_digest(&two_ab),
        entry_digest(&two_ba),
        "declared QoI order is semantic and moves it"
    );
}

#[test]
fn blank_authored_specs_have_no_wellformed_digest_and_refuse_admission() {
    let mut blank = pinned_probe();
    blank.deck = DeckPin::AuthoredSpec { spec: "   " };
    assert_eq!(blank.deck.digest(), None, "blank spec is not well-formed");
    assert!(canonical_row(&blank).contains("{\"authored\":null}"));
    assert!(matches!(
        validate_entry(&blank),
        Err(CitationRefusal::EmptyField {
            field: "deck.spec",
            ..
        })
    ));
}

#[test]
fn oversized_bead_ids_are_refused_before_the_copy() {
    let reg = registry();
    let receipt = reg.cite("g1-otto-cycle").expect("citable G1 seed");
    let oversized = " ".repeat(MAX_BEAD_ID_LEN + 1);
    assert_eq!(
        ConsumptionRecord::bind(&receipt, &oversized, ConsumptionStatus::Read),
        Err(ConsumptionRefusal::OversizedBead {
            len: MAX_BEAD_ID_LEN + 1
        })
    );
    let at_cap = "b".repeat(MAX_BEAD_ID_LEN);
    assert!(ConsumptionRecord::bind(&receipt, &at_cap, ConsumptionStatus::Read).is_ok());
}

#[test]
fn consumption_records_pin_the_exact_artifact_version() {
    let reg = registry();
    let receipt = reg.cite("g1-stefan-problem").expect("citable G1 seed");
    for status in ConsumptionStatus::all() {
        let record = ConsumptionRecord::bind(&receipt, "frankensim-example-bead", status)
            .expect("non-blank bead binds");
        assert_eq!(record.entry_id, "g1-stefan-problem");
        assert_eq!(record.entry_digest, receipt.entry_digest());
        assert_eq!(record.registry_version, VVREG_VERSION);
        assert_eq!(record.status, status);
    }
    let tags: Vec<&str> = ConsumptionStatus::all().iter().map(|s| s.tag()).collect();
    assert_eq!(
        tags,
        [
            "unread",
            "read",
            "derived",
            "reproduced",
            "independently_falsified"
        ]
    );
    assert_eq!(
        ConsumptionRecord::bind(&receipt, "   ", ConsumptionStatus::Read),
        Err(ConsumptionRefusal::EmptyBead)
    );
}

#[test]
fn color_rule_caps_never_launder_publisher_authority() {
    let reg = registry();
    let receipt = reg.cite("g1-otto-cycle").expect("citable G1 seed");
    assert_eq!(receipt.numerical_claim_cap(), ColorRank::Verified);
    // The physical cap is unconditionally Estimated in this slice: a bare
    // caller-asserted "held-out evidence" flag would be forgeable, so no
    // upgrade path is offered until a typed evidence binding exists.
    assert_eq!(receipt.physical_claim_cap(), ColorRank::Estimated);
    // The lattice ordering the caps rely on.
    assert!(ColorRank::Estimated < ColorRank::Validated);
    assert!(ColorRank::Validated < ColorRank::Verified);
}

/// The pinned (key, locator) identity of the full 30-reference seed:
/// omission, reordering, or citation-target drift fails here; the
/// remaining prose fields are covered by the per-field digest mutation
/// locks below plus the registry digest.
const EXPECTED_REFERENCES: [(&str, &str); 30] = [
    ("feec-stability-afw", "arXiv:0906.4325"),
    ("sheaf-spectra-hansen-ghrist", "arXiv:1808.01513"),
    ("sheaf-cosheaf-curry", "arXiv:1303.3255"),
    ("nasa9-thermo-mcbride", "NASA/TP-2002-211556"),
    ("entropy-stable-tadmor", "doi:10.1016/bs.hna.2016.09.006"),
    (
        "port-hamiltonian-dirac-cervera",
        "doi:10.1016/j.automatica.2006.08.014",
    ),
    ("rigidity-maxwell-calladine-rocks", "arXiv:2208.07419"),
    (
        "contact-ipc-li",
        "doi:10.1145/3386569.3392425; arXiv:2307.15908",
    ),
    ("codim-ipc-li", "arXiv:2012.04457"),
    ("nonsmooth-contact-acary", "arXiv:1410.2499"),
    ("validated-flowpipes-walawska-wilczak", "arXiv:1509.07388"),
    ("gear-te-athavale", "SAE 2001-01-1507"),
    ("wankel-seals-handschuh-owen", "NASA/TM-2010-216353"),
    (
        "iftomm-benchmark-library",
        "IFToMM benchmark library (exact artifact per entry)",
    ),
    (
        "team-em-benchmarks",
        "COMPUMAG TEAM suite (exact revision per entry)",
    ),
    (
        "nafems-thermal",
        "NAFEMS benchmark index (exact case ID/report/license per entry)",
    ),
    ("nasa9-cantera-oracle", "Cantera 3.2 release documentation"),
    (
        "nonholonomic-modin-verdier",
        "doi:10.1007/s00211-020-01126-y",
    ),
    ("hcurl-ams-hiptmair-xu", "doi:10.1137/060660588"),
    (
        "switched-descriptor-yildiz",
        "doi:10.1142/S0218126613500461",
    ),
    (
        "nonlinear-magnetic-energy-mandlmayr-egger",
        "arXiv:2311.02380",
    ),
    (
        "reacting-entropy-ching",
        "arXiv:2211.16254; arXiv:2211.16297",
    ),
    (
        "acoustics-burton-miller-fwh",
        "doi:10.1016/0022-247X(84)90146-X; doi:10.1098/rsta.1969.0031",
    ),
    ("tribology-ehl-hamrock-dowson", "NASA-TP-1342"),
    ("fatigue-astm-nasgro", "ASTM E466-21; ASTM E647-24"),
    (
        "assurance-iso12100-iec60034",
        "ISO 12100:2010; IEC 60034-1:2026; IEC 60034-2-1:2024",
    ),
    (
        "vvuq-asme-gum-nasa7009",
        "ASME VVUQ; JCGM GUM; NASA-STD-7009",
    ),
    (
        "combustion-sandia-ecn",
        "ECN data archive (exact configuration per entry)",
    ),
    ("interop-fmi-ssp", "FMI 3.0.2; SSP 2.0"),
    ("gear-iso6336", "ISO 6336-1:2019"),
];

#[test]
fn primary_reference_seed_is_complete_and_indexed() {
    let reg = registry();
    let actual: Vec<(u32, &str, &str)> = reg
        .references()
        .iter()
        .map(|r| (r.index, r.key, r.locator))
        .collect();
    let expected: Vec<(u32, &str, &str)> = EXPECTED_REFERENCES
        .iter()
        .enumerate()
        .map(|(i, (key, locator))| (u32::try_from(i).expect("30 fits") + 1, *key, *locator))
        .collect();
    assert_eq!(
        actual, expected,
        "all 30 bead references, keyed and located"
    );
    let feec = reg
        .reference("feec-stability-afw")
        .expect("reference lookup by key");
    assert_eq!(feec.index, 1);
    // The registry lint already proves no reference field is blank.
    assert!(registry().lint().integrity.is_empty());
}

#[test]
fn every_reference_field_moves_the_registry_digest() {
    let base_ref = dup_key_reference(1);
    let base = Registry::build(Vec::new(), vec![base_ref]).digest();
    let mut index = base_ref;
    index.index = 2;
    let mut key = base_ref;
    key.key = "other-key";
    let mut citation = base_ref;
    citation.citation = "c2";
    let mut locator = base_ref;
    locator.locator = "l2";
    let mut anchors = base_ref;
    anchors.anchors = "a2";
    let mut boundary = base_ref;
    boundary.boundary = "b2";
    for mutated in [index, key, citation, locator, anchors, boundary] {
        let digest = Registry::build(Vec::new(), vec![mutated]).digest();
        assert_ne!(digest, base, "field mutation must move digest: {mutated:?}");
    }
}
