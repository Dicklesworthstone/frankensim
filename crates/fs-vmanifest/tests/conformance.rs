//! Cross-instance conformance battery for the VerificationManifest schema and
//! the authored I01/I02/I03/I04/I08/I12/I15 drafts.
//!
//! Coverage maps to the bead's test plan: G0 schema/property gates
//! (required fields, identity stability, lattice separation, amendment
//! invalidation), G3 mutations (weakened hypotheses, swapped held-out
//! fixtures, loosened bands, production-oracle reuse — each fails closed
//! or moves the identity), chunked in-memory assembly identity as a G4
//! precursor, and G5 byte-stable canonicalization on the same ISA.

use fs_vmanifest::{
    Ambition, AmendmentRefusal, CampaignTier, ClaimPolarity, ClaimSpec, FiveExplicits, FixturePin,
    FixtureSource, FreezeRefusal, GauntletTier, MAX_CLAIMS, MAX_FIXTURES, MAX_MANIFEST_TEXT_BYTES,
    MAX_OBLIGATIONS, MAX_ROW_ITEMS, MAX_WAIVERS, ManifestDraft, ObligationRow, OracleRoute,
    Partition, ToleranceSemantics, VMANIFEST_SCHEMA_VERSION, Waiver, claim_digest, fixture_digest,
    i01_draft, obligation_digest, waiver_digest,
};

const PROBE_CLAIM: ClaimSpec = ClaimSpec {
    id: "probe-claim",
    ambition: Ambition::Solid,
    polarity: ClaimPolarity::Affirmative,
    statement: "probe statement",
    hypotheses: &["probe hypothesis"],
    qoi: "probe_qoi",
    unit: "1",
    tolerance: ToleranceSemantics::Exact,
    evidence_tier: GauntletTier::G0,
    oracle: OracleRoute {
        identity: "independent probe oracle",
        independent: true,
        tcb_overlap: "none",
    },
    activation: "immediately",
    kill: "never observed",
    fallback: "none needed",
    no_claim: "probe only",
};

const PROBE_FIXTURE: FixturePin = FixturePin {
    id: "probe-deck",
    source: FixtureSource::AuthoredSpec {
        spec: "PROBE: canonical probe deck bytes.",
    },
    partition: Partition::Development,
};

const PROBE_OBLIGATION: ObligationRow = ObligationRow {
    leaf: "probe-leaf",
    claims_covered: &["probe-claim"],
    unit_cases: &["happy", "empty", "error"],
    g0: "probe generators and laws",
    decks: &["probe-deck"],
    g3_relations: &["probe relabeling invariance"],
    g4_schedule: "cancel once at the only boundary",
    g5_matrix: "threads {1,2} x deterministic",
    entry_point: "scripts/e2e/leapfrog/probe.sh",
    tier: CampaignTier::Core,
    dsr_lane: "dsr quality --tool frankensim (probe)",
    obs_events: &["probe.event"],
    replay_command: "scripts/e2e/leapfrog/probe.sh --replay <artifact-id>",
};

fn probe_draft() -> ManifestDraft {
    ManifestDraft {
        initiative: "PROBE",
        title: "probe manifest",
        version: 1,
        explicits: FiveExplicits {
            units: "SI",
            seeds: "Philox probe streams",
            budgets: "1 s",
            versions: "schema v2",
            capabilities: "none",
        },
        claims: vec![PROBE_CLAIM],
        fixtures: vec![PROBE_FIXTURE],
        obligations: vec![PROBE_OBLIGATION],
        waivers: Vec::new(),
        amendment_rules: "successor versions only",
    }
}

fn shared_probe_draft() -> ManifestDraft {
    let mut sibling = PROBE_CLAIM;
    sibling.id = "probe-sibling";
    sibling.statement = "independent sibling statement";
    sibling.qoi = "probe_sibling_qoi";

    let mut shared = PROBE_OBLIGATION;
    shared.claims_covered = &["probe-claim", "probe-sibling"];

    let mut draft = probe_draft();
    draft.claims = vec![PROBE_CLAIM, sibling];
    draft.obligations = vec![shared];
    draft
}

const fn probe_waiver(subject: &'static str) -> Waiver {
    Waiver {
        subject,
        reason: "probe waiver reason",
        owner: "probe waiver owner",
        predicate: "probe waiver retirement predicate",
        expiry: "next probe review",
        promotion_effect: "affected probe evidence cannot promote",
    }
}

#[test]
fn i01_draft_freezes_with_the_declared_lattice_and_coverage() {
    let frozen = i01_draft().freeze().expect("the I01 seed must freeze");
    assert_eq!(frozen.initiative(), "I01");
    assert_eq!(frozen.version(), 1);
    assert_eq!(frozen.claims().len(), 9);
    let solid = frozen
        .claims()
        .iter()
        .filter(|c| c.ambition == Ambition::Solid)
        .count();
    let frontier = frozen
        .claims()
        .iter()
        .filter(|c| c.ambition == Ambition::Frontier)
        .count();
    let moonshot = frozen
        .claims()
        .iter()
        .filter(|c| c.ambition == Ambition::Moonshot)
        .count();
    assert_eq!(
        (solid, frontier, moonshot),
        (5, 2, 2),
        "baseline and maximal lattice elements stay separate"
    );
    assert_eq!(frozen.fixtures().len(), 6);
    let held_out = frozen
        .fixtures()
        .iter()
        .filter(|x| x.partition == Partition::HeldOut)
        .count();
    assert_eq!(held_out, 2, "held-out partitions are frozen up front");
    assert_eq!(frozen.obligations().len(), 6);
    assert_eq!(frozen.waivers().len(), 1);
    // The falsifier lane is preregistered with refutation polarity.
    let completeness = frozen
        .claim("i01-generated-law-completeness")
        .expect("moonshot completeness claim");
    assert_eq!(completeness.polarity, ClaimPolarity::Refutation);
    // Identity is stable and non-trivial.
    let again = i01_draft().freeze().expect("refreeze");
    assert_eq!(frozen.digest(), again.digest());
    assert_ne!(frozen.digest().as_bytes(), &[0_u8; 32]);
}

#[test]
fn authored_five_explicits_avoid_known_legacy_instance_revision_spellings() {
    assert_eq!(VMANIFEST_SCHEMA_VERSION, 2);
    for draft in [
        fs_vmanifest::i01_draft(),
        fs_vmanifest::i02_draft(),
        fs_vmanifest::i03_draft(),
        fs_vmanifest::i04_draft(),
        fs_vmanifest::i08_draft(),
        fs_vmanifest::i12_draft(),
        fs_vmanifest::i15_draft(),
    ] {
        let versions = draft.explicits.versions.to_ascii_lowercase();
        assert!(
            versions.contains("fs-vmanifest schema v2"),
            "{} does not pin canonical schema v2: {versions}",
            draft.initiative
        );
        let fields: Vec<_> = versions.split(';').map(str::trim).collect();
        for forbidden in [
            format!("manifest version {}", draft.version),
            format!("manifest revision {}", draft.version),
            format!("manifest instance version {}", draft.version),
            format!("manifest v{}", draft.version),
        ] {
            assert!(
                !fields.contains(&forbidden.as_str()),
                "{} uses known legacy numeric-revision spelling '{forbidden}' in FiveExplicits: {versions}",
                draft.initiative
            );
        }
    }
}

// The gate battery is one long deliberate sequence: each probe documents
// the order relative to its neighbors.
#[allow(clippy::too_many_lines)]
#[test]
fn freeze_gates_fire_in_the_documented_order() {
    // Cap check precedes every deep scan: an over-cap draft full of blank
    // fields refuses on the cap, not on the blanks.
    let mut over_cap = probe_draft();
    let mut blank_claim = PROBE_CLAIM;
    blank_claim.statement = "";
    over_cap.claims = vec![blank_claim; MAX_CLAIMS + 1];
    assert!(matches!(
        over_cap.freeze(),
        Err(FreezeRefusal::OverCap {
            what: "claims",
            len,
            cap,
        }) if len == MAX_CLAIMS + 1 && cap == MAX_CLAIMS
    ));

    let mut fixture_over_cap = probe_draft();
    fixture_over_cap.version = 0;
    fixture_over_cap.fixtures = vec![PROBE_FIXTURE; MAX_FIXTURES + 1];
    assert!(matches!(
        fixture_over_cap.freeze(),
        Err(FreezeRefusal::OverCap {
            what: "fixtures",
            len,
            cap,
        }) if len == MAX_FIXTURES + 1 && cap == MAX_FIXTURES
    ));

    let mut obligation_over_cap = probe_draft();
    obligation_over_cap.version = 0;
    obligation_over_cap.obligations = vec![PROBE_OBLIGATION; MAX_OBLIGATIONS + 1];
    assert!(matches!(
        obligation_over_cap.freeze(),
        Err(FreezeRefusal::OverCap {
            what: "obligations",
            len,
            cap,
        }) if len == MAX_OBLIGATIONS + 1 && cap == MAX_OBLIGATIONS
    ));

    let mut waiver_over_cap = probe_draft();
    waiver_over_cap.version = 0;
    waiver_over_cap.waivers = vec![probe_waiver("probe-claim"); MAX_WAIVERS + 1];
    assert!(matches!(
        waiver_over_cap.freeze(),
        Err(FreezeRefusal::OverCap {
            what: "waivers",
            len,
            cap,
        }) if len == MAX_WAIVERS + 1 && cap == MAX_WAIVERS
    ));

    let mut list_over_cap_and_zero_version = probe_draft();
    list_over_cap_and_zero_version.version = 0;
    let mut claim = PROBE_CLAIM;
    claim.hypotheses = vec!["bounded"; MAX_ROW_ITEMS + 1].leak();
    list_over_cap_and_zero_version.claims = vec![claim];
    assert!(matches!(
        list_over_cap_and_zero_version.freeze(),
        Err(FreezeRefusal::OverCap {
            what: "claim.hypotheses",
            len,
            cap,
        }) if len == MAX_ROW_ITEMS + 1 && cap == MAX_ROW_ITEMS
    ));

    let overlong_row: &'static [&'static str] =
        Box::leak(vec!["bounded"; MAX_ROW_ITEMS + 1].into_boxed_slice());
    for field in [
        "obligation.claims_covered",
        "obligation.unit_cases",
        "obligation.decks",
        "obligation.g3_relations",
        "obligation.obs_events",
    ] {
        let mut draft = probe_draft();
        draft.version = 0;
        match field {
            "obligation.claims_covered" => draft.obligations[0].claims_covered = overlong_row,
            "obligation.unit_cases" => draft.obligations[0].unit_cases = overlong_row,
            "obligation.decks" => draft.obligations[0].decks = overlong_row,
            "obligation.g3_relations" => draft.obligations[0].g3_relations = overlong_row,
            "obligation.obs_events" => draft.obligations[0].obs_events = overlong_row,
            _ => unreachable!("closed cap-test field table"),
        }
        assert!(
            matches!(
                draft.freeze(),
                Err(FreezeRefusal::OverCap { what, len, cap })
                    if what == field && len == MAX_ROW_ITEMS + 1 && cap == MAX_ROW_ITEMS
            ),
            "missing cap branch for {field}"
        );
    }

    // The cumulative byte cap is part of the same first gate. Reusing one
    // leaked 1 MiB slice in 17 hypothesis slots exercises cumulative text-byte
    // accounting without allocating a 17 MiB fixture solely for this test.
    let one_mebibyte: &'static str = Box::leak("x".repeat(1 << 20).into_boxed_str());
    let oversized_hypotheses: &'static [&'static str] =
        Box::leak(vec![one_mebibyte; 17].into_boxed_slice());
    let mut text_over_cap_and_zero_version = probe_draft();
    text_over_cap_and_zero_version.version = 0;
    text_over_cap_and_zero_version.claims[0].hypotheses = oversized_hypotheses;
    assert!(matches!(
        text_over_cap_and_zero_version.freeze(),
        Err(FreezeRefusal::OverCap {
            what: "manifest text bytes",
            len,
            cap,
        }) if len > cap && cap == MAX_MANIFEST_TEXT_BYTES
    ));

    let mut zero_version = probe_draft();
    zero_version.version = 0;
    assert!(matches!(
        zero_version.freeze(),
        Err(FreezeRefusal::ZeroVersion)
    ));

    let mut blank_title = probe_draft();
    blank_title.title = "  ";
    assert!(matches!(
        blank_title.freeze(),
        Err(FreezeRefusal::BlankField { field: "title", .. })
    ));

    let mut no_claims = probe_draft();
    no_claims.claims.clear();
    assert!(matches!(
        no_claims.freeze(),
        Err(FreezeRefusal::EmptyCollection { what: "claims" })
    ));

    let mut blank_statement = probe_draft();
    let mut claim = PROBE_CLAIM;
    claim.statement = " ";
    blank_statement.claims = vec![claim];
    assert!(matches!(
        blank_statement.freeze(),
        Err(FreezeRefusal::BlankField {
            field: "claim.statement",
            ..
        })
    ));

    let mut no_hypotheses = probe_draft();
    let mut claim = PROBE_CLAIM;
    claim.hypotheses = &[];
    no_hypotheses.claims = vec![claim];
    assert!(matches!(
        no_hypotheses.freeze(),
        Err(FreezeRefusal::BlankField {
            field: "claim.hypotheses",
            ..
        })
    ));

    let blank_list_cases: [(&str, fn(&mut ObligationRow)); 5] = [
        ("obligation.claims_covered", |row| {
            row.claims_covered = &[" "];
        }),
        ("obligation.unit_cases", |row| {
            row.unit_cases = &[" "];
        }),
        ("obligation.decks", |row| {
            row.decks = &[" "];
        }),
        ("obligation.g3_relations", |row| {
            row.g3_relations = &[" "];
        }),
        ("obligation.obs_events", |row| {
            row.obs_events = &[" "];
        }),
    ];
    for (expected_field, mutate) in blank_list_cases {
        let mut draft = probe_draft();
        mutate(&mut draft.obligations[0]);
        assert!(matches!(
            draft.freeze(),
            Err(FreezeRefusal::BlankField { field, .. }) if field == expected_field
        ));
    }

    let empty_list_cases: [(&str, fn(&mut ObligationRow)); 5] = [
        ("obligation.claims_covered", |row| {
            row.claims_covered = &[];
        }),
        ("obligation.unit_cases", |row| {
            row.unit_cases = &[];
        }),
        ("obligation.decks", |row| {
            row.decks = &[];
        }),
        ("obligation.g3_relations", |row| {
            row.g3_relations = &[];
        }),
        ("obligation.obs_events", |row| {
            row.obs_events = &[];
        }),
    ];
    for (expected_field, mutate) in empty_list_cases {
        let mut draft = probe_draft();
        mutate(&mut draft.obligations[0]);
        assert!(matches!(
            draft.freeze(),
            Err(FreezeRefusal::BlankField { field, .. }) if field == expected_field
        ));
    }

    let mut duplicate = probe_draft();
    duplicate.claims = vec![PROBE_CLAIM, PROBE_CLAIM];
    match duplicate.freeze() {
        Err(FreezeRefusal::DuplicateId { kind: "claim", id }) => {
            assert_eq!(id, "probe-claim");
        }
        other => panic!("duplicate claim must fail closed, got {other:?}"),
    }

    let mut ambiguous_evidence_id = probe_draft();
    ambiguous_evidence_id.obligations[0].leaf = "probe-claim";
    match ambiguous_evidence_id.freeze() {
        Err(FreezeRefusal::DuplicateId {
            kind: "claim/obligation evidence",
            id,
        }) => assert_eq!(id, "probe-claim"),
        other => panic!("ambiguous evidence id must fail closed, got {other:?}"),
    }

    let mut duplicate_mapping = probe_draft();
    duplicate_mapping.obligations[0].claims_covered = &["probe-claim", "probe-claim"];
    match duplicate_mapping.freeze() {
        Err(FreezeRefusal::DuplicateId {
            kind: "obligation claim mapping",
            id,
        }) => assert_eq!(id, "probe-claim"),
        other => panic!("duplicate claim mapping must fail closed, got {other:?}"),
    }

    let duplicate_set_cases: [(&str, fn(&mut ObligationRow)); 4] = [
        ("obligation unit case", |row| {
            row.unit_cases = &["happy", "happy"];
        }),
        ("obligation deck", |row| {
            row.decks = &["probe-deck", "probe-deck"];
        }),
        ("obligation G3 relation", |row| {
            row.g3_relations = &["same", "same"];
        }),
        ("obligation observation event", |row| {
            row.obs_events = &["probe.event", "probe.event"];
        }),
    ];
    for (expected_kind, mutate) in duplicate_set_cases {
        let mut draft = probe_draft();
        mutate(&mut draft.obligations[0]);
        match draft.freeze() {
            Err(FreezeRefusal::DuplicateId { kind, id }) => {
                assert_eq!(kind, expected_kind);
                assert!(!id.is_empty());
            }
            other => panic!("duplicate {expected_kind} must fail closed, got {other:?}"),
        }
    }

    // Reusing the production path as its own oracle fails closed.
    let mut reused_oracle = probe_draft();
    let mut claim = PROBE_CLAIM;
    claim.oracle = OracleRoute {
        identity: "the production path itself",
        independent: false,
        tcb_overlap: "total",
    };
    reused_oracle.claims = vec![claim];
    match reused_oracle.freeze() {
        Err(FreezeRefusal::ProductionOracleReuse { claim }) => {
            assert_eq!(claim, "probe-claim");
        }
        other => panic!("oracle reuse must fail closed, got {other:?}"),
    }

    // Oracle independence is a whole-manifest gate: a later production-path
    // reuse must beat an earlier claim's invalid tolerance.
    let mut mixed_claim_failures = shared_probe_draft();
    mixed_claim_failures.claims[0].tolerance = ToleranceSemantics::Absolute { atol: 0.0 };
    mixed_claim_failures.claims[1].oracle.independent = false;
    assert!(matches!(
        mixed_claim_failures.clone().freeze(),
        Err(FreezeRefusal::ProductionOracleReuse { claim }) if claim == "probe-sibling"
    ));
    mixed_claim_failures.claims.reverse();
    assert!(matches!(
        mixed_claim_failures.freeze(),
        Err(FreezeRefusal::ProductionOracleReuse { claim }) if claim == "probe-sibling"
    ));

    let tolerance_cases: [(ToleranceSemantics, &str); 4] = [
        (
            ToleranceSemantics::Absolute { atol: 0.0 },
            "absolute tolerance must be finite and > 0",
        ),
        (
            ToleranceSemantics::Relative { rtol: f64::NAN },
            "relative tolerance must be finite and > 0",
        ),
        (
            ToleranceSemantics::AbsRel {
                atol: 0.0,
                rtol: 0.0,
            },
            "zero-width tolerance (use Exact for bit verdicts)",
        ),
        (
            ToleranceSemantics::Interval { lo: 2.0, hi: 1.0 },
            "inverted interval",
        ),
    ];
    for (tolerance, expected) in tolerance_cases {
        let mut draft = probe_draft();
        let mut claim = PROBE_CLAIM;
        claim.tolerance = tolerance;
        draft.claims = vec![claim];
        match draft.freeze() {
            Err(FreezeRefusal::InvalidTolerance { reason, .. }) => {
                assert_eq!(reason, expected);
            }
            other => panic!("tolerance {tolerance:?} must refuse, got {other:?}"),
        }
    }

    let mut blank_spec = probe_draft();
    blank_spec.fixtures = vec![FixturePin {
        id: "probe-deck",
        source: FixtureSource::AuthoredSpec { spec: "   " },
        partition: Partition::Development,
    }];
    assert!(matches!(
        blank_spec.freeze(),
        Err(FreezeRefusal::MalformedFixture { .. })
    ));

    let mut bad_hex = probe_draft();
    bad_hex.fixtures = vec![FixturePin {
        id: "probe-deck",
        source: FixtureSource::External {
            digest_hex: "not-hex",
        },
        partition: Partition::Development,
    }];
    assert!(matches!(
        bad_hex.freeze(),
        Err(FreezeRefusal::MalformedFixture { .. })
    ));

    let mut orphan_claim_ref = probe_draft();
    let mut row = PROBE_OBLIGATION;
    row.claims_covered = &["no-such-claim"];
    orphan_claim_ref.obligations = vec![row];
    match orphan_claim_ref.freeze() {
        // The probe claim also becomes uncovered, but the orphan reference
        // gate fires first (documented order).
        Err(FreezeRefusal::OrphanClaimRef { claim, .. }) => {
            assert_eq!(claim, "no-such-claim");
        }
        other => panic!("orphan claim ref must fail closed, got {other:?}"),
    }

    let mut orphan_deck = probe_draft();
    let mut row = PROBE_OBLIGATION;
    row.decks = &["no-such-deck"];
    orphan_deck.obligations = vec![row];
    match orphan_deck.freeze() {
        Err(FreezeRefusal::OrphanDeck { deck, .. }) => assert_eq!(deck, "no-such-deck"),
        other => panic!("orphan deck must fail closed, got {other:?}"),
    }

    // Claim-reference resolution is a whole-manifest gate and therefore
    // precedes deck resolution even when the orphan deck occurs in row zero.
    let mut mixed_reference_failures = shared_probe_draft();
    mixed_reference_failures.obligations[0].claims_covered = &["probe-claim"];
    mixed_reference_failures.obligations[0].decks = &["no-such-deck"];
    let mut later_orphan_claim = PROBE_OBLIGATION;
    later_orphan_claim.leaf = "probe-second-leaf";
    later_orphan_claim.claims_covered = &["no-such-claim"];
    mixed_reference_failures
        .obligations
        .push(later_orphan_claim);
    assert!(matches!(
        mixed_reference_failures.clone().freeze(),
        Err(FreezeRefusal::OrphanClaimRef { leaf, claim })
            if leaf == "probe-second-leaf" && claim == "no-such-claim"
    ));
    mixed_reference_failures.obligations.reverse();
    assert!(matches!(
        mixed_reference_failures.freeze(),
        Err(FreezeRefusal::OrphanClaimRef { leaf, claim })
            if leaf == "probe-second-leaf" && claim == "no-such-claim"
    ));

    let mut uncovered = probe_draft();
    uncovered.obligations = Vec::new();
    match uncovered.freeze() {
        Err(FreezeRefusal::UncoveredClaim { claim }) => assert_eq!(claim, "probe-claim"),
        other => panic!("uncovered claim must fail closed, got {other:?}"),
    }

    // A waiver is an acceptable coverage route — but only a complete one.
    let mut waived = probe_draft();
    waived.obligations = Vec::new();
    waived.waivers = vec![Waiver {
        subject: "probe-claim",
        reason: "obligation lands with the first implementation bead",
        owner: "probe owner",
        predicate: "obligation row added via amendment",
        expiry: "next close burst",
        promotion_effect: "claim cannot promote while waived",
    }];
    assert!(waived.freeze().is_ok());

    let mut orphan_waiver = probe_draft();
    orphan_waiver.waivers = vec![probe_waiver("unused-waiver-subject")];
    match orphan_waiver.freeze() {
        Err(FreezeRefusal::OrphanWaiver { subject }) => {
            assert_eq!(subject, "unused-waiver-subject");
        }
        other => panic!("orphan waiver must fail closed, got {other:?}"),
    }

    let mut ambiguous_waiver = shared_probe_draft();
    ambiguous_waiver.obligations[0].decks = &["probe-deck", "probe-claim"];
    ambiguous_waiver.waivers = vec![probe_waiver("probe-claim")];
    match ambiguous_waiver.freeze() {
        Err(FreezeRefusal::AmbiguousWaiverSubject { subject }) => {
            assert_eq!(subject, "probe-claim");
        }
        other => panic!("ambiguous waiver must fail closed, got {other:?}"),
    }
}

#[test]
fn canonical_identity_is_input_order_invariant() {
    // G5: assembly order cannot move the digest.
    let baseline = i01_draft().freeze().expect("freeze").digest();
    let mut reordered = i01_draft();
    reordered.claims.reverse();
    reordered.fixtures.reverse();
    reordered.obligations.reverse();
    reordered.waivers.reverse();
    let reordered_digest = reordered.freeze().expect("freeze").digest();
    assert_eq!(reordered_digest, baseline);
}

#[test]
fn chunked_in_memory_assembly_matches_one_shot() {
    // G4 precursor at this layer: a draft assembled in chunks freezes to the
    // same identity as one-shot assembly. This does not encode a checkpoint,
    // restart a process, detect corruption, or exercise cancellation.
    let one_shot = i01_draft();
    let mut staged = ManifestDraft {
        initiative: one_shot.initiative,
        title: one_shot.title,
        version: one_shot.version,
        explicits: one_shot.explicits,
        claims: Vec::new(),
        fixtures: Vec::new(),
        obligations: Vec::new(),
        waivers: Vec::new(),
        amendment_rules: one_shot.amendment_rules,
    };
    let split = one_shot.claims.len() / 2;
    staged.claims.extend_from_slice(&one_shot.claims[..split]);
    // ... interruption boundary: everything above could be dropped and
    // rebuilt; the digest depends only on frozen content ...
    staged.claims.extend_from_slice(&one_shot.claims[split..]);
    staged.fixtures.extend_from_slice(&one_shot.fixtures);
    staged.obligations.extend_from_slice(&one_shot.obligations);
    staged.waivers.extend_from_slice(&one_shot.waivers);
    let a = one_shot.freeze().expect("freeze").digest();
    let b = staged.freeze().expect("freeze").digest();
    assert_eq!(a, b);
}

#[test]
fn g3_mutations_fail_closed_or_move_the_identity() {
    let baseline = i01_draft().freeze().expect("freeze").digest();

    // Loosen a band: identity must move.
    let mut loosened = i01_draft();
    for claim in &mut loosened.claims {
        if claim.id == "i01-residual-jvp-vjp-actions" {
            claim.tolerance = ToleranceSemantics::AbsRel {
                atol: 1e-6,
                rtol: 1e-3,
            };
        }
    }
    assert_ne!(loosened.freeze().expect("freeze").digest(), baseline);

    // Weaken hypotheses: identity must move.
    let mut weakened = i01_draft();
    for claim in &mut weakened.claims {
        if claim.id == "i01-residual-jvp-vjp-actions" {
            claim.hypotheses = &["operands drawn from the pinned coupled MMS deck family"];
        }
    }
    assert_ne!(weakened.freeze().expect("freeze").digest(), baseline);

    // Swap a held-out fixture: identity must move.
    let mut swapped = i01_draft();
    for fixture in &mut swapped.fixtures {
        if fixture.id == "i01-coupled-mms-holdout" {
            fixture.source = FixtureSource::AuthoredSpec {
                spec: "swapped generator text pretending to be the held-out family",
            };
        }
    }
    assert_ne!(swapped.freeze().expect("freeze").digest(), baseline);

    // Reuse a production oracle: fails closed (not merely a new identity).
    let mut reused = i01_draft();
    for claim in &mut reused.claims {
        if claim.id == "i01-conservation-exact-sequence" {
            claim.oracle = OracleRoute {
                identity: "the production assembly path",
                independent: false,
                tcb_overlap: "total",
            };
        }
    }
    assert!(matches!(
        reused.freeze(),
        Err(FreezeRefusal::ProductionOracleReuse { .. })
    ));

    // Altering flags AFTER freeze has no code path at all: FrozenManifest
    // exposes no mutating API, so the "alter after freeze" mutation is
    // impossible by construction; change goes through amend() below.
}

#[test]
#[allow(clippy::too_many_lines)]
fn per_component_identity_is_mutation_sensitive() {
    // Exhaustive destructuring is a compile-time maintenance tripwire: adding
    // a semantic field requires this battery to make an explicit decision.
    let ClaimSpec {
        id: _,
        ambition: _,
        polarity: _,
        statement: _,
        hypotheses: _,
        qoi: _,
        unit: _,
        tolerance: _,
        evidence_tier: _,
        oracle: _,
        activation: _,
        kill: _,
        fallback: _,
        no_claim: _,
    } = PROBE_CLAIM;
    let OracleRoute {
        identity: _,
        independent: _,
        tcb_overlap: _,
    } = PROBE_CLAIM.oracle;
    let FixturePin {
        id: _,
        source: _,
        partition: _,
    } = PROBE_FIXTURE;
    let ObligationRow {
        leaf: _,
        claims_covered: _,
        unit_cases: _,
        g0: _,
        decks: _,
        g3_relations: _,
        g4_schedule: _,
        g5_matrix: _,
        entry_point: _,
        tier: _,
        dsr_lane: _,
        obs_events: _,
        replay_command: _,
    } = PROBE_OBLIGATION;
    let waiver = probe_waiver("probe-claim");
    let Waiver {
        subject: _,
        reason: _,
        owner: _,
        predicate: _,
        expiry: _,
        promotion_effect: _,
    } = waiver;

    macro_rules! claim_field_changes_digest {
        ($field:ident, $value:expr) => {{
            let mut candidate = PROBE_CLAIM;
            candidate.$field = $value;
            assert_ne!(candidate, PROBE_CLAIM, stringify!($field));
            assert_ne!(
                claim_digest(&candidate),
                claim_digest(&PROBE_CLAIM),
                stringify!($field)
            );
        }};
    }
    claim_field_changes_digest!(id, "different-claim-id");
    claim_field_changes_digest!(ambition, Ambition::Moonshot);
    claim_field_changes_digest!(polarity, ClaimPolarity::Refutation);
    claim_field_changes_digest!(statement, "different statement");
    claim_field_changes_digest!(hypotheses, &["different hypothesis"]);
    claim_field_changes_digest!(qoi, "different_qoi");
    claim_field_changes_digest!(unit, "J");
    claim_field_changes_digest!(tolerance, ToleranceSemantics::Absolute { atol: 0.25 });
    claim_field_changes_digest!(evidence_tier, GauntletTier::G5);
    claim_field_changes_digest!(activation, "different activation");
    claim_field_changes_digest!(kill, "different kill");
    claim_field_changes_digest!(fallback, "different fallback");
    claim_field_changes_digest!(no_claim, "different no-claim boundary");

    for (field, oracle) in [
        (
            "oracle.identity",
            OracleRoute {
                identity: "different independent oracle",
                ..PROBE_CLAIM.oracle
            },
        ),
        (
            "oracle.independent",
            OracleRoute {
                independent: false,
                ..PROBE_CLAIM.oracle
            },
        ),
        (
            "oracle.tcb_overlap",
            OracleRoute {
                tcb_overlap: "different overlap",
                ..PROBE_CLAIM.oracle
            },
        ),
    ] {
        let mut candidate = PROBE_CLAIM;
        candidate.oracle = oracle;
        assert_ne!(candidate, PROBE_CLAIM, "{field}");
        assert_ne!(
            claim_digest(&candidate),
            claim_digest(&PROBE_CLAIM),
            "{field}"
        );
    }

    macro_rules! fixture_field_changes_digest {
        ($field:ident, $value:expr) => {{
            let mut candidate = PROBE_FIXTURE;
            candidate.$field = $value;
            assert_ne!(candidate, PROBE_FIXTURE, stringify!($field));
            assert_ne!(
                fixture_digest(&candidate),
                fixture_digest(&PROBE_FIXTURE),
                stringify!($field)
            );
        }};
    }
    fixture_field_changes_digest!(id, "different-fixture-id");
    fixture_field_changes_digest!(
        source,
        FixtureSource::AuthoredSpec {
            spec: "different fixture bytes"
        }
    );
    fixture_field_changes_digest!(
        source,
        FixtureSource::External {
            digest_hex: "0000000000000000000000000000000000000000000000000000000000000001"
        }
    );
    fixture_field_changes_digest!(partition, Partition::HeldOut);

    macro_rules! obligation_field_changes_digest {
        ($field:ident, $value:expr) => {{
            let mut candidate = PROBE_OBLIGATION;
            candidate.$field = $value;
            assert_ne!(candidate, PROBE_OBLIGATION, stringify!($field));
            assert_ne!(
                obligation_digest(&candidate),
                obligation_digest(&PROBE_OBLIGATION),
                stringify!($field)
            );
        }};
    }
    obligation_field_changes_digest!(leaf, "different-leaf");
    obligation_field_changes_digest!(claims_covered, &["different-claim"]);
    obligation_field_changes_digest!(unit_cases, &["different-case"]);
    obligation_field_changes_digest!(g0, "different G0 contract");
    obligation_field_changes_digest!(decks, &["different-deck"]);
    obligation_field_changes_digest!(g3_relations, &["different G3 relation"]);
    obligation_field_changes_digest!(g4_schedule, "different G4 schedule");
    obligation_field_changes_digest!(g5_matrix, "different G5 matrix");
    obligation_field_changes_digest!(entry_point, "different-entry-point");
    obligation_field_changes_digest!(tier, CampaignTier::Max);
    obligation_field_changes_digest!(dsr_lane, "different DSR lane");
    obligation_field_changes_digest!(obs_events, &["different.event"]);
    obligation_field_changes_digest!(replay_command, "different replay command");

    macro_rules! waiver_field_changes_digest {
        ($field:ident, $value:expr) => {{
            let mut candidate = waiver;
            candidate.$field = $value;
            assert_ne!(candidate, waiver, stringify!($field));
            assert_ne!(
                waiver_digest(&candidate),
                waiver_digest(&waiver),
                stringify!($field)
            );
        }};
    }
    waiver_field_changes_digest!(subject, "different-subject");
    waiver_field_changes_digest!(reason, "different reason");
    waiver_field_changes_digest!(owner, "different owner");
    waiver_field_changes_digest!(predicate, "different predicate");
    waiver_field_changes_digest!(expiry, "different expiry");
    waiver_field_changes_digest!(promotion_effect, "different effect");

    // External hex case is normalized: one canonical identity.
    let upper = FixturePin {
        id: "x",
        source: FixtureSource::External {
            digest_hex: "ABCDEFABCDEFABCDEFABCDEFABCDEFABCDEFABCDEFABCDEFABCDEFABCDEF1234",
        },
        partition: Partition::Development,
    };
    let lower = FixturePin {
        id: "x",
        source: FixtureSource::External {
            digest_hex: "abcdefabcdefabcdefabcdefabcdefabcdefabcdefabcdefabcdefabcdef1234",
        },
        partition: Partition::Development,
    };
    assert_eq!(upper.source, lower.source);
    assert_eq!(upper, lower);
    assert_eq!(fixture_digest(&upper), fixture_digest(&lower));

    for (field, left_tolerance, right_tolerance) in [
        (
            "Absolute.atol",
            ToleranceSemantics::Absolute { atol: 1.0 },
            ToleranceSemantics::Absolute { atol: 2.0 },
        ),
        (
            "Relative.rtol",
            ToleranceSemantics::Relative { rtol: 1.0 },
            ToleranceSemantics::Relative { rtol: 2.0 },
        ),
        (
            "AbsRel.atol",
            ToleranceSemantics::AbsRel {
                atol: 1.0,
                rtol: 2.0,
            },
            ToleranceSemantics::AbsRel {
                atol: 1.5,
                rtol: 2.0,
            },
        ),
        (
            "AbsRel.rtol",
            ToleranceSemantics::AbsRel {
                atol: 1.0,
                rtol: 2.0,
            },
            ToleranceSemantics::AbsRel {
                atol: 1.0,
                rtol: 2.5,
            },
        ),
        (
            "Interval.lo",
            ToleranceSemantics::Interval { lo: 1.0, hi: 2.0 },
            ToleranceSemantics::Interval { lo: 1.5, hi: 2.0 },
        ),
        (
            "Interval.hi",
            ToleranceSemantics::Interval { lo: 1.0, hi: 2.0 },
            ToleranceSemantics::Interval { lo: 1.0, hi: 2.5 },
        ),
        (
            "variant tag",
            ToleranceSemantics::Exact,
            ToleranceSemantics::Absolute { atol: 1.0 },
        ),
    ] {
        let mut left = PROBE_CLAIM;
        left.tolerance = left_tolerance;
        let mut right = PROBE_CLAIM;
        right.tolerance = right_tolerance;
        assert_ne!(left_tolerance, right_tolerance, "{field}");
        assert_ne!(claim_digest(&left), claim_digest(&right), "{field}");
    }
    let tag_cases = [
        ToleranceSemantics::Exact,
        ToleranceSemantics::Absolute { atol: 1.0 },
        ToleranceSemantics::Relative { rtol: 1.0 },
        ToleranceSemantics::AbsRel {
            atol: 1.0,
            rtol: 1.0,
        },
        ToleranceSemantics::Interval { lo: 1.0, hi: 1.0 },
    ];
    for (left_index, &left_tolerance) in tag_cases.iter().enumerate() {
        for &right_tolerance in &tag_cases[left_index + 1..] {
            let mut left = PROBE_CLAIM;
            left.tolerance = left_tolerance;
            let mut right = PROBE_CLAIM;
            right.tolerance = right_tolerance;
            assert_ne!(left_tolerance, right_tolerance, "variant tags");
            assert_ne!(claim_digest(&left), claim_digest(&right), "variant tags");
        }
    }

    // Valid signed zero bounds remain distinct authored identities because
    // canonical tolerance bytes use exact IEEE-754 encodings.
    let mut negative_zero = PROBE_CLAIM;
    negative_zero.tolerance = ToleranceSemantics::Interval { lo: -0.0, hi: 1.0 };
    let mut positive_zero = negative_zero;
    positive_zero.tolerance = ToleranceSemantics::Interval { lo: 0.0, hi: 1.0 };
    assert_ne!(negative_zero.tolerance, positive_zero.tolerance);
    assert_ne!(negative_zero, positive_zero);
    assert_ne!(claim_digest(&negative_zero), claim_digest(&positive_zero));
}

#[test]
fn canonical_set_and_external_hex_presentations_share_identity_and_equality() {
    let mut baseline = shared_probe_draft();
    baseline.fixtures.push(FixturePin {
        id: "probe-deck-2",
        source: FixtureSource::AuthoredSpec {
            spec: "PROBE: second canonical deck.",
        },
        partition: Partition::Development,
    });
    baseline.obligations[0].decks = &["probe-deck", "probe-deck-2"];
    baseline.obligations[0].g3_relations =
        &["probe relabeling invariance", "probe scaling invariance"];
    baseline.obligations[0].obs_events = &["probe.event", "probe.receipt"];

    let mut reordered = baseline.clone();
    reordered.obligations[0].claims_covered = &["probe-sibling", "probe-claim"];
    reordered.obligations[0].unit_cases = &["error", "empty", "happy"];
    reordered.obligations[0].decks = &["probe-deck-2", "probe-deck"];
    reordered.obligations[0].g3_relations =
        &["probe scaling invariance", "probe relabeling invariance"];
    reordered.obligations[0].obs_events = &["probe.receipt", "probe.event"];

    for row in [&baseline.obligations[0], &reordered.obligations[0]] {
        assert_eq!(
            row.canonical_claims_covered(),
            vec!["probe-claim", "probe-sibling"]
        );
        assert_eq!(row.canonical_unit_cases(), vec!["empty", "error", "happy"]);
        assert_eq!(row.canonical_decks(), vec!["probe-deck", "probe-deck-2"]);
        assert_eq!(
            row.canonical_g3_relations(),
            vec!["probe relabeling invariance", "probe scaling invariance"]
        );
        assert_eq!(
            row.canonical_obs_events(),
            vec!["probe.event", "probe.receipt"]
        );
    }

    assert_eq!(baseline.obligations[0], reordered.obligations[0]);
    assert_eq!(
        obligation_digest(&baseline.obligations[0]),
        obligation_digest(&reordered.obligations[0])
    );
    let predecessor = baseline.clone().freeze().expect("baseline freeze");
    let equivalent = reordered.clone().freeze().expect("reordered-set freeze");
    assert_eq!(predecessor.digest(), equivalent.digest());
    assert_eq!(predecessor, equivalent);
    assert_eq!(predecessor.obligations(), equivalent.obligations());
    let canonical = &predecessor.obligations()[0];
    assert_eq!(
        canonical.digest(),
        obligation_digest(&baseline.obligations[0])
    );
    assert_eq!(
        canonical.claims_covered(),
        &["probe-claim", "probe-sibling"]
    );
    assert_eq!(canonical.unit_cases(), &["empty", "error", "happy"]);
    assert_eq!(canonical.decks(), &["probe-deck", "probe-deck-2"]);
    assert_eq!(
        canonical.g3_relations(),
        &["probe relabeling invariance", "probe scaling invariance"]
    );
    assert_eq!(canonical.obs_events(), &["probe.event", "probe.receipt"]);

    reordered.version = 2;
    let (_, record) = predecessor
        .amend(reordered)
        .expect("presentation-only set reorder");
    assert!(record.invalidated.is_empty());

    const UPPER: &str = "ABCDEFABCDEFABCDEFABCDEFABCDEFABCDEFABCDEFABCDEFABCDEFABCDEF1234";
    const LOWER: &str = "abcdefabcdefabcdefabcdefabcdefabcdefabcdefabcdefabcdefabcdef1234";
    let mut upper = probe_draft();
    upper.fixtures[0].source = FixtureSource::External { digest_hex: UPPER };
    let mut lower = probe_draft();
    lower.fixtures[0].source = FixtureSource::External { digest_hex: LOWER };
    let upper = upper.freeze().expect("uppercase external digest");
    let lower = lower.freeze().expect("lowercase external digest");
    assert_eq!(upper.digest(), lower.digest());
    assert_eq!(upper, lower);
}

#[test]
fn amendment_refuses_identity_changes_version_errors_and_defective_successors() {
    let frozen = i01_draft().freeze().expect("freeze");

    let mut other_initiative = i01_draft();
    other_initiative.initiative = "I99";
    other_initiative.version = 2;
    match frozen.amend(other_initiative) {
        Err(AmendmentRefusal::InitiativeChanged { expected, offered }) => {
            assert_eq!(expected, "I01");
            assert_eq!(offered, "I99");
        }
        other => panic!("cross-initiative amendment must refuse, got {other:?}"),
    }

    // Wrong version: refused.
    let same_version = i01_draft();
    assert!(matches!(
        frozen.amend(same_version),
        Err(AmendmentRefusal::WrongVersion {
            expected: 2,
            offered: 1
        })
    ));

    // A defective successor is refused with its own gate.
    let mut broken = i01_draft();
    broken.version = 2;
    broken.title = " ";
    assert!(matches!(
        frozen.amend(broken),
        Err(AmendmentRefusal::SuccessorRefused(
            FreezeRefusal::BlankField { field: "title", .. }
        ))
    ));

    let mut max_version = probe_draft();
    max_version.version = u32::MAX;
    let max_version = max_version.freeze().expect("maximum version freezes");
    assert!(matches!(
        max_version.amend(probe_draft()),
        Err(AmendmentRefusal::VersionExhausted { version: u32::MAX })
    ));
}

#[test]
fn amendment_refuses_cross_version_evidence_kind_aliases_in_both_directions() {
    let predecessor = probe_draft().freeze().expect("predecessor freeze");

    let mut claim_to_leaf = probe_draft();
    claim_to_leaf.version = 2;
    claim_to_leaf.claims[0].id = "probe-renamed-claim";
    claim_to_leaf.obligations[0].leaf = "probe-claim";
    claim_to_leaf.obligations[0].claims_covered = &["probe-renamed-claim"];
    assert!(matches!(
        predecessor.amend(claim_to_leaf),
        Err(AmendmentRefusal::EvidenceKindChanged {
            id,
            from_kind: "claim",
            to_kind: "obligation leaf",
        }) if id == "probe-claim"
    ));

    let mut leaf_to_claim = probe_draft();
    leaf_to_claim.version = 2;
    leaf_to_claim.claims[0].id = "probe-leaf";
    leaf_to_claim.obligations[0].leaf = "probe-renamed-leaf";
    leaf_to_claim.obligations[0].claims_covered = &["probe-leaf"];
    assert!(matches!(
        predecessor.amend(leaf_to_claim),
        Err(AmendmentRefusal::EvidenceKindChanged {
            id,
            from_kind: "obligation leaf",
            to_kind: "claim",
        }) if id == "probe-leaf"
    ));
}

#[test]
fn version_only_and_new_authority_do_not_revoke_unchanged_predecessor_evidence() {
    let frozen = probe_draft().freeze().expect("freeze");

    let mut version_only = probe_draft();
    version_only.version = 2;
    let (successor, version_record) = frozen.amend(version_only).expect("version-only amendment");
    assert!(version_record.invalidated.is_empty());
    assert_eq!(
        (version_record.from_version, version_record.to_version),
        (1, 2)
    );
    assert_eq!(version_record.from_digest, frozen.digest());
    assert_eq!(version_record.to_digest, successor.digest());
    assert_eq!(successor.version(), 2);
    assert_ne!(version_record.from_digest, version_record.to_digest);

    let mut added_fixture = probe_draft();
    added_fixture.version = 2;
    let mut unrelated = PROBE_FIXTURE;
    unrelated.id = "probe-unreferenced-deck";
    unrelated.source = FixtureSource::AuthoredSpec {
        spec: "PROBE: unrelated future deck bytes.",
    };
    added_fixture.fixtures.push(unrelated);
    let (_, fixture_record) = frozen
        .amend(added_fixture)
        .expect("unrelated fixture addition");
    assert!(fixture_record.invalidated.is_empty());

    let mut added_authority = probe_draft();
    added_authority.version = 2;
    let mut new_claim = PROBE_CLAIM;
    new_claim.id = "probe-new-claim";
    new_claim.statement = "new claim with no predecessor authority";
    new_claim.qoi = "probe_new_qoi";
    added_authority.claims.push(new_claim);
    let mut new_leaf = PROBE_OBLIGATION;
    new_leaf.leaf = "probe-new-leaf";
    new_leaf.claims_covered = &["probe-new-claim"];
    added_authority.obligations.push(new_leaf);
    let (_, added_record) = frozen
        .amend(added_authority)
        .expect("new claim and leaf addition");
    assert!(added_record.invalidated.is_empty());

    let shared = shared_probe_draft().freeze().expect("shared freeze");
    let mut reordered_mapping = shared_probe_draft();
    reordered_mapping.version = 2;
    reordered_mapping.obligations[0].claims_covered = &["probe-sibling", "probe-claim"];
    let (_, reorder_record) = shared
        .amend(reordered_mapping)
        .expect("set-order-only mapping amendment");
    assert!(
        reorder_record.invalidated.is_empty(),
        "claim mapping presentation order is not execution semantics"
    );
}

#[test]
fn claim_changes_invalidate_the_claim_and_producer_leaf_but_not_siblings() {
    let frozen = shared_probe_draft().freeze().expect("freeze");
    let mut successor = shared_probe_draft();
    successor.version = 2;
    successor
        .claims
        .iter_mut()
        .find(|claim| claim.id == "probe-claim")
        .expect("probe claim")
        .statement = "changed probe statement";

    let (_, record) = frozen.amend(successor).expect("valid amendment");
    assert_eq!(record.invalidated, vec!["probe-claim", "probe-leaf"]);

    // Removing the claim also removes its row mapping, but the sibling
    // continues to consume byte-identical execution evidence and must not
    // be invalidated merely because the claims_covered list got shorter.
    let mut removed = shared_probe_draft();
    removed.version = 2;
    removed.claims.retain(|claim| claim.id != "probe-claim");
    removed.obligations[0].claims_covered = &["probe-sibling"];
    let (_, removed_record) = frozen.amend(removed).expect("claim removal amendment");
    assert_eq!(
        removed_record.invalidated,
        vec!["probe-claim", "probe-leaf"]
    );
}

#[test]
fn mapping_only_rewire_invalidates_removed_and_new_consumers_but_not_stable_siblings() {
    let mut predecessor = shared_probe_draft();
    let mut third = PROBE_CLAIM;
    third.id = "probe-third";
    third.statement = "third predecessor claim";
    third.qoi = "probe_third_qoi";
    predecessor.claims.push(third);
    let mut third_leaf = PROBE_OBLIGATION;
    third_leaf.leaf = "probe-third-leaf";
    third_leaf.claims_covered = &["probe-third"];
    predecessor.obligations.push(third_leaf);
    let frozen = predecessor.clone().freeze().expect("freeze");

    let mut successor = predecessor;
    successor.version = 2;
    successor.obligations[0].claims_covered = &["probe-sibling", "probe-third"];
    let (_, record) = frozen.amend(successor).expect("mapping-only rewire");
    assert_eq!(
        record.invalidated,
        vec!["probe-claim", "probe-leaf", "probe-third"]
    );
}

#[test]
fn fixture_and_obligation_changes_follow_reverse_dependencies() {
    let frozen = shared_probe_draft().freeze().expect("freeze");

    let mut changed_fixture = shared_probe_draft();
    changed_fixture.version = 2;
    changed_fixture.fixtures[0].partition = Partition::HeldOut;
    let (_, fixture_record) = frozen.amend(changed_fixture).expect("fixture amendment");
    assert_eq!(
        fixture_record.invalidated,
        vec!["probe-claim", "probe-leaf", "probe-sibling"]
    );

    // Filling a previously waived deck slot is also a semantic evidence
    // change even if the obligation and waiver bytes stay identical.
    let mut waived_slot = shared_probe_draft();
    waived_slot.fixtures.clear();
    waived_slot.waivers.push(probe_waiver("probe-deck"));
    let waived_slot = waived_slot.freeze().expect("waived deck predecessor");
    let mut filled_slot = shared_probe_draft();
    filled_slot.version = 2;
    filled_slot.waivers.push(probe_waiver("probe-deck"));
    let (_, filled_record) = waived_slot
        .amend(filled_slot)
        .expect("filled waived-deck amendment");
    assert_eq!(
        filled_record.invalidated,
        vec!["probe-claim", "probe-leaf", "probe-sibling"]
    );

    // The ordinary discharge transition adds the formerly waived deck and
    // removes its waiver in one successor. Both reverse-dependency walks hit
    // the same predecessor authority, and the set accumulator must name each
    // id exactly once.
    let mut discharged_slot = shared_probe_draft();
    discharged_slot.version = 2;
    let (_, discharged_record) = waived_slot
        .amend(discharged_slot)
        .expect("waiver discharge amendment");
    assert_eq!(
        discharged_record.invalidated,
        vec!["probe-claim", "probe-leaf", "probe-sibling"]
    );

    let mut changed_obligation = shared_probe_draft();
    changed_obligation.version = 2;
    changed_obligation.obligations[0].g0 = "changed generator and law contract";
    let (_, obligation_record) = frozen
        .amend(changed_obligation)
        .expect("obligation amendment");
    assert_eq!(
        obligation_record.invalidated,
        vec!["probe-claim", "probe-leaf", "probe-sibling"]
    );

    let mut added_obligation = shared_probe_draft();
    added_obligation.version = 2;
    let mut new_leaf = PROBE_OBLIGATION;
    new_leaf.leaf = "probe-added-leaf";
    new_leaf.claims_covered = &["probe-sibling"];
    added_obligation.obligations.push(new_leaf);
    let (_, added_record) = frozen
        .amend(added_obligation)
        .expect("added obligation amendment");
    assert_eq!(added_record.invalidated, vec!["probe-sibling"]);
}

#[test]
fn waiver_mutations_route_subjects_through_the_predecessor_graph() {
    let frozen = shared_probe_draft().freeze().expect("freeze");

    let mut added_claim_waiver = shared_probe_draft();
    added_claim_waiver.version = 2;
    added_claim_waiver.waivers.push(probe_waiver("probe-claim"));
    let (_, added_record) = frozen
        .amend(added_claim_waiver)
        .expect("added waiver amendment");
    assert_eq!(added_record.invalidated, vec!["probe-claim", "probe-leaf"]);

    let mut with_deck_waiver = shared_probe_draft();
    with_deck_waiver.waivers.push(probe_waiver("probe-deck"));
    let with_deck_waiver = with_deck_waiver.freeze().expect("waived predecessor");
    let mut changed_deck_waiver = shared_probe_draft();
    changed_deck_waiver.version = 2;
    let mut waiver = probe_waiver("probe-deck");
    waiver.reason = "changed probe waiver reason";
    changed_deck_waiver.waivers.push(waiver);
    let (_, changed_record) = with_deck_waiver
        .amend(changed_deck_waiver)
        .expect("changed waiver amendment");
    assert_eq!(
        changed_record.invalidated,
        vec!["probe-claim", "probe-leaf", "probe-sibling"]
    );

    let mut with_leaf_waiver = shared_probe_draft();
    with_leaf_waiver.waivers.push(probe_waiver("probe-leaf"));
    let with_leaf_waiver = with_leaf_waiver.freeze().expect("waived predecessor");
    let mut removed_leaf_waiver = shared_probe_draft();
    removed_leaf_waiver.version = 2;
    let (_, removed_record) = with_leaf_waiver
        .amend(removed_leaf_waiver)
        .expect("removed waiver amendment");
    assert_eq!(
        removed_record.invalidated,
        vec!["probe-claim", "probe-leaf", "probe-sibling"]
    );
}

#[test]
fn global_campaign_policy_changes_invalidate_all_predecessor_evidence() {
    let frozen = shared_probe_draft().freeze().expect("freeze");
    let expected = vec!["probe-claim", "probe-leaf", "probe-sibling"];

    let mut changed_explicits = shared_probe_draft();
    changed_explicits.version = 2;
    changed_explicits.explicits.capabilities = "changed probe capabilities";
    let (_, explicits_record) = frozen
        .amend(changed_explicits)
        .expect("Five Explicits amendment");
    assert_eq!(explicits_record.invalidated, expected);

    let mut changed_rules = shared_probe_draft();
    changed_rules.version = 2;
    changed_rules.amendment_rules = "changed successor-only amendment policy";
    let (_, rules_record) = frozen
        .amend(changed_rules)
        .expect("amendment-rules amendment");
    assert_eq!(rules_record.invalidated, expected);

    let mut changed_title = shared_probe_draft();
    changed_title.version = 2;
    changed_title.title = "changed campaign authority title";
    let (_, title_record) = frozen.amend(changed_title).expect("title amendment");
    assert_eq!(title_record.invalidated, expected);
}

#[test]
fn i01_amendment_names_exact_reverse_dependency_descendants() {
    let frozen = i01_draft().freeze().expect("freeze");

    // A claim change reaches its producer leaf but not the other claims
    // covered by that leaf.
    let mut successor = i01_draft();
    successor.version = 2;
    for claim in &mut successor.claims {
        if claim.id == "i01-solver-hint-optimality" {
            claim.tolerance = ToleranceSemantics::Interval { lo: 1.0, hi: 1.5 };
        }
    }
    let (amended, record) = frozen.amend(successor).expect("valid amendment");
    assert_eq!(amended.version(), 2);
    assert_eq!(record.from_version, 1);
    assert_eq!(record.to_version, 2);
    assert_eq!(record.from_digest, frozen.digest());
    assert_eq!(record.to_digest, amended.digest());
    assert_ne!(record.from_digest, record.to_digest);
    assert_eq!(
        record.invalidated,
        vec!["i01-maximal-lanes", "i01-solver-hint-optimality"]
    );

    // Removing an obligation invalidates the predecessor leaf and the
    // claim whose evidence it produced.
    let mut dropped = i01_draft();
    dropped.version = 2;
    dropped
        .obligations
        .retain(|row| row.leaf != "i01-replay-cancellation");
    // The claim it covered must be re-covered or waived for the successor
    // to freeze at all — waive it, honestly.
    dropped.waivers.push(Waiver {
        subject: "i01-deterministic-cancellation-replay",
        reason: "replay obligation withdrawn pending executor rework",
        owner: "I01 implementation beads",
        predicate: "obligation row restored via amendment",
        expiry: "next close burst",
        promotion_effect: "replay claim cannot promote while waived",
    });
    let (_, record) = frozen.amend(dropped).expect("valid amendment");
    assert_eq!(
        record.invalidated,
        vec![
            "i01-deterministic-cancellation-replay",
            "i01-replay-cancellation"
        ]
    );
}

#[test]
fn i02_draft_freezes_with_the_declared_lattice_and_coverage() {
    let frozen = fs_vmanifest::i02_draft()
        .freeze()
        .expect("the I02 seed must freeze");
    assert_eq!(frozen.initiative(), "I02");
    assert_eq!(frozen.version(), 1);
    assert_eq!(frozen.claims().len(), 9);
    let solid = frozen
        .claims()
        .iter()
        .filter(|c| c.ambition == Ambition::Solid)
        .count();
    let frontier = frozen
        .claims()
        .iter()
        .filter(|c| c.ambition == Ambition::Frontier)
        .count();
    let moonshot = frozen
        .claims()
        .iter()
        .filter(|c| c.ambition == Ambition::Moonshot)
        .count();
    assert_eq!((solid, frontier, moonshot), (5, 2, 2));
    assert_eq!(frozen.fixtures().len(), 7);
    let held_out = frozen
        .fixtures()
        .iter()
        .filter(|x| x.partition == Partition::HeldOut)
        .count();
    assert_eq!(held_out, 2, "DAE and adversarial held-out partitions");
    assert_eq!(frozen.obligations().len(), 6);
    assert_eq!(frozen.waivers().len(), 1);
    // The tearing-optimality falsifier lane is refutation-polarity.
    let tearing = frozen
        .claim("i02-globally-optimal-tearing")
        .expect("moonshot tearing claim");
    assert_eq!(tearing.polarity, ClaimPolarity::Refutation);
    // Distinct initiatives have distinct identities; both are stable.
    let i01 = i01_draft().freeze().expect("freeze");
    assert_ne!(frozen.digest(), i01.digest());
    let again = fs_vmanifest::i02_draft().freeze().expect("refreeze");
    assert_eq!(frozen.digest(), again.digest());
    // Input-order invariance holds for I02 as well (G5).
    let mut reordered = fs_vmanifest::i02_draft();
    reordered.claims.reverse();
    reordered.fixtures.reverse();
    reordered.obligations.reverse();
    assert_eq!(
        reordered.freeze().expect("freeze").digest(),
        frozen.digest()
    );
}

#[test]
fn i04_draft_freezes_with_the_declared_lattice_and_coverage() {
    let frozen = fs_vmanifest::i04_draft()
        .freeze()
        .expect("the I04 seed must freeze");
    assert_eq!(frozen.initiative(), "I04");
    assert_eq!(frozen.claims().len(), 8);
    let solid = frozen
        .claims()
        .iter()
        .filter(|c| c.ambition == Ambition::Solid)
        .count();
    let frontier = frozen
        .claims()
        .iter()
        .filter(|c| c.ambition == Ambition::Frontier)
        .count();
    let moonshot = frozen
        .claims()
        .iter()
        .filter(|c| c.ambition == Ambition::Moonshot)
        .count();
    assert_eq!((solid, frontier, moonshot), (4, 2, 2));
    assert_eq!(frozen.fixtures().len(), 8);
    let held_out = frozen
        .fixtures()
        .iter()
        .filter(|x| x.partition == Partition::HeldOut)
        .count();
    assert_eq!(held_out, 2, "defect and chart-pair held-out partitions");
    assert_eq!(frozen.obligations().len(), 5);
    assert_eq!(frozen.waivers().len(), 1);
    let completeness = frozen
        .claim("i04-counterfactual-completeness")
        .expect("moonshot completeness claim");
    assert_eq!(completeness.polarity, ClaimPolarity::Refutation);
    // Distinct from both prior initiatives; stable; order-invariant.
    let i01 = i01_draft().freeze().expect("freeze");
    let i02 = fs_vmanifest::i02_draft().freeze().expect("freeze");
    assert_ne!(frozen.digest(), i01.digest());
    assert_ne!(frozen.digest(), i02.digest());
    let mut reordered = fs_vmanifest::i04_draft();
    reordered.claims.reverse();
    reordered.fixtures.reverse();
    reordered.obligations.reverse();
    assert_eq!(
        reordered.freeze().expect("freeze").digest(),
        frozen.digest()
    );
}

#[test]
fn i08_draft_freezes_with_the_declared_lattice_and_coverage() {
    let frozen = fs_vmanifest::i08_draft()
        .freeze()
        .expect("the I08 seed must freeze");
    assert_eq!(frozen.initiative(), "I08");
    assert_eq!(frozen.claims().len(), 8);
    let solid = frozen
        .claims()
        .iter()
        .filter(|c| c.ambition == Ambition::Solid)
        .count();
    let frontier = frozen
        .claims()
        .iter()
        .filter(|c| c.ambition == Ambition::Frontier)
        .count();
    let moonshot = frozen
        .claims()
        .iter()
        .filter(|c| c.ambition == Ambition::Moonshot)
        .count();
    assert_eq!((solid, frontier, moonshot), (5, 2, 1));
    assert_eq!(frozen.fixtures().len(), 9);
    let held_out = frozen
        .fixtures()
        .iter()
        .filter(|x| x.partition == Partition::HeldOut)
        .count();
    assert_eq!(held_out, 2, "portfolio and Goodhart held-out partitions");
    assert_eq!(frozen.obligations().len(), 6);
    assert_eq!(frozen.waivers().len(), 1);
    let optimality = frozen
        .claim("i08-robust-multihorizon-optimality")
        .expect("moonshot optimality claim");
    assert_eq!(optimality.polarity, ClaimPolarity::Refutation);
    // Distinct from all prior initiatives; stable; order-invariant.
    for other in [
        i01_draft().freeze().expect("freeze").digest(),
        fs_vmanifest::i02_draft().freeze().expect("freeze").digest(),
        fs_vmanifest::i04_draft().freeze().expect("freeze").digest(),
    ] {
        assert_ne!(frozen.digest(), other);
    }
    let mut reordered = fs_vmanifest::i08_draft();
    reordered.claims.reverse();
    reordered.fixtures.reverse();
    reordered.obligations.reverse();
    assert_eq!(
        reordered.freeze().expect("freeze").digest(),
        frozen.digest()
    );
}

#[test]
fn i15_draft_freezes_with_the_declared_lattice_and_coverage() {
    let frozen = fs_vmanifest::i15_draft()
        .freeze()
        .expect("the I15 seed must freeze");
    assert_eq!(frozen.initiative(), "I15");
    assert_eq!(frozen.claims().len(), 9);
    let solid = frozen
        .claims()
        .iter()
        .filter(|c| c.ambition == Ambition::Solid)
        .count();
    let frontier = frozen
        .claims()
        .iter()
        .filter(|c| c.ambition == Ambition::Frontier)
        .count();
    let moonshot = frozen
        .claims()
        .iter()
        .filter(|c| c.ambition == Ambition::Moonshot)
        .count();
    assert_eq!((solid, frontier, moonshot), (5, 2, 2));
    assert_eq!(frozen.fixtures().len(), 8);
    let held_out = frozen
        .fixtures()
        .iter()
        .filter(|x| x.partition == Partition::HeldOut)
        .count();
    assert_eq!(
        held_out, 2,
        "edition and judgment-chain held-out partitions"
    );
    assert_eq!(frozen.obligations().len(), 6);
    assert_eq!(frozen.waivers().len(), 1);
    let impact = frozen
        .claim("i15-transitive-impact-completeness")
        .expect("moonshot impact claim");
    assert_eq!(impact.polarity, ClaimPolarity::Refutation);
    // Distinct from all prior initiatives; stable; order-invariant.
    for other in [
        i01_draft().freeze().expect("freeze").digest(),
        fs_vmanifest::i02_draft().freeze().expect("freeze").digest(),
        fs_vmanifest::i04_draft().freeze().expect("freeze").digest(),
        fs_vmanifest::i08_draft().freeze().expect("freeze").digest(),
    ] {
        assert_ne!(frozen.digest(), other);
    }
    let mut reordered = fs_vmanifest::i15_draft();
    reordered.claims.reverse();
    reordered.fixtures.reverse();
    reordered.obligations.reverse();
    assert_eq!(
        reordered.freeze().expect("freeze").digest(),
        frozen.digest()
    );
}

#[test]
fn i12_draft_freezes_with_the_declared_lattice_and_coverage() {
    let frozen = fs_vmanifest::i12_draft()
        .freeze()
        .expect("the I12 seed must freeze");
    assert_eq!(frozen.initiative(), "I12");
    assert_eq!(frozen.claims().len(), 10);
    let solid = frozen
        .claims()
        .iter()
        .filter(|c| c.ambition == Ambition::Solid)
        .count();
    let frontier = frozen
        .claims()
        .iter()
        .filter(|c| c.ambition == Ambition::Frontier)
        .count();
    let moonshot = frozen
        .claims()
        .iter()
        .filter(|c| c.ambition == Ambition::Moonshot)
        .count();
    assert_eq!((solid, frontier, moonshot), (6, 2, 2));
    assert_eq!(frozen.fixtures().len(), 7);
    let held_out = frozen
        .fixtures()
        .iter()
        .filter(|x| x.partition == Partition::HeldOut)
        .count();
    assert_eq!(
        held_out, 3,
        "rank-changing DAE, topology-reset, and tangency batteries stay held out"
    );
    assert_eq!(frozen.obligations().len(), 6);
    assert_eq!(frozen.waivers().len(), 1);
    let tripwire = frozen
        .claim("i12-grazing-false-certificate-falsifier")
        .expect("the Sev-0 grazing tripwire");
    assert_eq!(tripwire.polarity, ClaimPolarity::Refutation);
    assert_eq!(tripwire.ambition, Ambition::Solid);
    // Distinct from prior initiatives; stable; order-invariant.
    for other in [
        i01_draft().freeze().expect("freeze").digest(),
        fs_vmanifest::i02_draft().freeze().expect("freeze").digest(),
        fs_vmanifest::i15_draft().freeze().expect("freeze").digest(),
    ] {
        assert_ne!(frozen.digest(), other);
    }
    let mut reordered = fs_vmanifest::i12_draft();
    reordered.claims.reverse();
    reordered.fixtures.reverse();
    reordered.obligations.reverse();
    assert_eq!(
        reordered.freeze().expect("freeze").digest(),
        frozen.digest()
    );
}
