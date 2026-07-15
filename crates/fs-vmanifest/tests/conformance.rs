//! Conformance battery for the VerificationManifest schema and the I01
//! instance (bead frankensim-leapfrog-2026-program-i94v.1.1.8.1).
//!
//! Coverage maps to the bead's test plan: G0 schema/property gates
//! (required fields, identity stability, lattice separation, amendment
//! invalidation), G3 mutations (weakened hypotheses, swapped held-out
//! fixtures, loosened bands, production-oracle reuse — each fails closed
//! or moves the identity), G4 interruption/resume equivalence of manifest
//! assembly, and G5 byte-stable canonicalization on the same ISA.

use fs_vmanifest::{
    Ambition, AmendmentRefusal, CampaignTier, ClaimPolarity, ClaimSpec, FiveExplicits, FixturePin,
    FixtureSource, FreezeRefusal, GauntletTier, MAX_CLAIMS, ManifestDraft, ObligationRow,
    OracleRoute, Partition, ToleranceSemantics, Waiver, claim_digest, fixture_digest, i01_draft,
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
            versions: "schema v1",
            capabilities: "none",
        },
        claims: vec![PROBE_CLAIM],
        fixtures: vec![PROBE_FIXTURE],
        obligations: vec![PROBE_OBLIGATION],
        waivers: Vec::new(),
        amendment_rules: "successor versions only",
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
        Err(FreezeRefusal::OverCap { what: "claims", .. })
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

    let mut duplicate = probe_draft();
    duplicate.claims = vec![PROBE_CLAIM, PROBE_CLAIM];
    match duplicate.freeze() {
        Err(FreezeRefusal::DuplicateId { kind: "claim", id }) => {
            assert_eq!(id, "probe-claim");
        }
        other => panic!("duplicate claim must fail closed, got {other:?}"),
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
fn interrupted_two_stage_assembly_matches_one_shot() {
    // G4 analogue at this layer: a draft assembled in stages (as if
    // interrupted and resumed) freezes to the same identity as the
    // one-shot assembly.
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
fn per_component_identity_is_mutation_sensitive() {
    let mut statement = PROBE_CLAIM;
    statement.statement = "different statement";
    assert_ne!(claim_digest(&statement), claim_digest(&PROBE_CLAIM));

    let mut tier = PROBE_CLAIM;
    tier.evidence_tier = GauntletTier::G5;
    assert_ne!(claim_digest(&tier), claim_digest(&PROBE_CLAIM));

    let mut polarity = PROBE_CLAIM;
    polarity.polarity = ClaimPolarity::Refutation;
    assert_ne!(claim_digest(&polarity), claim_digest(&PROBE_CLAIM));

    let mut partition = PROBE_FIXTURE;
    partition.partition = Partition::HeldOut;
    assert_ne!(fixture_digest(&partition), fixture_digest(&PROBE_FIXTURE));

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
    assert_eq!(fixture_digest(&upper), fixture_digest(&lower));
}

#[test]
fn amendment_requires_successor_version_and_names_invalidated_descendants() {
    let frozen = i01_draft().freeze().expect("freeze");

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

    // A valid amendment touching exactly one claim invalidates exactly
    // that claim and nothing else.
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
    assert_eq!(record.invalidated, vec!["i01-solver-hint-optimality"]);

    // Removing an obligation invalidates its leaf.
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
    assert!(
        record
            .invalidated
            .contains(&"i01-replay-cancellation".to_string()),
        "dropped obligation leaves are named: {record:?}"
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
