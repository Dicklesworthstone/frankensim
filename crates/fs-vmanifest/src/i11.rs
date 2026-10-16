//! The I11 (experimental-campaign and DAQ compiler) VerificationManifest
//! draft (bead frankensim-leapfrog-2026-program-i94v.3.4.7.1).
//!
//! Baseline lattice ([S]): typed ContextOfUse-to-campaign IR, typed
//! signal chains (sensors/clocks/transfer functions/covariance),
//! deterministic safe excitation with interlock/abort/recovery
//! semantics, bounded DAQ packs with flagged-never-interpolated
//! dropout, calibration-chain integrity, and a standing
//! false-provenance refutation lane against stream tampering. Maximal
//! lattice ([F]/[M]): authenticated resumable raw streams, blind
//! partition discipline with receipted one-way unblinding, then
//! proof-carrying experiment execution with end-to-end uncertainty
//! attribution and adaptive safe excitation. A weaker receipt closes
//! its own element and is never relabeled as the stronger theorem;
//! version-1 prose cards mint no proof authority.

use crate::{
    Ambition, CampaignTier, ClaimPolarity, ClaimSpec, FiveExplicits, FixturePin, FixtureSource,
    GauntletTier, ManifestDraft, ObligationRow, OracleRoute, Partition, ToleranceSemantics, Waiver,
};

/// Build the I11 draft. Consumers freeze it themselves; the conformance
/// battery proves it freezes.
#[must_use]
pub fn i11_draft() -> ManifestDraft {
    ManifestDraft {
        initiative: "I11",
        title: "Experimental-campaign and DAQ compiler gate: ContextOfUse to typed \
                signal chains, safe excitation, bounded authenticated acquisition, \
                calibration-chained provenance, and blind-partition adjudication",
        version: 1,
        explicits: FiveExplicits {
            units: "SI base units throughout; sample times in seconds, rates in Hz, \
                    channel quantities in the per-fixture declared unit with typed \
                    transfer-function dimensions; exact bitwise/boolean verdicts use \
                    unit 'bit'; covariance entries in squared channel units",
            seeds: "Philox 4x32-10 counter streams keyed 'i11/<fixture-id>/<case-index>'; \
                    development indices 0..=16383; core held-out indices 65536..=81919; \
                    maximal held-out indices 131072..=147455 (disjoint by construction, \
                    split frozen here; per-holdout subranges pinned in the fixture \
                    specs)",
            budgets: "smoke tier <= 60 s on one host; core tier <= 30 min; max tier <= 8 h \
                      on a quiet perf host; <= 16 GiB memory per lane; DAQ pack rate/ \
                      size/duration caps are per-fixture declarations whose violation \
                      is a typed refusal; accuracy budgets are the per-claim tolerance \
                      fields",
            versions: "fs-vmanifest schema v2; toolchain pinned by \
                       rust-toolchain.toml; sibling pins by constellation.lock",
            capabilities: "no network; no FFI; deterministic mode mandatory for every G5 \
                           row; frontier/moonshot lanes stay behind feature flags; \
                           physical-rig and durable-receipt dependencies are waivered \
                           until their packs/beads land",
        },
        claims: i11_claims(),
        fixtures: i11_fixtures(),
        obligations: i11_obligations(),
        waivers: i11_waivers(),
        amendment_rules: "Any change is a successor version through FrozenManifest::amend; \
                          the amendment record names every invalidated claim/obligation \
                          descendant; an amendment after campaign start invalidates the \
                          affected evidence, which must be re-earned; there is no in-place \
                          edit path in the type system",
    }
}

#[allow(clippy::too_many_lines)]
fn i11_claims() -> Vec<ClaimSpec> {
    vec![
        ClaimSpec {
            id: "i11-typed-campaign-ir",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "The ContextOfUse-to-campaign IR (objectives, channels, \
                        excitation plans, acquisition windows, adjudication rules) \
                        round-trips its canonical encoding bit-stably and refuses \
                        every enumerated malformed class with a typed refusal naming \
                        the violated rule",
            hypotheses: &[
                "campaigns drawn from the pinned fixture families",
                "malformed classes are the enumerated ones: dangling channel \
                 references, unitless quantities, overlapping exclusive windows, \
                 adjudication rules naming absent partitions",
            ],
            qoi: "roundtrip_and_refusal_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G0,
            oracle: OracleRoute {
                identity: "independent decoder plus the fixture generator's own \
                           validity certificate (emitted during generation, not by \
                           the compiler)",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "first implementation bead of the I11 campaign-IR leaf opens",
            kill: "any accepted malformed campaign or refused valid campaign returns \
                   the IR to design review with the disagreeing fixture as receipt",
            fallback: "explicit user-supplied campaign normal form with per-field \
                       validation",
            no_claim: "no claim that an admitted campaign is scientifically adequate; \
                       admission is typing, not experiment design review",
        },
        ClaimSpec {
            id: "i11-typed-signal-chain",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Typed sensor/clock/transfer-function/covariance chains are \
                        dimension-checked end to end at admission, covariance \
                        declarations are validated symmetric positive-semidefinite, \
                        and every chain element carries its declared bandwidth and \
                        latency bounds",
            hypotheses: &[
                "chains drawn from the pinned families incl. multi-rate and \
                 cross-clock variants",
                "PSD validation at the declared numerical tolerance of the fixture",
            ],
            qoi: "chain_admission_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G0,
            oracle: OracleRoute {
                identity: "independent dimensional-analysis walker plus a separate \
                           eigenvalue-based PSD check on the fixture's own covariance \
                           certificate",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "campaign-IR leaf green at smoke tier",
            kill: "any dimension mismatch or indefinite covariance reaching execution \
                   blocks the admission leaf",
            fallback: "chains refuse with the first violated dimension/PSD diagnostic",
            no_claim: "no claim that declared covariances match physical reality; \
                       declaration validity only",
        },
        ClaimSpec {
            id: "i11-safe-excitation-interlocks",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Generated excitation never leaves the declared safety \
                        envelope (amplitude, rate, energy), interlock trips abort \
                        through request->drain->finalize into the declared safe \
                        state, and every abort/recovery carries a typed receipt — \
                        normal completion, safe-state completion, and containment \
                        failure are distinct recorded outcomes",
            hypotheses: &[
                "envelopes and interlock conditions from the pinned rig families",
                "synthetic interlock models (the physical-rig pack is the waivered \
                 dependency)",
            ],
            qoi: "envelope_and_interlock_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G2,
            oracle: OracleRoute {
                identity: "independent envelope monitor replaying the excitation \
                           stream against the declared limits (separate code path \
                           from the generator)",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "signal-chain leaf green at smoke tier",
            kill: "any envelope escape or undrained interlock abort is Sev-0 for the \
                   excitation lane",
            fallback: "excitation refuses to start without a validated envelope and \
                       interlock declaration",
            no_claim: "no physical-rig safety authority (waivered); synthetic-model \
                       semantics only; not IEC/ISO regulatory certification",
        },
        ClaimSpec {
            id: "i11-bounded-daq-packs",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "DAQ packs respect their declared rate/size/duration caps with \
                        typed refusals at cap violations, and saturation, dropout, \
                        and clock-skew artifacts are flagged in-band — a dropout is \
                        flagged, never interpolated silently",
            hypotheses: &[
                "packs from the pinned families incl. the saturation/dropout \
                 adversarial holdout",
                "caps are per-fixture declarations, not runtime discoveries",
            ],
            qoi: "pack_bounds_and_flagging_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G0,
            oracle: OracleRoute {
                identity: "independent pack auditor recomputing caps and artifact \
                           flags from the raw bytes (not the acquisition path)",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "signal-chain leaf green at smoke tier",
            kill: "any silent interpolation or unflagged saturation is Sev-0 for the \
                   acquisition lane",
            fallback: "packs exceeding caps refuse before acquisition starts",
            no_claim: "no claim about sensor physical validity; artifact flagging \
                       fidelity only",
        },
        ClaimSpec {
            id: "i11-calibration-chain-integrity",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Every channel's transfer function carries an unbroken, dated \
                        calibration lineage to a declared reference; a broken, \
                        expired, or circular link demotes the channel's evidence \
                        color with a typed diagnostic and never silently passes",
            hypotheses: &[
                "calibration graphs from the pinned families incl. the fault \
                 holdout's broken/expired/circular cases",
                "reference standards are fixture declarations",
            ],
            qoi: "lineage_integrity_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G0,
            oracle: OracleRoute {
                identity: "independent lineage walker over the fixture's own \
                           calibration-graph certificate",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "signal-chain leaf green at smoke tier",
            kill: "any broken lineage passing undemoted blocks the calibration lane",
            fallback: "channels without valid lineage carry Estimated color with the \
                       lineage gap named",
            no_claim: "no metrological accuracy claim about the reference standards \
                       themselves; lineage structure only",
        },
        ClaimSpec {
            id: "i11-false-provenance-falsifier",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Refutation,
            statement: "Adversarial stream families (splices, reorders, truncations, \
                        bit flips, replayed segments, forged sequence numbers) \
                        attempt to pass raw-stream authentication; any acceptance \
                        refutes the provenance lane's discipline",
            hypotheses: &[
                "tampering constructed by the fixture generator with ground-truth \
                 tamper maps in-spec",
                "authentication runs at production settings",
            ],
            qoi: "false_acceptance_count",
            unit: "1",
            tolerance: ToleranceSemantics::Interval { lo: 0.0, hi: 0.0 },
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "the fixture generator's ground-truth tamper map (not the \
                           authenticator)",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "runs in every campaign tier from the first authentication \
                         commit",
            kill: "this lane is never killed; it is the standing tripwire",
            fallback: "none: a nonzero count is a release blocker",
            no_claim: "absence of refutation is not completeness and cannot prove \
                       authenticator soundness; the lane only falsifies",
        },
        ClaimSpec {
            id: "i11-authenticated-raw-streams",
            ambition: Ambition::Frontier,
            polarity: ClaimPolarity::Affirmative,
            statement: "Raw acquisition streams are content-addressed and \
                        sequence-authenticated with resumable manifests: truncation, \
                        reordering, and replacement are detected on resume, and a \
                        verified manifest replays byte-identically",
            hypotheses: &[
                "streams from the pinned families incl. multi-device interleaving",
                "authentication metadata travels with the stream, never in a side \
                 channel",
            ],
            qoi: "stream_authentication_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G2,
            oracle: OracleRoute {
                identity: "independent verifier re-deriving content addresses from \
                           raw bytes (separate from the writer path)",
                independent: true,
                tcb_overlap: "shares the hash primitive (fs-blake3) only",
            },
            activation: "DAQ-pack and calibration leaves green at core tier",
            kill: "any undetected tamper class from the falsifier battery kills the \
                   lane",
            fallback: "streams carry Estimated provenance with the missing \
                       authentication step named",
            no_claim: "authentication proves byte identity and ordering, not sensor \
                       truth; a perfectly authenticated stream can still carry a \
                       miscalibrated signal",
        },
        ClaimSpec {
            id: "i11-blind-partition-discipline",
            ambition: Ambition::Frontier,
            polarity: ClaimPolarity::Affirmative,
            statement: "Blinded validation partitions are capability-separated: \
                        analysis paths cannot read blind labels before adjudication, \
                        every access attempt is receipted, and unblinding is a \
                        receipted one-way event that invalidates further blind-tier \
                        claims for the campaign",
            hypotheses: &[
                "blind assignments from the pinned adjudication holdout",
                "capability separation is the typed API boundary, not a convention",
            ],
            qoi: "blinding_discipline_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G2,
            oracle: OracleRoute {
                identity: "access-receipt audit trail replayed against the declared \
                           adjudication schedule (independent of the analysis code)",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "stream-authentication leaf green at core tier",
            kill: "any unreceipted blind-label access is Sev-0 for the adjudication \
                   lane",
            fallback: "campaigns without blind partitions run unblinded and say so in \
                       every report",
            no_claim: "blinding prevents information flow through the typed API only; \
                       out-of-band leakage is a process problem the manifest cannot \
                       certify",
        },
        ClaimSpec {
            id: "i11-proof-carrying-execution",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Affirmative,
            statement: "Experiment execution is proof-carrying: every reported \
                        uncertainty decomposes into attributed sensor, calibration, \
                        and processing terms whose composition is machine-checkable, \
                        and any unattributed residual is a named Unknown term, never \
                        absorbed silently",
            hypotheses: &[
                "campaigns from the pinned families with declared uncertainty \
                 budgets",
                "durable receipt schemas from the ledger migration are available \
                 (waivered dependency)",
            ],
            qoi: "attribution_completeness_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "independent recomposition checker summing attributed \
                           terms against the reported totals",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "behind the i11-moonshot feature flag on a pre-proof \
                         successor version; frontier lanes green",
            kill: "any silently absorbed residual on a declared-budget fixture kills \
                   the lane",
            fallback: "reports carry unattributed totals labeled Estimated",
            no_claim: "attribution correctness is relative to the declared budget \
                       model; version-1 prose mints no attribution authority",
        },
        ClaimSpec {
            id: "i11-adaptive-safe-excitation",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Affirmative,
            statement: "Adaptive excitation planning stays inside the declared safety \
                        envelope with certified margins at every adaptation step, and \
                        proposed evidence-gap closures are typed suggestions that \
                        never self-execute",
            hypotheses: &[
                "adaptation policies from the pinned families with declared margin \
                 requirements",
                "the safe-excitation [S] lane is green at core tier",
            ],
            qoi: "adaptive_margin_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "independent envelope monitor with margin recomputation \
                           (same oracle family as the [S] lane, distinct instance)",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "behind the i11-moonshot feature flag on a pre-proof \
                         successor version; requires the physical-rig waiver status \
                         reviewed",
            kill: "any adaptation step with an uncertified margin is Sev-0 for this \
                   lane",
            fallback: "static preregistered excitation plans only",
            no_claim: "no autonomous experiment authority: gap-closure proposals \
                       require human or governed approval; version-1 prose mints no \
                       theorem authority",
        },
    ]
}

const POLICY_SPEC: &str = "I11_CAMPAIGN_POLICY_V1
CAMPAIGN_IR= typed ContextOfUse-to-campaign lowering; admission is typing, never \
experiment-design adequacy; adjudication rules must name existing partitions
SIGNAL_CHAIN= sensors/clocks/transfer functions/covariance are dimension-checked end \
to end; covariance declarations validated PSD; declared bandwidth/latency travel with \
every element
EXCITATION_SAFETY= excitation never leaves the declared envelope; interlock trips \
abort through request->drain->finalize into the declared safe state; normal \
completion, safe-state completion, containment failure are distinct recorded outcomes
DAQ_PACKS= rate/size/duration caps are declarations with typed refusals; a dropout is \
flagged, never interpolated silently; saturation and clock skew are in-band artifacts
CALIBRATION= every channel carries an unbroken dated lineage to a declared reference; \
broken/expired/circular links demote evidence color, never silently pass
STREAM_AUTHENTICITY= raw streams are content-addressed and sequence-authenticated \
with resumable manifests; authentication proves byte identity and ordering, not \
sensor truth
BLINDING= blind labels are capability-separated from analysis paths; every access is \
receipted; unblinding is a receipted one-way event
THEOREM_AUTHORITY= version 1 has prose cards only and mints no proof; moonshot lanes \
activate only on a pre-proof successor version with machine-readable obligations
EVIDENCE_STATES= Verified/Validated/Estimated/Failed/Refuted/Unknown/no-claim are \
distinct; one axis never substitutes for another; a weaker receipt closes its own \
lattice element only
HOLDOUT= held-out fixtures adjudicate; development never touches held-out indices; \
each holdout has exactly one stage-local consumer row
LIFECYCLE= request->drain->finalize with checkpoint boundaries; cancellation is \
drained, never dropped; partial success cannot publish normal authority
LOGGING= structured fs-obs events per obligation row incl. the six lifecycle kinds; \
budgets, seeds, versions, capabilities logged at run start
RETENTION= receipts, failure bundles, and adjudication records are content-addressed \
and retained for replay; raw counters are diagnostics only
FAILURE_BUNDLE= every red or Unknown outcome retains a bounded reproducible bundle \
naming fixture, seed, budget, and disposition
PROMOTION= claims promote only through their preregistered obligation rows on frozen \
manifests; amendment invalidates affected descendants which must re-earn evidence
LEAF_REQUIREMENT= every execution leaf maps to exactly one obligation row; there are \
no unnamed skips";

fn i11_fixtures() -> Vec<FixturePin> {
    vec![
        FixturePin {
            id: "i11-campaign-policy-v1",
            source: FixtureSource::AuthoredSpec { spec: POLICY_SPEC },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i11-benchmark-campaign-family",
            source: FixtureSource::AuthoredSpec {
                spec: "i11 fixture v1: pinned benchmark campaigns — frequency-response \
                       identification of a two-mass oscillator and a thermal step \
                       test; channels {force N, velocity m/s, temperature K}; \
                       campaign IR incl. objectives, windows, adjudication rules; \
                       generator emits validity certificates; development indices \
                       0..=16383; seeds key 'i11/benchmark/<k>'",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i11-motor-rig-family",
            source: FixtureSource::AuthoredSpec {
                spec: "i11 fixture v1: synthetic motor-rig models — torque/speed/ \
                       current channels, declared excitation envelopes (amplitude, \
                       slew, energy), interlock conditions (overspeed, overcurrent, \
                       thermal), safe-state definitions, and abort/recovery \
                       schedules; development indices 0..=16383; seeds key \
                       'i11/motor/<k>'",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i11-clock-drift-family",
            source: FixtureSource::AuthoredSpec {
                spec: "i11 fixture v1: multi-device clock models with declared drift, \
                       jitter, and offset envelopes; cross-clock alignment cases and \
                       resync events; multi-rate interleaved streams; development \
                       indices 0..=16383; seeds key 'i11/clock/<k>'",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i11-saturation-dropout-core-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "i11 fixture v1 HOLDOUT: adversarial acquisition artifacts — \
                       saturation plateaus, dropout gaps, duplicated samples, and \
                       cap-boundary packs with ground-truth artifact maps in-spec; \
                       core held-out indices 65536..=69631; one I11.G3 consumer \
                       (i11-daq-packs); seeds key 'i11/satdrop/<k>'",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "i11-calibration-fault-core-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "i11 fixture v1 HOLDOUT: calibration-graph faults — broken \
                       links, expired certificates, circular lineages, and reference \
                       mismatches with ground-truth fault maps in-spec; core held-out \
                       indices 69632..=73727; one I11.G3 consumer \
                       (i11-calibration-chains); seeds key 'i11/calfault/<k>'",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "i11-blinded-adjudication-max-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "i11 fixture v1 HOLDOUT: blinded validation campaigns with \
                       sealed label assignments, declared adjudication schedules, \
                       and attribution-budget declarations for the proof-carrying \
                       lane; maximal held-out indices 131072..=135167; one I11.G3 \
                       consumer (i11-moonshot-proof-adaptive); seeds key \
                       'i11/blind/<k>'",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "i11-tamper-stream-family",
            source: FixtureSource::AuthoredSpec {
                spec: "i11 fixture v1: adversarial stream battery — splices, \
                       reorders, truncations, bit flips, replayed segments, forged \
                       sequence numbers, each with a ground-truth tamper map; used \
                       by the falsifier and authentication lanes; development \
                       indices 0..=16383; seeds key 'i11/tamper/<k>'",
            },
            partition: Partition::Development,
        },
    ]
}

#[allow(clippy::too_many_lines)]
fn i11_obligations() -> Vec<ObligationRow> {
    const UNIT_CASES: &[&str] = &[
        "happy",
        "empty",
        "boundary",
        "max",
        "error",
        "unit-dimension",
        "tie-break",
        "cancellation",
        "migration",
    ];
    vec![
        ObligationRow {
            leaf: "i11-campaign-admission",
            claims_covered: &["i11-typed-campaign-ir", "i11-typed-signal-chain"],
            unit_cases: UNIT_CASES,
            g0: "generators: campaign and signal-chain families incl. enumerated \
                 malformed classes; validity predicate: generator certificate \
                 agreement; laws: canonical round-trip bit-stability, refusal \
                 totality, dimension-check completeness, PSD validation; shrinker: \
                 channel removal preserving the violated rule; replay seeds per \
                 explicits",
            decks: &[
                "i11-benchmark-campaign-family",
                "i11-campaign-policy-v1",
                "i11-clock-drift-family",
            ],
            g3_relations: &[
                "channel relabeling invariance",
                "unit-rescaling covariance of admitted chains",
            ],
            g4_schedule: "request->drain->finalize injection between decode, \
                          dimension walk, and covariance validation; checkpoint at \
                          each phase boundary; drained cancellation leaves no \
                          partially admitted campaign",
            g5_matrix: "threads {1,2,7} x deterministic mode x same-ISA replay",
            entry_point: "scripts/e2e/leapfrog/i11_admission.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i11-admission slice)",
            obs_events: &[
                "request.received",
                "cancel.requested",
                "drain.completed",
                "finalize.completed",
                "failure_bundle.retained",
                "adjudication.receipt",
                "campaign.admitted",
                "campaign.refused",
            ],
            replay_command: "scripts/e2e/leapfrog/i11_admission.sh --manifest \
                             <manifest-id> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i11-excitation-safety",
            claims_covered: &["i11-safe-excitation-interlocks"],
            unit_cases: UNIT_CASES,
            g0: "generators: motor-rig envelope/interlock families; laws: envelope \
                 containment of every generated sample, interlock trip totality, \
                 distinct outcome recording (normal/safe-state/containment-failure); \
                 replay seeds per explicits",
            decks: &["i11-campaign-policy-v1", "i11-motor-rig-family"],
            g3_relations: &["envelope-tightening monotonicity of refusals"],
            g4_schedule: "the core of this lane IS G4: request->drain->finalize \
                          interlock aborts at every excitation phase; checkpoint at \
                          plan-step boundaries; injected trips must drain into the \
                          declared safe state with receipts",
            g5_matrix: "threads {1,2,7} x deterministic mode x same-ISA replay of \
                        excitation streams and abort receipts",
            entry_point: "scripts/e2e/leapfrog/i11_excitation.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i11-excitation slice)",
            obs_events: &[
                "request.received",
                "cancel.requested",
                "drain.completed",
                "finalize.completed",
                "failure_bundle.retained",
                "adjudication.receipt",
                "excitation.emitted",
                "interlock.tripped",
                "safestate.reached",
            ],
            replay_command: "scripts/e2e/leapfrog/i11_excitation.sh --manifest \
                             <manifest-id> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i11-daq-packs",
            claims_covered: &["i11-bounded-daq-packs"],
            unit_cases: UNIT_CASES,
            g0: "generators: pack families at and beyond declared caps plus the \
                 adversarial artifact holdout; laws: cap refusal exactness (boundary \
                 and one-beyond), artifact flag completeness against ground-truth \
                 maps, no silent interpolation; replay seeds per explicits",
            decks: &[
                "i11-campaign-policy-v1",
                "i11-clock-drift-family",
                "i11-saturation-dropout-core-holdout",
            ],
            g3_relations: &[
                "pack-splitting invariance of artifact flags",
                "rate-rescaling covariance of cap refusals",
            ],
            g4_schedule: "request->drain->finalize injection mid-acquisition; \
                          checkpoint at pack boundaries; drained cancellation retains \
                          complete packs only, never a torn pack",
            g5_matrix: "threads {1,2,7} x deterministic mode x same-ISA replay of \
                        packs and artifact flags",
            entry_point: "scripts/e2e/leapfrog/i11_daq.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i11-daq slice)",
            obs_events: &[
                "request.received",
                "cancel.requested",
                "drain.completed",
                "finalize.completed",
                "failure_bundle.retained",
                "adjudication.receipt",
                "artifact.flagged",
                "pack.refused",
                "pack.sealed",
            ],
            replay_command: "scripts/e2e/leapfrog/i11_daq.sh --manifest \
                             <manifest-id> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i11-calibration-chains",
            claims_covered: &["i11-calibration-chain-integrity"],
            unit_cases: UNIT_CASES,
            g0: "generators: calibration graphs incl. the fault holdout's broken/ \
                 expired/circular cases; laws: lineage walk totality, demotion \
                 exactness against ground-truth fault maps, dated-link ordering; \
                 replay seeds per explicits",
            decks: &[
                "i11-benchmark-campaign-family",
                "i11-calibration-fault-core-holdout",
                "i11-campaign-policy-v1",
            ],
            g3_relations: &["graph relabeling invariance of lineage verdicts"],
            g4_schedule: "request->drain->finalize injection mid-walk; checkpoint at \
                          per-channel boundaries; drained cancellation reports \
                          walked channels only",
            g5_matrix: "threads {1,2,7} x deterministic mode x same-ISA replay of \
                        lineage verdicts and demotions",
            entry_point: "scripts/e2e/leapfrog/i11_calibration.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i11-calibration slice)",
            obs_events: &[
                "request.received",
                "cancel.requested",
                "drain.completed",
                "finalize.completed",
                "failure_bundle.retained",
                "adjudication.receipt",
                "channel.demoted",
                "lineage.verified",
            ],
            replay_command: "scripts/e2e/leapfrog/i11_calibration.sh --manifest \
                             <manifest-id> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i11-stream-authentication",
            claims_covered: &[
                "i11-authenticated-raw-streams",
                "i11-blind-partition-discipline",
                "i11-false-provenance-falsifier",
            ],
            unit_cases: UNIT_CASES,
            g0: "generators: authenticated stream families plus the full tamper \
                 battery; laws: byte-identical verified replay, tamper detection \
                 completeness against ground-truth maps, zero false acceptances, \
                 receipted blind-label access, one-way unblinding; replay seeds per \
                 explicits",
            decks: &[
                "i11-campaign-policy-v1",
                "i11-clock-drift-family",
                "i11-tamper-stream-family",
            ],
            g3_relations: &[
                "chunking invariance of stream authentication",
                "device-interleaving invariance of sequence verdicts",
            ],
            g4_schedule: "request->drain->finalize injection mid-stream and \
                          mid-verification; checkpoint at manifest boundaries; \
                          resumed verification must agree bitwise with \
                          uninterrupted verification",
            g5_matrix: "threads {1,2,7} x deterministic mode x same-ISA replay of \
                        authentication verdicts and access receipts",
            entry_point: "scripts/e2e/leapfrog/i11_streams.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i11-streams slice)",
            obs_events: &[
                "request.received",
                "cancel.requested",
                "drain.completed",
                "finalize.completed",
                "failure_bundle.retained",
                "adjudication.receipt",
                "blind.access",
                "stream.verified",
                "tamper.detected",
            ],
            replay_command: "scripts/e2e/leapfrog/i11_streams.sh --manifest \
                             <manifest-id> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i11-moonshot-proof-adaptive",
            claims_covered: &[
                "i11-adaptive-safe-excitation",
                "i11-proof-carrying-execution",
            ],
            unit_cases: UNIT_CASES,
            g0: "generators: declared-budget campaigns and adaptive policies over the \
                 blinded holdout; laws: attribution recomposition exactness, named \
                 Unknown residuals, certified margins at every adaptation step, \
                 proposals never self-execute; replay seeds per explicits",
            decks: &[
                "i11-blinded-adjudication-max-holdout",
                "i11-campaign-policy-v1",
                "i11-motor-rig-family",
            ],
            g3_relations: &["budget-refinement monotonicity of attribution terms"],
            g4_schedule: "the core of this lane IS G4: request->drain->finalize \
                          injection with budget exhaustion, timeout, and interlock \
                          trips inside adaptation; checkpoint at adaptation-step \
                          boundaries; each outcome drained and typed; BudgetExhausted \
                          stays Unknown",
            g5_matrix: "threads {1,2,7} x deterministic mode x same-ISA replay of \
                        attribution and adaptation receipts",
            entry_point: "scripts/e2e/leapfrog/i11_moonshot.sh",
            tier: CampaignTier::Max,
            dsr_lane: "dsr quality --tool frankensim (i11-moonshot slice)",
            obs_events: &[
                "request.received",
                "cancel.requested",
                "drain.completed",
                "finalize.completed",
                "failure_bundle.retained",
                "adjudication.receipt",
                "adaptation.margin",
                "attribution.recomposed",
                "budget.receipt",
            ],
            replay_command: "scripts/e2e/leapfrog/i11_moonshot.sh --manifest \
                             <manifest-id> --replay <artifact-id>",
        },
    ]
}

fn i11_waivers() -> Vec<Waiver> {
    vec![
        Waiver {
            subject: "i11-excitation-safety",
            reason: "no physical rig pack exists in CI; interlock and envelope \
                     semantics are exercised against synthetic rig models only, so \
                     physical actuation authority cannot be earned in this campaign \
                     version",
            owner: "initiative-11 lane owner",
            predicate: "a governed physical motor-rig pack with its own safety \
                        sign-off lands and is pinned as a successor-version fixture",
            expiry: "first I11 campaign review after a rig pack lands; re-justify or \
                     retire at every manifest amendment",
            promotion_effect: "the [S] excitation claim closes on synthetic-model \
                               semantics only; no physical-rig safety authority is \
                               minted while the waiver is live",
        },
        Waiver {
            subject: "i11-moonshot-proof-adaptive",
            reason: "proof-carrying execution needs the durable receipt schemas from \
                     the ledger/package migration (bead \
                     frankensim-ext-ledger-package-migration-h61n, still open) and \
                     the adaptive lane needs the [S] excitation lane green at core \
                     tier; neither dependency has landed central proof",
            owner: "initiative-11 lane owner",
            predicate: "bead frankensim-ext-ledger-package-migration-h61n closes \
                        with green central proof AND the I11 [S] obligations are \
                        green at core tier",
            expiry: "first I11 campaign review after both dependencies land; \
                     re-justify or retire at every manifest amendment",
            promotion_effect: "the [M] claims stay Unknown and cannot close; [S]/[F] \
                               promotion is unaffected because their obligations run \
                               on synthetic fixtures with in-crate receipts",
        },
    ]
}
