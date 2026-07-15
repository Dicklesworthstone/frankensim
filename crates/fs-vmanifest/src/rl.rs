//! The RL reality-loop portfolio VerificationManifest draft (bead
//! frankensim-leapfrog-2026-program-i94v.3.5.1).
//!
//! The first CP/EM/RL/PD portfolio aggregate. RL freezes the SEAMS of the
//! lot-to-experiment-to-deployment loop — the composition of the I05
//! deployment-twin, I06 material-lot, I10 identifiability/coupon-design and
//! I11 campaign/DAQ instance manifests — not their internals: one identity
//! spine across stage boundaries; end-to-end blind-partition custody;
//! exactly-once uncertainty ownership; a frozen calibration/measurement
//! graph; deterministic whole-loop replay; selective reproof after mid-loop
//! change; and weakest-wins deployment gating. Stronger lanes target the
//! anytime-valid adaptive OED loop, target-exact HIL/timing composition, a
//! machine-checked end-to-end composition theorem, and physical-campaign
//! receipt parity. A refutation lane attacks forged cross-stage
//! certificates. Preregistration is governance, not proof of any lot,
//! model, campaign, target, or deployed answer.

use crate::{
    Ambition, CampaignTier, ClaimPolarity, ClaimSpec, FiveExplicits, FixturePin, FixtureSource,
    GauntletTier, ManifestDraft, ObligationRow, OracleRoute, Partition, ToleranceSemantics, Waiver,
};

const CAMPAIGN_POLICY_FIXTURE: &str = "rl-campaign-policy-v1";

/// Build the RL draft. Consumers freeze it themselves; focused
/// conformance tests prove that the authored seed satisfies schema v2.
#[must_use]
pub fn rl_draft() -> ManifestDraft {
    ManifestDraft {
        initiative: "RL",
        title: "Reality-loop portfolio gate, lot evidence to blind campaign to deployed twin: \
                one identity spine, sealed partitions, owned uncertainty, replayable loop, \
                selective reproof, and certified composition",
        version: 1,
        explicits: FiveExplicits {
            units: "SI coherent values throughout; stage-boundary quantities keep their \
                    instance-manifest units and a boundary re-VALIDATES units/frames/domains but never \
                    re-derives or re-converts them; timing/fault quantities in s with \
                    target-exact scope; information-gain and false-certificate counts \
                    dimensionless; exact/refusal/composition verdicts use unit 'bit'",
            seeds: "Philox 4x32-10 streams keyed by BLAKE3 domain \
                    org.frankensim.rl.fixture-stream.v1 over exact fixture, loop-world, stage \
                    and case ids; development indices 0..=16383, core held-out 65536..=77823, \
                    maximal held-out 131072..=143359; instance-manifest streams stay in their \
                    own i05/i06/i10/i11 domains and are never reused as loop seeds; duplicate \
                    logical ids or cross-stage reads are IntegrityFailed",
            budgets: "smoke <= 120 s and 2 GiB; core <= 90 min and 16 GiB; maximal <= 24 h and \
                      64 GiB; <= 10^3 loop worlds, 10^5 stage transitions per world, 10^6 \
                      retention receipts, 10^5 posterior draws per fit, 10^4 adaptive \
                      iterations and 10^7 mutant evaluations; caps are preflighted and \
                      exhaustion is BudgetExhausted/Unknown, never composition or deployment \
                      authority",
            versions: "fs-vmanifest schema v2; LoopStageGraph v1; EvidenceRetentionReceipt v1; \
                       CalibrationGraph v1; UncertaintyLedger v1; DeploymentGate v1; the \
                       consumed I05/I06/I10/I11 v1 manifest digests are receipt-bound at every \
                       stage boundary; rust-toolchain.toml and constellation.lock receipt-bound",
            capabilities: "safe Rust; no network or FFI in adjudication; deterministic mode; \
                           content-addressed receipt store mandatory; no stage skipping, no \
                           implicit identity bridging, no silent calibration bypass, no \
                           cross-stage holdout access, no favorable missing-receipt fill and \
                           no automatic physical actuation; maximal theorem/physical lanes \
                           feature-gated",
        },
        claims: rl_claims(),
        fixtures: rl_fixtures(),
        obligations: rl_obligations(),
        waivers: rl_waivers(),
        amendment_rules: "Every changed stage/identity/partition/ledger/calibration/threshold/ \
                          receipt/gate/seed/checker/profile field creates the exact next \
                          manifest version through FrozenManifest::amend. Candidate and blind \
                          partitions freeze before evidence. The amendment record names \
                          affected predecessor authority; no result edits this version, and \
                          unchanged evidence is rebound only through authenticated \
                          byte-identical lineage. Instance-manifest amendments propagate: a \
                          consumed I05/I06/I10/I11 successor version invalidates every RL \
                          receipt that bound the predecessor digest.",
    }
}

#[allow(clippy::too_many_lines)]
fn rl_claims() -> Vec<ClaimSpec> {
    vec![
        ClaimSpec {
            id: "rl-crossstage-identity-continuity",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "One identity spine spans the loop: every lot, specimen, coupon, \
                        instrument, calibration, campaign, model, target and receipt identity \
                        referenced at a stage boundary resolves to the same typed authority the \
                        producing stage minted, bound by exact content address; conflicting, \
                        re-minted, spliced or dangling boundary references fail closed",
            hypotheses: &[
                "stage boundaries consume predecessor artifacts only through \
                 EvidenceRetentionReceipt content addresses; display names, file paths and \
                 timestamps carry no identity authority",
                "identity semantics inside a stage belong to its instance manifest; the portfolio \
                 claim is continuity of the SAME authority across boundaries, never a re-derivation",
            ],
            qoi: "crossstage_identity_continuity_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G0,
            oracle: OracleRoute {
                identity: "independent loop-graph decoder and boundary-binding checker over \
                           separately generated valid/spliced/re-minted loop corpora",
                independent: true,
                tcb_overlap: "shares content-address primitives and public receipt schemas, not \
                              production stage assembly or boundary resolution",
            },
            activation: "the LoopStageGraph v1 implementation leaf opens",
            kill: "one boundary alias, accepted splice, silently re-minted identity or dangling \
                   consumed receipt blocks every downstream RL claim",
            fallback: "quarantine the affected loop world with typed BrokenSpine; no stage \
                       downstream of the break may publish loop authority",
            no_claim: "identity continuity does not prove any stage's internal correctness, \
                       physical custody, supplier honesty or instrument health — those belong to \
                       the instance manifests",
        },
        ClaimSpec {
            id: "rl-blind-partition-custody",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Blind partitions frozen at campaign-compile time survive the whole loop: \
                        no stage — lot selection, coupon fabrication records, DAQ streams, model \
                        fitting, refinement, deployment — can read another stage's held-out data \
                        or labels; custody transitions are typed and logged; any leak, premature \
                        read or post-result repartition is IntegrityFailed",
            hypotheses: &[
                "partition membership, label stores and truth parameters live in sealed custody \
                 with named stage-local consumers; the loop graph declares every authorized \
                 access edge before evidence exists",
                "side channels are in scope: fabrication metadata, adaptive design feedback, \
                 diagnostic logs and failure bundles must not encode held-out labels",
            ],
            qoi: "endtoend_blindness_custody_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "independent custody-transition replayer and leak detector over loop \
                           corpora with seeded overt and side-channel leaks",
                independent: true,
                tcb_overlap: "shares sealed-store schemas only; production stage code and logging \
                              paths are not reused",
            },
            activation: "identity continuity is green for the loop world",
            kill: "one undetected seeded leak, unauthorized custody edge or label-bearing side \
                   channel refutes loop blindness and blocks every gain/discrimination claim \
                   downstream",
            fallback: "burn the contaminated partition permanently, mark dependent claims \
                       EvidencePartial, and re-freeze successors on fresh ranges",
            no_claim: "custody discipline does not prove statistical validity of any analysis \
                       nor that synthetic blinding matches a real laboratory's information flow",
        },
        ClaimSpec {
            id: "rl-uncertainty-ownership-conservation",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Every uncertainty component in the loop — lot/spatial/aging aleatory, \
                        parameter epistemic, measurement, model discrepancy, timing — has exactly \
                        one owner at every stage; boundary handoffs neither drop nor double-count \
                        a component; the deployed twin's reported uncertainty decomposes exactly \
                        into upstream owned components under the frozen interval arithmetic",
            hypotheses: &[
                "the UncertaintyLedger binds each component to its owning stage and \
                 instance-manifest authority; a component with zero or two owners at any boundary is a \
                 frozen refusal, not a warning",
                "decomposition exactness is claimed under outward-rounded interval arithmetic on \
                 synthetic loop worlds with retained truth, not as a statement about real physics",
            ],
            qoi: "uncertainty_ownership_conservation_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G2,
            oracle: OracleRoute {
                identity: "independent ledger decoder plus exact interval recomposition checker \
                           against retained synthetic truth decompositions",
                independent: true,
                tcb_overlap: "shares frozen ledger bytes only; propagation and recomposition are \
                              independently implemented",
            },
            activation: "identity continuity and the ledger corpus are frozen",
            kill: "one dropped, doubly-owned or silently collapsed component, or a twin interval \
                   that fails exact recomposition, refutes ownership conservation",
            fallback: "widen the affected twin outputs to the conservative envelope of all \
                       candidate decompositions and mark them OwnershipUnknown",
            no_claim: "conserved ownership is bookkeeping soundness; it does not validate any \
                       component's magnitude, distribution family or physical origin",
        },
        ClaimSpec {
            id: "rl-calibration-measurement-graph",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "The calibration/measurement graph from raw DAQ counts to SI property \
                        observations is frozen, versioned and replayable; every deployed-twin QoI \
                        traces through concrete graph edges; a broken, expired, superseded or \
                        cycle-bearing calibration edge quarantines exactly its downstream \
                        evidence cone",
            hypotheses: &[
                "graph edges bind instrument, calibration certificate, transfer function, \
                 covariance and validity interval identities from the I11 and I06 authorities; \
                 edge validity is time-and-conditions qualified",
                "quarantine is cone-exact relative to the declared graph; graph-world completeness \
                 is not claimed at this tier",
            ],
            qoi: "calibration_graph_trace_and_quarantine_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G2,
            oracle: OracleRoute {
                identity: "independent graph decoder, edge-validity checker and exhaustive \
                           quarantine-cone oracle on bounded mutation graphs",
                independent: true,
                tcb_overlap: "shares typed graph bytes, not production traversal or cache code",
            },
            activation: "the CalibrationGraph v1 corpus and loop worlds are frozen",
            kill: "one QoI served without a complete valid trace, an under-quarantined cone or a \
                   stale edge accepted without authenticated lineage refutes the graph authority",
            fallback: "quarantine the whole loop world's measurement side conservatively and \
                       serve typed StaleCalibration refusals",
            no_claim: "a frozen graph does not calibrate any real instrument, prove traceability \
                       to a national standard or validate a transfer function's physics",
        },
        ClaimSpec {
            id: "rl-core-loop-replay",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "The complete core loop — lot passport, coupon design, campaign execution, \
                        model fit, twin update, deployment gate — replays deterministically from \
                        its EvidenceRetentionReceipt chain to bit-identical terminal states, \
                        receipts and served QoI enclosures on the same ISA/toolchain fingerprint",
            hypotheses: &[
                "every stage records exact command, seeds, budgets, versions, capability \
                 fingerprint and consumed receipt addresses; replay consumes only the chain",
                "inaccessible, redacted or expired inputs constrain replay explicitly and are \
                 named in the replay verdict rather than silently substituted",
            ],
            qoi: "whole_loop_bitwise_replay_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G5,
            oracle: OracleRoute {
                identity: "independent replay driver and receipt-chain verifier comparing full \
                           terminal artifacts bitwise",
                independent: true,
                tcb_overlap: "shares the retention store format only; no production stage \
                              orchestration code",
            },
            activation: "all four consumed instance-manifest digests are bound and a loop world \
                         completes end to end",
            kill: "any bit divergence on identical fingerprints, hidden input, or replay that \
                   silently proceeds past a missing artifact refutes loop replay",
            fallback: "demote the loop world's authority to EvidencePartial and retain the first \
                       divergence as a FailureBundle",
            no_claim: "bitwise replay proves process determinism, not correctness; cross-ISA \
                       behavior is governed by the instance manifests' own determinism classes",
        },
        ClaimSpec {
            id: "rl-substitution-selective-reproof",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "A mid-loop change — new lot, recalibrated sensor, model version bump, \
                        target or toolchain change — invalidates exactly the declared-graph \
                        reachable downstream receipts; unaffected byte-identical authority is \
                        reusable only through authenticated lineage; the twin cannot serve a QoI \
                        whose support chain contains an invalidated receipt",
            hypotheses: &[
                "the loop dependency graph types every receipt edge with semantic role, \
                 activation predicate and version; exactness is relative to this declared graph \
                 per the I06 doctrine",
                "reproof after invalidation re-enters at the earliest affected stage and rebinds \
                 the spine; a receipt is never revalidated because its path or name is unchanged",
            ],
            qoi: "selective_reproof_reachability_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "independent graph decoder and exhaustive reachability oracle on \
                           bounded loop mutations plus before/after authority digests",
                independent: true,
                tcb_overlap: "shares typed graph bytes, not production invalidation or cache code",
            },
            activation: "a loop world with complete declared dependency graph exists",
            kill: "one under-invalidated reachable receipt, an overclaimed reusable changed \
                   receipt or a served QoI with invalidated support refutes selective reproof",
            fallback: "invalidate the entire loop-world evidence closure conservatively and rerun \
                       from the changed stage",
            no_claim: "declared-graph exactness is not proof that every real organizational or \
                       physical dependency was modeled — that completeness is the theorem lane's \
                       separately gated target",
        },
        ClaimSpec {
            id: "rl-deployment-gate-composition",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "The deployed twin serves a QoI only when every receipt in its support \
                        chain is green, in-domain, unexpired and version-bound: composition is \
                        weakest-wins across the chain, and anything less yields a typed refusal \
                        — StaleEvidence, OutOfDomain, InvalidatedSupport, OwnershipUnknown or \
                        Unknown — never a silently degraded answer",
            hypotheses: &[
                "the DeploymentGate evaluates the full support chain at serve time against the \
                 frozen manifest digests, domains and expiry policies; gate decisions are \
                 receipts themselves and replay exactly",
                "refusal taxonomy is closed and typed; new failure modes require a manifest \
                 amendment, not an ad-hoc string",
            ],
            qoi: "deployment_gate_weakest_wins_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G2,
            oracle: OracleRoute {
                identity: "independent gate re-evaluator over loop worlds with seeded stale/ \
                           out-of-domain/invalidated support chains",
                independent: true,
                tcb_overlap: "shares gate receipt schemas only; the production gate evaluator is \
                              not reused",
            },
            activation: "loop replay and selective reproof are green for the world",
            kill: "one served QoI with a non-green support link, an untyped refusal or a gate \
                   decision that does not replay refutes deployment gating",
            fallback: "fail the gate closed for the affected QoI set and serve typed refusals \
                       with the ranked missing-evidence plan",
            no_claim: "a green gate is chain bookkeeping, not physical validity, fitness for any \
                       safety function or accreditation of laboratories, suppliers, devices, \
                       targets or models outside the declared ContextOfUse",
        },
        ClaimSpec {
            id: "rl-adaptive-oed-loop",
            ambition: Ambition::Frontier,
            polarity: ClaimPolarity::Affirmative,
            statement: "The adaptive loop — fit, design the next coupon/campaign, execute, refit \
                        — preserves anytime-valid statistical guarantees and partition custody \
                        across iterations under optional stopping and design adaptation; \
                        information gain accrues to preregistered per-iteration ledgers and the \
                        null false-certification bound holds",
            hypotheses: &[
                "iteration filtration, adaptation policy, stopping rules and gain ledgers freeze \
                 before the first adaptive iteration; the I10 discrimination and I11 campaign \
                 authorities are consumed, not re-derived",
                "design feedback is a declared custody edge: adaptive designs may depend on \
                 development data and prior posteriors, never on held-out labels",
            ],
            qoi: "anytime_valid_upper_bound_on_null_false_certification_probability",
            unit: "1",
            tolerance: ToleranceSemantics::Interval { lo: 0.0, hi: 0.05 },
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "independent e-process implementation and null loop simulations with \
                           adversarial stopping and adaptation",
                independent: true,
                tcb_overlap: "shares retained canonical iteration records and frozen policy bytes \
                              only",
            },
            activation: "feature rl-adaptive-loop; baseline custody, ownership and replay lanes \
                         are green",
            kill: "null false-certification above 0.05, a label-dependent design, edited stopping \
                   policy or double-counted iteration refutes the adaptive loop",
            fallback: "freeze adaptation and fall back to the preregistered static campaign \
                       sequence with descriptive-only adaptive summaries",
            no_claim: "anytime validity is under the frozen null and filtration; it neither \
                       proves the adapted designs are good nor that the model family contains \
                       the truth",
        },
        ClaimSpec {
            id: "rl-hil-timing-composition",
            ambition: Ambition::Frontier,
            polarity: ClaimPolarity::Affirmative,
            statement: "HIL, timing and fault-injection receipts from the I05 target profiles \
                        join the twin's support chain with target-exact scope: a receipt binds \
                        its exact binary, target, clock/power state and interference assumptions, \
                        composes with the numeric loop through the same weakest-wins gate, and \
                        never widens to another target, profile or binary",
            hypotheses: &[
                "target fingerprints are part of receipt identity; the gate treats a fingerprint \
                 mismatch as OutOfDomain, not as a warning",
                "timing evidence kinds keep their I05 semantics — a measured sample maximum never \
                 promotes to a WCET bound anywhere in the chain",
            ],
            qoi: "target_exact_timing_composition_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "independent chain checker over loop worlds with seeded cross-target \
                           and cross-profile receipt substitutions",
                independent: true,
                tcb_overlap: "shares receipt schemas and target fingerprint format only",
            },
            activation: "feature rl-hil-composition; deployment gating is green and an I05 \
                         profile is bound",
            kill: "one accepted cross-target substitution, promoted sample maximum or widened \
                   fault receipt refutes timing composition",
            fallback: "serve the affected QoIs with the numeric-only chain and a typed \
                       TimingEvidenceMissing refusal for timing-qualified requests",
            no_claim: "composition discipline mints no new timing bounds and no safety \
                       authority; per-target validity remains entirely the I05 manifests'",
        },
        ClaimSpec {
            id: "rl-physical-campaign-parity",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Affirmative,
            statement: "With the governed physical packs admitted, the full loop runs on real \
                        lots, coupons and rigs with the SAME receipt, blindness, ownership, \
                        replay and gating authority as synthetic worlds — receipt-machinery \
                        parity is exact, and synthetic-to-physical transfer quality is measured \
                        against preregistered ledgers, never assumed",
            hypotheses: &[
                "physical-world stages emit the same typed receipts through the same boundaries; \
                 any physical-only escape hatch, manual bypass or untyped human step is a parity \
                 break, not local color",
                "transfer metrics (calibration residuals, coverage deltas, gain shortfalls) are \
                 preregistered descriptive ledgers; no synthetic result is re-labeled as physical \
                 evidence",
            ],
            qoi: "physical_receipt_parity_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G5,
            oracle: OracleRoute {
                identity: "independent receipt-diff checker between synthetic and physical loop \
                           worlds plus a governed-pack admission auditor",
                independent: true,
                tcb_overlap: "shares receipt schemas and the retention store format only",
            },
            activation: "feature rl-physical-loop; the rl-governed-physical-loop-pack waiver is \
                         discharged and every baseline lane is independently green",
            kill: "one physical-only bypass, unreceipted manual step, parity-breaking schema \
                   fork or synthetic result presented as physical refutes parity",
            fallback: "keep physical campaigns descriptive-only; loop authority remains synthetic \
                       with the waiver's promotion effect in force",
            no_claim: "receipt parity does not make physical results correct, safe or \
                       representative; laboratory accreditation, device qualification and \
                       supplier trust stay out of scope",
        },
        ClaimSpec {
            id: "rl-endtoend-composition-theorem",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Affirmative,
            statement: "A machine-checked theorem composes the frozen stage contracts: if every \
                        receipt in a served QoI's support chain is green under its instance \
                        manifest's semantics, then the end-to-end property — the served \
                        enclosure's validity under the conjunction of declared stage premises — \
                        holds; premises, adapters and open-world boundaries are explicit and \
                        kernel replay is retained",
            hypotheses: &[
                "a pre-proof successor freezes stage-contract, receipt, chain and property ASTs \
                 with canonical bytes, total runtime-premise mapping and deterministic proof \
                 translation",
                "exact axiom allowlist and transitive closure, independent model decoder and \
                 nonvacuity witnesses are receipt-bound; manual and physical steps are explicit \
                 open-world boundaries, never silently assumed away",
            ],
            qoi: "kernel_checked_composition_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G5,
            oracle: OracleRoute {
                identity: "pinned Lean kernel replay plus independent canonical-model translation \
                           and finite counterexample checker over bounded loop worlds",
                independent: true,
                tcb_overlap: "kernel/axioms and the frozen semantic model are declared TCB; no \
                              production gate or chain code",
            },
            activation: "feature rl-composition-theorem; a pre-proof successor freezes all \
                         machine artifacts and the baseline lanes are independently green",
            kill: "an open-world dependency without a boundary, missing premise, disallowed \
                   axiom, translation mismatch, kernel rejection or genuine counterexample \
                   leaves this claim Unknown/Refuted and blocks composition wording",
            fallback: "retain weakest-wins gate bookkeeping as the only composition authority, \
                       exactly as in the baseline lane",
            no_claim: "version-1 prose mints no theorem authority; a composed guarantee is \
                       conditional on every stage premise and is not validity of any premise \
                       itself",
        },
        ClaimSpec {
            id: "rl-forged-loop-certificate-falsifier",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Refutation,
            statement: "Hidden cross-stage mutants try to obtain a green end-to-end certificate \
                        for a known invalid loop: any accepted mutant refutes the loop certifier \
                        and is retained as a minimized counterexample",
            hypotheses: &[
                "mutants include stage-skip splices, receipt replay across lots and worlds, \
                 partition leaks through fabrication metadata and diagnostic side channels, \
                 calibration-edge bypass, invalidated-support serving, cross-target timing \
                 substitution, adaptive-stopping edits, ledger double-ownership and \
                 retention-store truncation presented as complete",
                "each mutant has an independently authored lineage, custody, reachability, \
                 fingerprint or arithmetic witness hidden until maximal adjudication",
            ],
            qoi: "false_certificate_count",
            unit: "1",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "independent mutation generator and witness checker with no production \
                           stage, gate, graph or replay code",
                independent: true,
                tcb_overlap: "public receipt schemas and content-address primitives only",
            },
            activation: "a maximal certifier candidate and all checker/version identities are \
                         sealed before mutant reveal",
            kill: "the intended lane succeeds when a false certificate is accepted: mark the \
                   certifier Refuted and block dependent promotion; zero accepted mutants is \
                   finite falsification evidence only",
            fallback: "disable the affected certifier and fall back to the strongest \
                       independently surviving baseline gate with typed refusals",
            no_claim: "a finite adversarial corpus cannot prove absence of forgery, leakage, \
                       bookkeeping error or composition unsoundness in general",
        },
    ]
}

#[allow(clippy::too_many_lines)]
fn rl_fixtures() -> Vec<FixturePin> {
    vec![
        FixturePin {
            id: CAMPAIGN_POLICY_FIXTURE,
            source: FixtureSource::AuthoredSpec {
                spec: "RL_CAMPAIGN_POLICY_V1\n\
IDENTITY=Lot,Specimen,Coupon,Instrument,Calibration,Campaign,Model,Target,Receipt,LoopWorld,Gate and Checker identities have distinct typed domains, canonical preimages, versions and validity; stage boundaries bind exact content addresses; display text, paths and timestamps carry no identity authority; a signature authenticates a receipt, not the loop that produced it\n\
STAGES=the loop graph declares every stage, artifact, receipt and authorized access edge before evidence; stage internals answer to their instance manifests (I05 deployment/timing, I06 material-lot, I10 identifiability/design, I11 campaign/DAQ); the portfolio owns only the seams; no stage skipping or implicit identity bridging\n\
BLINDNESS=partitions frozen at campaign-compile time survive the whole loop; no stage may read another stage's holdout; fabrication metadata, adaptive feedback, diagnostics and failure bundles are leak surfaces in scope; a leak burns the partition permanently\n\
UNCERTAINTY=every component {lot/spatial/aging aleatory,parameter epistemic,measurement,model discrepancy,timing} has exactly one owner at every stage; handoffs neither drop nor double-count; the twin's reported uncertainty decomposes exactly into upstream owned components under outward-rounded interval arithmetic\n\
CALIBRATION=raw DAQ counts reach SI observations only through the frozen versioned calibration/measurement graph; broken,expired,superseded or cyclic edges quarantine exactly their downstream cone; no silent calibration bypass\n\
REPLAY=the whole core loop replays bitwise from its EvidenceRetentionReceipt chain on the same ISA/toolchain fingerprint; inaccessible,redacted or expired inputs constrain replay explicitly and are named in the verdict\n\
REPROOF=a mid-loop change invalidates exactly the declared-graph reachable receipts; unchanged authority is reusable only through authenticated byte-identical lineage; reproof re-enters at the earliest affected stage; a consumed instance-manifest successor version invalidates every receipt bound to the predecessor digest\n\
DEPLOYMENT=the twin serves a QoI only on an all-green,in-domain,unexpired,version-bound support chain; composition is weakest-wins; refusals are typed {StaleEvidence,OutOfDomain,InvalidatedSupport,OwnershipUnknown,Unknown}; no silently degraded answer; no automatic physical actuation\n\
ADAPTIVE=iteration filtration,adaptation policy,stopping rules and gain ledgers freeze before the first adaptive iteration; adaptive designs may consume development data and prior posteriors, never held-out labels; anytime-valid bounds hold under optional stopping\n\
TIMING=HIL/timing/fault receipts are target-exact: exact binary,target,clock/power state and interference assumptions; fingerprint mismatch is OutOfDomain; a measured sample maximum never promotes to a WCET bound anywhere in the chain\n\
THEOREM_AUTHORITY=version 1 prose mints no composition or parity authority; pre-proof successors freeze stage-contract/receipt/chain/property ASTs, semantics, total runtime premises, deterministic AST-to-Lean translation, exact axiom allowlist {propext,Quot.sound,Classical.choice} and transitive closure; sorryAx, custom postulates and native/unsafe oracle shortcuts are IntegrityFailed; manual and physical steps are explicit open-world boundaries\n\
EVIDENCE_STATES=Execution{Succeeded,Failed,Cancelled,TimedOut,BudgetExhausted,InfrastructureFailed,IntegrityFailed}; Predicate{Accepted,Refuted,Unknown,NotEvaluated}; Claim{NoClaim,EvidencePartial,Reproduced,Refuted,Proved}; Gate{Served,StaleEvidence,OutOfDomain,InvalidatedSupport,OwnershipUnknown,Unknown}; Promotion{Blocked,CoreEligible,MaxEligible}; axes never substitute\n\
HOLDOUT=loop worlds,truth parameters,labels,seeds,checker/policy versions and mutant corpora freeze before access; each Core/Max heldout fixture has one named consumer stage and disjoint range; premature/cross-stage read,label leak,replacement,retry-after-result or post-result tuning is IntegrityFailed\n\
LIFECYCLE=request->drain->finalize; poll at stage,boundary,receipt,iteration,gate and proof/search tile boundaries; checkpoints bind manifest/candidate/evidence cutoff and completed logical ids; resume/fork cannot duplicate stage executions,iterations,wealth or served answers\n\
LOGGING=bounded schema-versioned fs-obs JSONL with stable run/world/stage/case/claim/fixture/leaf/receipt/gate/checker/attempt ids, exact units,seeds,budgets,versions,capabilities and redaction reason; stdout is not evidence\n\
RETENTION=manifest,exact command,code/contract/schema/toolchain/registry hashes,consumed instance-manifest digests,receipt chains,ledgers,graphs,gate decisions,wealth paths,all terminal states,minimized counterexamples and replay-verifier result; promotion/refutation evidence durable; inaccessible/redacted/expired inputs constrain replay and never become silently verified\n\
FAILURE_BUNDLE=first boundary/custody/ownership/calibration/replay/gate divergence plus bounded context,identity chain,consumed receipts,expected/actual state,causal predecessors,minimized witness,checker disagreement,terminal state and replay command; partial success cannot publish normal authority\n\
PROMOTION=baseline claims only RL.G4 after G2 reproduction/G3 falsification; maximal claims only RL.G7 after G5 reproduction/G6 red-team; missing,stale,waived,integrity-failed or inaccessible evidence cannot promote\n\
LEAF_REQUIREMENT=every obligation references this policy fixture,all nine unit classes,smoke/core/max tier,DSR lane,events,replay,G4 drain/checkpoint,G5 matrix,independent adjudication,performance envelope,accessibility/agent parity and consuming gate",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "rl-loop-stage-graph-corpus",
            source: FixtureSource::AuthoredSpec {
                spec: "RL_FIXTURE_V1 typed loop-graph corpus: valid lot-to-deployment stage \
                       graphs across the four instance authorities with exact boundary bindings, \
                       plus invalid classes — stage-skip splice, re-minted identity, dangling \
                       consumed receipt, cross-world receipt replay, unauthorized custody edge, \
                       undeclared adaptive feedback edge and retention-store truncation. Graphs \
                       up to 10^4 receipts with exhaustive independent boundary/reachability \
                       witnesses.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "rl-synthetic-loop-worlds",
            source: FixtureSource::AuthoredSpec {
                spec: "RL_FIXTURE_V1 end-to-end synthetic loop worlds with retained truth: small \
                       material families (2..5 lots, 1..4 suppliers) feeding coupon designs, \
                       simulated campaigns with DAQ noise/nuisance, hierarchical fits, twin \
                       updates and served QoI sets; every stage emits full \
                       EvidenceRetentionReceipts; worlds include mid-loop changes (new lot, \
                       recalibration, model bump, target change) with exact expected \
                       invalidation cones and gate outcomes. Development seeds 0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "rl-calibration-chain-corpus",
            source: FixtureSource::AuthoredSpec {
                spec: "RL_FIXTURE_V1 calibration/measurement graphs: DAQ count to SI observation \
                       chains with instruments, certificates, transfer functions, covariance and \
                       validity intervals; valid chains plus broken, expired, superseded, cyclic, \
                       condition-violated and bypass-attempt variants, each with an exhaustive \
                       expected quarantine cone and an independent validity witness. Development \
                       seeds 4096..=8191.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "rl-uncertainty-ledger-corpus",
            source: FixtureSource::AuthoredSpec {
                spec: "RL_FIXTURE_V1 uncertainty-ownership ledgers over loop worlds: complete \
                       single-owner assignments with exact retained truth decompositions under \
                       outward-rounded interval arithmetic, plus trap classes — dropped \
                       component at a boundary, double ownership across stages, silent collapse \
                       of discrepancy into measurement, timing uncertainty smuggled into a \
                       numeric component and recomposition off-by-rounding cases. Development \
                       seeds 8192..=12287.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "rl-composition-theorem-card",
            source: FixtureSource::AuthoredSpec {
                spec: "RL_THEOREM_CARD_V1 TARGET. Intended successor proposition composes frozen \
                       stage contracts: for every served QoI, an all-green support chain under \
                       each instance manifest's semantics entails the served enclosure's \
                       validity under the conjunction of declared stage premises, with explicit \
                       adapters and open-world boundaries for manual and physical steps. \
                       REQUIRED MACHINE ARTIFACTS: stage-contract/receipt/chain/property ASTs \
                       and canonical bytes, total premise map, deterministic Lean \
                       translation/roundtrip, exact axiom allowlist \
                       {propext,Quot.sound,Classical.choice}, transitive closure, nonvacuity \
                       witnesses and retained kernel replay. Prose grants no composition \
                       authority.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "rl-physical-protocol-card",
            source: FixtureSource::AuthoredSpec {
                spec: "RL_PROTOCOL_CARD_V1 TARGET. Intended governed physical-loop protocol: \
                       consent/license/redaction policy, custody and key management for real \
                       lots/coupons/rigs, receipt-parity requirements (same typed receipts \
                       through the same boundaries, no physical-only bypass, every manual step \
                       receipted), blinded physical partitions, preregistered \
                       synthetic-to-physical transfer ledgers (calibration residuals, coverage deltas, gain \
                       shortfalls), third-party replay constraints and retention classes. \
                       Admission discharges the rl-governed-physical-loop-pack waiver and \
                       activates parity adjudication; the card itself grants no physical \
                       authority.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "rl-identity-custody-core-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "RL_FIXTURE_V1 CORE HOLDOUT identity/calibration adversaries, indices \
                       65536..=69631: sealed loop worlds with boundary splices, re-minted \
                       identities, replayed receipts, dangling bindings, expired/superseded/ \
                       cyclic calibration edges, condition-violated transfer functions and \
                       bypass attempts; each carries an independent lineage or graph witness. \
                       One RL.G3 consumer.",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "rl-blind-loop-core-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "RL_FIXTURE_V1 CORE HOLDOUT blind-loop worlds, indices 69632..=73727: \
                       sealed end-to-end worlds with seeded overt and side-channel partition \
                       leaks (fabrication metadata, adaptive feedback, diagnostic logs, failure \
                       bundles), unauthorized custody edges, ownership drop/double-count traps \
                       and recomposition stressors; leak labels and truth decompositions are \
                       withheld until one-shot submission. One RL.G3 consumer.",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "rl-replay-invalidation-core-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "RL_FIXTURE_V1 CORE HOLDOUT replay/reproof/gate worlds, indices \
                       73728..=77823: sealed worlds with hidden replay divergence seeds, \
                       mid-loop changes whose exact invalidation cones are withheld, stale and \
                       invalidated-support serving traps, untyped-refusal bait and gate \
                       decisions that must replay bitwise; expected terminal states sealed \
                       until submission. One RL.G3 consumer.",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "rl-adaptive-hil-max-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "RL_FIXTURE_V1 MAX HOLDOUT adaptive/timing worlds, indices \
                       131072..=135167: adaptive loops under either-null with adversarial \
                       stopping, label-dependent-design bait and iteration double-count traps; \
                       plus timing-composition worlds with cross-target and cross-profile \
                       receipt substitutions, sample-maximum promotion bait and \
                       fingerprint-mismatch cases. Oracle certificates withheld until adjudication. One \
                       RL.G6 consumer.",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "rl-loop-certifier-mutants-max-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "RL_FIXTURE_V1 MAX HOLDOUT loop-certifier mutants, indices \
                       135168..=143359: stage-skip splices, cross-lot receipt replays, \
                       partition leaks via side channels, calibration bypass, \
                       invalidated-support serving, cross-target timing substitution, adaptive-stopping \
                       edits, ledger double-ownership, retention truncation presented as \
                       complete, physical-parity bypass and disallowed theorem axioms; each \
                       with an independent witness sealed until the certifier candidate is \
                       frozen. One RL.G6 consumer.",
            },
            partition: Partition::HeldOut,
        },
    ]
}

#[allow(clippy::too_many_lines)]
fn rl_obligations() -> Vec<ObligationRow> {
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
            leaf: "rl-identity-calibration-spine",
            claims_covered: &[
                "rl-calibration-measurement-graph",
                "rl-crossstage-identity-continuity",
            ],
            unit_cases: UNIT_CASES,
            g0: "generators: loop-graph and calibration-chain mutations with independent \
                 boundary/validity witnesses; predicates: content-address binding, boundary \
                 resolution, edge validity, quarantine-cone exactness; laws: roundtrip, no \
                 boundary alias, splice/replay refusal, cone monotonicity under edge removal \
                 and stale-edge non-acceptance without authenticated lineage; shrink \
                 graph/chain while preserving the witness",
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                "rl-calibration-chain-corpus",
                "rl-identity-custody-core-holdout",
                "rl-loop-stage-graph-corpus",
            ],
            g3_relations: &[
                "receipt and node id relabeling preserves boundary resolution and quarantine \
                 cones up to canonical rebinding",
                "removing a calibration edge can preserve or enlarge but cannot shrink its \
                 quarantine cone",
                "a spliced or replayed receipt fails at the boundary before any downstream \
                 semantic use",
                "byte-identical unrelated sibling receipts remain reusable only with \
                 authenticated lineage",
                "re-validating a boundary never re-derives units, frames or domains",
            ],
            g4_schedule: "cancel/panic/timeout/corrupt input at graph decode, boundary \
                          resolution, edge-validity check, cone computation and receipt \
                          binding; poll each boundary/edge tile; request->drain->finalize \
                          publishes complete spine/cone verdicts or nothing; checkpoint \
                          completed boundary ids and content roots; resume/fork cannot \
                          duplicate bindings; retain bounded FailureBundle and every terminal \
                          state",
            g5_matrix: "graph/chain shards {1,2,7} x input orders {capture,reverse,permuted} x \
                        deterministic mode on identical toolchain fingerprint; boundary \
                        verdicts, cones, diagnostics, event JSON and digests match bitwise",
            entry_point: "scripts/e2e/leapfrog/rl_identity_spine.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (rl-identity-spine slice)",
            obs_events: &[
                "request.received",
                "spine.boundary_bound",
                "spine.boundary_refused",
                "calibration.edge_validated",
                "calibration.cone_quarantined",
                "cancel.requested",
                "drain.completed",
                "finalize.completed",
                "failure_bundle.retained",
                "adjudication.receipt",
            ],
            replay_command: "scripts/e2e/leapfrog/rl_identity_spine.sh --manifest <manifest-id> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "rl-blind-custody-ownership",
            claims_covered: &[
                "rl-blind-partition-custody",
                "rl-uncertainty-ownership-conservation",
            ],
            unit_cases: UNIT_CASES,
            g0: "generators: loop worlds with seeded overt/side-channel leaks and ledger \
                 drop/double-count traps; predicates: custody-transition typing, \
                 authorized-edge closure, leak detection, single-owner closure and exact interval \
                 recomposition; laws: no stage reads another stage's holdout, leak burns the \
                 partition, ownership is exactly-once at every boundary and recomposition \
                 equals retained truth under outward rounding; shrink world while preserving \
                 the leak or ledger witness",
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                "rl-blind-loop-core-holdout",
                "rl-synthetic-loop-worlds",
                "rl-uncertainty-ledger-corpus",
            ],
            g3_relations: &[
                "stage and world id relabeling preserves custody verdicts and ledger \
                 decompositions up to relabeling",
                "adding an authorized custody edge cannot retroactively legalize a prior leak",
                "removing applicable evidence moves ownership verdicts only toward \
                 OwnershipUnknown, never toward green",
                "diagnostic and failure-bundle content never distinguishes held-out labels \
                 beyond the declared access edges",
                "splitting one owned component across stages without a declared handoff is \
                 refused, not averaged",
            ],
            g4_schedule: "cancel/panic/timeout at custody transition, leak scan, ledger join, \
                          recomposition and receipt binding; poll each stage/component tile; \
                          request->drain->finalize publishes one terminal custody+ownership \
                          state per world or nothing; checkpoint completed stage ids, ledger \
                          roots and scan positions; resume/fork equals one-shot; retain leak \
                          and recomposition FailureBundle",
            g5_matrix: "world/component shards {1,2,7} x workers {1,2,7} x stage orders \
                        {forward,reverse} x deterministic mode on identical ISA; custody \
                        verdicts, leak findings, ledger decompositions, tie breaks, events and \
                        receipts match bitwise",
            entry_point: "scripts/e2e/leapfrog/rl_blind_custody.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (rl-blind-custody slice)",
            obs_events: &[
                "request.received",
                "custody.transition",
                "custody.leak_detected",
                "custody.partition_burned",
                "ledger.component_owned",
                "ledger.recomposition_checked",
                "cancel.requested",
                "drain.completed",
                "finalize.completed",
                "failure_bundle.retained",
                "adjudication.receipt",
            ],
            replay_command: "scripts/e2e/leapfrog/rl_blind_custody.sh --manifest <manifest-id> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "rl-replay-reproof-gate",
            claims_covered: &[
                "rl-core-loop-replay",
                "rl-deployment-gate-composition",
                "rl-substitution-selective-reproof",
            ],
            unit_cases: UNIT_CASES,
            g0: "generators: complete loop worlds with receipt chains, mid-loop mutations with \
                 exact expected cones and gate scenarios with seeded non-green links; \
                 predicates: chain completeness, bitwise replay, reachability exactness, \
                 authenticated reuse, weakest-wins evaluation and refusal typing; laws: replay \
                 determinism on identical fingerprints, under/over-invalidation refusal, no \
                 serve on invalidated support, gate decisions replay and refusals are typed \
                 members of the closed taxonomy; shrink world/mutation while preserving the \
                 divergence witness",
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                "rl-loop-stage-graph-corpus",
                "rl-replay-invalidation-core-holdout",
                "rl-synthetic-loop-worlds",
            ],
            g3_relations: &[
                "receipt-chain reordering within declared concurrency preserves terminal states \
                 and served enclosures bitwise",
                "a mutation's invalidation cone is exactly the declared-graph reachable set — \
                 node relabeling conjugates it",
                "adding a dependency edge can preserve or enlarge but cannot shrink an \
                 invalidation cone",
                "an instance-manifest digest bump invalidates every receipt bound to the \
                 predecessor digest",
                "weakest-wins is monotone: upgrading one link never downgrades a served \
                 verdict, degrading one link never upgrades it",
            ],
            g4_schedule: "cancel/panic/timeout/corrupt checkpoint at chain decode, stage replay, \
                          cone computation, reuse check, gate evaluation and receipt binding; \
                          poll each stage/frontier tile; request->drain->finalize publishes one \
                          terminal replay+reproof+gate state per world; checkpoint bound \
                          receipt roots, frontier and gate decisions; resume/fork equals \
                          one-shot; retain first-divergence and false-serve FailureBundle",
            g5_matrix: "world/chain shards {1,2,7} x workers {1,2,7} x replay orders \
                        {chain,reverse,permuted-within-concurrency} x deterministic mode on \
                        identical ISA/toolchain; terminal states, cones, gate decisions, \
                        refusal types, events and receipts match bitwise",
            entry_point: "scripts/e2e/leapfrog/rl_replay_reproof.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (rl-replay-reproof slice)",
            obs_events: &[
                "request.received",
                "replay.stage_completed",
                "replay.divergence",
                "reproof.receipt_invalidated",
                "reproof.authority_reused",
                "gate.served",
                "gate.refused",
                "cancel.requested",
                "drain.completed",
                "finalize.completed",
                "failure_bundle.retained",
                "adjudication.receipt",
            ],
            replay_command: "scripts/e2e/leapfrog/rl_replay_reproof.sh --manifest <manifest-id> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "rl-adaptive-hil-loop",
            claims_covered: &["rl-adaptive-oed-loop", "rl-hil-timing-composition"],
            unit_cases: UNIT_CASES,
            g0: "generators: adaptive loop iterations under either-null with stopping \
                 adversaries plus timing-composition worlds with cross-target substitutions; \
                 predicates: filtration/policy identity, nonnegative wealth, custody of design \
                 feedback, gain-ledger arithmetic, target-fingerprint binding and evidence-kind \
                 preservation; laws: label-independent adaptation, null stopping preserves the \
                 bound, iteration exactly-once, fingerprint mismatch is OutOfDomain and sample \
                 maxima never promote; shrink iteration/world while preserving the error",
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                "rl-adaptive-hil-max-holdout",
                "rl-synthetic-loop-worlds",
                "rl-uncertainty-ledger-corpus",
            ],
            g3_relations: &[
                "adversarial optional stopping under either null preserves the \
                 false-certification bound <= 0.05",
                "an adaptive design computed with held-out labels is detected and refused \
                 before execution",
                "iteration batching under the frozen filtration preserves exact wealth and \
                 ledger state",
                "substituting a timing receipt from another target, profile or binary is \
                 OutOfDomain before composition",
                "a measured sample maximum composed anywhere in the chain never upgrades to a \
                 WCET-bearing verdict",
            ],
            g4_schedule: "whole-campaign preflight freezes iteration/scenario budgets before max \
                          holdout access; cancel/panic/timeout at fit, design, execution, \
                          refit, wealth update, chain join and receipt binding; \
                          request->drain->finalize persists exact iteration and e-process state \
                          and one terminal state; checkpoint filtration position, wealth, \
                          ledger and chain roots; resume cannot double-count or edit stopping \
                          policy; BudgetExhausted stays Unknown; retain false-certification and \
                          promotion-bait FailureBundle",
            g5_matrix: "iteration/world shards {1,2,7,31} x workers {1,2,7} x orders \
                        {lex,reverse,permuted} x deterministic mode on identical ISA; wealth \
                        paths, ledgers, compositions, terminal states, events and receipts \
                        match bitwise",
            entry_point: "scripts/e2e/leapfrog/rl_adaptive_hil.sh",
            tier: CampaignTier::Max,
            dsr_lane: "dsr quality --tool frankensim (rl-adaptive-hil isolated lane)",
            obs_events: &[
                "request.received",
                "adaptive.iteration_completed",
                "adaptive.wealth_updated",
                "adaptive.design_refused",
                "timing.receipt_composed",
                "timing.out_of_domain",
                "cancel.requested",
                "drain.completed",
                "finalize.completed",
                "failure_bundle.retained",
                "adjudication.receipt",
            ],
            replay_command: "scripts/e2e/leapfrog/rl_adaptive_hil.sh --manifest <manifest-id> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "rl-loop-certifier",
            claims_covered: &[
                "rl-endtoend-composition-theorem",
                "rl-forged-loop-certificate-falsifier",
                "rl-physical-campaign-parity",
            ],
            unit_cases: UNIT_CASES,
            g0: "generators: stage-contract/receipt/chain ASTs, bounded loop counterexample \
                 searches, receipt-parity diffs and hidden cross-stage mutants; predicates: AST \
                 canonicality, premise/adapter/boundary closure, translation roundtrip, axiom \
                 closure, kernel replay, parity-diff emptiness and mutant witnesses; laws: \
                 alpha-renaming preserves the composition theorem up to rebinding, an \
                 open-world step without a boundary fails closed, parity breaks refuse and \
                 every known false certificate is refused; shrink model/mutant while \
                 preserving the witness",
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                "rl-composition-theorem-card",
                "rl-governed-physical-loop-pack",
                "rl-loop-certifier-mutants-max-holdout",
                "rl-physical-protocol-card",
            ],
            g3_relations: &[
                "stage/receipt alpha-renaming preserves the composition verdict up to canonical \
                 rebinding",
                "a manual or physical step without an explicit open-world boundary fails closed \
                 before kernel submission",
                "finite countermodels force Unknown for chains whose premises do not entail the \
                 served property",
                "a physical-only bypass or unreceipted manual step breaks parity before any \
                 transfer ledger is read",
                "every accepted splice/replay/leak/bypass/substitution/edit mutant refutes its \
                 certifier",
            ],
            g4_schedule: "whole-campaign preflight freezes proof/search/mutant budgets before \
                          reveal; cancel/panic/timeout/corrupt checkpoint at AST decode, \
                          premise mapping, translation, kernel replay, parity diff, mutant \
                          adjudication and receipt binding; request->drain->finalize drains \
                          children and publishes only complete terminal states; checkpoint \
                          frontier/proof/model/diff roots; BudgetExhausted stays Unknown; \
                          retain countermodel, parity-break and false-certificate \
                          FailureBundle",
            g5_matrix: "model/proof/mutant shards {1,2,7,31} x workers {1,2,7} x orders \
                        {lex,reverse,permuted} x deterministic mode on identical \
                        kernel/toolchain; canonical ASTs, premise closures, kernel results, \
                        parity diffs, mutant verdicts, events and receipt roots match bitwise",
            entry_point: "scripts/e2e/leapfrog/rl_loop_certifier.sh",
            tier: CampaignTier::Max,
            dsr_lane: "dsr quality --tool frankensim (rl-loop-certifier isolated lane)",
            obs_events: &[
                "request.received",
                "composition.premises_mapped",
                "composition.kernel_checked",
                "composition.unknown",
                "parity.diff_empty",
                "parity.break_detected",
                "certifier.mutant_refused",
                "certifier.refuted",
                "cancel.requested",
                "drain.completed",
                "finalize.completed",
                "failure_bundle.retained",
                "adjudication.receipt",
            ],
            replay_command: "scripts/e2e/leapfrog/rl_loop_certifier.sh --manifest <manifest-id> --replay <artifact-id>",
        },
    ]
}

fn rl_waivers() -> Vec<Waiver> {
    vec![Waiver {
        subject: "rl-governed-physical-loop-pack",
        reason: "all committed RL fixtures are synthetic loop worlds; independently authorized \
                 physical campaigns — real lots, coupons, rigs, DAQ streams, calibration \
                 certificates and blind lab truth spanning the whole loop — are not yet \
                 licensed, consent-cleared and content-addressed for replay",
        owner: "RL portfolio governance with the I05/I06/I10/I11 laboratory, metrology and \
                V&V owners",
        predicate: "a governed physical-loop pack is admitted under the \
                    rl-physical-protocol-card: exact consent/license/redaction policy, custody \
                    and key management, receipt-parity requirements with no physical-only \
                    bypass, blinded physical \
                    partitions, preregistered transfer ledgers, checker ownership, retention \
                    and third-party replay constraints",
        expiry: "before the first RL maximal physical-campaign submission; review on every \
                 manifest amendment, consumed instance-manifest amendment and laboratory/ \
                 license/calibration/rig revision",
        promotion_effect: "synthetic loop engineering may proceed, but no physical-loop \
                           authority, real-lot deployment claim, lab-validated composition or \
                           maximal RL promotion may be claimed while this waiver is live; \
                           baseline synthetic authority remains separately eligible",
    }]
}
