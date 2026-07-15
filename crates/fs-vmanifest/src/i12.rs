//! The I12 (hybrid mode-automaton compiler) VerificationManifest draft
//! (bead frankensim-leapfrog-2026-program-i94v.1.5.8.1).
//!
//! Baseline lattice ([S]): typed modes/guards/resets/invariants,
//! per-mode structural consistency, class-qualified complete guard-root
//! accounting on prescribed segments, deterministic superdense-time
//! sequencing, and event-time/generalized derivative protocols. Maximal
//! lattice ([F]/[M]): set-valued simultaneous-order semantics,
//! topology-changing resets, complete bounded-horizon mode coverage
//! with localized Unknown, and certified Zeno/accumulation handling.
//! One preregistered refutation lane hunts the doctrine-D3 Sev-0
//! failure (a fake finite certificate at a grazing tangency). A weaker
//! receipt closes its own element and is never relabeled as the
//! stronger theorem.

use crate::{
    Ambition, CampaignTier, ClaimPolarity, ClaimSpec, FiveExplicits, FixturePin, FixtureSource,
    GauntletTier, ManifestDraft, ObligationRow, OracleRoute, Partition, ToleranceSemantics, Waiver,
};

/// Build the I12 draft. Consumers freeze it themselves; the conformance
/// battery proves it freezes.
#[must_use]
pub fn i12_draft() -> ManifestDraft {
    ManifestDraft {
        initiative: "I12",
        title: "Hybrid mode-automaton compiler gate: typed automata to deterministic \
                superdense runtime with class-qualified event completeness, derivative \
                protocols, and set-valued nonsmooth outcomes",
        version: 1,
        explicits: FiveExplicits {
            units: "SI base units throughout; times in seconds, guard values in the \
                    declared per-fixture unit; dimensionless QoIs use unit '1'; exact \
                    bitwise/boolean/count verdicts use unit 'bit' or '1' as declared \
                    per claim",
            seeds: "Philox 4x32-10 counter streams keyed 'i12/<fixture-id>/<case-index>'; \
                    development cases use seed indices 0..=4095; held-out cases use \
                    65536..=69631 (disjoint by construction, split frozen here)",
            budgets: "smoke tier <= 60 s on one host; core tier <= 30 min; max tier <= 8 h \
                      on a quiet perf host; <= 16 GiB memory per lane; event-scan \
                      subdivision and Zeno budgets are declared per obligation and \
                      exhaustion is a typed set-valued outcome, never a silent stop; \
                      accuracy budgets are the per-claim tolerance fields",
            versions: "fs-vmanifest schema v1; manifest version 1; toolchain pinned by \
                       rust-toolchain.toml; sibling pins by constellation.lock",
            capabilities: "no network; no FFI; deterministic mode mandatory for every G5 \
                           row; frontier/moonshot lanes stay behind feature flags; the \
                           true-flow (ValidatedStep) dependency is waivered until its \
                           bead lands",
        },
        claims: i12_claims(),
        fixtures: i12_fixtures(),
        obligations: i12_obligations(),
        waivers: i12_waivers(),
        amendment_rules: "Any change is a successor version through FrozenManifest::amend; \
                          the amendment record names every invalidated claim/obligation \
                          descendant; an amendment after campaign start invalidates the \
                          affected evidence, which must be re-earned; there is no in-place \
                          edit path in the type system",
    }
}

#[allow(clippy::too_many_lines)]
fn i12_claims() -> Vec<ClaimSpec> {
    vec![
        ClaimSpec {
            id: "i12-typed-automaton-ir",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "The hybrid automaton IR (modes, guards, resets, invariants, \
                        superdense-time semantics) round-trips its canonical encoding \
                        bit-stably and rejects every malformed-automaton class in the \
                        pinned corpus with a typed refusal naming the violated rule",
            hypotheses: &[
                "automata drawn from the pinned fixture families",
                "malformed classes are the enumerated ones: dangling mode/guard/reset \
                 references, unit-inconsistent guards, non-total reset maps, invariant/ \
                 guard contradictions declared detectable at admission",
            ],
            qoi: "roundtrip_and_refusal_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G0,
            oracle: OracleRoute {
                identity: "independent decoder plus the fixture generator's own validity \
                           certificate (emitted during generation, not by the compiler)",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "first implementation bead of the I12 IR leaf opens",
            kill: "any accepted malformed automaton or refused valid automaton returns \
                   the IR to design review with the disagreeing fixture as receipt",
            fallback: "explicit user-supplied automaton normal form with per-field \
                       validation",
            no_claim: "no claim about semantic equivalence of distinct encodings beyond \
                       the canonical form; no reachability claims at admission",
        },
        ClaimSpec {
            id: "i12-per-mode-structural-consistency",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Every mode's continuous dynamics passes the I02-class structural \
                        consistency check (matching completeness, index classification, \
                        consistent-initialization solvability) before the mode is \
                        admitted into the runtime automaton",
            hypotheses: &[
                "per-mode dynamics drawn from the pinned families including the \
                 rank-changing DAE fixture",
                "structural checks are per mode; cross-mode coupling enters only through \
                 resets",
            ],
            qoi: "structural_admission_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G0,
            oracle: OracleRoute {
                identity: "I02 certificate checkers (augmenting-path, DM partition, BLT \
                           topological validation) reused as independent validators of \
                           the per-mode admission verdict",
                independent: true,
                tcb_overlap: "shares the I02 checker library, which is independent of \
                              the I12 compiler under test",
            },
            activation: "IR leaf green at smoke tier",
            kill: "an admitted structurally deficient mode (checker refutation) blocks \
                   the admission leaf",
            fallback: "modes with failed checks are refused with the I02 diagnostic; no \
                       degraded admission",
            no_claim: "no numerical solvability claim; structure only, per the I02 \
                       boundary",
        },
        ClaimSpec {
            id: "i12-guard-root-accounting",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "For the declared isolable guard class (Taylor-pair guards with \
                        certified derivative bands over prescribed analytic segments), \
                        event scanning classifies every subinterval as root-free, \
                        certified-unique-crossing, or an explicit possible-event window, \
                        and the true root count lies in [confirmed, confirmed+possible]",
            hypotheses: &[
                "guards belong to the declared Taylor-pair class (guard and true- \
                 derivative bands describing one real function)",
                "segments are prescribed analytic motor/DAE trajectories, not simulated \
                 true flows (the true-flow rung is the waivered ValidatedStep lane)",
            ],
            qoi: "accounting_completeness_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "100x-resolution dense sign-scan falsifier plus analytic \
                           root-count fixtures with closed-form event times",
                independent: true,
                tcb_overlap: "none (the falsifier evaluates guards pointwise through \
                              the fixture's closed form, not the scanner)",
            },
            activation: "IR and structural leaves green at smoke tier",
            kill: "any dense-scan sign change outside every certified/possible window \
                   is Sev-0 and halts the event lane",
            fallback: "the classical dense-output lane as an explicitly Estimated \
                       baseline",
            no_claim: "no completeness claim for black-box guards or simulated flows; \
                       roots at segment junctions may surface as possible windows",
        },
        ClaimSpec {
            id: "i12-deterministic-superdense-runtime",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Superdense-time execution (mode sequence, event indices, reset \
                        applications, mode-local step outputs) replays bit-identically \
                        across runs, worker counts, and checkpoint/resume splits on the \
                        same ISA in deterministic mode",
            hypotheses: &[
                "deterministic mode; same ISA; the declared tie-break order for \
                 simultaneous admitted events",
            ],
            qoi: "bitwise_replay_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G5,
            oracle: OracleRoute {
                identity: "bitwise transcript comparison harness over independently \
                           serialized run transcripts",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "event-accounting leaf green at smoke tier",
            kill: "any cross-worker or split-run divergence is a release blocker",
            fallback: "none: determinism is a house contract, not a feature",
            no_claim: "no cross-ISA bitwise claim in this campaign; cross-ISA goes \
                       through the golden-divergence classification lane",
        },
        ClaimSpec {
            id: "i12-event-derivative-protocol",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Event-time and saltation derivative records for regular isolated \
                        transitions agree with central finite differences of the whole \
                        hybrid map on the pinned smooth-transition fixtures",
            hypotheses: &[
                "transitions are regular and isolated (transversal guard crossing, \
                 locally unique event time)",
                "the FD stencil width is inside the declared regularity neighborhood",
            ],
            qoi: "relative_derivative_error",
            unit: "1",
            tolerance: ToleranceSemantics::Relative { rtol: 1e-6 },
            evidence_tier: GauntletTier::G1,
            oracle: OracleRoute {
                identity: "central finite differences of the composed hybrid map \
                           (independent of the derivative-record construction)",
                independent: true,
                tcb_overlap: "shares the primal simulator; derivative paths disjoint",
            },
            activation: "superdense runtime leaf green at smoke tier",
            kill: "gate failure on any smooth fixture kills the derivative-record lane \
                   until rederived from the actual discrete residual",
            fallback: "FD-only sensitivities labeled Estimated",
            no_claim: "no derivative claim at simultaneous, grazing, or active-set- \
                       change events: those produce set-valued records or explicit \
                       no-claim outcomes",
        },
        ClaimSpec {
            id: "i12-set-valued-simultaneous-orders",
            ambition: Ambition::Frontier,
            polarity: ClaimPolarity::Affirmative,
            statement: "Simultaneous or window-overlapping admitted events produce the \
                        complete enumeration of admissible orders (or an explicitly \
                        unordered set-valued outcome above the enumeration cap); no \
                        order is silently picked",
            hypotheses: &[
                "events admitted by the accounting leaf with overlapping certified \
                 windows",
                "enumeration cap per the declared budget; above it the outcome is \
                 explicitly unordered",
            ],
            qoi: "admissible_order_set_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "symmetric analytic fixtures whose admissible order sets are \
                           known by construction (permutation envelope checker)",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "event-accounting leaf green at core tier",
            kill: "any silently ordered simultaneous pair kills the sequencing lane",
            fallback: "refuse-and-report: the run halts at the simultaneous window with \
                       a set-valued receipt",
            no_claim: "no probabilistic order selection; no claim about which order the \
                       physical system takes",
        },
        ClaimSpec {
            id: "i12-topology-changing-resets",
            ambition: Ambition::Frontier,
            polarity: ClaimPolarity::Affirmative,
            statement: "Reset maps that change contact/constraint topology re-admit the \
                        target mode through the full structural pipeline and carry a \
                        typed receipt naming conserved and non-conserved quantities \
                        across the reset",
            hypotheses: &[
                "resets drawn from the pinned topology-reset fixture family",
                "conservation declarations are per-fixture inputs, not inferred",
            ],
            qoi: "reset_readmission_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G2,
            oracle: OracleRoute {
                identity: "independent post-reset structural check plus fixture-declared \
                           conservation ledgers",
                independent: true,
                tcb_overlap: "shares the I02 checker library only",
            },
            activation: "structural-consistency leaf green at core tier",
            kill: "an un-readmitted topology change reaching execution kills the reset \
                   lane",
            fallback: "topology-changing resets refuse at admission (static automata \
                       only)",
            no_claim: "no impact-law physics claim; the mechanical correctness of a \
                       specific restitution model is a physics-campaign question",
        },
        ClaimSpec {
            id: "i12-bounded-horizon-mode-coverage",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Affirmative,
            statement: "Bounded-horizon mode-reachability coverage is complete for the \
                        declared bounded mode universe: every mode sequence not proven \
                        unreachable is either explored or reported as a localized \
                        Unknown window with its budget receipt",
            hypotheses: &[
                "bounded mode universe and horizon declared per fixture",
                "guard accounting from the [S] lane supplies the event windows",
            ],
            qoi: "coverage_completeness_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "exhaustive small-universe enumeration on fixtures small \
                           enough to enumerate exactly",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "behind the i12-moonshot feature flag; requires the ValidatedStep \
                         waiver retired and the frontier lanes green",
            kill: "any unreported unreachable-unproven sequence on an enumerable \
                   fixture kills the coverage lane",
            fallback: "coverage restricted to development-declared mode sequences, \
                       labeled as such",
            no_claim: "no unbounded-horizon claim; no claim when the mode universe \
                       declaration is violated at runtime",
        },
        ClaimSpec {
            id: "i12-certified-zeno-handling",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Affirmative,
            statement: "Zeno/accumulation candidates are detected within the declared \
                        event budget and produce set-valued continuation outcomes \
                        (accumulation-point enclosure plus post-Zeno state set or an \
                        explicit refusal), never a silent truncation of the event \
                        sequence",
            hypotheses: &[
                "fixtures with analytically known Zeno points (chatter family)",
                "event budgets per the explicits; exhaustion is typed",
            ],
            qoi: "zeno_outcome_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G4,
            oracle: OracleRoute {
                identity: "closed-form Zeno accumulation times of the chatter fixture \
                           family",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "behind the i12-moonshot feature flag; frontier lanes green",
            kill: "a silently truncated event sequence on a known-Zeno fixture is Sev-0 \
                   for this lane",
            fallback: "runs refusing at the Zeno budget with the accumulation enclosure \
                       as the terminal receipt",
            no_claim: "no claim about physically correct post-Zeno dynamics selection; \
                       the outcome is set-valued by design",
        },
        ClaimSpec {
            id: "i12-grazing-false-certificate-falsifier",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Refutation,
            statement: "Adversarial tangency families (grazing guards, even-order touch \
                        points, near-tangency continua) attempt to extract a certified \
                        unique-crossing verdict where Unknown is the correct answer; \
                        any success refutes the event lane's certificate discipline",
            hypotheses: &[
                "fixtures constructed so the analytic verdict (touch without crossing) \
                 is provable by hand",
                "the scanner runs at its production settings, not loosened for the test",
            ],
            qoi: "false_certificate_count",
            unit: "1",
            tolerance: ToleranceSemantics::Interval { lo: 0.0, hi: 0.0 },
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "hand-proved analytic tangency constructions (the fixture \
                           spec carries the proof sketch)",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "runs in every campaign tier from the first event-lane commit",
            kill: "this lane is never killed; it is the standing Sev-0 tripwire",
            fallback: "none: a nonzero count is a release blocker, not a degradable \
                       result",
            no_claim: "absence of refutation is not a proof of completeness; the lane \
                       only falsifies",
        },
    ]
}

fn i12_fixtures() -> Vec<FixturePin> {
    vec![
        FixturePin {
            id: "i12-gear-backlash-family",
            source: FixtureSource::AuthoredSpec {
                spec: "i12 fixture v1: gear pair with backlash automaton; modes \
                       {drive+, lash, drive-}; guards on relative angle crossing \
                       +/- half-lash (rad); resets exchange contact constraint; \
                       parameters: lash in {1e-3, 1e-2} rad, inertia ratio in \
                       {0.5, 1, 2}, drive torque profiles {constant, sinusoid \
                       1 Hz}; horizon 10 s; superdense semantics; seeds per \
                       explicits key 'i12/gear-backlash/<k>'",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i12-gas-valve-0d-family",
            source: FixtureSource::AuthoredSpec {
                spec: "i12 fixture v1: 0-D gas control volume with poppet valve; \
                       modes {closed, open, choked}; guards on pressure ratio \
                       crossing crack/reseat and critical ratio (dimensionless); \
                       reset re-admits flow constraint; parameters: crack 1.2, \
                       reseat 1.1, critical 1.893 (diatomic), volumes {1e-4, \
                       1e-3} m^3; horizon 2 s; seeds key 'i12/gas-valve/<k>'",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i12-chatter-zeno-family",
            source: FixtureSource::AuthoredSpec {
                spec: "i12 fixture v1: bouncing-contact chatter with restitution \
                       e in {0.5, 0.8, 0.9}; closed-form Zeno accumulation time \
                       t_z = t_0 + (2 v_0 / g) * e/(1-e); modes {flight, \
                       contact}; guard height crossing zero (m); horizon beyond \
                       t_z by construction; seeds key 'i12/chatter/<k>'",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i12-simultaneous-symmetric-family",
            source: FixtureSource::AuthoredSpec {
                spec: "i12 fixture v1: symmetric double-guard constructions whose \
                       two guards cross at the same analytic instant (mirror \
                       rotors, equal-phase valves); admissible order set = both \
                       orders by symmetry proof in-spec; window-overlap variants \
                       with analytically separated near-coincidences; seeds key \
                       'i12/simultaneous/<k>'",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i12-rank-changing-dae-family",
            source: FixtureSource::AuthoredSpec {
                spec: "i12 fixture v1 (HELD-OUT): mode pair whose constraint \
                       Jacobian drops rank across the reset (pin joint engaging \
                       a slider); per-mode index classes {1, 2}; consistent \
                       initialization defined per mode; adjudication uses the \
                       I02 checker certificates; seeds key 'i12/rank-dae/<k>'",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "i12-topology-reset-family",
            source: FixtureSource::AuthoredSpec {
                spec: "i12 fixture v1 (HELD-OUT): contact make/break resets that \
                       change constraint topology (latch engage/release, clutch \
                       lock/slip); fixture declares conserved (momentum map \
                       components) and non-conserved (energy at impact) \
                       quantities per reset; seeds key 'i12/topology-reset/<k>'",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "i12-grazing-tangency-family",
            source: FixtureSource::AuthoredSpec {
                spec: "i12 fixture v1 (HELD-OUT): hand-proved tangency battery; \
                       guard families g(t) = 1 - cos(w t) - c at c = 0 (touch \
                       without crossing at t = 2 pi k / w; proof: g >= 0 with \
                       equality exactly at the lattice), even-order touch g = \
                       (t - t_0)^{2m} h(t) with h > 0, and near-tangency c = \
                       +/- 1e-12 controls; correct verdict is Unknown/possible- \
                       event at the touch points, certified crossing for the \
                       controls; seeds key 'i12/grazing/<k>'",
            },
            partition: Partition::HeldOut,
        },
    ]
}

#[allow(clippy::too_many_lines)]
fn i12_obligations() -> Vec<ObligationRow> {
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
            leaf: "i12-ir-admission",
            claims_covered: &[
                "i12-typed-automaton-ir",
                "i12-per-mode-structural-consistency",
            ],
            unit_cases: UNIT_CASES,
            g0: "generators: automaton families incl. enumerated malformed classes; \
                 validity predicate: generator certificate agreement; laws: canonical \
                 round-trip bit-stability, refusal totality over malformed classes, \
                 admission idempotence; shrinker: mode/edge removal preserving the \
                 violated rule; replay seeds per explicits",
            decks: &["i12-gear-backlash-family", "i12-gas-valve-0d-family"],
            g3_relations: &[
                "mode relabeling invariance",
                "guard unit-rescaling covariance",
            ],
            g4_schedule: "cancel admission between IR decode, per-mode structural check, \
                          and automaton assembly; drained cancellation leaves no \
                          partially admitted automaton",
            g5_matrix: "threads {1,2,7} x mode {deterministic} x same-ISA replay",
            entry_point: "scripts/e2e/leapfrog/i12_admission.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i12-admission slice)",
            obs_events: &[
                "automaton.admitted",
                "automaton.refused",
                "admission.cancelled",
            ],
            replay_command: "scripts/e2e/leapfrog/i12_admission.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i12-event-accounting",
            claims_covered: &[
                "i12-guard-root-accounting",
                "i12-grazing-false-certificate-falsifier",
            ],
            unit_cases: UNIT_CASES,
            g0: "generators: Taylor-pair guard families over prescribed segments; laws: \
                 every leaf classified, count interval validity, certified windows \
                 contain analytic roots, refinement never weakens a certificate; \
                 shrinker: segment bisection; replay seeds per explicits",
            decks: &[
                "i12-gear-backlash-family",
                "i12-simultaneous-symmetric-family",
                "i12-grazing-tangency-family",
            ],
            g3_relations: &[
                "time-translation covariance of event windows",
                "guard sign-flip maps rising to falling crossings",
                "dense-scan falsifier containment at 100x resolution",
            ],
            g4_schedule: "exhaust subdivision and Zeno budgets deliberately; every \
                          exhaustion surfaces as typed possible-event windows; cancel \
                          mid-scan at interval-pop boundaries with drained state",
            g5_matrix: "threads {1,2,7} x mode {deterministic} x same-ISA replay of \
                        windows, counts, and receipts",
            entry_point: "scripts/e2e/leapfrog/i12_events.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i12-events slice)",
            obs_events: &[
                "event.certified",
                "event.possible",
                "scan.budget",
                "scan.cancelled",
            ],
            replay_command: "scripts/e2e/leapfrog/i12_events.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i12-superdense-runtime",
            claims_covered: &["i12-deterministic-superdense-runtime"],
            unit_cases: UNIT_CASES,
            g0: "generators: admitted automata with event schedules from the accounting \
                 leaf; laws: transcript bit-equality across runs and splits, superdense \
                 index monotonicity, reset atomicity; replay seeds per explicits",
            decks: &["i12-gear-backlash-family", "i12-gas-valve-0d-family"],
            g3_relations: &["checkpoint/resume split-run equality at every event index"],
            g4_schedule: "inject cancellation at mode-step, event-resolution, and reset \
                          boundaries; verify request-drain-finalize with no partial \
                          transcript acceptance",
            g5_matrix: "threads {1,2,7} x mode {deterministic} x same-ISA transcript \
                        bit-replay",
            entry_point: "scripts/e2e/leapfrog/i12_runtime.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i12-runtime slice)",
            obs_events: &[
                "mode.entered",
                "event.resolved",
                "reset.applied",
                "run.cancelled",
            ],
            replay_command: "scripts/e2e/leapfrog/i12_runtime.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i12-derivative-protocols",
            claims_covered: &["i12-event-derivative-protocol"],
            unit_cases: UNIT_CASES,
            g0: "generators: regular isolated transitions from the smooth subfamilies; \
                 laws: saltation record vs central FD within rtol, record refusal at \
                 declared non-regular events; replay seeds per explicits",
            decks: &["i12-gear-backlash-family", "i12-gas-valve-0d-family"],
            g3_relations: &["parameter-shift covariance of event-time derivatives"],
            g4_schedule: "cancel between primal event resolution and derivative-record \
                          assembly; a drained cancellation yields no partial record",
            g5_matrix: "threads {1,2,7} x mode {deterministic} x same-ISA replay of \
                        derivative records",
            entry_point: "scripts/e2e/leapfrog/i12_derivatives.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i12-derivatives slice)",
            obs_events: &["saltation.recorded", "saltation.refused"],
            replay_command: "scripts/e2e/leapfrog/i12_derivatives.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i12-nonsmooth-outcomes",
            claims_covered: &[
                "i12-set-valued-simultaneous-orders",
                "i12-topology-changing-resets",
            ],
            unit_cases: UNIT_CASES,
            g0: "generators: symmetric simultaneous families and topology-reset \
                 families; laws: admissible-order completeness on symmetric fixtures, \
                 cap-exceeded outcomes explicitly unordered, post-reset re-admission \
                 totality, conservation-ledger agreement; replay seeds per explicits",
            decks: &[
                "i12-simultaneous-symmetric-family",
                "i12-topology-reset-family",
                "i12-rank-changing-dae-family",
            ],
            g3_relations: &[
                "symmetry-group action permutes admissible orders",
                "reset conjugation by mode relabeling",
            ],
            g4_schedule: "cancel inside order enumeration and inside reset re-admission; \
                          drained state names the completed prefix only",
            g5_matrix: "threads {1,2,7} x mode {deterministic} x same-ISA replay of \
                        order sets and reset receipts",
            entry_point: "scripts/e2e/leapfrog/i12_nonsmooth.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i12-nonsmooth slice)",
            obs_events: &["orders.enumerated", "orders.capped", "reset.readmitted"],
            replay_command: "scripts/e2e/leapfrog/i12_nonsmooth.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i12-moonshot-coverage-zeno",
            claims_covered: &[
                "i12-bounded-horizon-mode-coverage",
                "i12-certified-zeno-handling",
            ],
            unit_cases: UNIT_CASES,
            g0: "generators: enumerable small mode universes and the chatter family; \
                 laws: coverage completeness vs exhaustive enumeration, Zeno enclosure \
                 contains the closed-form accumulation time, budget exhaustion typed; \
                 replay seeds per explicits",
            decks: &["i12-chatter-zeno-family", "i12-gear-backlash-family"],
            g3_relations: &["horizon-extension monotonicity of the covered set"],
            g4_schedule: "the core of this lane IS G4: Zeno budget exhaustion, timeout, \
                          allocation-failure injection at accumulation windows, and \
                          cancellation inside coverage search, each with drained typed \
                          outcomes",
            g5_matrix: "threads {1,2,7} x mode {deterministic} x same-ISA replay of \
                        coverage sets and Zeno receipts",
            entry_point: "scripts/e2e/leapfrog/i12_moonshot.sh",
            tier: CampaignTier::Max,
            dsr_lane: "dsr quality --tool frankensim (i12-moonshot slice)",
            obs_events: &["coverage.window", "zeno.enclosed", "budget.receipt"],
            replay_command: "scripts/e2e/leapfrog/i12_moonshot.sh --replay <artifact-id>",
        },
    ]
}

fn i12_waivers() -> Vec<Waiver> {
    vec![Waiver {
        subject: "i12-moonshot-coverage-zeno",
        reason: "the true-flow enclosure dependency (ValidatedStep tubes over simulated \
                 dynamics) is unimplemented; only prescribed-path guard accounting \
                 exists today, so the moonshot coverage/Zeno lane cannot activate \
                 against simulated flows",
        owner: "initiative-12 lane owner",
        predicate: "bead frankensim-ext-time-validated-step-ow2o closes with green \
                    central proof and its event-location rung is wired into the I12 \
                    runtime",
        expiry: "first I12 campaign review after ow2o closes; re-justify or retire at \
                 every manifest amendment",
        promotion_effect: "the [M] claims stay Unknown and cannot close; [S]/[F] \
                           promotion is unaffected because their obligations run on \
                           prescribed-path fixtures only",
    }]
}
