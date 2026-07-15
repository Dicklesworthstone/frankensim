//! The I02 (machine causalization and structural-index compiler)
//! VerificationManifest draft (bead
//! frankensim-leapfrog-2026-program-i94v.1.2.8.1).
//!
//! Baseline lattice ([S]): equation-variable incidence, deterministic
//! matching/DM/SCC/BLT, scoped index reduction, consistent
//! initialization, executable block plans with repair witnesses. Maximal
//! lattice ([F]/[M]): minimality and presentation-invariance of causal
//! witnesses, hybrid-mode-wide structural completeness, certified
//! hidden-constraint discovery, and globally optimal tearing (falsifier
//! lane). A weaker receipt closes its own element and is never relabeled
//! as the stronger theorem.

use crate::{
    Ambition, CampaignTier, ClaimPolarity, ClaimSpec, FiveExplicits, FixturePin, FixtureSource,
    GauntletTier, ManifestDraft, ObligationRow, OracleRoute, Partition, ToleranceSemantics, Waiver,
};

/// Build the I02 draft. Consumers freeze it themselves; the conformance
/// battery proves it freezes.
#[must_use]
pub fn i02_draft() -> ManifestDraft {
    ManifestDraft {
        initiative: "I02",
        title: "Machine causalization and structural-index compiler gate: incidence to \
                 executable block plans with reduction, initialization, and repair \
                 obligations",
        version: 1,
        explicits: FiveExplicits {
            units: "SI base units throughout; dimensionless QoIs use unit '1'; exact \
                    bitwise/boolean/count verdicts use unit 'bit' or '1' as declared \
                    per claim",
            seeds: "Philox 4x32-10 counter streams keyed 'i02/<fixture-id>/<case-index>'; \
                    development cases use seed indices 0..=4095; held-out cases use \
                    65536..=69631 (disjoint by construction, split frozen here)",
            budgets: "smoke tier <= 60 s on one host; core tier <= 30 min; max tier <= 8 h \
                      on a quiet perf host; <= 16 GiB memory per lane; accuracy budgets \
                      are the per-claim tolerance fields",
            versions: "fs-vmanifest schema v1; manifest version 1; toolchain pinned by \
                       rust-toolchain.toml; sibling pins by constellation.lock",
            capabilities: "no network; no FFI; deterministic mode mandatory for every G5 \
                           row; frontier/moonshot lanes stay behind feature flags",
        },
        claims: i02_claims(),
        fixtures: i02_fixtures(),
        obligations: i02_obligations(),
        waivers: i02_waivers(),
        amendment_rules: "Any change is a successor version through FrozenManifest::amend; \
                          the amendment record names every invalidated claim/obligation \
                          descendant; an amendment after campaign start invalidates the \
                          affected evidence, which must be re-earned; there is no in-place \
                          edit path in the type system",
    }
}

#[allow(clippy::too_many_lines)]
fn i02_claims() -> Vec<ClaimSpec> {
    vec![
        ClaimSpec {
            id: "i02-equation-variable-incidence",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "The extracted equation-variable incidence structure matches an \
                        independent extraction from the typed source on every pinned \
                        model, including differentiated-variable and mode-conditional \
                        incidences",
            hypotheses: &[
                "models drawn from the pinned toy-machine and DAE fixture families",
                "incidence is structural (nonzero pattern), not numerical",
            ],
            qoi: "incidence_agreement_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G0,
            oracle: OracleRoute {
                identity: "independent walker over the fixture generator's own incidence \
                           certificate (emitted during generation, not by the compiler)",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "first implementation bead of the I02 incidence leaf opens",
            kill: "any incidence disagreement the certificate cannot explain returns the \
                   extractor to design review with the disagreeing model as receipt",
            fallback: "explicit user-declared incidence annotations",
            no_claim: "no claim about numerical nonzero cancellation; structure only",
        },
        ClaimSpec {
            id: "i02-deterministic-matching-dm-scc-blt",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Maximum matching, Dulmage-Mendelsohn decomposition, SCC \
                        condensation, and BLT ordering are deterministic across runs and \
                        worker counts, and each is verified by its independent checker: \
                        no augmenting path exists, DM blocks partition correctly, and the \
                        BLT order is a valid topological order of the condensation",
            hypotheses: &[
                "graphs from the pinned fixture families including deficient graphs",
                "tie-breaking is the declared deterministic order (lexicographic by \
                 stable id)",
            ],
            qoi: "checker_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G0,
            oracle: OracleRoute {
                identity: "certificate checkers: augmenting-path search for matching \
                           maximality, partition-property checks for DM, topological \
                           validation for BLT (each independent of the construction \
                           algorithms)",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "incidence leaf green at smoke tier",
            kill: "any nondeterminism across worker counts is a release blocker; any \
                   checker refutation kills the affected construction algorithm",
            fallback: "reference (slow) constructions with the same checkers",
            no_claim: "no claim of optimal asymptotic complexity; correctness and \
                       determinism only",
        },
        ClaimSpec {
            id: "i02-scoped-index-reduction",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Structural index reduction lowers every pinned high-index DAE \
                        deck to structural index <= 1 with exactly the deck's expected \
                        per-equation differentiation counts, and the reduced system's \
                        incidence admits a complete matching",
            hypotheses: &[
                "decks from the pinned high-index DAE family with generator-computed \
                 expected differentiation counts",
                "reduction is structural (Pantelides-class); numerical regularity of \
                 the reduced system is a separate claim",
            ],
            qoi: "differentiation_count_agreement",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G1,
            oracle: OracleRoute {
                identity: "the fixture generator's differentiation-count certificate \
                           plus an independent complete-matching check on the reduced \
                           incidence",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "matching/decomposition leaf green at smoke tier",
            kill: "a deck reduced to a system whose incidence still lacks a complete \
                   matching kills the reduction route for that structure class",
            fallback: "user-declared index annotations with runtime constraint \
                       monitoring",
            no_claim: "structural index only: no claim that the reduced system is \
                       numerically nonsingular along trajectories",
        },
        ClaimSpec {
            id: "i02-consistent-initialization",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Computed consistent initial values satisfy every original and \
                        hidden (differentiated) constraint of the pinned decks within \
                        the stated band",
            hypotheses: &[
                "decks from the pinned high-index DAE family with known consistent \
                 initial manifolds",
                "initialization uses the declared deterministic solve schedule",
            ],
            qoi: "max_constraint_residual",
            unit: "1",
            tolerance: ToleranceSemantics::AbsRel {
                atol: 1e-10,
                rtol: 1e-8,
            },
            evidence_tier: GauntletTier::G1,
            oracle: OracleRoute {
                identity: "direct residual evaluation of the ORIGINAL (pre-reduction) \
                           constraint set at the computed initial point, evaluated \
                           without the initialization solver",
                independent: true,
                tcb_overlap: "shares fs-math deterministic scalar kernels",
            },
            activation: "index-reduction leaf green at smoke tier",
            kill: "persistent hidden-constraint violation beyond band kills the \
                   initialization schedule for that deck class",
            fallback: "user-supplied initial values with residual reporting",
            no_claim: "no claim of initialization uniqueness; one consistent point \
                       within band, on the pinned decks only",
        },
        ClaimSpec {
            id: "i02-block-plans-repair-witnesses",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Executable block plans schedule every BLT block exactly once in \
                        dependency order, and structurally deficient inputs yield typed \
                        repair witnesses naming the exact DM fine-structure deficiency \
                        (under/overdetermined parts) rather than a generic failure",
            hypotheses: &[
                "well-posed graphs and engineered deficient graphs from the pinned \
                 under/overdetermined family with generator-recorded deficiencies",
            ],
            qoi: "plan_and_witness_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G0,
            oracle: OracleRoute {
                identity: "schedule validator (each block once, dependencies before \
                           dependents) plus the generator's recorded deficiency \
                           certificate for the witness side",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "matching/decomposition leaf green at smoke tier",
            kill: "a deficient input that produces a plan instead of a witness (or a \
                   witness naming the wrong part) kills silent-repair behavior",
            fallback: "refuse all deficient inputs with the raw DM decomposition \
                       attached",
            no_claim: "witnesses name structural deficiencies; no claim that suggested \
                       repairs are physically meaningful",
        },
        ClaimSpec {
            id: "i02-causal-witness-minimality-invariance",
            ambition: Ambition::Frontier,
            polarity: ClaimPolarity::Affirmative,
            statement: "Causal witnesses are minimal (no proper sub-witness certifies \
                        the same causalization) and invariant across the pinned \
                        adversarial equivalent presentations of the same model",
            hypotheses: &[
                "presentation pairs carry the fixture generator's equivalence \
                 certificate",
                "minimality is with respect to the declared witness order frozen with \
                 the fixture family",
            ],
            qoi: "minimality_and_invariance_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "sub-witness enumerator (bounded by the frozen witness order) \
                           plus the generator's equivalence certificate",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "all five baseline claims green at core tier",
            kill: "a certified-equivalent pair with divergent witnesses kills \
                   presentation-dependent causalization",
            fallback: "canonical-form pre-normalization before causalization",
            no_claim: "invariance is claimed only for pairs carrying the generator's \
                       certificate, not for arbitrary user rewrites",
        },
        ClaimSpec {
            id: "i02-hybrid-mode-structural-completeness",
            ambition: Ambition::Frontier,
            polarity: ClaimPolarity::Affirmative,
            statement: "Structural analysis covers every reachable mode of the pinned \
                        hybrid mode families: each mode gets its own valid matching, \
                        BLT order, and index report, with mode-transition consistency \
                        recorded",
            hypotheses: &[
                "mode families from the pinned hybrid fixture set with \
                 generator-enumerated reachable modes",
                "mode reachability is the generator's declared enumeration, not \
                 discovered by the compiler",
            ],
            qoi: "uncovered_mode_count",
            unit: "1",
            tolerance: ToleranceSemantics::Interval { lo: 0.0, hi: 0.0 },
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "per-mode independent checker sweep (same checkers as the \
                           baseline claims) over the generator's mode enumeration",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "all five baseline claims green at core tier",
            kill: "any enumerated mode the analysis cannot causalize (without a typed \
                   witness) kills mode-generic analysis for that family",
            fallback: "per-mode explicit analysis with mode-set restrictions recorded",
            no_claim: "completeness is relative to the generator's mode enumeration; \
                       no claim about modes outside it",
        },
        ClaimSpec {
            id: "i02-certified-hidden-constraint-discovery",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Affirmative,
            statement: "Hidden-constraint discovery emits, for every pinned deck, a \
                        certificate deriving each hidden constraint from the original \
                        equations by the recorded differentiation/elimination steps, \
                        checkable without re-running discovery",
            hypotheses: &[
                "decks from the pinned high-index families (development and held-out)",
                "certificate checking is symbolic step replay, not numerical agreement",
            ],
            qoi: "certificate_replay_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "independent certificate replayer executing the recorded \
                           derivation steps against the original equations",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "maximal lane opens only after baseline closes at core tier",
            kill: "a certificate that fails replay on any deck kills certified \
                   discovery (discovery demotes to heuristic with runtime monitoring)",
            fallback: "runtime constraint monitoring without discovery certificates",
            no_claim: "certificates prove derivability of the emitted constraints; no \
                       claim that the emitted set is exhaustive (exhaustiveness is the \
                       falsifier lane's job)",
        },
        ClaimSpec {
            id: "i02-globally-optimal-tearing",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Refutation,
            statement: "Falsifier lane: search for a tearing of any pinned block whose \
                        declared cost beats the compiler's certified-optimal tearing; \
                        the optimality claim survives only while this search comes up \
                        empty",
            hypotheses: &[
                "cost model is the frozen per-fixture tearing cost declared with the \
                 deck",
                "search budget is the max-tier budget in explicits; search space is \
                 the deck's declared finite candidate enumeration",
            ],
            qoi: "counterexample_count",
            unit: "1",
            tolerance: ToleranceSemantics::Interval { lo: 0.0, hi: 0.0 },
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "independent exhaustive/branch-and-bound tearing enumerator \
                           over the frozen candidate space and cost model",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "maximal lane opens only after baseline closes at core tier",
            kill: "first verified counterexample kills the global-optimality claim; \
                   the counterexample tearing is retained as the receipt",
            fallback: "optimality restated as within-factor bound via amendment, or \
                       tearing demotes to advisory",
            no_claim: "an empty search is Estimated evidence bounded by the search \
                       budget and the frozen candidate space; it is never an \
                       optimality proof",
        },
    ]
}

fn i02_fixtures() -> Vec<FixturePin> {
    vec![
        FixturePin {
            id: "i02-toy-machine-slice",
            source: FixtureSource::AuthoredSpec {
                spec: "GENERATOR: a toy electromechanical machine vertical slice \
(voltage source, motor with electrical and mechanical ports, shaft \
compliance, inertial load, speed sensor) parameterized over \
component counts 1..=8 per stage; each instance emitted WITH its \
incidence certificate, expected matching, DM blocks, BLT order, \
and causal witness recorded at generation time.\n\
SEEDS: Philox stream 'i02/toy-machine/<k>', k in 0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i02-high-index-dae-family",
            source: FixtureSource::AuthoredSpec {
                spec: "GENERATOR: high-index DAE decks built from mechanism and \
constraint archetypes (pendulum-class index-3 loops, slider \
chains, coupled constraint rings) with generator-computed \
expected structural index, per-equation differentiation counts, \
hidden-constraint sets, and a known consistent initial manifold \
per deck.\n\
SEEDS: Philox stream 'i02/dae/<k>', k in 0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i02-high-index-dae-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "Same generator and certificate discipline as \
i02-high-index-dae-family; SEEDS: Philox stream 'i02/dae/<k>', \
k in 65536..=69631; withheld until adjudication.",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "i02-under-overdetermined-graphs",
            source: FixtureSource::AuthoredSpec {
                spec: "GENERATOR: structurally deficient equation graphs with \
engineered under- and overdetermined parts (dangling variables, \
redundant equations, mixed deficiencies) and the exact expected \
Dulmage-Mendelsohn fine-structure diagnosis plus minimal repair \
witness recorded for each.\n\
SEEDS: Philox stream 'i02/deficient/<k>', k in 0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i02-hybrid-mode-families",
            source: FixtureSource::AuthoredSpec {
                spec: "GENERATOR: hybrid machine families with enumerated discrete \
modes (clutch engage/disengage, contact on/off, saturation \
regions), each mode's expected structure recorded, plus \
mode-transition consistency pairs.\n\
SEEDS: Philox stream 'i02/hybrid/<k>', k in 0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i02-adversarial-presentations",
            source: FixtureSource::AuthoredSpec {
                spec: "GENERATOR: adversarial equivalent-presentation pairs of the \
same model (equation reordering, variable renaming, algebraic \
rewrites within the declared rewrite system, alias insertion) \
each carrying an equivalence certificate constructed during \
generation.\n\
SEEDS: Philox stream 'i02/adversarial/<k>', k in 0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i02-adversarial-presentations-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "Same generator as i02-adversarial-presentations; SEEDS: Philox \
stream 'i02/adversarial/<k>', k in 65536..=69631; withheld until \
adjudication.",
            },
            partition: Partition::HeldOut,
        },
    ]
}

#[allow(clippy::too_many_lines)]
fn i02_obligations() -> Vec<ObligationRow> {
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
            leaf: "i02-incidence-extraction",
            claims_covered: &["i02-equation-variable-incidence"],
            unit_cases: UNIT_CASES,
            g0: "generators: toy-machine + DAE families; validity predicate: incidence \
                 certificate agreement; laws: extraction idempotence, mode-conditional \
                 incidence stability; shrinker: equation removal preserving certificate; \
                 replay seeds per explicits",
            decks: &["i02-toy-machine-slice", "i02-high-index-dae-family"],
            g3_relations: &[
                "equation reordering invariance",
                "variable renaming invariance",
            ],
            g4_schedule: "cancel extraction at every tile boundary class once; verify \
                          drain leaves no partial incidence artifact",
            g5_matrix: "threads {1,2,7} x mode {deterministic} x same-ISA replay",
            entry_point: "scripts/e2e/leapfrog/i02_incidence.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i02-incidence slice)",
            obs_events: &["incidence.extracted", "incidence.cancelled"],
            replay_command: "scripts/e2e/leapfrog/i02_incidence.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i02-matching-decomposition",
            claims_covered: &["i02-deterministic-matching-dm-scc-blt"],
            unit_cases: UNIT_CASES,
            g0: "generators: all graph families; laws: matching maximality (no \
                 augmenting path), DM partition properties, BLT topological validity, \
                 deterministic tie-breaks; shrinker: vertex removal; replay seeds per \
                 explicits",
            decks: &[
                "i02-toy-machine-slice",
                "i02-high-index-dae-family",
                "i02-under-overdetermined-graphs",
            ],
            g3_relations: &["graph isomorphism covariance of DM blocks"],
            g4_schedule: "cancel between matching, DM, SCC, and BLT stages; each drained \
                          stage leaves no partial decomposition",
            g5_matrix: "threads {1,2,7} x mode {deterministic} x same-ISA replay",
            entry_point: "scripts/e2e/leapfrog/i02_decomposition.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i02-decomposition slice)",
            obs_events: &["matching.done", "dm.blocks", "blt.order", "stage.cancelled"],
            replay_command: "scripts/e2e/leapfrog/i02_decomposition.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i02-index-reduction",
            claims_covered: &["i02-scoped-index-reduction"],
            unit_cases: UNIT_CASES,
            g0: "generators: high-index DAE family; laws: differentiation-count \
                 certificate agreement, complete matching on the reduced incidence, \
                 reduction idempotence on index-1 inputs; replay seeds per explicits",
            decks: &["i02-high-index-dae-family"],
            g3_relations: &["equation scaling invariance of structural reduction"],
            g4_schedule: "cancel mid-reduction at differentiation steps; resume must \
                          reproduce uninterrupted counts",
            g5_matrix: "threads {1,2,7} x mode {deterministic} x same-ISA replay",
            entry_point: "scripts/e2e/leapfrog/i02_reduction.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i02-reduction slice)",
            obs_events: &["reduction.step", "reduction.done", "reduction.cancelled"],
            replay_command: "scripts/e2e/leapfrog/i02_reduction.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i02-initialization",
            claims_covered: &["i02-consistent-initialization"],
            unit_cases: UNIT_CASES,
            g0: "generators: high-index DAE family with known consistent manifolds; \
                 laws: residuals of original + hidden constraints within band at the \
                 computed point; replay seeds per explicits",
            decks: &["i02-high-index-dae-family"],
            g3_relations: &["unit-rescaling covariance of initialization residuals"],
            g4_schedule: "cancel during the initialization solve; drained solve leaves \
                          no partial initial-value artifact",
            g5_matrix: "threads {1,2,7} x mode {deterministic} x same-ISA replay",
            entry_point: "scripts/e2e/leapfrog/i02_initialization.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i02-initialization slice)",
            obs_events: &["init.solved", "init.residuals", "init.cancelled"],
            replay_command: "scripts/e2e/leapfrog/i02_initialization.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i02-block-plans-repair",
            claims_covered: &["i02-block-plans-repair-witnesses"],
            unit_cases: UNIT_CASES,
            g0: "generators: well-posed + deficient graph families; laws: every block \
                 scheduled exactly once in dependency order; every engineered deficiency \
                 yields its exact recorded witness; replay seeds per explicits",
            decks: &["i02-toy-machine-slice", "i02-under-overdetermined-graphs"],
            g3_relations: &["block permutation invariance of the witness set"],
            g4_schedule: "cancel plan construction between blocks; drain leaves no \
                          partial plan",
            g5_matrix: "threads {1,2,7} x mode {deterministic} x same-ISA replay",
            entry_point: "scripts/e2e/leapfrog/i02_plans.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i02-plans slice)",
            obs_events: &["plan.block", "repair.witness", "plan.cancelled"],
            replay_command: "scripts/e2e/leapfrog/i02_plans.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i02-maximal-lanes",
            claims_covered: &[
                "i02-causal-witness-minimality-invariance",
                "i02-hybrid-mode-structural-completeness",
                "i02-certified-hidden-constraint-discovery",
                "i02-globally-optimal-tearing",
            ],
            unit_cases: UNIT_CASES,
            g0: "generators: hybrid + adversarial families and held-out partitions; \
                 laws: witness minimality under the frozen order, per-mode checker \
                 sweeps, certificate replay, tearing enumeration replay; replay seeds \
                 per explicits",
            decks: &[
                "i02-hybrid-mode-families",
                "i02-adversarial-presentations",
                "i02-adversarial-presentations-holdout",
                "i02-high-index-dae-holdout",
                "i02-industrial-machine-benchmark",
            ],
            g3_relations: &[
                "certified presentation-pair equivalence",
                "held-out vs development distribution parity checks",
            ],
            g4_schedule: "cancel each maximal lane once mid-campaign; a resumed lane \
                          must report identical verdicts",
            g5_matrix: "threads {1,2,7} x mode {deterministic} x same-ISA replay",
            entry_point: "scripts/e2e/leapfrog/i02_maximal.sh",
            tier: CampaignTier::Max,
            dsr_lane: "dsr quality --tool frankensim (i02-maximal slice)",
            obs_events: &[
                "witness.minimality",
                "hybrid.mode",
                "hidden.certificate",
                "tearing.counterexample",
            ],
            replay_command: "scripts/e2e/leapfrog/i02_maximal.sh --replay <artifact-id>",
        },
    ]
}

fn i02_waivers() -> Vec<Waiver> {
    vec![Waiver {
        subject: "i02-industrial-machine-benchmark",
        reason: "the external industrial machine benchmark deck (IFToMM-class multibody \
                 or vendor machine model) is not yet licensed or content-addressed; \
                 only its slot is preregistered",
        owner: "I02 implementation beads (frankensim-leapfrog-2026-program-i94v.1.2.x)",
        predicate: "deck bytes pinned with exact edition, license, digest, QoIs, and \
                    acceptance envelopes through the fs-vvreg registry discipline",
        expiry: "before the first Max-tier campaign run; review at each Phase-2 close \
                 burst",
        promotion_effect: "maximal claims may run without this deck but cannot cite it; \
                           any result depending on it stays Estimated until the pin \
                           lands",
    }]
}
