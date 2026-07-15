//! The I10 constitutive identifiability and optimal coupon design
//! VerificationManifest draft (bead
//! frankensim-leapfrog-2026-program-i94v.3.3.7.1).
//!
//! The baseline freezes non-confusable law/experiment/nuisance/discrepancy
//! schema identities; adjoint sensitivities published only inside independent
//! outward-rounded finite-difference enclosures; three-valued gauge-witnessed
//! structural identifiability; profile-likelihood/sloppiness practical
//! identifiability with preregistered coverage; conservative three-valued
//! manufacturable coupon/environment/sensor design admission; and blind
//! held-out design gain. Stronger lanes target discrepancy-aware robust
//! design, anytime-valid adaptive law discrimination, a machine-checked
//! identifiable-combination quotient theorem, and certified global design
//! optimality over a frozen finite grammar. A refutation lane attacks false
//! identifiability and design certificates. Preregistration is governance,
//! not proof of a law, parameter, design, theorem, or lab campaign.

use crate::{
    Ambition, CampaignTier, ClaimPolarity, ClaimSpec, FiveExplicits, FixturePin, FixtureSource,
    GauntletTier, ManifestDraft, ObligationRow, OracleRoute, Partition, ToleranceSemantics, Waiver,
};

const CAMPAIGN_POLICY_FIXTURE: &str = "i10-campaign-policy-v1";

/// Build the I10 draft. Consumers freeze it themselves; focused
/// conformance tests prove that the authored seed satisfies schema v2.
#[must_use]
pub fn i10_draft() -> ManifestDraft {
    ManifestDraft {
        initiative: "I10",
        title: "Constitutive identifiability and optimal coupon design gate: gauge-aware \
                verdicts, adjoint-checked sensitivities, manufacturable designs, blind held-out \
                gain, and certified design optimality",
        version: 1,
        explicits: FiveExplicits {
            units: "SI coherent values with exact source-unit records; stress Pa, temperature K, \
                    force N, displacement m, time s, rate 1/s; strain, coverage gaps and gain \
                    ratios dimensionless; parameter blocks bind per-component units, tensor order \
                    and frames; FIM entries carry inverse-product parameter units; \
                    exact/refusal/coverage verdicts use unit 'bit'; no unit-erased parameter \
                    vector enters analysis",
            seeds: "Philox 4x32-10 streams keyed by BLAKE3 domain \
                    org.frankensim.i10.fixture-stream.v1 over exact fixture, law-family, stage \
                    and case ids; development indices 0..=16383, core held-out 65536..=77823, \
                    maximal held-out 131072..=147455; retained truth draws and blind partition \
                    assignments are sealed; duplicate logical ids or cross-stage reads are \
                    IntegrityFailed",
            budgets: "smoke <= 90 s and 2 GiB; core <= 60 min and 16 GiB; maximal <= 12 h and \
                      32 GiB; <= 10^4 law families, 10^6 forward solves, 10^5 adjoint solves, \
                      10^7 FIM entries, 10^6 candidate designs, 10^5 posterior draws per fit and \
                      10^7 scenario evaluations; caps are preflighted and exhaustion is \
                      BudgetExhausted/Unknown, never identifiability or optimality",
            versions: "fs-vmanifest schema v2; ConstitutiveLawSchema v1; ExperimentDesign v1; \
                       SensitivityReceipt v1; IdentifiabilityVerdict v1; CouponDesign v1; \
                       GainFunctional v1; DiscriminationPolicy v1; rust-toolchain.toml and \
                       constellation.lock receipt-bound",
            capabilities: "safe Rust; no network or FFI in adjudication; deterministic mode; \
                           exact unit/frame/domain registry and constraint grammar; \
                           content-addressed truth retention and redaction reasons mandatory; no \
                           implicit unit conversion, silent gauge fixing, unverified symbolic \
                           simplification, favorable missing-scenario fill or automatic lab \
                           actuation; maximal theorem/design lanes feature-gated",
        },
        claims: i10_claims(),
        fixtures: i10_fixtures(),
        obligations: i10_obligations(),
        waivers: i10_waivers(),
        amendment_rules: "Every changed law/experiment/nuisance/discrepancy/parameter/constraint/ \
                          design/sensor/truth/threshold/seed/checker/grammar field creates the \
                          exact next manifest version through FrozenManifest::amend. Candidate \
                          and blind partitions freeze before evidence. The amendment record \
                          names affected predecessor authority; no result edits this version, \
                          and unchanged evidence is rebound only through authenticated \
                          byte-identical lineage.",
    }
}

#[allow(clippy::too_many_lines)]
fn i10_claims() -> Vec<ClaimSpec> {
    vec![
        ClaimSpec {
            id: "i10-law-experiment-schema-identity",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "ConstitutiveLawSchema, parameter-block, ExperimentDesign, loading-program, \
                        sensor/DAQ-channel, nuisance, discrepancy and campaign identities \
                        round-trip canonically; conflicting, aliased, unit-inconsistent or \
                        domain-violating declarations fail closed without merging distinct \
                        modeling authority",
            hypotheses: &[
                "every schema has a typed issuer namespace, canonical preimage, version and explicit \
                 parameter/unit/frame/domain blocks; display names never alias authority",
                "objectivity and material-symmetry metadata (isotropy class, invariant basis, reference \
                 frame) is declared, not inferred; thermodynamic admissibility flags are recorded claims, \
                 not ingestion-verified theorems",
            ],
            qoi: "schema_identity_roundtrip_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G0,
            oracle: OracleRoute {
                identity: "independent canonical decoder plus unit/frame/domain containment checker over \
                           a separately generated schema-mutation corpus",
                independent: true,
                tcb_overlap: "shares public schema bytes and the frozen unit registry, not production \
                              parsing, assembly or conflict resolution",
            },
            activation: "the ConstitutiveLawSchema v1 implementation leaf opens",
            kill: "one authority alias, silent unit coercion, accepted domain violation or lost \
                   symmetry/frame declaration blocks every downstream I10 claim",
            fallback: "quarantine the record as typed UnresolvedSchema; preserve raw bytes and permit \
                       no sensitivity, identifiability or design use of that branch",
            no_claim: "canonical schemas do not prove a law's physical validity, thermodynamic \
                       admissibility, fit quality or relevance to any real material",
        },
        ClaimSpec {
            id: "i10-sensitivity-adjoint-consistency",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "For every catalog law/experiment pair, adjoint parameter gradients of each \
                        declared QoI lie inside independently computed outward-rounded \
                        finite-difference enclosures under the frozen step/widening schedule, and \
                        forward/adjoint sensitivities agree exactly under frozen reparameterization \
                        chain-rule maps",
            hypotheses: &[
                "finite-difference step schedules, enclosure widening and near-kink exclusion windows \
                 freeze before any gradient is published",
                "nuisance and discrepancy parameters are carried through both routes; a gradient with a \
                 missing block is refused, not zero-filled",
            ],
            qoi: "adjoint_gradient_enclosure_containment_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G2,
            oracle: OracleRoute {
                identity: "independent directed-rounding finite-difference enclosure oracle and \
                           reparameterization chain-rule replayer",
                independent: true,
                tcb_overlap: "shares frozen law/experiment bytes only; forward, adjoint and enclosure \
                              implementations are disjoint",
            },
            activation: "schema identity passes for the referenced law/experiment pair",
            kill: "one adjoint component escaping its enclosure, a silently zero-filled nuisance block \
                   or a chain-rule mismatch quarantines the pair and its dependent analyses",
            fallback: "publish forward-only sensitivities marked Unverified with the enclosure gap; \
                       no identifiability or design authority consumes an unverified gradient",
            no_claim: "gradient agreement is not identifiability, convexity, well-posedness or design \
                       optimality, and holds only inside the declared parameter/loading domains",
        },
        ClaimSpec {
            id: "i10-structural-identifiability-verdicts",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Structural identifiability analysis over frozen law/experiment families \
                        returns three-valued verdicts — Identifiable, NonIdentifiable with an \
                        explicit replayable gauge witness (a parameter-space curve or finite \
                        reparameterization leaving every declared observable invariant), or \
                        Unknown — and never mints Identifiable from numerics alone",
            hypotheses: &[
                "the symbolic derivation method (differential-algebra/generating-series identity), \
                 observable set and input class freeze before any verdict; a full-rank numeric FIM alone \
                 never upgrades Unknown to Identifiable",
                "gauge-equivalent parameter sets are one physical authority; every NonIdentifiable \
                 witness replays through the forward model with observables invariant to the declared \
                 enclosure",
                "known confoundings are seeded facts: kinematic and isotropic hardening under monotone \
                 proportional loading, Prony relabeling and coincident-time-constant merging, and \
                 modulus-thickness products under single-axis membrane loading",
            ],
            qoi: "structural_identifiability_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "independent symbolic derivation checker plus forward-simulation witness \
                           replayer over families with sealed identifiability labels",
                independent: true,
                tcb_overlap: "shares canonical law/observable bytes, not production elimination, \
                              simplification or ranking code",
            },
            activation: "schema identity and adjoint consistency are green for the family",
            kill: "one false Identifiable on a known gauge family, an unreplayable witness or a \
                   numerics-only upgrade refutes the analyzer and blocks design-gain claims",
            fallback: "downgrade the family to Unknown with the failed derivation retained; design \
                       lanes treat Unknown parameters as unidentified and design for combinations only",
            no_claim: "structural verdicts hold for exact model structure and noise-free observables; \
                       they do not certify practical recoverability under noise or model discrepancy",
        },
        ClaimSpec {
            id: "i10-practical-identifiability-sloppiness",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Practical identifiability reports FIM spectra, sloppy directions and \
                        profile-likelihood intervals under the frozen noise/nuisance model; on \
                        synthetic truth, profile-likelihood interval coverage meets its \
                        preregistered band and eigen-direction ownership separates data-supported \
                        stiffness from prior/constraint-induced stiffness",
            hypotheses: &[
                "noise model, nuisance treatment, profile grid/termination and coverage scoring freeze \
                 before held-out access; acceptance uses exact rank arithmetic or a preregistered \
                 one-sided confidence bound, never a raw empirical proportion alone",
                "a small FIM eigenvalue is evidence of sloppiness, not a certificate of structural \
                 non-identifiability; the two claims carry distinct evidence and wording",
            ],
            qoi: "worst_profile_likelihood_coverage_gap",
            unit: "1",
            tolerance: ToleranceSemantics::Absolute { atol: 0.05 },
            evidence_tier: GauntletTier::G2,
            oracle: OracleRoute {
                identity: "exact linear/Gaussian microcase oracle with closed-form FIM/profile \
                           intervals plus a retained simulation-truth coverage adjudicator",
                independent: true,
                tcb_overlap: "shares frozen model/data bytes only; profiling, eigensolvers and coverage \
                              scoring are independently implemented",
            },
            activation: "structural verdicts exist for the family and sensitivities are enclosure-checked",
            kill: "coverage gap above 0.05, spectrum mismatch on closed-form microcases, unowned \
                   eigen-direction or profile truncation presented as a bounded interval refutes this lane",
            fallback: "widen to conservative interval reports marked PracticalUnknown; never narrow an \
                       interval by dropping a nuisance dimension",
            no_claim: "calibration on synthetic families is not truth of the noise model, coverage for \
                       real laboratories or authority outside the declared parameter/loading domains",
        },
        ClaimSpec {
            id: "i10-manufacturable-coupon-admission",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Coupon/environment/sensor design admission is a three-valued predicate over \
                        the frozen manufacturability, fixture, loading, environment and sensor \
                        constraint grammar: Feasible only when every hard constraint is \
                        conservatively satisfied, Infeasible on a witnessed violation, and Unknown \
                        for missing or incomparable declarations; feasibility precedes any \
                        information criterion",
            hypotheses: &[
                "constraints declare direction, inclusive/exclusive bounds, joint envelopes, units/frames \
                 and conservative evaluation semantics; plate-stock orientation, minimum radii, tolerance \
                 classes, grip/fixture envelopes, load-frame limits, chamber ranges and per-sensor \
                 range/resolution/bandwidth are all hard-typed",
                "no information-criterion value, however favorable, can compensate a hard-constraint \
                 violation or an Unknown feasibility input",
            ],
            qoi: "coupon_design_admission_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "independent exact-rational constraint evaluator over a manually adjudicated \
                           valid/invalid design corpus with feasibility witnesses",
                independent: true,
                tcb_overlap: "shares the frozen constraint grammar bytes only; evaluation and search \
                              code are disjoint",
            },
            activation: "the CouponDesign v1 grammar and constraint corpus are frozen",
            kill: "one false Feasible (unmanufacturable radius, over-range sensor, fixture collision, \
                   unreachable temperature-rate pair) blocks the design engine",
            fallback: "return Infeasible with a minimized witness or Unknown with ranked missing \
                       declarations; no design is emitted for fabrication",
            no_claim: "Feasible is not fabrication approval, lab scheduling, cost authority or a claim \
                       that the design is informative",
        },
        ClaimSpec {
            id: "i10-heldout-design-gain",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Frozen optimality-designed coupon campaigns, executed one-shot on blind \
                        held-out synthetic experiments, reduce posterior uncertainty of the \
                        declared identifiable combinations relative to frozen reference designs by \
                        at least the preregistered margin, with shortfall measured by the frozen \
                        gain functional",
            hypotheses: &[
                "candidate designs, reference designs, gain functional, identifiable-combination basis \
                 and scoring arithmetic freeze before held-out access",
                "gain is claimed only for identifiable combinations; unidentified directions report \
                 unchanged-by-construction rather than counting toward gain",
            ],
            qoi: "heldout_information_gain_shortfall",
            unit: "1",
            tolerance: ToleranceSemantics::Interval { lo: 0.0, hi: 0.05 },
            evidence_tier: GauntletTier::G2,
            oracle: OracleRoute {
                identity: "independent blind-campaign executor and gain scorer over sealed truth \
                           parameters and partitions",
                independent: true,
                tcb_overlap: "shares frozen design/campaign bytes only; execution, inference and \
                              scoring are independently implemented",
            },
            activation: "identifiability verdicts, coupon admission and the gain functional are frozen",
            kill: "shortfall above 0.05, gain claimed on an unidentified direction, reference-design \
                   replacement after results or partition leakage refutes the design-gain claim",
            fallback: "report the honest shortfall with per-combination decomposition; the reference \
                       design remains the recommended campaign",
            no_claim: "synthetic held-out gain is not real-laboratory gain, cost-effectiveness, \
                       robustness to unmodeled physics or superiority over undeclared designs",
        },
        ClaimSpec {
            id: "i10-discrepancy-aware-robust-design",
            ambition: Ambition::Frontier,
            polarity: ClaimPolarity::Affirmative,
            statement: "Under the declared model-discrepancy ambiguity set, returned robust designs \
                        keep the frozen worst-case information criterion within the preregistered \
                        regret band of the scenario-exhaustive oracle on bounded microcases, and \
                        report Unknown rather than a favorable criterion when any scenario \
                        evaluation is missing",
            hypotheses: &[
                "the ambiguity set, scenario grid, criterion, and regret normalization freeze before \
                 evaluation; widening the ambiguity set cannot improve the certified worst case",
                "every scenario evaluation is interval-valued or Unknown; missing physics never \
                 contributes a favorable point value",
            ],
            qoi: "worst_scenario_robust_criterion_regret",
            unit: "1",
            tolerance: ToleranceSemantics::Absolute { atol: 0.05 },
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "scenario-exhaustive oracle on bounded design/discrepancy microcases plus an \
                           independent regret checker",
                independent: true,
                tcb_overlap: "shares frozen scenario/design bytes, not production search or surrogate \
                              code",
            },
            activation: "feature i10-robust-design; baseline admission and gain lanes are green",
            kill: "regret above 0.05, a favorable criterion published with a missing scenario or an \
                   ambiguity-widening improvement refutes the robust designer",
            fallback: "fall back to the best verified nominal-design authority with discrepancy \
                       sensitivity reported descriptively",
            no_claim: "robustness is only over the declared ambiguity set and criterion; it is not \
                       validity of any discrepancy model or protection against unmodeled physics",
        },
        ClaimSpec {
            id: "i10-anytime-valid-adaptive-discrimination",
            ambition: Ambition::Frontier,
            polarity: ClaimPolarity::Affirmative,
            statement: "Sequential adaptive experiment selection between rival constitutive laws \
                        uses preregistered nonnegative e-processes valid under the declared \
                        filtration and optional stopping; under either null the false-selection \
                        probability bound holds, and Selected, Undecided, DataInvalid and \
                        ModelUnknown remain distinct",
            hypotheses: &[
                "logical arrival order, filtration, rival pair, betting functions, truncation, design- \
                 adaptation policy and multiplicity family freeze before monitoring",
                "certifier validation uses exhaustive finite microcases plus a preregistered \
                 time-uniform e-bound for simulation families; the upper bound, not the raw observed \
                 false-selection fraction, is the acceptance QoI",
            ],
            qoi: "anytime_valid_upper_bound_on_null_false_selection_probability",
            unit: "1",
            tolerance: ToleranceSemantics::Interval { lo: 0.0, hi: 0.05 },
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "independent exact/log-domain e-process implementation and either-null \
                           simulations with adversarial stopping and adaptive designs",
                independent: true,
                tcb_overlap: "shares only retained canonical streams and frozen policy bytes",
            },
            activation: "feature i10-adaptive-discrimination; rival laws and the filtration are frozen",
            kill: "null false-selection above 0.05, invalid supermartingale arithmetic, post-hoc \
                   family editing or a duplicated stream record blocks discrimination authority",
            fallback: "freeze automated selection and surface Undecided/DataInvalid with diagnostics; \
                       fixed-horizon summaries remain descriptive only",
            no_claim: "Undecided never proves model equivalence; Selected localizes statistical \
                       evidence under the frozen design policy, not physical truth of the chosen law",
        },
        ClaimSpec {
            id: "i10-gauge-quotient-completeness-theorem",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Affirmative,
            statement: "A machine-checked theorem constructs the identifiable-combination quotient \
                        for each frozen law/experiment family: every declared gauge action leaves \
                        the observables invariant, the quotient coordinates are complete (equal \
                        coordinates imply identical observables on the declared domain) and \
                        independent (no further collapse), with kernel replay retained",
            hypotheses: &[
                "a pre-proof successor freezes law/observable/group ASTs and canonical bytes, quotient \
                 semantics, total runtime-premise mapping and deterministic proof translation",
                "exact axiom allowlist and transitive closure, independent model decoder and nonvacuity \
                 witnesses are receipt-bound; domain restrictions are explicit, never silently widened",
            ],
            qoi: "kernel_checked_quotient_completeness_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G5,
            oracle: OracleRoute {
                identity: "pinned Lean kernel replay plus independent canonical-model translation and \
                           finite countermodel checker for incomplete or dependent coordinates",
                independent: true,
                tcb_overlap: "kernel/axioms and the frozen semantic model are declared TCB; no \
                              production quotient construction or proof search",
            },
            activation: "feature i10-quotient-theorem; a pre-proof successor freezes all machine \
                         artifacts and the baseline identifiability lanes are independently green",
            kill: "invariance-breaking counterexample, incomplete/dependent coordinates, missing \
                   premise, disallowed axiom, translation mismatch or kernel rejection leaves this \
                   claim Unknown/Refuted and blocks completeness wording",
            fallback: "retain three-valued structural verdicts with replayable witnesses; design lanes \
                       keep treating non-quotient directions as unidentified",
            no_claim: "version-1 prose mints no theorem authority; a quotient over declared gauges is \
                       not exhaustiveness over undeclared symmetries, real materials or discrepancy",
        },
        ClaimSpec {
            id: "i10-global-design-optimality",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Affirmative,
            statement: "Over a frozen finite coupon/environment/sensor design grammar and declared \
                        discrepancy ambiguity set, the engine returns a feasible design with \
                        checker-verified coincident global lower and feasible upper bounds on the \
                        frozen robust information criterion — a certified zero normalized \
                        optimality gap",
            hypotheses: &[
                "a pre-candidate successor freezes the finite grammar, validity and exclusion order, \
                 enumeration or verified relaxation bounds, rank/unrank/sharding, independent \
                 decoder/checker, preflight and completeness root",
                "all bound and criterion arithmetic and the frozen positive reporting scale are exact \
                 or outward-rounded so a zero normalized gap cannot hide rounding",
                "every feasibility and criterion evaluation is bound to the same identifiability and \
                 constraint evidence and returns interval/Unknown rather than a favorable point",
            ],
            qoi: "certified_normalized_design_optimality_gap",
            unit: "1",
            tolerance: ToleranceSemantics::Interval { lo: 0.0, hi: 0.0 },
            evidence_tier: GauntletTier::G5,
            oracle: OracleRoute {
                identity: "independent finite-grammar exhaustive oracle on microcases plus \
                           primal-feasibility and dual/lower-bound checker",
                independent: true,
                tcb_overlap: "shares frozen grammar/evaluation bytes, not production search, surrogate \
                              or bound implementation",
            },
            activation: "feature i10-global-design; the finite grammar and every consumed evidence \
                         authority are frozen",
            kill: "infeasible recommendation, missing design, underestimated uncertainty, inconsistent \
                   lower/upper bounds, any positive certified gap or an independently better \
                   admissible design refutes global optimality",
            fallback: "return the best verified feasible design with its honest gap/Unknown state and \
                       ranked evidence plan; no optimality wording",
            no_claim: "optimality is only over the finite encoded grammar and declared criterion; it \
                       cannot compensate a misspecified ambiguity set, constraint grammar, gain \
                       functional or invalid identifiability evidence",
        },
        ClaimSpec {
            id: "i10-false-identifiability-design-certificate-falsifier",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Refutation,
            statement: "Hidden schema, gauge, statistical, design and certificate mutants try to \
                        obtain a green identifiability/design/optimality certificate for a known \
                        invalid case; any accepted mutant refutes the corresponding certifier and \
                        is retained as a minimized counterexample",
            hypotheses: &[
                "mutants include gauge-equivalent reparameterizations presented as identifiable, \
                 monotone-loading hardening spoofs, Prony relabelings, unit-rescaled parameter aliases, \
                 noise-model swaps, profile-interval truncations, forged design feasibility, \
                 sensor-range spoofs, holdout label leaks, e-process stopping edits, dropped grammar \
                 shards and unsound lower bounds",
                "each mutant has an independently authored invariance, exact-rational, enclosure, \
                 countermodel or finite-grammar witness hidden until maximal adjudication",
            ],
            qoi: "false_certificate_count",
            unit: "1",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "independent mutation generator and witness checker with no production \
                           analyzer, profiler, designer or optimizer code",
                independent: true,
                tcb_overlap: "public schemas and the frozen unit registry only",
            },
            activation: "a maximal certifier candidate and all checker/version identities are sealed \
                         before mutant reveal",
            kill: "the intended lane succeeds when a false certificate is accepted: mark the certifier \
                   Refuted and block dependent promotion; zero accepted mutants is finite \
                   falsification evidence only",
            fallback: "disable the affected certifier and use the strongest independently surviving \
                       baseline with Unknown/no-design",
            no_claim: "a finite adversarial corpus cannot prove absence of gauge confusion, statistical \
                       invalidity, constraint forgery, grammar incompleteness or optimizer unsoundness",
        },
    ]
}

#[allow(clippy::too_many_lines)]
fn i10_fixtures() -> Vec<FixturePin> {
    vec![
        FixturePin {
            id: CAMPAIGN_POLICY_FIXTURE,
            source: FixtureSource::AuthoredSpec {
                spec: "I10_CAMPAIGN_POLICY_V1\n\
IDENTITY=ConstitutiveLawSchema,ParameterBlock,ExperimentDesign,LoadingProgram,SensorChannel,NuisanceBlock,DiscrepancyBlock,CouponDesign,Campaign,GainFunctional,Checker and Decision identities have distinct typed domains, canonical preimages, versions and validity; same display text never aliases authority\n\
LAW=every law binds parameter units/frames/domains, objectivity/material-symmetry metadata, admissible ranges and declared observables; thermodynamic admissibility flags are recorded claims, not ingestion-verified theorems; no implicit unit conversion or silent domain widening\n\
SENSITIVITY=adjoint gradients are published only inside independently computed outward-rounded finite-difference enclosures under frozen step/widening schedules; reparameterizations use frozen chain-rule maps; an escaped component quarantines the pair\n\
IDENTIFIABILITY=verdicts {Identifiable,NonIdentifiable,Unknown,IntegrityFailed} are three-valued and witness-bearing; gauge-equivalent parameter sets are one physical authority; a full-rank numeric FIM alone never mints Identifiable; kinematic and isotropic hardening are confounded under monotone proportional loading until a discriminating path is added\n\
SLOPPINESS=FIM spectra, sloppy directions and profile-likelihood intervals carry noise/nuisance-model identity and eigen-direction ownership; a small eigenvalue is evidence of sloppiness, not a certificate of structural non-identifiability; data-supported and prior/constraint-induced stiffness remain distinct\n\
DESIGN=coupon/environment/sensor feasibility {Feasible,Infeasible,Unknown,IntegrityFailed} is evaluated conservatively over the frozen constraint grammar; design feasibility precedes information optimality; no automatic lab actuation\n\
GAIN=held-out gain is scored one-shot on blind campaigns against frozen reference designs with the frozen gain functional over declared identifiable combinations; truth, labels and partitions stay sealed until submission; unidentified directions never count toward gain\n\
DISCRIMINATION=rival-law selection uses preregistered nonnegative e-processes under the declared filtration with frozen stopping/adaptation/multiplicity policy; Selected,Undecided,DataInvalid,ModelUnknown remain distinct; Undecided never proves model equivalence\n\
THEOREM_AUTHORITY=version 1 prose mints no quotient-completeness or global-design-optimality authority; pre-proof successors freeze law/observable/group ASTs, semantics, total runtime premises, deterministic AST-to-Lean translation, exact axiom allowlist {propext,Quot.sound,Classical.choice} and transitive closure; sorryAx, custom postulates and native/unsafe oracle shortcuts are IntegrityFailed; pre-candidate successors freeze finite grammar,validity,enumeration/bounds,rank-unrank/sharding and completeness root\n\
EVIDENCE_STATES=Execution{Succeeded,Failed,Cancelled,TimedOut,BudgetExhausted,InfrastructureFailed,IntegrityFailed}; Predicate{Accepted,Refuted,Unknown,NotEvaluated}; Claim{NoClaim,EvidencePartial,Reproduced,Refuted,Proved}; Design{Feasible,Infeasible,Unknown,NoDesign}; Promotion{Blocked,CoreEligible,MaxEligible}; axes never substitute\n\
HOLDOUT=laws,experiments,designs,truth parameters,labels,seeds,checker/policy versions and grammars freeze before access; each Core/Max heldout fixture has one named consumer stage and disjoint range; premature/cross-stage read,label leak,replacement,retry-after-result or post-result tuning is IntegrityFailed\n\
LIFECYCLE=request->drain->finalize; poll at schema,solve,witness,profile,design,campaign,scenario and proof/search tile boundaries; checkpoints bind manifest/candidate/evidence cutoff and completed logical ids; resume/fork cannot duplicate solves,draws,wealth or published designs\n\
LOGGING=bounded schema-versioned fs-obs JSONL with stable run/case/claim/fixture/leaf/law/experiment/design/sensor/checker/attempt ids, exact units,seeds,budgets,versions,capabilities and redaction reason; stdout is not evidence\n\
RETENTION=manifest,exact command,code/contract/schema/toolchain/registry hashes,law/experiment/design bytes,solver and adjoint receipts,witnesses,spectra,intervals,wealth paths,bounds,all terminal states,minimized counterexamples and replay-verifier result; promotion/refutation evidence durable; inaccessible/redacted/expired inputs constrain replay and never become silently verified\n\
FAILURE_BUNDLE=first semantic/sensitivity/identifiability/design/certificate divergence plus bounded context,schema identity,expected/actual unit/value/enclosure/verdict,causal predecessors,minimized witness,checker disagreement,terminal state and replay command; partial success cannot publish normal authority\n\
PROMOTION=baseline claims only I10.G4 after G2 reproduction/G3 falsification; maximal claims only I10.G7 after G5 reproduction/G6 red-team; missing,stale,waived,integrity-failed or inaccessible evidence cannot promote\n\
LEAF_REQUIREMENT=every obligation references this policy fixture,all nine unit classes,smoke/core/max tier,DSR lane,events,replay,G4 drain/checkpoint,G5 matrix,independent adjudication,performance envelope,accessibility/agent parity and consuming gate",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i10-constitutive-law-catalog",
            source: FixtureSource::AuthoredSpec {
                spec: "I10_FIXTURE_V1 canonical constitutive families with exact parameter/unit/frame \
                       blocks: linear isotropic and orthotropic elasticity; Neo-Hookean, Mooney-Rivlin \
                       and one-, two- and three-term Ogden hyperelasticity; Prony-series linear \
                       viscoelasticity with temperature shift; J2 elastoplasticity with Voce isotropic, \
                       Armstrong-Frederick kinematic and combined hardening; Norton-Bailey creep; \
                       scalar damage; thermoelastic coupling. Each family pins reference parameter \
                       sets, admissible domains, loading programs (monotone, cyclic, non-proportional, \
                       relaxation, creep, thermal), analytic or manufactured observable solutions and \
                       known gauge facts (Ogden term permutation, Prony relabeling, \
                       kinematic-versus-isotropic monotone confounding, membrane modulus-thickness \
                       product).",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i10-gauge-symmetry-counterexamples",
            source: FixtureSource::AuthoredSpec {
                spec: "I10_FIXTURE_V1 symmetry/gauge counterexample corpus: parameter-space curves and \
                       finite reparameterization groups leaving every declared observable invariant — \
                       Ogden exponent/modulus permutations, Prony term relabeling and merging at \
                       coincident time constants, modulus-thickness and modulus-area products under \
                       single-axis loading, hardening decompositions invariant under monotone \
                       proportional paths that split under cyclic or non-proportional paths, \
                       unit-rescaling gauges and near-gauge sloppy valleys with pinned FIM \
                       null/near-null spaces. Each entry carries an independent invariance witness and \
                       the discriminating experiment that breaks it when one exists.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i10-multisensor-coupon-corpus",
            source: FixtureSource::AuthoredSpec {
                spec: "I10_FIXTURE_V1 synthetic multi-sensor coupon campaigns: dogbone, notched, \
                       cruciform biaxial, Iosipescu shear, torsion and bending coupons instrumented \
                       with load cells, extensometers, strain gauges, full-field DIC displacement maps \
                       and thermocouples; heteroscedastic and spatially correlated noise, calibration \
                       drift, saturation/dropout, misalignment, grip slippage and machine-compliance \
                       nuisances. Every channel carries units, frames, covariance, sampling schedule \
                       and validity windows; raw and independently normalized streams are both pinned. \
                       Development seeds 4096..=8191.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i10-environment-loading-constraints",
            source: FixtureSource::AuthoredSpec {
                spec: "I10_FIXTURE_V1 manufacturability/environment/loading constraint grammar: \
                       plate-stock orientation cuts, minimum fillet radii, thickness and tolerance \
                       classes, machining and surface-state limits, grip and fixture envelopes, \
                       load-frame force/displacement/rate limits, chamber temperature/humidity ranges, \
                       sensor mounting exclusions and per-sensor range/resolution/bandwidth. Valid and \
                       invalid design corpus with independent feasibility witnesses: unmanufacturable \
                       radius, over-range sensor, fixture collision, unreachable temperature-rate pair \
                       and joint-envelope violations masked by favorable marginals.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i10-sloppy-model-families",
            source: FixtureSource::AuthoredSpec {
                spec: "I10_FIXTURE_V1 sloppy-model families with pinned spectra: Prony ladders with \
                       eigenvalue decades, multi-exponential relaxation chains, hyperelastic families \
                       near incompressibility and hardening families near monotone confounding. Exact \
                       linear/Gaussian microcases with closed-form FIM and profile intervals; \
                       simulation-truth cases with retained generating parameters, noise and nuisance \
                       draws; identical-observable pseudo-identifiability traps. Development seeds \
                       0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i10-quotient-theorem-card",
            source: FixtureSource::AuthoredSpec {
                spec: "I10_THEOREM_CARD_V1 TARGET. Intended successor proposition constructs, for each \
                       frozen law/experiment family, the identifiable-combination quotient: declared \
                       gauge group action, observable invariance, completeness (equal quotient \
                       coordinates imply identical observables on the declared domain) and \
                       independence (no further collapse), with explicit domain restrictions. \
                       REQUIRED MACHINE ARTIFACTS: law/observable/group ASTs and canonical bytes, \
                       total premise map, deterministic Lean translation/roundtrip, exact axiom \
                       allowlist {propext,Quot.sound,Classical.choice}, transitive closure, \
                       nonvacuity witnesses and retained kernel replay. Prose grants no theorem or \
                       completeness authority.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i10-design-grammar-card",
            source: FixtureSource::AuthoredSpec {
                spec: "I10_EXHAUSTIVENESS_CARD_V1 TARGET. Intended successor freezes the finite \
                       coupon/environment/sensor/loading design grammar, canonical design state, \
                       feasibility and interval information-criterion semantics, discrepancy ambiguity \
                       set, validity and exclusion order, exact enumeration or verified relaxation \
                       bounds, rank/unrank/sharding, independent decoder/checker, preflight and \
                       completeness root. Returned feasible upper bound, global lower bound and \
                       normalized gap are separately replayed. Prose grants no global-optimality or \
                       design authority.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i10-schema-adjoint-adversaries-core-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "I10_FIXTURE_V1 CORE HOLDOUT schema/adjoint adversaries, indices 65536..=69631: \
                       aliased law names with divergent parameter blocks, unit-inconsistent \
                       parameter/observable pairs, frame-swapped anisotropy, domain-violating \
                       reference sets, adjoint traps near loading-program kinks and activation \
                       boundaries, step-size-sensitive finite-difference cases and reparameterization \
                       chain-rule stressors. Each carries an independent decoder or enclosure witness. \
                       One I10.G3 consumer.",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "i10-identifiability-core-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "I10_FIXTURE_V1 CORE HOLDOUT identifiability families, indices 69632..=73727: \
                       blind law/experiment pairs with sealed Identifiable/NonIdentifiable/Unknown \
                       labels and hidden gauge witnesses — hardening splits revealed only under cyclic \
                       paths, Prony merges at coincident time constants, product gauges broken by an \
                       added sensor channel, sloppy-but-structurally-identifiable spectra and \
                       profile-coverage stressors. Labels and witnesses are inaccessible until \
                       one-shot submission. One I10.G3 consumer.",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "i10-blind-experiment-gain-core-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "I10_FIXTURE_V1 CORE HOLDOUT blind gain campaigns, indices 73728..=77823: \
                       sealed synthetic lab campaigns executing frozen candidate and reference designs \
                       on hidden truth parameters with realistic noise/nuisance draws; the gain \
                       functional over declared identifiable combinations is scored one-shot; truth, \
                       labels and partitions withheld until submission; includes feasibility traps \
                       whose favorable information values must not rescue infeasible designs. One \
                       I10.G3 consumer.",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "i10-discrimination-max-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "I10_FIXTURE_V1 MAX HOLDOUT robust-design/discrimination cases, indices \
                       131072..=135167: rival-law streams under either-null and alternative regimes \
                       with adversarial optional stopping, drift and nuisance contamination, adaptive \
                       design feedback, delayed/revised/duplicate records; plus scenario-deleted \
                       robust-design microcases whose scenario-exhaustive oracle certificates are \
                       withheld. One I10.G6 consumer.",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "i10-gauge-quotient-max-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "I10_FIXTURE_V1 MAX HOLDOUT quotient certifier stress family, indices \
                       135168..=139263: declared groups with invariance-breaking mutations, incomplete \
                       quotient coordinate sets, dependent coordinates, domain-restriction traps, \
                       translation and axiom adversaries and kernel-replay corruption cases; finite \
                       countermodels and invariance witnesses sealed until adjudication. One I10.G6 \
                       consumer.",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "i10-design-certifier-mutants-max-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "I10_FIXTURE_V1 MAX HOLDOUT certifier mutants, indices 139264..=147455: forged \
                       design feasibility, sensor-range spoof, dropped grammar shard, unsound \
                       relaxation/lower bound, favorable Unknown coercion, gauge-alias \
                       identifiability spoof, monotone-loading hardening spoof, unit-rescaled \
                       parameter alias, profile-interval truncation, noise-model swap, e-process \
                       stopping edit, holdout label leak and disallowed theorem axiom. Each mutant \
                       carries an independent witness sealed until the certifier candidate is frozen. \
                       One I10.G6 consumer.",
            },
            partition: Partition::HeldOut,
        },
    ]
}

#[allow(clippy::too_many_lines)]
fn i10_obligations() -> Vec<ObligationRow> {
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
            leaf: "i10-law-schema-adjoint",
            claims_covered: &[
                "i10-law-experiment-schema-identity",
                "i10-sensitivity-adjoint-consistency",
            ],
            unit_cases: UNIT_CASES,
            g0: "generators: law/experiment/nuisance/discrepancy schema mutations and catalog QoI/ \
                 adjoint pairs; predicates: canonical identity, unit/frame/domain containment, \
                 parameter-block typing, finite-difference enclosure containment and \
                 reparameterization chain rule; laws: roundtrip, no authority alias, unit-coercion \
                 refusal, enclosure directedness, gradient invariance under frozen \
                 reparameterization and refused zero-fill of missing blocks; shrink schema/loading \
                 program while preserving the independent witness",
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                "i10-constitutive-law-catalog",
                "i10-multisensor-coupon-corpus",
                "i10-schema-adjoint-adversaries-core-holdout",
            ],
            g3_relations: &[
                "canonical schema field order and id relabeling preserve identity and observables",
                "an exact compatible unit/frame transform of parameters and observables conjugates \
                 gradients exactly",
                "tightening a finite-difference step schedule under the frozen widening policy cannot \
                 evict a contained adjoint gradient",
                "a schema alias or domain-violating reference set fails before any sensitivity is \
                 published",
                "an observable-preserving loading-program reparameterization preserves QoI gradients",
            ],
            g4_schedule: "cancel/panic/timeout/corrupt input after schema decode, unit/frame check, \
                          forward solve, adjoint solve, enclosure assembly and receipt binding; poll \
                          each law/QoI pair; request->drain->finalize publishes matched \
                          forward/adjoint/enclosure artifacts or none; checkpoint completed pair ids \
                          and content roots; resume/fork cannot duplicate pairs; retain bounded \
                          FailureBundle and every terminal state",
            g5_matrix: "law/pair shards {1,2,7} x input orders {catalog,reverse,permuted} x \
                        deterministic mode on identical toolchain fingerprint; canonical schemas, \
                        gradients, enclosures, verdicts, event JSON and digests match bitwise",
            entry_point: "scripts/e2e/leapfrog/i10_law_schema.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i10-law-schema slice)",
            obs_events: &[
                "request.received",
                "schema.admitted",
                "schema.refused",
                "sensitivity.forward_completed",
                "sensitivity.adjoint_completed",
                "sensitivity.enclosure_checked",
                "cancel.requested",
                "drain.completed",
                "finalize.completed",
                "failure_bundle.retained",
                "adjudication.receipt",
            ],
            replay_command: "scripts/e2e/leapfrog/i10_law_schema.sh --manifest <manifest-id> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i10-identifiability-analysis",
            claims_covered: &[
                "i10-practical-identifiability-sloppiness",
                "i10-structural-identifiability-verdicts",
            ],
            unit_cases: UNIT_CASES,
            g0: "generators: law/experiment pairs with known gauge groups, sloppy spectra and \
                 analytic FIM/profile microcases; predicates: three-valued verdict typing, \
                 invariance-witness forward replay, FIM/profile arithmetic, eigen-direction \
                 ownership and coverage scoring; laws: witness replay invariance, no Identifiable \
                 from a full-rank numeric FIM alone, verdict monotonicity under experiment \
                 enrichment, profile-interval nesting and structural-versus-practical wording \
                 separation; shrink family while preserving the gauge witness",
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                "i10-constitutive-law-catalog",
                "i10-gauge-symmetry-counterexamples",
                "i10-identifiability-core-holdout",
                "i10-sloppy-model-families",
            ],
            g3_relations: &[
                "adding an observable channel or loading segment never turns Identifiable into \
                 NonIdentifiable for the same parameters",
                "deleting the discriminating cyclic segment restores the monotone hardening \
                 confounding verdict",
                "a pinned gauge curve replayed through the forward model leaves every declared \
                 observable invariant",
                "coverage scoring is invariant to holdout case order and external-id reassignment",
                "parameter relabeling and unit rescaling conjugate verdicts, spectra and witnesses \
                 exactly",
            ],
            g4_schedule: "cancel/panic/timeout at symbolic derivation, witness replay, FIM assembly, \
                          eigendecomposition, profile sweep and verdict binding; poll each \
                          pair/profile tile; request->drain->finalize publishes one terminal verdict \
                          per pair plus witnesses or nothing; checkpoint completed pair/profile ids \
                          and content roots; resume/fork equals one-shot; retain false-Identifiable \
                          and coverage FailureBundle",
            g5_matrix: "pair/profile shards {1,2,7} x workers {1,2,7} x case orders {forward,reverse} \
                        x deterministic mode on identical ISA; verdicts, witnesses, spectra, \
                        intervals, tie breaks, events and receipts match bitwise",
            entry_point: "scripts/e2e/leapfrog/i10_identifiability.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i10-identifiability slice)",
            obs_events: &[
                "request.received",
                "identifiability.identifiable",
                "identifiability.nonidentifiable",
                "identifiability.unknown",
                "identifiability.witness_replayed",
                "sloppiness.spectrum_reported",
                "cancel.requested",
                "drain.completed",
                "finalize.completed",
                "failure_bundle.retained",
                "adjudication.receipt",
            ],
            replay_command: "scripts/e2e/leapfrog/i10_identifiability.sh --manifest <manifest-id> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i10-coupon-design-gain",
            claims_covered: &[
                "i10-heldout-design-gain",
                "i10-manufacturable-coupon-admission",
            ],
            unit_cases: UNIT_CASES,
            g0: "generators: design-grammar states crossed with the constraint corpus and synthetic \
                 campaigns with retained truth; predicates: hard-constraint three-valued evaluation, \
                 conservative feasibility arithmetic, gain-functional binding to declared \
                 identifiable combinations and blind-partition integrity; laws: no false Feasible on \
                 the independent evaluator, Unknown monotonicity under declaration removal, \
                 reference-design pinning, gain-score campaign-order invariance and no gain on \
                 unidentified directions; shrink design/campaign while preserving the violation \
                 witness",
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                "i10-blind-experiment-gain-core-holdout",
                "i10-environment-loading-constraints",
                "i10-multisensor-coupon-corpus",
            ],
            g3_relations: &[
                "tightening any hard manufacturability/environment/sensor bound cannot convert \
                 Infeasible or Unknown to Feasible",
                "removing a constraint declaration moves feasibility only toward Unknown, never \
                 silently toward Feasible",
                "design and sensor-channel id permutation preserves admission verdicts and gain \
                 scores",
                "unit/frame transforms of constraint and design quantities preserve every verdict \
                 exactly",
                "replacing the frozen reference design after results is refused as an integrity \
                 failure",
            ],
            g4_schedule: "cancel/panic/timeout at grammar decode, constraint evaluation, campaign \
                          assembly, blind scoring and receipt binding; whole-campaign preflight \
                          precedes holdout access; poll each design/campaign tile; \
                          request->drain->finalize publishes one terminal \
                          Feasible/Infeasible/Unknown or scored-gain state; checkpoint \
                          frontier/visited/scored ids and content roots; resume equals one-shot; \
                          retain false-Feasible and gain-shortfall FailureBundle",
            g5_matrix: "design/campaign shards {1,2,7} x workers {1,2,7} x orders \
                        {lex,reverse,permuted} x deterministic mode on identical ISA; verdicts, \
                        witnesses, gain scores, tie breaks, events and receipts match bitwise",
            entry_point: "scripts/e2e/leapfrog/i10_coupon_design.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i10-coupon-design slice)",
            obs_events: &[
                "request.received",
                "design.admitted",
                "design.infeasible",
                "design.unknown",
                "gain.campaign_scored",
                "gain.shortfall",
                "cancel.requested",
                "drain.completed",
                "finalize.completed",
                "failure_bundle.retained",
                "adjudication.receipt",
            ],
            replay_command: "scripts/e2e/leapfrog/i10_coupon_design.sh --manifest <manifest-id> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i10-robust-discrimination-designer",
            claims_covered: &[
                "i10-anytime-valid-adaptive-discrimination",
                "i10-discrepancy-aware-robust-design",
            ],
            unit_cases: UNIT_CASES,
            g0: "generators: discrepancy ambiguity scenarios, rival-law streams under either null \
                 with stopping adversaries and bounded robust-design microcases with exhaustive \
                 oracles; predicates: scenario-complete evaluation, regret arithmetic, nonnegative \
                 wealth, filtration/bet identity and terminal-state separation; laws: ambiguity \
                 widening cannot improve the certified worst case, a missing scenario forces \
                 Unknown, null stopping preserves the false-selection bound, batching under the \
                 frozen filtration preserves wealth and Selected/Undecided/DataInvalid/ModelUnknown \
                 stay distinct; shrink scenario/stream while preserving the error",
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                "i10-constitutive-law-catalog",
                "i10-discrimination-max-holdout",
                "i10-environment-loading-constraints",
                "i10-sloppy-model-families",
            ],
            g3_relations: &[
                "widening the declared discrepancy ambiguity set cannot improve the certified \
                 worst-case criterion",
                "deleting a scenario evaluation forces Unknown rather than a favorable robust \
                 criterion",
                "adversarial optional stopping under either null preserves the false-selection bound \
                 <= 0.05",
                "rival-law relabeling conjugates selections and wealth paths exactly",
                "duplicate or revised stream records follow the frozen policy and cannot be counted \
                 twice favorably",
            ],
            g4_schedule: "whole-campaign preflight freezes scenario/stream budgets before max holdout \
                          access; cancel/panic/timeout at scenario evaluation, design update, wealth \
                          update, selection and receipt binding; request->drain->finalize persists \
                          exact e-process state and one terminal state; checkpoint filtration \
                          position, logical ids, wealth and scenario roots; resume cannot \
                          double-count or edit the stopping policy; BudgetExhausted stays Unknown; \
                          retain false-selection and regret FailureBundle",
            g5_matrix: "scenario/stream shards {1,2,7,31} x workers {1,2,7} x orders \
                        {lex,reverse,permuted} x deterministic mode on identical ISA; regrets, \
                        wealth paths, selections, terminal states, events and receipts match bitwise",
            entry_point: "scripts/e2e/leapfrog/i10_discrimination.sh",
            tier: CampaignTier::Max,
            dsr_lane: "dsr quality --tool frankensim (i10-discrimination isolated lane)",
            obs_events: &[
                "request.received",
                "robust.scenario_evaluated",
                "robust.design_updated",
                "robust.unknown",
                "discrimination.wealth_updated",
                "discrimination.selected",
                "discrimination.undecided",
                "cancel.requested",
                "drain.completed",
                "finalize.completed",
                "failure_bundle.retained",
                "adjudication.receipt",
            ],
            replay_command: "scripts/e2e/leapfrog/i10_discrimination.sh --manifest <manifest-id> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i10-quotient-theorem-certifier",
            claims_covered: &["i10-gauge-quotient-completeness-theorem"],
            unit_cases: UNIT_CASES,
            g0: "generators: law/observable/group ASTs, gauge actions, quotient candidates and \
                 kernel-replay adversaries; predicates: AST canonicality, group-action invariance, \
                 completeness and independence witnesses, translation roundtrip, axiom closure and \
                 kernel replay; laws: alpha-renaming preserves the quotient up to canonical \
                 rebinding, invariance-breaking mutation fails closed, incomplete or dependent \
                 coordinates are refused and a domain restriction never silently widens; shrink \
                 family/proof while preserving the counterexample",
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                "i10-gauge-quotient-max-holdout",
                "i10-gauge-symmetry-counterexamples",
                "i10-quotient-theorem-card",
            ],
            g3_relations: &[
                "law/parameter alpha-renaming preserves the quotient theorem up to canonical \
                 rebinding",
                "an invariance-breaking group mutation fails closed before kernel submission",
                "finite countermodels force Unknown for incomplete or dependent quotient coordinates",
                "a disallowed axiom, missing premise or translation mutation fails closed",
                "declared-family identifiability verdicts remain available when the maximal theorem \
                 is Unknown/Refuted",
            ],
            g4_schedule: "whole-campaign preflight precedes max holdout; cancel/panic/timeout/corrupt \
                          checkpoint at AST decode, invariance check, completeness search, \
                          translation, kernel replay and witness checking; poll per family/proof \
                          tile; request->drain->finalize drains children and preserves the terminal \
                          state; checkpoint frontier/derivation/proof/model roots; BudgetExhausted \
                          stays Unknown; retain countermodel and kernel FailureBundle",
            g5_matrix: "family/proof shards {1,2,7,31} x workers {1,2,7} x orders \
                        {lex,reverse,permuted} x deterministic mode on identical kernel/toolchain; \
                        canonical ASTs, invariance records, completeness verdicts, kernel results, \
                        witnesses, events and receipt roots match bitwise",
            entry_point: "scripts/e2e/leapfrog/i10_quotient.sh",
            tier: CampaignTier::Max,
            dsr_lane: "dsr quality --tool frankensim (i10-quotient isolated lane)",
            obs_events: &[
                "request.received",
                "quotient.family_admitted",
                "quotient.invariance_checked",
                "quotient.completeness_checked",
                "quotient.kernel_checked",
                "quotient.unknown",
                "cancel.requested",
                "drain.completed",
                "finalize.completed",
                "failure_bundle.retained",
                "adjudication.receipt",
            ],
            replay_command: "scripts/e2e/leapfrog/i10_quotient.sh --manifest <manifest-id> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i10-global-design-certifier",
            claims_covered: &[
                "i10-false-identifiability-design-certificate-falsifier",
                "i10-global-design-optimality",
            ],
            unit_cases: UNIT_CASES,
            g0: "generators: finite design grammars, interval feasibility/criterion evaluations, \
                 exact microcase optima, primal/dual bounds and hidden certificate mutants; \
                 predicates: grammar decode/validity/completeness, evidence binding, hard \
                 feasibility, interval criterion, lower <= feasible upper, normalized gap, \
                 rank/unrank/shard coverage and mutant witnesses; laws: design relabeling, \
                 exhaustive microcase agreement, adding a feasible design cannot worsen the \
                 optimum, ambiguity widening cannot improve the robust criterion and every known \
                 false certificate is refused; shrink grammar/mutant while preserving the witness",
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                "i10-design-certifier-mutants-max-holdout",
                "i10-design-grammar-card",
                "i10-environment-loading-constraints",
                "i10-instrumented-lab-campaign-pack",
            ],
            g3_relations: &[
                "design/scenario id permutation preserves the optimal design, bounds and gap up to \
                 canonical rebinding",
                "adding a feasible design can preserve or improve but cannot worsen the certified \
                 optimum",
                "widening ambiguity or evaluation intervals cannot improve a worst-case criterion \
                 without proof",
                "shard/frontier order preserves the completeness root, bounds and deterministic \
                 argmin tie break",
                "every accepted schema/gauge/statistical/design/certificate mutant refutes its \
                 certifier",
            ],
            g4_schedule: "whole-campaign preflight freezes grammar/evaluation budgets before mutant \
                          reveal; cancel/panic/timeout/corrupt checkpoint at design generation, \
                          feasibility, criterion evaluation, bound update, shard join, completeness \
                          and mutant adjudication; request->drain->finalize drains children and \
                          publishes only complete bound/terminal states; checkpoint \
                          frontier/visited/primal/dual roots; BudgetExhausted stays Unknown; retain \
                          the best feasible incumbent and every minimized false-certificate \
                          FailureBundle",
            g5_matrix: "design/scenario shards {1,2,7,31} x workers {1,2,7} x frontier orders \
                        {lex,reverse,permuted} x deterministic mode on identical toolchain; \
                        valid/excluded/evaluated counts, completeness roots, feasible incumbent, \
                        lower/upper/gap, optimal tie, mutant verdicts, events and receipts match \
                        bitwise",
            entry_point: "scripts/e2e/leapfrog/i10_global_design.sh",
            tier: CampaignTier::Max,
            dsr_lane: "dsr quality --tool frankensim (i10-global-design isolated lane)",
            obs_events: &[
                "request.received",
                "design.candidate_generated",
                "design.candidate_infeasible",
                "design.incumbent_updated",
                "design.bound_updated",
                "design.completeness_checked",
                "certifier.mutant_refused",
                "certifier.refuted",
                "cancel.requested",
                "drain.completed",
                "finalize.completed",
                "failure_bundle.retained",
                "adjudication.receipt",
            ],
            replay_command: "scripts/e2e/leapfrog/i10_global_design.sh --manifest <manifest-id> --replay <artifact-id>",
        },
    ]
}

fn i10_waivers() -> Vec<Waiver> {
    vec![Waiver {
        subject: "i10-instrumented-lab-campaign-pack",
        reason: "all committed I10 fixtures are synthetic coupon-shaped data; independently \
                 authorized physical lab campaigns — real coupons, DIC/load-frame raw streams, \
                 calibration certificates, machine/fixture records and blind lab truth — are not \
                 yet licensed, consent-cleared and content-addressed for replay",
        owner: "I10 laboratory governance, metrology and V&V owners",
        predicate: "a governed instrumented-lab pack is admitted with exact consent/license/ \
                    redaction policy, raw authorized streams and calibrations, blinded partitions, \
                    design-execution records, checker ownership, retention and third-party replay \
                    constraints",
        expiry: "before the first I10 maximal design-campaign submission; review on every manifest \
                 amendment and laboratory/license/calibration/rig revision",
        promotion_effect: "synthetic schema/identifiability/design engineering may proceed, but no \
                           physical-lab identifiability, real-coupon design gain, lab-validated \
                           discrimination or maximal I10 promotion may be claimed while this waiver \
                           is live; baseline synthetic authority remains separately eligible",
    }]
}
