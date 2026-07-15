//! The I06 material-lot passport and substitution-impact engine
//! VerificationManifest draft (bead
//! frankensim-leapfrog-2026-program-i94v.3.2.7.1).
//!
//! The baseline freezes non-confusable material/lot/specimen/process and
//! custody identities; method- and condition-qualified property ingestion;
//! hierarchical lot/spatial/aging uncertainty; ContextOfUse-specific
//! substitution decisions; exact dependency invalidation; anytime-valid
//! supplier drift; and replayable decision receipts. Stronger lanes target
//! coupled multi-property posteriors, identifiable causal transport,
//! transitively complete substitution impact, and robustly optimal material
//! replacement. A refutation lane attacks forged provenance and false impact
//! certificates. Preregistration is governance, not proof of a material,
//! supplier, posterior, causal graph, or replacement decision.

use crate::{
    Ambition, CampaignTier, ClaimPolarity, ClaimSpec, FiveExplicits, FixturePin, FixtureSource,
    GauntletTier, ManifestDraft, ObligationRow, OracleRoute, Partition, ToleranceSemantics, Waiver,
};

const CAMPAIGN_POLICY_FIXTURE: &str = "i06-campaign-policy-v1";

/// Build the I06 draft. Consumers freeze it themselves; focused
/// conformance tests prove that the authored seed satisfies schema v2.
#[must_use]
pub fn i06_draft() -> ManifestDraft {
    ManifestDraft {
        initiative: "I06",
        title: "Material-lot passport and substitution-impact gate: method-qualified properties, \
                hierarchical evidence, causal deltas, selective reproof, and adversarial provenance",
        version: 1,
        explicits: FiveExplicits {
            units: "SI coherent values plus an exact source-unit/scale/offset record; each property binds \
                    semantic quantity kind, tensor/order, frame/orientation, basis and method/condition domain. \
                    Hardness scales and categorical grades are typed non-convertible ordinals; contact angle \
                    distinguishes advancing/receding/static; exact/refusal/coverage verdicts use unit 'bit'",
            seeds: "Philox 4x32-10 streams keyed by BLAKE3 domain \
                    org.frankensim.i06.fixture-stream.v1 over exact fixture, lot-family, stage and case ids; \
                    development indices 0..=16383, core held-out 65536..=81919, maximal held-out \
                    131072..=147455; authenticated raw draws and blind partition assignments are retained; \
                    duplicate logical ids or cross-stage reads are IntegrityFailed",
            budgets: "smoke <= 90 s and 2 GiB; core <= 60 min and 16 GiB; maximal <= 16 h and 32 GiB; \
                      <= 10^6 passport nodes, 10^7 property observations, 10^8 dependency edges, 10^5 \
                      posterior draws per fit, and 10^7 candidate-action/scenario evaluations; caps are \
                      preflighted and exhaustion is BudgetExhausted/Unknown, never compatibility or proof",
            versions: "fs-vmanifest schema v2; MaterialPassport v1; PropertyObservation v1; \
                       ContextOfUse v1; UncertaintyOwnership v1; SubstitutionDeltaGraph v1; \
                       SelectiveReproofReceipt v1; rust-toolchain.toml and constellation.lock receipt-bound",
            capabilities: "safe Rust; no network or FFI in adjudication; deterministic mode; exact unit/ \
                           quantity-kind registry, method/calibration/custody schemas, signed artifact support, \
                           content-addressed raw-data retention and redaction reasons mandatory; no implicit \
                           property conversion, supplier trust, grade equivalence, missing-data fill, causal \
                           transport, or substitution approval; maximal theorem/decision lanes feature-gated",
        },
        claims: i06_claims(),
        fixtures: i06_fixtures(),
        obligations: i06_obligations(),
        waivers: i06_waivers(),
        amendment_rules: "Every changed material/spec/lot/specimen/process/custody/property/method/calibration/ \
                          context/model/prior/dependency/threshold/seed/checker/decision field creates the exact \
                          next manifest version through FrozenManifest::amend. Candidate and blind partitions \
                          freeze before evidence. The amendment record plus the typed substitution graph names \
                          affected predecessor authority; no result edits this version, and unchanged evidence is \
                          rebound only through authenticated byte-identical lineage.",
    }
}

#[allow(clippy::too_many_lines)]
fn i06_claims() -> Vec<ClaimSpec> {
    vec![
        ClaimSpec {
            id: "i06-passport-identity-custody",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "MaterialDefinition, specification revision, supplier site, heat/batch/lot, lot position, \
                        specimen/coupon, process step, storage/aging interval, custody transfer, calibration and raw \
                        artifact identities round-trip canonically; conflicting, missing, replayed or impossible \
                        lineage fails closed without merging distinct physical authority",
            hypotheses: &[
                "every identity has a typed issuer namespace, canonical preimage, validity interval and supersession/ \
                 split/merge semantics; physical labels and electronic records remain separately observed",
                "signatures authenticate bytes/key possession only; physical truth, supplier honesty and label-to- \
                 specimen correspondence require independent observations",
            ],
            qoi: "passport_identity_and_custody_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G0,
            oracle: OracleRoute {
                identity: "independent canonical decoder, custody-state-machine checker and signature/artifact binder \
                           over a separately generated lineage corpus",
                independent: true,
                tcb_overlap: "shares cryptographic primitives and public schemas, not production passport assembly \
                              or conflict resolution",
            },
            activation: "the MaterialPassport v1 implementation leaf opens",
            kill: "one authority alias, silently merged conflict, accepted impossible custody transition, stale key/ \
                   artifact binder, or loss of source bytes blocks every downstream I06 claim",
            fallback: "quarantine the affected record/lot branch with typed UnresolvedIdentity; preserve raw bytes and \
                       permit no substitution decision from that branch",
            no_claim: "canonical custody records do not prove physical composition, absence of counterfeiting, legal \
                       title, supplier qualification or conformance to a material specification",
        },
        ClaimSpec {
            id: "i06-property-observation-semantics",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Every COA/metrology observation binds an unambiguous property kind, value/censoring/ \
                        detection-limit semantics, exact source and SI conversion, uncertainty/covariance, method, \
                        calibration chain, specimen geometry/location/orientation, temperature, pressure, frequency, \
                        strain/shear rate, humidity, surface/process/history state and validity domain",
            hypotheses: &[
                "property registry distinguishes electrical/thermal conductivity, intrinsic permeability/gas \
                 permeability/permeance/diffusivity, magnetic moment/magnetization/permeability/coercivity, \
                 dynamic/kinematic/non-Newtonian viscosity, latent heat per mass/volume, density, ductility and \
                 hardness scale",
                "surface observations distinguish roughness/contamination and advancing/receding/static contact angle; \
                 anisotropic tensors retain frame, symmetry and component covariance",
                "missing, below/above-detection, interval-censored, qualitative and ordinal values are not coerced to \
                 exact numeric zeros or midpoints",
            ],
            qoi: "typed_observation_roundtrip_and_conversion_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G0,
            oracle: OracleRoute {
                identity: "independent dimension/quantity-kind/method registry checker plus exact-rational source-unit \
                           conversion and calibration-lineage replay",
                independent: true,
                tcb_overlap: "shares the frozen registry bytes but no production parser, conversion or ingestion path",
            },
            activation: "passport identity/custody passes for the referenced specimen and artifacts",
            kill: "quantity-kind aliasing, unit/frame loss, invalid calibration, censored-data coercion, or a converted \
                   interval failing exact-rational containment quarantines the observation and dependent decisions",
            fallback: "retain the raw signed record as Uninterpreted/Unsupported with a required mapping task; do not \
                       guess a unit, method, scale or property equivalence",
            no_claim: "semantic ingestion does not endorse the method, calibrate the instrument, prove traceability, \
                       establish representativeness or imply a material-law parameter outside its validity domain",
        },
        ClaimSpec {
            id: "i06-hierarchical-lot-posterior",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "For the frozen hierarchical model family, inference keeps supplier/site, heat/lot, spatial \
                        position, specimen, process, aging, environment, measurement error and model discrepancy at \
                        their declared levels; generated posteriors reproduce pinned exact/analytic moments and achieve \
                        preregistered held-out interval coverage without pseudoreplication",
            hypotheses: &[
                "likelihood, priors/hyperpriors, covariance parameterization, censoring/missingness mechanism and \
                 exchangeability/partial-pooling groups freeze before held-out access",
                "aleatory lot/spatial variation, epistemic parameter uncertainty, measurement uncertainty and model \
                 discrepancy have explicit owners and are not summed twice or silently collapsed",
                "coverage is property/ContextOfUse specific; optional model selection and stopping use the frozen \
                 multiplicity/e-process policy; acceptance uses exact finite-sample rank arithmetic or a preregistered \
                 one-sided confidence/e-value upper bound on coverage error, never a raw empirical proportion alone",
            ],
            qoi: "worst_heldout_interval_coverage_gap",
            unit: "1",
            tolerance: ToleranceSemantics::Absolute { atol: 0.03 },
            evidence_tier: GauntletTier::G2,
            oracle: OracleRoute {
                identity: "independent exact-enumeration/conjugate-model oracle and retained simulation-based- \
                           calibration plus held-out coverage adjudicator",
                independent: true,
                tcb_overlap: "shares frozen model/data bytes only; inference, random streams and coverage code are \
                              independently implemented",
            },
            activation: "typed observation ingestion is green and uncertainty ownership is complete",
            kill: "posterior moment mismatch outside its directed band, coverage gap above 0.03, unowned/double-counted \
                   uncertainty or lot/specimen pseudoreplication refutes this model receipt",
            fallback: "widen to a conservative nonparametric/interval lot envelope or report PosteriorUnknown; never \
                       borrow evidence across a rejected exchangeability boundary",
            no_claim: "calibration on synthetic/pinned families is not truth of a chosen prior/model, universal \
                       frequentist coverage, causal attribution or authority for an unobserved supplier/process regime",
        },
        ClaimSpec {
            id: "i06-contextual-substitution-admission",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "A substitution verdict is a three-valued ContextOfUse-specific predicate over exact source/ \
                        candidate lot evidence, requirements and validity domains: Compatible only when every hard \
                        requirement is conservatively satisfied, Incompatible on a witnessed violation, and Unknown \
                        for missing, incomparable, out-of-domain or insufficient evidence",
            hypotheses: &[
                "ContextOfUse freezes part/function, geometry/surface state, process route, load/environment/time/frequency \
                 domains, property transforms, reliability/safety margins, standards/spec revisions and decision date",
                "requirements declare direction, inclusive/exclusive bounds, joint constraints, units/property kinds, \
                 confidence/coverage semantics and whether nominal, worst-case, quantile or posterior probability applies",
                "grade/spec-name equality, visual similarity and overlapping nominal ranges grant no compatibility by name",
            ],
            qoi: "substitution_admission_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G2,
            oracle: OracleRoute {
                identity: "independent exact-rational requirement evaluator and manually adjudicated cross-property \
                           ContextOfUse corpus",
                independent: true,
                tcb_overlap: "shares frozen requirement/property semantics only; graph traversal and decision code are disjoint",
            },
            activation: "passport, observation and applicable posterior authority are green or explicitly Unknown",
            kill: "one false Compatible, hidden hard requirement, favorable missing-value coercion, unit/kind mismatch or \
                   out-of-domain extrapolation blocks the substitution engine",
            fallback: "return Unknown with ranked missing evidence/tests or Incompatible with a minimized witness; retain \
                       the incumbent material and make no automatic procurement/manufacturing change",
            no_claim: "Compatibility is not supplier approval, regulatory/standards certification, process qualification, \
                       universal interchangeability, or a claim that the candidate is optimal",
        },
        ClaimSpec {
            id: "i06-selective-impact-invalidation",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Given the frozen typed dependency graph, a material/lot/property/process change computes the exact \
                        reachable affected set of requirements, laws/parameters, geometry/process constraints, analyses, \
                        controllers, hazards, decisions and evidence; affected receipts are invalidated and unaffected \
                        byte-identical authority is reusable only through authenticated lineage",
            hypotheses: &[
                "every graph edge has typed source/target authority, semantic role, direction, activation predicate, version \
                 and independent construction/check receipt; dynamic/reflection/manual dependencies are declared",
                "exactness is relative to this declared graph; graph-world completeness is a separate maximal claim",
            ],
            qoi: "declared_graph_reachable_impact_set_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "independent graph decoder and exhaustive reachability oracle on bounded mutation graphs plus \
                           retained before/after authority digests",
                independent: true,
                tcb_overlap: "shares typed graph bytes, not production traversal, invalidation or cache code",
            },
            activation: "a ContextOfUse substitution/change proposal and complete declared dependency graph exist",
            kill: "one under-invalidated reachable node, overclaimed reusable changed node, nondeterministic impact set, or \
                   stale receipt accepted without authenticated lineage refutes selective reproof",
            fallback: "invalidate the entire ContextOfUse evidence closure conservatively and rerun; never preserve a receipt \
                       merely because its file path or display name is unchanged",
            no_claim: "exact declared-graph reachability is not proof that every real causal/organizational dependency was modeled",
        },
        ClaimSpec {
            id: "i06-anytime-valid-supplier-drift",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Sequential supplier/lot/process monitoring uses preregistered nonnegative e-processes whose null \
                        expectation is at most one under the declared filtration; optional stopping and multiple property/ \
                        supplier surveillance use frozen e-BH/family allocation, with Alarm, NoAlarm, DataInvalid and \
                        ModelUnknown remaining distinct",
            hypotheses: &[
                "logical arrival order, filtration, null/alternative, betting functions, truncation, reset/change policy, \
                 missing/censored/late/revised record treatment and multiplicity family freeze before monitoring",
                "selection into tested lots and supplier reporting delay satisfy the declared missingness/selection model or \
                 the lane returns ModelUnknown",
                "certifier validation uses exhaustive finite microcases plus a preregistered time-uniform binomial e-bound \
                 for simulation families; its upper bound, not the raw observed false-alarm fraction, is the acceptance QoI",
            ],
            qoi: "anytime_valid_upper_bound_on_null_false_alarm_probability",
            unit: "1",
            tolerance: ToleranceSemantics::Interval { lo: 0.0, hi: 0.05 },
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "independent exact/log-domain e-process implementation and null simulations with adversarial \
                           stopping, delayed/censored observations and multiplicity",
                independent: true,
                tcb_overlap: "shares only retained canonical observations and frozen policy bytes",
            },
            activation: "typed observations arrive under a frozen monitoring policy and filtration",
            kill: "null false-alarm rate above 0.05, invalid supermartingale arithmetic, post-hoc family/reset editing, \
                   duplicated observation or selection-model violation blocks drift authority",
            fallback: "freeze automated alarms and surface ModelUnknown/DataInvalid with a diagnostic; conventional fixed- \
                       horizon summaries remain descriptive only",
            no_claim: "NoAlarm proves neither no physical drift nor supplier conformance; Alarm localizes statistical \
                       evidence, not physical cause or fraud",
        },
        ClaimSpec {
            id: "i06-replayable-substitution-decision",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Every user-facing substitution result carries an immutable decision receipt binding candidate/ \
                        incumbent lots, ContextOfUse, evidence cutoff, model/requirements/utility versions, uncertainty \
                        owners, hard-constraint witnesses, changed/invalidated/reused authority, alternatives, fallback and \
                        exact replay; deterministic tie-breaking and ranked missing-evidence actions are agent/human parity",
            hypotheses: &[
                "the decision consumes only evidence authorized at its cutoff and records inaccessible/redacted inputs",
                "manual overrides are signed separate decisions with reason, owner, scope and expiry; they do not rewrite \
                 analytic compatibility or posterior receipts",
            ],
            qoi: "decision_receipt_replay_and_explanation_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G5,
            oracle: OracleRoute {
                identity: "independent receipt decoder/replayer and authority-lineage checker with representative human/agent \
                           explanation tasks",
                independent: true,
                tcb_overlap: "shares frozen artifacts and public schema, not production planner or UI rendering",
            },
            activation: "contextual admission and impact invalidation return terminal receipts",
            kill: "non-replay, hidden evidence, unstable tie, stale authority, unexplained override or divergent human/agent \
                   semantics blocks decision publication",
            fallback: "publish an explicit NoDecision/Unknown receipt and the ranked evidence plan; no hidden default material",
            no_claim: "a transparent decision is not necessarily correct, economically optimal, approved by an engineer, \
                       supplier-qualified or safe outside the exact ContextOfUse",
        },
        ClaimSpec {
            id: "i06-coupled-property-posterior",
            ambition: Ambition::Frontier,
            polarity: ClaimPolarity::Affirmative,
            statement: "A coupled posterior jointly models declared mechanical, electrical, magnetic, thermal, transport, \
                        tribological, wetting and aging properties with cross-property covariance, tensor/frame constraints, \
                        positive/monotone thermodynamic admissibility and method-specific observation operators, and achieves \
                        preregistered multivariate held-out calibration",
            hypotheses: &[
                "the finite property vector, transformations, structural zeros, covariance prior, physical constraints and \
                 method/discrepancy operators freeze before held-out access",
                "joint identifiability diagnostics distinguish data-supported correlation from prior/constraint-induced coupling",
                "calibration acceptance uses exact multivariate rank tests where enumerable or a preregistered one-sided \
                 confidence/e-value bound; an unadjusted heldout score alone carries no population-calibration authority",
            ],
            qoi: "worst_multivariate_calibration_and_constraint_violation",
            unit: "1",
            tolerance: ToleranceSemantics::Interval { lo: 0.0, hi: 0.05 },
            evidence_tier: GauntletTier::G2,
            oracle: OracleRoute {
                identity: "independent constrained synthetic generator, exact low-dimensional posterior oracle and multivariate \
                           rank/energy calibration checker",
                independent: true,
                tcb_overlap: "shares frozen property definitions only; inference and calibration implementations are disjoint",
            },
            activation: "feature i06-coupled-posterior; scalar/hierarchical baseline and property semantics are green",
            kill: "calibration/constraint violation above 0.05, covariance non-PSD, frame inconsistency, hidden method pooling \
                   or unsupported identifiable-correlation claim refutes the coupled model",
            fallback: "factor into independently calibrated marginal/conditional envelopes and mark cross-property decisions Unknown",
            no_claim: "calibrated correlation is not causation, universal material law, extrapolation beyond validity domains or \
                       proof that a finite property vector captures every substitution-relevant behavior",
        },
        ClaimSpec {
            id: "i06-identifiable-causal-transport",
            ambition: Ambition::Frontier,
            polarity: ClaimPolarity::Affirmative,
            statement: "For frozen structural causal/selection diagrams and explicitly admissible interventions, the engine \
                        derives only identifiable source-to-target lot/process effects, emits a symbolic identifying functional \
                        with positivity/overlap and transport premises, and returns CausalUnknown when no valid derivation exists",
            hypotheses: &[
                "variables, time/order, structural equations or graph Markov assumptions, latent/confounded edges, selection \
                 nodes, environments, interventions, consistency/SUTVA scope and positivity domains freeze before data",
                "graph discovery output is never silently treated as the causal graph; sensitivity bounds cover declared \
                 unmeasured-confounding alternatives",
            ],
            qoi: "identification_and_transport_derivation_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "independent do-calculus/selection-diagram derivation checker plus finite SCM countermodel pairs for \
                           non-identifiability",
                independent: true,
                tcb_overlap: "shares canonical graph/query bytes but no production identification or simplification code",
            },
            activation: "feature i06-causal-transport; baseline observations/posteriors and target ContextOfUse are frozen",
            kill: "a claimed identifying functional disagrees on a countermodel pair with identical observed distribution, or \
                   omits positivity/selection/consistency premises",
            fallback: "return CausalUnknown with experimental/OED requirements and conservative associational sensitivity bounds",
            no_claim: "identifiability under an assumed graph does not validate that graph, prove effect homogeneity, identify an \
                       individual counterfactual or authorize transport outside the named source/target environments",
        },
        ClaimSpec {
            id: "i06-transitive-impact-completeness",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Affirmative,
            statement: "A machine-checked theorem relates a complete versioned machine/requirement/evidence dependency-and-causal \
                        model to substitution deltas and proves that the returned impact closure contains every authority whose \
                        truth or applicability can change, while every retained authority has a preserved semantic witness",
            hypotheses: &[
                "a pre-proof successor freezes canonical graph/SCM/authority/proposition/definition ASTs, semantics of change, \
                 applicability and preservation, total runtime-premise mapping, and deterministic proof translation",
                "reflection, generated code, dynamic queries, external/manual decisions and physical dependencies have complete \
                 reified adapters or are explicit open-world boundaries",
                "exact axiom allowlist and transitive closure, independent model decoder and nonvacuity witnesses are receipt-bound",
            ],
            qoi: "kernel_checked_impact_completeness_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G5,
            oracle: OracleRoute {
                identity: "pinned Lean kernel replay plus independent canonical-model/proposition translation and runtime-premise checker",
                independent: true,
                tcb_overlap: "kernel/axioms and frozen semantic model are declared TCB; no production impact traversal/proof search",
            },
            activation: "feature i06-impact-completeness-theorem; a pre-proof successor freezes all machine artifacts and \
                         baseline/frontier premises are independently green",
            kill: "open-world dependency, missing adapter/premise, translation mismatch, disallowed axiom, kernel rejection or \
                   genuine counterexample leaves this claim Unknown/Refuted and blocks completeness wording",
            fallback: "retain exact reachability relative to the declared graph and conservatively invalidate the whole ContextOfUse closure",
            no_claim: "version-1 prose mints no theorem or real-world completeness authority; complete graph-world impact is not \
                       omniscience about unknown physics, organizations, fraud or future requirements",
        },
        ClaimSpec {
            id: "i06-robust-decision-optimality",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Affirmative,
            statement: "Over a frozen finite candidate/action/test plan and ambiguity set, the engine returns a globally optimal \
                        feasible replacement/no-replacement/evidence-acquisition policy for declared lifecycle cost, performance, \
                        reliability, supply and safety utility, with coincident checker-verified global lower and feasible upper bounds",
            hypotheses: &[
                "candidate set, incumbent, ContextOfUse, hard constraints, utility/loss, time horizon, discounting, information \
                 actions, ambiguity/posterior set, scenario probabilities and decision-maker risk functional freeze before results",
                "every feasibility and utility evaluation is bound to the same property/impact evidence and returns interval/ \
                 Unknown rather than a favorable point for missing physics",
                "globality is only over the finite encoded decision grammar; a pre-candidate successor freezes enumeration/ \
                 relaxation, validity, bound and completeness-check semantics; all bound/utility arithmetic and the frozen \
                 positive reporting scale are exact or outward-rounded so a zero normalized gap cannot hide rounding",
            ],
            qoi: "certified_normalized_global_optimality_gap",
            unit: "1",
            tolerance: ToleranceSemantics::Interval { lo: 0.0, hi: 0.0 },
            evidence_tier: GauntletTier::G5,
            oracle: OracleRoute {
                identity: "independent finite-action exhaustive oracle on microcases plus primal-feasibility and dual/lower-bound checker",
                independent: true,
                tcb_overlap: "shares frozen decision/evaluation bytes, not production search, surrogate or bound implementation",
            },
            activation: "feature i06-robust-decision; finite grammar and every consumed evidence authority are frozen",
            kill: "infeasible recommendation, missing action, underestimated uncertainty, inconsistent lower/upper bounds, \
                   any positive certified gap or independently better admissible policy refutes global optimality",
            fallback: "return the best verified feasible incumbent with its honest gap/Unknown state and ranked evidence plan; no auto-substitution",
            no_claim: "optimality is not universal, value-neutral or regulatory approval and cannot compensate for a misspecified \
                       utility, ambiguity set, candidate grammar, causal graph, missing hazard or invalid evidence",
        },
        ClaimSpec {
            id: "i06-false-provenance-impact-certificate-falsifier",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Refutation,
            statement: "Hidden provenance, semantic, statistical, graph, causal and decision mutants try to obtain a green \
                        passport/substitution/impact/optimality certificate for a known invalid case; any accepted mutant refutes \
                        the corresponding certifier and is retained as a minimized counterexample",
            hypotheses: &[
                "mutants include forged/replayed/stale signatures, specimen/lot relabeling, duplicate raw data, unit/property/ \
                 hardness-scale spoofing, censoring loss, calibration break, selection leak, pseudoreplication, dependency-edge \
                 omission/double ownership, stale evidence reuse, nonidentifiable causal query and hidden infeasible action",
                "each mutant has an independently authored cryptographic, exact-rational, lineage, graph-reachability, SCM \
                 countermodel or finite-decision witness hidden until maximal adjudication",
            ],
            qoi: "false_certificate_count",
            unit: "1",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "independent mutation generator and witness checker with no production parser, inference, graph or optimizer code",
                independent: true,
                tcb_overlap: "public schemas and cryptographic primitives only",
            },
            activation: "a maximal certifier candidate and all checker/version identities are sealed before mutant reveal",
            kill: "the intended lane succeeds when a false certificate is accepted: mark the certifier Refuted and block \
                   dependent promotion; zero accepted mutants is finite falsification evidence only",
            fallback: "disable the affected certifier and use the strongest independently surviving baseline with Unknown/no-decision",
            no_claim: "a finite adversarial corpus cannot prove absence of forgery, data poisoning, semantic confusion, omitted \
                       dependency, causal misspecification or optimizer unsoundness",
        },
    ]
}

#[allow(clippy::too_many_lines)]
fn i06_fixtures() -> Vec<FixturePin> {
    vec![
        FixturePin {
            id: CAMPAIGN_POLICY_FIXTURE,
            source: FixtureSource::AuthoredSpec {
                spec: "I06_CAMPAIGN_POLICY_V1\n\
IDENTITY=MaterialDefinition,SpecRevision,SupplierSite,Heat,Batch,Lot,LotPosition,Specimen,Coupon,ProcessStep,StorageInterval,CustodyTransfer,Instrument,Calibration,RawArtifact and Decision identities have distinct typed domains, canonical preimages, validity/supersession and split/merge rules; same display text never aliases authority\n\
AUTHENTICITY=a valid signature proves byte/key binding only; physical label, specimen correspondence, composition, method execution and supplier truth require independent evidence; stale/revoked/unauthorized keys, replay, inconsistent custody and conflicting signed claims are explicit IntegrityFailed/Conflict states\n\
PROPERTY=bind quantity kind,tensor/order,frame/orientation,basis,value/censoring/detection limits,source/SI units,uncertainty/covariance,method,calibration,specimen geometry/location,temperature,pressure,frequency,strain/shear rate,humidity,surface/process/history and validity domain; electrical vs thermal conductivity, intrinsic/gas permeability vs permeance/diffusivity, magnetic moment/magnetization/permeability/coercivity, dynamic/kinematic/non-Newtonian viscosity, mass/volume latent heat, density, ductility, ordinal hardness scales and advancing/receding/static contact angle are non-confusable\n\
UNCERTAINTY={aleatory lot/spatial/process/aging,epistemic parameter,measurement,model discrepancy,scenario}; each has one owner and covariance/conditioning rule; missing/censored/qualitative observations remain typed; exchangeability, selection, priors, likelihood and optional model/stopping policy freeze before holdout\n\
SUBSTITUTION=verdict {Compatible,Incompatible,Unknown,IntegrityFailed} is only for exact ContextOfUse, evidence cutoff and requirement versions; hard constraints precede utility; direction/inclusivity/joint probability semantics are explicit; grade/spec-name equality or nominal overlap grants no equivalence\n\
IMPACT=declared graph edges bind typed source/target authority, semantic role, activation and version; changed/reachable/reused/inapplicable/unknown states are distinct; exact selective invalidation is relative to the graph, while real-world graph completeness is a separately gated theorem target\n\
DRIFT=logical filtration, null/alternative, e-process bets/truncation, missing/censoring/revision, reset/change policy, supplier/property family and e-BH allocation freeze first; Alarm,NoAlarm,DataInvalid,ModelUnknown remain distinct; NoAlarm never proves no drift\n\
CAUSAL=structural/selection graph,variables,time,latent edges,environments,interventions,query,consistency/SUTVA,positivity and sensitivity alternatives freeze first; learned association is not a causal graph; unidentified queries return CausalUnknown\n\
DECISION=finite candidates/actions/tests,incumbent,hard constraints,utility/loss,risk functional,horizon,discount,ambiguity/posterior/scenarios and tie break freeze first; feasible incumbent,lower bound,upper bound,gap and Unknown are retained; no automatic procurement/manufacturing action\n\
THEOREM_AUTHORITY=version 1 prose mints no completeness,causal or global-optimality authority; pre-proof successors freeze canonical graph/SCM/authority/proposition/definition ASTs, semantics, total runtime premises, deterministic AST-to-Lean translation, exact axiom allowlist {propext,Quot.sound,Classical.choice} and transitive closure; sorryAx, custom postulates and native/unsafe oracle shortcuts are IntegrityFailed; pre-search successors freeze finite grammar,validity,enumeration/bounds,rank-unrank/sharding and completeness root\n\
EVIDENCE_STATES=Execution{Succeeded,Failed,Cancelled,TimedOut,BudgetExhausted,InfrastructureFailed,IntegrityFailed}; Predicate{Accepted,Refuted,Unknown,NotEvaluated}; Claim{NoClaim,EvidencePartial,Reproduced,Refuted,Proved}; Decision{Compatible,Incompatible,Unknown,NoDecision}; Promotion{Blocked,CoreEligible,MaxEligible}; axes never substitute\n\
HOLDOUT=candidate,models,graphs,targets,thresholds,seeds,checker/policy versions and decision grammar freeze before access; each Core/Max heldout fixture has one named consumer stage and disjoint range; premature/cross-stage read,label leak,replacement,retry-after-result or post-result tuning is IntegrityFailed\n\
LIFECYCLE=request->drain->finalize; poll at record,lot,posterior-chain,graph-frontier,monitor-update,scenario and proof/search tile boundaries; checkpoints bind manifest/candidate/evidence cutoff and completed logical ids; resume/fork cannot duplicate observations,posterior draws,e-process wealth or published decisions\n\
LOGGING=bounded schema-versioned fs-obs JSONL with stable run/case/claim/fixture/leaf/material/spec/lot/specimen/process/property/method/calibration/context/model/graph/decision/checker/attempt ids, exact units,seeds,budgets,versions,capabilities and redaction reason; stdout is not evidence\n\
RETENTION=manifest,exact command,code/contract/schema/toolchain/registry hashes,raw authorized bytes/signatures/calibrations,canonical observations,model/prior/graph/context/requirement/decision bytes,draws/diagnostics,all terminal states,minimized counterexamples and replay-verifier result; promotion/refutation/decision evidence durable; inaccessible/redacted/expired inputs constrain replay and never become silently verified\n\
FAILURE_BUNDLE=first semantic/inference/impact/decision divergence plus bounded context,identity/custody chain,raw artifact,expected/actual unit/value/interval/state,causal predecessors,minimized witness,checker disagreement,terminal state and replay command; partial success cannot publish normal authority\n\
PROMOTION=baseline claims only I06.G4 after G2 reproduction/G3 falsification; maximal claims only I06.G7 after G5 reproduction/G6 red-team; missing,stale,waived,integrity-failed or inaccessible evidence cannot promote\n\
LEAF_REQUIREMENT=every obligation references this policy fixture,all nine unit classes,smoke/core/max tier,DSR lane,events,replay,G4 drain/checkpoint,G5 matrix,independent adjudication,performance envelope,accessibility/agent parity and consuming gate",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i06-passport-custody-corpus",
            source: FixtureSource::AuthoredSpec {
                spec: "I06_FIXTURE_V1 passport/custody state machines: supplier sites, heats, batches, lots, \
                       locations, cuts/splits/merges, specimens/coupons, process/storage/aging intervals, custody \
                       transfers, instruments/calibrations and raw artifacts. Valid cases include partial use and \
                       recall/supersession. Invalid classes include namespace/display collisions, impossible time/ \
                       mass ancestry, duplicate specimen, orphan split, cycle, overlapping exclusive custody, stale/ \
                       revoked/wrong-scope key, replayed signature/artifact, conflicting signed COAs and electronic/ \
                       physical-label disagreement. Generator emits an independent validity/conflict witness.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i06-property-ontology-v1",
            source: FixtureSource::AuthoredSpec {
                spec: "I06_FIXTURE_V1 property ontology and conversion KATs. KINDS: density; heat capacity; \
                       thermal/electrical conductivity; resistivity; dielectric permittivity/loss; magnetic moment, \
                       magnetization, permeability, remanence and coercivity; elastic/plastic/creep/fatigue/fracture \
                       parameters; tensile elongation/ductility; Rockwell/Vickers/Brinell hardness as distinct \
                       method-qualified scales; dynamic/kinematic and shear-rate-dependent viscosity; intrinsic porous \
                       permeability, gas permeability, permeance, diffusivity; vapor pressure; latent heat per \
                       mass/volume; surface energy and advancing/receding/static contact angle. Tensor/frame, temperature, \
                       pressure, frequency, rate, humidity, roughness, contamination, process/history domains explicit. \
                       Exact-rational affine/multiplicative conversions plus forbidden cross-kind conversion corpus.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i06-metrology-ingestion-corpus",
            source: FixtureSource::AuthoredSpec {
                spec: "I06_FIXTURE_V1 synthetic COA/LIMS/metrology records in bounded canonical CSV/JSON/binary forms: \
                       exact, repeated, interval, left/right-censored, below detection, qualitative and ordinal values; \
                       significant digits; correlated calibration covariance; instrument drift; specimen geometry and \
                       spatial coordinates; method/standard edition; unit/frame/orientation. Includes decimal/binary \
                       locale traps, offset temperatures, percent/fraction, mass/volume normalization, duplicate raw \
                       rows, missingness, transposed tensor basis, hardness-scale spoof and calibration discontinuity. \
                       Raw bytes and independently authored normalized records are both pinned.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i06-hierarchical-synthetic-families",
            source: FixtureSource::AuthoredSpec {
                spec: "I06_FIXTURE_V1 hierarchical generators with analytic/conjugate microcases and nonconjugate \
                       simulation truth: suppliers 2..5, sites 1..3, heats 2..20, lots/heat 1..8, spatial positions \
                       3..50, specimens/position 1..5, aging times 0..8 and instruments 1..4; crossed/nested effects, \
                       heteroscedastic/correlated/censored measurements, selection and missingness. Truth separately \
                       owns lot/spatial/process/aging/measurement/discrepancy variance. Includes identical-row \
                       pseudoreplication traps, exchangeability breaks, nonidentifiable variance allocation and prior- \
                       dominated small-lot cases. Development seeds 0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i06-material-family-matrix",
            source: FixtureSource::AuthoredSpec {
                spec: "I06_FIXTURE_V1 synthetic material-shaped families: carburized gear steels with heat treatment, \
                       hardness depth, fatigue/fracture and viscosity-dependent lubricant interface; NdFeB/ferrite \
                       magnet lots with remanence/coercivity/permeability vs temperature and coating/corrosion; polymer/ \
                       insulation lots with dielectric, thermal, moisture, gas permeability/diffusivity, aging and \
                       partial-discharge context; cold-plate alloys with conductivity, strength, corrosion and process \
                       state; porous seals/membranes with intrinsic permeability versus permeance; wetting surfaces with \
                       roughness/contamination and advancing/receding angle; phase-change media with density and latent \
                       heat per mass/volume. Every number is synthetic and method/domain-qualified; no supplier/spec claim.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i06-dependency-delta-graphs",
            source: FixtureSource::AuthoredSpec {
                spec: "I06_FIXTURE_V1 bounded typed graphs with material/lot/property/process nodes feeding requirements, \
                       constitutive laws, geometry/surface/process constraints, meshes/solvers, controllers, hazards, \
                       standards clauses, decisions and evidence. Graphs 1..10000 nodes and up to 100000 edges include \
                       diamonds, cycles condensed by typed SCC rules, conditional edges, shared evidence, unchanged-byte \
                       lineage and manual/external open-world boundaries. Mutations have exact exhaustive reachability \
                       sets. Invalid graphs cover missing/double owner, dangling authority, wrong direction/type, stale \
                       activation and hidden dynamic/reflection dependency declarations.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i06-drift-null-alternative-families",
            source: FixtureSource::AuthoredSpec {
                spec: "I06_FIXTURE_V1 sequential supplier/property streams with exact/log-domain e-process KATs: IID \
                       and conditionally valid nulls; gradual/abrupt mean, variance, tail, covariance and censoring drift; \
                       delayed/revised/duplicate/missing/selected reports; suppliers/properties entering/leaving families; \
                       adversarial stopping and reset attempts. Logical ids/order, filtration, bets, truncation, family/ \
                       e-BH allocation and change policy pinned. Null simulations target alpha=0.05; alternatives measure \
                       delay descriptively. Development seeds 4096..=8191.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i06-impact-causal-theorem-card",
            source: FixtureSource::AuthoredSpec {
                spec: "I06_THEOREM_CARD_V1 TARGET. Intended successor proposition relates typed versioned authority graph \
                       G, causal applicability model C, semantic change delta D and impact closure I(D), proving every node \
                       whose truth/applicability differs is in I and every retained node has a preservation witness. It \
                       also freezes source/target selection diagrams and identifiable transport query semantics. REQUIRED \
                       MACHINE ARTIFACTS: graph/SCM/authority/proposition/definition ASTs and canonical bytes, complete \
                       adapters/open-world boundaries, total premise map, deterministic Lean translation/roundtrip, exact \
                       axiom allowlist {propext,Quot.sound,Classical.choice}, transitive closure, nonvacuity witnesses and \
                       retained kernel replay. Prose grants no theorem, graph-completeness or causal authority.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i06-decision-grammar-card",
            source: FixtureSource::AuthoredSpec {
                spec: "I06_EXHAUSTIVENESS_CARD_V1 TARGET. Intended successor freezes finite incumbent/candidate/ \
                       substitute/no-substitute/test/acquire-information policy grammar, canonical action/scenario state, \
                       feasibility and interval utility semantics, ambiguity/posterior set, risk functional, validity and \
                       exclusion order, exact enumeration or verified relaxation bounds, rank/unrank/sharding, independent \
                       decoder/checker, preflight and completeness root. Returned feasible upper bound, global lower bound \
                       and normalized gap are separately replayed. Prose grants no global-optimality or decision authority.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i06-provenance-adversaries-core-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "I06_FIXTURE_V1 CORE HOLDOUT provenance/semantic adversaries, indices 65536..=69631: valid-looking \
                       signature with wrong issuer scope, revoked/stale key, byte replay under new lot, specimen swap, \
                       impossible split/merge/custody time, duplicate raw data, source/SI unit spoof, tensor-frame swap, \
                       hardness-scale alias, censoring/detection-limit loss, calibration break and signed COA conflict. \
                       Each has independent cryptographic/lineage/exact-rational witness. One I06.G3 consumer.",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "i06-posterior-coverage-core-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "I06_FIXTURE_V1 CORE HOLDOUT hierarchical families, indices 69632..=73727: unseen supplier/site/ \
                       lot/spatial/aging combinations, censoring/missingness, correlated instruments, exchangeability \
                       breaks and prior-dominated sparse lots. Truth retains all uncertainty-owner components; labels and \
                       blind partition remain inaccessible until one-shot posterior submission. Scores exact/conjugate \
                       moments, SBC/rank and property/context interval coverage. One I06.G3 consumer.",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "i06-substitution-impact-core-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "I06_FIXTURE_V1 CORE HOLDOUT substitution cases, indices 73728..=77823: gear steel, magnet, \
                       insulation polymer, cold-plate alloy, porous seal and wetting-surface candidates crossed with \
                       exact ContextOfUse and dependency graphs. Includes compatible, witnessed-incompatible and Unknown \
                       cases from missing domain/method/joint evidence; property-kind/unit/grade traps; conditional graph \
                       edges and byte-identical reusable siblings. Independent evaluator and exhaustive graph reachability \
                       sets withheld until final submission. One I06.G3 consumer.",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "i06-supplier-drift-core-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "I06_FIXTURE_V1 CORE HOLDOUT supplier streams, indices 77824..=81919: conditionally valid nulls \
                       with adversarial stopping plus gradual/abrupt property, tail, covariance, censoring and selection \
                       changes; delayed/revised/duplicate/missing reports and multiplicity-family churn. Logical order, \
                       stopping adversary and labels are sealed; one-shot output includes wealth path, alarm time/state, \
                       family decisions and DataInvalid/ModelUnknown. One I06.G3 consumer.",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "i06-coupled-properties-max-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "I06_FIXTURE_V1 MAX HOLDOUT coupled-property families, indices 131072..=135167: anisotropic \
                       mechanical/electrical/magnetic/thermal/transport/tribological/wetting/aging vectors with known \
                       structural zeros, PSD covariance, positive/monotone admissibility, method discrepancy and \
                       identifiable versus prior-induced correlation. Exact low-dimensional and simulation-truth cases \
                       score multivariate rank/energy calibration and constraint violations. One I06.G6 consumer.",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "i06-causal-impact-max-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "I06_FIXTURE_V1 MAX HOLDOUT causal/impact cases, indices 135168..=139263: finite SCM and selection- \
                       diagram source/target environments with identifiable, transportable, nonidentifiable and positivity- \
                       violating queries; latent confounding, selection, process interventions and sensitivity pairs; typed \
                       authority graphs with complete adapters and deliberate open-world boundaries. Countermodel pairs \
                       share observed distributions for nonidentifiable queries; graph delta closures and preservation \
                       witnesses are exact. One I06.G6 consumer.",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "i06-decision-certifier-mutants-max-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "I06_FIXTURE_V1 MAX HOLDOUT decision/certifier mutants, indices 139264..=147455: stale/replayed/ \
                       forged provenance, lot/specimen relabel, unit/kind/censoring/covariance mutation, pseudoreplication, \
                       e-process optional-stopping edit, omitted/double dependency, nonidentifiable causal functional, \
                       hidden hard constraint/action/scenario, favorable Unknown point, unsound relaxation/lower bound, \
                       incomplete shard and disallowed theorem axiom. Each mutant has an independent witness and is hidden \
                       until the sealed candidate. One I06.G6 consumer.",
            },
            partition: Partition::HeldOut,
        },
    ]
}

#[allow(clippy::too_many_lines)]
fn i06_obligations() -> Vec<ObligationRow> {
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
            leaf: "i06-passport-observation-ingestion",
            claims_covered: &[
                "i06-passport-identity-custody",
                "i06-property-observation-semantics",
            ],
            unit_cases: UNIT_CASES,
            g0: "generators: passport/custody state machines, signed artifact mutations and typed property/metrology \
                 records; predicates: canonical identity, namespace/issuer/key scope, lineage/custody validity, raw-artifact \
                 binding, quantity-kind/unit/frame/method/domain and calibration containment; laws: roundtrip, no authority \
                 alias, signature replay/staleness refusal, exact-rational conversion, censoring preservation, tensor basis \
                 covariance and conflict non-erasure; shrink lineage/record fields while preserving independent witness",
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                "i06-passport-custody-corpus",
                "i06-property-ontology-v1",
                "i06-metrology-ingestion-corpus",
                "i06-provenance-adversaries-core-holdout",
            ],
            g3_relations: &[
                "canonical id/record order relabeling preserves identity graph and normalized observation",
                "exact compatible unit/frame transform and inverse preserves physical value/covariance",
                "signature/key/lot/specimen/artifact substitution fails before semantic use",
                "duplicating raw bytes or normalized rows cannot manufacture independent observations",
                "censored/ordinal/qualitative data never become an exact zero or midpoint under format migration",
            ],
            g4_schedule: "cancel/panic/timeout/corrupt input after raw retention, decode, signature/key check, custody join, \
                          unit/method/calibration normalization and canonical assembly; poll each record/lineage edge; \
                          request->drain->finalize publishes all bound raw+normalized+conflict artifacts or none; checkpoint \
                          completed logical ids and content roots; resume/fork cannot duplicate records or erase conflicts; \
                          retain bounded FailureBundle and every terminal state",
            g5_matrix: "record/lineage shards {1,2,7} x input orders {capture,reverse,permuted} x deterministic mode on \
                        identical registry/crypto/toolchain fingerprint; canonical passports, observation/covariance bytes, \
                        conflicts, diagnostics, event JSON and digests match bitwise",
            entry_point: "scripts/e2e/leapfrog/i06_ingestion.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i06-ingestion slice)",
            obs_events: &[
                "request.received",
                "passport.record_admitted",
                "passport.conflict",
                "custody.transition",
                "observation.normalized",
                "observation.quarantined",
                "cancel.requested",
                "drain.completed",
                "finalize.completed",
                "failure_bundle.retained",
                "adjudication.receipt",
            ],
            replay_command: "scripts/e2e/leapfrog/i06_ingestion.sh --manifest <manifest-id> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i06-hierarchical-posterior",
            claims_covered: &["i06-hierarchical-lot-posterior"],
            unit_cases: UNIT_CASES,
            g0: "generators: analytic/conjugate and simulation-truth supplier/site/lot/spatial/specimen/process/aging/ \
                 instrument hierarchies with censored/missing/selected observations; predicates: frozen model/prior/ \
                 likelihood/group/covariance/missingness identity, uncertainty-owner completeness, exact-moment and SBC/ \
                 heldout coverage arithmetic; laws: observation permutation, unit covariance, duplicate-id refusal, \
                 group-relabel conjugacy, posterior draw logical-id determinism and no pseudoreplication; shrink hierarchy \
                 while retaining moment/coverage/ownership failure",
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                "i06-hierarchical-synthetic-families",
                "i06-material-family-matrix",
                "i06-posterior-coverage-core-holdout",
            ],
            g3_relations: &[
                "consistent affine unit transform maps posterior moments/intervals/covariance exactly",
                "lot/specimen/spatial id relabeling preserves exchangeability-group posterior up to relabeling",
                "duplicating a measurement with the same logical id is refused rather than narrowing uncertainty",
                "adding an unowned uncertainty source cannot improve coverage or compatibility authority",
                "heldout order and external-id reassignment preserve per-case scores and aggregate coverage",
            ],
            g4_schedule: "cancel/panic/timeout at model build, chain initialization, draw tile, diagnostics, heldout prediction \
                          and reduction; request->drain->finalize publishes no posterior without all chains/owners/diagnostics \
                          and terminal states; checkpoint binds logical chain/draw ids, RNG counters and model/data roots; \
                          resume/fork equals one-shot without duplicated draws; retain nonconvergence/coverage FailureBundle",
            g5_matrix: "chains {1,2,4} x draw shards {1,2,7} x workers {1,2,7} x observation orders {forward,reverse} x \
                        deterministic mode on identical ISA; exact logical RNG draws, summaries, intervals, diagnostics, \
                        heldout scores, tie breaks, events and receipt roots match bitwise",
            entry_point: "scripts/e2e/leapfrog/i06_posterior.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i06-posterior slice)",
            obs_events: &[
                "request.received",
                "posterior.model_admitted",
                "posterior.draw_completed",
                "posterior.diagnostic",
                "posterior.coverage_scored",
                "posterior.unknown",
                "cancel.requested",
                "drain.completed",
                "finalize.completed",
                "failure_bundle.retained",
                "adjudication.receipt",
            ],
            replay_command: "scripts/e2e/leapfrog/i06_posterior.sh --manifest <manifest-id> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i06-substitution-impact",
            claims_covered: &[
                "i06-contextual-substitution-admission",
                "i06-selective-impact-invalidation",
            ],
            unit_cases: UNIT_CASES,
            g0: "generators: material-family incumbents/candidates x exact ContextOfUse/requirement semantics x typed \
                 dependency graphs and deltas; predicates: hard-requirement three-valued evaluation, validity-domain/ \
                 evidence-cutoff binding, exact graph reachability, changed-byte authority and authenticated reuse; laws: \
                 false Compatible impossible on independent evaluator, Unknown monotonic under evidence removal, graph id \
                 relabeling, reachable-set exactness, unchanged sibling preservation, conditional-edge activation and stale \
                 receipt refusal; shrink to minimal requirement/property/edge witness",
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                "i06-property-ontology-v1",
                "i06-material-family-matrix",
                "i06-dependency-delta-graphs",
                "i06-substitution-impact-core-holdout",
            ],
            g3_relations: &[
                "consistent unit/frame/property transform preserves requirement verdict and impact set",
                "tightening a hard admissible bound cannot convert Incompatible/Unknown to Compatible",
                "removing applicable evidence cannot improve a three-valued compatibility verdict",
                "node-id permutation conjugates the exact reachable affected set and deterministic tie breaks",
                "adding a dependency edge can preserve or enlarge but cannot shrink impact closure",
                "byte-identical unrelated sibling authority remains reusable only with authenticated lineage",
            ],
            g4_schedule: "cancel/panic/timeout during requirement evaluation, graph decode, frontier expansion, invalidation, \
                          reuse check and receipt assembly; poll each requirement/frontier tile; request->drain->finalize \
                          publishes one terminal Compatible/Incompatible/Unknown/IntegrityFailed plus complete impact state; \
                          checkpoint binds frontier/visited/activation and content roots; resume/fork equals one-shot; retain \
                          minimized false-admission/under-or-over-invalidation FailureBundle",
            g5_matrix: "graph shards {1,2,7} x workers {1,2,7} x node/edge/requirement orders {forward,reverse,permuted} x \
                        deterministic mode; verdict, witnesses, changed/reachable/reused/unknown sets, ranked fixes, events \
                        and canonical receipt match bitwise across the same ISA/toolchain",
            entry_point: "scripts/e2e/leapfrog/i06_substitution_impact.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i06-substitution-impact slice)",
            obs_events: &[
                "request.received",
                "substitution.context_bound",
                "substitution.compatible",
                "substitution.incompatible",
                "substitution.unknown",
                "impact.node_invalidated",
                "impact.authority_reused",
                "cancel.requested",
                "drain.completed",
                "finalize.completed",
                "failure_bundle.retained",
                "adjudication.receipt",
            ],
            replay_command: "scripts/e2e/leapfrog/i06_substitution_impact.sh --manifest <manifest-id> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i06-drift-decision-replay",
            claims_covered: &[
                "i06-anytime-valid-supplier-drift",
                "i06-replayable-substitution-decision",
            ],
            unit_cases: UNIT_CASES,
            g0: "generators: exact/log e-process streams under null/alternative with stopping/missing/revision/selection \
                 adversaries plus substitution decision receipts/overrides; predicates: filtration/bet/family/reset identity, \
                 nonnegative wealth, adversarial-stopping error, evidence cutoff, hard witnesses, authority lineage, override \
                 separation and replay parity; laws: logical-id/order determinism, no duplicate wealth update, e-BH allocation, \
                 NoAlarm non-authority, stable tie/ranked fixes and decision replay; shrink stream/decision while preserving error",
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                "i06-drift-null-alternative-families",
                "i06-material-family-matrix",
                "i06-supplier-drift-core-holdout",
            ],
            g3_relations: &[
                "null adversarial stopping preserves family type-I error <= 0.05",
                "stream batching/order under frozen logical filtration preserves exact wealth and alarm state",
                "duplicate/revised/late records follow frozen policy and cannot be counted twice favorably",
                "candidate/order/id relabeling preserves deterministic decision and explanation up to relabeling",
                "signed manual override changes decision authority but not analytic compatibility/posterior receipts",
            ],
            g4_schedule: "cancel/panic/timeout at observation admission, wealth update, family decision, receipt join, \
                          explanation and publication; request->drain->finalize persists exact e-process state and one decision \
                          terminal state; checkpoint binds filtration position, logical ids, wealth and evidence cutoff; resume \
                          cannot double-count or edit stopping policy; retain alarm/false-alarm/non-replay FailureBundle",
            g5_matrix: "stream/decision shards {1,2,7} x workers {1,2,7} x authorized batchings {single,chunked} x \
                        deterministic mode on identical ISA; wealth path, family decision, alarm time/state, evidence cutoff, \
                        tie, explanation, events and receipt bytes match bitwise",
            entry_point: "scripts/e2e/leapfrog/i06_drift_decision.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i06-drift-decision slice)",
            obs_events: &[
                "request.received",
                "drift.wealth_updated",
                "drift.alarm",
                "drift.no_alarm",
                "drift.model_unknown",
                "decision.receipt_built",
                "decision.override_recorded",
                "cancel.requested",
                "drain.completed",
                "finalize.completed",
                "failure_bundle.retained",
                "adjudication.receipt",
            ],
            replay_command: "scripts/e2e/leapfrog/i06_drift_decision.sh --manifest <manifest-id> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i06-coupled-property-model",
            claims_covered: &["i06-coupled-property-posterior"],
            unit_cases: UNIT_CASES,
            g0: "generators: constrained multivariate property/tensor/method/aging families with exact low-dimensional and \
                 simulation truth; predicates: registry/frame/method consistency, PSD covariance, positive/monotone/ \
                 thermodynamic constraints, uncertainty ownership, joint identifiability and multivariate calibration; laws: \
                 property/vector permutation covariance, frame rotation, unit transform, structural-zero preservation, \
                 marginalization consistency and no prior-induced correlation overclaim; shrink property graph/dataset",
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                "i06-property-ontology-v1",
                "i06-material-family-matrix",
                "i06-coupled-properties-max-holdout",
            ],
            g3_relations: &[
                "orthogonal frame change conjugates tensors/covariance and preserves scalar decision QoIs",
                "consistent per-property unit transforms preserve joint density/calibration after Jacobian accounting",
                "property permutation conjugates covariance and leaves invariant constraints/scores unchanged",
                "adding an unidentified correlation cannot acquire data-supported authority from a physical constraint alone",
                "heldout order/id relabel preserves multivariate rank/energy and constraint violation scores",
            ],
            g4_schedule: "preflight dimensions/draw budget before heldout access; cancel/panic/timeout at constrained model build, \
                          factorization, draw, projection, diagnostic and heldout scoring tiles; request->drain->finalize publishes \
                          no joint posterior without all owners/constraints/terminal states; checkpoints bind chains/draws/model/ \
                          data roots; resume equals one-shot; retain PSD/constraint/calibration FailureBundle",
            g5_matrix: "chains {1,2,4} x property/draw shards {1,2,7} x workers {1,2,7} x deterministic mode on identical \
                        ISA; logical draws, covariance/constraint records, summaries, calibration scores, tie breaks, events and \
                        receipt root match bitwise",
            entry_point: "scripts/e2e/leapfrog/i06_coupled_properties.sh",
            tier: CampaignTier::Max,
            dsr_lane: "dsr quality --tool frankensim (i06-coupled-properties isolated lane)",
            obs_events: &[
                "request.received",
                "coupled_model.admitted",
                "coupled_model.constraint_checked",
                "coupled_model.draw_completed",
                "coupled_model.calibration_scored",
                "coupled_model.refuted",
                "cancel.requested",
                "drain.completed",
                "finalize.completed",
                "failure_bundle.retained",
                "adjudication.receipt",
            ],
            replay_command: "scripts/e2e/leapfrog/i06_coupled_properties.sh --manifest <manifest-id> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i06-causal-impact-certifier",
            claims_covered: &[
                "i06-identifiable-causal-transport",
                "i06-transitive-impact-completeness",
            ],
            unit_cases: UNIT_CASES,
            g0: "generators: finite SCM/selection diagrams, causal queries/countermodel pairs and typed authority/impact \
                 models plus proposition/definition ASTs; predicates: graph/query canonicality, do/selection derivation, positivity/ \
                 premise binding, nonidentifiability witness, adapter/open-world closure, translation roundtrip, axiom closure, \
                 kernel replay and preservation/nonvacuity witnesses; laws: graph-id alpha-renaming, derivation soundness on \
                 finite SCMs, countermodel refusal, impact theorem specialization and no open-world completeness; shrink model/query/proof",
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                "i06-dependency-delta-graphs",
                "i06-impact-causal-theorem-card",
                "i06-causal-impact-max-holdout",
            ],
            g3_relations: &[
                "graph/variable alpha-renaming preserves identification and impact theorem up to canonical rebinding",
                "source/target environment identity substitution invalidates transport before data comparison",
                "finite countermodel pairs force CausalUnknown for nonidentifiable queries",
                "adding an open-world boundary prevents but cannot strengthen completeness authority",
                "disallowed axiom, missing premise/adapter or semantic/translation mutation fails closed",
                "declared-graph reachability remains available if maximal completeness is Unknown/Refuted",
            ],
            g4_schedule: "whole-campaign preflight precedes max holdout; cancel/panic/timeout/corrupt checkpoint at model decode, \
                          causal derivation, countermodel search, AST translation, kernel replay and witness checking; poll per graph/ \
                          proof tile; request->drain->finalize drains children and preserves terminal state; checkpoints bind frontier/ \
                          derivation/proof/model roots; BudgetExhausted stays Unknown; retain countermodel/kernel/impact FailureBundle",
            g5_matrix: "model/proof shards {1,2,7,31} x workers {1,2,7} x graph/frontier orders {lex,reverse,permuted} x \
                        deterministic mode on identical kernel/toolchain; canonical graphs/queries, derivations, countermodels, AST/ \
                        premise/axiom closure, kernel result, witnesses, events and receipt roots match bitwise",
            entry_point: "scripts/e2e/leapfrog/i06_causal_impact.sh",
            tier: CampaignTier::Max,
            dsr_lane: "dsr quality --tool frankensim (i06-causal-impact isolated lane)",
            obs_events: &[
                "request.received",
                "causal.query_admitted",
                "causal.identified",
                "causal.unknown",
                "causal.countermodel",
                "impact.ast_checked",
                "impact.kernel_checked",
                "impact.completeness_unknown",
                "cancel.requested",
                "drain.completed",
                "finalize.completed",
                "failure_bundle.retained",
                "adjudication.receipt",
            ],
            replay_command: "scripts/e2e/leapfrog/i06_causal_impact.sh --manifest <manifest-id> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i06-robust-decision-certifier",
            claims_covered: &[
                "i06-robust-decision-optimality",
                "i06-false-provenance-impact-certificate-falsifier",
            ],
            unit_cases: UNIT_CASES,
            g0: "generators: finite candidate/action/test/scenario grammars, interval feasibility/utility, ambiguity/risk \
                 models, exact microcase optima, primal/dual bounds and hidden certificate mutants; predicates: grammar decode/ \
                 validity/completeness, evidence/ContextOfUse binding, hard feasibility, interval utility, lower<=feasible upper, \
                 normalized gap, rank/unrank/shard coverage, mutant witness; laws: action/scenario relabeling, exhaustive microcase \
                 agreement, adding feasible actions cannot worsen optimum, uncertainty widening cannot improve robust value, every \
                 known false certificate refused; shrink grammar/scenario/mutant while preserving witness",
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                "i06-material-family-matrix",
                "i06-decision-grammar-card",
                "i06-decision-certifier-mutants-max-holdout",
                "i06-governed-industrial-lot-pack",
            ],
            g3_relations: &[
                "candidate/action/scenario id permutation preserves optimal policy/value/gap up to canonical rebinding",
                "adding a feasible action can preserve or improve but cannot worsen global optimum",
                "widening ambiguity/evaluation intervals cannot improve a worst-case utility bound without proof",
                "hard-constraint, evidence-cutoff or ContextOfUse mutation moves authority before optimization",
                "shard/frontier order preserves completeness root, bounds and deterministic argmin tie break",
                "every accepted provenance/semantic/statistical/graph/causal/decision mutant refutes its certifier",
            ],
            g4_schedule: "whole-campaign preflight freezes grammar/evaluation budget before mutant reveal; cancel/panic/timeout/ \
                          corrupt checkpoint at action generation, feasibility, scenario evaluation, bound update, shard join, \
                          completeness and mutant adjudication; request->drain->finalize drains children and publishes only complete \
                          bound/terminal states; checkpoints bind frontier/visited/primal/dual/roots; BudgetExhausted stays Unknown; \
                          retain best feasible incumbent and every minimized false-certificate FailureBundle",
            g5_matrix: "action/scenario shards {1,2,7,31} x workers {1,2,7} x frontier orders {lex,reverse,permuted} x \
                        deterministic mode on identical toolchain; valid/excluded/evaluated counts, completeness roots, feasible \
                        incumbent, lower/upper/gap, optimal tie, mutant verdicts, events and receipts match bitwise",
            entry_point: "scripts/e2e/leapfrog/i06_robust_decision.sh",
            tier: CampaignTier::Max,
            dsr_lane: "dsr quality --tool frankensim (i06-robust-decision isolated lane)",
            obs_events: &[
                "request.received",
                "decision.action_generated",
                "decision.action_infeasible",
                "decision.incumbent_updated",
                "decision.bound_updated",
                "decision.completeness_checked",
                "certifier.mutant_refused",
                "certifier.refuted",
                "cancel.requested",
                "drain.completed",
                "finalize.completed",
                "failure_bundle.retained",
                "adjudication.receipt",
            ],
            replay_command: "scripts/e2e/leapfrog/i06_robust_decision.sh --manifest <manifest-id> --replay <artifact-id>",
        },
    ]
}

fn i06_waivers() -> Vec<Waiver> {
    vec![Waiver {
        subject: "i06-governed-industrial-lot-pack",
        reason: "all committed I06 fixtures are synthetic material-shaped data; independently authorized supplier COAs, \
                 raw metrology, custody attestations, process histories, blind lot labels, commercial property methods and \
                 real substitution outcomes are not yet licensed, privacy-cleared and content-addressed for replay",
        owner: "I06 material-data governance, laboratory V&V and supplier-registry owners",
        predicate: "a multi-supplier governed pack is admitted with exact consent/license/redaction, key and custody policy, \
                    raw authorized measurements/calibrations/methods, blinded partitions, ContextOfUse/substitution outcomes, \
                    checker ownership, retention and third-party replay constraints",
        expiry: "before the first I06 maximal decision campaign submission; review on every manifest amendment, supplier/ \
                 laboratory/license/key/spec/method revision and evidence-access change",
        promotion_effect: "synthetic schema/statistical/graph/theorem engineering may proceed, but no real supplier/lot \
                           qualification, industrial substitution performance, causally complete field impact, decision- \
                           optimal replacement or maximal I06 promotion may be claimed while this waiver is live; baseline \
                           synthetic authority remains separately eligible",
    }]
}
