//! The I03 (electrostatic/electroquasistatic physics) VerificationManifest
//! draft (bead frankensim-leapfrog-2026-program-i94v.2.1.8.1).
//!
//! Baseline lattice ([S]): exact cochain identities and independently measured
//! FEEC/CutFEM convergence, gauge and floating-conductor discipline,
//! capacitance extraction, locally conservative conduction, class-specific
//! dielectric admission, quantitative EQS routing, descriptor-circuit power,
//! and held-variable force/adjoint closure. Maximal lattice ([F]/[M]): partial
//! discharge and breakdown, space charge and aging, electrostriction, exact
//! variational naturality, certified refinement defects, topology-event jump
//! theorems, and an independent falsifier lane. A weaker receipt closes only
//! its own lattice element and is never relabeled as a stronger theorem.

use crate::{
    Ambition, CampaignTier, ClaimPolarity, ClaimSpec, FiveExplicits, FixturePin, FixtureSource,
    GauntletTier, ManifestDraft, ObligationRow, OracleRoute, Partition, ToleranceSemantics, Waiver,
};

const CAMPAIGN_POLICY_FIXTURE: &str = "i03-campaign-policy-v1";

/// Build the I03 draft. Consumers freeze it themselves; the conformance
/// battery proves the draft freezes and that its scientific seams remain
/// separately addressable.
#[must_use]
pub fn i03_draft() -> ManifestDraft {
    ManifestDraft {
        initiative: "I03",
        title: "Electrostatic/EQS physics gate: exact-sequence fields, conductors, \
                dielectrics, circuits, forces, discharge, aging, electrostriction, \
                and topology-aware force theorems",
        version: 1,
        explicits: FiveExplicits {
            units: "six-base quantity system L M T I Theta N; volts, amperes, \
                    coulombs, farads, siemens, joules, watts, pascals, and \
                    volts-per-metre are derived quantities; radians are \
                    dimensionless with semantic angle kind; every normalized QoI \
                    declares its dimensional denominator",
            seeds: "Philox 4x32-10 uses each deterministic fixture-declared alias 'i03/<stream>'; \
                    BLAKE3 derive-key domain org.frankensim.i03.fixture-stream.v1 \
                    hashes the exact UTF-8 alias and digest bytes 0..8 little-endian \
                    form the 64-bit Philox key; counter low64 is case index and \
                    counter high64 is output-block ordinal; development indices \
                    0..=4095, core held-out 65536..=69631, maximal held-out \
                    131072..=135167; paired development/holdout fixtures reuse an \
                    alias but are separated by disjoint inclusive ranges, while \
                    distinct aliases must have distinct derived keys. Statistical \
                    IID heldout authority is the separately governed raw-block \
                    commit/reveal protocol in i03-campaign-policy, conditional on \
                    its explicit at-least-one-honest uniform-source assumption and \
                    never a public Philox range or short-seed expansion",
            budgets: "smoke <= 90 s on one host; core <= 45 min and <= 24 GiB; \
                      max <= 12 h on a quiet campaign host and <= 64 GiB; each \
                      cancellation poll interval <= one declared operator tile; \
                      numerical accuracy is fixed by the per-claim tolerance",
            versions: "fs-vmanifest schema v2; fs-qty six-base \
                       schema; material/passivity-card schema v1; circuit-port \
                       ownership schema v1; theorem-card schema v1; toolchain pinned \
                       by rust-toolchain.toml and sibling revisions by \
                       constellation.lock",
            capabilities: "baseline: electrostatic, steady-conduction, \
                           eqs-regime-router, eqs-descriptor-circuit, \
                           electrostatic-force, \
                           exact-discrete-adjoint; maximal feature gates: \
                           partial-discharge, space-charge-aging, electrostriction, \
                           topology-force-theorem; no network or FFI; deterministic \
                           mode mandatory for G5",
        },
        claims: i03_claims(),
        fixtures: i03_fixtures(),
        obligations: i03_obligations(),
        waivers: i03_waivers(),
        amendment_rules: "After campaign start every semantic change forks a \
                          successor through FrozenManifest::amend. A changed \
                          hypothesis, material/regime pin, held-variable convention, \
                          fixture partition, tolerance, oracle, operational policy, or \
                          capability invalidates the affected descendants; no result \
                          may edit this version in place",
    }
}

#[allow(clippy::too_many_lines)]
fn i03_claims() -> Vec<ClaimSpec> {
    vec![
        ClaimSpec {
            id: "i03-electrostatic-exact-sequence",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "The admitted scalar-potential and mixed-field electrostatic \
                        complexes satisfy their integer cochain identities exactly: \
                        d1*d0=0, relative-boundary restrictions commute, gauge kernels \
                        are classified per unanchored connected component, and any mixed \
                        curl-free formulation explicitly represents H1 while true \
                        electrostatic closed-loop periods are zero; relative terminal \
                        voltages remain boundary/port data rather than loop periods",
            hypotheses: &[
                "the finite complex and relative boundary subcomplex are admitted and their incidence maps are exact integers",
                "scalar-potential gauges use relative H0: one constant mode per unanchored connected component, with pure-Neumann compatibility checked",
                "a direct or mixed field formulation explicitly carries an H1 basis, proves every physical electrostatic loop period is zero absent a separately modeled non-electrostatic EMF, and keeps relative terminal voltages in trace/port data",
            ],
            qoi: "exact_incidence_boundary_gauge_and_period_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G0,
            oracle: OracleRoute {
                identity: "i03.oracle.exact_sequence.v1 at fs-vmanifest-oracles/i03/exact_sequence.rs::check_incidence_and_gauge",
                independent: true,
                tcb_overlap: "shares canonical complex bytes only; independently multiplies integer incidence maps and computes component/cohomology ranks",
            },
            activation: "electrostatic complex, relative-boundary, gauge, and mixed-period admission are implemented",
            kill: "one nonzero exact composition, misclassified gauge mode, missing harmonic sector, nonzero electrostatic loop period without an admitted EMF, or rank mismatch kills this algebraic claim",
            fallback: "return Unsupported with the exact complex/boundary/gauge/period blocker; never regularize silently",
            no_claim: "exact cochain algebra alone says nothing about metric accuracy, stability, convergence, or full-wave validity",
        },
        ClaimSpec {
            id: "i03-electrostatic-convergence-gauss",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Under the pinned regularity, coercivity, commuting-projection, \
                        geometry, cut-quadrature, and asymptotic-family hypotheses, \
                        gauge-fixed FEEC/CutFEM electrostatics satisfies manufactured \
                        Gauss balance and the separately pinned energy-, flux-, and \
                        goal-functional convergence floors",
            hypotheses: &[
                "permittivity is uniformly bounded and positive on the solved quotient; boundary data satisfy compatibility and the admitted problem is coercive modulo its declared gauge",
                "the FEEC family has the pinned bounded commuting projection and the manufactured solution has the regularity required by each stated order",
                "body-fitted or certified-cut geometry, stabilization, quadrature, mesh sequence, norms, fit window, and pre-asymptotic rejection rule are frozen by the fixture",
            ],
            qoi: "maximum_preregistered_normalized_gauss_error_and_directed_order_shortfall",
            unit: "1",
            tolerance: ToleranceSemantics::Interval { lo: 0.0, hi: 1.0 },
            evidence_tier: GauntletTier::G1,
            oracle: OracleRoute {
                identity: "i03.oracle.electrostatic_mms.v1 at fs-vmanifest-oracles/i03/electrostatic_mms.rs::adjudicate",
                independent: true,
                tcb_overlap: "shares exact fixture equations and units; independently evaluates fields, fluxes, Gauss integrals, regression windows, and directed floors",
            },
            activation: "the exact-sequence claim is green and certified geometry/quadrature receipts exist",
            kill: "one normalized balance score above one, invalid convergence window, or missed directed floor kills this numerical claim",
            fallback: "return the admitted field with an Unknown accuracy receipt and the failed norm, level, cut class, or regularity hypothesis named",
            no_claim: "no rate outside the pinned regularity and asymptotic regime; a non-nested CutFEM goal error need not decrease at every refinement step",
        },
        ClaimSpec {
            id: "i03-floating-conductor-capacitance",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "With material state, free volume charge, and permanent \
                        polarization frozen, linear terminal charge is affine, \
                        q(V)=q0+C V. For reciprocal uniformly coercive permittivity, C \
                        is the Hessian of electrostatic energy. In a connected closed \
                        system containing every conductor including the enclosure, C is \
                        symmetric positive semidefinite with C*1=0 and ker(C)=span{1}; \
                        after a declared ground, infinity condition, or quotient \
                        projection, the corresponding coercive reduced operator is \
                        positive definite",
            hypotheses: &[
                "terminal charge orientation is fixed; the closed all-conductor q includes the enclosure and obeys 1^T q0+Q_vol=0, while an open/truncated q_terminal excludes the exterior boundary and separately obeys 1^T q_terminal+Q_vol-Q_exterior=0; these ledgers are distinct typed contracts",
                "C is d q/d V while free charge, permanent polarization, geometry, constitutive branch, and material state remain fixed; no state-changing voltage secant is relabeled as capacitance",
                "the closed theorem includes the enclosing conductor and a connected dielectric quotient; grounded, infinity-referenced, and explicitly projected matrices are separately typed contracts",
                "permittivity is reciprocal and uniformly positive on the admitted quotient; nonreciprocal cards and singular or disconnected quotients route to narrower predicates",
                "for every reciprocal closed-system case, an exact self-adjoint assembly receipt proves symmetry of the admitted bilinear form and adjointness of terminal lift/trace maps on the gauge quotient, hence C=C^T exactly; exact oriented columnwise Gauss incidence proves 1^T C=0 and exact common-shift incidence proves C*1=0 before any floating-point structural score is evaluated",
            ],
            qoi: "maximum_preregistered_normalized_charge_and_matrix_structure_defect",
            unit: "1",
            tolerance: ToleranceSemantics::Interval { lo: 0.0, hi: 1.0 },
            evidence_tier: GauntletTier::G2,
            oracle: OracleRoute {
                identity: "i03.oracle.capacitance.v1 at fs-vmanifest-oracles/i03/capacitance.rs::solve_and_check",
                independent: true,
                tcb_overlap: "shares canonical geometry, source, and material bytes; analytic coax/sphere, certified source-free boundary-integral, and independently equilibrated volume/source routes consume neither the production field, flux, nor matrix",
            },
            activation: "electrostatic exact-sequence and convergence claims are green at smoke tier",
            kill: "a source-changing column solve, missing conductor/reference, post-hoc symmetrization or row-sum projection, global-charge failure, wrong quotient/kernel dimension, gauge-dependent observables, or failed reciprocal/coercive structure blocks promotion",
            fallback: "return the admitted terminal solution and an Unknown capacitance receipt with the failed structural property named",
            no_claim: "a high-precision solve or extrapolation trend alone is not an enclosure; symmetry is not claimed for nonreciprocal laws, and zero row sums, a constant null mode, or strict positivity are asserted only for their typed closed or coercive-reduced contracts",
        },
        ClaimSpec {
            id: "i03-steady-conduction-conservation",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Steady conduction reconstructs current in an admitted H(div) or \
                        dual-cochain flux space, conserves charge cellwise and globally, \
                        assigns every terminal-current sign once, and reports \
                        non-negative Joule loss under the pinned passive-tensor contract",
            hypotheses: &[
                "the bounded conductivity tensor has non-negative symmetric part and is uniformly positive on the solved quotient whenever uniqueness is claimed",
                "a locally conservative H(div)/dual-cochain current or certified equilibrated-flux reconstruction is part of the result",
                "electrode/current boundary ownership and outward-normal convention are frozen",
                "sources and sinks carry compatible charge-per-time units; pure-Neumann data satisfy global compatibility and use an explicit gauge",
            ],
            qoi: "maximum_preregistered_normalized_local_global_balance_and_passivity_defect",
            unit: "1",
            tolerance: ToleranceSemantics::Interval { lo: 0.0, hi: 1.0 },
            evidence_tier: GauntletTier::G1,
            oracle: OracleRoute {
                identity: "i03.oracle.conduction.v1 at fs-vmanifest-oracles/i03/conduction.rs::check_flux_balance",
                independent: true,
                tcb_overlap: "shares material and boundary bytes; independently reconstructs face fluxes and evaluates analytic network reductions",
            },
            activation: "steady-conduction constitutive and port schemas are admitted",
            kill: "a missing current owner, absent conservative flux receipt, negative passive loss beyond its normalized budget, or local/global balance score above one kills the claim",
            fallback: "refuse decision-grade terminal currents and retain the minimized conservation defect",
            no_claim: "a primal scalar FE gradient alone is not claimed cellwise conservative; no electrochemical reaction, ballistic carrier transport, or unpinned contact resistance",
        },
        ClaimSpec {
            id: "i03-declared-dielectric-laws",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Every linear, nonlinear, anisotropic, and dispersive dielectric \
                        law is routed by a versioned material-class card: LTI rational \
                        laws prove proper pole stability, realness/causality, and the \
                        pinned matrix-valued positive-real supply convention including \
                        feedthrough; nonlinear/history laws bind their evolution \
                        equation or maximal-monotone inclusion to a storage/supply/ \
                        dissipation inequality, or supply an explicitly scoped \
                        incremental-passivity theorem; empirical cards receive no \
                        thermodynamic authority beyond their validation domain",
            hypotheses: &[
                "constitutive parameters, transform-sign convention, temperature/frequency domain, internal state, and history initialization are explicit",
                "the card declares LTI-rational, nonlinear-storage, incremental/IQC, or empirical class and routes to its class-specific independent checker",
                "passivity and reciprocity predicates are properties of the pinned card version, electrical effort/flow convention, boundary/history terms, state ensemble, and evolution law, not inferred from convex functions or a type name",
            ],
            qoi: "material_regime_admission_and_passivity_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G0,
            oracle: OracleRoute {
                identity: "i03.oracle.dielectric_cards.v1 at fs-vmanifest-oracles/i03/dielectric.rs::dispatch_class_checker",
                independent: true,
                tcb_overlap: "shares canonical card bytes/parser; independent rational, storage-inequality, IQC, and domain checkers do not share production update kernels",
            },
            activation: "the dielectric-card registry and history-state schema are implemented",
            kill: "one class/checker mismatch, unsupported state admitted as Supported, failed storage/dissipation inequality, or reciprocity claim without its hypothesis returns the card family to review",
            fallback: "route to a narrower admitted law or return Unknown with the missing regime/state receipt",
            no_claim: "positive-realness of an unspecified transfer function is not a nonlinear-passivity proof; no universal dielectric model, gyrotropic reciprocity, or history-free dispersive solve",
        },
        ClaimSpec {
            id: "i03-eqs-regime-routing",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "A quantitative electroquasistatic regime card admits the \
                        curl-free approximation only through one of two certified routes: \
                        a separable complete modal expansion with an enclosed excitation \
                        tail, or a lifted full-Maxwell residual paired with certified \
                        stability and dual-weighted QoI remainder bounds. Propagation, \
                        induction, omitted magnetic energy, radiation, and requested-QoI \
                        monitors must meet their pinned budgets; skin-depth and charge- \
                        relaxation monitors classify only their validated conductive/ \
                        displacement subcards. Terminal current always includes J_cond+ \
                        dD/dt, and a failed or indeterminate certificate escalates or \
                        returns Unknown/Unsupported with no fabricated EQS result",
            hypotheses: &[
                "geometry, stationary or moving port class, waveform spectrum, passive possibly anisotropic/dispersive epsilon/mu/sigma bounds, norms, requested QoI, and positive error-budget allocation are explicit",
                "the modal route is used only where a branch-certified separable Maxwell pencil, complete excited-mode set, and discarded-tail enclosure exist; general heterogeneous anisotropic geometry uses the lifted operator-residual route with an enclosed residual-norm upper bound, a strictly positive enclosed inf-sup lower bound, the resulting rigorous field-error quotient, and a dual-weighted QoI bound",
                "the regime card pins exact norms/projections and candidate-independent positive scale formulas for propagation, induction-to-conservative-field, omitted magnetic-energy, gross-input radiation, and dimensionless direct QoI-remainder upper bounds plus a boundary band around every threshold; signed net work is never a denominator, and skin-depth/charge-relaxation quantities only select a validated static/ohmic/displacement subcard",
                "the independent analytic or separately implemented Maxwell solution adjudicates the bounds but never supplies a missing runtime certificate; stationary terminal current integrates J_cond+dD/dt with one orientation, while moving ports use the separate pullback/GCL contract",
            ],
            qoi: "eqs_admit_escalate_unknown_and_total_current_acceptance_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G2,
            oracle: OracleRoute {
                identity: "i03.oracle.eqs_regime.v1 at fs-vmanifest-oracles/i03/eqs_regime.rs::adjudicate",
                independent: true,
                tcb_overlap: "shares canonical geometry/material/waveform bytes; analytic and reference-Maxwell monitors do not share the production EQS operator",
            },
            activation: "electrostatic convergence, dielectric admission, and total-current port ownership are green",
            kill: "using a bulk scalar propagation constant for an unseparated heterogeneous problem, omitting modal completeness/tail or residual stability/dual bounds, an undefined norm/projection, a signed or zero energy denominator, a unit-bearing unnormalized QoI comparison, admitting from a sampled reference error, one missing displacement-current contribution, or one silent failure to escalate kills the router claim",
            fallback: "invoke the separately admitted Maxwell capability when available; otherwise return Unknown or Unsupported with the violated monitor and required capability",
            no_claim: "small frequency, large skin depth, or short charge-relaxation time alone is not an EQS theorem; admission is fixture-, norm-, port-, and QoI-specific and does not claim magnetic-force or radiation accuracy",
        },
        ClaimSpec {
            id: "i03-field-circuit-power-closure",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Electroquasistatic field/descriptor-circuit coupling owns each \
                        terminal voltage, total current, and orientation exactly once. \
                        With W_fc positive into the field, three separate identities close: \
                        W_fc+W_volume,in=Delta(W_field+W_material)+D_field+W_mech,out+ \
                        D_num,field; W_external-W_fc=Delta(W_circuit)+D_circuit+ \
                        D_num,circuit; and their sum W_external+W_volume,in= \
                        Delta(W_field+W_material+W_circuit)+D_field+D_circuit+ \
                        W_mech,out+D_num,field+D_num,circuit, where internal W_fc cancels",
            hypotheses: &[
                "descriptor circuit is admitted with an explicit gauge and consistent initial state",
                "total terminal current is the primary oriented port flow. With dielectric-outward n and current positive into the field, a fixed port satisfies I_total=-integral_Gamma(J_f+partial_t D) dot n=dot(Q_e_free)-dot(Q_e_transfer), where Q_e_free=-integral_Gamma D dot n and Q_e_transfer=integral_0^t integral_Gamma J_f dot n are disjoint owners; I_total=dot(Q_e_free) only for a blocking port with J_f dot n=0. Moving ports use the corresponding pinned ALE/Reynolds pullback of relative carrier flux plus D, swept/geometric terms, and the discrete geometric-conservation law; every field, circuit, external/volume source, mechanical, storage, and dissipation term has exactly one owner in each declared control volume",
                "exact equality is limited to cards with a verified storage/dissipation realization and the same pinned discrete time/work quadrature; every D_num term is a preregistered method-derived quantity, never a fitted closure slack; other passive cards receive only their declared inequality",
            ],
            qoi: "maximum_preregistered_normalized_total_current_and_discrete_work_balance_defect",
            unit: "1",
            tolerance: ToleranceSemantics::Interval { lo: 0.0, hi: 1.0 },
            evidence_tier: GauntletTier::G1,
            oracle: OracleRoute {
                identity: "i03.oracle.field_circuit_power.v1 at fs-vmanifest-oracles/i03/field_circuit_power.rs::integrate_and_balance",
                independent: true,
                tcb_overlap: "shares port schema, card-defined storage functions, and canonical quadrature declaration; uses an independent high-precision small-DAE integrator and ledger evaluator",
            },
            activation: "EQS routing, capacitance/conduction, storage-realization, and descriptor-circuit port ownership are green",
            kill: "duplicate/missing term ownership, mixing subsystem and coupled ledgers, source-sign ambiguity, fitted numerical slack, inconsistent DAE state, an equality asserted for an inequality-only card, or normalized closure score above one kills the coupled claim",
            fallback: "uncoupled field solve with explicit prescribed terminals; coupled result remains Unknown",
            no_claim: "no equality for undeclared material memory, missing magnetic/circuit storage, or moving geometry without a certified ALE/geometric-conservation receipt; no full-wave/radiation claim",
        },
        ClaimSpec {
            id: "i03-force-adjoint-held-variable-closure",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Electrostatic generalized force, virtual work, and exact- \
                        discrete adjoints agree under the frozen constrained \
                        thermodynamic ensemble. An admitted scalar-potential card gives \
                        the full electrical coenergy W*(x,V,z,history) with Q=d_V W*; \
                        its strictly convex charge-side Legendre dual governs fixed Q. \
                        Mixed/floating, nonintegrable, nonreciprocal, and hysteretic \
                        cards instead differentiate their complete constrained field/ \
                        circuit action without inventing a scalar coenergy. The familiar \
                        quadratic C(x) formulas are corollaries only for reciprocal \
                        linear zero-offset Q=C(x)V cards; derivatives freeze topology, \
                        state/history branch, and event mode",
            hypotheses: &[
                "the equilibrium shape path is smooth, orientation-preserving, and remains on one regular/stable branch with fixed topology, cut classification, terminal ownership, active material state, and event mode",
                "a scalar-coenergy card proves that the terminal-charge one-form is closed and exact on its declared same-branch voltage domain; the radial formula is allowed only when 0 belongs to that domain and every segment {lambda V | 0<=lambda<=1} stays inside it, while a merely simply connected domain uses a separately content-pinned in-domain path after path independence is proved",
                "for a closed all-conductor system the fixed-Q route uses Vbar=R^m/span{1}, Qbar={Q | 1^T Q=-Q_vol}, qhat=Q-q0 in 1-perp, and the gauge-invariant reduced coenergy Wbar*([V])=W*(V)-<V,q0>; its Legendre-Fenchel dual pairs [V] only with qhat, proves gauge-section independence, and requires strict convexity/invertibility only on Vbar; grounded/open contracts use their own coercive spaces and incompatible Q or an unproved affine offset/source reconstruction is Unsupported",
                "the complete fixed-V, fixed-Q, mixed, or floating constraint ensemble, charge offset, internal/history state, source terms, and corresponding signed coenergy/Legendre/constrained action are explicit fixture fields; nonintegrable cards route directly to the action formulation",
                "interface traction is derived from the same total electrical/material energy; complex-step is used only when the entire composed geometry, solve, state update, and QoI path is holomorphic, otherwise an independent interval-controlled symmetric/Ridders difference is used",
            ],
            qoi: "maximum_preregistered_normalized_generalized_force_covector_virtual_work_and_adjoint_discrepancy",
            unit: "1",
            tolerance: ToleranceSemantics::Interval { lo: 0.0, hi: 1.0 },
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "i03.oracle.force_adjoint.v1 at fs-vmanifest-oracles/i03/force_adjoint.rs::check_ensemble_derivative",
                independent: true,
                tcb_overlap: "shares frozen energy values, ensemble card, and geometry path; derivative stencil, truncation/roundoff enclosure, and constrained functional assembly are independent",
            },
            activation: "field/circuit power closure and shape-path identity are green",
            kill: "a missing integrability/Legendre certificate, incompatible affine charge, gauge-dependent Legendre value, sign/ensemble ambiguity, invented scalar coenergy, non-holomorphic complex-step use, branch leakage, or normalized adjoint/virtual-work score above one blocks force and optimization promotion",
            fallback: "report force/gradient Unknown and use derivative-free bounded design exploration",
            no_claim: "no scalar coenergy for a nonintegrable terminal-charge one-form and no differentiability through breakdown, contact/event switching, material branch changes, or topology changes",
        },
        ClaimSpec {
            id: "i03-partial-discharge-breakdown-routing",
            ambition: Ambition::Frontier,
            polarity: ClaimPolarity::Affirmative,
            statement: "A certified regime router separates field-only stress, partial- \
                        discharge inception, and breakdown lanes. On the synthetic \
                        held-out family, simultaneous familywise-99% exact one-sided \
                        lower confidence bounds over independent material-lot clusters \
                        must each be at least 0.90 for every preregistered fixture-family/ \
                        outcome cell while cellwise median and 90th-percentile normalized \
                        interval widths meet their pinned sharpness caps; censoring is \
                        scored conservatively and out-of-domain \
                        inputs become Unknown/Unsupported, never fabricated predictions",
            hypotheses: &[
                "insulation geometry/defects, gas pressure/composition, waveform and measurement bandwidth, apparent-charge/event threshold, temperature, humidity, and material-lot distribution are pinned",
                "one preregistered outcome per fixture-family and endpoint is formed for each material lot; within-lot dependence and censoring are arbitrary but retained, while exact Clopper-Pearson authority is conditional on an independently governed audited receipt that the pinned finite material-lot law has exact probabilities n_i/2^256 and that 512 lot blocks plus 512 disjoint candidate-seed blocks are jointly IID uniform 256-bit values",
                "the frozen custodian roster, verification-key digests, commitment encoding/domain, and 3-of-3 reveal rule are admission-bound, but candidate/model/toolchain and every population, censoring, normalization, and stopping semantic are irrevocably committed first. Only afterward do the custodians generate and commit exactly 1024 indexed 256-bit blocks; exact authority assumes and independently audits that at least one custodian generated its complete block vector information-theoretically uniform and independent of every other source and the already-fixed candidate. The commitment transcript stays outside the candidate sandbox; sorted reveal blocks are combined only by frozen coordinate-wise XOR, never by deterministic short-seed expansion. A separately bound future 256-bit beacon challenges evaluation order but is not an entropy input to any lot/seed block, so its distribution or commitment dependence cannot perturb IID sampling; a missing/late/duplicate reveal, beacon substitution, abort, retry, or resampling is IntegrityFailed rather than an omitted input or a new holdout",
                "the governed evaluator applies one identical measurable candidate map to every lot. Its complete visible input is exactly that lot's features derived from U_(LOT,l), its disjoint IID U_(CANDIDATE_SEED,l), and frozen campaign constants; lot ordinal, stable case/receipt id, processing order, shard, worker, nonce, beacon, other blocks, and logging metadata remain outside the candidate sandbox. Label/prediction/feature feedback, online updates, cache or RNG carryover, and cross-lot stopping are forbidden, and both secret order permutation and fresh random reassignment of every external receipt id must leave each lot's prediction byte-identical before identically distributed outcome/Clopper-Pearson authority exists",
                "coverage and sharpness use the exact compatible-set containment and full-domain normalized-width formulas and are evaluated once on the untouched maximal holdout; synthetic generator conformance is not experimental material validation",
            ],
            qoi: "coverage_confidence_sharpness_and_regime_routing_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G2,
            oracle: OracleRoute {
                identity: "i03.oracle.discharge_coverage.v1 at fs-vmanifest-oracles/i03/discharge_coverage.rs::adjudicate_clusters",
                independent: true,
                tcb_overlap: "shares sealed labels, normalization, and censoring declaration; independently computes cluster outcomes, exact confidence bounds, width statistics, and regime decisions",
            },
            activation: "all baseline field/material claims are green and the partial-discharge feature is enabled",
            kill: "coverage lower-bound or sharpness failure, missing/compromised IID sampling receipt or commit/reveal access control, cross-lot state/adaptation/order dependence, holdout leakage, silent extrapolation, or Supported outside the declared regime kills promotion of the model family",
            fallback: "field-stress envelope plus Unknown discharge/breakdown disposition and ranked missing evidence",
            no_claim: "synthetic holdout evidence validates algorithmic calibration only; independent experimental records are required for physical/industrial promotion, and no universal inception, deterministic breakdown, or service-life law is claimed",
        },
        ClaimSpec {
            id: "i03-space-charge-aging-singular-routing",
            ambition: Ambition::Frontier,
            polarity: ClaimPolarity::Affirmative,
            statement: "Space-charge transport conserves a typed augmented charge ledger \
                        across dimensionally distinct bulk densities, interface/surface \
                        densities, lumped electrode charges, and external port transfers. \
                        Dimensionally homogeneous closed bulk reactions satisfy \
                        q_bulk^T S_bulk=0; free bulk charge obeys continuity, \
                        polarization supplies distributional bound charge/current exactly \
                        once, carrier transfer is separated from electrode D-trace free \
                        charge and external total current, and the integrated control- \
                        volume ledger reconciles all owners. Aging states enforce \
                        positivity, capacity, and damage bounds, while every conductor/ \
                        insulator/high-contrast limit has a pinned scaling, limiting \
                        problem, topology, uniform estimate, and refusal boundary",
            hypotheses: &[
                "mobility, trapping, injection, recombination, electrode transfer, thermal, and damage stoichiometry are versioned material cards, and every bulk, surface, electrode, or port reservoir has its own unit, support dimension, sign, capacity, ownership, initial state, and volume, surface, lumped, or boundary evolution law",
                "q_bulk and S_bulk act only on dimensionally homogeneous bulk concentration columns; surface/electrode/global charges enter after the pinned measure pairing through a typed control-volume incidence ledger, and any reduced non-charge-null bulk reaction has an explicit compensating surface/electrode/port transfer",
                "the free-charge/displacement and total-charge/polarization formulations are related by a pinned commuting-incidence equivalence proof with rho_b=-div P in the distributional sense; for a fixed interface with n from minus to plus, sigma_b=n dot (P_minus-P_plus) and dot(sigma_b)+n dot (J_b_plus-J_b_minus)=0 for J_b=partial_t P, while moving interfaces use a separately certified spacetime transport/geometric term",
                "with dielectric-outward electrode normal and external total current positive into the field, each fixed-port ledger proves I_total=-integral_Gamma(J_f+partial_t D) dot n=dot(Q_e_free)-dot(Q_e_transfer); carrier transfer, D-trace free charge, and total current remain separate owners, and moving ports use the correspondingly pulled-back relative-carrier plus displacement/geometric/GCL identity",
                "energy equality is asserted only for thermodynamically generated cards; empirical aging cards carry validation-domain and service-life no-claim receipts",
                "each singular family freezes its nondimensional scaling, limiting DAE/PDE and boundary conditions, norm/topology, uniform error estimate, and index-change escalation",
            ],
            qoi: "maximum_preregistered_normalized_charge_state_energy_and_limit_route_defect",
            unit: "1",
            tolerance: ToleranceSemantics::Interval { lo: 0.0, hi: 1.0 },
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "i03.oracle.space_charge_aging.v1 at fs-vmanifest-oracles/i03/space_charge_aging.rs::check_ledger_and_limits",
                independent: true,
                tcb_overlap: "shares material/stoichiometry bytes; uses an independent conservative finite-volume reference and analytic limiting problems",
            },
            activation: "baseline conduction/dielectric claims are green and the independent space-charge-aging feature is enabled",
            kill: "a dimensionally mixed reaction vector, omitted charged reservoir, q_bulk^T S_bulk failure without a compensating owner, double-counted or sign-wrong polarization surface term, missing moving-interface transport, electrode carrier/free-charge/total-current mismatch, state/capacity violation, thermodynamic overclaim, wrong limiting problem, or finite answer after a required refusal kills the lane",
            fallback: "bounded static field envelope with frozen material state and explicit aging no-claim",
            no_claim: "no service-life extrapolation beyond pinned horizons/lots, no energy authority for empirical damage cards, and no convergence claim for an unclassified reservoir, moving-boundary flux, DAE-index transition, or singular limit",
        },
        ClaimSpec {
            id: "i03-electrostriction-energy-interface-closure",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Affirmative,
            statement: "For the pinned finite-strain electrostrictive energy class and \
                        electrical ensemble, numerical bulk stress, total Piola/Cauchy \
                        and configurational-interface traction, virtual work, and exact- \
                        discrete adjoints close as derivatives of one total potential on \
                        regular equilibrium branches and admitted singular-limit families",
            hypotheses: &[
                "the objective free-energy density is expressed through C=F^T F, referential electric variables and reference structural tensors (or an equivalent frame-indifferent invariant basis); voltage/charge ensemble, Legendre transform, pullback/pushforward, and interface normal are fixed by the material card",
                "the deformation selects one exact invertibility route: a C1 global-homeomorphism theorem with positive Jacobian, proper degree-one local diffeomorphism, and certified boundary embedding, or a W^{1,p}, p>3, Ciarlet-Nečas a.e.-injectivity theorem with bounded Lipschitz domain, continuous/Lusin-N representative, positive determinant a.e., area formula, trace/degree premises, and independently enclosed image measure satisfying integral det(F)<=measure(y(Omega)); the receipt states whether global or only a.e. noninterpenetration was proved",
                "the equilibrium branch carries either an exact reduced-field-elimination receipt with a field-block isomorphism and coercive reduced Schur Hessian on the mechanical constraint tangent, or an exact mixed bordered/KKT-Hessian isomorphism with uniform inf-sup constant and fixed Morse index; coercivity of the full fixed-voltage saddle Hessian is never inferred",
                "topology and constitutive branch remain fixed; a singular-limit lane supplies a uniform lower bound on the selected reduced-coercivity or mixed inf-sup constant, compactness, convergence topology, and a bound justifying interchange of differentiation and limit",
            ],
            qoi: "maximum_preregistered_normalized_energy_stress_traction_virtual_work_and_adjoint_defect",
            unit: "1",
            tolerance: ToleranceSemantics::Interval { lo: 0.0, hi: 1.0 },
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "i03.oracle.electrostriction_numeric.v1 at fs-vmanifest-oracles/i03/electrostriction_numeric.rs::enclose_closure",
                independent: true,
                tcb_overlap: "shares exact energy/theorem-card bytes; independent symbolic derivatives, interval virtual work, and interface-balance evaluation do not share production kernels",
            },
            activation: "baseline force/dielectric closure is green and the electrostriction feature plus numerical theorem-card domain are enabled",
            kill: "normalized closure above one, loss of injectivity/stability, traction/sign mismatch, or unclassified singular route kills numerical promotion",
            fallback: "one-way electrostatic-to-mechanical loading with electrostriction and theorem authority disabled",
            no_claim: "a.e. injectivity is not relabeled as a global homeomorphism; no blanket Maxwell-stress equivalence or singular-limit interchange outside the pinned energy, ensemble, regularity, interface, topology, and theorem domain",
        },
        ClaimSpec {
            id: "i03-electrostriction-interface-theorem",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Affirmative,
            statement: "A machine-checked theorem card derives the total bulk and \
                        configurational-interface electrostrictive stress, virtual-work \
                        identity, held-ensemble signs, and exact-discrete adjoint from one \
                        pinned free energy, and proves its explicitly conditioned \
                        singular-limit corollary when differentiation commutes with limit",
            hypotheses: &[
                "the exact theorem-card bytes fix free energy, independent electrical/mechanical state, ensemble Legendre transform, frames, traces, normals, and interface jump convention",
                "the declaration selects and binds either the smooth proper-degree-one global-homeomorphism premises or the precise W^{1,p} Ciarlet-Nečas a.e.-injectivity premises, including domain, representative, Lusin-N/area-formula, determinant, trace/degree, and image-measure hypotheses; branch regularity/stability and any uniform singular-limit hypotheses are explicit",
                "equilibrium stability selects either reduced electrical/internal elimination with an exact coercive mechanical Schur Hessian, or a mixed bordered/KKT-Hessian isomorphism with certified inf-sup constant and fixed Morse index; the full fixed-voltage saddle Hessian is never silently called coercive, and singular-limit bounds are uniform in the selected route",
                "a canonical binding receipt proves byte/semantic equivalence among the manifest claim, serialized theorem card, generated Lean proposition, elaborated checked proposition, and exported declaration; the pinned independent kernel verifies the proof term, complete transitive declaration/environment axiom closure, and exact i03.lean-axioms.v1 policy whose only permitted foundational names are propext, Quot.sound, and Classical.choice before any theorem color is consumed",
            ],
            qoi: "independent_electrostriction_theorem_checker_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "i03.oracle.electrostriction_lean.v1 at proofs/i03/Electrostriction.lean::electrostrictionInterfaceClosure checked by pinned Lean4 kernel receipt",
                independent: true,
                tcb_overlap: "shares theorem-card serialization and rational constitutive expressions only; proof kernel and numerical production path are disjoint",
            },
            activation: "the target theorem card is frozen, and a pre-candidate manifest successor has frozen the complete machine theorem AST, symbol definitions, runtime-premise schema, and total AST-to-Lean translation required by FORMAL_PROJECTION; theorem checking is independent of partial-discharge and space-charge campaign outcomes",
            kill: "kernel or declaration-binding rejection, sorryAx, a custom postulate/theorem-equivalent axiom, native-oracle proof authority outside the admitted kernel/checker TCB, any transitive axiom outside exactly {propext, Quot.sound, Classical.choice}, or one independently verified premise-satisfying counterexample refutes this manifest version's theorem",
            fallback: "retain numerical closure receipts without formal theorem authority",
            no_claim: "checker acceptance proves only the exact frozen theorem; it does not establish that a runtime case satisfies the premises or that an unpinned singular limit is valid",
        },
        ClaimSpec {
            id: "i03-cohomology-force-naturality-theorem",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Affirmative,
            statement: "A machine-checked variational-sheaf stationary-naturality theorem \
                        derives, rather than assumes, an equivalence of gauge-reduced \
                        equilibrium solution groupoids and exact pullback covariance of the \
                        resulting generalized-force one-form. Representations may use \
                        different covers and complexes when coherent descent and a \
                        differentiable restriction-compatible chain equivalence are lifted \
                        either to a certified equivariant stationary condensation with a \
                        unique smooth vertical stationary lift and exact reduced-action/ \
                        Schur identity, or to complete filtered pronilpotent cyclic \
                        L_infinity/BV data with convergent homotopy inverse that proves \
                        equivalence of the selected derived critical groupoids. If phi maps \
                        configuration charts, the conclusion is f_A=phi^*f_B; a vector- \
                        force statement follows only when the pinned configuration metrics \
                        and Riesz maps are themselves natural",
            hypotheses: &[
                "the representations are pinned variational sheaves on finite relative-complex covers; restriction maps, overlap cocycles, Cech descent, chain maps P_x and Q_x, and homotopies Q_x P_x~id and P_x Q_x~id are certified and differentiable on the declared configuration domain",
                "the stationary-condensation route pins a restriction-compatible vertical/horizontal splitting, contractible vertical fibers, a unique smooth equivariant stationary lift, uniformly invertible or coercive vertical Hessian, and exact reduced-action plus Schur-complement identity; the cyclic route instead uses complete descending filtered pronilpotent cyclic L_infinity/BV data, filtration-preserving maps and homotopies, a nondegenerate pairing and master equation, termination or certified convergence of every higher series, and a filtration-preserving homotopy inverse proving equivalence of the card-selected Deligne-Getzler or Hinich derived critical groupoids; a finite-polynomial alternative must supply explicit mutually inverse functors and checked natural equivalences rather than invoke a generic quasi-isomorphism slogan",
                "after gauge reduction the complete scalar ensemble actions agree through the selected variational-transfer route under phi and P_x, including geometry, weighted Hodge/material energy and state, loads, traces, ports, source/coenergy terms, constraints, and perturbation transport; the proof retains total derivatives of phi, P_x, the stationary lift or higher maps, restriction maps, and every shape-dependent primitive",
                "gauge and constraint quotients are regular and stationary classes exist on the compared branch; ordinary chain-homotopy equivalence without one of the two variational-transfer receipts has no equilibrium authority",
                "cohomology bases and coordinates transform contragrediently; compensated potential/exact-component changes preserve the assembled physical cochain, and harmonic Gram/Hodge derivatives are retained",
                "solution correspondence and force/QoI equality are conclusions, never premises; an additional QoI requires its own naturality square or appears as an explicit transfer defect",
                "the canonical theorem-binding receipt matches the manifest claim bytes, theorem-card AST, fully qualified declaration, elaborated type, environment, proof term, and complete transitive axiom closure against exact policy i03.lean-axioms.v1={propext, Quot.sound, Classical.choice} before kernel success is consumable",
            ],
            qoi: "independent_variational_sheaf_force_pullback_checker_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "i03.oracle.topology_force_lean.v1 at proofs/i03/TopologyForce.lean::representationNaturalityForcePullback checked by pinned Lean4 kernel receipt",
                independent: true,
                tcb_overlap: "shares the bound theorem-card AST, incidence integers, rational variational data, and commuting/homotopy diagrams only; proof kernel and production force path are disjoint",
            },
            activation: "baseline electrostatic/force closure is green, the topology-force target card is frozen, and a pre-candidate manifest successor has frozen the complete naturality machine AST, symbol definitions, runtime-premise schema, and total AST-to-Lean translation required by FORMAL_PROJECTION; electrostriction, discharge, and aging outcomes are not prerequisites",
            kill: "checker or declaration-binding rejection, a failed descent/homotopy/action-naturality premise, a failed stationary-lift/Hessian/Schur receipt, a failed cyclic-L_infinity/BV critical-groupoid receipt, or one independently verified in-domain counterexample refutes exactly the bound theorem declaration",
            fallback: "retain per-chart numerical force receipts with no topology-invariance transfer",
            no_claim: "chain-homotopy or cohomology equivalence alone preserves neither stationary points nor force; a vector representative is not natural without a natural metric/Riesz map; missing variational-transfer, QoI, action, descent, material, load, trace, port, or source naturality routes to the defect/event theorem rather than being silently inferred",
        },
        ClaimSpec {
            id: "i03-refinement-force-defect-enclosure",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Affirmative,
            statement: "For an admitted non-isomorphic representation transfer or \
                        refinement, a proof-carrying averaged-adjoint identity expresses \
                        the pulled-back target generalized-force one-form minus the source \
                        one-form as a signed algebraic sum of primal/dual residual, \
                        ensemble-transfer, QoI-transfer, chain-homotopy/descent, Hodge/material, \
                        geometry, load, trace, port, gauge/constraint, quadrature, solver, \
                        state/perturbation-transfer, endpoint, and nonlinear path- \
                        remainder contributions. A dependency-preserving joint enclosure \
                        retains signed cancellation. Componentwise zero is sufficient for \
                        exact covariance, but a nontrivial formally certified cancellation \
                        of the total sum may also prove exact covariance",
            hypotheses: &[
                "the source and target variational sheaves/actions, configuration map, prolongation/restriction maps, chain homotopies, descent data, source and target QoI/ensemble definitions, states, averaged path adjoint, loads, traces, ports, gauges, constraints, and perturbation fields are content-pinned",
                "a smooth or explicitly piecewise-smooth common-path construction is admitted; every path stratum, endpoint term, noncommuting square, QoI mismatch, and homotopy remainder is named, while a genuine topology/event crossing routes to the event theorem",
                "stability, regularity, geometry/quadrature consistency, and every saturation assumption used by the enclosure are independently checked and named in the receipt",
                "the same physical held-variable protocol, material branch, relative cohomology sector, and intended observable are compared; representational differences are charged to their transfer defects",
                "the signed-identity declaration, the exact configuration transfer T_x whose cotangent pullback compares force one-forms, and every generated defect label are bound byte/semantically to the theorem card and checked proposition before an enclosure receives theorem authority",
                "the Lean identity receipt binds the complete transitive axiom closure and accepts exactly policy i03.lean-axioms.v1={propext, Quot.sound, Classical.choice}; sorryAx, custom postulates, theorem-equivalent axioms, and native-oracle shortcuts outside the admitted kernel/checker TCB are IntegrityFailed",
            ],
            qoi: "normalized_pulled_generalized_force_outside_dependency_preserving_signed_defect_enclosure",
            unit: "1",
            tolerance: ToleranceSemantics::Interval { lo: 0.0, hi: 1.0 },
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "i03.oracle.refinement_force_defect.v2 composite proofs/i03/TopologyForce.lean::refinementForceDefectIdentity plus fs-vmanifest-oracles/i03/refinement_force.rs::enclose_signed_pullback_defect",
                independent: true,
                tcb_overlap: "shares the bound theorem-card AST, frozen states, actions, QoIs, maps, and homotopies; the Lean kernel checks the signed identity while an independent numerical adjudicator recomputes residuals, preserves interval dependencies, and evaluates high-precision dual-work differences",
            },
            activation: "baseline force closure and the source/target representation certificates are green, and a pre-candidate manifest successor has frozen the complete signed-defect machine AST, symbol definitions, runtime-premise schema, and total AST-to-Lean translation required by FORMAL_PROJECTION",
            kill: "kernel/declaration-binding rejection, an observed pulled-back force difference outside the joint signed enclosure, an unnamed ensemble/QoI/homotopy/noncommuting contribution, double-counting, or certified exact covariance without a proof that the total signed defect is zero refutes this manifest version's refinement theorem",
            fallback: "report chart-specific forces and the unresolved defect budget without transfer authority",
            no_claim: "componentwise nonzero defects do not refute exact covariance because they may cancel; no stepwise monotonic force or goal-error decrease is claimed, and convergence requires the pinned stability, regularity, consistency, transfer, and asymptotic hypotheses",
        },
        ClaimSpec {
            id: "i03-topology-event-jump-theorem",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Affirmative,
            statement: "For a pinned regularized topology-changing family on a \
                        stratified spacetime/configuration cobordism, a machine-checked \
                        weak-limit theorem passes the exact finite-regularization cochain \
                        charge and energy-work balances to a declared distribution or \
                        Radon-measure topology. Bulk, boundary, interface, field, material, \
                        circuit, source, contact/Joule, heat/radiation, and mechanical \
                        terms have disjoint owners and every internal transfer cancels \
                        once by orientation. A finite force/work atom or mechanical \
                        impulse is concluded only under its pinned tightness, \
                        integrability, parameter, and renormalization premises",
            hypotheses: &[
                "the event card fixes pre/post complexes and port spaces, stratified cobordism, event parameter and units, regularization/scaling family, test-function space, limiting topology, charge-transfer and circuit/source maps, material-state map, and pre/post trace ownership",
                "finite-regularization models satisfy exact typed charge and energy-work ledgers; uniform bounds, tightness/compactness, and convergence of every owned term are proved or the term is explicitly DeclaredDivergent, and a certified epsilon-to-zero enclosure selects a unique limit",
                "every subtraction or finite-part renormalization has a content-pinned counterterm owner and theorem stating scheme dependence or independence; voltage-current work is defined as a limit of finite pairings, never as an undefined product of distributions",
                "mechanical impulse is named only when the event parameter is physical time and momentum balance is modeled; otherwise the atomic conclusion is generalized work or a force measure",
                "the runtime event satisfies the checker-visible premises and crosses no undeclared contact, breakdown, source, or constitutive branch",
                "the manifest event claim, event-card AST, fully qualified declaration, elaborated type, ownership schema, proof term, and complete transitive axiom closure are linked by the canonical binding receipt and checked against exact policy i03.lean-axioms.v1={propext, Quot.sound, Classical.choice}",
            ],
            qoi: "independent_topology_event_jump_checker_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "i03.oracle.topology_event_lean.v1 at proofs/i03/TopologyEventJump.lean::distributionalEventBalance checked by pinned Lean4 kernel receipt",
                independent: true,
                tcb_overlap: "shares the bound event-card AST, cobordism, typed ownership partition, and rational weak-balance terms only; proof kernel is independent",
            },
            activation: "the event-jump target card is frozen, and a pre-candidate manifest successor has frozen the complete event machine AST, symbol definitions, runtime-premise schema, and total AST-to-Lean translation required by FORMAL_PROJECTION; it is separate from the class-preserving naturality theorem",
            kill: "checker or declaration-binding rejection, a missing/duplicate owner, failed internal-transfer cancellation, undefined pairing, unproved weak limit/renormalization, or one premise-satisfying counterexample refutes the bound event theorem",
            fallback: "stop at the event boundary, retain the complete pre-limit family and ownership ledger, and return claim adjudication Unknown with observable disposition DeclaredDivergent for every divergent/unproved post-event observable (unless the requested predicate itself is divergence)",
            no_claim: "cohomology/class change alone predicts neither a finite force jump nor impulse; divergent limits, scheme-dependent finite parts, and unpinned contact/breakdown regularizations remain explicit outcomes rather than being coerced finite",
        },
        ClaimSpec {
            id: "i03-topology-force-counterexample-search",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Refutation,
            statement: "Falsifier lane: exhaust a cardinality-audited canonical microgrammar, \
                        formally quantify over the full bounded parameter grammar, and \
                        search a separately seeded larger grammar for adversarial \
                        variational sheaves, weighted actions, transports, refinements, \
                        homotopies, harmonic-coordinate changes, and regularized topology \
                        events inside each bound theorem domain. Promotion requires every \
                        preregistered nonvacuity/coverage floor plus zero independently \
                        verified in-domain counterexamples; every valid candidate and \
                        minimized counterexample is retained",
            hypotheses: &[
                "the topology-force adversary fixture pins a cardinality-audited exhaustive microgrammar, symbolic full-domain theorem parameters, and a separately seeded larger grammar, with finite dimensions, cell/conductor/event counts, numerator/denominator/norm/sparsity/degree bounds, differentiable maps, refinement depth, evaluation budgets, and full decorated-object canonicalization",
                "an independently checked rank/unrank enumeration-completeness certificate binds the microgrammar cardinality and complete decorated canonical quotient; the formal declarations quantify over the full bounded parameter grammar, and candidate validity plus exact theorem-domain membership are checked independently before counting",
                "the receipt reports nonzero nontrivial premise-satisfying witness counts for local-to-global naturality with d_x P_x nonzero, genuine refinement defects, and each topology-event regularization, so zero counterexamples cannot green a vacuous domain",
                "each candidate carries the exact manifest/card/declaration digest tuple used by the premise checker and independent minimizer",
            ],
            qoi: "exact_nonvacuity_coverage_and_zero_verified_in_domain_counterexample_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "i03.oracle.topology_force_falsifier.v1 at fs-vmanifest-oracles/i03/topology_force_falsifier.rs::verify_membership_cover_and_minimize",
                independent: true,
                tcb_overlap: "shares canonical candidate bytes and bound declaration digests; independent finite enumerator, premise checker, formal kernels, high-precision evaluator, coverage adjudicator, and minimizer are version-pinned in each receipt",
            },
            activation: "all declaration names, elaborated theorem types, manifest bindings, and theorem-card digests are frozen, and a pre-candidate FrozenManifest::amend successor has replaced M0 target prose with the complete executable grammar/encoding/predicate/tag/parameter/rank-unrank AST and receipts required by M0_FORMALIZATION before search or promotion",
            kill: "the first independently verified in-domain counterexample refutes exactly its bound declaration; an empty theorem domain, missed nonvacuity floor, grammar escape, rank/unrank or decorated-canonicalization defect, budget-incomplete microgrammar, or declaration-digest mismatch fails the campaign rather than passing it",
            fallback: "restrict or replace the theorem through an authenticated amendment while retaining every candidate, coverage failure, and counterexample",
            no_claim: "a bounded counterexample search is not a proof even when complete, nonvacuous, and empty; proof-kernel acceptance of the exactly bound declaration remains separately mandatory",
        },
    ]
}

#[allow(clippy::too_many_lines)]
fn i03_fixtures() -> Vec<FixturePin> {
    vec![
        FixturePin {
            id: "i03-parallel-plate-mms",
            source: FixtureSource::AuthoredSpec {
                spec: "I03_FIXTURE_V1 parallel-plate-mms. EQUATIONS: -div(epsilon grad \
phi)=rho_vol on [0,Lx]x[0,Ly]x[0,Lz], Dirichlet electrodes x=0,Lx; phi is the \
rational-coefficient trigonometric polynomial selected by Philox and rho_vol is symbolic \
substitution inside each material. For every discontinuous-epsilon interface Gamma with \
normal n from minus to plus, the manufactured free surface charge is explicitly \
sigma_Gamma=n dot (D_plus-D_minus) in C/m2 with D=-epsilon grad(phi); smooth-epsilon \
twins set sigma_Gamma=0. GRID: L_i in {1/2,1,2} m, epsilon_r layers in {1,2,7}, polynomial \
degree p in {1,2,3}, h levels {1/4,1/8,1/16,1/32}; body-fitted and planar cut \
offsets {1/7,2/7}. ACCEPTANCE: exact d1*d0; per cell s_G=|integral_boundary(D dot n)-integral_cell(rho_vol)-integral_interface(sigma_Gamma)|/\
(1e-12 C+5e-7 Q_scale), maximum s_G<=1; last-three-level energy order >=p-0.20, \
reconstructed-flux order >=p-0.35, and stored-energy goal-functional order >=2p-0.50, \
each with regression R2>=0.98. SCALE_BINDING sets \
Q_scale=max(1e-12 C,max_cell |Q_boundary,ref|,|Q_volume,ref|,|Q_interface,ref|) from the \
exact symbolic MMS on the identical oriented cell/interface supports; formula/unit/source \
digest and exact quadrature are fixed before candidate execution. SERIALIZATION: \
canonical little-endian IEEE-754 inputs and sorted cell ids. SEEDS: Philox alias \
'i03/parallel-plate', k=0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i03-coax-sphere-harmonic",
            source: FixtureSource::AuthoredSpec {
                spec: "I03_FIXTURE_V1 coax-sphere-harmonic. GEOMETRIES: coax radii \
(a,b) in {(1,2),(1,4),(2,5)} cm and spherical radii in the same ratios; analytic \
phi,E,D,C,W are evaluated at 256-bit precision. TOPOLOGY: scalar H0 gauges are \
constants per unanchored component; mixed-field annulus/torus cases carry an explicit \
integer H1 basis and a physically zero closed-loop electrostatic period vector; relative \
terminal voltages are separate trace data. NEGATIVE TWINS: omit one gauge constraint or \
harmonic basis vector, inject a nonzero loop period without EMF, or corrupt one relative- \
boundary incidence. MESH: p={1,2,3}, four uniform levels. ACCEPTANCE: component/rank \
identities exact and max(s_phi,s_E,s_D,s_C,s_W)<=1, where s_phi=|Delta phi|/(1e-11 V+1e-8 phi_scale), s_E=||Delta E||/(1e-9 V/m+1e-8 E_scale), \
s_D=||Delta D||/(1e-20 C/m2+1e-8 D_scale), s_C=|Delta C|/(1e-18 F+1e-8 C_scale), \
and s_W=|Delta W|/(1e-18 J+1e-8 W_scale). SCALE_BINDING sets each named scale to the \
maximum of its written absolute floor and the matching 256-bit analytic reference norm: \
L_infinity for phi, pinned L2 measures for E and D, absolute terminal C, and absolute stored \
W; norm/support/quadrature, formula id, unit id, exact value bits, and analytic-source digest \
are fixed before candidate execution. SEEDS: Philox alias 'i03/harmonic', k=0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i03-floating-conductors",
            source: FixtureSource::AuthoredSpec {
                spec: "I03_FIXTURE_V1 floating-conductors development. CASES: 2..8 \
disjoint ellipsoid conductors inside an explicit enclosing conductor, plus separate \
fixed-ground and infinity-boundary variants; axes {0.5,1,2} cm, gaps {0.2,0.5,1} cm, \
epsilon tensors with eigenvalues {1,2,5}epsilon0. SOURCE-FREE reciprocal cases test the \
all-conductor indefinite matrix. OFFSET cases freeze rho_f and permanent P, compute \
q0=q(V=0), and form every column C_:j=[q(Vstar*e_j)-q0]/Vstar without changing material \
state; plus/minus-Vstar agreement is checked. Nonreciprocal twins suppress symmetry \
authority. TERMINALS: fixed, floating, and charge constrained with separately typed \
closed, grounded, infinity, and quotient projections. ORACLE ROUTES: analytic coax/sphere \
fields; a 256-bit interval piecewise-linear Galerkin BEM only for its admitted source-free \
boundary class on nested 512m,2048m,8192m-panel meshes for m conductors, with enclosed \
geometry, order-8 singular quadrature, residual norm, coercivity lower bound, and observed \
energy-order interval lower bound 0.75; and an independent equilibrated volume/interface \
route for volume-source or polarization cases. Extrapolation is diagnostic, never a \
substitute for the residual enclosure; a case with no certified route is Unsupported. \
ORIENTATION: n_Omega points out of the dielectric and \
q_i=-integral_Gamma_i D dot n_Omega. With Q_vol=integral_Omega rho_f, the closed \
all-conductor vector includes the enclosure and uses 1^T q+Q_vol=0. An open/truncated \
vector excludes Gamma_ext, defines Q_exterior=integral_Gamma_ext D dot n_Omega, and \
separately uses 1^T q_terminal+Q_vol-Q_exterior=0. No case may count Gamma_ext in both q \
and Q_exterior. ACCEPTANCE routes predicates by contract. Every case requires its typed \
s_global=|closed or open ledger residual|/(1e-12 C+1e-8 Qstar)<=1 and \
s_C=||C-C_oracle||_F/(1e-18 F+1e-8 Cstar)<=1. SCALE_BINDING fixes \
Qstar=max(1e-12 C,||q_oracle||_2,|Q_vol|,|Q_exterior|) and \
Cstar=max(1e-18 F,||C_oracle||_F), with exact norms and \
independent-oracle source digests fixed before candidate execution. Closed reciprocal cases \
additionally require an exact self-adjoint bilinear-form/terminal-adjointness receipt proving \
C=C^T, an exact oriented column-ledger bit proving 1^T C=0, and an exact common-shift \
incidence bit proving C*1=0. The interval encoding shares one variable for each reciprocal \
pair C_ij=C_ji; s_sym=||C-C^T||_F/(1e-18 F+1e-8 Cstar)<=1 and \
s_rows=||C*1||_2/(1e-18 F+1e-8 Cstar)<=1 remain independent numerical corroboration, \
never theorem substitutes. For the exact full-column difference basis B spanning 1-perp, \
interval LDL and the generalized eigenproblem (B^T C B,B^T B) prove \
lambda_min>=1e-10 Cstar. Together with the exact symmetry and null receipts these certify \
interval PSD and ker(C)=span{1}. Every off-diagonal interval has upper endpoint \
<=0. Reduction contracts are distinct: grounding compares C_grounded to the declared \
principal minor C_uu; quotient coordinates compare C_quot to B^T C B; eliminating declared \
floating terminals compares C_eff to C_aa-C_af C_ff^{-1} C_fa using an interval-certified \
invertible C_ff; and an infinity/open-boundary operator compares only to its independently \
constructed open-domain oracle, never to an unnamed Schur complement. Each applicable \
reduction requires s_reduction=||C_reported-C_declared_reduction||_F/ \
(1e-18 F+1e-8 Cstar)<=1 and s_ref=||R_a q_a-R_b q_b||_2/ \
(1e-12 C+1e-8 Qstar)<=1 under pinned terminal reconstruction maps; interval directional \
derivatives certify declared terminal-monotonicity signs. Reduced cases have no zero-row- \
sum/kernel predicate. Nonreciprocal cases require \
conservation and oracle agreement but no symmetry predicate. Malformed negative twins must \
refuse; post-hoc structural projection is forbidden. \
SEEDS: Philox alias 'i03/floating', k=0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i03-floating-conductors-core-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "I03_FIXTURE_V1 floating-conductors CORE HOLDOUT. The exact \
generator/acceptance protocol is i03-floating-conductors, with unseen ellipsoid poses, \
boundary variants, charge offsets, and material tensors. SEEDS: Philox alias \
'i03/floating', k=65536..=69631. SEALED until the sole consumer \
i03-conductors-capacitance enters core adjudication; any earlier read is IntegrityFailed.",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "i03-resistor-conduction-mms",
            source: FixtureSource::AuthoredSpec {
                spec: "I03_FIXTURE_V1 resistor-conduction-mms. EQUATION: \
-div(sigma grad phi)=s on unit boxes/layered prisms; sigma symmetric-part eigenvalues \
in {0,0.1,1,10} S/m with skew twins and coercive-quotient flag. BOUNDARIES: voltage, \
current, mixed, and compatible/incompatible pure Neumann. FLUX: Raviart-Thomas or \
dual-cochain reconstruction at p={0,1,2}; exact polynomial MMS and series/parallel \
linear network references. ACCEPTANCE: max(s_cell,s_global,s_I,s_P)<=1, where each current \
balance score is |Delta I|/(1e-12 A+1e-8 I_scale), terminal-current score s_I=||Delta I_terminal||/(1e-12 A+1e-8 I_scale), and passive-loss score s_P=max(0,-P_Joule)/(1e-15 W+1e-8 P_scale); missing local-flux receipt is failure. \
SCALE_BINDING sets I_scale=max(1e-12 A,||I_network,ref||_2,max_face|I_MMS,ref|) and \
P_scale=max(1e-15 W,P_Joule,ref) from the exact MMS/series-parallel network on identical \
oriented terminals/faces; norm/support, formula/unit, quadrature, exact value bits, and \
reference-source digest are fixed before candidate execution. \
SEEDS: Philox alias 'i03/conduction', k=0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i03-dispersive-dielectric-cards",
            source: FixtureSource::AuthoredSpec {
                spec: "I03_FIXTURE_V1 dielectric-cards. LTI RATIONAL: Debye relaxation \
times tau={1e-9,1e-6,1e-3}s have poles p=-1/tau; damped Lorentz, proper direct- \
feedthrough, lossless-boundary-pole, active-residue, and wrong-sign twins are included. \
With exp(+i omega t), the physical boundary is s=0+i omega. Exact rational state-space/ \
KYP or equivalent algebraic certificates establish analyticity for Re(s)>0, realness, \
and Hermitian(Yp(s))>=0 for matrix Yp=s*Chi in S/m, including feedthrough and boundary- \
pole residues; stable poles alone never pass. NONLINEAR/HISTORY: each valid card pins an \
energy U(D,z), E=d_D U, state evolution 0 in partial_(dot z)Phi+partial_z U, state and \
boundary/history terms, and proves U(t1)-U(t0)<=integral E dot dot(D) dt-integral \
Diss dt with Diss>=0. Equivalent E-coenergy cards carry the exact Legendre certificate. \
Invalid twins preserve convex-looking syntax while breaking the constitutive power \
inequality. Incremental/IQC twins pin compared trajectories, filter state, initialization, \
and horizon and receive no absolute-passivity authority. EMPIRICAL cards pin a \
temperature/frequency/state hull and carry no energy theorem. Reciprocity is checked \
independently of passivity. SEEDS: Philox alias 'i03/dielectric', k=0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i03-eqs-regime-boundaries",
            source: FixtureSource::AuthoredSpec {
                spec: "I03_FIXTURE_V1 eqs-regime development. CASES: coax, parallel \
plate, RC dielectric, and conductive slab with L={1e-4,1e-2,1}m, log-spaced \
f=10^k Hz for integer k=1..11, omega=2*pi*f, epsilon_r={1,4,100}, mu_r={1,10}, \
sigma={0,1e-4,1e4}S/m. NORMS: the card freezes the relative-boundary Maxwell graph norm \
||.||_X, electric-energy norm ||.||_epsilon, projection onto the nonconservative/solenoidal \
electric component, positive gross injected/incident energy W_in,+, and requested-QoI \
dual norm. SCALE_BINDING provides exact positive fixture-only S_X,S_E,S_WE,S_Q and authored \
absolute floors X_abs,E_abs,W_abs,Q_abs with their unit, formula, and source digests. ROUTE M is \
allowed only for separable cases: a branch-certified modal Maxwell pencil defines \
eta_route,M=max over excitation-certified retained modes of sup_omega \
L_m|gamma_m(omega)|+U_tail/(X_abs+S_X), where U_tail is the discarded-mode bound in \
||.||_X. ROUTE R lifts the EQS field into the full-Maxwell complex, computes its exact \
residual R_M(E_EQS) in X*, encloses an upper bound R_hi>=||R_M(E_EQS)||_(X*) and a \
strictly positive lower bound beta_lo<=beta_R for the applicable full-Maxwell inf-sup \
constant, defines U_X,R:=R_hi/beta_lo, and proves \
||E_Maxwell-lift(E_EQS)||_X<=U_X,R before defining \
eta_route,R=U_X,R/(X_abs+S_X). \
For the selected route, the same certified enclosure supplies U_ind bounding the nonconservative electric component in \
||.||_epsilon, U_B bound omitted magnetic energy, U_rad bound outward radiated energy, \
and U_Q bound dual-weighted requested-QoI error. Define eta_ind=U_ind/(E_abs+S_E), \
eta_B=U_B/(W_abs+S_WE), eta_rad=U_rad/(W_abs+W_in,+), and \
eta_Q=U_Q/(Q_abs+S_Q). Signed net source work is never a denominator; zero/nonfinite \
gross input, an absent norm/projection, or an absent scale routes Unknown. Admission \
requires eta_route,M<=1e-2 or eta_route,R<=1e-2 as selected, eta_ind<=1e-2, \
eta_B<=1e-3, eta_rad<=1e-4, and \
eta_Q<=min(1e-2,0.25*b_Q), where b_Q>0 is the normalized requested-QoI budget. Missing \
modal completeness, tail, stability, or dual bound routes Unknown/escalate. Directional \
L/skin_depth and omega*tau_charge select only the pinned static [0,0.1], transition \
[0.1,10], or displacement [10,infinity] subcard; sigma=0 is a capacitive classifier, \
not rejection. For every upper-bound monitor u and threshold beta, boundary multipliers \
{0.5,0.9,1.0,1.1,2.0} are exercised: an enclosure with upper endpoint <=beta may Admit, \
one with lower endpoint >beta Escalates, and one straddling beta is Unknown/escalated. \
Stationary current integrates J_cond+dD/dt. Analytic/full-Maxwell solutions independently \
adjudicate field, current, and QoI error but never supply a missing admission certificate. \
SEEDS: Philox alias 'i03/eqs-regime', k=0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i03-eqs-regime-core-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "I03_FIXTURE_V1 eqs-regime CORE HOLDOUT. Uses the exact monitors, \
threshold multipliers, total-current definition, and reference routes of \
i03-eqs-regime-boundaries on unseen geometry/material/waveform tuples. SEEDS: Philox alias \
'i03/eqs-regime', k=65536..=69631. SEALED until sole consumer \
i03-field-circuit-force-adjoint enters core adjudication; any earlier read is \
IntegrityFailed.",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "i03-field-circuit-transients",
            source: FixtureSource::AuthoredSpec {
                spec: "I03_FIXTURE_V1 field-circuit-transients. NETWORKS: 1..8 field \
ports coupled to R/C, ideal voltage/current sources, floating nodes, and declared circuit- \
owned inductors. GEOMETRY: stationary capacitors plus prescribed translating comb and \
deforming-gap ports with analytic reference maps and zero-velocity twins. DAE index is \
{1,2} with consistent/inconsistent initial states; backward-Euler and midpoint use dt/T \
in {1/32,1/64,1/128}. PORT FLOW: with dielectric-outward n and I_k positive into the \
field, stationary ports prove I_k=-integral_Gamma_k(J_f+partial_t D) dot n \
=dot(Q_e_free,k)-dot(Q_e_transfer,k), where Q_e_free,k=-integral_Gamma_k D dot n and \
Q_e_transfer,k=integral_0^t integral_Gamma_k J_f dot n are disjoint owners. Only a blocking \
port with J_f dot n=0 identifies I_k with dot(Q_e_free,k). Moving ports pull back relative \
carrier flux and D to the reference port and include every swept/geometric term. The moving \
formula must satisfy terminal continuity and the discrete GCL \
exactly and reduce to the stationary formula when velocity is zero. SIGNS: I_k enters the \
field, V_k is terminal-minus-reference, Wfc=integral sum_k V_k I_k dt, source work is \
positive into its owner, and traction-times-boundary-velocity work is positive outward. \
FIELD LEDGER: Wfc+Wvolume=Delta(Wfield+Wmaterial)+Dfield+Wmech_out+Wnum_field. \
CIRCUIT LEDGER: Wexternal-Wfc=Delta(Wcircuit)+Dcircuit+Wnum_circuit. COMBINED LEDGER is \
their algebraic sum with no Wfc. Every Wnum is computed independently from the declared \
update identity before ledger comparison, never residual-fitted. ACCEPTANCE: \
max(s_I,s_charge,s_GCL,s_field,s_circuit,s_combined)<=1; current scores use \
1e-12 A+1e-8 Istar, charge/GCL scores 1e-12 C+1e-8 Qstar, and work scores \
1e-15 J+1e-8 Wstar. SCALE_BINDING fixes Istar,Qstar,Wstar as maxima of their stated \
positive floors and the matching independent reference-trajectory norms or gross owned-work \
magnitudes, with exact orientation, quadrature, formula, unit, and source digest before \
candidate execution. REFERENCE: analytic moving- \
capacitor identities plus independent 256-bit descriptor integration. SEEDS: Philox \
alias 'i03/circuit', k=0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i03-force-held-variable-benchmarks",
            source: FixtureSource::AuthoredSpec {
                spec: "I03_FIXTURE_V1 force-held-variable development. CONFIGURATIONS: \
parallel plate, comb drive, dielectric interface, and 2..4 floating conductors; admitted \
translation, rotation, and shape-coordinate charts with their tangent spaces, dual \
pairings, coordinate units, constraint reactions, and pinned configuration metrics; \
paths are smooth with fixed topology/cut class. ENSEMBLES: a scalar-coenergy card pins \
Wstar with q=d_V Wstar on its stated same-branch voltage domain and certifies the closed, \
exact terminal-charge one-form and path independence. The radial formula \
Wstar(x,V,z,h)-Wstar(x,0,z,h)=integral_0^1 V^T q(x,lambda V,z,h) d lambda is permitted \
only when 0 and the entire segment {lambda V | 0<=lambda<=1} lie in that domain; a merely \
simply connected domain uses a separately content-pinned in-domain path. For a closed \
all-conductor system, Vbar=R^m/span{1}, Qbar={Q | 1^T Q=-Q_vol}, q0 is in Qbar, and \
qhat=Q-q0 is in 1-perp. Define Wbar*([V])=Wstar(V)-<V,q0> with its exact source/offset \
term, and W_Q(x,Q,z,h)=sup_[V]in Vbar {<[V],qhat>-Wbar*([V])}. The declared gauge-section implementation \
proves section independence, affine offset/source covariance, reduced strict convexity, \
and quotient reconstruction. Grounded/open cards use their own coercive spaces; an \
incompatible Q is Unsupported. Nonintegrable, nonreciprocal, history-switching, and \
mixed/floating cases use the complete constrained action A without a scalar-coenergy claim. \
Generalized force is +d_x Wstar at fixed V, -d_x W at fixed Q, or d_x A for the declared \
mixed action. Only reciprocal linear \
zero-offset q=C V reduces to the two quadratic capacitance formulas. STENCILS: complex- \
step only when the entire composed geometry/solve/state/QoI path is holomorphic; otherwise \
8th-order symmetric/Ridders differences carry interval truncation+roundoff enclosures at \
steps 2^-8..2^-20. INTERFACES: total-energy traction, normal, admissible virtual- \
displacement trace, integrated resultant, and moment arm are pinned. ACCEPTANCE: for the \
frozen unit-normalized admissible test-variation set Vunit, max(s_dual,s_VW,s_adj)<=1. \
s_dual=sup_v |(G_energy-G_traction)[v]|/(1e-15 J+1e-7 Wstar_scale); \
s_VW=|Delta A-integral_path G(q)[dq]|/(1e-15 J+1e-7 Wstar_scale) with certified path- \
quadrature remainder; s_adj uses the same dual-work pairing. SCALE_BINDING fixes \
Wstar_scale=max(1e-15 J,sup_{v in Vunit}|G_oracle[v]|,|Delta A_oracle|) in the \
unit-normalized dual-work pairing from fixture-only independent reference actions before \
candidate execution. Coordinate reports use \
1e-12 N for translations, 1e-15 J/rad for rotations, and a declared unit-bearing \
absolute floor for every other derivative coordinate; every scale is positive and frozen \
before execution. SEEDS: Philox alias 'i03/force', k=0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i03-force-held-variable-core-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "I03_FIXTURE_V1 force-held-variable CORE HOLDOUT. Uses the exact \
ensemble formulas, regularity predicates, stencils/enclosures, normalization, and \
interface conventions of i03-force-held-variable-benchmarks on unseen off-grid \
paths. SEEDS: Philox alias 'i03/force', k=65536..=69631. SEALED until sole consumer \
i03-field-circuit-force-adjoint enters core adjudication; any earlier read is \
IntegrityFailed.",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "i03-insulation-adversaries",
            source: FixtureSource::AuthoredSpec {
                spec: "I03_FIXTURE_V1 insulation-adversaries development. FAMILIES: \
void, needle, triple junction, surface discharge, and contaminated interface. FACTORS: \
gas {air,N2,SF6-surrogate} with pressure {0.5,1,2}bar, waveform {DC,50Hz,1kHz,pulse}, \
temperature {253,293,333}K, relative humidity {0.1,0.5,0.9}, and right censoring. \
DEVELOPMENT CLUSTERS: 512 deterministic material-lot-shaped clusters; each consumes exactly eight consecutive case \
indices and evaluates all five fixture families and both endpoints on every repeat. For a \
partition beginning at seed k0, cluster=floor((k-k0)/8) and repeat=(k-k0) mod 8. For \
each family/endpoint cell, the cluster outcome is the worst-covered repeat, so within-lot \
dependence is arbitrary. This public generator validates algorithms and replay but mints no \
IID or confidence authority; maximal promotion uses the separately governed sampled-lot \
receipt. POPULATION RNG: the output-block ordinal reserves its high 16 bits for factor_id \
and low 48 bits for that factor's local block, so rejection in one factor cannot shift any \
other factor. Lot-level factors use case index k0+8*lot and ids {1:gas,2:pressure, \
3:waveform,4:temperature,5:humidity,6:defect-scale,7:material-threshold}; repeat-level \
measurement factors use the actual k and ids {8:calibration-gain,9:apparent-charge-offset, \
10:voltage-noise,11:stress-noise,12:right-censoring}. For ids 1..5 the local draw starts \
at j=0; for family-specific ids 6..7 it is j=family_index; for ids 8..11 it is \
j=2*family_index+endpoint_index; and id 12 owns the pair \
j=2*(2*family_index+endpoint_index)+{0,1}. Gas, pressure, waveform, temperature, and \
humidity use rejection-sampled uniform choices in their own lanes; defect scale and material \
threshold are independent log-affine draws over [0.5,2] times their family nominal values; \
calibration gain is uniform [0.98,1.02], apparent-charge offset is uniform [-0.5,0.5] pC, \
and voltage and stress noise are two independent uniforms on [-0.005,0.005] of their \
separately declared domain spans. The first id-12 draw censors with probability 0.20; when \
censored, its second draw sets the apparatus limit L_c=(0.5+0.5*u)*nominal-family-span \
inside the pinned finite label domain, and only a true label above L_c is emitted as the \
right-censored compatible set [L_c,2*nominal-family-span]. Otherwise the exact synthetic \
label is emitted. Censoring is therefore repeat/family/endpoint owned, never inferred after \
coverage. The exact analytic family law and a separate 256-bit \
interval label evaluator are serialized per case; no learned candidate generates labels. \
EVENTS: partial discharge is apparent charge >=10 pC within the \
100kHz..20MHz band. AC breakdown is a conductive transition sustained for >=5 drive \
periods; pulse breakdown uses >=5 explicitly declared repetition periods; DC breakdown \
requires current density >=1 A/m2 and field collapse >=80% continuously for the pinned \
absolute dwell 10 ms. WIDTHS: for family f and endpoint e, freeze S_fe>0 and label domain \
D_fe=[0,2*S_fe]. Every predicted P_lrfe is a nonempty closed interval in D_fe; the compatible \
observation C_lrfe is {y} for an exact label or [L_c,2*S_fe] for a censored label. Coverage \
is exactly C_lrfe subseteq P_lrfe. Define H_lrfe=conv(P_lrfe union C_lrfe) and \
w_lrfe=length(H_lrfe)/(2*S_fe); empty, reversed, nonfinite, or out-of-domain intervals fail. \
The lot-cluster coverage bit is the conjunction over eight repeats and its width is \
max_r w_lrfe. STATISTICS: ten cells=five families x {PD inception,breakdown}; \
nominal directed coverage 0.90; each cell uses a one-sided exact Clopper-Pearson lower \
bound at alpha=0.001, so Bonferroni gives simultaneous familywise confidence >=99%, \
and every lower bound must be >=0.90. Exact confidence authority is conditional on the \
independently audited IID sampling receipt; deterministic development outputs receive only \
descriptive coverage. Sharpness is explicitly empirical on the \
512-cluster holdout: after sorting normalized cluster widths nondecreasing, the upper \
median w_(257) must be <=0.20 and nearest-rank p90 w_(ceil(0.90*512))=w_(461) <=0.35 \
in every cell; no population-quantile claim is inferred. Out-of-hull routes Unknown/ \
Unsupported. SEEDS: Philox alias \
'i03/insulation', k=0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i03-insulation-adversaries-max-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "I03_FIXTURE_V1 insulation-adversaries MAX HOLDOUT. Uses the exact \
families, eight-repeat within-lot construction, conservative censoring, width formula, event \
labels, nominal coverage, ten cellwise alpha=0.001 exact lower bounds with simultaneous \
familywise confidence >=99%, and per-cell width caps of i03-insulation-adversaries. \
STATISTICAL AUTHORITY: the roster is exactly \
{i03-stat-custodian-a,i03-stat-custodian-b,i03-stat-custodian-c}, requires 3-of-3 reveals, \
and binds each verification-key digest plus the commitment scheme/encoding in the admission \
receipt. Candidate/model/toolchain, manifest, exact finite population atom table, censoring \
rule, normalization, and stopping rule are irrevocably content-committed first. Only afterward \
does each custodian generate its entropy vector and submit a commitment into a sealed \
simultaneous phase in which every other custodian's commitment bytes are withheld. Each \
custodian irrevocably signs before any of the three commitment records are released; only the \
all-three-fixed phase receipt releases them in custodian-id order, after which the beacon round \
may occur. The authenticated transcript schema is I03StatHoldoutTranscriptV1, an exact tagged \
union. Every record encodes u32-big-endian-length-framed domain, u8 stage tag \
{CANDIDATE_FIXED=0,COMMIT=1,COMMIT_SET_FIXED=2,BEACON=3,REVEAL=4,FINAL=5}, u64-big-endian \
phase_seq, 32-byte link digest, fixed 32-byte manifest, candidate, model, toolchain, \
population, and campaign digests, u32-framed actor id, exact 32-byte actor-key digest, then only \
its stage payload and 64-byte signature; there are no optional/common commitment fields or \
zero sentinels. The exact domain value is org.frankensim.i03.stat-holdout.transcript.v1. \
Signing bytes are every preceding record byte through the payload, excluding \
only the signature. A record digest is BLAKE3 derive-key domain \
org.frankensim.i03.stat-holdout.record.v1 over the entire encoded record including signature. \
The only valid automaton has ten records and contiguous phase_seq 0..9: one \
CANDIDATE_FIXED by actor i03-stat-governance; exactly three COMMIT records in written roster \
order; one COMMIT_SET_FIXED by governance whose payload is the three custodian-id-ordered \
COMMIT record digests; one BEACON by governance; exactly three REVEAL records in the same \
roster order; and one FINAL by governance. CANDIDATE_FIXED payload is the candidate-semantics \
root, defined by BLAKE3 derive-key domain \
org.frankensim.i03.stat-holdout.candidate-root.v1 over the six common manifest/candidate/ \
model/toolchain/population/campaign digests in that written order; \
COMMIT payload is the 32-byte entropy commitment; COMMIT_SET_FIXED payload is exactly \
the three 32-byte COMMIT digests; BEACON payload is u32-framed exact source UTF-8, round u64, \
32-byte B, and 32-byte source-authentication-receipt digest; REVEAL payload is its own \
32-byte COMMIT record digest, exact 32768-byte vector, and 32-byte salt. Link semantics form \
an exact join DAG compatible with peer-withheld commitments: record 0 uses 32 zero bytes; \
every COMMIT links to CANDIDATE_FIXED; COMMIT_SET_FIXED also links to CANDIDATE_FIXED and \
joins the three COMMIT digests; BEACON links to COMMIT_SET_FIXED; every REVEAL links to \
BEACON and names its own COMMIT digest; FINAL links to BEACON. FINAL payload contains the \
three roster-ordered REVEAL digests, pre-FINAL transcript root, and sampler-output root. The \
pre-FINAL root is BLAKE3 derive-key domain \
org.frankensim.i03.stat-holdout.transcript-root.v1 over the ordered u32-length-framed nine \
record digests for phase_seq 0..8, each already covering its signature. The sampler-output \
root is BLAKE3 derive-key domain org.frankensim.i03.stat-holdout.sampler-root.v1 over 512 \
lot-index-ordered u32-length-framed SamplerLotV1 records. Each record is exactly u16-big- \
endian lot index, 32-byte U_LOT, 32-byte U_CANDIDATE_SEED, u32-big-endian-length-framed UTF-8 \
selected atom id, 32-byte candidate-output digest, then ten CellResultV1 values in outer \
family order {void,needle,triple-junction,surface-discharge,contaminated-interface} and inner \
endpoint order {PD_INCEPTION,BREAKDOWN}. CellResultV1 is coverage:u8 in {0,1} followed by the \
correctly rounded roundTiesToEven IEEE-754 binary64 value of the exact normalized cluster- \
width, encoded as u64 big-endian bits; mathematical zero must be canonical +0.0 bits=0 and \
-0.0 is rejected. Noncanonical atom ids, nonfinite/out-of-[0,1] widths, wrong count/order, or \
trailing bytes are IntegrityFailed. \
FINAL is therefore \
nonrecursive, and openings cannot replay against another beacon. No extra, missing, reordered, \
wrong-link, wrong-actor, or wrong-payload record is decodable. Custodians sign COMMIT/REVEAL; the admission-bound governance \
key signs CANDIDATE_FIXED/COMMIT_SET_FIXED/BEACON/FINAL. Ed25519 uses RFC8032 PureEdDSA \
with SHA-512, empty context, and no prehash, \
canonical point/scalar encodings, S<L, nonidentity A/R, exact prime-subgroup checks [L]A=0 \
and [L]R=0, k=SHA512(R||A||M) mod L, and the uncofactored equation [S]B_base=R+[k]A; \
noncanonical, off-curve, torsion/mixed-order, or equation-failing inputs are IntegrityFailed. \
Public-key digests use BLAKE3 derive-key domain \
org.frankensim.i03.stat-holdout.key.v1 over the exact 32 key bytes. The entropy commitment is BLAKE3 derive-key \
domain org.frankensim.i03.stat-holdout.commit.v1 over canonical length-framed custodian id, \
verification-key digest, exact 32768-byte vector, and exact 32-byte hiding salt, with the \
signed commit/reveal receipts fixing field order and lengths and the reveal opening both vector \
and salt. The transcript is never candidate- \
visible and exact IID authority does not rely on computational hiding because the candidate is \
already fixed. Each commitment \
opens exactly 32768 bytes parsed as 1024 indexed 256-bit blocks: roles LOT and CANDIDATE_SEED \
crossed with lot indices 0..511, with no spare or variable-length entropy. For custodian c, \
r=LOT=0 or CANDIDATE_SEED=1, and l in 0..511, R_c(r,l) is exactly bytes \
[32*(2*l+r),32*(2*l+r+1)); this lot-major layout admits no padding or alternate role-major \
parse. Exact statistical \
authority is conditional on the explicit, independently audited assumption that at least one \
named custodian generated its entire 1024-block vector information-theoretically IID uniform \
and independent of the other custodians, every adversarial mask and commitment, the already- \
fixed candidate, and population construction; no independence from the order-only beacon is \
required. The sealed simultaneous phase makes \
later commitment-byte adaptation an IntegrityFailed transcript. The other committed vectors \
may be adversarial. Population probabilities are \
canonical integers n_i/2^256, n_i>=0 and sum_i n_i=2^256. Atom order is frozen; \
cumulative endpoints are exact unsigned 257-bit integers canonically encoded as 33-byte \
big-endian values: c_0=0, c_(i+1)=c_i+n_i, and c_N=2^256. A draw u in \
{0,...,2^256-1} selects the unique half-open interval [c_i,c_(i+1)); n_i=0 atoms are \
unreachable. Endpoints must be nondecreasing with no overflow, overlap, or gap. A future \
unpredictable exactly 256-bit public beacon B identified by round/source is bound only after \
COMMIT_SET_FIXED. For role r in {LOT=0,CANDIDATE_SEED=1} and lot l, \
U_(r,l)=R_a(r,l) xor R_b(r,l) xor R_c(r,l). Every reveal/U block is exactly 32-byte \
big-endian; roles use the written numeric tags and l is canonical u16 big-endian. This \
coordinate-wise bijection preserves joint uniformity \
supplied by the honest vector. Interpreting U_(LOT,l) as an unsigned big-endian integer and \
selecting the unique frozen cumulative interval records exactly 512 IID material-lot atoms; \
U_(CANDIDATE_SEED,l) is the disjoint IID per-lot candidate seed. No hash, XOF, PRG, rejection \
stream, short mixed seed, or repeated draw may mint lot-level independence; no candidate or \
single evaluator can choose or grind the realized holdout. B is not a sampling input: its low \
9 bits select an exact cyclic shift of lot evaluation order, its next bit selects direction, and \
the remaining bits are retained as an unused challenge nonce; order-invariance is mandatory, \
so no assumption about B's distribution can alter the lot law or exact confidence authority. \
Every committed reveal is mandatory: a missing/late/duplicate reveal, wrong key, beacon \
substitution, abort, retry, or resampling is IntegrityFailed with no replacement holdout. \
within-lot repeats may be arbitrarily dependent. A secret campaign nonce and raw lot bytes \
are commitment-bound and inaccessible to the candidate until one-shot adjudication is final, \
then revealed and retained for deterministic replay. EXECUTION applies one identical measurable \
candidate map in a fresh stateless scope for every lot. Its complete candidate-visible input is \
exactly {features derived from U_(LOT,l), U_(CANDIDATE_SEED,l), frozen campaign constants}; \
lot ordinal, stable case/receipt id, processing order, shard, worker, campaign nonce, beacon, \
other entropy blocks, and logging metadata remain outside the candidate sandbox. No prior \
label, prediction, feature, cache, model update, RNG state, stopping decision, or lot order may \
affect another lot. Candidate randomness is a \
deterministic function only of the disjoint U_(CANDIDATE_SEED,l) block and frozen artifact and \
is included in each lot receipt; any within-lot expansion is merely part of that deterministic \
randomized-algorithm map and mints no additional independent draw. A secret order permutation \
and a fresh random reassignment of all external case/receipt ids must each produce byte-identical \
per-lot predictions before identically distributed outcome-level IID authority exists. \
Philox may expand retained within-lot inputs into repeat-local simulation noise but mints no \
lot-level IID premise. Exact Clopper-Pearson authority is conditional on the IID receipt and \
the at-least-one-honest uniform-source assumption; absent or compromised entropy, precommitment, access \
control, or one-shot discipline is IntegrityFailed and yields at most descriptive empirical \
coverage. The protocol is public and preregistered; only the committed nonce, raw lots, \
features, predictions, and labels are SEALED until sole consumer i03-breakdown-routing \
enters maximal adjudication, and any earlier realization-byte read is IntegrityFailed. \
Synthetic labels are not experimental validation.",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "i03-space-charge-aging-limits",
            source: FixtureSource::AuthoredSpec {
                spec: "I03_FIXTURE_V1 space-charge-aging development. TYPED RESERVOIRS: \
nonnegative mobile, trapped, and immobile bulk concentrations c_bulk have declared \
number-density units, signed per-particle charges q_bulk in C, and integer closed-reaction \
matrix S_bulk satisfying q_bulk^T S_bulk=0 exactly for every dimensionally homogeneous \
reaction column. Interface populations carry surface-density units, electrodes carry lumped \
C, and imposed sources are oriented port transfers; they are never concatenated into \
c_bulk. Negative twins delete one bulk product, surface transfer, or electrode owner. Define \
rho_f=q_bulk^T c_bulk plus separately declared fixed bulk charge, and prove \
d_t rho_f+div J_f=0 with oriented surface/port transfers through the integrated control- \
volume incidence ledger. POLARIZATION: globally rho_b=-div P and J_b=partial_t P as \
distributions. For each fixed interface with n from minus to plus, \
sigma_b=n dot (P_minus-P_plus) in C/m2 and \
dot(sigma_b)+n dot (J_b_plus-J_b_minus)=0 in A/m2. Moving interfaces use the pinned \
spacetime pullback, relative flux, and geometric-conservation term. Shared incidence proves \
d_t rho_b+div J_b=0 and div(epsilon0 E)=rho_f+rho_b without evolving any bound \
quantity independently. ELECTRODES: with dielectric-outward n, carrier transfer is \
Q_e_transfer=integral_0^t integral_Gamma_e J_f dot n; electrode free charge is separately \
Q_e_free=-integral_Gamma_e D dot n. With external total current positive into the field, a \
fixed port proves I_total=-integral_Gamma_e(J_f+partial_t D) dot n \
=dot(Q_e_free)-dot(Q_e_transfer); only a blocking port with J_f dot n=0 reduces this to \
dot(Q_e_free). Moving ports use the corresponding pulled-back relative-carrier plus D, \
geometric, and GCL ledger. The integrated bulk-plus-surface-plus-electrode-plus-port ledger includes each \
source exactly once. TRAP CAPACITY: each material card declares either one shared trapped- \
species group or disjoint per-species groups, positive dimensionless rational site weights \
w_i, and rho_trap_occ=sum_(i in group) w_i*|q_i|*c_trapped,i in C/m3. Each trapped species \
belongs to exactly one group. Its charge-density capacity rho_trap_cap is in C/m3; no implicit \
conversion to number density or signed-charge cancellation is allowed. PARAMETERS: mobility \
{1e-14,1e-10,1e-6}m2/Vs, rho_trap_cap {1e-6,1e-3}C/m3, \
injection/recombination Damkohler numbers {1e-3,1,1e3}. STATES: positivity, capacity, \
thermal {253..373}K, damage [0,1], and polarization history. CARDS: gradient-flow \
storage/dissipation and empirical-aging no-energy-authority twins. LIMITS: sigma or \
mobility lambda in {1e-6..1e6}; each case pins nondimensional scaling, limiting PDE/DAE, \
boundary conditions, L2/graph norm, index change, uniform estimate, and refusal band. \
ACCEPTANCE: every support keeps its own unit. The maximum of s_Q,s_bulk_cont, \
s_surface_free_cont,s_gauss,s_bound_bulk,s_bound_surface,s_bound_surface_cont,s_moving, \
s_electrode,s_state_bulk,s_state_surface,s_state_electrode,s_capacity,s_thermal,s_damage,s_energy, \
and s_limit must be <=1. Global/electrode charge scores use 1e-12 C+1e-7 Q_scale, and \
electrode-current scores use 1e-12 A+1e-7 I_scale. Bulk free continuity uses \
||d_t rho_f+div J_f||/(1e-12 A/m3+1e-7 jdiv_scale); free surface continuity and the \
fixed-interface bound continuity each use 1e-12 A/m2+1e-7 jsurf_scale. Gauss and \
rho_b=-div P bulk equivalence use 1e-12 C/m3+1e-7 rho_scale; \
sigma_b=n dot(P_minus-P_plus) uses 1e-12 C/m2+1e-7 sigma_scale. Each moving-interface \
spacetime weak-balance pairing against a frozen dimensionless test uses \
1e-15 C+1e-7 Q_test_scale. State violations are support-specific: bulk number \
concentration uses 1 m^-3+1e-7 c_bulk_scale, surface number concentration uses \
1 m^-2+1e-7 c_surface_scale, electrode state uses 1e-12 C+1e-7 Q_scale, temperature \
uses 1e-6 K+1e-7 T_scale, trap-capacity excess uses \
max(0,rho_trap_occ-rho_trap_cap)/(1e-12 C/m3+1e-7 rho_trap_cap_scale), damage uses \
1e-12+1e-7 d_scale, and polarization uses \
1e-12 C/m2+1e-7 P_scale. Energy uses 1e-15 J+1e-7 W_scale; limiting- \
norm excess uses its preregistered unit-bearing positive enclosure width. Exact q_bulk^T \
S_bulk, support-dimension typing, and incidence-complex identities are additional bit \
predicates. SCALE_BINDING defines every named scale as the maximum of its declared positive \
absolute floor and the matching independent reference norm on the identical bulk, surface, \
electrode, port, or test-function support; exact measures/quadrature, formula ids, unit ids, \
IEEE bits, and source digests are fixed before candidate or holdout access. SEEDS: Philox \
alias 'i03/space-charge', k=0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i03-space-charge-aging-max-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "I03_FIXTURE_V1 space-charge-aging MAX HOLDOUT. Uses the exact \
species/stoichiometry, state constraints, card classes, singular scalings, limiting \
problems/topologies, and normalized acceptance of i03-space-charge-aging-limits on \
unseen parameter tuples. SEEDS: Philox alias 'i03/space-charge', k=131072..=135167. \
SEALED until sole consumer i03-space-charge-aging enters maximal adjudication; any \
earlier read is IntegrityFailed.",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "i03-electrostriction-finite-strain",
            source: FixtureSource::AuthoredSpec {
                spec: "I03_FIXTURE_V1 electrostriction development. DEFORMATION: \
orientation-preserving affine/cubic W^{1,p}, p>3, maps with singular values in [0.5,2] \
on bounded connected Lipschitz reference domains. ROUTE H (GLOBAL HOMEOMORPHISM) requires a \
C1 map on the closure, interval Bernstein det(F)>=0.2, a proper degree-one local- \
diffeomorphism certificate, and an interval boundary-degree/self-separation proof that the \
trace is an embedding. ROUTE CN (A.E. INJECTIVITY) binds the continuous/Lusin-N \
representative, positive determinant a.e., area formula, trace/degree data, and independently \
computes the image measure from certified boundary winding/degree cells before checking \
integral_Omega det(F)<=measure(y(Omega)); the image computation may not assume injectivity. \
The exact pinned Ciarlet-Necas theorem then establishes only a.e. injectivity. Every receipt \
states H or CN, and sampling alone has no authority. ENERGY: \
frame-indifferent rational polynomials Psi(C,E_ref,z,A_ref), C=F^T F, built from \
objective contractions with pinned reference structural tensors A_ref, plus the D_ref \
fixed-charge Legendre dual. Every scalar, referential/spatial vector, and second/fourth- \
order internal component of z has an explicit material-frame transformation; superposed \
spatial rotation sends F to R F and spatial vector/tensor state by R while leaving the \
referential invariant representation unchanged. Deliberately componentwise-F \
nonobjective and unstable cards are invalid twins. STABILITY is independently typed from \
Route H/CN: REDUCED certifies electrical/internal block elimination and a coercive \
mechanical Schur Hessian, while MIXED certifies the bordered/KKT-Hessian isomorphism, \
inf-sup beta>0, and fixed Morse index without a full-Hessian positivity claim; all four \
invertibility/stability cross-products are represented. \
Coefficients are {1/2,1,2}. INTERFACES: two-phase planar/curved patches with fixed \
referential/spatial normal, total Piola/Cauchy and configurational stress; each card \
either pins zero surface energy or includes its surface-energy first variation exactly. \
LIMITS: \
contrast lambda={1e-4..1e4}, uniform selected-route stability (REDUCED Schur coercivity or \
MIXED inf-sup/bordered isomorphism) and nonuniform negative twins. ORACLE: \
exact symbolic derivatives plus 256-bit interval virtual work. ACCEPTANCE: \
max(s_energy,s_stress,s_traction,s_resultant,s_VW,s_adj)<=1; energy/work use \
1e-15 J+1e-7 Wstar, stress and pointwise traction use 1e-9 Pa+1e-7 Sstar, integrated \
resultant uses 1e-12 N+1e-7 Fstar, integrated moment uses 1e-15 J/rad+1e-7 Mstar, \
and each adjoint component uses its declared derivative unit with an explicit positive \
absolute floor plus 1e-7 relative scale. SCALE_BINDING defines Wstar,Sstar,Fstar,Mstar and \
each component scale as the maximum of the stated positive floor and the matching \
independent symbolic/interval reference norm on the identical support, measure, and \
quadrature; candidate outputs cannot set them. SEEDS: Philox alias \
'i03/electrostriction', k=0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i03-electrostriction-max-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "I03_FIXTURE_V1 electrostriction MAX HOLDOUT. Uses the exact \
deformation, ensemble, energy, interface, regularity/stability, contrast-limit, and \
normalized-acceptance protocol of i03-electrostriction-finite-strain on unseen \
coefficients and curved patches. SEEDS: Philox alias 'i03/electrostriction', \
k=131072..=135167. SEALED until sole consumer i03-electrostriction-theorem enters \
maximal adjudication; any earlier read is IntegrityFailed.",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "i03-electrostriction-theorem-card",
            source: FixtureSource::AuthoredSpec {
                spec: "I03_THEOREM_CARD_V1 electrostriction. VARIABLES: C=F^T F, E_ref,z \
and reference structural tensors for the objective fixed-voltage potential, with \
F,D_ref,z for its fixed-charge Legendre dual and the exact referential/spatial electric \
transforms. PREMISES: Psi is C2 on the declared open state set and frame-indifferent. \
INVERTIBILITY ROUTE H binds the C1-closure, positive-Jacobian, properness, degree-one, and \
boundary-embedding hypotheses of the named global-homeomorphism theorem. ROUTE CN binds \
the bounded-Lipschitz domain, W^{1,p} p>3 continuous/Lusin-N representative, positive \
determinant a.e., area formula, trace/degree, independently certified image measure, and \
integral det(F)<=image-measure hypotheses of the named Ciarlet-Necas a.e.-injectivity \
theorem. The conclusion records global homeomorphism for H and only a.e. injectivity for CN. \
STABILITY ROUTE: exactly one checker-visible route is selected. REDUCED uniquely solves \
electrical/internal stationary variables on their gauge/constraint quotient with a certified \
field-block isomorphism, forms the exact reduced Schur Hessian, and proves it coercive on \
the mechanical constraint tangent. MIXED retains saddle variables and proves the bordered/KKT \
Hessian is an isomorphism with certified inf-sup beta>0 and fixed Morse index; it makes no \
positive-definite full-Hessian claim. The chosen route proves local uniqueness on the \
declared branch and \
implicit-function regularity. A singular-limit corollary requires a uniform lower bound on \
the same reduced-coercivity or mixed inf-sup constant. Every \
internal-state objectivity law, trace, interface normal, and zero/nonzero surface-energy \
owner is fixed; topology/branch fixed. CONCLUSIONS: first variation equals total Piola \
stress and configurational interface jump including any surface energy; spatial pushforward \
gives Cauchy traction; constrained virtual work and exact-discrete adjoint are identical. LIMIT \
COROLLARY additionally requires uniform selected-route stability (reduced-Schur coercivity \
or mixed inf-sup/bordered isomorphism), compactness, strong state and weak \
flux convergence, and uniform derivative domination. Any absent premise yields \
Unknown, never theorem success. FORMAL PROJECTION GATE: this manifest-version-1 authored \
card freezes the theorem target and scientific no-claim boundary but contains no complete \
machine proposition AST, symbol-definition closure, or total translation, so it mints no \
theorem promotion authority. Before any proof/candidate bytes exist, a FrozenManifest::amend \
successor must replace this card with canonical machine AST bytes, every referenced definition \
or content digest, a total runtime-premise schema mapping each formal hypothesis to evidence, \
and the exact version/source digest of a deterministic AST-to-Lean translator; generated and \
elaborated types must structurally round-trip to that AST. A later admission receipt may not \
choose or weaken the projection for this version. AXIOM POLICY: exact policy \
i03.lean-axioms.v1 permits only the fully qualified names {propext,Quot.sound,Classical.choice}; \
its digest is BLAKE3 derive-key domain org.frankensim.i03.lean-axioms.v1 over the lexically \
sorted sequence of u32-little-endian byte length plus exact UTF-8 name. Every receipt binds \
that digest and the complete transitive declaration/environment axiom closure; sorryAx, any \
custom/theorem-equivalent postulate, and native-oracle proof authority outside the admitted kernel/checker TCB are \
IntegrityFailed. BINDING: a canonical table binds \
i03-electrostriction-interface-theorem to its schema-v2 manifest canonical claim_digest \
using the retained claim.v1 component domain, covering every identity-bearing ClaimSpec \
field, and separately hashes the \
explicit theorem projection. The serialized card, generated Lean proposition, elaborated \
checked proposition, and exported receipt carry mutually verified semantic/byte digests; \
declaration, elaborated-type projection, environment, axiom report, Lean4 toolchain, and \
proof-term digests are mandatory, and any mismatch is IntegrityFailed. The sibling numerical \
closure claim is bound to its independent interval-oracle receipt and gains no theorem color \
from this export. EXPORT: \
proofs/i03/Electrostriction.lean::electrostrictionInterfaceClosure.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i03-topology-force-adversaries",
            source: FixtureSource::AuthoredSpec {
                spec: "I03_FIXTURE_V1 topology-force development. TYPED FINITE \
SUPERGRAMMAR: oriented abstract simplicial pairs (K,L) of dimension 1..3 with 1..=8 vertices, \
1..=12 simplices, L a subcomplex, <=4 conductor labels, and Betti numbers <=3; \
orientation/permutation is quotiented by canonical lexicographic labeling and every \
boundary-square and relative-subcomplex identity is checked exactly. CANONICAL OBJECT means \
the entire typed decorated tuple, not its bare complex: an isomorphism transports vertex/ \
simplex order and orientation signs, relative/conductor/support labels, cover nerve and \
restrictions, cochain bases, Hodge/material matrices by signed-permutation congruence, \
actions/state/QoIs by coordinate substitution, ports/traces/gauges/constraints, P/Q/maps/ \
homotopies by conjugation, and every theorem/declaration binding. The canonicalizer emits \
the lexicographically minimal full encoding plus an exact isomorphism witness; idempotence, \
orbit equality, automorphism stabilizers, and rank/unrank completeness are independently \
checked. Topology-only canonicalization has no completeness authority. COVERS: 1..=4 \
cover elements encoded as subcomplexes with their open stars, closed under every represented \
nonempty intersection; induced restrictions and exact Cech descent are checked through \
degree 3 and through every higher nerve/cochain degree consumed by the bound action. \
M0_TARGET_BEGIN. The intended exhaustive microgrammar has bounds dimension<=2, vertices<=5, \
simplices<=8, state dimension<=3, matrix dimension<=4 with <=8 nonzeros, <=3 cover \
elements, <=1 refinement, ordinary-polynomial degree<=2 with <=6 monomials, and ordinary \
reduced coefficients p/q with |p|<=2,1<=q<=2. Its written-order 16 strata are {stationary- \
naturality,filtered-cyclic-naturality,nonidentity-P,nonzero-H1,nonzero-d_x-P,overlap-descent, \
ensemble-defect,QoI-defect,homotopy-defect,signed-cancellation,dielectric-gap-event, \
conductive-neck-event,conductivity-bridge-event,birth,merge,DeclaredDivergent-control}; its \
written-order 16 checker tags are {identity,orientation,restriction,descent,chain,homotopy, \
action,gauge,ensemble,QoI,Hodge,material,load-trace-port,stability,event-ownership,binding}. \
The target factorization is exactly 16 strata x 16 distinct bases x 16 tags x 16 exact \
parameter tuples, hence N_micro=65536 and case rank \
(((stratum*16+base)*16+tag)*16+parameter). M0_TARGET_END. This target prose is not a \
machine grammar AST and mints no exhaustive authority. Before any candidate/search bytes \
exist, a FrozenManifest::amend successor must replace this fixture with: (1) a total typed \
BaseRecord schema covering the entire decorated object listed above; (2) canonical byte \
encodings and finite domains for every field; (3) formal, executable validity and each of the \
16 stratum-membership predicates; (4) exact executable semantics and expected verdict for \
every tag; (5) the arity, types, and explicit values of all 16 parameter tuples; (6) total \
candidate enumeration and lexical order with cross-stratum exclusion; and (7) exact rank/ \
unrank and shard functions. The successor binds the literal AST bytes, parser/checker source \
digests, independent schema-decoder receipt, derived N_micro=65536 proof, rank/unrank \
bijection receipt, worst-case cost/preflight, and complete shard/Merkle root before the \
candidate digest. A later admission receipt cannot choose any missing grammar semantics. \
Candidate bytes are forbidden inputs to parsing, inclusion, ordering, rank/unrank, or shard \
assignment; absent/mismatched machine bytes are IntegrityFailed and exhaustive authority is \
Unknown, never inferred from this target prose. The target digest is BLAKE3 derive-key domain \
org.frankensim.i03.microgrammar-target.v1 over the exact UTF-8 bytes from M0_TARGET_BEGIN \
through M0_TARGET_END and binds successor intent without pretending to be the future AST. \
RATIONALS: every scalar is unique reduced p/q with \
|p|<=16 and 1<=q<=16; matrices have dimension <=12, <=48 nonzeros, entry magnitude<=16, \
and SPD candidates pass exact principal-minor checks. SCALAR FUNCTION GRAMMAR: each action, \
load, trace, port, constraint, gauge, state, and QoI scalar is a reduced-rational polynomial \
in x and <=12 scalar state coordinates with total degree<=4, <=64 nonzero monomials, and \
coefficients from the bounded rational grammar; vector/tensor objects have dimension<=12 and \
<=128 total nonzero scalar-polynomial components. WEIGHTS/DATA: rational SPD Hodge/material \
matrices and only those finite scalar/vector/tensor objects. MAP FAMILIES: x belongs to the closed interval [-1,1]; entries of \
P_x,Q_x, restriction transformations, configuration map phi_x, and homotopies are \
rational polynomials of degree<=2 with the same coefficient grammar. Exact polynomial \
chain, inverse-homotopy, descent, action-naturality, and derivative identities certify each \
family; possibly rectangular chain/restriction maps use exact rank and required-minor profiles, \
homotopies use their exact typed identities, and only square configuration/chart maps that \
must be locally invertible use interval Jacobian-determinant or singular-value exclusion of \
zero. VARIATIONAL TRANSFER: condensation cases pin rational vertical/horizontal splittings, \
contractible fibers, interval-coercive vertical Hessians, unique stationary lifts, and exact \
reduced-action/Schur identities. Cyclic cases pin bounded-arity (<=3) rational L_infinity/BV \
higher maps on dimensions<=12, with <=128 nonzero bounded-rational tensor coefficients per \
arity. CYCLIC FILTER: complete descending pronilpotent filtrations, rational filtration \
degrees, nilpotence/completion depth, filtration-preserving maps/homotopies, nondegenerate \
cyclic pairings, master-equation identities, and terminating or convergence-certified higher \
series are bounded and checked. A homotopy inverse proves equivalence of the selected derived \
critical groupoids; the finite-polynomial alternative supplies explicit mutually inverse \
functors and checked natural equivalences. Invalid noncomplete/nonconvergent \
quasi-isomorphism twins are mandatory. REFINEMENTS: <=3 canonical \
stellar or barycentric subdivisions with induced transfers and bounded homotopies. True \
gauge/compensated exact changes leave the physical cochain identical; GL(b,Z) bases have \
entries |a_ij|<=3 and contragredient coordinates; uncompensated shifts are negative twins. \
EVENTS: <=2 elementary conductor births/merges and exactly three content-pinned primitive \
tokens EVENT(kind,k), encoded as bytes {kind:u8,k:u8} with kind GAP=0, ROUNDED_NECK=1, \
CONDUCTIVITY_BRIDGE=2 and k in 4..16. Each token owns the exact dyadic epsilon=(1,2^k); \
this pair is not an ordinary rational coefficient and is exempt only from the written q bounds. \
GAP is the signed parallel-boundary separation epsilon, ROUNDED_NECK uses the exact algebraic \
circle/sphere radius equation sum_i x_i^2=epsilon^2, and CONDUCTIVITY_BRIDGE uses the exact \
C1 cubic smoothstep h(t)=3*t^2-2*t^3 on t in [0,1]. These named primitive constructors and \
their fixed cubic coefficients are exempt only from ordinary polynomial degree/coefficient \
bounds, including M0's; all composed state/action/QoI and test polynomials remain inside their \
declared grammar. Charge/source/circuit/material ownership and \
distributional test polynomials of degree<=4 use the bounded rational grammar; each family \
also carries an interval epsilon-to-zero remainder enclosure. PROOF/FALSIFIER SPLIT: formal \
declarations quantify symbolically over the full bounded parameter supergrammar and never \
rely on enumerating its Cartesian product. A separate exhaustive microgrammar is admitted \
only with exact cardinality N_micro, rank/unrank bijection proof, worst-case per-case cost, \
shard partition, and a whole-campaign preflight proof covering exhaustive enumeration, formal \
checking, adversarial search, counterexample minimization, Merkle construction, checkpointing, \
and retention inside the frozen 12 h/64 GiB envelope. The committed campaign requires \
N_micro=65536, visits every index 0..N_micro-1 exactly once, and binds the independent shard/ \
Merkle completeness digest. The full 12-coordinate/degree-4 supergrammar receives exactly \
4096 deterministic search trajectories indexed by Philox case k=0..4095, each with exactly \
4096 candidate evaluations j=0..4095, for exactly 16,777,216 adversarial evaluations total; \
it is never called exhaustive. The preflight accounts separately for all 65,536 M0 cases, \
all 16,777,216 adversarial evaluations, and every named overhead. Budget exhaustion yields \
execution disposition BudgetExhausted with evidence completeness PartialEvidence, never a \
zero-counterexample pass. NONVACUITY: >=64 admitted naturality cases, >=64 refinement cases, >=32 \
birth and >=32 merge cases; coverage includes >=16 each for stationary-condensation and \
cyclic-L_infinity/BV transfer, nonidentity P, H1 nonzero, d_x P nonzero, nontrivial overlap \
descent, nonzero ensemble-transfer defect, nonzero QoI-transfer defect, nonzero homotopy \
defect, certified signed cancellation, nonzero event atom, and DeclaredDivergent control; \
at least one single-premise negative twin exists per premise class. ACCEPTANCE: all \
coverage floors hold, formal checker bits are exact, verified in-domain counterexample \
count is zero, and every numerical pulled-force discrepancy lies inside its dependency- \
preserving joint enclosure under the theorem-card dual-work normalization. SEEDS: Philox \
alias 'i03/topology-force', k=0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i03-topology-force-max-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "I03_FIXTURE_V1 topology-force MAX HOLDOUT. Uses the exact complex, \
weight/load/trace/port grammar, map/refinement, harmonic-coordinate, event/gap/contact, \
search-budget, full decorated canonicalization, minimization, and numerical acceptance \
protocol of i03-topology-force-adversaries for unseen cases beyond the exhaustive \
microgrammar. \
SEEDS: Philox alias 'i03/topology-force', k=131072..=135167. SEALED until sole consumer \
i03-topology-force-theorem-falsifier enters maximal adjudication; any earlier read is \
IntegrityFailed.",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "i03-topology-force-theorem-card",
            source: FixtureSource::AuthoredSpec {
                spec: "I03_THEOREM_CARD_V1 topology-force. FORMAL PROJECTION GATE: this \
manifest-version-1 authored card freezes the naturality, signed-refinement-defect, and event- \
jump theorem targets and no-claim boundaries but contains no complete machine proposition \
AST, definition closure, runtime-premise schema, or total translator, so it mints no theorem \
promotion authority. Before any proof/candidate bytes exist, a FrozenManifest::amend successor \
must replace it with canonical machine AST bytes for all three targets, every referenced \
definition or content digest, total schemas mapping each formal hypothesis to runtime evidence, \
and the exact version/source digest of a deterministic AST-to-Lean translator; generated and \
elaborated types must structurally round-trip to the corresponding AST. A later admission \
receipt may not choose or weaken a projection for this version. AXIOM POLICY: exact policy \
i03.lean-axioms.v1 permits only {propext,Quot.sound,Classical.choice}; its digest is BLAKE3 \
derive-key domain org.frankensim.i03.lean-axioms.v1 over the lexically sorted sequence of \
u32-little-endian byte length plus exact UTF-8 name. Every receipt binds that digest and the \
complete transitive declaration/environment axiom closure; sorryAx, any custom/theorem- \
equivalent postulate, and native-oracle proof authority outside the admitted kernel/checker TCB are IntegrityFailed. BINDING: \
a canonical table maps \
each theorem-authority manifest claim id to one fully qualified Lean declaration, binds the \
schema-v2 manifest canonical claim_digest using the retained claim.v1 component domain and \
covering every identity-bearing ClaimSpec field, \
and separately hashes the explicit theorem projection and elaborated declaration type. The \
checker receipt carries manifest/card/claim digests, \
declaration name, elaborated-type digest, complete import/environment digest, kernel/ \
toolchain digest, proof-term digest, exact axiom-policy digest, and complete transitive axiom \
closure; any mismatch or axiom outside the exact policy is IntegrityFailed, and runtime premise instantiation is a separate \
content-addressed receipt. The numerical refutation claim \
i03-topology-force-counterexample-search instead binds its canonical claim_digest to the \
independent exhaustive-search/minimization oracle and receives no Lean declaration by alias. \
NATURALITY: representations are variational sheaves over \
pinned finite relative-complex covers. A differentiable restriction-compatible chain \
equivalence P_x,Q_x with coherent descent and certified inverse homotopies receives \
variational authority only through either an equivariant stationary condensation with \
unique smooth vertical stationary lift, coercive/invertible vertical Hessian, and exact \
reduced-action/Schur identity, or complete descending filtered pronilpotent cyclic \
L_infinity/BV data with filtration-preserving morphism/homotopy inverse, fixed completion \
and pairing, master equation, terminating or convergence-certified higher series, and a \
proof of equivalence of the explicitly selected Deligne-Getzler or Hinich derived critical \
groupoids. A finite- \
polynomial alternative supplies checked mutually inverse groupoid functors and natural \
equivalences rather than invoking generic quasi-isomorphism invariance. Each route transports the complete gauge- \
reduced ensemble action including every shape derivative. The declaration derives an \
equivalence of stationary solution groupoids and f_A=phi^*f_B for the \
generalized-force one-form; neither solution nor force/QoI equality is a premise, and \
vector covariance additionally requires a natural metric/Riesz map. REFINEMENT: the \
declaration defines the pinned configuration transfer T_x and gives a signed averaged- \
adjoint identity T_x^*f_h'-f_h=D_primal+D_dual+D_ensemble+D_QoI_transfer+ \
D_homotopy_descent+D_Hodge+D_material+D_geom+D_load+D_trace+D_port+ \
D_gauge+D_constraint+D_quad+D_solver+D_state+D_perturbation+D_endpoint+D_nonlinear. \
Its joint enclosure preserves dependency and cancellation. Componentwise zero is \
sufficient, not necessary; exact covariance requires a proof that the total signed sum is \
zero. EVENT: a pinned regularized stratified cobordism passes exact finite-epsilon charge \
and energy-work ledgers to a declared distribution/Radon-measure limit with disjoint \
ownership and pairwise cancellation of internal transfers. Finite jump, work atom, or \
time-parameterized mechanical impulse requires separate tightness, integrability, and \
renormalization premises; divergence remains explicit. EXPORTS: \
proofs/i03/TopologyForce.lean::representationNaturalityForcePullback, \
proofs/i03/TopologyForce.lean::refinementForceDefectIdentity, and \
proofs/i03/TopologyEventJump.lean::distributionalEventBalance. NUMERICAL ACCEPTANCE: for \
pulled force-covector discrepancy delta_f and certified joint enclosure E, take the \
supremum over the frozen unit-normalized admissible variation set of the distance between \
dual-work pairings divided by 1e-15 J+1e-7 A_dualwork. SCALE_BINDING fixes \
A_dualwork=max(1e-15 J,sup_{v in Vunit}|f_h,oracle[v]|, \
sup_{v in Vunit}|T_x^*f_h',oracle[v]|) from fixture-only independent reference actions, \
with exact norm/pairing, positive floor, formula id, unit id, IEEE bits, and source digest \
before any candidate result. It \
exists for scalar-coenergy and complete constrained-action routes alike; score<=1, every nonvacuity floor holds, \
and verified in-domain counterexample count equals zero. Missing premises yield Unknown; \
binding, ownership, or integrity defects yield IntegrityFailed.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: CAMPAIGN_POLICY_FIXTURE,
            source: FixtureSource::AuthoredSpec {
                spec: "I03_CAMPAIGN_POLICY_V1\n\
EXECUTION_DISPOSITIONS=Completed|Cancelled|TimedOut|BudgetExhausted|InfrastructureFailed\n\
PREDICATE_OUTCOMES=Satisfied|Violated|Indeterminate\n\
CLAIM_ADJUDICATIONS=Supported|Failed|Refuted|Unknown\n\
EVIDENCE_COMPLETENESS=CompleteEvidence|PartialEvidence|NoEvidence\n\
EVIDENCE_INTEGRITY=IntegrityVerified|IntegrityFailed\n\
OPERATIONAL_SUPPORT=SupportedOperation|UnsupportedOperation\n\
OBSERVABLE_DISPOSITIONS=FiniteObservable|DeclaredDivergent\n\
PROMOTION_EFFECTS=Promotes|BlocksPromotion|NoPromotionAuthority\n\
NO_COLLAPSE=execution disposition, requested predicate outcome, claim adjudication, evidence completeness, evidence integrity, operational support, observable disposition, and promotion effect remain separate; DeclaredDivergent is an owned observable/payload disposition while claim adjudication remains Unknown unless the requested predicate is divergence\n\
LOGGING=fs-obs schema-versioned bounded JSONL; stable case/leaf/claim/fixture/oracle/attempt ids; exact units,seeds,budgets,versions,capabilities; I03_FIXTURE_V1 hashes the exact UTF-8 alias with BLAKE3 derive-key domain org.frankensim.i03.fixture-stream.v1 and interprets digest bytes [0,8) little-endian as the Philox 4x32-10 64-bit key; counter words c0,c1 are case-index low/high32 and c2,c3 are output-block low/high32; ordinary draws occupy factor_id=0, encode floor(j/4) in output-block bits 0..47, and use lane j mod 4; declared factor substreams require factor_id in 1..=65535, encode it in output-block bits 48..63 and floor(j/4) in bits 0..47, and use factor-local lane j mod 4, so rejection cannot shift another factor; every local block must be <2^48, and factor-zero mixing, overflow, or counter aliasing is IntegrityFailed; an n-option choice requires 1<=n<=2^32, computes t=floor(2^32/n)*n in exact/u64 arithmetic, accepts u64(x)<t, returns u64(x) mod n, and advances j on rejection; n=2^32 has t=2^32, accepts every x, and returns x, while n=0 or wrapped/u32 threshold arithmetic is IntegrityFailed; canonical floats use (x+0.5)/2^32; the finite alias table must have no derived-key collision before generation; options and factor ids retain written order and inputs serialize canonical little-endian; human/JSONL semantic parity; large fields content-addressed; material/fixture/license secrets redacted by policy\n\
HELDOUT_COMMIT_REVEAL=the named roster {i03-stat-custodian-a,i03-stat-custodian-b,i03-stat-custodian-c}, custodian/governance key digests, 3-of-3 rule, exact ten-record tagged-union/join-DAG I03StatHoldoutTranscriptV1 stages {CANDIDATE_FIXED,COMMIT,COMMIT_SET_FIXED,BEACON,REVEAL,FINAL}, strict RFC8032 SHA-512 Ed25519 authentication with canonical encodings, prime-subgroup checks and uncofactored verification equation, and BLAKE3 record/root/key/commitment domains are admission-bound; BEACON binds source, round, exact 32-byte value, and source-authentication digest while FINAL binds the three reveal digests, nonrecursive pre-FINAL root over signed records 0..8, and sampler-output root; candidate/model/toolchain, the finite atom law with probabilities n_i/2^256, and all campaign semantics become irrevocable before the custodians generate/commit entropy; commitments are sibling records linked to CANDIDATE_FIXED and use a sealed simultaneous all-three-fixed phase that withholds peer bytes until every signed COMMIT record is immutable, then COMMIT_SET_FIXED joins their roster-ordered digests before the beacon and signed REVEAL records; each opening is exactly 1024 indexed 256-bit blocks in lot-major offset [32*(2*l+r),32*(2*l+r+1)) for r=LOT=0/CANDIDATE_SEED=1 and l=0..511, and exact authority is conditional on an independently audited at-least-one-honest custodian whose complete vector is information-theoretically IID uniform and independent of the already-fixed candidate and every adversarial mask/commitment; the transcript remains outside the candidate sandbox and exactness does not rely on computational hiding; the three per-index custodian blocks are combined only by frozen coordinate-wise XOR, with exact 257-bit cumulative endpoints and unique half-open [c_i,c_(i+1)) dyadic-interval sampling and no hash/XOF/PRG/rejection expansion minting lot independence; the separately bound 256-bit beacon is never a sampling input and only selects a cyclic shift/direction order challenge, so its distribution cannot perturb the lot law; one identical candidate map sees only its current lot features, disjoint candidate-seed block, and frozen constants, while ordinals, ids, order, shard/worker, nonce/beacon, other blocks, and logging metadata remain outside its sandbox; secret order permutation and fresh external-id reassignment must leave per-lot predictions byte-identical; missing/late/duplicate reveal, wrong length/key, beacon substitution, abort, retry, or resampling is IntegrityFailed with no replacement holdout; access log exposes no labels, nonce, or raw draws until one-shot final submission, then releases retained replay bytes; any locally regenerable public statistical max-holdout has no untouched or IID authority\n\
THEOREM_AXIOMS=policy i03.lean-axioms.v1 permits exactly {propext,Quot.sound,Classical.choice}; its digest uses BLAKE3 derive-key domain org.frankensim.i03.lean-axioms.v1 over lexically sorted u32-little-endian-length-framed UTF-8 names; every theorem receipt binds this digest and the complete transitive declaration/environment axiom closure; sorryAx, custom or theorem-equivalent postulates, any extra axiom, and native-oracle proof authority outside the admitted kernel/checker TCB are IntegrityFailed\n\
FORMAL_PROJECTION=manifest version 1 freezes ambitious theorem targets but its prose cards mint no theorem promotion authority; before proof/candidate bytes exist, a FrozenManifest::amend successor must freeze canonical machine proposition AST bytes, all referenced definition bytes/digests, a total formal-hypothesis-to-runtime-evidence schema, and the exact source/version digest of a deterministic total AST-to-Lean translator, with structural generated/elaborated-type round trips; later admission receipts cannot choose or weaken the proposition\n\
M0_FORMALIZATION=manifest version 1 freezes the exact 16x16x16x16, N_micro=65536 exhaustive-search target but its prose mints no exhaustive authority; before search/candidate bytes exist, a FrozenManifest::amend successor must freeze the total decorated BaseRecord schema, canonical field encodings/domains, all validity and 16 stratum predicates, exact semantics/verdicts for 16 tags, explicit typed values of 16 parameter tuples, event primitive encodings, total enumeration/order/exclusion and rank/unrank/shard algorithms, parser/checker source digests, independent decoder and bijection proofs, cost preflight, and Merkle root; later admission receipts cannot choose missing grammar semantics\n\
SCALE_BINDING=every unit-bearing residual/error score not separately assigned an exact dimensionless formula has form ||Delta||/(a_abs+a_rel*S), with unit-bearing a_abs>0, dimensionless a_rel>0, finite S>0, norm, support, aggregation, formula id, unit id, IEEE bits, and source digest frozen in fixture bytes; every dimensionless ratio, empirical width/order statistic, and exact bit instead freezes its complete formula and all candidate-independent inputs; S and every candidate-independent exception input are computed/frozen before candidate execution solely from authored parameters and independent fixture/oracle reference data within declared upper/lower formulas; candidate outputs may enter only the frozen numerator/statistic/bit evaluation and cannot set or enlarge a reference scale, threshold, support, or aggregation; nonfinite, missing, or mismatched scales/inputs are IntegrityFailed\n\
RETENTION=promotion and independent manifest-adjudication receipts, theorem cards, refutations, minimized failures, checkpoint/fork lineage, amendment lineage, and integrity failures durable; smoke telemetry bounded/expirable; each non-success emits a replayable FailureBundle or an authenticated reason it cannot\n\
ACCESSIBILITY_AGENT_PARITY=every public campaign action and diagnostic reachable without pointer-only interaction; stable non-TTY JSONL surface carries the same choices, state, remediation, and replay identity as the human surface\n\
PERFORMANCE=smoke/core/max wall,memory,cancellation-latency envelopes come from FiveExplicits; regressions are evidence failures, never silently widened thresholds\n\
PROMOTION=baseline claims consumed only by I03.G4 core adjudication after independent G2 reproduction and G3 falsification; maximal claims consumed only by I03.G7 after independent G5 reproduction and G6 red-team; each heldout fixture has one named stage-local consumer and any premature/cross-stage read is IntegrityFailed; missing/stale/waived/integrity-failed evidence cannot promote\n\
AMENDMENT=ManifestDraft.version is the sole machine-interpreted instance-revision authority; successor version only; exact affected descendants invalidated; prior receipts remain bound to their original manifest and unchanged component evidence may be rebound only through authenticated amendment lineage plus identical component digest\n\
LEAF_REQUIREMENT=every I03 obligation row references this exact fixture id; declares smoke/core/max tier, DSR lane, events, replay command, request-drain-finalize plus checkpoint/resume/fork schedule, determinism matrix, cross-crate/IR/API roundtrip, and independent adjudication receipt; and covers claims that separately declare independent oracle and no-claim effect.",
            },
            partition: Partition::Development,
        },
    ]
}

#[allow(clippy::too_many_lines)]
fn i03_obligations() -> Vec<ObligationRow> {
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
            leaf: "i03-feec-cutfem-electrostatics",
            claims_covered: &[
                "i03-electrostatic-exact-sequence",
                "i03-electrostatic-convergence-gauss",
            ],
            unit_cases: UNIT_CASES,
            g0: "generators: parallel-plate/layered scalar H0 and mixed H1 harmonic \
                 complexes with valid and omitted-gauge/period twins; predicates: exact \
                 relative complex, bounded commuting projection, coercive quotient, \
                 regularity, cut/quadrature receipt; laws: d1*d0=0 exactly, correct H0 \
                 gauge, zero physical electrostatic H1 periods, gauge invariance, unit covariance, Gauss \
                 balance and directed convergence floors; shrinkers retain topology/cut/ \
                 regularity defect; cross-crate/IR/API roundtrip preserves all \
                 units and receipt ids; replay seeds per FiveExplicits",
            decks: &[
                "i03-parallel-plate-mms",
                "i03-coax-sphere-harmonic",
                CAMPAIGN_POLICY_FIXTURE,
            ],
            g3_relations: &[
                "rigid chart transform leaves normalized observables invariant",
                "under the pinned nested/stable/consistent asymptotic family the certified upper envelope converges at its directed floor, without demanding stepwise raw goal-error monotonicity",
                "a common gauge shift of every conductor and its reference changes potential but not field, charge, energy, or force",
            ],
            g4_schedule: "request-drain-finalize: inject cancellation before classification, after each assembly \
                          tile, during solve/finalize, and around content-addressed \
                          checkpoint save/resume/fork; request-drain-finalize publishes no \
                          partial field; resumed/forked deterministic state equals the \
                          uninterrupted prefix; corrupt gauge, topology, cut, material, \
                          checkpoint, and adjudication receipts independently; retain the \
                          minimized FailureBundle",
            g5_matrix: "threads {1,2,7} x shards {1,3} x deterministic scheduler \
                        permutations x ISA families {Apple-aarch64,x86_64}; bitwise \
                        comparison is only within an identical ISA fingerprint; exact \
                        topology, reduction tree, and content identity must match",
            entry_point: "scripts/e2e/leapfrog/i03_electrostatics.sh",
            tier: CampaignTier::Core,
            dsr_lane: "env FRANKENSIM_VMANIFEST_LEAF=i03-feec-cutfem-electrostatics dsr quality --tool frankensim",
            obs_events: &[
                "electrostatic.admission",
                "electrostatic.gauss_balance",
                "electrostatic.convergence",
                "electrostatic.cancelled",
                "execution.cancelled",
                "checkpoint.saved",
                "checkpoint.resumed",
                "checkpoint.forked",
                "evidence.adjudication_receipt",
                "evidence.failure_bundle",
            ],
            replay_command: "scripts/e2e/leapfrog/i03_electrostatics.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i03-conductors-capacitance",
            claims_covered: &["i03-floating-conductor-capacitance"],
            unit_cases: UNIT_CASES,
            g0: "generators: fixed/floating/charge-constrained conductor ensembles and \
                 closed-all-conductor versus fixed-ground/infinity formulations, q0+C*V, \
                 reciprocal/nonreciprocal and disconnected negative twins; predicates: \
                 explicit reference model, connected source-free theorem domain, gauge, \
                 terminal ownership, admitted material; laws: charge conservation; \
                 PSD/kernel span{1}/zero row sums only for the closed indefinite matrix; \
                 structure on C not q0; reduced SPD under coercivity; cross-crate/IR/API \
                 roundtrip preserves terminal authority; shrinkers preserve the boundary/reference defect",
            decks: &[
                "i03-floating-conductors",
                "i03-floating-conductors-core-holdout",
                "i03-coax-sphere-harmonic",
                CAMPAIGN_POLICY_FIXTURE,
            ],
            g3_relations: &[
                "terminal permutation conjugates the capacitance matrix",
                "a common shift of all conductors including the enclosure/reference changes no terminal charge",
                "reference-terminal choice preserves reduced physical predictions",
            ],
            g4_schedule: "request-drain-finalize: cancel after each terminal solve/matrix column and around \
                          checkpoint save/resume/fork; drain publishes the complete \
                          authenticated matrix or none; resumed column accumulation equals \
                          uninterrupted order; inject duplicated owners, wrong reference \
                          model, stale gauge/holdout/checkpoint, and non-finite flux; retain \
                          the first minimized structural counterexample",
            g5_matrix: "terminal orders forward/reverse x threads {1,2,7} x shards \
                        {1,4} x deterministic mode x ISA families \
                        {Apple-aarch64,x86_64}; bitwise comparison is only within an \
                        identical ISA fingerprint; bit-identical canonical matrix and \
                        terminal receipt identity",
            entry_point: "scripts/e2e/leapfrog/i03_capacitance.sh",
            tier: CampaignTier::Core,
            dsr_lane: "env FRANKENSIM_VMANIFEST_LEAF=i03-conductors-capacitance dsr quality --tool frankensim",
            obs_events: &[
                "terminal.charge_balance",
                "capacitance.structure",
                "capacitance.oracle",
                "capacitance.cancelled",
                "execution.cancelled",
                "checkpoint.saved",
                "checkpoint.resumed",
                "checkpoint.forked",
                "evidence.adjudication_receipt",
                "evidence.failure_bundle",
            ],
            replay_command: "scripts/e2e/leapfrog/i03_capacitance.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i03-conduction-dielectric-laws",
            claims_covered: &[
                "i03-steady-conduction-conservation",
                "i03-declared-dielectric-laws",
            ],
            unit_cases: UNIT_CASES,
            g0: "generators: resistor MMS and dielectric cards with passive, active, \
                 reciprocal, gyrotropic, incomplete-history, empirical, and out-of-regime \
                 twins; predicates: H(div)/dual conservative flux, Neumann compatibility, \
                 class-specific LTI/storage/IQC/empirical card, state/history, units and \
                 domain; laws: local/global charge, nonnegative symmetric-part loss, \
                 exp(+i omega t) uses causal s=0+i omega, positive-real Yp only for LTI, \
                 nonlinear storage inequality, declared- \
                 only reciprocity, unit covariance, and cross-crate/IR/API roundtrip; \
                 shrinkers preserve the first card/flux/balance defect",
            decks: &[
                "i03-resistor-conduction-mms",
                "i03-dispersive-dielectric-cards",
                CAMPAIGN_POLICY_FIXTURE,
            ],
            g3_relations: &[
                "series/parallel reduction preserves terminal response only for the declared linear lumped/network subfixtures",
                "for autonomous laws, translating forcing, time origin, and complete material history together preserves response",
                "reciprocal-axis rotation covariantly rotates flux without changing loss",
            ],
            g4_schedule: "request-drain-finalize: cancel during card admission/history initialization, each \
                          conduction tile/loss accumulation, and checkpoint save/resume/ \
                          fork; drain leaves no half-admitted state; resumed/forked material \
                          history and flux are deterministic; inject missing history, wrong \
                          checker class, negative passive loss, non-finite parameters, port \
                          sign flips, and corrupt checkpoints; retain bounded redacted \
                          material evidence",
            g5_matrix: "material-card order permutations x threads {1,2,7} x shards \
                        {1,3} x deterministic mode x ISA families \
                        {Apple-aarch64,x86_64}; bitwise comparison is only within an \
                        identical ISA fingerprint; identical admission, state, current, \
                        loss, and refusal receipts",
            entry_point: "scripts/e2e/leapfrog/i03_material_conduction.sh",
            tier: CampaignTier::Core,
            dsr_lane: "env FRANKENSIM_VMANIFEST_LEAF=i03-conduction-dielectric-laws dsr quality --tool frankensim",
            obs_events: &[
                "material.regime",
                "material.passivity",
                "conduction.balance",
                "conduction.loss",
                "material_conduction.cancelled",
                "execution.cancelled",
                "checkpoint.saved",
                "checkpoint.resumed",
                "checkpoint.forked",
                "evidence.adjudication_receipt",
                "evidence.failure_bundle",
            ],
            replay_command: "scripts/e2e/leapfrog/i03_material_conduction.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i03-field-circuit-force-adjoint",
            claims_covered: &[
                "i03-eqs-regime-routing",
                "i03-field-circuit-power-closure",
                "i03-force-adjoint-held-variable-closure",
            ],
            unit_cases: UNIT_CASES,
            g0: "generators: EQS/full-wave boundary sweeps, descriptor-circuit transients, \
                 and fixed-V/fixed-Q/mixed/floating force paths with sign/control/branch \
                 twins; predicates: quantitative EQS card, stationary total current or \
                 moving-port pullback/GCL receipt, consistent DAE, disjoint term owners, \
                 storage realization, and whole-composition holomorphic eligibility or \
                 interval finite-difference route; scalar coenergy additionally requires a \
                 certified integrable/path-independent terminal-charge one-form, otherwise \
                 the complete constrained action is mandatory; laws: exact regime routing, terminal- \
                 charge continuity, zero-velocity current reduction, subsystem/combined \
                 work identity or inequality, complete coenergy/Legendre functional, \
                 generalized virtual-work/adjoint pairing, and cross-crate/ \
                 IR/API roundtrip; shrinkers preserve the first monitor/port/ \
                 ensemble defect",
            decks: &[
                "i03-eqs-regime-boundaries",
                "i03-eqs-regime-core-holdout",
                "i03-field-circuit-transients",
                "i03-force-held-variable-benchmarks",
                "i03-force-held-variable-core-holdout",
                CAMPAIGN_POLICY_FIXTURE,
            ],
            g3_relations: &[
                "port orientation reversal flips effort-flow signs but preserves total power",
                "fixed-V, fixed-Q, mixed, and floating ensembles agree only after applying their exact constrained Legendre functionals; no generic opposite-force-sign rule is assumed",
                "an orientation-preserving reparameterization of the same smooth stable equilibrium branch leaves integrated virtual work invariant",
                "deforming a voltage-space integration path preserves scalar coenergy only for an integrable card; a nonzero closed-loop charge integral routes to the complete constrained action",
            ],
            g4_schedule: "request-drain-finalize: cancel/fault at regime admission, DAE causalization, field solve, \
                          port exchange, energy ledger, adjoint reverse tiles, final force \
                          receipt, and checkpoint save/resume/fork; losers drain; resumed/ \
                          forked DAE+material+adjoint state matches uninterrupted replay; \
                          inconsistent state, duplicate work ownership, event/topology/ \
                          branch change, or non-holomorphic complex-step must refuse; \
                          retain exact closure ledger and minimized witness",
            g5_matrix: "threads {1,2,7} x shards {1,3} x port insertion permutations x \
                        deterministic forward/adjoint replay x ISA families \
                        {Apple-aarch64,x86_64}; bitwise comparison is only within an \
                        identical ISA fingerprint; bit-identical plan, field, port, \
                        energy, force, and gradient receipts",
            entry_point: "scripts/e2e/leapfrog/i03_field_circuit_force.sh",
            tier: CampaignTier::Core,
            dsr_lane: "env FRANKENSIM_VMANIFEST_LEAF=i03-field-circuit-force-adjoint dsr quality --tool frankensim",
            obs_events: &[
                "field_circuit.causalized",
                "field_circuit.power_closure",
                "force.held_variable",
                "adjoint.identity",
                "eqs.regime",
                "field_circuit_force.cancelled",
                "execution.cancelled",
                "checkpoint.saved",
                "checkpoint.resumed",
                "checkpoint.forked",
                "evidence.adjudication_receipt",
                "evidence.failure_bundle",
            ],
            replay_command: "scripts/e2e/leapfrog/i03_field_circuit_force.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i03-breakdown-routing",
            claims_covered: &["i03-partial-discharge-breakdown-routing"],
            unit_cases: UNIT_CASES,
            g0: "generators: insulation adversaries with censored, biased, outlier, \
                 within-lot correlated, holdout-leak, wide-interval, and out-of-domain \
                 twins; predicates: geometry/gas/waveform/measurement/lot/event/censoring \
                 declaration complete, candidate-first commit/reveal intact, and an \
                 independently governed audited receipt establishes heldout lot clusters \
                 IID from the pinned population; laws: exact \
                 Bonferroni-simultaneous familywise-99% one-sided per-family/per-endpoint \
                 cluster bounds at alpha=0.001, each >=0.90, conservative censoring, \
                 factor-local Philox lanes with explicit repeat/family/endpoint censoring ownership \
                 only for deterministic within-lot expansion, never for the IID premise, \
                 cellwise median/p90 width caps, explicit Supported/Unknown/Unsupported, synthetic-versus-\
                 experimental authority and cross-crate/IR/API roundtrip; shrinkers \
                 retain censoring/dependence/sharpness defect",
            decks: &[
                "i03-insulation-adversaries",
                "i03-insulation-adversaries-max-holdout",
                "i03-external-hv-industrial-pack",
                CAMPAIGN_POLICY_FIXTURE,
            ],
            g3_relations: &[
                "unit-rescaled equivalent waveforms preserve normalized onset coverage",
                "permuting the eight within-lot repeats or the IID lot-cluster order leaves each preregistered worst-repeat cell outcome and exact confidence bound unchanged",
                "uniformly widening intervals cannot rescue a sharpness failure even if raw coverage increases",
                "rejection or censoring in one factor-local substream leaves every other factor and cell stream bit-identical",
            ],
            g4_schedule: "request-drain-finalize: cancel under adversarial stopping, during label/cluster \
                          formation, before confidence/sharpness finalization, and around \
                          checkpoint save/resume/fork; drain retains censoring and never \
                          publishes a partial physical claim; inject missing lot/event \
                          definitions, correlated clusters, stale/premature holdout, \
                          oracle misuse, and budget exhaustion; retain refutations and \
                          minimized integrity failures durably",
            g5_matrix: "threads {1,2,7} x shards {1,4} x deterministic stochastic streams \
                        x held-out replay x ISA families {Apple-aarch64,x86_64}; bitwise \
                        comparison is only within an identical ISA fingerprint; exact \
                        sample/cluster membership, labels, route, coverage bound, \
                        sharpness, and evidence identity",
            entry_point: "scripts/e2e/leapfrog/i03_breakdown.sh",
            tier: CampaignTier::Max,
            dsr_lane: "env FRANKENSIM_VMANIFEST_LEAF=i03-breakdown-routing dsr quality --tool frankensim",
            obs_events: &[
                "insulation.regime",
                "discharge.coverage",
                "discharge.sharpness",
                "breakdown.cancelled",
                "execution.cancelled",
                "checkpoint.saved",
                "checkpoint.resumed",
                "checkpoint.forked",
                "evidence.adjudication_receipt",
                "evidence.failure_bundle",
            ],
            replay_command: "scripts/e2e/leapfrog/i03_breakdown.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i03-space-charge-aging",
            claims_covered: &["i03-space-charge-aging-singular-routing"],
            unit_cases: UNIT_CASES,
            g0: "generators: conservative mobile/trapped free species, derived polarization \
                 bound volume/interface charge/current, separately oriented carrier- \
                 transfer/electrode-free-charge/total-current ledgers, charge-null \
                 reactions, thermodynamic and empirical aging cards, and scaled conductor/ \
                 insulator/index-changing limit twins; predicates: stoichiometry, all \
                 owners, positivity/capacity/history, card authority, scaling, limiting \
                 problem, topology/norm and uniform estimate complete; laws: q_bulk^T S_bulk=0, \
                 support-dimension typing with distinct bulk/surface/electrode/port scores, \
                 free-charge continuity, div D=rho_f, \
                 polarization-derived total-charge equivalence without double counting, \
                 fixed/moving-interface weak balance, electrode/circuit reconciliation, state \
                 invariance, declared storage/dissipation only, convergence \
                 or explicit refusal, and cross-crate/IR/API roundtrip; shrinkers \
                 retain the first species, authority, index, or limit defect",
            decks: &[
                "i03-space-charge-aging-limits",
                "i03-space-charge-aging-max-holdout",
                CAMPAIGN_POLICY_FIXTURE,
            ],
            g3_relations: &[
                "refining a charge-source partition preserves the stoichiometric global ledger",
                "a units-and-state rescaling derived from the pinned nondimensionalization preserves the limiting route",
                "approaching a pinned singular family converges in its declared topology or crosses its explicit refusal boundary",
                "Eulerian and ALE moving-interface charge ledgers agree only after the pinned spacetime pullback and geometric-conservation term",
            ],
            g4_schedule: "request-drain-finalize: cancel mid-transport, reaction, material/damage update, limit \
                          estimator, and checkpoint save/resume/fork; drain publishes no \
                          partial service life; resumed/forked species/history/ledger state \
                          equals uninterrupted replay; inject missing owners/history, \
                          capacity violations, empirical energy overclaim, wrong limiting \
                          problem, index change, stale holdout, and budget exhaustion",
            g5_matrix: "threads {1,2,7} x shards {1,4} x species/card permutations x \
                        deterministic streams x held-out replay x ISA families \
                        {Apple-aarch64,x86_64}; bitwise comparison is only within an \
                        identical ISA fingerprint; exact state, ledger, route, enclosure, \
                        checkpoint, and evidence identity",
            entry_point: "scripts/e2e/leapfrog/i03_space_charge_aging.sh",
            tier: CampaignTier::Max,
            dsr_lane: "env FRANKENSIM_VMANIFEST_LEAF=i03-space-charge-aging dsr quality --tool frankensim",
            obs_events: &[
                "transport.charge_balance",
                "transport.state_admissibility",
                "aging.authority",
                "singular_limit.route",
                "space_charge_aging.cancelled",
                "execution.cancelled",
                "checkpoint.saved",
                "checkpoint.resumed",
                "checkpoint.forked",
                "evidence.adjudication_receipt",
                "evidence.failure_bundle",
            ],
            replay_command: "scripts/e2e/leapfrog/i03_space_charge_aging.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i03-electrostriction-theorem",
            claims_covered: &[
                "i03-electrostriction-energy-interface-closure",
                "i03-electrostriction-interface-theorem",
            ],
            unit_cases: UNIT_CASES,
            g0: "generators: objective-invariant energies and orientation-preserving deformations, material states, \
                 electrical ensembles, interfaces, contrast limits, regular/unstable \
                 branches, injective/interpenetrating and theorem-valid/invalid twins; \
                 predicates: exact canonical-claim-digest/theorem/declaration binding, independent state+Legendre transform, \
                 interval-global det(F) floor, boundary one-to-one/Ciarlet-Necas \
                 injectivity with Route H global and Route CN only a.e. authority; \
                 Route H/CN invertibility is orthogonal to REDUCED/MIXED equilibrium \
                 stability and every cross-product is typed; exact reduced-Schur \
                 coercivity or mixed bordered-Hessian inf-sup, frame/interface, fixed topology/ \
                 branch, uniform singular hypotheses; laws: formal checker bit separate \
                 from normalized numerical frame/energy/total stress/traction/virtual- \
                 work/adjoint closure and cross-crate/IR/API roundtrip; \
                 shrinkers preserve the failed premise/interface patch",
            decks: &[
                "i03-electrostriction-finite-strain",
                "i03-electrostriction-max-holdout",
                "i03-electrostriction-theorem-card",
                CAMPAIGN_POLICY_FIXTURE,
            ],
            g3_relations: &[
                "superposed rigid motion preserves energy and covariantly rotates stress",
                "interface orientation reversal swaps side/normal signs without changing balance",
                "admitted contrast-limit sequences preserve the theorem error enclosure",
                "Route H global-homeomorphism authority may be downcast to Route CN a.e.-injectivity authority, never promoted in reverse",
            ],
            g4_schedule: "request-drain-finalize: cancel formal export/checking, interval subdivision, nonlinear \
                          solve, adjoint, receipt finalization, and checkpoint save/resume/ \
                          fork; corrupt each theorem premise, ensemble, branch, material/ \
                          interface convention, proof term, and checkpoint; drain publishes \
                          no theorem color before checker+adjudicator completion; resumed \
                          state is equivalent; proof rejection and minimized traction \
                          counterexamples are durable",
            g5_matrix: "threads {1,2,7} x interval shards {1,4} x chart/material order \
                        permutations x deterministic mode x ISA families \
                        {Apple-aarch64,x86_64}; bitwise comparison is only within an \
                        identical ISA fingerprint; identical theorem bytes, checker \
                        receipt, enclosures, and numerical witness",
            entry_point: "scripts/e2e/leapfrog/i03_electrostriction.sh",
            tier: CampaignTier::Max,
            dsr_lane: "env FRANKENSIM_VMANIFEST_LEAF=i03-electrostriction-theorem dsr quality --tool frankensim",
            obs_events: &[
                "theorem.exported",
                "theorem.checked",
                "electrostriction.energy_closure",
                "electrostriction.traction",
                "electrostriction.cancelled",
                "execution.cancelled",
                "checkpoint.saved",
                "checkpoint.resumed",
                "checkpoint.forked",
                "evidence.adjudication_receipt",
                "evidence.failure_bundle",
            ],
            replay_command: "scripts/e2e/leapfrog/i03_electrostriction.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i03-topology-force-theorem-falsifier",
            claims_covered: &[
                "i03-cohomology-force-naturality-theorem",
                "i03-refinement-force-defect-enclosure",
                "i03-topology-event-jump-theorem",
                "i03-topology-force-counterexample-search",
            ],
            unit_cases: UNIT_CASES,
            g0: "generators: canonical finite complexes, relative subcomplexes, \
                 variational-sheaf covers/descent, weighted actions/Hodge/material/load/ \
                 trace/port diagrams, harmonic coordinates, chain-equivalence homotopies, \
                 chart/refinement maps, signed defect decompositions, cobordisms, gap/contact \
                 limits, charge/circuit/source laws, and valid/invalid theorem twins; \
                 predicates: pre-candidate M0_FORMALIZATION successor gate, exact bound \
                 theorem-domain membership, full decorated-object canonicalization, \
                 cardinality/rank-unrank/budget proof for the exhaustive \
                 microgrammar, stationary condensation or complete filtered pronilpotent \
                 cyclic-L_infinity/BV variational-transfer receipt with terminating/ \
                 convergent higher series and derived-groupoid equivalence, and nonvacuity \
                 floors independently checked; laws: variational-sheaf descent and complete-action \
                 naturality (not chain/cohomology equivalence alone), derived equilibrium-groupoid equivalence \
                 and generalized-force-one-form pullback, metric-conditioned vector \
                 covariance, contragredient harmonic coordinates, signed averaged-adjoint \
                 identity under the exact T_x pullback with D_ensemble plus QoI-transfer/ \
                 homotopy/descent and dependency-preserving \
                 cancellation, typed distributional event ownership/weak balance, \
                 canonical claim_digest/declaration/proof binding, proof-kernel/adjudicator agreement, cross-crate/ \
                 IR/API roundtrip, and counterexample minimization; shrinkers retain the \
                 in-domain violation",
            decks: &[
                "i03-topology-force-adversaries",
                "i03-topology-force-max-holdout",
                "i03-topology-force-theorem-card",
                "i03-coax-sphere-harmonic",
                CAMPAIGN_POLICY_FIXTURE,
            ],
            g3_relations: &[
                "a complete variational-sheaf transfer carrying either the stationary-condensation or cyclic-L_infinity/BV critical-groupoid receipt pulls the generalized-force one-form back exactly; bare chain/cohomology equivalence is an invalid negative twin",
                "cohomology basis changes with contragredient coordinates and true gauge or compensated exact/potential changes that leave the total physical cochain identical preserve the observable including Gram/Hodge derivatives; an uncompensated exact-cochain shift is a negative twin",
                "general refinement satisfies its signed dependency-preserving enclosure; componentwise-zero defects are sufficient, while a formally certified nontrivial cancellation may also yield exact covariance",
                "a QoI-definition change or square commuting only up to homotopy changes its explicitly owned transfer/homotopy defect; omitting either is an invalid negative twin",
                "a premise-satisfying regularized cobordism preserves typed distributional charge and energy-work balances; a finite jump or impulse needs its additional tightness/integrability premises while class change alone predicts nothing",
                "canonical relabeling preserves theorem/falsifier verdict",
                "mutating any bound claim, card, theorem projection, declaration type, environment, axiom report, or proof identity invalidates theorem-receipt reuse",
            ],
            g4_schedule: "request-drain-finalize: cancel enumeration, optimization, proof checking, minimization, \
                          retention, and checkpoint save/resume/fork independently; \
                          request-drain-finalize loses no discovered candidate and resumed \
                          canonical search frontier equals uninterrupted replay; inject \
                          theorem-byte, premise, Hodge/load/trace map, event regularization, \
                          checker, canonicalization, checkpoint, adjudication, and artifact- \
                          integrity faults; every verified counterexample is durable and \
                          cannot be waived away",
            g5_matrix: "threads {1,2,7} x search shards {1,3,8} x enumeration orders x \
                        deterministic mode x ISA families {Apple-aarch64,x86_64}; bitwise \
                        comparison is only within an identical ISA fingerprint; exact \
                        candidate set, canonical minima, theorem bytes, checker receipts, \
                        and terminal disposition",
            entry_point: "scripts/e2e/leapfrog/i03_topology_force_theorem.sh",
            tier: CampaignTier::Max,
            dsr_lane: "env FRANKENSIM_VMANIFEST_LEAF=i03-topology-force-theorem-falsifier dsr quality --tool frankensim",
            obs_events: &[
                "cohomology.class",
                "topology_force.theorem_checked",
                "topology_force.candidate",
                "topology_force.counterexample",
                "topology_force.cancelled",
                "execution.cancelled",
                "checkpoint.saved",
                "checkpoint.resumed",
                "checkpoint.forked",
                "evidence.adjudication_receipt",
                "evidence.failure_bundle",
            ],
            replay_command: "scripts/e2e/leapfrog/i03_topology_force_theorem.sh --replay <artifact-id>",
        },
    ]
}

fn i03_waivers() -> Vec<Waiver> {
    vec![Waiver {
        subject: "i03-external-hv-industrial-pack",
        reason: "the external high-voltage partial-discharge/breakdown campaign pack \
                 is license-restricted and has not yet been registered with edition, \
                 raw-byte digest, normalization, censoring/dependence declaration, \
                 and redistribution policy",
        owner: "I03 implementation and V&V registry beads",
        predicate: "fs-vvreg admits the exact licensed bytes, edition, normalization, \
                    oracle ownership, directed acceptance arithmetic, and \
                    license/export controls",
        expiry: "before the first maximal I03 campaign submission; review at every \
                 Phase-2 close burst",
        promotion_effect: "baseline remains independently promotable, but partial- \
                           discharge/breakdown evidence cannot reach maximal promotion \
                           or cite the industrial pack while this waiver is live",
    }]
}
