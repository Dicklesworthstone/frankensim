//! The I09 (periodic-machine orbit/Floquet/adjoint/continuation
//! service) VerificationManifest draft (bead
//! frankensim-leapfrog-2026-program-i94v.4.1.7.1).
//!
//! Baseline lattice ([S]): typed orbit problems, cross-method orbit
//! agreement (shooting/collocation/harmonic balance), consistent-
//! tangent monodromy, descriptor-infinity separation, exact-discrete
//! adjoints, and a standing false-stability refutation lane for
//! nonnormal/near-defective monodromy. Maximal lattice ([F]/[M]):
//! certified multiplier intervals and phase-robust continuation with
//! branch-point receipts, then bounded-box global branch topology
//! with localized Unknown and hybrid/grazing continuation. A weaker
//! receipt closes its own element and is never relabeled as the
//! stronger theorem; version-1 prose cards mint no proof authority.

use crate::{
    Ambition, CampaignTier, ClaimPolarity, ClaimSpec, FiveExplicits, FixturePin, FixtureSource,
    GauntletTier, ManifestDraft, ObligationRow, OracleRoute, Partition, ToleranceSemantics, Waiver,
};

/// Build the I09 draft. Consumers freeze it themselves; the conformance
/// battery proves it freezes.
#[must_use]
pub fn i09_draft() -> ManifestDraft {
    ManifestDraft {
        initiative: "I09",
        title: "Periodic-machine orbit gate: typed orbit problems to cross-verified \
                orbits, consistent-tangent monodromy, Floquet classification with \
                certified multiplier intervals, exact-discrete adjoints, and \
                receipt-carrying continuation",
        version: 1,
        explicits: FiveExplicits {
            units: "SI base units throughout; periods in seconds, multipliers and \
                    phase conditions dimensionless (unit '1'); exact bitwise/boolean \
                    verdicts use unit 'bit'; residual norms in the per-fixture \
                    declared unit",
            seeds: "Philox 4x32-10 counter streams keyed 'i09/<fixture-id>/<case-index>'; \
                    development indices 0..=16383; core held-out indices 65536..=81919; \
                    maximal held-out indices 131072..=147455 (disjoint by construction, \
                    split frozen here; per-holdout subranges pinned in the fixture \
                    specs)",
            budgets: "smoke tier <= 60 s on one host; core tier <= 30 min; max tier <= 8 h \
                      on a quiet perf host; <= 16 GiB memory per lane; eigensolver \
                      service ticks and continuation steps carry declared budgets whose \
                      exhaustion is a typed resumable outcome; accuracy budgets are the \
                      per-claim tolerance fields",
            versions: "fs-vmanifest schema v2; toolchain pinned by \
                       rust-toolchain.toml; sibling pins by constellation.lock",
            capabilities: "no network; no FFI; deterministic mode mandatory for every G5 \
                           row; frontier/moonshot lanes stay behind feature flags; the \
                           degenerate-cluster and hybrid-continuation dependencies are \
                           waivered until their beads land",
        },
        claims: i09_claims(),
        fixtures: i09_fixtures(),
        obligations: i09_obligations(),
        waivers: i09_waivers(),
        amendment_rules: "Any change is a successor version through FrozenManifest::amend; \
                          the amendment record names every invalidated claim/obligation \
                          descendant; an amendment after campaign start invalidates the \
                          affected evidence, which must be re-earned; there is no in-place \
                          edit path in the type system",
    }
}

#[allow(clippy::too_many_lines)]
fn i09_claims() -> Vec<ClaimSpec> {
    vec![
        ClaimSpec {
            id: "i09-typed-orbit-problems",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "The typed orbit-problem IR (periodic, relative-periodic, and \
                        hybrid orbit families with phase conditions, symmetry \
                        declarations, and descriptor structure) round-trips its \
                        canonical encoding bit-stably and refuses every enumerated \
                        malformed class with a typed refusal naming the violated rule",
            hypotheses: &[
                "problems drawn from the pinned fixture families",
                "malformed classes are the enumerated ones: missing/inconsistent phase \
                 condition, non-positive declared period, symmetry group mismatch, \
                 descriptor mass matrix shape violations",
            ],
            qoi: "roundtrip_and_refusal_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G0,
            oracle: OracleRoute {
                identity: "independent decoder plus the fixture generator's own validity \
                           certificate (emitted during generation, not by the service)",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "first implementation bead of the I09 problem-IR leaf opens",
            kill: "any accepted malformed problem or refused valid problem returns the \
                   IR to design review with the disagreeing fixture as receipt",
            fallback: "explicit user-supplied normal form with per-field validation",
            no_claim: "no existence claim for orbits of an admitted problem; admission \
                       is typing, not analysis",
        },
        ClaimSpec {
            id: "i09-cross-method-orbit-agreement",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Shooting, collocation, and harmonic-balance solutions of the \
                        same admitted smooth orbit problem agree on orbit state and \
                        period within the declared band on every pinned smooth fixture, \
                        and each reports its own residual honestly",
            hypotheses: &[
                "smooth (non-hybrid) fixtures with isolated orbits",
                "each method converged under its own declared criterion; agreement is \
                 checked only between converged runs",
            ],
            qoi: "max_cross_method_orbit_discrepancy",
            unit: "1",
            tolerance: ToleranceSemantics::AbsRel {
                atol: 1e-8,
                rtol: 1e-6,
            },
            evidence_tier: GauntletTier::G2,
            oracle: OracleRoute {
                identity: "closed-form Mathieu/Hill references where available; \
                           otherwise the three methods are mutually independent \
                           routes and the discrepancy IS the QoI",
                independent: true,
                tcb_overlap: "methods share the fixture right-hand side only",
            },
            activation: "problem-IR leaf green at smoke tier",
            kill: "persistent cross-method disagreement beyond the band on a smooth \
                   fixture kills promotion of every involved method until root-caused",
            fallback: "single-method results labeled with method identity and residual, \
                       no cross-verified tag",
            no_claim: "no uniqueness claim; distinct orbits found by distinct methods \
                       are a reportable outcome, not a defect",
        },
        ClaimSpec {
            id: "i09-consistent-tangent-monodromy",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "The monodromy action assembled from the ACTUAL discrete \
                        tangent of the orbit discretization agrees with central finite \
                        differences of the period map on the pinned smooth fixtures",
            hypotheses: &[
                "smooth fixtures with isolated orbits and regular period maps",
                "FD stencil width inside the declared regularity neighborhood",
            ],
            qoi: "relative_monodromy_action_error",
            unit: "1",
            tolerance: ToleranceSemantics::Relative { rtol: 1e-6 },
            evidence_tier: GauntletTier::G1,
            oracle: OracleRoute {
                identity: "central finite differences of the composed period map \
                           (independent of the tangent assembly)",
                independent: true,
                tcb_overlap: "shares the primal integrator; tangent paths disjoint",
            },
            activation: "cross-method leaf green at smoke tier",
            kill: "gate failure on any smooth fixture kills the monodromy lane until \
                   rederived from the actual discrete residual",
            fallback: "FD monodromy labeled Estimated",
            no_claim: "no claim at non-smooth events; the hybrid rung enters through \
                       the waivered continuation lane",
        },
        ClaimSpec {
            id: "i09-descriptor-infinity-separation",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "For descriptor (DAE) orbit problems, spectral contributions of \
                        the constraint infinity modes are separated from finite Floquet \
                        multipliers with a typed classification, and no infinity mode \
                        is ever reported as a finite stability multiplier",
            hypotheses: &[
                "descriptor fixtures with declared index and regular pencil structure",
                "the declared index matches the fixture certificate",
            ],
            qoi: "classification_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G0,
            oracle: OracleRoute {
                identity: "analytic multiplier sets of the pinned constrained-DAE \
                           fixtures (constructed with known finite spectra)",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "monodromy leaf green at smoke tier",
            kill: "any infinity mode reported finite is a release blocker for the \
                   descriptor lane",
            fallback: "descriptor problems refuse Floquet classification and report \
                       raw pencil data",
            no_claim: "no claim for singular pencils or index drift at runtime; those \
                       refuse with the pencil diagnostic",
        },
        ClaimSpec {
            id: "i09-exact-discrete-adjoints",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Adjoint gradients of periodic-orbit QoIs (period, phase-locked \
                        state functionals) derived from the ACTUAL discrete orbit \
                        residual agree with central finite differences on the pinned \
                        smooth fixtures",
            hypotheses: &[
                "smooth fixtures; QoIs differentiable at the solution",
                "the orbit solve converged to its declared tolerance before \
                 differentiation",
            ],
            qoi: "relative_gradient_error",
            unit: "1",
            tolerance: ToleranceSemantics::Relative { rtol: 1e-6 },
            evidence_tier: GauntletTier::G1,
            oracle: OracleRoute {
                identity: "central finite differences of the full solve-then-evaluate \
                           pipeline (independent of the adjoint assembly)",
                independent: true,
                tcb_overlap: "shares the primal solver; adjoint path disjoint",
            },
            activation: "monodromy leaf green at smoke tier",
            kill: "gate failure kills the adjoint lane until rederived; no silent \
                   fallback to FD inside the adjoint API",
            fallback: "FD sensitivities labeled Estimated",
            no_claim: "no adjoint claim across bifurcation points or non-smooth \
                       events; those return typed no-claim records",
        },
        ClaimSpec {
            id: "i09-false-stability-falsifier",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Refutation,
            statement: "Adversarial nonnormal and near-defective monodromy families \
                        (large transient growth, multiplier clusters near the unit \
                        circle, Jordan-block limits) attempt to extract a certified \
                        asymptotic-stability verdict where Unknown or set-valued is \
                        the correct answer; any success refutes the classification \
                        lane's certificate discipline, and no measured transient is \
                        ever relabeled as asymptotic stability",
            hypotheses: &[
                "fixtures constructed so the true classification is provable by hand \
                 (spec carries the construction and proof sketch)",
                "the classifier runs at production settings",
            ],
            qoi: "false_certificate_count",
            unit: "1",
            tolerance: ToleranceSemantics::Interval { lo: 0.0, hi: 0.0 },
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "hand-proved analytic monodromy constructions with known \
                           defective/nonnormal structure",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "runs in every campaign tier from the first classification \
                         commit",
            kill: "this lane is never killed; it is the standing tripwire",
            fallback: "none: a nonzero count is a release blocker",
            no_claim: "absence of refutation is not completeness and cannot prove \
                       classifier soundness; the lane only falsifies",
        },
        ClaimSpec {
            id: "i09-certified-multiplier-intervals",
            ambition: Ambition::Frontier,
            polarity: ClaimPolarity::Affirmative,
            statement: "Floquet multiplier magnitudes are reported with certified \
                        existence intervals derived from residual-backed eigensolver \
                        receipts (Weyl-style containment through the fs-spectral \
                        service), and local bifurcation candidates (fold, flip, \
                        Neimark-Sacker) are flagged exactly when a certified interval \
                        crosses the unit circle; existence intervals never mint \
                        multiplicity",
            hypotheses: &[
                "monodromy actions from the consistent-tangent lane",
                "eigensolver residuals converged to the declared tolerance; intervals \
                 certify existence, not multiplicity",
            ],
            qoi: "interval_containment_and_flag_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G2,
            oracle: OracleRoute {
                identity: "analytic multiplier values of Mathieu/Hill and constructed \
                           rotor fixtures (closed-form or independently derived)",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "monodromy and descriptor leaves green at core tier",
            kill: "an analytic multiplier escaping every certified interval kills the \
                   classification lane",
            fallback: "point estimates labeled Estimated with residual context, no \
                       certified flags",
            no_claim: "no multiplicity claims (existence intervals only); degenerate \
                       cluster qualification is the waivered lane",
        },
        ClaimSpec {
            id: "i09-phase-robust-continuation",
            ambition: Ambition::Frontier,
            polarity: ClaimPolarity::Affirmative,
            statement: "Pseudo-arclength continuation of orbit branches under the \
                        declared phase conditions traverses the pinned parameter \
                        ranges, and every step carries a receipt naming step size, \
                        residuals, phase-condition conditioning, and detected \
                        branch-point candidates with typed classifications",
            hypotheses: &[
                "branches within the pinned parameter boxes of the fixture families",
                "step budgets per the explicits; exhaustion is a typed resumable \
                 outcome",
            ],
            qoi: "branch_traversal_receipt_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G2,
            oracle: OracleRoute {
                identity: "independently computed branch samples at pinned parameter \
                           values (fresh solves, not continuation-predicted)",
                independent: true,
                tcb_overlap: "shares the orbit solvers; continuation logic disjoint",
            },
            activation: "multiplier-interval leaf green at core tier",
            kill: "a silently skipped branch point or unreceipted step kills the \
                   continuation lane",
            fallback: "parameter sweeps of independent solves, no branch structure \
                       claims",
            no_claim: "no global completeness claim (that is the moonshot lane); \
                       between-sample branch structure is interpolated context only",
        },
        ClaimSpec {
            id: "i09-global-branch-topology",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Affirmative,
            statement: "Within declared bounded parameter boxes, the reported branch \
                        topology (components, folds, period-doubling cascades within \
                        resolution) is complete relative to the declared orbit class: \
                        every candidate region not covered by a certified branch or \
                        exclusion is a localized Unknown window with its budget \
                        receipt",
            hypotheses: &[
                "bounded parameter boxes and orbit class declared per fixture",
                "certified multiplier intervals from the frontier lane",
            ],
            qoi: "topology_completeness_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "exhaustive fine-grid enumeration on fixtures small enough \
                           to enumerate; hand-derived cascade structure for the \
                           period-doubling family",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "behind the i09-moonshot feature flag on a pre-proof successor \
                         version; frontier lanes green",
            kill: "any unreported uncovered region on an enumerable fixture kills the \
                   lane",
            fallback: "continuation-lane receipts without topology claims",
            no_claim: "no claim outside the declared boxes or orbit class; chaotic \
                       regions report as Unknown windows, not branches; version-1 \
                       prose mints no completeness authority",
        },
        ClaimSpec {
            id: "i09-hybrid-grazing-continuation",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Affirmative,
            statement: "Continuation through hybrid-event boundaries (impact orbits, \
                        grazing transitions) either carries a theorem-qualified \
                        crossing receipt (regular transversal case) or returns a \
                        set-valued/Unknown outcome at grazing and corner cases — \
                        never a silently smoothed branch",
            hypotheses: &[
                "hybrid fixtures from the pinned impact family with known grazing \
                 parameter values",
                "event accounting supplied by the I12 machinery for the prescribed \
                 rung; true-flow coverage is the waivered dependency",
            ],
            qoi: "hybrid_continuation_outcome_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "closed-form grazing parameter values of the impact-orbit \
                           family",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "behind the i09-moonshot feature flag on a pre-proof successor \
                         version; requires the hybrid waiver retired",
            kill: "a smoothed-through grazing point on a known fixture is Sev-0 for \
                   this lane",
            fallback: "continuation halts at hybrid boundaries with the boundary \
                       receipt",
            no_claim: "no claim about physical branch selection at grazing; outcomes \
                       are set-valued by design; version-1 prose mints no theorem \
                       authority",
        },
    ]
}

const POLICY_SPEC: &str = "I09_CAMPAIGN_POLICY_V1
ORBIT_PROBLEMS= typed periodic/relative-periodic/hybrid families; admission is typing, \
never existence; descriptor structure declared with index certificates
METHODS= shooting, collocation, harmonic balance are mutually independent routes; \
agreement is checked only between converged runs; one method never adjudicates \
another's convergence
MONODROMY= consistent discrete tangents of the actual discretization; finite \
differences are the independent falsifier route, never the production path
MULTIPLIER_EVIDENCE= existence intervals from residual-backed eigensolver receipts; \
existence intervals never mint multiplicity; no measured transient is ever relabeled \
as asymptotic stability; defective/nonnormal cases stay Unknown or set-valued
ADJOINTS= derived from the actual discrete residual; FD-verified; no silent FD \
fallback inside the adjoint API; non-smooth events return typed no-claim records
CONTINUATION= every step carries a receipt; budget exhaustion is typed and resumable; \
branch points are classified or explicitly Unknown, never silently skipped
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

fn i09_fixtures() -> Vec<FixturePin> {
    vec![
        FixturePin {
            id: "i09-campaign-policy-v1",
            source: FixtureSource::AuthoredSpec { spec: POLICY_SPEC },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i09-mathieu-hill-family",
            source: FixtureSource::AuthoredSpec {
                spec: "i09 fixture v1: Mathieu/Hill equations x'' + (delta + eps \
                       cos t) x = 0 with (delta, eps) grids straddling the first \
                       two stability tongues; closed-form/independently tabulated \
                       multipliers and tongue boundaries in-spec; periods 2 pi and \
                       pi families; development indices 0..=16383; seeds key \
                       'i09/mathieu/<k>'",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i09-gear-cycle-family",
            source: FixtureSource::AuthoredSpec {
                spec: "i09 fixture v1: periodically forced gear-mesh stiffness \
                       cycle (parametric stiffness square-ish wave smoothed to \
                       declared C2), torque forcing at mesh frequency ratios \
                       {1, 2, 3}; smooth subfamily only (backlash variants live \
                       in the hybrid family); development indices 0..=16383; \
                       seeds key 'i09/gear/<k>'",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i09-rotor-whirl-family",
            source: FixtureSource::AuthoredSpec {
                spec: "i09 fixture v1: Jeffcott-class rotor with anisotropic \
                       supports and unbalance forcing; whirl orbits across spin \
                       speeds bracketing the critical; gyroscopic variant with \
                       declared skew term; independently derived multiplier \
                       magnitudes at pinned speeds; development indices \
                       0..=16383; seeds key 'i09/rotor/<k>'",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i09-wankel-cycle-family",
            source: FixtureSource::AuthoredSpec {
                spec: "i09 fixture v1: Wankel-style periodic machine cycle — \
                       crank-periodic pressure/inertia forcing on the rotor pose \
                       path (prescribed kinematics per the fs-motion wankel_tube \
                       semantics), QoIs are cycle-averaged loads and phase-locked \
                       state functionals; smooth lane only; development indices \
                       0..=16383; seeds key 'i09/wankel/<k>'",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i09-constrained-dae-orbit-core-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "i09 fixture v1 HOLDOUT: index-2/index-3 constrained \
                       mechanisms (crank-slider loop, four-bar) with periodic \
                       forcing; constructed so the finite multiplier set is known \
                       and infinity modes are countable from the declared index; \
                       adjudication compares typed classifications; core held-out \
                       indices 65536..=69631; one I09.G3 consumer \
                       (i09-problem-admission); seeds key 'i09/dae-orbit/<k>'",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "i09-defective-monodromy-core-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "i09 fixture v1 HOLDOUT: hand-proved adversarial monodromy \
                       battery — Jordan-limit families M(eps) = [[1, 1],[eps^2, \
                       1]]-type with defectiveness at eps = 0, nonnormal \
                       transient-growth pairs with spectral radius < 1 but \
                       transient norms >> 1, and unit-circle clusters at distance \
                       1e-12; correct verdicts (Unknown/set-valued) proved \
                       in-spec; core held-out indices 69632..=73727; one I09.G3 \
                       consumer (i09-floquet-classification); seeds key \
                       'i09/defective/<k>'",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "i09-hybrid-impact-max-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "i09 fixture v1 HOLDOUT: impact oscillator orbits with \
                       restitution e in {0.7, 0.9}; closed-form grazing parameter \
                       values in-spec; period-1 and period-2 impacting branches \
                       plus the grazing boundary; maximal held-out indices \
                       131072..=135167; one I09.G3 consumer \
                       (i09-moonshot-topology-hybrid); seeds key \
                       'i09/impact/<k>'",
            },
            partition: Partition::HeldOut,
        },
    ]
}

#[allow(clippy::too_many_lines)]
fn i09_obligations() -> Vec<ObligationRow> {
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
            leaf: "i09-problem-admission",
            claims_covered: &[
                "i09-typed-orbit-problems",
                "i09-descriptor-infinity-separation",
            ],
            unit_cases: UNIT_CASES,
            g0: "generators: orbit-problem families incl. enumerated malformed classes \
                 and descriptor variants; validity predicate: generator certificate \
                 agreement; laws: canonical round-trip bit-stability, refusal totality, \
                 infinity-mode count agreement with declared index; shrinker: \
                 constraint removal preserving the violated rule; replay seeds per \
                 explicits",
            decks: &[
                "i09-campaign-policy-v1",
                "i09-constrained-dae-orbit-core-holdout",
                "i09-mathieu-hill-family",
            ],
            g3_relations: &[
                "state relabeling invariance",
                "unit-rescaling covariance of admitted problems",
            ],
            g4_schedule: "request->drain->finalize injection between decode, pencil \
                          analysis, and classification; checkpoint at each phase \
                          boundary; drained cancellation leaves no partially admitted \
                          problem",
            g5_matrix: "threads {1,2,7} x deterministic mode x same-ISA replay",
            entry_point: "scripts/e2e/leapfrog/i09_admission.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i09-admission slice)",
            obs_events: &[
                "request.received",
                "cancel.requested",
                "drain.completed",
                "finalize.completed",
                "failure_bundle.retained",
                "adjudication.receipt",
                "orbitproblem.admitted",
                "orbitproblem.refused",
            ],
            replay_command: "scripts/e2e/leapfrog/i09_admission.sh --manifest \
                             <manifest-id> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i09-cross-method-orbits",
            claims_covered: &["i09-cross-method-orbit-agreement"],
            unit_cases: UNIT_CASES,
            g0: "generators: smooth fixture families; laws: per-method residual \
                 honesty, convergence-criterion enforcement before comparison, \
                 discrepancy symmetry; replay seeds per explicits",
            decks: &[
                "i09-campaign-policy-v1",
                "i09-gear-cycle-family",
                "i09-mathieu-hill-family",
                "i09-rotor-whirl-family",
                "i09-wankel-cycle-family",
            ],
            g3_relations: &[
                "phase-condition change covariance",
                "time-translation invariance of orbit discrepancy",
            ],
            g4_schedule: "request->drain->finalize injection inside each method's \
                          solve; checkpoint at iteration boundaries; drained \
                          cancellation yields a typed resumable state, never a \
                          comparison against a partial solve",
            g5_matrix: "threads {1,2,7} x deterministic mode x same-ISA replay of \
                        orbits and discrepancies",
            entry_point: "scripts/e2e/leapfrog/i09_orbits.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i09-orbits slice)",
            obs_events: &[
                "request.received",
                "cancel.requested",
                "drain.completed",
                "finalize.completed",
                "failure_bundle.retained",
                "adjudication.receipt",
                "orbit.converged",
                "orbit.discrepancy",
            ],
            replay_command: "scripts/e2e/leapfrog/i09_orbits.sh --manifest \
                             <manifest-id> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i09-monodromy-adjoint",
            claims_covered: &[
                "i09-consistent-tangent-monodromy",
                "i09-exact-discrete-adjoints",
            ],
            unit_cases: UNIT_CASES,
            g0: "generators: smooth fixtures with regular period maps; laws: tangent \
                 record vs FD within rtol, adjoint-tangent duality (inner-product \
                 identity), refusal at declared non-differentiable QoIs; replay seeds \
                 per explicits",
            decks: &[
                "i09-campaign-policy-v1",
                "i09-mathieu-hill-family",
                "i09-rotor-whirl-family",
            ],
            g3_relations: &["parameter-shift covariance of period sensitivities"],
            g4_schedule: "request->drain->finalize injection between primal solve, \
                          tangent assembly, and adjoint sweep; checkpoint between \
                          phases; a drained cancellation yields no partial gradient \
                          record",
            g5_matrix: "threads {1,2,7} x deterministic mode x same-ISA replay of \
                        monodromy actions and gradients",
            entry_point: "scripts/e2e/leapfrog/i09_monodromy.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i09-monodromy slice)",
            obs_events: &[
                "request.received",
                "cancel.requested",
                "drain.completed",
                "finalize.completed",
                "failure_bundle.retained",
                "adjudication.receipt",
                "adjoint.recorded",
                "adjoint.refused",
                "monodromy.assembled",
            ],
            replay_command: "scripts/e2e/leapfrog/i09_monodromy.sh --manifest \
                             <manifest-id> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i09-floquet-classification",
            claims_covered: &[
                "i09-certified-multiplier-intervals",
                "i09-false-stability-falsifier",
            ],
            unit_cases: UNIT_CASES,
            g0: "generators: monodromy actions from smooth and descriptor lanes plus \
                 the adversarial battery; laws: interval containment of analytic \
                 multipliers, flag exactly-when-crossing semantics, Unknown at \
                 defective/nonnormal constructions, zero false certificates; replay \
                 seeds per explicits",
            decks: &[
                "i09-campaign-policy-v1",
                "i09-defective-monodromy-core-holdout",
                "i09-mathieu-hill-family",
                "i09-rotor-whirl-family",
            ],
            g3_relations: &[
                "multiplier conjugation symmetry for real monodromy",
                "similarity-transform invariance of certified magnitudes",
            ],
            g4_schedule: "request->drain->finalize injection with deliberate \
                          eigensolver tick-budget exhaustion (typed resumable \
                          outcomes); checkpoint between ticks; cancel mid-tick with \
                          rollback verified",
            g5_matrix: "threads {1,2,7} x deterministic mode x same-ISA replay of \
                        intervals, flags, and receipts",
            entry_point: "scripts/e2e/leapfrog/i09_floquet.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i09-floquet slice)",
            obs_events: &[
                "request.received",
                "cancel.requested",
                "drain.completed",
                "finalize.completed",
                "failure_bundle.retained",
                "adjudication.receipt",
                "multiplier.certified",
                "stability.unknown",
                "tick.budget",
            ],
            replay_command: "scripts/e2e/leapfrog/i09_floquet.sh --manifest \
                             <manifest-id> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i09-continuation",
            claims_covered: &["i09-phase-robust-continuation"],
            unit_cases: UNIT_CASES,
            g0: "generators: branch segments in pinned parameter boxes; laws: per-step \
                 receipt completeness, independent-sample agreement at pinned \
                 parameters, typed budget exhaustion, branch-point candidate \
                 classification totality; replay seeds per explicits",
            decks: &[
                "i09-campaign-policy-v1",
                "i09-gear-cycle-family",
                "i09-mathieu-hill-family",
                "i09-rotor-whirl-family",
            ],
            g3_relations: &["parameter reparameterization covariance of branches"],
            g4_schedule: "request->drain->finalize injection inside predictor, \
                          corrector, and branch-point processing; checkpoint at step \
                          boundaries; resume must reproduce the uninterrupted branch \
                          bitwise",
            g5_matrix: "threads {1,2,7} x deterministic mode x same-ISA replay of \
                        branch receipts",
            entry_point: "scripts/e2e/leapfrog/i09_continuation.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i09-continuation slice)",
            obs_events: &[
                "request.received",
                "cancel.requested",
                "drain.completed",
                "finalize.completed",
                "failure_bundle.retained",
                "adjudication.receipt",
                "branch.step",
                "branchpoint.candidate",
                "budget.receipt",
            ],
            replay_command: "scripts/e2e/leapfrog/i09_continuation.sh --manifest \
                             <manifest-id> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i09-moonshot-topology-hybrid",
            claims_covered: &[
                "i09-global-branch-topology",
                "i09-hybrid-grazing-continuation",
            ],
            unit_cases: UNIT_CASES,
            g0: "generators: enumerable parameter boxes and the impact family; laws: \
                 topology completeness vs exhaustive enumeration, grazing outcomes \
                 set-valued at closed-form grazing parameters, localized Unknown \
                 windows carry budget receipts; replay seeds per explicits",
            decks: &[
                "i09-campaign-policy-v1",
                "i09-hybrid-impact-max-holdout",
                "i09-mathieu-hill-family",
            ],
            g3_relations: &["box-refinement monotonicity of the covered set"],
            g4_schedule: "the core of this lane IS G4: request->drain->finalize \
                          injection with budget exhaustion, timeout, and cancellation \
                          inside topology search and hybrid crossing processing; \
                          checkpoint at window boundaries; each outcome drained and \
                          typed; BudgetExhausted stays Unknown",
            g5_matrix: "threads {1,2,7} x deterministic mode x same-ISA replay of \
                        coverage sets and crossing receipts",
            entry_point: "scripts/e2e/leapfrog/i09_moonshot.sh",
            tier: CampaignTier::Max,
            dsr_lane: "dsr quality --tool frankensim (i09-moonshot slice)",
            obs_events: &[
                "request.received",
                "cancel.requested",
                "drain.completed",
                "finalize.completed",
                "failure_bundle.retained",
                "adjudication.receipt",
                "budget.receipt",
                "grazing.setvalued",
                "topology.window",
            ],
            replay_command: "scripts/e2e/leapfrog/i09_moonshot.sh --manifest \
                             <manifest-id> --replay <artifact-id>",
        },
    ]
}

fn i09_waivers() -> Vec<Waiver> {
    vec![
        Waiver {
            subject: "i09-floquet-classification",
            reason: "theorem-qualified Floquet clusters under degeneracy need the \
                     fs-spectral eigensolver service's central proof to land (bead \
                     bfid is in_progress: correction pass plus batch verification \
                     pending); until then interval certificates are existence-only \
                     and degenerate clusters stay Unknown",
            owner: "initiative-09 lane owner",
            predicate: "bead frankensim-ext-spectral-eigensolver-service-bfid closes \
                        with green central proof",
            expiry: "first I09 campaign review after bfid closes; re-justify or \
                     retire at every manifest amendment",
            promotion_effect: "the [F] multiplier-interval claim can close on \
                               existence-only semantics; any degenerate-cluster \
                               qualification stays Unknown and unclosable while the \
                               waiver is live",
        },
        Waiver {
            subject: "i09-moonshot-topology-hybrid",
            reason: "hybrid/grazing continuation requires the I12 event machinery's \
                     [S] lanes green and the true-flow ValidatedStep rung (bead \
                     ow2o) for simulated dynamics; neither dependency has landed \
                     central proof",
            owner: "initiative-09 lane owner",
            predicate: "I12 [S] obligations green at core tier AND bead \
                        frankensim-ext-time-validated-step-ow2o closes with green \
                        central proof",
            expiry: "first I09 campaign review after both dependencies land; \
                     re-justify or retire at every manifest amendment",
            promotion_effect: "the [M] claims stay Unknown and cannot close; [S]/[F] \
                               promotion is unaffected because their obligations run \
                               on smooth and descriptor fixtures only",
        },
    ]
}
