//! Focused I03 electrostatic/EQS VerificationManifest conformance.
//!
//! This target proves the I03-specific lattice, stage-separated held-out
//! partitions, leaf/policy mapping, mutation sensitivity, canonical assembly,
//! and transitive amendment invalidation. Generic schema-gate ordering remains
//! in `tests/conformance.rs`.

use fs_blake3::hash_domain;
use fs_vmanifest::{
    Ambition, CampaignTier, ClaimPolarity, FixturePin, FixtureSource, FreezeRefusal, GauntletTier,
    ManifestDraft, Partition, ToleranceSemantics, i03_draft,
};
use std::collections::{BTreeMap, BTreeSet};

const CAMPAIGN_POLICY_FIXTURE: &str = "i03-campaign-policy-v1";

const SOLID_CLAIMS: &[&str] = &[
    "i03-declared-dielectric-laws",
    "i03-electrostatic-convergence-gauss",
    "i03-electrostatic-exact-sequence",
    "i03-eqs-regime-routing",
    "i03-field-circuit-power-closure",
    "i03-floating-conductor-capacitance",
    "i03-force-adjoint-held-variable-closure",
    "i03-steady-conduction-conservation",
];

const FRONTIER_CLAIMS: &[&str] = &[
    "i03-partial-discharge-breakdown-routing",
    "i03-space-charge-aging-singular-routing",
];

const MOONSHOT_CLAIMS: &[&str] = &[
    "i03-cohomology-force-naturality-theorem",
    "i03-electrostriction-energy-interface-closure",
    "i03-electrostriction-interface-theorem",
    "i03-refinement-force-defect-enclosure",
    "i03-topology-event-jump-theorem",
    "i03-topology-force-counterexample-search",
];

type ClaimAuthority = (
    &'static str,
    GauntletTier,
    &'static str,
    &'static str,
    ToleranceSemantics,
    &'static str,
);

const CLAIM_AUTHORITY: &[ClaimAuthority] = &[
    (
        "i03-electrostatic-exact-sequence",
        GauntletTier::G0,
        "exact_incidence_boundary_gauge_and_period_verdict",
        "bit",
        ToleranceSemantics::Exact,
        "i03.oracle.exact_sequence.v1 at fs-vmanifest-oracles/i03/exact_sequence.rs::check_incidence_and_gauge",
    ),
    (
        "i03-electrostatic-convergence-gauss",
        GauntletTier::G1,
        "maximum_preregistered_normalized_gauss_error_and_directed_order_shortfall",
        "1",
        ToleranceSemantics::Interval { lo: 0.0, hi: 1.0 },
        "i03.oracle.electrostatic_mms.v1 at fs-vmanifest-oracles/i03/electrostatic_mms.rs::adjudicate",
    ),
    (
        "i03-floating-conductor-capacitance",
        GauntletTier::G2,
        "maximum_preregistered_normalized_charge_and_matrix_structure_defect",
        "1",
        ToleranceSemantics::Interval { lo: 0.0, hi: 1.0 },
        "i03.oracle.capacitance.v1 at fs-vmanifest-oracles/i03/capacitance.rs::solve_and_check",
    ),
    (
        "i03-steady-conduction-conservation",
        GauntletTier::G1,
        "maximum_preregistered_normalized_local_global_balance_and_passivity_defect",
        "1",
        ToleranceSemantics::Interval { lo: 0.0, hi: 1.0 },
        "i03.oracle.conduction.v1 at fs-vmanifest-oracles/i03/conduction.rs::check_flux_balance",
    ),
    (
        "i03-declared-dielectric-laws",
        GauntletTier::G0,
        "material_regime_admission_and_passivity_verdict",
        "bit",
        ToleranceSemantics::Exact,
        "i03.oracle.dielectric_cards.v1 at fs-vmanifest-oracles/i03/dielectric.rs::dispatch_class_checker",
    ),
    (
        "i03-eqs-regime-routing",
        GauntletTier::G2,
        "eqs_admit_escalate_unknown_and_total_current_acceptance_verdict",
        "bit",
        ToleranceSemantics::Exact,
        "i03.oracle.eqs_regime.v1 at fs-vmanifest-oracles/i03/eqs_regime.rs::adjudicate",
    ),
    (
        "i03-field-circuit-power-closure",
        GauntletTier::G1,
        "maximum_preregistered_normalized_total_current_and_discrete_work_balance_defect",
        "1",
        ToleranceSemantics::Interval { lo: 0.0, hi: 1.0 },
        "i03.oracle.field_circuit_power.v1 at fs-vmanifest-oracles/i03/field_circuit_power.rs::integrate_and_balance",
    ),
    (
        "i03-force-adjoint-held-variable-closure",
        GauntletTier::G3,
        "maximum_preregistered_normalized_generalized_force_covector_virtual_work_and_adjoint_discrepancy",
        "1",
        ToleranceSemantics::Interval { lo: 0.0, hi: 1.0 },
        "i03.oracle.force_adjoint.v1 at fs-vmanifest-oracles/i03/force_adjoint.rs::check_ensemble_derivative",
    ),
    (
        "i03-partial-discharge-breakdown-routing",
        GauntletTier::G2,
        "coverage_confidence_sharpness_and_regime_routing_verdict",
        "bit",
        ToleranceSemantics::Exact,
        "i03.oracle.discharge_coverage.v1 at fs-vmanifest-oracles/i03/discharge_coverage.rs::adjudicate_clusters",
    ),
    (
        "i03-space-charge-aging-singular-routing",
        GauntletTier::G3,
        "maximum_preregistered_normalized_charge_state_energy_and_limit_route_defect",
        "1",
        ToleranceSemantics::Interval { lo: 0.0, hi: 1.0 },
        "i03.oracle.space_charge_aging.v1 at fs-vmanifest-oracles/i03/space_charge_aging.rs::check_ledger_and_limits",
    ),
    (
        "i03-electrostriction-energy-interface-closure",
        GauntletTier::G3,
        "maximum_preregistered_normalized_energy_stress_traction_virtual_work_and_adjoint_defect",
        "1",
        ToleranceSemantics::Interval { lo: 0.0, hi: 1.0 },
        "i03.oracle.electrostriction_numeric.v1 at fs-vmanifest-oracles/i03/electrostriction_numeric.rs::enclose_closure",
    ),
    (
        "i03-electrostriction-interface-theorem",
        GauntletTier::G3,
        "independent_electrostriction_theorem_checker_verdict",
        "bit",
        ToleranceSemantics::Exact,
        "i03.oracle.electrostriction_lean.v1 at proofs/i03/Electrostriction.lean::electrostrictionInterfaceClosure checked by pinned Lean4 kernel receipt",
    ),
    (
        "i03-cohomology-force-naturality-theorem",
        GauntletTier::G3,
        "independent_variational_sheaf_force_pullback_checker_verdict",
        "bit",
        ToleranceSemantics::Exact,
        "i03.oracle.topology_force_lean.v1 at proofs/i03/TopologyForce.lean::representationNaturalityForcePullback checked by pinned Lean4 kernel receipt",
    ),
    (
        "i03-refinement-force-defect-enclosure",
        GauntletTier::G3,
        "normalized_pulled_generalized_force_outside_dependency_preserving_signed_defect_enclosure",
        "1",
        ToleranceSemantics::Interval { lo: 0.0, hi: 1.0 },
        "i03.oracle.refinement_force_defect.v2 composite proofs/i03/TopologyForce.lean::refinementForceDefectIdentity plus fs-vmanifest-oracles/i03/refinement_force.rs::enclose_signed_pullback_defect",
    ),
    (
        "i03-topology-event-jump-theorem",
        GauntletTier::G3,
        "independent_topology_event_jump_checker_verdict",
        "bit",
        ToleranceSemantics::Exact,
        "i03.oracle.topology_event_lean.v1 at proofs/i03/TopologyEventJump.lean::distributionalEventBalance checked by pinned Lean4 kernel receipt",
    ),
    (
        "i03-topology-force-counterexample-search",
        GauntletTier::G3,
        "exact_nonvacuity_coverage_and_zero_verified_in_domain_counterexample_verdict",
        "bit",
        ToleranceSemantics::Exact,
        "i03.oracle.topology_force_falsifier.v1 at fs-vmanifest-oracles/i03/topology_force_falsifier.rs::verify_membership_cover_and_minimize",
    ),
];

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum ClaimField {
    Activation,
    Statement,
    Hypotheses,
    OracleTcbOverlap,
}

type ClaimClauseAuthority = (&'static str, ClaimField, &'static [&'static str]);

const CLAIM_CLAUSE_AUTHORITY: &[ClaimClauseAuthority] = &[
    (
        "i03-electrostatic-exact-sequence",
        ClaimField::Statement,
        &["d1*d0=0", "closed-loop periods"],
    ),
    (
        "i03-electrostatic-convergence-gauss",
        ClaimField::Statement,
        &["convergence floors"],
    ),
    (
        "i03-electrostatic-convergence-gauss",
        ClaimField::OracleTcbOverlap,
        &["independently evaluates"],
    ),
    (
        "i03-floating-conductor-capacitance",
        ClaimField::Statement,
        &["q(V)=q0+C V"],
    ),
    (
        "i03-floating-conductor-capacitance",
        ClaimField::Hypotheses,
        &["state-changing voltage secant", "C=C^T exactly"],
    ),
    (
        "i03-steady-conduction-conservation",
        ClaimField::Statement,
        &["cellwise and globally", "H(div)"],
    ),
    (
        "i03-declared-dielectric-laws",
        ClaimField::Statement,
        &["evolution equation", "positive-real"],
    ),
    (
        "i03-eqs-regime-routing",
        ClaimField::Statement,
        &["Unknown/Unsupported"],
    ),
    (
        "i03-eqs-regime-routing",
        ClaimField::Hypotheses,
        &["branch-certified", "signed net work is never a denominator"],
    ),
    (
        "i03-field-circuit-power-closure",
        ClaimField::Statement,
        &["internal W_fc cancels"],
    ),
    (
        "i03-field-circuit-power-closure",
        ClaimField::Hypotheses,
        &[
            "I_total=-integral_Gamma(J_f+partial_t D) dot n",
            "I_total=dot(Q_e_free) only for a blocking port",
            "ALE/Reynolds",
        ],
    ),
    (
        "i03-force-adjoint-held-variable-closure",
        ClaimField::Statement,
        &[
            "full electrical coenergy",
            "quadratic C(x) formulas are corollaries only",
        ],
    ),
    (
        "i03-force-adjoint-held-variable-closure",
        ClaimField::Hypotheses,
        &["Qbar={Q | 1^T Q=-Q_vol}", "radial formula"],
    ),
    (
        "i03-partial-discharge-breakdown-routing",
        ClaimField::Statement,
        &["familywise-99%"],
    ),
    (
        "i03-partial-discharge-breakdown-routing",
        ClaimField::Hypotheses,
        &[
            "synthetic generator conformance",
            "independently governed audited receipt",
            "fresh random reassignment of every external receipt id",
            "prediction byte-identical",
        ],
    ),
    (
        "i03-space-charge-aging-singular-routing",
        ClaimField::Statement,
        &["carrier transfer is separated", "q_bulk^T S_bulk=0"],
    ),
    (
        "i03-electrostriction-energy-interface-closure",
        ClaimField::Hypotheses,
        &["Ciarlet", "mixed bordered/KKT-Hessian"],
    ),
    (
        "i03-electrostriction-interface-theorem",
        ClaimField::Hypotheses,
        &[
            "byte/semantic equivalence",
            "canonical binding receipt",
            "i03.lean-axioms.v1",
            "propext, Quot.sound, and Classical.choice",
        ],
    ),
    (
        "i03-electrostriction-interface-theorem",
        ClaimField::Activation,
        &["pre-candidate manifest successor", "FORMAL_PROJECTION"],
    ),
    (
        "i03-cohomology-force-naturality-theorem",
        ClaimField::Statement,
        &[
            "variational-sheaf",
            "stationary condensation",
            "f_A=phi^*f_B",
        ],
    ),
    (
        "i03-cohomology-force-naturality-theorem",
        ClaimField::Hypotheses,
        &[
            "complete descending filtered pronilpotent",
            "certified convergence",
        ],
    ),
    (
        "i03-cohomology-force-naturality-theorem",
        ClaimField::Activation,
        &["pre-candidate manifest successor", "FORMAL_PROJECTION"],
    ),
    (
        "i03-refinement-force-defect-enclosure",
        ClaimField::Statement,
        &[
            "ensemble-transfer",
            "dependency-preserving",
            "formally certified cancellation",
        ],
    ),
    (
        "i03-refinement-force-defect-enclosure",
        ClaimField::Hypotheses,
        &["exact configuration transfer T_x"],
    ),
    (
        "i03-refinement-force-defect-enclosure",
        ClaimField::Activation,
        &["pre-candidate manifest successor", "FORMAL_PROJECTION"],
    ),
    (
        "i03-topology-event-jump-theorem",
        ClaimField::Statement,
        &["Radon-measure"],
    ),
    (
        "i03-topology-event-jump-theorem",
        ClaimField::Hypotheses,
        &["undefined product of distributions"],
    ),
    (
        "i03-topology-event-jump-theorem",
        ClaimField::Activation,
        &["pre-candidate manifest successor", "FORMAL_PROJECTION"],
    ),
    (
        "i03-topology-force-counterexample-search",
        ClaimField::Statement,
        &["nonvacuity"],
    ),
    (
        "i03-topology-force-counterexample-search",
        ClaimField::Hypotheses,
        &[
            "rank/unrank enumeration-completeness",
            "full decorated-object canonicalization",
        ],
    ),
    (
        "i03-topology-force-counterexample-search",
        ClaimField::Activation,
        &[
            "pre-candidate FrozenManifest::amend successor",
            "M0_FORMALIZATION",
        ],
    ),
];

const UNIT_CASES: &[&str] = &[
    "boundary",
    "cancellation",
    "empty",
    "error",
    "happy",
    "max",
    "migration",
    "tie-break",
    "unit-dimension",
];

const LEAF_MAP: &[(&str, CampaignTier, &[&str])] = &[
    (
        "i03-feec-cutfem-electrostatics",
        CampaignTier::Core,
        &[
            "i03-electrostatic-convergence-gauss",
            "i03-electrostatic-exact-sequence",
        ],
    ),
    (
        "i03-conductors-capacitance",
        CampaignTier::Core,
        &["i03-floating-conductor-capacitance"],
    ),
    (
        "i03-conduction-dielectric-laws",
        CampaignTier::Core,
        &[
            "i03-declared-dielectric-laws",
            "i03-steady-conduction-conservation",
        ],
    ),
    (
        "i03-field-circuit-force-adjoint",
        CampaignTier::Core,
        &[
            "i03-eqs-regime-routing",
            "i03-field-circuit-power-closure",
            "i03-force-adjoint-held-variable-closure",
        ],
    ),
    (
        "i03-breakdown-routing",
        CampaignTier::Max,
        &["i03-partial-discharge-breakdown-routing"],
    ),
    (
        "i03-space-charge-aging",
        CampaignTier::Max,
        &["i03-space-charge-aging-singular-routing"],
    ),
    (
        "i03-electrostriction-theorem",
        CampaignTier::Max,
        &[
            "i03-electrostriction-energy-interface-closure",
            "i03-electrostriction-interface-theorem",
        ],
    ),
    (
        "i03-topology-force-theorem-falsifier",
        CampaignTier::Max,
        &[
            "i03-cohomology-force-naturality-theorem",
            "i03-refinement-force-defect-enclosure",
            "i03-topology-event-jump-theorem",
            "i03-topology-force-counterexample-search",
        ],
    ),
];

const G0_CLAUSE_AUTHORITY: &[(&str, &[&str])] = &[
    (
        "i03-feec-cutfem-electrostatics",
        &[
            "directed convergence floors",
            "cross-crate/IR/API roundtrip",
        ],
    ),
    (
        "i03-conductors-capacitance",
        &["q0+C*V", "PSD/kernel span{1}/zero row sums only"],
    ),
    (
        "i03-conduction-dielectric-laws",
        &[
            "exp(+i omega t) uses causal s=0+i omega",
            "positive-real Yp only for LTI",
            "nonlinear storage inequality",
        ],
    ),
    (
        "i03-field-circuit-force-adjoint",
        &[
            "moving-port pullback/GCL receipt",
            "integrable/path-independent terminal-charge one-form",
            "complete coenergy/Legendre functional",
        ],
    ),
    (
        "i03-breakdown-routing",
        &[
            "Bonferroni-simultaneous familywise-99%",
            "factor-local Philox lanes",
            "synthetic-versus-experimental authority",
        ],
    ),
    (
        "i03-space-charge-aging",
        &[
            "q_bulk^T S_bulk=0",
            "distinct bulk/surface/electrode/port scores",
            "polarization-derived total-charge equivalence without double counting",
        ],
    ),
    (
        "i03-electrostriction-theorem",
        &[
            "Route H global and Route CN only a.e. authority",
            "Ciarlet-Necas",
            "formal checker bit separate",
        ],
    ),
    (
        "i03-topology-force-theorem-falsifier",
        &[
            "complete-action naturality (not chain/cohomology equivalence alone)",
            "dependency-preserving cancellation",
            "exact T_x pullback with D_ensemble",
            "canonical claim_digest/declaration/proof binding",
            "typed distributional event ownership/weak balance",
        ],
    ),
];

const G4_CLAUSE_AUTHORITY: &[(&str, &[&str])] = &[
    (
        "i03-feec-cutfem-electrostatics",
        &[
            "publishes no partial field",
            "corrupt gauge, topology, cut",
            "minimized FailureBundle",
        ],
    ),
    (
        "i03-conductors-capacitance",
        &[
            "complete authenticated matrix or none",
            "stale gauge/holdout/checkpoint",
            "first minimized structural counterexample",
        ],
    ),
    (
        "i03-conduction-dielectric-laws",
        &[
            "no half-admitted state",
            "missing history",
            "bounded redacted material evidence",
        ],
    ),
    (
        "i03-field-circuit-force-adjoint",
        &[
            "losers drain",
            "duplicate work ownership",
            "exact closure ledger",
        ],
    ),
    (
        "i03-breakdown-routing",
        &[
            "never publishes a partial physical claim",
            "correlated clusters",
            "refutations and minimized integrity failures durably",
        ],
    ),
    (
        "i03-space-charge-aging",
        &[
            "no partial service life",
            "missing owners/history",
            "budget exhaustion",
        ],
    ),
    (
        "i03-electrostriction-theorem",
        &[
            "no theorem color before checker+adjudicator completion",
            "proof rejection",
            "counterexamples are durable",
        ],
    ),
    (
        "i03-topology-force-theorem-falsifier",
        &[
            "loses no discovered candidate",
            "every verified counterexample is durable",
            "cannot be waived away",
        ],
    ),
];

const G5_CLAUSE_AUTHORITY: &[(&str, &str)] = &[
    (
        "i03-feec-cutfem-electrostatics",
        "exact topology, reduction tree, and content identity must match",
    ),
    (
        "i03-conductors-capacitance",
        "bit-identical canonical matrix and terminal receipt identity",
    ),
    (
        "i03-conduction-dielectric-laws",
        "identical admission, state, current, loss, and refusal receipts",
    ),
    (
        "i03-field-circuit-force-adjoint",
        "bit-identical plan, field, port, energy, force, and gradient receipts",
    ),
    (
        "i03-breakdown-routing",
        "exact sample/cluster membership, labels, route, coverage bound, sharpness, and evidence identity",
    ),
    (
        "i03-space-charge-aging",
        "exact state, ledger, route, enclosure, checkpoint, and evidence identity",
    ),
    (
        "i03-electrostriction-theorem",
        "identical theorem bytes, checker receipt, enclosures, and numerical witness",
    ),
    (
        "i03-topology-force-theorem-falsifier",
        "exact candidate set, canonical minima, theorem bytes, checker receipts, and terminal disposition",
    ),
];

const G3_RELATION_AUTHORITY: &[(&str, &[&str])] = &[
    (
        "i03-feec-cutfem-electrostatics",
        &[
            "rigid chart transform leaves normalized observables invariant",
            "under the pinned nested/stable/consistent asymptotic family the certified upper envelope converges at its directed floor, without demanding stepwise raw goal-error monotonicity",
            "a common gauge shift of every conductor and its reference changes potential but not field, charge, energy, or force",
        ],
    ),
    (
        "i03-conductors-capacitance",
        &[
            "terminal permutation conjugates the capacitance matrix",
            "a common shift of all conductors including the enclosure/reference changes no terminal charge",
            "reference-terminal choice preserves reduced physical predictions",
        ],
    ),
    (
        "i03-conduction-dielectric-laws",
        &[
            "series/parallel reduction preserves terminal response only for the declared linear lumped/network subfixtures",
            "for autonomous laws, translating forcing, time origin, and complete material history together preserves response",
            "reciprocal-axis rotation covariantly rotates flux without changing loss",
        ],
    ),
    (
        "i03-field-circuit-force-adjoint",
        &[
            "port orientation reversal flips effort-flow signs but preserves total power",
            "fixed-V, fixed-Q, mixed, and floating ensembles agree only after applying their exact constrained Legendre functionals; no generic opposite-force-sign rule is assumed",
            "an orientation-preserving reparameterization of the same smooth stable equilibrium branch leaves integrated virtual work invariant",
            "deforming a voltage-space integration path preserves scalar coenergy only for an integrable card; a nonzero closed-loop charge integral routes to the complete constrained action",
        ],
    ),
    (
        "i03-breakdown-routing",
        &[
            "unit-rescaled equivalent waveforms preserve normalized onset coverage",
            "permuting the eight within-lot repeats or the IID lot-cluster order leaves each preregistered worst-repeat cell outcome and exact confidence bound unchanged",
            "uniformly widening intervals cannot rescue a sharpness failure even if raw coverage increases",
            "rejection or censoring in one factor-local substream leaves every other factor and cell stream bit-identical",
        ],
    ),
    (
        "i03-space-charge-aging",
        &[
            "refining a charge-source partition preserves the stoichiometric global ledger",
            "a units-and-state rescaling derived from the pinned nondimensionalization preserves the limiting route",
            "approaching a pinned singular family converges in its declared topology or crosses its explicit refusal boundary",
            "Eulerian and ALE moving-interface charge ledgers agree only after the pinned spacetime pullback and geometric-conservation term",
        ],
    ),
    (
        "i03-electrostriction-theorem",
        &[
            "superposed rigid motion preserves energy and covariantly rotates stress",
            "interface orientation reversal swaps side/normal signs without changing balance",
            "admitted contrast-limit sequences preserve the theorem error enclosure",
            "Route H global-homeomorphism authority may be downcast to Route CN a.e.-injectivity authority, never promoted in reverse",
        ],
    ),
    (
        "i03-topology-force-theorem-falsifier",
        &[
            "a complete variational-sheaf transfer carrying either the stationary-condensation or cyclic-L_infinity/BV critical-groupoid receipt pulls the generalized-force one-form back exactly; bare chain/cohomology equivalence is an invalid negative twin",
            "cohomology basis changes with contragredient coordinates and true gauge or compensated exact/potential changes that leave the total physical cochain identical preserve the observable including Gram/Hodge derivatives; an uncompensated exact-cochain shift is a negative twin",
            "general refinement satisfies its signed dependency-preserving enclosure; componentwise-zero defects are sufficient, while a formally certified nontrivial cancellation may also yield exact covariance",
            "a QoI-definition change or square commuting only up to homotopy changes its explicitly owned transfer/homotopy defect; omitting either is an invalid negative twin",
            "a premise-satisfying regularized cobordism preserves typed distributional charge and energy-work balances; a finite jump or impulse needs its additional tightness/integrability premises while class change alone predicts nothing",
            "canonical relabeling preserves theorem/falsifier verdict",
            "mutating any bound claim, card, theorem projection, declaration type, environment, axiom report, or proof identity invalidates theorem-receipt reuse",
        ],
    ),
];

const HOLDOUT_MAP: &[(&str, CampaignTier, &str)] = &[
    (
        "i03-floating-conductors-core-holdout",
        CampaignTier::Core,
        "i03-conductors-capacitance",
    ),
    (
        "i03-eqs-regime-core-holdout",
        CampaignTier::Core,
        "i03-field-circuit-force-adjoint",
    ),
    (
        "i03-force-held-variable-core-holdout",
        CampaignTier::Core,
        "i03-field-circuit-force-adjoint",
    ),
    (
        "i03-insulation-adversaries-max-holdout",
        CampaignTier::Max,
        "i03-breakdown-routing",
    ),
    (
        "i03-space-charge-aging-max-holdout",
        CampaignTier::Max,
        "i03-space-charge-aging",
    ),
    (
        "i03-electrostriction-max-holdout",
        CampaignTier::Max,
        "i03-electrostriction-theorem",
    ),
    (
        "i03-topology-force-max-holdout",
        CampaignTier::Max,
        "i03-topology-force-theorem-falsifier",
    ),
];

const SEEDED_FIXTURE_MAP: &[(&str, Partition, &str, u32, u32)] = &[
    (
        "i03-parallel-plate-mms",
        Partition::Development,
        "i03/parallel-plate",
        0,
        4_095,
    ),
    (
        "i03-coax-sphere-harmonic",
        Partition::Development,
        "i03/harmonic",
        0,
        4_095,
    ),
    (
        "i03-floating-conductors",
        Partition::Development,
        "i03/floating",
        0,
        4_095,
    ),
    (
        "i03-floating-conductors-core-holdout",
        Partition::HeldOut,
        "i03/floating",
        65_536,
        69_631,
    ),
    (
        "i03-resistor-conduction-mms",
        Partition::Development,
        "i03/conduction",
        0,
        4_095,
    ),
    (
        "i03-dispersive-dielectric-cards",
        Partition::Development,
        "i03/dielectric",
        0,
        4_095,
    ),
    (
        "i03-eqs-regime-boundaries",
        Partition::Development,
        "i03/eqs-regime",
        0,
        4_095,
    ),
    (
        "i03-eqs-regime-core-holdout",
        Partition::HeldOut,
        "i03/eqs-regime",
        65_536,
        69_631,
    ),
    (
        "i03-field-circuit-transients",
        Partition::Development,
        "i03/circuit",
        0,
        4_095,
    ),
    (
        "i03-force-held-variable-benchmarks",
        Partition::Development,
        "i03/force",
        0,
        4_095,
    ),
    (
        "i03-force-held-variable-core-holdout",
        Partition::HeldOut,
        "i03/force",
        65_536,
        69_631,
    ),
    (
        "i03-insulation-adversaries",
        Partition::Development,
        "i03/insulation",
        0,
        4_095,
    ),
    (
        "i03-space-charge-aging-limits",
        Partition::Development,
        "i03/space-charge",
        0,
        4_095,
    ),
    (
        "i03-space-charge-aging-max-holdout",
        Partition::HeldOut,
        "i03/space-charge",
        131_072,
        135_167,
    ),
    (
        "i03-electrostriction-finite-strain",
        Partition::Development,
        "i03/electrostriction",
        0,
        4_095,
    ),
    (
        "i03-electrostriction-max-holdout",
        Partition::HeldOut,
        "i03/electrostriction",
        131_072,
        135_167,
    ),
    (
        "i03-topology-force-adversaries",
        Partition::Development,
        "i03/topology-force",
        0,
        4_095,
    ),
    (
        "i03-topology-force-max-holdout",
        Partition::HeldOut,
        "i03/topology-force",
        131_072,
        135_167,
    ),
];

const NO_PUBLIC_SEED_DECLARATION_MAP: &[(&str, Partition)] = &[
    ("i03-insulation-adversaries-max-holdout", Partition::HeldOut),
    ("i03-electrostriction-theorem-card", Partition::Development),
    ("i03-topology-force-theorem-card", Partition::Development),
    (CAMPAIGN_POLICY_FIXTURE, Partition::Development),
];

const FIXTURE_CLAUSE_AUTHORITY: &[(&str, &[&str])] = &[
    (
        "i03-parallel-plate-mms",
        &[
            "manufactured free surface charge",
            "last-three-level energy order >=p-0.20",
            "stored-energy goal-functional order >=2p-0.50",
            "Q_scale=max(1e-12 C",
        ],
    ),
    (
        "i03-coax-sphere-harmonic",
        &[
            "physically zero closed-loop electrostatic period vector",
            "inject a nonzero loop period without EMF",
            "256-bit precision",
            "maximum of its written absolute floor",
        ],
    ),
    (
        "i03-floating-conductors",
        &[
            "q0=q(V=0)",
            "separately typed closed, grounded, infinity, and quotient projections",
            "closed all-conductor vector includes the enclosure",
            "grounding compares C_grounded to the declared principal minor C_uu",
            "exact self-adjoint bilinear-form/terminal-adjointness receipt",
            "Extrapolation is diagnostic",
            "post-hoc structural projection is forbidden",
        ],
    ),
    (
        "i03-floating-conductors-core-holdout",
        &[
            "exact generator/acceptance protocol is i03-floating-conductors",
            "SEALED until the sole consumer",
            "any earlier read is IntegrityFailed",
        ],
    ),
    (
        "i03-resistor-conduction-mms",
        &[
            "dual-cochain reconstruction",
            "compatible/incompatible pure Neumann",
            "missing local-flux receipt is failure",
            "I_scale=max(1e-12 A",
        ],
    ),
    (
        "i03-dispersive-dielectric-cards",
        &[
            "matrix Yp=s*Chi in S/m",
            "stable poles alone never pass",
            "receive no absolute-passivity authority",
        ],
    ),
    (
        "i03-eqs-regime-boundaries",
        &[
            "omega=2*pi*f",
            "ROUTE M is allowed only for separable cases",
            "ROUTE R lifts the EQS field",
            "R_hi>=||R_M(E_EQS)||_(X*)",
            "beta_lo<=beta_R",
            "U_X,R:=R_hi/beta_lo",
            "||E_Maxwell-lift(E_EQS)||_X<=U_X,R",
            "eta_route,R=U_X,R/(X_abs+S_X)",
            "Signed net source work is never a denominator",
            "straddling beta is Unknown/escalated",
        ],
    ),
    (
        "i03-eqs-regime-core-holdout",
        &[
            "exact monitors, threshold multipliers, total-current definition",
            "SEALED until sole consumer",
            "any earlier read is IntegrityFailed",
        ],
    ),
    (
        "i03-field-circuit-transients",
        &[
            "I_k=-integral_Gamma_k(J_f+partial_t D) dot n",
            "Q_e_free,k=-integral_Gamma_k D dot n",
            "Only a blocking port",
            "relative carrier flux and D",
            "discrete GCL exactly",
            "COMBINED LEDGER is their algebraic sum with no Wfc",
            "never residual-fitted",
        ],
    ),
    (
        "i03-force-held-variable-benchmarks",
        &[
            "integral_0^1 V^T q",
            "qhat=Q-q0",
            "Wbar*([V])",
            "Wstar_scale=max(1e-15 J",
            "Only reciprocal linear zero-offset",
            "entire composed geometry/solve/state/QoI path is holomorphic",
            "J/rad for rotations",
        ],
    ),
    (
        "i03-force-held-variable-core-holdout",
        &[
            "exact ensemble formulas, regularity predicates, stencils/enclosures",
            "SEALED until sole consumer",
            "any earlier read is IntegrityFailed",
        ],
    ),
    (
        "i03-insulation-adversaries",
        &[
            "512 deterministic material-lot-shaped clusters",
            "mints no IID or confidence authority",
            "Censoring is therefore repeat/family/endpoint owned",
            "separate 256-bit interval label evaluator",
            "w_lrfe=length(H_lrfe)/(2*S_fe)",
            "alpha=0.001",
            "w_(461) <=0.35",
        ],
    ),
    (
        "i03-insulation-adversaries-max-holdout",
        &[
            "eight-repeat within-lot construction, conservative censoring",
            "simultaneous familywise confidence >=99%",
            "exactly 32768 bytes parsed as 1024 indexed 256-bit blocks",
            "I03StatHoldoutTranscriptV1",
            "an exact tagged union",
            "{CANDIDATE_FIXED=0,COMMIT=1,COMMIT_SET_FIXED=2,BEACON=3,REVEAL=4,FINAL=5}",
            "ten records and contiguous phase_seq 0..9",
            "exact domain value is org.frankensim.i03.stat-holdout.transcript.v1",
            "org.frankensim.i03.stat-holdout.record.v1",
            "record including signature",
            "every COMMIT links to CANDIDATE_FIXED",
            "COMMIT_SET_FIXED also links to CANDIDATE_FIXED",
            "BEACON payload is u32-framed exact source UTF-8, round u64",
            "pre-FINAL root is BLAKE3 derive-key domain",
            "org.frankensim.i03.stat-holdout.transcript-root.v1",
            "org.frankensim.i03.stat-holdout.sampler-root.v1",
            "SamplerLotV1",
            "CellResultV1",
            "mathematical zero must be canonical +0.0 bits=0",
            "-0.0 is rejected",
            "FINAL is therefore nonrecursive",
            "Custodians sign COMMIT/REVEAL",
            "governance key signs CANDIDATE_FIXED/COMMIT_SET_FIXED/BEACON/FINAL",
            "Ed25519 uses RFC8032 PureEdDSA",
            "with SHA-512, empty context, and no prehash",
            "[L]A=0 and [L]R=0",
            "uncofactored equation [S]B_base=R+[k]A",
            "sealed simultaneous phase",
            "all-three-fixed phase receipt",
            "at least one named custodian",
            "probabilities are canonical integers n_i/2^256",
            "33-byte big-endian values",
            "unique half-open interval",
            "bytes [32*(2*l+r),32*(2*l+r+1))",
            "coordinate-wise bijection preserves joint uniformity",
            "No hash, XOF, PRG, rejection stream, short mixed seed",
            "one identical measurable candidate map",
            "stable case/receipt id",
            "fresh random reassignment of all external case/receipt ids",
            "Philox may expand retained within-lot inputs",
            "Synthetic labels are not experimental validation",
        ],
    ),
    (
        "i03-space-charge-aging-limits",
        &[
            "q_bulk^T S_bulk=0 exactly",
            "Q_e_free=-integral_Gamma_e D dot n",
            "I_total=-integral_Gamma_e(J_f+partial_t D) dot n",
            "rho_trap_occ=sum_(i in group)",
            "trap-capacity excess",
            "without evolving any bound quantity independently",
            "every support keeps its own unit",
            "free surface continuity",
            "uniform estimate, and refusal band",
        ],
    ),
    (
        "i03-space-charge-aging-max-holdout",
        &[
            "exact species/stoichiometry, state constraints, card classes",
            "SEALED until sole consumer",
            "any earlier read is IntegrityFailed",
        ],
    ),
    (
        "i03-electrostriction-finite-strain",
        &[
            "W^{1,p}, p>3",
            "interval boundary-degree/self-separation proof",
            "Ciarlet-Necas theorem",
            "MIXED certifies the bordered/KKT-Hessian isomorphism",
            "surface-energy first variation exactly",
        ],
    ),
    (
        "i03-electrostriction-max-holdout",
        &[
            "exact deformation, ensemble, energy, interface",
            "SEALED until sole consumer",
            "any earlier read is IntegrityFailed",
        ],
    ),
    (
        "i03-electrostriction-theorem-card",
        &[
            "generated Lean proposition",
            "schema-v2 manifest canonical claim_digest",
            "retained claim.v1 component domain",
            "Ciarlet-Necas a.e.-injectivity theorem",
            "uniform selected-route stability",
            "FORMAL PROJECTION GATE",
            "mints no theorem promotion authority",
            "i03.lean-axioms.v1",
            "sorryAx",
            "proof-term digests",
        ],
    ),
    (
        "i03-topology-force-adversaries",
        &[
            "TYPED FINITE SUPERGRAMMAR",
            "entire typed decorated tuple",
            "M0_TARGET_BEGIN",
            "M0_TARGET_END",
            "hence N_micro=65536",
            "This target prose is not a machine grammar AST and mints no exhaustive authority",
            "FrozenManifest::amend successor",
            "derived N_micro=65536 proof",
            "org.frankensim.i03.microgrammar-target.v1",
            "exact UTF-8 bytes from M0_TARGET_BEGIN through M0_TARGET_END",
            "EVENT(kind,k)",
            "CONDUCTIVITY_BRIDGE uses the exact C1 cubic smoothstep",
            "whole-campaign preflight proof",
            "exactly 16,777,216 adversarial evaluations total",
            "all 65,536 M0 cases",
            "execution disposition BudgetExhausted with evidence completeness PartialEvidence",
            "SCALAR FUNCTION GRAMMAR",
            "<=128 nonzero bounded-rational tensor coefficients per arity",
            "never called exhaustive",
            "d_x P nonzero",
            "DeclaredDivergent control",
        ],
    ),
    (
        "i03-topology-force-max-holdout",
        &[
            "exact complex, weight/load/trace/port grammar",
            "unseen cases beyond the exhaustive microgrammar",
            "any earlier read is IntegrityFailed",
        ],
    ),
    (
        "i03-topology-force-theorem-card",
        &[
            "each theorem-authority manifest claim id",
            "counterexample-search instead binds its canonical claim_digest",
            "FORMAL PROJECTION GATE",
            "mints no theorem promotion authority",
            "i03.lean-axioms.v1",
            "sorryAx",
            "elaborated-type digest",
            "complete descending filtered pronilpotent cyclic",
            "f_A=phi^*f_B",
            "Componentwise zero is sufficient, not necessary",
            "distribution/Radon-measure limit",
            "SCALE_BINDING fixes",
            "A_dualwork=max(1e-15 J",
        ],
    ),
    (
        CAMPAIGN_POLICY_FIXTURE,
        &[
            "NO_COLLAPSE=",
            "HELDOUT_COMMIT_REVEAL=",
            "THEOREM_AXIOMS=",
            "FORMAL_PROJECTION=",
            "M0_FORMALIZATION=",
            "SCALE_BINDING=",
            "every dimensionless ratio, empirical width/order statistic, and exact bit",
            "PROMOTION=baseline claims consumed only by I03.G4",
            "LEAF_REQUIREMENT=every I03 obligation row references this exact fixture id",
        ],
    ),
];

const DECK_MAP: &[(&str, &[&str])] = &[
    (
        "i03-feec-cutfem-electrostatics",
        &[
            "i03-parallel-plate-mms",
            "i03-coax-sphere-harmonic",
            CAMPAIGN_POLICY_FIXTURE,
        ],
    ),
    (
        "i03-conductors-capacitance",
        &[
            "i03-floating-conductors",
            "i03-floating-conductors-core-holdout",
            "i03-coax-sphere-harmonic",
            CAMPAIGN_POLICY_FIXTURE,
        ],
    ),
    (
        "i03-conduction-dielectric-laws",
        &[
            "i03-resistor-conduction-mms",
            "i03-dispersive-dielectric-cards",
            CAMPAIGN_POLICY_FIXTURE,
        ],
    ),
    (
        "i03-field-circuit-force-adjoint",
        &[
            "i03-eqs-regime-boundaries",
            "i03-eqs-regime-core-holdout",
            "i03-field-circuit-transients",
            "i03-force-held-variable-benchmarks",
            "i03-force-held-variable-core-holdout",
            CAMPAIGN_POLICY_FIXTURE,
        ],
    ),
    (
        "i03-breakdown-routing",
        &[
            "i03-insulation-adversaries",
            "i03-insulation-adversaries-max-holdout",
            "i03-external-hv-industrial-pack",
            CAMPAIGN_POLICY_FIXTURE,
        ],
    ),
    (
        "i03-space-charge-aging",
        &[
            "i03-space-charge-aging-limits",
            "i03-space-charge-aging-max-holdout",
            CAMPAIGN_POLICY_FIXTURE,
        ],
    ),
    (
        "i03-electrostriction-theorem",
        &[
            "i03-electrostriction-finite-strain",
            "i03-electrostriction-max-holdout",
            "i03-electrostriction-theorem-card",
            CAMPAIGN_POLICY_FIXTURE,
        ],
    ),
    (
        "i03-topology-force-theorem-falsifier",
        &[
            "i03-topology-force-adversaries",
            "i03-topology-force-max-holdout",
            "i03-topology-force-theorem-card",
            "i03-coax-sphere-harmonic",
            CAMPAIGN_POLICY_FIXTURE,
        ],
    ),
];

type OperationalSpec = (
    &'static str,
    &'static str,
    &'static str,
    &'static str,
    &'static [&'static str],
);

const OPERATIONAL_MAP: &[OperationalSpec] = &[
    (
        "i03-feec-cutfem-electrostatics",
        "scripts/e2e/leapfrog/i03_electrostatics.sh",
        "env FRANKENSIM_VMANIFEST_LEAF=i03-feec-cutfem-electrostatics dsr quality --tool frankensim",
        "scripts/e2e/leapfrog/i03_electrostatics.sh --replay <artifact-id>",
        &[
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
    ),
    (
        "i03-conductors-capacitance",
        "scripts/e2e/leapfrog/i03_capacitance.sh",
        "env FRANKENSIM_VMANIFEST_LEAF=i03-conductors-capacitance dsr quality --tool frankensim",
        "scripts/e2e/leapfrog/i03_capacitance.sh --replay <artifact-id>",
        &[
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
    ),
    (
        "i03-conduction-dielectric-laws",
        "scripts/e2e/leapfrog/i03_material_conduction.sh",
        "env FRANKENSIM_VMANIFEST_LEAF=i03-conduction-dielectric-laws dsr quality --tool frankensim",
        "scripts/e2e/leapfrog/i03_material_conduction.sh --replay <artifact-id>",
        &[
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
    ),
    (
        "i03-field-circuit-force-adjoint",
        "scripts/e2e/leapfrog/i03_field_circuit_force.sh",
        "env FRANKENSIM_VMANIFEST_LEAF=i03-field-circuit-force-adjoint dsr quality --tool frankensim",
        "scripts/e2e/leapfrog/i03_field_circuit_force.sh --replay <artifact-id>",
        &[
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
    ),
    (
        "i03-breakdown-routing",
        "scripts/e2e/leapfrog/i03_breakdown.sh",
        "env FRANKENSIM_VMANIFEST_LEAF=i03-breakdown-routing dsr quality --tool frankensim",
        "scripts/e2e/leapfrog/i03_breakdown.sh --replay <artifact-id>",
        &[
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
    ),
    (
        "i03-space-charge-aging",
        "scripts/e2e/leapfrog/i03_space_charge_aging.sh",
        "env FRANKENSIM_VMANIFEST_LEAF=i03-space-charge-aging dsr quality --tool frankensim",
        "scripts/e2e/leapfrog/i03_space_charge_aging.sh --replay <artifact-id>",
        &[
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
    ),
    (
        "i03-electrostriction-theorem",
        "scripts/e2e/leapfrog/i03_electrostriction.sh",
        "env FRANKENSIM_VMANIFEST_LEAF=i03-electrostriction-theorem dsr quality --tool frankensim",
        "scripts/e2e/leapfrog/i03_electrostriction.sh --replay <artifact-id>",
        &[
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
    ),
    (
        "i03-topology-force-theorem-falsifier",
        "scripts/e2e/leapfrog/i03_topology_force_theorem.sh",
        "env FRANKENSIM_VMANIFEST_LEAF=i03-topology-force-theorem-falsifier dsr quality --tool frankensim",
        "scripts/e2e/leapfrog/i03_topology_force_theorem.sh --replay <artifact-id>",
        &[
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
    ),
];

fn authored_spec(fixture: &FixturePin) -> &'static str {
    match fixture.source {
        FixtureSource::AuthoredSpec { spec } => spec,
        FixtureSource::External { .. } => {
            panic!(
                "fixture '{}' must be an authored I03 specification",
                fixture.id
            )
        }
    }
}

fn fixture<'a>(fixtures: &'a [FixturePin], id: &str) -> &'a FixturePin {
    fixtures
        .iter()
        .find(|candidate| candidate.id == id)
        .unwrap_or_else(|| panic!("missing I03 fixture '{id}'"))
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct SeedDeclaration<'a> {
    alias: &'a str,
    start: u32,
    end: u32,
}

fn seed_declaration(spec: &str) -> Option<SeedDeclaration<'_>> {
    const LABEL: &str = "SEEDS:";
    const MARKER: &str = "SEEDS: Philox alias '";
    let mut labels = spec.match_indices(LABEL);
    let (label_index, _) = labels.next()?;
    let prefix = spec.get(..label_index)?.trim_end();
    if labels.next().is_some()
        || (!prefix.is_empty() && !prefix.ends_with('.'))
        || spec.match_indices("Philox alias").count() != 1
    {
        return None;
    }
    let suffix = spec.get(label_index..)?.strip_prefix(MARKER)?;
    let (alias, suffix) = suffix.split_once("', k=")?;
    let alias_tail = alias.strip_prefix("i03/")?;
    if alias_tail.is_empty()
        || !alias_tail
            .bytes()
            .all(|byte| byte.is_ascii_lowercase() || byte.is_ascii_digit() || byte == b'-')
    {
        return None;
    }
    let range_token = suffix.split_ascii_whitespace().next()?;
    let range = range_token.strip_suffix('.')?;
    let trailing = suffix.strip_prefix(range_token)?;
    if trailing.contains("k=") || trailing.contains("..=") {
        return None;
    }
    let (start, end) = range.split_once("..=")?;
    let canonical_u32 = |token: &str| {
        if token.is_empty()
            || !token.bytes().all(|byte| byte.is_ascii_digit())
            || (token.len() > 1 && token.starts_with('0'))
        {
            return None;
        }
        token.parse::<u32>().ok()
    };
    let start = canonical_u32(start)?;
    let end = canonical_u32(end)?;
    if start > end {
        return None;
    }
    Some(SeedDeclaration { alias, start, end })
}

fn stream_key(alias: &str) -> u64 {
    let digest = hash_domain("org.frankensim.i03.fixture-stream.v1", alias.as_bytes());
    u64::from_le_bytes(
        digest.as_bytes()[..8]
            .try_into()
            .expect("an fs-blake3 digest always has eight prefix bytes"),
    )
}

fn all_claim_ids() -> BTreeSet<&'static str> {
    SOLID_CLAIMS
        .iter()
        .chain(FRONTIER_CLAIMS)
        .chain(MOONSHOT_CLAIMS)
        .copied()
        .collect()
}

fn all_authority_ids() -> BTreeSet<&'static str> {
    all_claim_ids()
        .into_iter()
        .chain(LEAF_MAP.iter().map(|(leaf, _, _)| *leaf))
        .collect()
}

fn clone_in_memory_assembly_state(draft: &ManifestDraft) -> ManifestDraft {
    draft.clone()
}

fn policy_spec() -> &'static str {
    let draft = i03_draft();
    authored_spec(fixture(&draft.fixtures, CAMPAIGN_POLICY_FIXTURE))
}

#[test]
#[allow(clippy::too_many_lines)]
fn i03_freezes_with_exact_lattice_leaf_and_policy_mapping() {
    let frozen = i03_draft().freeze().expect("the I03 seed must freeze");
    assert_eq!(frozen.initiative(), "I03");
    assert_eq!(frozen.version(), 1);
    assert_eq!(frozen.claims().len(), 16);
    assert_eq!(frozen.fixtures().len(), 22);
    assert_eq!(frozen.obligations().len(), 8);
    assert_eq!(frozen.waivers().len(), 1);
    let versions = frozen.explicits().versions.to_ascii_lowercase();
    assert!(versions.contains("fs-vmanifest schema v2"));
    let version_fields: Vec<_> = versions.split(';').map(str::trim).collect();
    for forbidden in [
        "manifest version 1",
        "manifest revision 1",
        "manifest instance version 1",
        "manifest v1",
    ] {
        assert!(
            !version_fields.contains(&forbidden),
            "known legacy numeric-revision field '{forbidden}' duplicates ManifestDraft.version"
        );
    }
    let waiver = frozen.waivers()[0];
    assert_eq!(waiver.subject, "i03-external-hv-industrial-pack");
    assert_eq!(waiver.owner, "I03 implementation and V&V registry beads");
    assert!(waiver.reason.contains("license-restricted"));
    assert!(
        waiver
            .predicate
            .contains("fs-vvreg admits the exact licensed bytes")
    );
    assert!(
        waiver
            .expiry
            .contains("before the first maximal I03 campaign")
    );
    assert!(
        waiver
            .promotion_effect
            .contains("cannot reach maximal promotion")
    );

    let expected_fixtures: BTreeSet<_> = SEEDED_FIXTURE_MAP
        .iter()
        .map(|(id, _, _, _, _)| *id)
        .chain(NO_PUBLIC_SEED_DECLARATION_MAP.iter().map(|(id, _)| *id))
        .collect();
    assert_eq!(
        expected_fixtures.len(),
        SEEDED_FIXTURE_MAP.len() + NO_PUBLIC_SEED_DECLARATION_MAP.len(),
        "fixture authority tables contain duplicate or overlapping ids"
    );
    assert_eq!(
        frozen
            .fixtures()
            .iter()
            .map(|pin| pin.id)
            .collect::<BTreeSet<_>>(),
        expected_fixtures,
        "the complete fixture identity map is authority, not just its cardinality"
    );
    let fixture_clause_ids: BTreeSet<_> =
        FIXTURE_CLAUSE_AUTHORITY.iter().map(|(id, _)| *id).collect();
    assert_eq!(fixture_clause_ids, expected_fixtures);
    assert_eq!(
        FIXTURE_CLAUSE_AUTHORITY.len(),
        fixture_clause_ids.len(),
        "fixture clause authority contains duplicate ids"
    );
    for &(id, fragments) in FIXTURE_CLAUSE_AUTHORITY {
        let spec = authored_spec(fixture(frozen.fixtures(), id));
        for &fragment in fragments {
            assert!(
                spec.contains(fragment),
                "fixture '{id}' lost load-bearing clause '{fragment}'"
            );
        }
    }

    let claims_for = |ambition| {
        frozen
            .claims()
            .iter()
            .filter(|claim| claim.ambition == ambition)
            .map(|claim| claim.id)
            .collect::<BTreeSet<_>>()
    };
    assert_eq!(
        claims_for(Ambition::Solid),
        SOLID_CLAIMS.iter().copied().collect()
    );
    assert_eq!(
        claims_for(Ambition::Frontier),
        FRONTIER_CLAIMS.iter().copied().collect()
    );
    assert_eq!(
        claims_for(Ambition::Moonshot),
        MOONSHOT_CLAIMS.iter().copied().collect()
    );
    let refutations: BTreeSet<_> = frozen
        .claims()
        .iter()
        .filter(|claim| claim.polarity == ClaimPolarity::Refutation)
        .map(|claim| claim.id)
        .collect();
    assert_eq!(
        refutations,
        BTreeSet::from(["i03-topology-force-counterexample-search"])
    );

    let claim_authority_ids: BTreeSet<_> = CLAIM_AUTHORITY
        .iter()
        .map(|(id, _, _, _, _, _)| *id)
        .collect();
    assert_eq!(claim_authority_ids, all_claim_ids());
    assert_eq!(CLAIM_AUTHORITY.len(), claim_authority_ids.len());
    let mut oracle_ids = BTreeSet::new();
    for &(id, tier, qoi, unit, tolerance, oracle) in CLAIM_AUTHORITY {
        let claim = frozen
            .claim(id)
            .unwrap_or_else(|| panic!("missing claim authority '{id}'"));
        assert_eq!(claim.evidence_tier, tier, "wrong Gauntlet tier for {id}");
        assert_eq!(claim.qoi, qoi, "wrong QoI for {id}");
        assert_eq!(claim.unit, unit, "wrong QoI unit for {id}");
        assert_eq!(claim.tolerance, tolerance, "wrong tolerance for {id}");
        assert_eq!(claim.oracle.identity, oracle, "wrong oracle route for {id}");
        assert!(
            oracle_ids.insert(claim.oracle.identity),
            "oracle route reused by {id}"
        );
    }
    assert_eq!(oracle_ids.len(), frozen.claims().len());

    let clause_ids: BTreeSet<_> = CLAIM_CLAUSE_AUTHORITY
        .iter()
        .map(|(id, _, _)| *id)
        .collect();
    assert_eq!(clause_ids, all_claim_ids());
    let clause_fields: BTreeSet<_> = CLAIM_CLAUSE_AUTHORITY
        .iter()
        .map(|(id, field, _)| (*id, *field))
        .collect();
    assert_eq!(
        CLAIM_CLAUSE_AUTHORITY.len(),
        clause_fields.len(),
        "duplicate claim/field clause-authority row"
    );
    for &(id, field, fragments) in CLAIM_CLAUSE_AUTHORITY {
        let claim = frozen.claim(id).expect("clause-pinned claim");
        assert!(
            !fragments.is_empty(),
            "empty clause authority for {id:?}/{field:?}"
        );
        assert_eq!(
            fragments.iter().copied().collect::<BTreeSet<_>>().len(),
            fragments.len(),
            "duplicate clause fragment for {id:?}/{field:?}"
        );
        for &fragment in fragments {
            assert!(
                !fragment.trim().is_empty(),
                "blank clause for {id:?}/{field:?}"
            );
            let present = match field {
                ClaimField::Activation => claim.activation.contains(fragment),
                ClaimField::Statement => claim.statement.contains(fragment),
                ClaimField::Hypotheses => claim
                    .hypotheses
                    .iter()
                    .any(|hypothesis| hypothesis.contains(fragment)),
                ClaimField::OracleTcbOverlap => claim.oracle.tcb_overlap.contains(fragment),
            };
            assert!(
                present,
                "claim '{id}' lost load-bearing {field:?} clause '{fragment}'"
            );
        }
    }
    assert_eq!(
        frozen
            .claim("i03-cohomology-force-naturality-theorem")
            .expect("topology-force naturality theorem")
            .tolerance,
        ToleranceSemantics::Exact
    );

    let expected_units: BTreeSet<_> = UNIT_CASES.iter().copied().collect();
    let mut coverage = BTreeMap::<&str, usize>::new();
    let actual_leaves: BTreeSet<_> = frozen.obligations().iter().map(|row| row.leaf()).collect();
    assert_eq!(
        actual_leaves,
        LEAF_MAP.iter().map(|(leaf, _, _)| *leaf).collect()
    );
    let deck_keys: BTreeSet<_> = DECK_MAP.iter().map(|(leaf, _)| *leaf).collect();
    assert_eq!(deck_keys, actual_leaves);
    assert_eq!(deck_keys.len(), DECK_MAP.len(), "duplicate DECK_MAP key");
    let operational_keys: BTreeSet<_> = OPERATIONAL_MAP
        .iter()
        .map(|(leaf, _, _, _, _)| *leaf)
        .collect();
    assert_eq!(operational_keys, actual_leaves);
    assert_eq!(
        operational_keys.len(),
        OPERATIONAL_MAP.len(),
        "duplicate OPERATIONAL_MAP key"
    );
    let g0_clause_keys: BTreeSet<_> = G0_CLAUSE_AUTHORITY.iter().map(|(leaf, _)| *leaf).collect();
    assert_eq!(g0_clause_keys, actual_leaves);
    assert_eq!(
        g0_clause_keys.len(),
        G0_CLAUSE_AUTHORITY.len(),
        "duplicate G0_CLAUSE_AUTHORITY key"
    );
    let g4_clause_keys: BTreeSet<_> = G4_CLAUSE_AUTHORITY.iter().map(|(leaf, _)| *leaf).collect();
    assert_eq!(g4_clause_keys, actual_leaves);
    assert_eq!(
        g4_clause_keys.len(),
        G4_CLAUSE_AUTHORITY.len(),
        "duplicate G4_CLAUSE_AUTHORITY key"
    );
    let g5_clause_keys: BTreeSet<_> = G5_CLAUSE_AUTHORITY.iter().map(|(leaf, _)| *leaf).collect();
    assert_eq!(g5_clause_keys, actual_leaves);
    assert_eq!(
        g5_clause_keys.len(),
        G5_CLAUSE_AUTHORITY.len(),
        "duplicate G5_CLAUSE_AUTHORITY key"
    );
    let g3_authority_keys: BTreeSet<_> = G3_RELATION_AUTHORITY
        .iter()
        .map(|(leaf, _)| *leaf)
        .collect();
    assert_eq!(g3_authority_keys, actual_leaves);
    assert_eq!(
        g3_authority_keys.len(),
        G3_RELATION_AUTHORITY.len(),
        "duplicate G3_RELATION_AUTHORITY key"
    );
    for (leaf, tier, expected_claims) in LEAF_MAP {
        let row = frozen
            .obligations()
            .iter()
            .find(|row| row.leaf() == *leaf)
            .unwrap_or_else(|| panic!("missing I03 obligation '{leaf}'"));
        assert_eq!(row.tier(), *tier, "wrong campaign tier for {leaf}");
        assert_eq!(
            row.claims_covered()
                .iter()
                .copied()
                .collect::<BTreeSet<_>>(),
            expected_claims.iter().copied().collect(),
            "wrong claim map for {leaf}"
        );
        for &claim in row.claims_covered() {
            *coverage.entry(claim).or_default() += 1;
        }
        let unit_set: BTreeSet<_> = row.unit_cases().iter().copied().collect();
        assert_eq!(unit_set, expected_units, "wrong unit cases for {leaf}");
        assert_eq!(
            row.unit_cases().len(),
            unit_set.len(),
            "duplicate unit case for {leaf}"
        );
        let expected_decks = DECK_MAP
            .iter()
            .find(|(candidate, _)| candidate == leaf)
            .map(|(_, decks)| *decks)
            .unwrap_or_else(|| panic!("missing deck authority for {leaf}"));
        assert_eq!(
            row.decks().iter().copied().collect::<BTreeSet<_>>(),
            expected_decks.iter().copied().collect::<BTreeSet<_>>(),
            "wrong exact deck set for {leaf}"
        );
        assert_eq!(
            row.decks().iter().copied().collect::<BTreeSet<_>>().len(),
            row.decks().len(),
            "duplicate deck for {leaf}"
        );

        let (_, entry, dsr, replay, expected_events) = OPERATIONAL_MAP
            .iter()
            .find(|(candidate, _, _, _, _)| candidate == leaf)
            .copied()
            .unwrap_or_else(|| panic!("missing operational authority for {leaf}"));
        assert_eq!(row.entry_point(), entry, "wrong entry point for {leaf}");
        assert_eq!(
            row.dsr_lane(),
            dsr,
            "wrong executable DSR command for {leaf}"
        );
        assert_eq!(
            row.replay_command(),
            replay,
            "wrong replay command for {leaf}"
        );
        assert_eq!(
            row.obs_events().iter().copied().collect::<BTreeSet<_>>(),
            expected_events.iter().copied().collect::<BTreeSet<_>>(),
            "wrong exact event set for {leaf}"
        );
        assert_eq!(
            row.obs_events().len(),
            expected_events.len(),
            "duplicate event for {leaf}"
        );
        for lifecycle_event in [
            "execution.cancelled",
            "checkpoint.saved",
            "checkpoint.resumed",
            "checkpoint.forked",
        ] {
            assert!(
                row.obs_events().contains(&lifecycle_event),
                "{leaf} omits common lifecycle event {lifecycle_event}"
            );
        }
        let expected_g3 = G3_RELATION_AUTHORITY
            .iter()
            .find(|(candidate, _)| candidate == leaf)
            .map(|(_, relations)| *relations)
            .expect("G3 relation authority");
        assert_eq!(
            row.g3_relations().iter().copied().collect::<BTreeSet<_>>(),
            expected_g3.iter().copied().collect::<BTreeSet<_>>(),
            "wrong exact G3 relation set for {leaf}"
        );
        assert_eq!(
            row.g3_relations()
                .iter()
                .copied()
                .collect::<BTreeSet<_>>()
                .len(),
            row.g3_relations().len(),
            "duplicate G3 relation for {leaf}"
        );
        assert_eq!(
            row.g3_relations().len(),
            expected_g3.len(),
            "G3 relation authority for {leaf} contains duplicates"
        );
        assert!(
            row.g3_relations()
                .iter()
                .all(|relation| !relation.trim().is_empty()),
            "blank G3 relation for {leaf}"
        );
        for required in ["request-drain-finalize", "checkpoint", "resume", "fork"] {
            assert!(
                row.g4_schedule().contains(required),
                "{leaf} G4 schedule omits {required}"
            );
        }
        let required_g4_clauses = G4_CLAUSE_AUTHORITY
            .iter()
            .find(|(candidate, _)| candidate == leaf)
            .map(|(_, clauses)| *clauses)
            .expect("G4 clause authority");
        for clause in required_g4_clauses {
            assert!(
                row.g4_schedule().contains(clause),
                "{leaf} G4 schedule lost load-bearing clause '{clause}'"
            );
        }
        for required in [
            "generators:",
            "predicates:",
            "laws:",
            "shrinkers",
            "cross-crate",
            "IR/API roundtrip",
        ] {
            assert!(row.g0().contains(required), "{leaf} G0 omits {required}");
        }
        let required_g0_clauses = G0_CLAUSE_AUTHORITY
            .iter()
            .find(|(candidate, _)| candidate == leaf)
            .map(|(_, clauses)| *clauses)
            .expect("G0 clause authority");
        for clause in required_g0_clauses {
            assert!(
                row.g0().contains(clause),
                "{leaf} G0 lost load-bearing clause '{clause}'"
            );
        }
        for required in [
            "threads",
            "shards",
            "deterministic",
            "ISA families {Apple-aarch64,x86_64}",
            "bitwise comparison is only within an identical ISA fingerprint",
        ] {
            assert!(
                row.g5_matrix().contains(required),
                "{leaf} G5 matrix omits {required}"
            );
        }
        assert!(
            ["exact", "identical", "bit-identical"]
                .iter()
                .any(|required| row.g5_matrix().contains(required)),
            "{leaf} G5 matrix omits an exact-output identity requirement"
        );
        let required_g5_clause = G5_CLAUSE_AUTHORITY
            .iter()
            .find(|(candidate, _)| candidate == leaf)
            .map(|(_, clause)| *clause)
            .expect("G5 clause authority");
        assert!(
            row.g5_matrix().contains(required_g5_clause),
            "{leaf} G5 matrix lost exact-output authority '{required_g5_clause}'"
        );
    }
    assert_eq!(coverage.len(), 16);
    for claim in all_claim_ids() {
        assert_eq!(
            coverage.get(claim),
            Some(&1),
            "claim '{claim}' must be covered exactly once"
        );
    }

    let policy = policy_spec();
    assert_eq!(policy.lines().next(), Some("I03_CAMPAIGN_POLICY_V1"));
    let mut policy_map = BTreeMap::new();
    for line in policy.lines().skip(1) {
        let (key, value) = line
            .trim()
            .split_once('=')
            .unwrap_or_else(|| panic!("policy line has no key: {line}"));
        assert!(
            policy_map.insert(key, value).is_none(),
            "duplicate policy key {key}"
        );
    }
    let expected_policy_keys: BTreeSet<_> = [
        "ACCESSIBILITY_AGENT_PARITY",
        "AMENDMENT",
        "CLAIM_ADJUDICATIONS",
        "EVIDENCE_COMPLETENESS",
        "EVIDENCE_INTEGRITY",
        "EXECUTION_DISPOSITIONS",
        "FORMAL_PROJECTION",
        "HELDOUT_COMMIT_REVEAL",
        "LEAF_REQUIREMENT",
        "LOGGING",
        "M0_FORMALIZATION",
        "NO_COLLAPSE",
        "OBSERVABLE_DISPOSITIONS",
        "OPERATIONAL_SUPPORT",
        "PERFORMANCE",
        "PREDICATE_OUTCOMES",
        "PROMOTION",
        "PROMOTION_EFFECTS",
        "RETENTION",
        "SCALE_BINDING",
        "THEOREM_AXIOMS",
    ]
    .into_iter()
    .collect();
    assert_eq!(
        policy_map.keys().copied().collect::<BTreeSet<_>>(),
        expected_policy_keys
    );
    assert_eq!(
        policy_map.len(),
        expected_policy_keys.len(),
        "policy keys must be unique"
    );
    assert_eq!(
        policy_map["EXECUTION_DISPOSITIONS"],
        "Completed|Cancelled|TimedOut|BudgetExhausted|InfrastructureFailed"
    );
    assert_eq!(
        policy_map["PREDICATE_OUTCOMES"],
        "Satisfied|Violated|Indeterminate"
    );
    assert_eq!(
        policy_map["CLAIM_ADJUDICATIONS"],
        "Supported|Failed|Refuted|Unknown"
    );
    assert_eq!(
        policy_map["EVIDENCE_COMPLETENESS"],
        "CompleteEvidence|PartialEvidence|NoEvidence"
    );
    assert_eq!(
        policy_map["EVIDENCE_INTEGRITY"],
        "IntegrityVerified|IntegrityFailed"
    );
    assert_eq!(
        policy_map["OPERATIONAL_SUPPORT"],
        "SupportedOperation|UnsupportedOperation"
    );
    assert_eq!(
        policy_map["OBSERVABLE_DISPOSITIONS"],
        "FiniteObservable|DeclaredDivergent"
    );
    assert_eq!(
        policy_map["PROMOTION_EFFECTS"],
        "Promotes|BlocksPromotion|NoPromotionAuthority"
    );
    let require = |key: &str, clauses: &[&str]| {
        let value = policy_map
            .get(key)
            .unwrap_or_else(|| panic!("missing policy key {key}"));
        for clause in clauses {
            assert!(
                value.contains(clause),
                "policy {key} lost owning clause: {clause}"
            );
        }
    };
    require(
        "NO_COLLAPSE",
        &[
            "execution disposition",
            "requested predicate outcome",
            "claim adjudication",
            "evidence completeness",
            "evidence integrity",
            "operational support",
            "observable disposition",
            "DeclaredDivergent is an owned observable/payload disposition",
            "claim adjudication remains Unknown",
            "promotion effect",
        ],
    );
    require(
        "LOGGING",
        &[
            "fs-obs schema-versioned bounded JSONL",
            "BLAKE3 derive-key domain org.frankensim.i03.fixture-stream.v1",
            "ordinary draws occupy factor_id=0",
            "factor_id in 1..=65535",
            "encode it in output-block bits 48..63",
            "rejection cannot shift another factor",
            "local block must be <2^48",
            "requires 1<=n<=2^32",
            "exact/u64 arithmetic",
            "n=2^32 has t=2^32",
            "n=0 or wrapped/u32 threshold arithmetic is IntegrityFailed",
            "finite alias table must have no derived-key collision",
            "human/JSONL semantic parity",
        ],
    );
    require(
        "HELDOUT_COMMIT_REVEAL",
        &[
            "candidate/model/toolchain",
            "I03StatHoldoutTranscriptV1",
            "{CANDIDATE_FIXED,COMMIT,COMMIT_SET_FIXED,BEACON,REVEAL,FINAL}",
            "strict RFC8032 SHA-512 Ed25519",
            "prime-subgroup checks and uncofactored verification equation",
            "BEACON binds source, round, exact 32-byte value",
            "FINAL binds the three reveal digests",
            "nonrecursive pre-FINAL root over signed records 0..8",
            "sampler-output root",
            "sealed simultaneous all-three-fixed phase",
            "exactly 1024 indexed 256-bit blocks",
            "lot-major offset [32*(2*l+r),32*(2*l+r+1))",
            "at-least-one-honest custodian",
            "information-theoretically IID uniform",
            "probabilities n_i/2^256",
            "separately bound 256-bit beacon",
            "exact 257-bit cumulative endpoints",
            "unique half-open [c_i,c_(i+1))",
            "combined only by frozen coordinate-wise XOR",
            "beacon is never a sampling input",
            "cyclic shift/direction order challenge",
            "no hash/XOF/PRG/rejection expansion",
            "one identical candidate map sees only its current lot features",
            "fresh external-id reassignment",
            "abort, retry, or resampling is IntegrityFailed",
            "one-shot final submission",
            "locally regenerable public statistical max-holdout has no untouched or IID authority",
        ],
    );
    require(
        "THEOREM_AXIOMS",
        &[
            "i03.lean-axioms.v1",
            "{propext,Quot.sound,Classical.choice}",
            "complete transitive declaration/environment axiom closure",
            "sorryAx",
            "unsafe/native-oracle proof authority",
        ],
    );
    require(
        "FORMAL_PROJECTION",
        &[
            "prose cards mint no theorem promotion authority",
            "before proof/candidate bytes exist",
            "canonical machine proposition AST bytes",
            "formal-hypothesis-to-runtime-evidence schema",
            "deterministic total AST-to-Lean translator",
            "cannot choose or weaken the proposition",
        ],
    );
    require(
        "M0_FORMALIZATION",
        &[
            "16x16x16x16, N_micro=65536",
            "prose mints no exhaustive authority",
            "all validity and 16 stratum predicates",
            "explicit typed values of 16 parameter tuples",
            "rank/unrank/shard algorithms",
            "cannot choose missing grammar semantics",
        ],
    );
    require(
        "SCALE_BINDING",
        &[
            "||Delta||/(a_abs+a_rel*S)",
            "every candidate-independent exception input",
            "candidate outputs may enter only the frozen numerator/statistic/bit evaluation",
            "cannot set or enlarge a reference scale, threshold, support, or aggregation",
            "formula id, unit id, IEEE bits, and source digest",
            "IntegrityFailed",
        ],
    );
    require(
        "RETENTION",
        &[
            "promotion and independent manifest-adjudication receipts",
            "each non-success emits a replayable FailureBundle",
            "integrity failures durable",
        ],
    );
    require(
        "ACCESSIBILITY_AGENT_PARITY",
        &[
            "without pointer-only interaction",
            "stable non-TTY JSONL surface",
        ],
    );
    require(
        "PERFORMANCE",
        &[
            "come from FiveExplicits",
            "never silently widened thresholds",
        ],
    );
    require(
        "PROMOTION",
        &[
            "baseline claims consumed only by I03.G4",
            "maximal claims consumed only by I03.G7",
            "each heldout fixture has one named stage-local consumer",
            "IntegrityFailed",
        ],
    );
    require(
        "AMENDMENT",
        &[
            "ManifestDraft.version is the sole machine-interpreted instance-revision authority",
            "exact affected descendants invalidated",
            "authenticated amendment lineage plus identical component digest",
        ],
    );
    require(
        "LEAF_REQUIREMENT",
        &[
            "request-drain-finalize plus checkpoint/resume/fork",
            "cross-crate/IR/API roundtrip",
            "independent adjudication receipt",
        ],
    );
}

#[test]
#[allow(clippy::too_many_lines)]
fn i03_holdouts_are_stage_separated_and_statistical_holdout_is_not_publicly_regenerable() {
    let frozen = i03_draft().freeze().expect("freeze");
    let actual_holdouts: BTreeSet<_> = frozen
        .fixtures()
        .iter()
        .filter(|fixture| fixture.partition == Partition::HeldOut)
        .map(|fixture| fixture.id)
        .collect();
    assert_eq!(
        actual_holdouts,
        HOLDOUT_MAP.iter().map(|(id, _, _)| *id).collect()
    );
    assert_eq!(
        actual_holdouts.len(),
        HOLDOUT_MAP.len(),
        "HOLDOUT_MAP ids must be unique"
    );

    let mut derived_keys = BTreeMap::<u64, &str>::new();
    for (id, partition, alias, start, end) in SEEDED_FIXTURE_MAP {
        let pin = fixture(frozen.fixtures(), id);
        assert_eq!(pin.partition, *partition, "wrong partition for {id}");
        let declaration = seed_declaration(authored_spec(pin))
            .unwrap_or_else(|| panic!("{id} must have one exact Philox declaration"));
        assert_eq!(
            declaration,
            SeedDeclaration {
                alias: *alias,
                start: *start,
                end: *end,
            },
            "wrong exact seed authority for {id}"
        );
        let key = stream_key(alias);
        if let Some(previous_alias) = derived_keys.insert(key, *alias) {
            assert_eq!(
                previous_alias, *alias,
                "distinct fixture aliases collide on Philox key {key:#018x}"
            );
        }
    }
    assert_eq!(
        derived_keys.len(),
        SEEDED_FIXTURE_MAP
            .iter()
            .map(|(_, _, alias, _, _)| *alias)
            .collect::<BTreeSet<_>>()
            .len(),
        "each distinct alias must own one collision-free derived key"
    );
    for (id, partition) in NO_PUBLIC_SEED_DECLARATION_MAP {
        let pin = fixture(frozen.fixtures(), id);
        assert_eq!(pin.partition, *partition, "wrong partition for {id}");
        let spec = authored_spec(pin);
        assert!(
            !spec.contains("SEEDS:") && !spec.contains("Philox alias"),
            "fixture {id} must omit a public per-fixture seed declaration"
        );
        assert_eq!(
            seed_declaration(spec),
            None,
            "unseeded fixture {id} unexpectedly gained a seed authority"
        );
    }

    for (id, tier, sole_consumer) in HOLDOUT_MAP {
        let consumer_tier = LEAF_MAP
            .iter()
            .find(|(leaf, _, _)| leaf == sole_consumer)
            .map(|(_, tier, _)| *tier)
            .expect("holdout consumer is a declared leaf");
        assert_eq!(
            consumer_tier, *tier,
            "holdout {id} tier disagrees with consumer {sole_consumer}"
        );
        let spec = authored_spec(fixture(frozen.fixtures(), id));
        if *id == "i03-insulation-adversaries-max-holdout" {
            assert_eq!(
                seed_declaration(spec),
                None,
                "the governed statistical holdout must not be publicly regenerable"
            );
            for required in [
                "Candidate/model/toolchain",
                "{i03-stat-custodian-a,i03-stat-custodian-b,i03-stat-custodian-c}",
                "all-three-fixed phase receipt",
                "Ed25519 uses RFC8032 PureEdDSA",
                "with SHA-512, empty context, and no prehash",
                "[L]A=0 and [L]R=0",
                "BEACON payload is u32-framed exact source UTF-8",
                "pre-FINAL transcript root",
                "SamplerLotV1",
                "independent of the other custodians, every adversarial mask and commitment",
                "future unpredictable exactly 256-bit public beacon",
                "exactly 512 IID material-lot atoms",
                "information-theoretically IID uniform",
                "bytes [32*(2*l+r),32*(2*l+r+1))",
                "unique half-open interval",
                "No hash, XOF, PRG, rejection stream",
                "complete candidate-visible input is exactly",
                "fresh random reassignment of all external case/receipt ids",
                "secret campaign nonce",
                "Exact Clopper-Pearson authority is conditional",
            ] {
                assert!(
                    spec.contains(required),
                    "statistical holdout lost commit/reveal authority '{required}'"
                );
            }
        } else {
            let declaration = seed_declaration(spec).expect("deterministic holdout seed");
            let expected_range = match tier {
                CampaignTier::Core => (65_536, 69_631),
                CampaignTier::Max => (131_072, 135_167),
                CampaignTier::Smoke => panic!("I03 has no held-out smoke partition"),
            };
            assert_eq!((declaration.start, declaration.end), expected_range);
            assert!(
                SEEDED_FIXTURE_MAP
                    .iter()
                    .any(
                        |(_, partition, alias, start, end)| *partition == Partition::Development
                            && *alias == declaration.alias
                            && (*start, *end) == (0, 4_095)
                    ),
                "{id} needs a matching development fixture on alias '{}'",
                declaration.alias
            );
        }
        let consumers: BTreeSet<_> = frozen
            .obligations()
            .iter()
            .filter(|row| row.decks().contains(id))
            .map(|row| row.leaf())
            .collect();
        assert_eq!(
            consumers,
            BTreeSet::from([*sole_consumer]),
            "{id} must have exactly one stage-local consumer"
        );
    }

    let seeds = frozen.explicits().seeds;
    for required in [
        "fixture-declared alias 'i03/<stream>'",
        "BLAKE3 derive-key domain org.frankensim.i03.fixture-stream.v1",
        "counter low64 is case index",
        "counter high64 is output-block ordinal",
        "0..=4095",
        "65536..=69631",
        "131072..=135167",
        "distinct aliases must have distinct derived keys",
        "Statistical IID heldout authority",
        "never a public Philox range",
    ] {
        assert!(
            seeds.contains(required),
            "missing seed authority: {required}"
        );
    }
}

#[test]
fn seed_clause_parser_rejects_ambiguous_or_noncanonical_authority() {
    assert_eq!(
        seed_declaration("fixture. SEEDS: Philox alias 'i03/probe', k=0..=4095."),
        Some(SeedDeclaration {
            alias: "i03/probe",
            start: 0,
            end: 4_095,
        })
    );
    for malformed in [
        "fixture. NOT_SEEDS: Philox alias 'i03/probe', k=0..=4095.",
        "fixture. NOT SEEDS: Philox alias 'i03/probe', k=0..=4095.",
        "fixture. XSEEDS: Philox alias 'i03/probe', k=0..=4095.",
        "fixture. SEEDS: Philox alias 'i03/probe', k=0..=4095. SEEDS: alternative.",
        "fixture. SEEDS: Philox alias 'i03/probe', k=00..=4095.",
        "fixture. SEEDS: Philox alias 'i03/probe', k=+0..=4095.",
        "fixture. SEEDS: Philox alias 'i03/probe', k=4095..=0.",
        "fixture. SEEDS: Philox alias 'i03/Probe', k=0..=4095.",
        "fixture. SEEDS: Philox alias 'i03/probe_name', k=0..=4095.",
        "fixture. SEEDS: Philox alias 'other/probe', k=0..=4095.",
        "fixture. SEEDS: Philox alias 'i03/probe', k=0..=4095. trailing Philox alias 'i03/other'.",
    ] {
        assert_eq!(
            seed_declaration(malformed),
            None,
            "ambiguous/noncanonical seed clause parsed: {malformed}"
        );
    }
}

#[test]
fn i03_identity_is_stable_and_top_level_input_order_invariant() {
    let frozen = i03_draft().freeze().expect("freeze");
    assert_eq!(
        frozen.digest(),
        i03_draft().freeze().expect("refreeze").digest()
    );

    let mut reordered = i03_draft();
    reordered.claims.reverse();
    reordered.fixtures.reverse();
    reordered.obligations.reverse();
    reordered.waivers.reverse();
    assert_eq!(
        reordered.freeze().expect("reordered freeze").digest(),
        frozen.digest()
    );
}

#[test]
fn i03_chunked_in_memory_assembly_matches_one_shot_freeze() {
    // I03 G4 precursor: clone the partial in-memory draft at deterministic
    // chunk boundaries. Frozen identity must depend only on completed content.
    // This does not serialize a checkpoint, restart a process, detect
    // corruption, or exercise runtime cancellation/draining.
    let one_shot = i03_draft();
    let expected = one_shot.clone().freeze().expect("one-shot freeze");
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

    for chunk in one_shot.claims.chunks(3) {
        staged.claims.extend_from_slice(chunk);
        staged = clone_in_memory_assembly_state(&staged);
    }
    for chunk in one_shot.fixtures.chunks(4) {
        staged.fixtures.extend_from_slice(chunk);
        staged = clone_in_memory_assembly_state(&staged);
    }
    for chunk in one_shot.obligations.chunks(2) {
        staged.obligations.extend_from_slice(chunk);
        staged = clone_in_memory_assembly_state(&staged);
    }
    for chunk in one_shot.waivers.chunks(1) {
        staged.waivers.extend_from_slice(chunk);
        staged = clone_in_memory_assembly_state(&staged);
    }

    let chunked = staged.freeze().expect("chunked in-memory freeze");
    assert_eq!(chunked.digest(), expected.digest());
    assert_eq!(chunked, expected);
}

#[test]
fn i03_g3_mutations_refuse_or_move_the_frozen_authority() {
    let baseline = i03_draft().freeze().expect("freeze").digest();

    let mut weakened = i03_draft();
    weakened
        .claims
        .iter_mut()
        .find(|claim| claim.id == "i03-force-adjoint-held-variable-closure")
        .expect("force claim")
        .hypotheses = &["shape path remains smooth"];
    assert_ne!(
        weakened.freeze().expect("weakened successor").digest(),
        baseline
    );

    let mut missing_hypotheses = i03_draft();
    missing_hypotheses
        .claims
        .iter_mut()
        .find(|claim| claim.id == "i03-force-adjoint-held-variable-closure")
        .expect("force claim")
        .hypotheses = &[];
    assert!(matches!(
        missing_hypotheses.freeze(),
        Err(FreezeRefusal::BlankField {
            field: "claim.hypotheses",
            ..
        })
    ));

    let mut correlated_oracle = i03_draft();
    correlated_oracle
        .claims
        .iter_mut()
        .find(|claim| claim.id == "i03-electrostatic-exact-sequence")
        .expect("electrostatic claim")
        .oracle
        .independent = false;
    assert!(matches!(
        correlated_oracle.freeze(),
        Err(FreezeRefusal::ProductionOracleReuse { .. })
    ));

    let mut relaxed = i03_draft();
    relaxed
        .claims
        .iter_mut()
        .find(|claim| claim.id == "i03-force-adjoint-held-variable-closure")
        .expect("force claim")
        .tolerance = ToleranceSemantics::Interval { lo: 0.0, hi: 2.0 };
    assert_ne!(
        relaxed.freeze().expect("relaxed successor").digest(),
        baseline
    );

    let mut swapped_holdout = i03_draft();
    swapped_holdout
        .fixtures
        .iter_mut()
        .find(|fixture| fixture.id == "i03-insulation-adversaries-max-holdout")
        .expect("maximal insulation holdout")
        .source = FixtureSource::AuthoredSpec {
        spec: "unauthorized post-result held-out replacement",
    };
    assert_ne!(
        swapped_holdout
            .freeze()
            .expect("held-out replacement successor")
            .digest(),
        baseline
    );

    let mut missing_policy = i03_draft();
    missing_policy
        .fixtures
        .retain(|fixture| fixture.id != CAMPAIGN_POLICY_FIXTURE);
    assert!(matches!(
        missing_policy.freeze(),
        Err(FreezeRefusal::OrphanDeck { deck, .. }) if deck == CAMPAIGN_POLICY_FIXTURE
    ));
}

#[test]
#[allow(clippy::too_many_lines)]
fn i03_amendments_invalidate_exact_transitive_authority() {
    let frozen = i03_draft().freeze().expect("freeze");

    let mut version_only = i03_draft();
    version_only.version = 2;
    let (version_successor, version_record) =
        frozen.amend(version_only).expect("version-only amendment");
    assert!(
        version_record.invalidated.is_empty(),
        "the numeric revision itself does not falsify byte-identical component evidence"
    );
    assert_eq!(
        (version_record.from_version, version_record.to_version),
        (1, 2)
    );
    assert_eq!(version_record.from_digest, frozen.digest());
    assert_eq!(version_record.to_digest, version_successor.digest());
    assert_eq!(version_successor.version(), 2);
    assert_ne!(version_record.from_digest, version_record.to_digest);

    let mut changed_force = i03_draft();
    changed_force.version = 2;
    changed_force
        .claims
        .iter_mut()
        .find(|claim| claim.id == "i03-force-adjoint-held-variable-closure")
        .expect("force claim")
        .tolerance = ToleranceSemantics::Interval { lo: 0.0, hi: 2.0 };
    let (_, force_record) = frozen.amend(changed_force).expect("force amendment");
    assert_eq!(
        force_record.invalidated,
        vec![
            "i03-field-circuit-force-adjoint",
            "i03-force-adjoint-held-variable-closure",
        ]
    );

    let mut changed_holdout = i03_draft();
    changed_holdout.version = 2;
    changed_holdout
        .fixtures
        .iter_mut()
        .find(|fixture| fixture.id == "i03-insulation-adversaries-max-holdout")
        .expect("maximal insulation holdout")
        .source = FixtureSource::AuthoredSpec {
        spec: "replacement maximal insulation holdout with a new semantic identity",
    };
    let (_, holdout_record) = frozen
        .amend(changed_holdout)
        .expect("held-out fixture amendment");
    assert_eq!(
        holdout_record.invalidated,
        vec![
            "i03-breakdown-routing",
            "i03-partial-discharge-breakdown-routing",
        ]
    );

    let mut changed_policy = i03_draft();
    changed_policy.version = 2;
    changed_policy
        .fixtures
        .iter_mut()
        .find(|fixture| fixture.id == CAMPAIGN_POLICY_FIXTURE)
        .expect("policy fixture")
        .source = FixtureSource::AuthoredSpec {
        spec: "I03_CAMPAIGN_POLICY_V2 unauthorized threshold edit",
    };
    let (_, policy_record) = frozen.amend(changed_policy).expect("policy amendment");
    assert_eq!(
        policy_record
            .invalidated
            .iter()
            .map(String::as_str)
            .collect::<BTreeSet<_>>(),
        all_authority_ids(),
        "policy changes invalidate every claim and execution leaf, but not the fixture id"
    );
    assert_eq!(policy_record.invalidated.len(), 24);

    let mut changed_explicits = i03_draft();
    changed_explicits.version = 2;
    changed_explicits.explicits.capabilities =
        "test-only changed I03 capability authority; every descendant must invalidate";
    let (_, explicits_record) = frozen
        .amend(changed_explicits)
        .expect("Five Explicits amendment");
    assert_eq!(
        explicits_record
            .invalidated
            .iter()
            .map(String::as_str)
            .collect::<BTreeSet<_>>(),
        all_authority_ids()
    );
    assert_eq!(explicits_record.invalidated.len(), 24);

    let mut changed_title = i03_draft();
    changed_title.version = 2;
    changed_title.title = "test-only changed campaign semantics in the title";
    let (_, title_record) = frozen.amend(changed_title).expect("title amendment");
    assert_eq!(
        title_record
            .invalidated
            .iter()
            .map(String::as_str)
            .collect::<BTreeSet<_>>(),
        all_authority_ids(),
        "a title change is global campaign authority, not harmless metadata"
    );
}
