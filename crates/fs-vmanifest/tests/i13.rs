//! Focused I13 cohomology-preserving EM topology-optimization manifest
//! conformance.
//!
//! These tests pin the independently promotable core, maximal theorem
//! ratchets, stage-separated holdouts, campaign-policy semantics, exact leaf
//! ownership, no-false-authority boundaries, canonical identity, and
//! transitive amendment invalidation. They test manifest authority only; they
//! do not claim that any numerical solver, physical model, theorem, or
//! counterexample campaign exists or has passed.

use fs_vmanifest::{
    Ambition, CampaignTier, ClaimPolarity, FixturePin, FixtureSource, FreezeRefusal, GauntletTier,
    I13_FRESH_V2_AUTHORITY_DAG_EDGES, I13_FRESH_V2_AUTHORITY_TAGGED_SUM_NODES,
    I13_PROTOCOL_FIXED_TERMINAL_EVENTS_V2, ManifestDraft, Partition, ToleranceSemantics, i13_draft,
};
use std::collections::{BTreeMap, BTreeSet};

const CAMPAIGN_POLICY_FIXTURE: &str = "i13-campaign-policy-v1";

const SOLID_CLAIMS: &[&str] = &[
    "i13-design-field-identity",
    "i13-energy-consistent-material-interpolation",
    "i13-manufacturability-insulation-audit",
    "i13-route-netlist-realization-audit",
    "i13-terminal-relative-class-audit",
    "i13-topology-event-lineage-invalidation",
];

const FRONTIER_CLAIMS: &[&str] = &[
    "i13-fixed-topology-coupled-adjoint",
    "i13-fixed-topology-optimizer-improvement",
    "i13-governed-industrial-machine-validation",
    "i13-power-force-thermal-stress-closure",
    "i13-robust-guarded-machine-qois",
];

const MOONSHOT_CLAIMS: &[&str] = &[
    "i13-certified-topology-change-continuation",
    "i13-global-robust-manufacturable-optimum",
    "i13-integral-winding-class-synthesis",
    "i13-maximal-theorem-counterexample-search",
    "i13-relative-differential-character-invariants",
    "i13-topology-to-route-performance-theorem",
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
        "i13-design-field-identity",
        GauntletTier::G0,
        "design_field_roundtrip_admission_and_identity_verdict",
        "bit",
        ToleranceSemantics::Exact,
        "i13.oracle.design_field.v1 at fs-vmanifest-oracles/i13/design_field.rs::decode_and_adjudicate",
    ),
    (
        "i13-energy-consistent-material-interpolation",
        GauntletTier::G1,
        "maximum_preregistered_normalized_pure_limit_energy_tangent_objectivity_and_dissipation_defect",
        "1",
        ToleranceSemantics::Interval { lo: 0.0, hi: 1.0 },
        "i13.oracle.material_energy.v1 at fs-vmanifest-oracles/i13/material_energy.rs::differentiate_and_enclose",
    ),
    (
        "i13-terminal-relative-class-audit",
        GauntletTier::G2,
        "terminal_relative_incidence_integral_period_and_naturality_verdict",
        "bit",
        ToleranceSemantics::Exact,
        "i13.oracle.relative_topology.v1 at fs-vmanifest-oracles/i13/relative_topology.rs::rebuild_exact_pairing",
    ),
    (
        "i13-topology-event-lineage-invalidation",
        GauntletTier::G3,
        "topology_event_classification_transaction_and_invalidation_verdict",
        "bit",
        ToleranceSemantics::Exact,
        "i13.oracle.topology_event_lineage.v1 at fs-vmanifest-oracles/i13/topology_event.rs::diff_and_check_invalidation",
    ),
    (
        "i13-manufacturability-insulation-audit",
        GauntletTier::G2,
        "maximum_preregistered_normalized_manufacturing_insulation_thermal_and_stress_margin_violation",
        "1",
        ToleranceSemantics::Interval { lo: 0.0, hi: 1.0 },
        "i13.oracle.manufacturing.v1 at fs-vmanifest-oracles/i13/manufacturing.rs::reconstruct_and_measure",
    ),
    (
        "i13-route-netlist-realization-audit",
        GauntletTier::G2,
        "route_connectivity_isolation_integral_class_and_field_comparison_verdict",
        "bit",
        ToleranceSemantics::Exact,
        "i13.oracle.route_realization.v1 at fs-vmanifest-oracles/i13/route_realization.rs::rebuild_graph_class_and_fields",
    ),
    (
        "i13-fixed-topology-coupled-adjoint",
        GauntletTier::G3,
        "maximum_preregistered_normalized_weighted_adjoint_finite_difference_and_representation_gradient_defect",
        "1",
        ToleranceSemantics::Interval { lo: 0.0, hi: 1.0 },
        "i13.oracle.fixed_adjoint.v1 at fs-vmanifest-oracles/i13/fixed_adjoint.rs::perturb_and_crosscheck",
    ),
    (
        "i13-fixed-topology-optimizer-improvement",
        GauntletTier::G2,
        "fixed_topology_independent_improvement_and_guard_adjudication_verdict",
        "bit",
        ToleranceSemantics::Exact,
        "i13.oracle.optimizer_replay.v1 at fs-vmanifest-oracles/i13/optimizer_replay.rs::recompute_candidate",
    ),
    (
        "i13-robust-guarded-machine-qois",
        GauntletTier::G2,
        "robust_primary_improvement_and_all_guarded_qoi_decision_verdict",
        "bit",
        ToleranceSemantics::Exact,
        "i13.oracle.robust_holdout.v1 at fs-vmanifest-oracles/i13/robust_holdout.rs::adjudicate_one_shot",
    ),
    (
        "i13-power-force-thermal-stress-closure",
        GauntletTier::G2,
        "maximum_preregistered_normalized_electrical_field_mechanical_loss_thermal_and_stress_balance_defect",
        "1",
        ToleranceSemantics::Interval { lo: 0.0, hi: 1.0 },
        "i13.oracle.machine_balance.v1 at fs-vmanifest-oracles/i13/machine_balance.rs::integrate_disjoint_ledger",
    ),
    (
        "i13-governed-industrial-machine-validation",
        GauntletTier::G2,
        "governed_industrial_machine_validation_scope_and_simultaneous_interval_verdict",
        "bit",
        ToleranceSemantics::Exact,
        "i13.oracle.industrial_validation.v1 at fs-vmanifest-oracles/i13/industrial_validation.rs::reconstruct_and_adjudicate",
    ),
    (
        "i13-integral-winding-class-synthesis",
        GauntletTier::G2,
        "integral_terminal_balanced_winding_synthesis_and_bound_verdict",
        "bit",
        ToleranceSemantics::Exact,
        "i13.oracle.integer_synthesis.v1 at fs-vmanifest-oracles/i13/integer_synthesis.rs::enumerate_or_verify_bound",
    ),
    (
        "i13-certified-topology-change-continuation",
        GauntletTier::G3,
        "topology_event_class_state_work_restart_and_generalized_derivative_verdict",
        "bit",
        ToleranceSemantics::Exact,
        "i13.oracle.topology_continuation.v1 at fs-vmanifest-oracles/i13/topology_continuation.rs::rebuild_relation_transfer_and_limits",
    ),
    (
        "i13-global-robust-manufacturable-optimum",
        GauntletTier::G3,
        "certified_global_robust_objective_gap_and_all_machine_constraints_verdict",
        "bit",
        ToleranceSemantics::Exact,
        "i13.oracle.global_robust.v1 at fs-vmanifest-oracles/i13/global_robust.rs::check_partition_bounds_and_feasibility",
    ),
    (
        "i13-relative-differential-character-invariants",
        GauntletTier::G3,
        "independent_relative_differential_character_theorem_checker_verdict",
        "bit",
        ToleranceSemantics::Exact,
        "i13.oracle.diffchar_lean.v1 at proofs/i13/RelativeDifferentialCharacter.lean::electromechanicalNaturality checked by pinned Lean4 kernel receipt",
    ),
    (
        "i13-topology-to-route-performance-theorem",
        GauntletTier::G3,
        "independent_topology_to_route_existence_construction_and_performance_theorem_verdict",
        "bit",
        ToleranceSemantics::Exact,
        "i13.oracle.topology_route_lean.v1 at proofs/i13/TopologyRoute.lean::realizationAndPerformance checked by pinned Lean4 kernel receipt",
    ),
    (
        "i13-maximal-theorem-counterexample-search",
        GauntletTier::G3,
        "exact_nonvacuity_coverage_and_zero_genuine_admitted_countermodel_verdict",
        "bit",
        ToleranceSemantics::Exact,
        "i13.oracle.maximal_falsifier.v1 at fs-vmanifest-oracles/i13/maximal_falsifier.rs::verify_cover_admission_and_minimize",
    ),
];

const CLAIM_CLAUSES: &[(&str, &[&str])] = &[
    (
        "i13-design-field-identity",
        &[
            "canonical identity",
            "validity predicate",
            "not a physical microstructure",
        ],
    ),
    (
        "i13-governed-industrial-machine-validation",
        &[
            "AP242 product/configuration/occurrence",
            "jointly committed before any protected byte",
            "one governed pack is not AP242 conformance",
        ],
    ),
    (
        "i13-energy-consistent-material-interpolation",
        &[
            "free-energy",
            "arbitrary coefficient blending",
            "held-source convention",
        ],
    ),
    (
        "i13-terminal-relative-class-audit",
        &[
            "physical pair (K,L)",
            "abstract torsion-rich",
            "primitive",
            "self-field exclusion",
            "never from homology alone",
            "homologous finite-radius routes may have different",
        ],
    ),
    (
        "i13-topology-event-lineage-invalidation",
        &[
            "deterministically invalidated",
            "terminal/netlist",
            "do not prove a canonical class map",
        ],
    ),
    (
        "i13-manufacturability-insulation-audit",
        &[
            "clearance",
            "creepage",
            "current-density",
            "AP242 export/reimport errors",
            "current field-derived guard receipts",
            "NoPromotionAuthority",
            "universal manufacturability",
        ],
    ),
    (
        "i13-route-netlist-realization-audit",
        &[
            "finite-cross-section coil route",
            "semantic netlist",
            "does not synthesize",
        ],
    ),
    (
        "i13-fixed-topology-coupled-adjoint",
        &[
            "exact-discrete adjoint",
            "frozen topology",
            "no classical derivative crosses a topology",
        ],
    ),
    (
        "i13-fixed-topology-optimizer-improvement",
        &[
            "frozen fixed-topology",
            "catalog baseline",
            "UPSTREAM_BINDING receipts",
            "not local or global optimality",
        ],
    ),
    (
        "i13-robust-guarded-machine-qois",
        &[
            "frozen uncertainty law",
            "one-shot held-out",
            "physical population",
        ],
    ),
    (
        "i13-power-force-thermal-stress-closure",
        &[
            "one stable owner",
            "held variables",
            "physical validation remains separate",
        ],
    ),
    (
        "i13-integral-winding-class-synthesis",
        &[
            "terminal-balanced integral",
            "mixed-integer synthesizer",
            "not a connected finite-radius coil route",
        ],
    ),
    (
        "i13-certified-topology-change-continuation",
        &[
            "classical saltation",
            "Bouligand",
            "symmetric finite difference",
        ],
    ),
    (
        "i13-global-robust-manufacturable-optimum",
        &[
            "certified global",
            "frozen mixed discrete/continuous machine-design grammar",
            "globality is only inside the frozen model/grammar/uncertainty domain",
        ],
    ),
    (
        "i13-relative-differential-character-invariants",
        &[
            "real-valued lane; R/Lambda",
            "physical period lattice Lambda",
            "topology alone",
        ],
    ),
    (
        "i13-topology-to-route-performance-theorem",
        &[
            "necessary-and-sufficient",
            "semantic netlists/joints",
            "no universal route",
        ],
    ),
    (
        "i13-maximal-theorem-counterexample-search",
        &[
            "GenuineCountermodel",
            "finite survival never proves",
            "checker implementation defect",
        ],
    ),
];

const FIXTURE_AUTHORITY: &[(&str, Partition, &[&str])] = &[
    (
        "i13-abstract-integral-topology",
        Partition::Development,
        &["ABSTRACT_NOT_R3_PHYSICAL", "lens-space", "no material"],
    ),
    (
        CAMPAIGN_POLICY_FIXTURE,
        Partition::Development,
        &[
            "STATUS_ALGEBRA",
            "ACCEPTANCE_ARITHMETIC",
            "FORMAL_PROJECTION",
            "M0_FORMALIZATION",
            "UPSTREAM_BINDING",
            "16777216",
            "LEAF_REQUIREMENT",
        ],
    ),
    (
        "i13-differential-character-theorem-card",
        Partition::Development,
        &[
            "COEFFICIENT LANES",
            "R/Lambda",
            "NONVACUITY",
            "prose has no proof authority",
        ],
    ),
    (
        "i13-fixed-topology-adjoint",
        Partition::Development,
        &[
            "identical discrete residual",
            "symmetric directional difference",
            "OPTIMIZER CHECKS",
            "clean independent final-candidate replay",
            "topology event",
        ],
    ),
    (
        "i13-fixed-topology-adjoint-core-holdout",
        Partition::HeldOut,
        &[
            "sole consumer i13-adjoint-optimizer-core",
            "one-shot",
            "OPTIMIZER STRATA",
            "independently replayed objective reversals",
            "event/nonsmooth negatives",
        ],
    ),
    (
        "i13-global-robust-max",
        Partition::Development,
        &[
            "entire frozen grammar",
            "coverage tree",
            "local stationary point called global",
        ],
    ),
    (
        "i13-global-robust-max-holdout",
        Partition::HeldOut,
        &[
            "sole consumer i13-global-robust-max-leaf",
            "one-shot",
            "beyond the frozen finite grammar",
        ],
    ),
    (
        "i13-integer-winding-synthesis-max",
        Partition::Development,
        &[
            "terminal-balanced integral phase/branch class vector",
            "typed Infeasible/Unknown",
            "BudgetExhausted is Unknown",
            "not a connected finite-radius geometric coil route",
        ],
    ),
    (
        "i13-integer-winding-synthesis-max-holdout",
        Partition::HeldOut,
        &[
            "sole consumer i13-integer-synthesis-max-leaf",
            "one-shot",
            "heuristic failure is not infeasibility",
            "class success is not geometric routability",
        ],
    ),
    (
        "i13-manufacturing-insulation-core-deck",
        Partition::Development,
        &[
            "manufacture, insulation, and the declared current/thermal/stress feasibility guards only",
            "positive conductor/insulation section",
            "lamination direction/thickness",
            "current-density, thermal-hotspot, and stress margins",
            "AP242 product/configuration/occurrence identity",
            "ACCEPTANCE_ARITHMETIC",
            "candidate-independent scales",
            "no route is synthesized",
        ],
    ),
    (
        "i13-route-realization-core-deck",
        Partition::Development,
        &[
            "finite-cross-section conductor solids",
            "semantic phase, branch, parallel-path",
            "galvanic isolation",
            "never repairs or synthesizes",
        ],
    ),
    (
        "i13-manufacturing-insulation-core-holdout",
        Partition::HeldOut,
        &[
            "sole consumer i13-manufacturability-core",
            "inaccessible committed realization key",
            "current-density/thermal/stress guard boundaries",
            "AP242 occurrence/unit/frame",
            "process-invalid",
            "no route may be repaired or synthesized",
        ],
    ),
    (
        "i13-route-realization-core-holdout",
        Partition::HeldOut,
        &[
            "sole consumer i13-route-realization-core",
            "finite-cross-section coil routes",
            "semantic weld/splice/contact",
            "never synthesizes or repairs",
        ],
    ),
    (
        "i13-material-field-mms",
        Partition::Development,
        &[
            "pure endpoints",
            "hidden stateful hysteresis",
            "quotient basis",
        ],
    ),
    (
        "i13-material-field-core-holdout",
        Partition::HeldOut,
        &[
            "sole consumer i13-material-fields-energy-core",
            "one access",
            "no retuning",
        ],
    ),
    (
        "i13-maximal-adversary-grammar",
        Partition::Development,
        &[
            "65536 raw tuples",
            "rank/unrank/shard bijection",
            "16777216",
            "Search survival never proves",
        ],
    ),
    (
        "i13-maximal-adversary-max-holdout",
        Partition::HeldOut,
        &[
            "sole consumer i13-maximal-falsifier-max-leaf",
            "one-shot",
            "GenuineCountermodel",
        ],
    ),
    (
        "i13-multiphysics-closure-core-deck",
        Partition::Development,
        &[
            "PMSM/generator",
            "held variables",
            "one stable owner",
            "power/force/thermal/stress closure residuals",
        ],
    ),
    (
        "i13-robust-qoi-core-deck",
        Partition::Development,
        &[
            "frozen bounded manufacturing/operating uncertainty law",
            "optional-stopping rule",
            "candidate-independent scale",
            "population validity",
        ],
    ),
    (
        "i13-multiphysics-closure-core-holdout",
        Partition::HeldOut,
        &[
            "sole consumer i13-multiphysics-closure-core",
            "one-shot",
            "physical-population authority",
        ],
    ),
    (
        "i13-robust-qoi-core-holdout",
        Partition::HeldOut,
        &[
            "sole consumer i13-robust-qoi-core",
            "frozen uncertainty law",
            "jointly committed",
            "population-representativeness authority",
        ],
    ),
    (
        "i13-r3-relative-topology",
        Partition::Development,
        &[
            "embedded in declared subsets of R3",
            "torsion-free only under the terminal-subspace premise",
            "geometry contact without semantic joint",
            "circuit-work pairing exact",
        ],
    ),
    (
        "i13-r3-relative-topology-core-holdout",
        Partition::HeldOut,
        &[
            "sole consumer i13-relative-topology-lineage-core",
            "physical strata",
            "contact-without-joint",
        ],
    ),
    (
        "i13-topology-continuation-max",
        Partition::Development,
        &[
            "classical saltation",
            "Clarke/set-valued",
            "symmetric finite difference",
        ],
    ),
    (
        "i13-topology-continuation-max-holdout",
        Partition::HeldOut,
        &[
            "sole consumer i13-topology-continuation-max-leaf",
            "One-shot",
            "refuse classical saltation",
        ],
    ),
    (
        "i13-topology-events-core",
        Partition::Development,
        &[
            "detect, bracket",
            "refuse continuation",
            "geometry contact silently creates a joint",
        ],
    ),
    (
        "i13-topology-events-core-holdout",
        Partition::HeldOut,
        &[
            "sole consumer i13-relative-topology-lineage-core",
            "remesh-only",
            "never requires maximal continuation",
        ],
    ),
    (
        "i13-topology-route-theorem-card",
        Partition::Development,
        &[
            "Necessary, sufficient and performance",
            "NO UNIVERSALITY",
            "FORMAL_PROJECTION",
        ],
    ),
];

const SEEDED_FIXTURES: &[(&str, &str, &str)] = &[
    ("i13-material-field-mms", "i13/material-field", "0..=4095"),
    (
        "i13-material-field-core-holdout",
        "i13/material-field",
        "65536..=69631",
    ),
    ("i13-r3-relative-topology", "i13/r3-relative", "0..=4095"),
    (
        "i13-r3-relative-topology-core-holdout",
        "i13/r3-relative",
        "65536..=69631",
    ),
    (
        "i13-topology-events-core",
        "i13/topology-events",
        "0..=4095",
    ),
    (
        "i13-topology-events-core-holdout",
        "i13/topology-events",
        "65536..=69631",
    ),
    (
        "i13-fixed-topology-adjoint",
        "i13/fixed-adjoint",
        "0..=4095",
    ),
    (
        "i13-fixed-topology-adjoint-core-holdout",
        "i13/fixed-adjoint",
        "65536..=69631",
    ),
    (
        "i13-manufacturing-insulation-core-deck",
        "i13/manufacturing-insulation",
        "0..=4095",
    ),
    (
        "i13-manufacturing-insulation-core-holdout",
        "i13/manufacturing-insulation",
        "65536..=69631",
    ),
    (
        "i13-route-realization-core-deck",
        "i13/route-realization",
        "0..=4095",
    ),
    (
        "i13-route-realization-core-holdout",
        "i13/route-realization",
        "65536..=69631",
    ),
    ("i13-robust-qoi-core-deck", "i13/robust-qoi", "0..=4095"),
    (
        "i13-robust-qoi-core-holdout",
        "i13/robust-qoi",
        "65536..=69631",
    ),
    (
        "i13-multiphysics-closure-core-deck",
        "i13/multiphysics-closure",
        "0..=4095",
    ),
    (
        "i13-multiphysics-closure-core-holdout",
        "i13/multiphysics-closure",
        "65536..=69631",
    ),
    (
        "i13-integer-winding-synthesis-max",
        "i13/integer-synthesis",
        "0..=4095",
    ),
    (
        "i13-integer-winding-synthesis-max-holdout",
        "i13/integer-synthesis",
        "131072..=135167",
    ),
    (
        "i13-topology-continuation-max",
        "i13/topology-continuation",
        "0..=4095",
    ),
    (
        "i13-topology-continuation-max-holdout",
        "i13/topology-continuation",
        "131072..=135167",
    ),
    ("i13-global-robust-max", "i13/global-robust", "0..=4095"),
    (
        "i13-global-robust-max-holdout",
        "i13/global-robust",
        "131072..=135167",
    ),
    (
        "i13-maximal-adversary-grammar",
        "i13/maximal-adversary",
        "0..=4095",
    ),
    (
        "i13-maximal-adversary-max-holdout",
        "i13/maximal-adversary",
        "131072..=135167",
    ),
];

const UNSEEDED_GOVERNANCE_FIXTURES: &[&str] = &[
    CAMPAIGN_POLICY_FIXTURE,
    "i13-abstract-integral-topology",
    "i13-differential-character-theorem-card",
    "i13-topology-route-theorem-card",
];

const LEAF_MAP: &[(&str, CampaignTier, &[&str])] = &[
    (
        "i13-adjoint-optimizer-core",
        CampaignTier::Core,
        &[
            "i13-fixed-topology-coupled-adjoint",
            "i13-fixed-topology-optimizer-improvement",
        ],
    ),
    (
        "i13-differential-character-theorem-max",
        CampaignTier::Max,
        &["i13-relative-differential-character-invariants"],
    ),
    (
        "i13-global-robust-max-leaf",
        CampaignTier::Max,
        &["i13-global-robust-manufacturable-optimum"],
    ),
    (
        "i13-industrial-validation-max",
        CampaignTier::Max,
        &["i13-governed-industrial-machine-validation"],
    ),
    (
        "i13-integer-synthesis-max-leaf",
        CampaignTier::Max,
        &["i13-integral-winding-class-synthesis"],
    ),
    (
        "i13-manufacturability-core",
        CampaignTier::Core,
        &["i13-manufacturability-insulation-audit"],
    ),
    (
        "i13-route-realization-core",
        CampaignTier::Core,
        &["i13-route-netlist-realization-audit"],
    ),
    (
        "i13-robust-qoi-core",
        CampaignTier::Core,
        &["i13-robust-guarded-machine-qois"],
    ),
    (
        "i13-multiphysics-closure-core",
        CampaignTier::Core,
        &["i13-power-force-thermal-stress-closure"],
    ),
    (
        "i13-material-fields-energy-core",
        CampaignTier::Core,
        &[
            "i13-design-field-identity",
            "i13-energy-consistent-material-interpolation",
        ],
    ),
    (
        "i13-maximal-falsifier-max-leaf",
        CampaignTier::Max,
        &["i13-maximal-theorem-counterexample-search"],
    ),
    (
        "i13-relative-topology-lineage-core",
        CampaignTier::Core,
        &[
            "i13-terminal-relative-class-audit",
            "i13-topology-event-lineage-invalidation",
        ],
    ),
    (
        "i13-topology-continuation-max-leaf",
        CampaignTier::Max,
        &["i13-certified-topology-change-continuation"],
    ),
    (
        "i13-topology-route-theorem-max",
        CampaignTier::Max,
        &["i13-topology-to-route-performance-theorem"],
    ),
];

const LEAF_DECK_MAP: &[(&str, &[&str])] = &[
    (
        "i13-adjoint-optimizer-core",
        &[
            "i13-fixed-topology-adjoint",
            "i13-fixed-topology-adjoint-core-holdout",
            "i13-material-field-mms",
            "i13-multiphysics-closure-core-deck",
            CAMPAIGN_POLICY_FIXTURE,
        ],
    ),
    (
        "i13-differential-character-theorem-max",
        &[
            "i13-differential-character-theorem-card",
            "i13-r3-relative-topology",
            "i13-abstract-integral-topology",
            CAMPAIGN_POLICY_FIXTURE,
        ],
    ),
    (
        "i13-global-robust-max-leaf",
        &[
            "i13-global-robust-max",
            "i13-global-robust-max-holdout",
            "i13-manufacturing-insulation-core-deck",
            "i13-route-realization-core-deck",
            "i13-robust-qoi-core-deck",
            "i13-multiphysics-closure-core-deck",
            CAMPAIGN_POLICY_FIXTURE,
        ],
    ),
    (
        "i13-industrial-validation-max",
        &[
            "i13-manufacturing-insulation-core-deck",
            "i13-route-realization-core-deck",
            "i13-robust-qoi-core-deck",
            "i13-multiphysics-closure-core-deck",
            "i13-external-industrial-ap242-pack",
            CAMPAIGN_POLICY_FIXTURE,
        ],
    ),
    (
        "i13-integer-synthesis-max-leaf",
        &[
            "i13-integer-winding-synthesis-max",
            "i13-integer-winding-synthesis-max-holdout",
            "i13-r3-relative-topology",
            CAMPAIGN_POLICY_FIXTURE,
        ],
    ),
    (
        "i13-manufacturability-core",
        &[
            "i13-manufacturing-insulation-core-deck",
            "i13-manufacturing-insulation-core-holdout",
            CAMPAIGN_POLICY_FIXTURE,
        ],
    ),
    (
        "i13-material-fields-energy-core",
        &[
            "i13-material-field-mms",
            "i13-material-field-core-holdout",
            CAMPAIGN_POLICY_FIXTURE,
        ],
    ),
    (
        "i13-maximal-falsifier-max-leaf",
        &[
            "i13-maximal-adversary-grammar",
            "i13-maximal-adversary-max-holdout",
            "i13-integer-winding-synthesis-max",
            "i13-topology-continuation-max",
            "i13-global-robust-max",
            "i13-differential-character-theorem-card",
            "i13-topology-route-theorem-card",
            CAMPAIGN_POLICY_FIXTURE,
        ],
    ),
    (
        "i13-multiphysics-closure-core",
        &[
            "i13-multiphysics-closure-core-deck",
            "i13-multiphysics-closure-core-holdout",
            CAMPAIGN_POLICY_FIXTURE,
        ],
    ),
    (
        "i13-relative-topology-lineage-core",
        &[
            "i13-r3-relative-topology",
            "i13-abstract-integral-topology",
            "i13-r3-relative-topology-core-holdout",
            "i13-topology-events-core",
            "i13-topology-events-core-holdout",
            CAMPAIGN_POLICY_FIXTURE,
        ],
    ),
    (
        "i13-robust-qoi-core",
        &[
            "i13-robust-qoi-core-deck",
            "i13-robust-qoi-core-holdout",
            CAMPAIGN_POLICY_FIXTURE,
        ],
    ),
    (
        "i13-route-realization-core",
        &[
            "i13-route-realization-core-deck",
            "i13-route-realization-core-holdout",
            "i13-r3-relative-topology",
            CAMPAIGN_POLICY_FIXTURE,
        ],
    ),
    (
        "i13-topology-continuation-max-leaf",
        &[
            "i13-topology-continuation-max",
            "i13-topology-continuation-max-holdout",
            "i13-topology-events-core",
            CAMPAIGN_POLICY_FIXTURE,
        ],
    ),
    (
        "i13-topology-route-theorem-max",
        &[
            "i13-topology-route-theorem-card",
            "i13-r3-relative-topology",
            "i13-integer-winding-synthesis-max",
            "i13-manufacturing-insulation-core-deck",
            "i13-route-realization-core-deck",
            CAMPAIGN_POLICY_FIXTURE,
        ],
    ),
];

const CORE_HELDOUT_CONSUMERS: &[(&str, &str)] = &[
    (
        "i13-fixed-topology-adjoint-core-holdout",
        "i13-adjoint-optimizer-core",
    ),
    (
        "i13-manufacturing-insulation-core-holdout",
        "i13-manufacturability-core",
    ),
    (
        "i13-route-realization-core-holdout",
        "i13-route-realization-core",
    ),
    (
        "i13-material-field-core-holdout",
        "i13-material-fields-energy-core",
    ),
    (
        "i13-multiphysics-closure-core-holdout",
        "i13-multiphysics-closure-core",
    ),
    ("i13-robust-qoi-core-holdout", "i13-robust-qoi-core"),
    (
        "i13-r3-relative-topology-core-holdout",
        "i13-relative-topology-lineage-core",
    ),
    (
        "i13-topology-events-core-holdout",
        "i13-relative-topology-lineage-core",
    ),
];

const MAX_HELDOUT_CONSUMERS: &[(&str, &str)] = &[
    (
        "i13-global-robust-max-holdout",
        "i13-global-robust-max-leaf",
    ),
    (
        "i13-integer-winding-synthesis-max-holdout",
        "i13-integer-synthesis-max-leaf",
    ),
    (
        "i13-maximal-adversary-max-holdout",
        "i13-maximal-falsifier-max-leaf",
    ),
    (
        "i13-topology-continuation-max-holdout",
        "i13-topology-continuation-max-leaf",
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

const REQUIRED_LIFECYCLE_EVENTS: &[&str] = &[
    "execution.requested",
    "execution.admitted",
    "execution.started",
    "cancellation.requested",
    "cancellation.observed",
    "execution.cancelled",
    "execution.drain.started",
    "execution.drained",
    "execution.finalized",
    "evidence.atomic_publish_committed",
    "artifact.published",
    "checkpoint.saved",
    "checkpoint.resumed",
    "checkpoint.forked",
    "evidence.adjudication_receipt",
    "evidence.failure_bundle",
];

const REQUIRED_GOVERNANCE_EVENTS: &[&str] =
    &["governance.transition", "governance.protected_access"];

const GOVERNED_EVIDENCE_CONSUMERS: &[&str] = &[
    "i13-adjoint-optimizer-core",
    "i13-global-robust-max-leaf",
    "i13-industrial-validation-max",
    "i13-integer-synthesis-max-leaf",
    "i13-manufacturability-core",
    "i13-material-fields-energy-core",
    "i13-maximal-falsifier-max-leaf",
    "i13-multiphysics-closure-core",
    "i13-relative-topology-lineage-core",
    "i13-robust-qoi-core",
    "i13-route-realization-core",
    "i13-topology-continuation-max-leaf",
];

fn authored_spec(fixture: &FixturePin) -> &'static str {
    match fixture.source {
        FixtureSource::AuthoredSpec { spec } => spec,
        FixtureSource::External { .. } => {
            panic!(
                "fixture '{}' must be an authored I13 specification",
                fixture.id
            )
        }
    }
}

fn fixture<'a>(fixtures: &'a [FixturePin], id: &str) -> &'a FixturePin {
    fixtures
        .iter()
        .find(|candidate| candidate.id == id)
        .unwrap_or_else(|| panic!("missing I13 fixture '{id}'"))
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

fn policy_spec() -> &'static str {
    let draft = i13_draft();
    authored_spec(fixture(&draft.fixtures, CAMPAIGN_POLICY_FIXTURE))
}

#[test]
#[allow(clippy::too_many_lines)]
fn i13_freezes_with_exact_lattice_fixture_and_leaf_authority() {
    let frozen = i13_draft().freeze().expect("the I13 seed must freeze");
    assert_eq!(frozen.initiative(), "I13");
    assert_eq!(frozen.version(), 1);
    assert_eq!(frozen.claims().len(), 17);
    assert_eq!(frozen.fixtures().len(), 28);
    assert_eq!(frozen.obligations().len(), 14);
    assert_eq!(frozen.waivers().len(), 1);
    assert!(
        frozen
            .explicits()
            .versions
            .contains("fs-vmanifest schema v2")
    );
    assert!(frozen.explicits().capabilities.contains("core:"));
    assert!(
        frozen
            .explicits()
            .capabilities
            .contains("maximal feature gates:")
    );
    assert!(frozen.explicits().budgets.contains("BudgetExhausted"));
    for fragment in [
        "Philox 4x32-10",
        "exact case-sensitive UTF-8 fixture aliases",
        "without Unicode normalization",
        "ten-round Random123 word convention",
        "M0=0xD2511F53",
        "M1=0xCD9E8D57",
        "W0=0x9E3779B9",
        "W1=0xBB67AE85",
        "For round j=0..9",
        "(hi0,lo0)=mulhilo(M0,c0)",
        "(hi1,lo1)=mulhilo(M1,c2)",
        "c'=(hi1 xor c1 xor k0,lo1,hi0 xor c3 xor k1,lo0)",
        "After each round j<9",
        "k0=k0+W0 and k1=k1+W1 modulo 2^32",
        "realization_key_32=BLAKE3::derive_key('org.frankensim.i13.development-realization.v1', alias_utf8)",
        "independently sampled 32-byte custodian key",
        "d=BLAKE3::derive_key('org.frankensim.i13.fixture-stream.v2', U64LE(byte_len(alias_utf8))||alias_utf8||realization_key_32)",
        "byte_len counts the exact UTF-8 bytes",
        "U64LE is fixed-width unsigned little-endian",
        "k0=LE32(d[0..4])",
        "k1=LE32(d[4..8])",
        "c0=low32(case_index)",
        "c1=high32(case_index)",
        "c2=low32(output_block_ordinal)",
        "c3=high32(output_block_ordinal)",
        "output_block_ordinal starts at zero",
        "lane is exactly one of {0,1,2,3}",
        "never folded into the counter",
        "LE32(r0)||LE32(r1)||LE32(r2)||LE32(r3)",
        "byte offset n selects block floor(n/16)",
        "lane floor((n mod 16)/4)",
        "byte n mod 4",
        "increments after four output words",
        "KAT_GATE",
        "before any fixture realization or evidence generation",
        "at least one public-development vector and one heldout-format vector",
        "nonsecret fixed test realization key",
        "exact alias bytes, realization key, derived digest d",
        "k0, k1, case_index, output_block_ordinal, counter words, r0..r3",
        "serialized output bytes",
        "reproduced by two independent implementations",
        "nonzero high32 case/block counter words",
        "block 0-to-1 and byte offsets 15-to-16",
        "every output lane and exact LE bytes",
        "alias case and non-normalized UTF-8 twins",
        "cross-alias/key/counter collision refusal",
        "Separate KATs pin the holdout-realization commitment and governed transaction-intent encodings",
        "per-field mutation sensitivity",
        "tests encoding only and grants no custody",
        "actual secret custodian key is never exposed by a KAT",
        "Version-1 prose does not invent those computed words",
        "Native-endian casts are forbidden",
        "Development indices are 0..=4095",
        "core held-out 65536..=69631",
        "maximal held-out 131072..=135167",
        "committed inaccessible custodian key and realized artifact root",
    ] {
        assert!(
            frozen.explicits().seeds.contains(fragment),
            "seed algebra lost '{fragment}'"
        );
    }
    assert!(
        frozen
            .explicits()
            .seeds
            .contains("Physical-population validity never comes from Philox")
    );

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
    assert_eq!(
        frozen
            .claims()
            .iter()
            .filter(|claim| claim.polarity == ClaimPolarity::Refutation)
            .map(|claim| claim.id)
            .collect::<BTreeSet<_>>(),
        BTreeSet::from(["i13-maximal-theorem-counterexample-search"])
    );

    let authority_ids: BTreeSet<_> = CLAIM_AUTHORITY
        .iter()
        .map(|(id, _, _, _, _, _)| *id)
        .collect();
    assert_eq!(authority_ids, all_claim_ids());
    assert_eq!(authority_ids.len(), CLAIM_AUTHORITY.len());
    let mut oracle_ids = BTreeSet::new();
    for &(id, tier, qoi, unit, tolerance, oracle) in CLAIM_AUTHORITY {
        let claim = frozen
            .claim(id)
            .unwrap_or_else(|| panic!("missing claim '{id}'"));
        assert_eq!(claim.evidence_tier, tier, "wrong Gauntlet tier for {id}");
        assert_eq!(claim.qoi, qoi, "wrong QoI for {id}");
        assert_eq!(claim.unit, unit, "wrong unit for {id}");
        assert_eq!(claim.tolerance, tolerance, "wrong tolerance for {id}");
        assert_eq!(claim.oracle.identity, oracle, "wrong oracle for {id}");
        assert!(claim.oracle.independent, "production oracle reuse for {id}");
        assert!(
            oracle_ids.insert(claim.oracle.identity),
            "oracle route reused by {id}"
        );
        assert!(!claim.hypotheses.is_empty(), "claim {id} has no hypotheses");
        for required in [claim.activation, claim.kill, claim.fallback, claim.no_claim] {
            assert!(
                !required.trim().is_empty(),
                "claim {id} lost lifecycle authority"
            );
        }
    }

    assert_eq!(
        CLAIM_CLAUSES
            .iter()
            .map(|(id, _)| *id)
            .collect::<BTreeSet<_>>(),
        all_claim_ids()
    );
    assert_eq!(CLAIM_CLAUSES.len(), all_claim_ids().len());
    for &(id, fragments) in CLAIM_CLAUSES {
        let claim = frozen.claim(id).expect("clause-pinned claim");
        let authority = format!(
            "{}\n{}\n{}\n{}\n{}\n{}\n{}",
            claim.statement,
            claim.hypotheses.join("\n"),
            claim.activation,
            claim.kill,
            claim.fallback,
            claim.no_claim,
            claim.oracle.tcb_overlap
        );
        for &fragment in fragments {
            assert!(
                authority.contains(fragment),
                "claim '{id}' lost '{fragment}'"
            );
        }
    }

    assert_eq!(
        frozen
            .fixtures()
            .iter()
            .map(|pin| pin.id)
            .collect::<BTreeSet<_>>(),
        FIXTURE_AUTHORITY.iter().map(|(id, _, _)| *id).collect()
    );
    assert_eq!(FIXTURE_AUTHORITY.len(), frozen.fixtures().len());
    for &(id, partition, fragments) in FIXTURE_AUTHORITY {
        let pin = fixture(frozen.fixtures(), id);
        assert_eq!(pin.partition, partition, "wrong partition for {id}");
        assert!(pin.digest().is_some(), "fixture {id} has no identity");
        let spec = authored_spec(pin);
        for &fragment in fragments {
            assert!(spec.contains(fragment), "fixture '{id}' lost '{fragment}'");
        }
    }

    let seeded_ids: BTreeSet<_> = SEEDED_FIXTURES.iter().map(|(id, _, _)| *id).collect();
    let unseeded_ids: BTreeSet<_> = UNSEEDED_GOVERNANCE_FIXTURES.iter().copied().collect();
    assert_eq!(seeded_ids.len(), SEEDED_FIXTURES.len());
    assert_eq!(unseeded_ids.len(), UNSEEDED_GOVERNANCE_FIXTURES.len());
    assert!(seeded_ids.is_disjoint(&unseeded_ids));
    assert_eq!(
        seeded_ids
            .union(&unseeded_ids)
            .copied()
            .collect::<BTreeSet<_>>(),
        frozen.fixtures().iter().map(|pin| pin.id).collect()
    );
    for &(id, alias, case_range) in SEEDED_FIXTURES {
        let pin = fixture(frozen.fixtures(), id);
        let spec = authored_spec(pin);
        assert!(
            spec.contains(&format!("alias '{alias}'")),
            "fixture '{id}' lost seed alias"
        );
        assert!(
            spec.contains(&format!("case_index={case_range}")),
            "fixture '{id}' lost case-index range"
        );
        match pin.partition {
            Partition::Development => assert_eq!(case_range, "0..=4095"),
            Partition::HeldOut => {
                assert!(matches!(case_range, "65536..=69631" | "131072..=135167"))
            }
        }
    }
    let mut stream_pairs = BTreeMap::<&str, Vec<(&str, Partition, &str)>>::new();
    for &(id, alias, case_range) in SEEDED_FIXTURES {
        stream_pairs.entry(alias).or_default().push((
            id,
            fixture(frozen.fixtures(), id).partition,
            case_range,
        ));
    }
    assert_eq!(stream_pairs.len(), 12);
    for (alias, members) in stream_pairs {
        assert_eq!(
            members.len(),
            2,
            "stream alias '{alias}' must name exactly one development/holdout pair"
        );
        assert_eq!(
            members
                .iter()
                .filter(|(_, partition, range)| {
                    *partition == Partition::Development && *range == "0..=4095"
                })
                .count(),
            1,
            "stream alias '{alias}' lost its sole development range"
        );
        assert_eq!(
            members
                .iter()
                .filter(|(_, partition, range)| {
                    *partition == Partition::HeldOut
                        && matches!(*range, "65536..=69631" | "131072..=135167")
                })
                .count(),
            1,
            "stream alias '{alias}' lost its sole disjoint held-out range"
        );
    }

    let expected_units: BTreeSet<_> = UNIT_CASES.iter().copied().collect();
    let expected_leaves: BTreeSet<_> = LEAF_MAP.iter().map(|(leaf, _, _)| *leaf).collect();
    assert_eq!(
        LEAF_DECK_MAP
            .iter()
            .map(|(leaf, _)| *leaf)
            .collect::<BTreeSet<_>>(),
        expected_leaves,
        "every I13 leaf needs one exact dependency-deck declaration"
    );
    assert_eq!(LEAF_DECK_MAP.len(), expected_leaves.len());
    assert_eq!(
        REQUIRED_LIFECYCLE_EVENTS
            .iter()
            .copied()
            .collect::<BTreeSet<_>>()
            .len(),
        REQUIRED_LIFECYCLE_EVENTS.len(),
        "the I13 lifecycle event contract itself must not contain aliases"
    );
    assert_eq!(
        REQUIRED_GOVERNANCE_EVENTS
            .iter()
            .copied()
            .collect::<BTreeSet<_>>()
            .len(),
        REQUIRED_GOVERNANCE_EVENTS.len(),
        "the I13 governance event contract itself must not contain aliases"
    );
    let governed_consumers: BTreeSet<_> = GOVERNED_EVIDENCE_CONSUMERS.iter().copied().collect();
    assert_eq!(governed_consumers.len(), 12);
    let governed_from_decks: BTreeSet<_> = CORE_HELDOUT_CONSUMERS
        .iter()
        .chain(MAX_HELDOUT_CONSUMERS)
        .map(|(_, consumer)| *consumer)
        .chain(["i13-industrial-validation-max"])
        .collect();
    assert_eq!(governed_consumers, governed_from_decks);
    assert_eq!(
        frozen
            .obligations()
            .iter()
            .map(|row| row.leaf())
            .collect::<BTreeSet<_>>(),
        expected_leaves
    );
    let mut coverage = BTreeMap::<&str, usize>::new();
    for &(leaf, tier, claims) in LEAF_MAP {
        let row = frozen
            .obligations()
            .iter()
            .find(|row| row.leaf() == leaf)
            .unwrap_or_else(|| panic!("missing leaf '{leaf}'"));
        assert_eq!(row.tier(), tier, "wrong tier for {leaf}");
        assert_eq!(
            row.claims_covered()
                .iter()
                .copied()
                .collect::<BTreeSet<_>>(),
            claims.iter().copied().collect(),
            "wrong claim map for {leaf}"
        );
        for &claim in row.claims_covered() {
            *coverage.entry(claim).or_default() += 1;
        }
        assert_eq!(
            row.unit_cases().iter().copied().collect::<BTreeSet<_>>(),
            expected_units,
            "wrong unit cases for {leaf}"
        );
        assert_eq!(row.unit_cases().len(), expected_units.len());
        let expected_decks = LEAF_DECK_MAP
            .iter()
            .find_map(|(candidate, decks)| (*candidate == leaf).then_some(*decks))
            .unwrap_or_else(|| panic!("missing deck map for '{leaf}'"));
        assert_eq!(
            row.decks().iter().copied().collect::<BTreeSet<_>>(),
            expected_decks.iter().copied().collect(),
            "wrong dependency decks for {leaf}"
        );
        assert_eq!(row.decks().len(), expected_decks.len());
        assert!(row.decks().contains(&CAMPAIGN_POLICY_FIXTURE));
        for required in [
            "generators:",
            "predicates:",
            "laws:",
            "shrinkers",
            "cross-crate",
            "IR/API roundtrip",
        ] {
            assert!(row.g0().contains(required), "{leaf} G0 lost '{required}'");
        }
        for required in [
            "request-drain-finalize",
            "checkpoint save/resume/fork",
            "FailureBundle",
        ] {
            assert!(
                row.g4_schedule().contains(required),
                "{leaf} G4 lost '{required}'"
            );
        }
        for required in [
            "threads",
            "shards",
            "mode {deterministic}",
            "ISA families {Apple-aarch64,x86_64}",
            "bitwise comparison is only within an identical ISA fingerprint",
            "output artifact identity must match",
        ] {
            assert!(
                row.g5_matrix().contains(required),
                "{leaf} G5 lost '{required}'"
            );
        }
        assert_eq!(row.entry_point(), "scripts/e2e/leapfrog/i13-em-topology.sh");
        let expected_lane_selector = format!("FRANKENSIM_VMANIFEST_LEAF={leaf}");
        assert!(row.dsr_lane().contains(expected_lane_selector.as_str()));
        let expected_replay =
            format!("scripts/e2e/leapfrog/i13-em-topology.sh --leaf {leaf} --replay <artifact-id>");
        assert_eq!(row.replay_command(), expected_replay.as_str());
        for event in REQUIRED_LIFECYCLE_EVENTS {
            assert_eq!(
                row.obs_events()
                    .iter()
                    .filter(|observed| *observed == event)
                    .count(),
                1,
                "{leaf} must declare lifecycle event {event} exactly once"
            );
        }
        for event in REQUIRED_GOVERNANCE_EVENTS {
            let expected_count = usize::from(governed_consumers.contains(leaf));
            assert_eq!(
                row.obs_events()
                    .iter()
                    .filter(|observed| *observed == event)
                    .count(),
                expected_count,
                "{leaf} has the wrong governed-evidence declaration count for {event}"
            );
        }
        let governance_namespace = row
            .obs_events()
            .iter()
            .copied()
            .filter(|event| event.starts_with("governance."))
            .collect::<BTreeSet<_>>();
        let expected_governance_namespace = if governed_consumers.contains(leaf) {
            REQUIRED_GOVERNANCE_EVENTS.iter().copied().collect()
        } else {
            BTreeSet::new()
        };
        assert_eq!(
            governance_namespace, expected_governance_namespace,
            "{leaf} must use exactly the closed governance.* vocabulary"
        );
        assert!(!row.g3_relations().is_empty(), "{leaf} has no G3 relation");
    }
    assert_eq!(
        coverage.keys().copied().collect::<BTreeSet<_>>(),
        all_claim_ids()
    );
    assert!(
        coverage.values().all(|count| *count == 1),
        "each claim needs one evidence owner"
    );

    let waiver = frozen.waivers()[0];
    assert_eq!(waiver.subject, "i13-external-industrial-ap242-pack");
    assert!(waiver.reason.contains("license restricted"));
    assert!(
        waiver
            .owner
            .contains("frankensim-leapfrog-2026-program-i94v.2.3.8.1")
    );
    assert!(waiver.owner.contains("fs-vvreg custodian"));
    for required in [
        "I13IndustrialDischargeReceiptV1",
        "exact waiver subject and predicate",
        "one atomic FrozenManifest::amend authority transaction",
        "same-ID I13IndustrialDischargeEnvelopeV1 External root",
        "removes this Waiver row",
        "verifies the AmendmentRecord",
        "advances the authority head",
        "publisher, standard, application-protocol, Part-21 edition",
        "corrigenda, schema/profile/conformance",
        "exact licensed bytes",
        "raw and normalized roots",
        "license/export/access/custody",
        "AP242 semantic occurrence/configuration/netlist mapping",
        "material/process/heat/lot/as-built/test-cell lineage",
        "coordinate/unit/frame/time synchronization",
        "calibration certificates/ranges/validity windows/transfer functions",
        "uncertainty/covariance",
        "environment/control/source/held-variable state",
        "operating domain/interpolation/extrapolation",
        "censoring/missingness/failed-run/dependence/discrepancy policy",
        "candidate/model/toolchain/checker/AcceptanceCard commitment",
        "expected QoIs/intervals/bands/scales/simultaneous/multiplicity/stopping",
        "independent AP242/material/test-cell/oracle/reviewer/adjudicator ownership",
        "redaction/disclosure policy",
        "custodian/reviewer/adjudicator signatures",
        "revocation evidence",
        "complete access/attempt ledger",
        "predecessor manifest",
        "domain-separated transaction intent",
        "verified AmendmentRecord separately binds the final successor",
        "atomic transaction advances the authority head",
        "raw/all-zero digest",
        "fixture-only replacement",
        "receipt-only change",
        "waiver-only removal",
        "split transaction",
        "structurally frozen successor is not discharge",
    ] {
        assert!(
            waiver.predicate.contains(required),
            "industrial waiver predicate lost '{required}'"
        );
    }
    for required in [
        "discharge this waiver only in the atomic RealizationCommitted authority transaction",
        "membership fact, aggregate, derived statistic, commitment opening and side channel",
        "agents, tools, descendants and transitive capabilities",
        "explicitly named least-privilege custodian",
        "independent AP242/material/test-cell/calibration reviewers",
        "frozen GovernanceCommitted access/custody/redaction protocol",
        "durable write-ahead complete access/attempt/output ledger",
        "custodian realization begins only after CandidateFrozen",
        "ingest, canonicalize, map, calibrate, construct and independently verify",
        "may release only protocol-approved commitments and receipts",
        "No protected byte, label, membership fact",
        "candidate-side or adjudication-side principal",
        "grants no candidate, adjudication or promotion authority",
        "cannot be assumed by anyone with candidate influence or later independent-adjudication responsibility",
        "only the i13-industrial-validation-max/governance attempt may be admitted",
        "with no science, reveal, adjudication or promotion authority",
        "before i13-industrial-validation-max/science admission",
        "before one-shot RevealedForAdjudication",
        "adjudication starts only after the atomic RealizationCommitted transaction and complete governance-attempt drain/finalize reconciliation",
        "promotion only after Closed",
        "review again for every standard, application-protocol",
        "Part-21 edition, corrigendum, schema/profile/conformance",
        "license/export/access/custody",
        "corpus item/root",
        "semantic mapping",
        "occurrence/configuration/netlist",
        "material/process/heat/lot/as-built/test-cell link",
        "calibration/time synchronization",
        "environment/control/source/held-variable state",
        "operating domain/interpolation/extrapolation",
        "censoring/missingness/failed-run/dependence/discrepancy",
        "QoI/expected interval/band/scale/simultaneous/multiplicity/stopping decision rule",
        "candidate/model/toolchain/checker/AcceptanceCard",
        "oracle/reviewer/adjudicator/custodian",
        "redaction/disclosure",
        "access capability/attempt ledger",
        "envelope/receipt, transaction intent",
        "authority head, signature or revocation change",
    ] {
        assert!(
            waiver.expiry.contains(required),
            "industrial waiver expiry lost '{required}'"
        );
    }
    assert!(!waiver.predicate.contains("Phase-2"));
    assert!(!waiver.expiry.contains("Phase-2"));
    assert!(waiver.promotion_effect.contains("synthetic/analytic"));
    assert!(
        waiver
            .promotion_effect
            .contains("only i13-governed-industrial-machine-validation")
    );
    let waiver_consumers: Vec<_> = frozen
        .obligations()
        .iter()
        .filter(|row| row.decks().contains(&waiver.subject))
        .map(|row| row.leaf())
        .collect();
    assert_eq!(waiver_consumers, ["i13-industrial-validation-max"]);
}

#[test]
fn i13_holdouts_are_stage_separated_single_consumer_and_tier_local() {
    let frozen = i13_draft().freeze().expect("freeze");
    let assert_governance_events = |row: &fs_vmanifest::FrozenObligationRow| {
        for event in REQUIRED_GOVERNANCE_EVENTS {
            assert_eq!(
                row.obs_events()
                    .iter()
                    .filter(|observed| *observed == event)
                    .count(),
                1,
                "{} must declare governed event {event} exactly once",
                row.leaf()
            );
        }
    };
    let held_out: BTreeSet<_> = frozen
        .fixtures()
        .iter()
        .filter(|pin| pin.partition == Partition::HeldOut)
        .map(|pin| pin.id)
        .collect();
    let expected: BTreeSet<_> = CORE_HELDOUT_CONSUMERS
        .iter()
        .chain(MAX_HELDOUT_CONSUMERS)
        .map(|(deck, _)| *deck)
        .collect();
    assert_eq!(held_out, expected);
    assert_eq!(held_out.len(), 12);

    for &(deck, sole_consumer) in CORE_HELDOUT_CONSUMERS {
        let consumers: Vec<_> = frozen
            .obligations()
            .iter()
            .filter(|row| row.decks().contains(&deck))
            .collect();
        assert_eq!(consumers.len(), 1, "core holdout {deck} needs one consumer");
        assert_eq!(consumers[0].leaf(), sole_consumer);
        assert_eq!(consumers[0].tier(), CampaignTier::Core);
        assert_governance_events(consumers[0]);
        let spec = authored_spec(fixture(frozen.fixtures(), deck));
        assert!(spec.contains("sole consumer"));
        assert!(spec.contains("65536..=69631"));
    }
    for &(deck, sole_consumer) in MAX_HELDOUT_CONSUMERS {
        let consumers: Vec<_> = frozen
            .obligations()
            .iter()
            .filter(|row| row.decks().contains(&deck))
            .collect();
        assert_eq!(consumers.len(), 1, "max holdout {deck} needs one consumer");
        assert_eq!(consumers[0].leaf(), sole_consumer);
        assert_eq!(consumers[0].tier(), CampaignTier::Max);
        assert_governance_events(consumers[0]);
        let spec = authored_spec(fixture(frozen.fixtures(), deck));
        assert!(spec.contains("sole consumer"));
        assert!(spec.contains("131072..=135167"));
    }
    let industrial = frozen
        .obligations()
        .iter()
        .find(|row| row.leaf() == "i13-industrial-validation-max")
        .expect("industrial waiver consumer");
    assert_governance_events(industrial);
    let joint_relative_topology = frozen
        .obligations()
        .iter()
        .find(|row| row.leaf() == "i13-relative-topology-lineage-core")
        .expect("joint two-stratum relative-topology consumer");
    let joint_heldouts = joint_relative_topology
        .decks()
        .iter()
        .filter(|deck| held_out.contains(*deck))
        .copied()
        .collect::<BTreeSet<_>>();
    assert_eq!(
        joint_heldouts,
        BTreeSet::from([
            "i13-r3-relative-topology-core-holdout",
            "i13-topology-events-core-holdout",
        ]),
        "the joint consumer must retain both independently governed strata"
    );
    for required in [
        "claim_ids, leaf_ids, fixture_ids and stratum_ids are each nonempty, lexically raw-UTF-8-sorted and duplicate-free",
        "one group transition, never per-stratum transitions that permit split reveal",
    ] {
        assert!(
            policy_spec().contains(required),
            "joint multi-stratum governance lost '{required}'"
        );
    }

    let core_claims: BTreeSet<_> = SOLID_CLAIMS
        .iter()
        .chain(FRONTIER_CLAIMS)
        .copied()
        .filter(|id| *id != "i13-governed-industrial-machine-validation")
        .collect();
    let max_claims: BTreeSet<_> = MOONSHOT_CLAIMS
        .iter()
        .copied()
        .chain(["i13-governed-industrial-machine-validation"])
        .collect();
    for row in frozen.obligations() {
        let covered: BTreeSet<_> = row.claims_covered().iter().copied().collect();
        match row.tier() {
            CampaignTier::Core => {
                assert!(covered.is_subset(&core_claims));
                assert!(covered.is_disjoint(&max_claims));
                for &(deck, _) in MAX_HELDOUT_CONSUMERS {
                    assert!(
                        !row.decks().contains(&deck),
                        "core leaf consumed max holdout {deck}"
                    );
                }
            }
            CampaignTier::Max => {
                assert!(covered.is_subset(&max_claims));
                assert!(covered.is_disjoint(&core_claims));
                for &(deck, _) in CORE_HELDOUT_CONSUMERS {
                    assert!(
                        !row.decks().contains(&deck),
                        "max leaf consumed core holdout {deck}"
                    );
                }
            }
            CampaignTier::Smoke => panic!("I13 has no smoke adjudication owner"),
        }
    }
}

#[test]
#[allow(clippy::too_many_lines)]
fn i13_policy_keeps_execution_science_theorem_and_physical_axes_orthogonal() {
    let policy = policy_spec();
    for required in [
        "execution_disposition={Completed=0,Cancelled=1,BudgetExhausted=2,TimedOut=3,InfrastructureFailed=4}",
        "predicate_outcome={Satisfied,Violated,Indeterminate}",
        "claim_adjudication={Supported,Failed,Refuted,Unknown}",
        "evidence_completeness={CompleteEvidence,PartialEvidence,NoEvidence}",
        "evidence_integrity={IntegrityVerified,IntegrityFailed}",
        "domain_applicability={Admitted,OutOfDomain,Indeterminate,Malformed}",
        "operational_support={SupportedOperation,UnsupportedOperation}",
        "promotion_effect={Promotes,BlocksPromotion,NoPromotionAuthority}",
        "proof_status={Unattempted,Proved,Disproved,Unknown}",
        "physical_validation={NotAttempted,InDomainValidated,Failed,OutOfDomain,Unknown}",
        "falsifier_status={NotRun,Running,NoCountermodelWithinFrozenDomain,GenuineCountermodelFound,CampaignInvalid}",
        "No theorem or physical axis is inferred from another",
        "HOLDOUT_REALIZATION",
        "JOINT_REVEAL",
        "CANONICAL_GOVERNANCE_FRAMING_V1",
        "GOVERNED_GROUP_IDENTITY_V1",
        "GOVERNANCE_STAGES",
        "GOVERNANCE_GENESIS_V1",
        "GOVERNANCE_FIELD_INVENTORY_V1",
        "GOVERNANCE_VALUE_MATRIX_V1",
        "GOVERNANCE_OBSERVABILITY_V1",
        "INDUSTRIAL_TWO_PHASE_ADMISSION_V1",
        "GOVERNED_EXTERNAL_DISCHARGE",
        "TRANSACTION_INTENT_V1",
        "TRANSACTION_INTENT_ENCODING_V1",
        "GOVERNED_RECEIPT_AUTHORITY_V1",
        "MUTATION_IDEMPOTENCY_AND_FENCING",
        "MUTATION_FRONTIER_RECONCILIATION_V1",
        "DRAIN_TRIGGER_V1",
        "CANCELLATION_SCOPE_SEAL",
        "CANCELLATION_BOUNDS",
        "LIFECYCLE_EVENT_CONTRACT",
        "SCALE_BINDING",
        "ACCEPTANCE_ARITHMETIC",
        "UPSTREAM_BINDING",
        "INCONSISTENCY_RESOLUTION",
    ] {
        assert!(
            policy.contains(required),
            "policy lost orthogonal axis '{required}'"
        );
    }
    for required in [
        "Only GenuineCountermodel independently admitted against every exact frozen assumption refutes",
        "Cancelled, TimedOut, BudgetExhausted",
        "never map to Satisfied, Supported, or Refuted",
        "Philox provides replay only",
        "This numeric target is not an exhaustiveness claim",
        "semantic netlist and joint technology",
        "sole evidence-owning complete green",
        "A formal kernel proof is not inferred from held-out survival",
        "a held-out falsifier campaign is not a substitute for proof",
        "pre-access FrozenManifest::amend successor replaces the schema slot",
        "one atomic signed pre-access successor",
        "publishing a content hash alone grants no byte capability",
        "publicly derivable Philox stream is replay/development evidence only",
        "holdout-realization-commitment.v1",
        "U64LE(byte_len(alias_utf8))",
        "byte_len counts the exact case-sensitive UTF-8 bytes without normalization",
        "raw32(realization_key_32)",
        "raw32(generator_digest)",
        "raw32(case_order_digest)",
        "protocol_kind is the closed enum {HoldoutSlotRealization=1,ExternalDischarge=2}",
        "GovernanceCommitted=1, CandidateFrozen=2, RealizationCommitted=3, RevealedForAdjudication=4, and Closed=5",
        "protocol-specific transition grammars",
        "GOVERNANCE_VALUE_TABLE_V1",
        "claim_ids, leaf_ids, fixture_ids and stratum_ids are each nonempty, lexically raw-UTF-8-sorted and duplicate-free",
        "governed_group_root",
        "one group transition, never per-stratum transitions",
        "H3_RealizationCommitted",
        "E3_RealizationCommitted",
        "Every transition serializes all and only the twenty-five inventory fields in F01..F25 order",
        "D, P and N in each normative matrix row mean Digest, Pending and NotApplicable",
        "every protocol-relevant unresolved future output is P",
        "every protocol-inapplicable field is N from the first row onward",
        "A stage transition is authority-bearing only when committed by the named fs-vvreg authority",
        "GovernanceCommitted freezes the protected source universe",
        "CandidateFrozen binds the immutable candidate/model/toolchain/checker/AcceptanceCard",
        "For HoldoutSlotRealization, RealizationCommitted atomically installs every committed schema-slot successor",
        "For ExternalDischarge, RealizationCommitted is the single atomic envelope/waiver/protected-binding/amendment/authority-head transaction",
        "RevealedForAdjudication grants only its committed one-shot adjudication capability",
        "Closed follows independent adjudication",
        "explicit union {NotApplicable=0,Pending=1,Digest=2}",
        "Exact-key byte-identical replay is response behavior, not another stage outcome",
        "governance.transition and governance.protected_access are the entire stable governance.* namespace",
        "outcome={Committed,Refused,IntegrityFailed}",
        "A Committed transition event and authority-head CAS are one durable transaction",
        "replay_disposition=OriginalCommittedReceipt",
        "cannot encode IdempotentReplay as a transition outcome",
        "action={CapabilityIssue,AccessAttempt,ProtectedRead,DerivedOutput,Exclusion,Disclosure,CapabilityConsume,CapabilityRevoke}",
        "outcome={Prepared,Completed,Refused,Aborted}",
        "access_operation_id, idempotency key, predecessor access-chain receipt",
        "expected/observed/successor access-chain heads and generations",
        "terminal_effect_option={Absent=0,NotPerformed=1,Performed=2,MayHaveOccurred=3}",
        "Prepared is valid if and only if terminal_effect_option=Absent",
        "Prepared receipt identity is copied into exactly one terminal row",
        "atomically compare-and-swap the access-chain head and generation",
        "stale head, fork, ABA epoch or missing predecessor is IntegrityFailed",
        "Exact replay while Prepared returns the original Prepared receipt",
        "replay after any Completed, Refused or Aborted terminal returns that original terminal receipt",
        "crash-dangling Prepared is potential access and blocks authority",
        "Aborted with MayHaveOccurred",
        "Closed requires zero unresolved Prepared operations",
        "JSONL is a non-authoritative mirror",
        "I13IndustrialDischargeReceiptV1",
        "exact waiver subject and predicate",
        "publisher, standard, application-protocol, Part-21 edition",
        "failed-run/dependence/discrepancy rules",
        "simultaneous/multiplicity/stopping rule and adjudication code",
        "protocol-specific slot/envelope/waiver changes",
        "One atomic RealizationCommitted FrozenManifest::amend authority transaction",
        "compare-and-swap advances exactly the bound authority head",
        "org.frankensim.i13.governed-transaction-intent.v1",
        "canonical successor-intent projection for exactly one atomic protocol transaction",
        "predecessor governance-stage receipt",
        "immutable CandidateFrozen candidate/model/toolchain/checker/AcceptanceCard/QoI/scale/band/decision commitment",
        "role-addressed protected bindings",
        "independently sorted slot and root lists are forbidden",
        "initiative_id='I13'",
        "schema_identity='i13-governed-transaction-intent-v1'",
        "checked successor_version=predecessor_version+1",
        "governance_stage=RealizationCommitted",
        "prevents ABA replay",
        "retired waivers, protected roles and future artifacts form an exact role-addressed bijection",
        "pairwise distinct",
        "strict P -> authorization -> realized outputs -> successor -> AmendmentRecord -> NewAuthorityHead -> commit-receipt DAG",
        "25 exact ASCII bytes I13_TRANSACTION_INTENT_V1 followed by one trailing byte 0x00",
        "then U16LE(19)",
        "Digest is exactly 32 raw bytes and must be nonzero",
        "GovernanceCommitted=1,CandidateFrozen=2,RealizationCommitted=3,RevealedForAdjudication=4,Closed=5",
        "TRANSACTION_INTENT_COMPOUND_TYPES_V1",
        "Utf8OrEmpty=U64LE(byte_len)||exact case-sensitive bytes without normalization",
        "ProtectedBindingList=U64LE(count)||concat(FRAME_BYTES(ProtectedBinding))",
        "FutureArtifactList=U64LE(count)||concat(FRAME_BYTES(FutureArtifact))",
        "P is at most 16777216 bytes",
        "every Utf8/Utf8OrEmpty payload is at most 65536 bytes",
        "checked arithmetic against remaining input and these caps before allocation",
        "digest_union is exactly byte 0=Pending with no following bytes or byte 1=Digest followed by 32 raw bytes",
        "tags 0x000f..0x0011 require Pending",
        "Missing, extra, duplicate, out-of-order, unknown, empty-when-required, non-minimal or trailing bytes are refused",
        "two independent encoders must match a published 19-field intent KAT",
        "receipt_schema_set_root Digest",
        "ReceiptSchemaSetRoot=BLAKE3::derive_key('org.frankensim.i13.receipt-schema-set.v1'",
        "with exactly that case-sensitive role order and no alias",
        "RealizedOutputRoot=BLAKE3::derive_key('org.frankensim.i13.realized-output-set.v1'",
        "NewAuthorityHeadDigest=BLAKE3::derive_key('org.frankensim.i13.authority-head.v1'",
        "NewAuthorityHead=NewAuthorityHeadDigest||U64LE(checked_add(old_generation,1))",
        "atomically writes exactly the derived new head plus the non-embedded commit receipt",
        "unavailable/mismatched schema bytes, missing/extra/swapped role, output-set mismatch",
        "Every other mutating governance, reveal, adjudication, publication, resume or fork operation defines its own versioned domain-separated canonical request projection",
        "Exactly one fenced commit may win; losing branches are cancelled and fully drained",
        "prevents a content-hash cycle",
        "stale head/generation, permissive/unregistered protocol or output schema, receipt-role/schema mismatch, unavailable schema, role-swapped root, omitted CandidateFrozen input, raw or all-zero digest",
        "partial slot/waiver change",
        "split successor",
        "post-access amendment",
        "structurally frozen successor",
        "grants no discharge or promotion authority",
        "every submitted attempt emits exactly one execution.requested",
        "A refused admission emits no execution.admitted, execution.started, execution.drain.started or work event",
        "Every admitted attempt emits exactly one execution.admitted and then exactly one execution.started",
        "first canonically selected cancellation for one open scope/epoch is the primary request",
        "request-time observer-catalog root/count, descendant-frontier root/count, transition-frontier root/count, access-frontier root/count",
        "worst-case terminal row-and-encoded-byte reservation",
        "later distinct request for the sealed/draining/finalized epoch receives a deterministic AlreadySealed/AlreadyDraining/AlreadyFinalized response",
        "When CancellationSeal wins terminal-trigger selection, its atomically appended execution.drain.started initiates propagation exactly once",
        "Every sealed observer emits exactly one cancellation.observed binding the primary request event identity",
        "independently reconciled cardinality-zero root and exactly zero observations",
        "request-to-observation maximum is defined as zero at seal",
        "one supervisor monotonic clock",
        "conservative outward uncertainty and drift enclosure",
        "execution.cancelled is atomically appended by the terminal-selection transaction only when CancellationSeal wins",
        "CanonicalTerminalPath atomically selects exactly one DRAIN_TRIGGER_V1 trigger and emits exactly one execution.drain.started",
        "then emits exactly one execution.drained and performs exactly one indivisible TERMINAL_CONDITION_FRONTIER_V1 cutoff/execution.finalized transaction under TERMINAL_DEADLINE_GATE_V1",
        "TerminalPersistenceCatastrophe is permitted only when the reserved priority terminal sink, the independently reserved emergency terminal journal and every authenticated durable-prefix recovery route all fail",
        "the attempt is TerminalPersistencePending and is neither finalized nor catastrophic",
        "successful recovery enters CanonicalTerminalPath rather than inventing a second branch",
        "derives no canonical execution disposition or evidence-completeness axis including NoEvidence",
        "out-of-band explicitly non-authoritative incident diagnostic with zero scientific, theorem, optimization or promotion authority",
        "execution.finalized binds final execution/input/integrity/domain/support state, selected trigger, TerminalConditionBitsV1, TerminalConditionCutoffV1, TerminalConditionWitnessSetRoot/count, deadline-resolution root/count, derived final disposition and immutable evidence payload",
        "exactly one evidence.adjudication_receipt after finalization binds claim/evidence axes",
        "Successful adjudication requires exactly one evidence.atomic_publish_committed followed by exactly one artifact.published",
        "may publish at most one paired atomic terminal/failure evidence transaction",
        "with zero promotion authority",
        "partial/unpaired publication is IntegrityFailed",
        "Duplicate, missing, out-of-order or inapplicable lifecycle events",
        "work/cancellation/drain after execution.finalized",
        "A missing or late sealed observation contributes an SLO-overrun/TimedOut-bit witness",
        "a receipt-qualified observer, supervisor or sink loss contributes an InfrastructureFailed-bit witness",
        "simultaneous conditions contribute both",
        "yields evidence_integrity=IntegrityFailed",
        "On CanonicalTerminalPath, evidence completeness is PartialEvidence only when a durable authenticated prefix/FailureBundle exists and NoEvidence otherwise",
        "neither state can become scientific Failed or Refuted",
        "observation_latency_i=checked_sub(observation_interval_upper_i,request_interval_lower)",
        "max_i(observation_latency_i) is compared inclusively with 250000000 Core or 1000000000 Max",
        "complete observer-observation and descendant-exit reconciliation roots",
        "membership-plus-exit receipt",
        "calibration artifact and watchdog arrival remain telemetry",
        "equality_miss=abs(observed-expected)/S",
        "upper_miss=max(0,observed-upper)/S",
        "lower_miss=max(0,lower-observed)/S",
        "interval_miss=max(0,lo-observed,observed-hi)/S",
        "can never be imputed as zero or omitted",
        "canonical upstream set as exact ClaimSpec digests",
        "maximal falsifier additionally binds the canonical digest of every targeted ClaimSpec",
        "request observation <=250 ms",
        "IntegrityFailed contradiction",
    ] {
        assert!(
            policy.contains(required),
            "policy lost no-collapse rule '{required}'"
        );
    }
    let intent_fields = [
        ("0x0001", "initiative_id", "Utf8"),
        ("0x0002", "schema_identity", "Utf8"),
        ("0x0003", "protocol_kind", "U16LE"),
        ("0x0004", "successor_version", "U64LE"),
        ("0x0005", "predecessor_manifest_digest", "Digest"),
        ("0x0006", "expected_authority_head", "AuthorityHead"),
        ("0x0007", "predecessor_stage_receipt_digest", "Digest"),
        ("0x0008", "candidate_freeze_commitment_digest", "Digest"),
        ("0x0009", "retired_waiver_subjects", "Utf8List"),
        ("0x000a", "protected_bindings", "ProtectedBindingList"),
        ("0x000b", "coupled_transaction_group", "Utf8"),
        ("0x000c", "governance_stage", "U16LE"),
        ("0x000d", "authority_scope", "Utf8"),
        ("0x000e", "mutation_fence", "MutationFence"),
        ("0x000f", "successor_artifact_slots", "FutureArtifactList"),
        ("0x0010", "final_successor_slot", "FutureDigest"),
        ("0x0011", "amendment_record_slot", "FutureDigest"),
        ("0x0012", "receipt_schema_set_root", "Digest"),
        ("0x0013", "governance_protocol_schema_set_root", "Digest"),
    ];
    assert_eq!(intent_fields.len(), 19);
    let intent_header = b"I13_TRANSACTION_INTENT_V1\0";
    assert_eq!(intent_header.len(), 26);
    assert_eq!(intent_header.last(), Some(&0));
    assert!(!intent_header.ends_with(br"\0"));
    assert!(!policy.contains(r"I13_TRANSACTION_INTENT_V1\0"));
    let mut previous_mapping_end = 0;
    for (tag, name, field_type) in intent_fields {
        let exact_mapping = format!("{tag} {name} {field_type}");
        assert_eq!(
            policy.match_indices(exact_mapping.as_str()).count(),
            1,
            "transaction-intent mapping must occur exactly once: {exact_mapping}"
        );
        let mapping_start = policy
            .find(exact_mapping.as_str())
            .expect("exact transaction-intent mapping");
        assert!(
            mapping_start >= previous_mapping_end,
            "transaction-intent mappings must remain in increasing tag order: {exact_mapping}"
        );
        previous_mapping_end = mapping_start + exact_mapping.len();
    }
    let semantic_authority_header = b"I13_OUTPUT_SCHEMA_SEMANTIC_AUTHORITY_V1\0";
    let output_artifact_schema_header = b"I13_OUTPUT_ARTIFACT_SCHEMA_V1\0";
    let external_body_schema_header = b"I13_EXTERNAL_DISCHARGE_BODY_SCHEMA_V1\0";
    assert_eq!(semantic_authority_header.len(), 40);
    assert_eq!(output_artifact_schema_header.len(), 30);
    assert_eq!(external_body_schema_header.len(), 38);
    assert_eq!(40 + 2 + 2 * 4, 50);
    for header in [
        semantic_authority_header.as_slice(),
        output_artifact_schema_header.as_slice(),
        external_body_schema_header.as_slice(),
    ] {
        assert_eq!(header.last(), Some(&0));
        assert!(!header.ends_with(br"\0"));
    }
    let governed_output_fields = [
        (1, "schema_identity", "Utf8", "ConstantEnvelopeIdentity"),
        (2, "protocol_kind", "U16LE", "MatrixProtocolKind"),
        (3, "governance_stage", "U16LE", "MatrixGovernanceStage"),
        (4, "authority_scope", "Utf8", "MatrixAuthorityScope"),
        (5, "target_slot_id", "Utf8", "FutureTargetSlotId"),
        (6, "artifact_role", "Utf8", "DerivedArtifactRole"),
        (
            7,
            "related_waiver_subject",
            "Utf8OrEmpty",
            "FutureRelatedWaiverSubject",
        ),
        (
            8,
            "transaction_intent_digest",
            "Digest",
            "TransactionIntentDigest",
        ),
        (
            9,
            "protocol_authorization_receipt_digest",
            "Digest",
            "ProtocolAuthorizationReceiptDigest",
        ),
        (
            10,
            "governance_protocol_schema_set_root",
            "Digest",
            "GovernanceProtocolSchemaSetRoot",
        ),
        (
            11,
            "candidate_freeze_commitment_digest",
            "Digest",
            "IntentCandidateFreezeDigest",
        ),
        (12, "protected_root", "Digest", "ProtectedBindingRoot"),
        (
            13,
            "realized_payload_schema_digest",
            "Digest",
            "ProtocolPayloadSchemaDigest",
        ),
        (
            14,
            "realized_payload_digest",
            "Digest",
            "RehashedPayloadDigest",
        ),
        (
            15,
            "realized_payload_byte_count",
            "U64LE",
            "StreamedPayloadByteCount",
        ),
        (
            16,
            "realization_receipt_root",
            "Digest",
            "RoleSpecificRealizationReceiptRoot",
        ),
    ];
    assert_eq!(governed_output_fields.len(), 16);
    let governed_output_name_bytes: usize = governed_output_fields
        .iter()
        .map(|(_, name, _, _)| name.len())
        .sum();
    assert_eq!(governed_output_name_bytes, 357);
    assert_eq!(30 + 32 + 2 + 32 + 2 + 16 * 21 + 357, 791);
    for (ordinal, name, field_type, constraint) in governed_output_fields {
        let exact_row = format!("{ordinal}/{name}/{field_type}/{constraint}");
        assert!(
            policy.contains(exact_row.as_str()),
            "governed output schema lost exact row '{exact_row}'"
        );
    }
    let external_body_fields = [
        (1, "schema_identity", "Utf8"),
        (2, "waiver_subject", "Utf8"),
        (3, "waiver_predicate_digest", "Digest"),
        (4, "protected_binding_root", "Digest"),
        (5, "acquisition_access_custody_digest", "Digest"),
        (6, "publisher_standard_profile_digest", "Digest"),
        (7, "raw_evidence_byte_root", "Digest"),
        (8, "canonical_normalized_root", "Digest"),
        (9, "occurrence_configuration_netlist_root", "Digest"),
        (10, "material_process_lineage_root", "Digest"),
        (11, "coordinate_unit_frame_time_root", "Digest"),
        (12, "instrumentation_calibration_uncertainty_root", "Digest"),
        (13, "environment_control_source_state_root", "Digest"),
        (14, "operating_domain_rules_root", "Digest"),
        (15, "missingness_dependence_discrepancy_root", "Digest"),
        (16, "qoi_scale_band_decision_root", "Digest"),
        (
            17,
            "candidate_model_toolchain_checker_acceptance_root",
            "Digest",
        ),
        (18, "ownership_threat_graph_root", "Digest"),
        (19, "redaction_disclosure_policy_digest", "Digest"),
        (20, "review_signature_set_root", "Digest"),
        (21, "revocation_evidence_root", "Digest"),
        (22, "complete_access_attempt_ledger_root", "Digest"),
        (23, "predecessor_manifest_digest", "Digest"),
        (24, "transaction_intent_digest", "Digest"),
        (25, "authorization_receipt_digest", "Digest"),
        (26, "no_claim_boundary_digest", "Digest"),
        (27, "realized_evidence_set_root", "Digest"),
    ];
    assert_eq!(external_body_fields.len(), 27);
    let external_body_name_bytes: usize = external_body_fields
        .iter()
        .map(|(_, name, _)| name.len())
        .sum();
    assert_eq!(external_body_name_bytes, 783);
    assert_eq!(38 + 32 + 2 + 27 * 21 + 783, 1_422);
    let nested_semantic_authority_header = b"I13_NESTED_SEMANTIC_AUTHORITY_V1\0";
    assert_eq!(nested_semantic_authority_header.len(), 33);
    assert_eq!(33 + 2 + 20 * 129, 2_615);
    for (ordinal, name, field_type) in external_body_fields {
        let exact_row = format!("{ordinal}/{name}/{field_type}");
        assert!(
            policy.contains(exact_row.as_str()),
            "external discharge body schema lost exact row '{exact_row}'"
        );
    }
    let authorization_schema_fields = [
        (1, "receipt_role", "U8"),
        (2, "receipt_schema_digest", "Digest"),
        (3, "receipt_schema_set_root", "Digest"),
        (4, "receipt_schema_membership_ordinal", "U8"),
        (5, "transaction_intent_digest", "Digest"),
        (6, "predecessor_stage_receipt_digest", "Digest"),
        (7, "expected_authority_head_digest", "Digest"),
        (8, "expected_authority_head_generation", "U64LE"),
        (9, "governance_protocol_schema_set_root", "Digest"),
        (10, "protocol_kind", "U16LE"),
        (11, "governance_stage", "U16LE"),
        (12, "authority_scope", "Utf8"),
        (13, "future_artifact_count", "U64LE"),
        (14, "output_schema_membership_proof_set_root", "Digest"),
    ];
    let commit_schema_fields = [
        (1, "receipt_role", "U8"),
        (2, "receipt_schema_digest", "Digest"),
        (3, "receipt_schema_set_root", "Digest"),
        (4, "receipt_schema_membership_ordinal", "U8"),
        (5, "transaction_intent_digest", "Digest"),
        (6, "protocol_authorization_receipt_digest", "Digest"),
        (7, "realized_output_root", "Digest"),
        (8, "predecessor_manifest_digest", "Digest"),
        (9, "final_successor_digest", "Digest"),
        (10, "amendment_record_digest", "Digest"),
        (11, "governance_protocol_schema_set_root", "Digest"),
        (12, "observed_old_authority_head_digest", "Digest"),
        (13, "observed_old_authority_head_generation", "U64LE"),
        (14, "new_authority_head_digest", "Digest"),
        (15, "new_authority_head_generation", "U64LE"),
        (16, "cas_result", "U8"),
    ];
    assert_eq!(authorization_schema_fields.len(), 14);
    assert_eq!(commit_schema_fields.len(), 16);
    for (ordinal, name, field_type) in authorization_schema_fields
        .into_iter()
        .chain(commit_schema_fields)
    {
        let exact_row = format!("{ordinal}/{name}/{field_type}");
        assert!(
            policy.contains(exact_row.as_str()),
            "receipt schema lost exact row '{exact_row}'"
        );
    }
    let authorization_schema_name_bytes: usize = authorization_schema_fields
        .iter()
        .map(|(_, name, _)| name.len())
        .sum();
    let commit_schema_name_bytes: usize = commit_schema_fields
        .iter()
        .map(|(_, name, _)| name.len())
        .sum();
    assert_eq!(authorization_schema_name_bytes, 349);
    assert_eq!(commit_schema_name_bytes, 414);
    assert_eq!(16 + 36 + 2 + 14 * 21 + 349, 697);
    assert_eq!(16 + 37 + 2 + 16 * 21 + 414, 805);
    let holdout_realization_header = b"I13_HOLDOUT_REALIZATION_RECEIPT_V1\0";
    let external_realization_header = b"I13_EXTERNAL_REALIZATION_RECEIPT_V1\0";
    let authorization_receipt_header = b"I13_PROTOCOL_AUTHORIZATION_V1\0";
    let commit_receipt_header = b"I13_PROTOCOL_COMMIT_RECEIPT_V1\0";
    let governed_output_header = b"I13_GOVERNED_OUTPUT_ENVELOPE_V1\0";
    let external_body_header = b"I13_EXTERNAL_DISCHARGE_BODY_V1\0";
    let schema_proof_set_header = b"I13_OUTPUT_SCHEMA_PROOF_SET_V1\0";
    assert_eq!(holdout_realization_header.len(), 35);
    assert_eq!(external_realization_header.len(), 36);
    assert_eq!(authorization_receipt_header.len(), 30);
    assert_eq!(commit_receipt_header.len(), 31);
    assert_eq!(governed_output_header.len(), 32);
    assert_eq!(external_body_header.len(), 31);
    assert_eq!(schema_proof_set_header.len(), 31);
    assert_eq!(507 + 2 * 65_536, 131_579);
    assert_eq!(388 + 3 * 65_536, 196_996);
    assert_eq!(292 + 65_536, 65_828);
    assert_eq!(31 + 3 + 11 * 32 + 2 * 8, 402);
    assert_eq!(507 + 65_513 + 65_536, 131_556);
    assert_eq!(388 + 65_536 + 65_536 + 65_518, 196_978);
    assert_eq!(292 + 28, 320);
    assert_eq!(292 + 22, 314);
    assert_eq!(8 + (8 + 8 + 21 + 32) + (8 + 8 + 14 + 32), 139);
    assert_eq!(1 + 1 + 32, 34);
    for exact in [
        "GOVERNED_PROTOCOL_SCHEMA_AUTHORITY_V1",
        "33 exact ASCII bytes I13_GOVERNANCE_PROTOCOL_MATRIX_V1 followed by one byte 0x00",
        "U16LE(2), and exactly two rows in protocol-kind order",
        "MatrixRowV1=U16LE(protocol_kind)||U16LE(governance_stage)",
        "HoldoutSlotRealization=1,RealizationCommitted=3",
        "ExternalDischarge=2,RealizationCommitted=3",
        "GovernanceProtocolMatrixDigest=BLAKE3::derive_key('org.frankensim.i13.governance-protocol-matrix.v1'",
        "GovernanceProtocolSchemaMemberV1=U16LE(protocol_kind)",
        "GovernanceProtocolSchemaSetProjectionV1=raw32(GovernanceProtocolMatrixDigest)",
        "GovernanceProtocolSchemaMembershipProofV1=U64LE(ordinal)||FRAME_BYTES(GovernanceProtocolSchemaMemberV1)",
        "an isolated record is never a membership proof",
        "tag 0x0013 byte-equals GovernanceProtocolSchemaSetRoot",
        "never the set root itself",
        "matrix plus all fetched schema bytes and the complete set projection together are at most 67108864 bytes",
        "cumulative decoded constraint-AST nodes are at most 1048576",
        "depth 64 is admitted while 65 is refused before recursion/allocation",
        "ordinal/member/full-rehash membership",
        "permissive-schema mutation twins before governance authority",
        "GOVERNED_OUTPUT_SCHEMA_ROLE_CONFORMANCE_V1",
        "39 exact ASCII bytes I13_OUTPUT_SCHEMA_SEMANTIC_AUTHORITY_V1 followed by one byte 0x00",
        "the complete projection is exactly 50 bytes",
        "HoldoutSlotRealization=1/HoldoutRealizationEnvelope=1",
        "ExternalDischarge=2/ExternalDischargeEnvelope=2",
        "The semantic-authority row selected by protocol_kind supplies schema_kind; the protocol matrix alone does not",
        "29 ASCII bytes I13_OUTPUT_ARTIFACT_SCHEMA_V1",
        "raw32(required_protocol_payload_schema_digest)",
        "exactly 791 bytes",
        "U16LE(16), and sixteen SchemaFieldV1 records in ordinal order",
        "requiredness is exactly Required=1",
        "schema_identity is exactly i13-governed-output-envelope-v1",
        "fields 2..12 byte-equal their matrix, FutureArtifact, P, authorization receipt and matching ProtectedBinding authorities",
        "field 15 is in 1..=U64::MAX",
        "governance cannot invent, weaken or substitute that payload schema",
        "37 exact ASCII bytes I13_EXTERNAL_DISCHARGE_BODY_SCHEMA_V1 followed by one byte 0x00",
        "U16LE(27), and exactly these required ordinal/name/type semantic rows with no optional or extension field",
        "raw32(NestedSemanticAuthorityRoot), U16LE(27)",
        "exactly 1422 bytes under EXTERNAL_DISCHARGE_BODY_SCHEMA_ENCODING_V1",
        "exact 2615-byte nested authority projection/root",
        "every member and the set root transitively bind OutputSchemaSemanticAuthorityRoot, schema kind, the closed wrapper table and exact protocol payload-schema authority",
        "initially permissive, optional, omitted, renamed, retyped, extended",
        "initially-permissive mutation before governance authority",
        "EXTERNAL_DISCHARGE_BODY_SCHEMA_ENCODING_V1",
        "ExternalBodySchemaFieldV1=U16LE(ordinal)||FRAME_BYTES(Utf8(field_name))",
        "I13ExternalDischargeBodySchemaBytesV1 is the 38-byte header/NUL followed by raw32(NestedSemanticAuthorityRoot)",
        "exactly 1422 bytes",
        "27 fixed 21-byte descriptor overheads and 783 field-name bytes",
        "There is no separately caller-supplied nested-schema digest",
        "NESTED_SEMANTIC_AUTHORITY_V1",
        "exactly twenty 129-byte members",
        "raw32(independent_decoder_set_root)",
        "two distinct nonzero fingerprints in lexical order",
        "fixed by the predecessor manifest before GovernanceCommitted",
        "CandidateFrozen subsequently copies the identical root before candidate execution",
        "GOVERNED_REALIZATION_RECEIPTS_V1",
        "org.frankensim.i13.realized-payload.v1",
        "34 exact ASCII bytes I13_HOLDOUT_REALIZATION_RECEIPT_V1 followed by one byte 0x00",
        "at most 131579 bytes",
        "35 exact ASCII bytes I13_EXTERNAL_REALIZATION_RECEIPT_V1 followed by one byte 0x00",
        "at most checked_add(388,checked_mul(3,65536))=196996 bytes",
        "The semantic-authority row selected by protocol_kind—not the protocol matrix alone—selects schema kind",
        "GOVERNED_RECEIPT_SCHEMA_CONFORMANCE_V1",
        "ProtocolAuthorizationSchemaBytes is exactly FRAME_BYTES(Utf8('i13-protocol-authorization-schema-v1'))||U16LE(14)",
        "exactly 697 bytes",
        "29 exact ASCII bytes I13_PROTOCOL_AUTHORIZATION_V1 followed by one byte 0x00",
        "exactly checked_add(292,byte_len(authority_scope_utf8)) bytes and at most 65828 bytes",
        "ProtocolCommitReceiptSchemaBytes is exactly FRAME_BYTES(Utf8('i13-protocol-commit-receipt-schema-v1'))||U16LE(16)",
        "exactly 805 bytes",
        "30 exact ASCII bytes I13_PROTOCOL_COMMIT_RECEIPT_V1 followed by one byte 0x00",
        "exactly 402 bytes",
        "hash-correct initially permissive receipt schema",
        "GOVERNED_OUTPUT_INSTANCE_ENCODING_V1",
        "31 exact ASCII bytes I13_GOVERNED_OUTPUT_ENVELOPE_V1 followed by one byte 0x00",
        "governed_output_envelope_digest=BLAKE3::derive_key('org.frankensim.i13.governed-output-envelope.v1'",
        "sole realized_artifact_digest used by FutureDigest, the successor slot and RealizedOutputRoot",
        "30 exact ASCII bytes I13_EXTERNAL_DISCHARGE_BODY_V1 followed by one byte 0x00",
        "OUTPUT_SCHEMA_MEMBERSHIP_PROOF_SET_V1",
        "30 exact ASCII bytes I13_OUTPUT_SCHEMA_PROOF_SET_V1 followed by one byte 0x00",
        "output_schema_membership_proof_set_root=BLAKE3::derive_key('org.frankensim.i13.output-schema-proof-set.v1'",
        "GOVERNANCE_SCHEMA_SET_CLOSURE_V1",
        "lexicographically sorted, duplicate-free and all-and-only the derived role set",
        "schema_count is in 1..=4096 and byte-equals the protected-binding count, the future-artifact count",
        "GOVERNED_RECEIPT_IDENTITY_V1",
        "ReceiptSchemaSetProjectionV1=U64LE(2)",
        "exactly 139 bytes",
        "ReceiptSchemaMembershipProofV1 is exactly one byte ordinal||one byte receipt_role||raw32(receipt_schema_digest), 34 bytes",
        "protocol_authorization_receipt_digest=BLAKE3::derive_key('org.frankensim.i13.protocol-authorization-receipt.v1'",
        "protocol_commit_receipt_digest=BLAKE3::derive_key('org.frankensim.i13.protocol-commit-receipt.v1'",
        "ROLE_DERIVED_RECEIPT_BOUNDS_V1",
        "protocol-valid maxima are 131556 for HoldoutRealizationReceiptBytesV1, 196978 for ExternalRealizationReceiptBytesV1",
        "respectively 320 or 314 for ProtocolAuthorizationReceiptBytesV1",
    ] {
        assert!(
            policy.contains(exact),
            "governed protocol schema authority lost '{exact}'"
        );
    }
}

#[test]
#[allow(clippy::too_many_lines)]
fn i13_governance_schema_has_exact_framing_genesis_and_closed_matrix() {
    let policy = policy_spec();
    let fields = [
        ("F01", "source_universe"),
        ("F02", "acquisition_access_custody"),
        ("F03", "authorized_principals"),
        ("F04", "candidate_input_permissions"),
        ("F05", "selection_exclusion"),
        ("F06", "receipt_schema"),
        ("F07", "decision_rule"),
        ("F08", "candidate_freeze"),
        ("F09", "schema_slot_successors"),
        ("F10", "realized_root_commitments"),
        ("F11", "hidden_opening_commitments"),
        ("F12", "protected_bindings"),
        ("F13", "retired_waivers"),
        ("F14", "discharge_envelopes"),
        ("F15", "pre_reveal_access_chain_root"),
        ("F16", "amendment_record"),
        ("F17", "realization_authority_head"),
        ("F18", "one_shot_capability"),
        ("F19", "reveal_receipt"),
        ("F20", "reveal_access_chain_head"),
        ("F21", "adjudication_receipt"),
        ("F22", "final_access_chain_root"),
        ("F23", "action_outcome_effect_counts"),
        ("F24", "terminal_receipt"),
        ("F25", "failure_bundle_or_clean"),
    ];
    let all_fields: BTreeSet<_> = fields.iter().map(|(tag, _)| *tag).collect();
    assert_eq!(all_fields.len(), 25);

    let mut previous_mapping_end = 0;
    for (tag, name) in fields {
        let mapping = format!("{tag}={name}");
        assert_eq!(
            policy.match_indices(mapping.as_str()).count(),
            1,
            "governance field mapping must occur exactly once: {mapping}"
        );
        let mapping_start = policy.find(mapping.as_str()).expect("field mapping");
        assert!(
            mapping_start >= previous_mapping_end,
            "governance fields must remain in F01..F25 order: {mapping}"
        );
        previous_mapping_end = mapping_start + mapping.len();
    }

    type MatrixRow = (
        &'static str,
        &'static [&'static str],
        &'static [&'static str],
        &'static [&'static str],
    );
    let rows: [MatrixRow; 10] = [
        (
            "H1_GovernanceCommitted",
            &["F01", "F02", "F03", "F04", "F05", "F06", "F07"],
            &[
                "F08", "F09", "F10", "F11", "F12", "F15", "F16", "F17", "F18", "F19", "F20", "F21",
                "F22", "F23", "F24", "F25",
            ],
            &["F13", "F14"],
        ),
        (
            "H2_CandidateFrozen",
            &["F01", "F02", "F03", "F04", "F05", "F06", "F07", "F08"],
            &[
                "F09", "F10", "F11", "F12", "F15", "F16", "F17", "F18", "F19", "F20", "F21", "F22",
                "F23", "F24", "F25",
            ],
            &["F13", "F14"],
        ),
        (
            "H3_RealizationCommitted",
            &[
                "F01", "F02", "F03", "F04", "F05", "F06", "F07", "F08", "F09", "F10", "F11", "F12",
                "F15", "F16", "F17",
            ],
            &["F18", "F19", "F20", "F21", "F22", "F23", "F24", "F25"],
            &["F13", "F14"],
        ),
        (
            "H4_RevealedForAdjudication",
            &[
                "F01", "F02", "F03", "F04", "F05", "F06", "F07", "F08", "F09", "F10", "F11", "F12",
                "F15", "F16", "F17", "F18", "F19", "F20",
            ],
            &["F21", "F22", "F23", "F24", "F25"],
            &["F13", "F14"],
        ),
        (
            "H5_Closed",
            &[
                "F01", "F02", "F03", "F04", "F05", "F06", "F07", "F08", "F09", "F10", "F11", "F12",
                "F15", "F16", "F17", "F18", "F19", "F20", "F21", "F22", "F23", "F24", "F25",
            ],
            &[],
            &["F13", "F14"],
        ),
        (
            "E1_GovernanceCommitted",
            &["F01", "F02", "F03", "F04", "F05", "F06", "F07"],
            &[
                "F08", "F12", "F13", "F14", "F15", "F16", "F17", "F18", "F19", "F20", "F21", "F22",
                "F23", "F24", "F25",
            ],
            &["F09", "F10", "F11"],
        ),
        (
            "E2_CandidateFrozen",
            &["F01", "F02", "F03", "F04", "F05", "F06", "F07", "F08"],
            &[
                "F12", "F13", "F14", "F15", "F16", "F17", "F18", "F19", "F20", "F21", "F22", "F23",
                "F24", "F25",
            ],
            &["F09", "F10", "F11"],
        ),
        (
            "E3_RealizationCommitted",
            &[
                "F01", "F02", "F03", "F04", "F05", "F06", "F07", "F08", "F12", "F13", "F14", "F15",
                "F16", "F17",
            ],
            &["F18", "F19", "F20", "F21", "F22", "F23", "F24", "F25"],
            &["F09", "F10", "F11"],
        ),
        (
            "E4_RevealedForAdjudication",
            &[
                "F01", "F02", "F03", "F04", "F05", "F06", "F07", "F08", "F12", "F13", "F14", "F15",
                "F16", "F17", "F18", "F19", "F20",
            ],
            &["F21", "F22", "F23", "F24", "F25"],
            &["F09", "F10", "F11"],
        ),
        (
            "E5_Closed",
            &[
                "F01", "F02", "F03", "F04", "F05", "F06", "F07", "F08", "F12", "F13", "F14", "F15",
                "F16", "F17", "F18", "F19", "F20", "F21", "F22", "F23", "F24", "F25",
            ],
            &[],
            &["F09", "F10", "F11"],
        ),
    ];

    let mut prior_holdout_digests = BTreeSet::new();
    let mut prior_external_digests = BTreeSet::new();
    for (row_name, digests, pending, not_applicable) in rows {
        let expected = format!(
            "{row_name} D={{{}}};P={{{}}};N={{{}}}.",
            digests.join(","),
            pending.join(","),
            not_applicable.join(",")
        );
        assert_eq!(
            policy.match_indices(expected.as_str()).count(),
            1,
            "normative governance matrix row must occur exactly once: {expected}"
        );

        let digest_set: BTreeSet<_> = digests.iter().copied().collect();
        let pending_set: BTreeSet<_> = pending.iter().copied().collect();
        let not_applicable_set: BTreeSet<_> = not_applicable.iter().copied().collect();
        assert_eq!(digest_set.len(), digests.len(), "duplicate D in {row_name}");
        assert_eq!(
            pending_set.len(),
            pending.len(),
            "duplicate P in {row_name}"
        );
        assert_eq!(
            not_applicable_set.len(),
            not_applicable.len(),
            "duplicate N in {row_name}"
        );
        assert!(
            digest_set.is_disjoint(&pending_set),
            "D/P overlap in {row_name}"
        );
        assert!(
            digest_set.is_disjoint(&not_applicable_set),
            "D/N overlap in {row_name}"
        );
        assert!(
            pending_set.is_disjoint(&not_applicable_set),
            "P/N overlap in {row_name}"
        );
        let partition: BTreeSet<_> = digest_set
            .iter()
            .chain(&pending_set)
            .chain(&not_applicable_set)
            .copied()
            .collect();
        assert_eq!(
            partition, all_fields,
            "matrix row is not closed: {row_name}"
        );

        let prior = if row_name.starts_with('H') {
            &mut prior_holdout_digests
        } else {
            &mut prior_external_digests
        };
        assert!(
            prior.is_subset(&digest_set),
            "a prior digest regressed in {row_name}"
        );
        *prior = digest_set;
    }

    let group_header = b"I13_GOVERNED_GROUP_V1\0";
    assert_eq!(group_header.len(), 22);
    assert_eq!(group_header.last(), Some(&0));
    assert!(!group_header.ends_with(br"\0"));
    assert!(!policy.contains(r"I13_GOVERNED_GROUP_V1\0"));
    for exact in [
        "FRAME_BYTES(x_bytes)=U64LE(byte_len(x_bytes))||x_bytes",
        "FRAME_BYTES(Utf8(s)) intentionally carries an outer encoded-payload length and the inner UTF-8 byte length",
        "governed_group_root=BLAKE3::derive_key('org.frankensim.i13.governed-group-root.v1',G)",
        "Each IdSet has at most 4096 identifiers",
        "every identifier is nonempty",
        "each identifier has at most 65536 UTF-8 bytes",
        "G has at most 4194304 bytes under checked arithmetic",
        "two independently implemented encoders must match a published KAT containing exact G bytes and governed_group_root",
        "Uninitialized=0 exists only as the from_stage of the first authority transition",
        "transition_genesis_receipt=BLAKE3::derive_key('org.frankensim.i13.transition-genesis.v1'",
        "access_genesis_receipt=BLAKE3::derive_key('org.frankensim.i13.access-genesis.v1'",
        "advances the authority-head generation from 0 to 1",
        "advances access generation 0 to 1",
        "terminal_effect_option={Absent=0,NotPerformed=1,Performed=2,MayHaveOccurred=3}",
        "Prepared is valid if and only if terminal_effect_option=Absent",
        "Closed requires zero unresolved Prepared operations",
        "one governed_group_root, one genesis-bound capability epoch, and one monotone governance/authority/access lineage",
    ] {
        assert!(policy.contains(exact), "governance schema lost '{exact}'");
    }
}

#[test]
fn i13_drain_and_logging_arithmetic_is_exact_and_end_to_end() {
    let policy = policy_spec();
    let trigger_candidate_header = b"I13_TRIGGER_CANDIDATE_V1\0";
    assert_eq!(trigger_candidate_header.len(), 25);
    assert_eq!(trigger_candidate_header.last(), Some(&0));
    let cancellation_proposal_header = b"I13_CANCELLATION_PROPOSAL_V1\0";
    assert_eq!(cancellation_proposal_header.len(), 29);
    assert_eq!(cancellation_proposal_header.last(), Some(&0));
    let unavoidable_source_universe_header = b"I13_UNAVOIDABLE_SOURCE_UNIVERSE_V1\0";
    let candidate_source_class_header = b"I13_CANDIDATE_SOURCE_CLASS_MATRIX_V1\0";
    let candidate_breach_header = b"I13_CANDIDATE_CAPACITY_BREACH_V1\0";
    let candidate_breach_selection_header = b"I13_CAPACITY_BREACH_SELECTION_V1\0";
    let post_selection_class_header = b"I13_POST_SELECTION_SOURCE_CLASS_MATRIX_V1\0";
    let post_selection_breach_header = b"I13_POST_SELECTION_CAPACITY_BREACH_V1\0";
    assert_eq!(unavoidable_source_universe_header.len(), 35);
    assert_eq!(candidate_source_class_header.len(), 37);
    assert_eq!(candidate_breach_header.len(), 33);
    assert_eq!(candidate_breach_selection_header.len(), 33);
    assert_eq!(post_selection_class_header.len(), 42);
    assert_eq!(post_selection_breach_header.len(), 38);
    assert_eq!(35 + 2 + 256 * 43, 11_045);
    assert_eq!(276 + 1_048_000 + 76, 1_048_352);
    assert_eq!(493 + 8, 501);
    assert_eq!(322 + 1_048_000, 1_048_322);
    let terminal_condition_header = b"I13_TERMINAL_CONDITION_V1\0";
    assert_eq!(terminal_condition_header.len(), 26);
    assert_eq!(26 + 1 + 8 + 32, 67);
    let terminal_deadline_ticket_header = b"I13_TERMINAL_DEADLINE_TICKET_V1\0";
    assert_eq!(terminal_deadline_ticket_header.len(), 32);
    let terminal_deadline_resolution_header = b"I13_TERMINAL_DEADLINE_RESOLUTION_V1\0";
    assert_eq!(terminal_deadline_resolution_header.len(), 36);
    assert_eq!(32 + 1 + 8 + 8 + 32 + 32 + 8 + 9 + 32 + 8 + 32, 202);
    assert_eq!(32 + 1 + 8 + 8 + 32 + 32 + 8 + 17 + 32 + 8 + 32, 210);
    assert_eq!(36 + 1 + 1 + 32 + 32, 102);
    assert_eq!(32 + 32 + 8 + 2 * (8 + 102), 292);
    assert_eq!(16 + 84 * 65_536, 5_505_040);
    assert_eq!(5_505_040 + 1, 5_505_041);
    assert_eq!(40 + 84 * 65_535, 5_504_980);
    assert_eq!(5_504_980 + 1, 5_504_981);
    assert!(255_u64 < 256_u64);
    assert!(255_u64.to_le_bytes() > 256_u64.to_le_bytes());
    let earlier_time_key = (255_u64, 4_u8, 256_u64, [1_u8; 32]);
    let later_time_key = (256_u64, 0_u8, 255_u64, [0_u8; 32]);
    assert!(earlier_time_key < later_time_key);
    let earlier_sequence_key = (7_u64, 4_u8, 255_u64, [1_u8; 32]);
    let later_sequence_key = (7_u64, 4_u8, 256_u64, [0_u8; 32]);
    assert!(earlier_sequence_key < later_sequence_key);
    for exact in [
        "terminal_trigger is the closed enum {CancellationSeal=1,SuccessfulWorkBoundary=2,BudgetOrTimeoutBoundary=3,InfrastructureFailureBoundary=4}",
        "TriggerCandidateV1=(effective_interval_lower_ns,causal_rank,trigger_logical_sequence,stable_event_id,trigger_kind,terminal_subcause)",
        "terminal_subcause is the closed enum {NotApplicable=0,Timeout=1,BudgetExhausted=2}",
        "24 exact ASCII bytes I13_TRIGGER_CANDIDATE_V1 followed by one byte 0x00",
        "exactly 76 bytes total",
        "Valid kind/rank/subcause triples are exactly {(InfrastructureFailureBoundary,0,NotApplicable),(CancellationSeal,1,NotApplicable),(BudgetOrTimeoutBoundary,2,Timeout),(BudgetOrTimeoutBoundary,3,BudgetExhausted),(SuccessfulWorkBoundary,4,NotApplicable)}",
        "TriggerCandidateSetRoot=BLAKE3::derive_key('org.frankensim.i13.trigger-candidate-set.v1'",
        "Every durably eligible cause and its canonical candidate bytes commit in one attempt/epoch-scoped atomic registration",
        "The frontier is AUTHENTICATED_SET_FRONTIER_V1 with domain 'org.frankensim.i13.trigger-candidate-membership'",
        "candidate_membership_root_n=AuthenticatedSetRootV1 over exactly n unique TriggerCandidateBytes",
        "generation=count=n",
        "A cap+1 InternalSchedulable operation is UnsupportedOperation before it is armed, proof allocation or effect",
        "34 exact ASCII bytes I13_UNAVOIDABLE_SOURCE_UNIVERSE_V1 followed by one byte 0x00",
        "SourceAuthorityRowV1=one byte phase where PreSelection=1 or PostSelection=2",
        "exactly 43 bytes",
        "source_count is in 1..=256",
        "max_canonical_facts is exactly 1",
        "complete universe is at most 11045 bytes",
        "UnavoidableSourceUniverseRoot=BLAKE3::derive_key('org.frankensim.i13.unavoidable-source-universe.v1'",
        "binds its root into candidate/witness frontier contexts, both capacity-breach records and all primary/emergency reservations",
        "256/257 boundary",
        "An UnavoidableExternal eligible fact at full capacity enters CANDIDATE_CAPACITY_BREACH_V1",
        "36 exact ASCII bytes I13_CANDIDATE_SOURCE_CLASS_MATRIX_V1 followed by one byte 0x00",
        "SuccessfulWorkBoundary=1/InternalSchedulable=1",
        "CancellationProposal=3/UnavoidableExternal=2",
        "SupervisorTimeoutBoundary=4/UnavoidableExternal=2",
        "ReceiptQualifiedInfrastructureFailure=5/UnavoidableExternal=2",
        "CandidateSourceClassMatrixRoot=BLAKE3::derive_key('org.frankensim.i13.candidate-source-class-matrix.v1'",
        "context U64LE(attempt_id)||U64LE(epoch)||raw32(CandidateSourceClassMatrixRoot)||raw32(UnavoidableSourceUniverseRoot)",
        "32 exact ASCII bytes I13_CANDIDATE_CAPACITY_BREACH_V1 followed by one byte 0x00",
        "full_generation=full_count=65536",
        "breach_reason where UnavoidableExternal=1 or InternalAdmissionInvariantViolation=2",
        "whole breach record is at most 1048352 bytes",
        "candidate_capacity_breach_digest=BLAKE3::derive_key('org.frankensim.i13.candidate-capacity-breach.v1'",
        "selection_boundary_logical_sequence=checked_add(max(predecessor snapshot maximum logical sequence,breach source-event logical sequence),1)",
        "terminal_selection_transaction_max_logical_sequence is the maximum sequence of every event that batch emits",
        "sets causal_rank=0, stable_event_id=candidate_capacity_breach_digest, trigger_kind=InfrastructureFailureBoundary and terminal_subcause=NotApplicable",
        "32 exact ASCII bytes I13_CAPACITY_BREACH_SELECTION_V1 followed by one byte 0x00",
        "exactly 501 bytes",
        "never claims TriggerCandidateSetRoot, RejectedCandidateRoot or typed-minimum arbitration",
        "A successfully persisted breach is never catastrophe and cannot remain indefinitely unfinalized",
        "Terminal-selection precomputation snapshots the exact membership root/generation/count/authenticated traversal and exact terminal-event root/generation/authenticated traversal in one fenced read",
        "selection_boundary_logical_sequence=checked_add(max(logical_sequence over all snapshotted terminal-event records, with genesis maximum 0),1)",
        "every snapshotted candidate source sequence is strictly less than that boundary",
        "there is no independent candidate-registry alias",
        "Any candidate/noncandidate event insert or Open-to-Breached transition changes one expectation and aborts a racing normal selection without effect",
        "Candidate membership and terminal-event authority are both order-independent sparse-Merkle sets",
        "generation/count-versus-logical-sequence separation",
        "candidate_count is in 1..=65536",
        "candidate-set projection size is exactly checked_add(16,checked_mul(candidate_count,84))",
        "MAX_TRIGGER_CANDIDATE_SET_PROJECTION_BYTES=5505040",
        "count 65536 reaches that exact byte cap",
        "candidate-set byte-cap-plus-one mutation is the 5505040-byte maximum valid candidate-set projection plus one trailing byte",
        "rejected_count=checked_sub(candidate_count,1)",
        "RejectedCandidateRoot=BLAKE3::derive_key('org.frankensim.i13.trigger-rejected-set.v1'",
        "Nonwinning candidate bytes remain unique and lexicographically sorted",
        "rejected-set projection size is exactly checked_add(40,checked_mul(rejected_count,84))",
        "MAX_REJECTED_CANDIDATE_SET_PROJECTION_BYTES=5504980",
        "rejected-set byte-cap-plus-one mutation is the 5504980-byte maximum valid rejected-set projection plus one trailing byte",
        "The Arbitrated selection receipt binds boundary identity, closure epoch, candidate-set root/count, exact winning bytes, RejectedCandidateRoot/count and terminal_selection_transaction_max_logical_sequence",
        "minimum of decoded typed arbitration key K=(effective_interval_lower_ns:u64,causal_rank:u8,trigger_logical_sequence:u64,stable_event_id:[u8;32])",
        "integer fields compare as unsigned numeric values and only stable_event_id compares lexicographically by raw bytes",
        "This is not lexicographic comparison of TriggerCandidateBytes because U64LE does not preserve numeric order",
        "encoded-byte ordering has authority only for canonical candidate-set and rejected-set roots",
        "two non-byte-identical candidates with the same K are IntegrityFailed",
        "InfrastructureFailure=0,CancellationSeal=1,Timeout=2,BudgetExhausted=3,SuccessfulWork=4",
        "Arrival, worker, thread, shard and candidate presentation order have no authority",
        "Compare-and-swap arrival order is not arbitration",
        "two independent proposal/candidate/set/source-class/breach encoders, published exact-byte/root KATs",
        "simultaneous-candidate, identical-time/rank, candidate-permutation, thread/shard-permutation, zero/count-cap/count-cap-plus-one",
        "candidate-projection exact-cap/cap-plus-one, rejected-projection exact-cap/cap-plus-one",
        "numeric-255-versus-256 U64LE counterexamples in both time and sequence positions",
        "proposal UTF-8/projection exact-cap/cap-plus-one",
        "same-K/nonidentical-candidate, omitted/extra-candidate, candidate-head-race",
        "proposal-wins/proposal-loses/racing-proposal, selection-versus-breach, every partial-transaction crash and exact replay/conflict test",
        "accepted cancellation ingress atomically appends one durable nonterminal cancellation.proposed event",
        "28 exact ASCII bytes I13_CANCELLATION_PROPOSAL_V1 followed by one byte 0x00",
        "complete projection has at most 65709 bytes under checked arithmetic",
        "proposal_event_id=BLAKE3::derive_key('org.frankensim.i13.cancellation-proposal.v1'",
        "cancellation.proposed neither seals the scope nor performs cancellation and is not cancellation.requested",
        "a CancellationSeal trigger means this proposal event/interval, never the later request event",
        "Only if that proposal's candidate wins does the one atomic terminal-selection transaction commit its selection receipt",
        "freeze the observer catalog and close the descendant/transition/access frontiers",
        "A losing proposal remains committed under RejectedCandidateRoot",
        "A racing ingress is either fully appended and registered before closure and therefore included, or refused without an event/effect",
        "Exact-key replay of the proposal key always returns the byte-identical original proposal receipt without another event",
        "a separate non-authoritative status envelope reports Pending or the bound immutable winning/losing selection result",
        "Exact-key replay of the distinct terminal-selection operation key returns the byte-identical all-or-none terminal-selection transaction",
        "For every non-cancellation winner the same transaction appends no cancellation.requested, cancellation.observed or execution.cancelled event",
        "an explicit ZeroWorkSuccessfulBoundary receipt when the admitted work set is empty",
        "precede the selected terminal trigger—for CancellationSeal this means the durable cancellation.proposed event/interval, never the later cancellation.requested event",
        "No post-trigger effect may be performed",
        "The fenced Arbitrated transaction atomically validates the typed minimum",
        "compare-and-swaps the expected terminal head, candidate-registration-frontier head/generation/count and derived TriggerCandidateSetRoot/count",
        "observer-catalog, descendant-frontier, transition-frontier and access-frontier heads/root/counts and exact reservation",
        "Any expected-head, root, count or reservation change aborts without an event or effect and forces recomputation",
        "no durable state may contain a CancellationSeal selection without its request/seal or a request/seal without that exact selection",
        "A racing registration either commits before terminal-selection closure with deterministic catalog/frontier membership and reservation receipts, or is refused without catalog/frontier membership or effect",
        "observer-catalog, descendant-frontier, transition-frontier and access-frontier roots/counts",
        "observer-observation, descendant-exit, transition-reconciliation and access-reconciliation roots/counts",
        "drain_latency=checked_sub(drained_interval_upper,trigger_interval_lower)",
        "finalize_latency=checked_sub(finalized_interval_upper,drained_interval_lower)",
        "terminal_latency=checked_add(drain_latency,finalize_latency)",
        "{2000000000,500000000,2500000000} Core",
        "{8000000000,2000000000,10000000000} Max",
        "Equality at each cap passes and cap+1 fails",
        "priority_observation_latency_ns=checked_add(checked_add(worst_case_request_to_next_priority_poll_ns,worst_case_observer_reaction_and_enqueue_ns),checked_add(worst_case_priority_queue_delay_ns,checked_add(checked_mul(priority_segment_count,worst_case_durable_service_time_ns),clock_calibration_uncertainty_ns)))",
        "must be <=250000000 Core or <=1000000000 Max",
        "Equality passes; cap+1, asymmetric capacity, arithmetic failure, borrowing, priority starvation, stale service receipt or non-independent emergency service refuses admission/capability issue",
        "The drained/finalized path is included in CANCELLATION_BOUNDS and TERMINAL_PERSISTENCE_LATENCY_V1",
        "A persisted sink failure adds an InfrastructureFailure/bit4 witness",
        "While neither lane has canonically finalized but a reserved prefix/replay route remains active, the attempt is TerminalPersistencePending",
        "Only failure of both mirrored lanes and every authenticated replay route permits TerminalPersistenceCatastrophe",
        "with no canonical finalization, disposition, completeness, adjudication or publication",
        "later evidence-adjudication derivation maps every admissible SLO-overrun witness to claim_adjudication=Unknown",
        "No overrun ever becomes scientific Failed or Refuted",
        "Exact resource-ceiling exhaustion adds a BudgetExhausted/bit2 witness",
        "receipt-qualified supervisor or sink failure adds an InfrastructureFailed/bit4 witness",
        "Witnesses are append-only facts; none directly assigns execution_disposition or claim_adjudication",
        "SuccessfulWorkBoundary->Completed/bit0, CancellationSeal->Cancelled/bit1",
        "BudgetOrTimeoutBoundary/BudgetExhausted->BudgetExhausted/bit2",
        "BudgetOrTimeoutBoundary/Timeout->TimedOut/bit3",
        "InfrastructureFailureBoundary->InfrastructureFailed/bit4",
        "25 exact ASCII bytes I13_TERMINAL_CONDITION_V1 followed by one byte 0x00",
        "ResourceExhaustion=2,SloOverrun=3,InfrastructureFailure=4",
        "exactly 67 bytes",
        "logical sequence is strictly greater than the receipt-bound terminal_selection_transaction_max_logical_sequence and strictly less than terminal_condition_cutoff_logical_sequence",
        "TerminalSelectionReceiptV1 record is a transitive ancestor",
        "Witness bytes are unique and lexicographically sorted",
        "witness_count is in 0..=65536",
        "exact size checked_add(40,checked_mul(witness_count,75))",
        "MAX_TERMINAL_CONDITION_WITNESS_PROJECTION_BYTES=4915240",
        "cardinality-zero root remains present",
        "two independent witness encoders must reproduce a published exact-byte/root KAT",
        "derived-byte-cap/trailing-byte-cap-plus-one, unknown-kind, zero/reused/mismatched-digest",
        "selection-max-equal/cutoff-equal/post-cutoff sequence",
        "The terminal-selection transaction opens one attempt/epoch witness frontier bound to selection_receipt_digest and terminal_selection_transaction_max_logical_sequence",
        "It is AUTHENTICATED_SET_FRONTIER_V1 with domain 'org.frankensim.i13.terminal-condition-witness-frontier'",
        "Every eligible post-selection condition source and its receipt insert in one bounded sparse-Merkle compare-and-swap",
        "The source and receipt are both causally after the complete selection transaction and before cutoff",
        "A source at or before the selection boundary belonged in TriggerCandidateSet and cannot be laundered later",
        "POST_SELECTION_CAPACITY_BREACH_V1",
        "41 exact ASCII bytes I13_POST_SELECTION_SOURCE_CLASS_MATRIX_V1 followed by one byte 0x00",
        "ResourceExhaustionBoundary=1/InternalSchedulable=1",
        "SupervisorSloOverrun=2/UnavoidableExternal=2",
        "37 exact ASCII bytes I13_POST_SELECTION_CAPACITY_BREACH_V1 followed by one byte 0x00",
        "raw32(PostSelectionSourceClassMatrixRoot), raw32(UnavoidableSourceUniverseRoot)",
        "frontier_kind where WitnessFrontier=1 or TerminalEventFrontier=2",
        "breach_reason where UnavoidableExternal=1 or InternalAdmissionInvariantViolation=2",
        "complete record is at most 1048322 bytes",
        "reason UnavoidableExternal requires one of the two external matrix rows plus exact PostSelection universe membership",
        "InternalAdmissionInvariantViolation requires the internal ResourceExhaustionBoundary row",
        "WitnessFrontier requires full_generation=full_count=65536",
        "TerminalEventFrontier requires full_generation=full_count=177407",
        "PostSelectionCapacityBreachUnionV1={Absent=0,Present=1||raw32(digest)}",
        "It does not pretend the overflow source is a sparse-set member or terminal-event parent",
        "never waits for an impossible cap+1 insert",
        "forces evidence_integrity=IntegrityFailed and InfrastructureFailed/bit4",
        "exactly one non-promoting FailureBundle",
        "Finalization precomputation atomically snapshots either the Open witness root/generation/count or the Breached root/generation/count plus exact Present breach digest",
        "expected_witness_frontier_count=witness_count",
        "terminal_condition_cutoff_logical_sequence=checked_add(max(logical_sequence over all snapshotted terminal-event records, with empty maximum 0),checked_add(barrier_level_count,1))",
        "Under one valid CommitBeforeGrantV1, finalization consumes the permit",
        "TerminalConditionCutoffV1=(terminal_condition_cutoff_logical_sequence,TerminalConditionWitnessSetRoot,witness_count,deadline_resolution_root,deadline_resolution_count)",
        "There is no cutoff without finalization or finalization without the exact barrier/cutoff/breach union",
        "A split, stale core, expired/reused permit, root/generation/count/state/all-and-only mismatch",
        "drain_latency=checked_sub(drained_interval_upper,trigger_interval_lower)",
        "finalize_deadline=checked_add(drained_interval_lower,finalize_cap)",
        "remaining_total=checked_sub(total_cap,drain_latency)",
        "total_finalize_deadline=checked_add(drained_interval_lower,remaining_total)",
        "DeadlineV1 is the closed union AbsoluteDeadline=1||U64LE(ns) or AlreadyExceededByDrain=2||U64LE(drain_latency)||U64LE(total_cap)",
        "Finalize uses AbsoluteDeadline, while Total uses AbsoluteDeadline when checked subtraction succeeds and AlreadyExceededByDrain when it underflows",
        "latter immediately resolves Total Exceeded with an idempotent pre-finalization SloOverrun receipt and is not an arithmetic-integrity failure",
        "31 exact ASCII bytes I13_TERMINAL_DEADLINE_TICKET_V1 plus NUL",
        "exactly 202 bytes for AbsoluteDeadline or 210 bytes for AlreadyExceededByDrain",
        "ticket_digest=BLAKE3::derive_key('org.frankensim.i13.terminal-deadline-ticket.v1'",
        "Ticket fields are derived, never caller-selected",
        "Recompute drain_latency from the bound selection/drained receipt intervals and select the frozen Core or Max caps committed for that attempt",
        "Finalize requires outer cap=the frozen finalize_cap and exactly AbsoluteDeadline(checked_add(drained_interval_lower,finalize_cap))",
        "Total requires outer cap=the frozen total_cap",
        "exactly AbsoluteDeadline(checked_add(drained_interval_lower,remaining_total))",
        "requires exactly AlreadyExceededByDrain(drain_latency,total_cap)",
        "encoded union arm and every duplicated cap, latency and deadline value must byte-equal that recomputation",
        "attempt/epoch mismatch, unauthorized monitor, non-atomic or non-causal arm sequence or arithmetic substitution is IntegrityFailed",
        "exact-byte/root KAT twins mutate each derived field and cross-field relation independently",
        "An unresolved absolute ticket becomes OnTime only through the pre-existing CommitBeforeGrantV1",
        "If the bound misses a deadline, no grant is issued for that core",
        "No finalization may commit until each ticket is either already Exceeded in FinalizationCoreV1",
        "35 exact ASCII bytes I13_TERMINAL_DEADLINE_RESOLUTION_V1 plus NUL",
        "one kind byte Finalize=1 or Total=2, one outcome byte OnTime=1 or Exceeded=2",
        "exactly 102 bytes",
        "deadline_resolution_count is exactly 2",
        "deadline_resolution_root=BLAKE3::derive_key('org.frankensim.i13.terminal-deadline-resolution-set.v1'",
        "exact 292-byte root projection",
        "neither may bind execution.finalized or its digest",
        "both deadline-union arms, equality/cap-plus-one, both tickets exceeded, watchdog/finalizer permutations, crash/replay",
        "unknown kind/outcome, swapped/missing/duplicate records, cross-ticket receipt, zero digest, count/root mismatch",
        "no finalized-before-resolution and no cutoff/finalized split",
        "post-trigger witness set contains every and only unique receipt-qualified",
        "PostConditionMaskV1 sets bit2 if and only if",
        "post-trigger evidence can never set Completed/bit0 or Cancelled/bit1",
        "TerminalConditionBitsV1 must byte-equal base_bit OR PostConditionMaskV1",
        "every ordinary witness must contribute its bit, a Present breach must contribute bit4 and IntegrityFailed",
        "TerminalConditionWitnessSetRoot=BLAKE3::derive_key('org.frankensim.i13.terminal-condition-witness-set.v1'",
        "InfrastructureFailed>TimedOut>BudgetExhausted>Cancelled>Completed",
        "independently of witness registration, append, worker, thread, shard or serialization order",
        "there is no mutable execution-disposition setter",
        "this post-trigger join order is distinct from and never changes DRAIN_TRIGGER_V1 causal ranks",
        "The selected trigger remains separately immutable and is never rewritten",
        "execution.finalized binds the exact condition byte, TerminalConditionCutoffV1, witness-set root/count, breach union, deadline-resolution root/count and derived disposition/integrity",
        "An omitted/extra/duplicate/inadmissible witness, missing/forged breach, missing or unearned bit, unknown bit",
        "one-byte final_disposition tag is exactly Completed=0,Cancelled=1,BudgetExhausted=2,TimedOut=3,InfrastructureFailed=4",
        "byte-equals max{i in 0..=4 | TerminalConditionBitsV1 AND (1<<i) is nonzero}",
        "all 31 syntactically bounded nonzero five-bit inputs are exercised for the pure highest-set-bit decoder",
        "lifecycle acceptance separately enforces exactly one base bit plus only earned post bits",
        "exactly the derived tag accepted",
        "witness omission/addition/duplication, breach add/drop/swap, simultaneous additions and every witness/breach/condition arrival permutation",
        "simultaneous conditions contribute both and only FINAL_EXECUTION_DISPOSITION_JOIN_V1 derives the final disposition",
    ] {
        assert!(
            policy.contains(exact),
            "exact latency contract lost '{exact}'"
        );
    }
}

#[test]
fn i13_finalization_authority_is_acyclic_bounded_and_recoverable() {
    let policy = policy_spec();

    let finalization_core_header = b"I13_FINALIZATION_CORE_V1\0";
    let commit_before_grant_header = b"I13_COMMIT_BEFORE_GRANT_V1\0";
    let finalization_commit_receipt_header = b"I13_FINALIZATION_COMMIT_RECEIPT_V1\0";
    let emergency_terminal_slot_header = b"I13_EMERGENCY_TERMINAL_SLOT_V1\0";
    let emergency_journal_commit_header = b"I13_EMERGENCY_JOURNAL_COMMIT_V1\0";
    let ontime_receipt_header = b"I13_ONTIME_RECEIPT_V1\0";
    let finalization_envelope_header = b"I13_FINALIZATION_ENVELOPE_V1\0";
    let terminal_barrier_header = b"I13_TERMINAL_CAUSAL_BARRIER_V1\0";
    let lane_epoch_transition_header = b"I13_LANE_EPOCH_TRANSITION_V1\0";
    let lane_epoch_receipt_header = b"I13_LANE_EPOCH_TRANSITION_RECEIPT_V1\0";

    assert_eq!(finalization_core_header.len(), 25);
    assert_eq!(commit_before_grant_header.len(), 27);
    assert_eq!(finalization_commit_receipt_header.len(), 35);
    assert_eq!(emergency_terminal_slot_header.len(), 31);
    assert_eq!(emergency_journal_commit_header.len(), 32);
    assert_eq!(ontime_receipt_header.len(), 22);
    assert_eq!(finalization_envelope_header.len(), 29);
    assert_eq!(terminal_barrier_header.len(), 31);
    assert_eq!(lane_epoch_transition_header.len(), 29);
    assert_eq!(lane_epoch_receipt_header.len(), 37);
    assert_eq!(31 + 2 * 8 + 32 + 5 * 8 + 32, 151);
    assert_eq!(151 + 64 * 32, 2_199);
    assert_eq!(27 + 10 * 32 + 11 * 8 + 1, 436);
    assert_eq!(35 + 1 + 3 * 8 + 9 * 32 + 2 * 8 + 8 + 2 * 8, 388);
    assert_eq!(22 + 3 * 32 + 3 * 8, 142);
    assert_eq!(31 + 1 + 3 * 8 + 9 * 32 + 3 * 8, 368);
    assert_eq!(368 + 1_048_576, 1_048_944);
    assert_eq!(
        32 + 3 * 8 + 32 + 8 + 32 + 32 + 8 + 32 + 8 + 2 * 8 + 32 + 32,
        288
    );
    assert_eq!(1_048_944 + 288 + 388, 1_049_620);
    assert_eq!(4_194_304 - 1_049_620, 3_144_684);
    assert_eq!(29 + 2 * 8 + 32 + 2 * 8 + 2 + 32 + 32 + 8 + 4 * 32, 295);
    assert_eq!(37 + 5 * 32 + 4 * 8 + 1, 230);
    assert_eq!(256 * 32 + 2 * 8 + 8 + 76, 8_292);
    assert_eq!(256 * 32 + 2 * 8 + 8 + 67, 8_283);
    assert_eq!(8 + 8 + 32 + (8 + 8 + 256) + 32, 352);
    assert_eq!(256 * 32 + 2 * 8 + 8 + 352, 8_568);
    assert_eq!(8 + 64 * 32, 2_056);
    assert_eq!(8 + 16_384 * 66, 1_081_352);
    assert_eq!(180_224 * 512, 92_274_688);
    let artifact_bytes = |n: u64| n * 354 + (n - 1) * 66 + n * 16 + n.div_ceil(1_024) * 120;
    assert_eq!(artifact_bytes(1), 490);
    assert_eq!(artifact_bytes(1_024), 446_518);
    assert_eq!(artifact_bytes(1_025), 447_074);
    assert_eq!(artifact_bytes(180_224), 78_598_718);
    for n in [1, 2, 1_024, 1_025, 180_224] {
        assert!(artifact_bytes(n) <= n * 512);
    }
    assert_eq!(67_108_864 + 2 * 92_274_688 + 2 * 4_194_304, 260_046_848);
    assert_eq!(268_435_456 - 260_046_848, 8_388_608);
    assert_eq!(763 + 1 + 1_047_812, 1_048_576);
    assert_eq!(763 + 33 + 1_047_780, 1_048_576);
    assert_eq!(763 + 1 + 64 * 32 + 1_045_764, 1_048_576);
    let barrier_shape = |n: u64| {
        let mut width = n;
        let mut levels = 0_u64;
        let mut events = 0_u64;
        while width > 64 {
            width = width.div_ceil(64);
            levels += 1;
            events += width;
        }
        (levels, events, width, n + events + 1)
    };
    assert_eq!(barrier_shape(0), (0, 0, 0, 1));
    assert_eq!(barrier_shape(64), (0, 0, 64, 65));
    assert_eq!(barrier_shape(65), (1, 2, 2, 68));
    assert_eq!(barrier_shape(4_096), (1, 64, 64, 4_161));
    assert_eq!(barrier_shape(4_097), (2, 67, 2, 4_165));
    assert_eq!(barrier_shape(177_407), (2, 2_816, 44, 180_224));
    assert_eq!(barrier_shape(177_408), (2, 2_816, 44, 180_225));

    for exact in [
        "candidate_membership_root_n=AuthenticatedSetRootV1 over exactly n unique TriggerCandidateBytes",
        "AUTHENTICATED_SET_FRONTIER_V1",
        "SparseSetUpdateProofBytesV1=U64LE(old_count)||U64LE(new_count)||FRAME_BYTES(record)",
        "bit 7-(d mod 8) of record_key byte floor(d/8)",
        "exact maxima are 8292 bytes for a 76-byte candidate, 8283 for a 67-byte witness and 8568 for a terminal-event record",
        "A valid insertion proof reconstructs old leaf=empty_256 and new leaf=leaf for the same derived record_key",
        "deletion, replacement and occupied-leaf insertion are forbidden",
        "insertion never sorts or rehashes the existing set",
        "authenticated hash-key-order leaf traversal proves all-and-only leaves/count; it does not claim lexical record order",
        "bounded stable external radix sort then orders the exact full record bytes",
        "its exact genesis uses empty_0/count 0, and generation=count=n",
        "independently verified insertion root/n+1",
        "Exact-key replay returns the byte-identical receipt without insertion",
        "duplicate candidate, stale root/generation/count, malformed proof, overflow, cross-attempt/epoch/matrix/universe reuse or ABA is IntegrityFailed",
        "reversed multi-insert permutation, noncandidate-interleaving, exact-cap, pre-arm internal refusal, authenticated internal-invariant breach, unavoidable cap-plus-one, selection-versus-breach race and replay KATs",
        "both tickets must arm atomically at one sequence strictly after execution.drained",
        "Both tickets' attempt/epoch must byte-equal the attempt/epoch bound by the selection and drained receipts",
        "monitor identity must be authorized by the frozen attempt profile and capability epoch",
        "equal/reversed arm sequence, cross-attempt/epoch and unauthorized-monitor cases",
        "one non-borrowable FinalizationServiceReservationV1 from the primary priority sink and one from the independently serviced emergency terminal journal",
        "one exclusive 4194304-byte terminal segment",
        "24 exact ASCII bytes I13_FINALIZATION_CORE_V1 plus NUL",
        "The template has exactly two kind-ordered entries",
        "PendingOnTime=0 with no payload or Exceeded=1||raw32(ticket_digest)||raw32(SloOverrun_receipt_digest)",
        "The core contains no grant/permit digest, OnTime receipt, deadline-resolution root, execution.finalized digest or post-commit receipt",
        "CommitBeforeGrantV1 is an authenticated, single-use forward FinalizationServiceCertificateV1 issued before execution.finalized serialization",
        "exact 436-byte projection",
        "one byte lane||U64LE(lane_epoch)",
        "lane/lane_epoch must be exactly Primary/0 or Emergency/1",
        "certified_persistence_upper=checked_add(finalization_ready_interval_upper,checked_add(max_prepare_and_submit_ns,checked_add(max_queue_delay_ns,checked_add(max_durable_service_ns,clock_uncertainty_ns))))",
        "All inputs are pre-existing receipt-bound values, never an observation or prediction derived from this record's own eventual durability",
        "expiration_ns must be at least certified_persistence_upper",
        "A grant may leave an entry PendingOnTime only if certified_persistence_upper is no later than that ticket's inclusive absolute deadline",
        "OnTimeResolutionReceiptV1 is independently derived from the grant and ticket before the final envelope",
        "The grant binds no OnTime receipt, resolution root, final-envelope/execution.finalized digest or post-commit receipt",
        "either makes every derived barrier record, exact core-derived final envelope and terminal-root compare-and-swap durable by certified_persistence_upper or has zero canonical finalization effect",
        "Head/core/ticket/barrier/reservation/lane-epoch change, expiry, duplicate consume or primary-to-emergency replay aborts and invalidates the grant",
        "34 exact ASCII bytes I13_FINALIZATION_COMMIT_RECEIPT_V1 plus NUL",
        "exactly 388 bytes",
        "successor_terminal_generation=checked_add(predecessor_terminal_generation,checked_add(barrier_event_count,1))",
        "successor_terminal_head_digest is TERMINAL_EVENT_CHAIN_V1's independently recomputed order-independent root",
        "This receipt is downstream of and never embedded in the final envelope",
        "actual_durable_interval_upper<=certified_persistence_upper",
        "reservation -> FinalizationCoreV1 -> CommitBeforeGrantV1 -> OnTime receipts/resolution root/final envelope -> barrier/final terminal-event-set insertion -> FinalizationCommitReceiptV1 -> adjudication",
        "An Exceeded resolution sequence is strictly after ticket arming and strictly before the cutoff",
        "OnTime resolution sequence is the same atomic logical sequence as its grant-consuming cutoff/finalization",
        "equal/reversed arm and resolution sequences, grant/final-envelope cycles",
        "exactly one non-borrowable 4194304-byte EmergencyTerminalSegmentV1 for each attempt",
        "At most one EmergencyTerminalSlotV1 and one correlated EmergencyJournalCommitReceiptV1 may commit for one attempt/epoch",
        "30 exact ASCII bytes I13_EMERGENCY_TERMINAL_SLOT_V1 plus NUL",
        "The fixed overhead is exactly 368 bytes",
        "terminal_transaction_byte_count=exact_finalized_record_byte_count",
        "complete slot size is exactly checked_add(368,terminal_transaction_byte_count) and at most 1048944 bytes",
        "slot_digest=BLAKE3::derive_key('org.frankensim.i13.emergency-terminal-slot.v1'",
        "emergency_journal_head_0=BLAKE3::derive_key('org.frankensim.i13.emergency-journal-head.v1'",
        "successor_generation=checked_add(n,1)",
        "31 exact ASCII bytes I13_EMERGENCY_JOURNAL_COMMIT_V1 plus NUL",
        "exactly 288 bytes",
        "maximum 1048944-byte slot plus the 288-byte journal and 388-byte finalization receipts is 1049620 bytes, leaving 3144684 bytes",
        "run_reserve_i=checked_add(67108864,checked_add(checked_mul(terminal_event_artifact_reserve_i,2),checked_mul(checked_mul(p_i,2),4194304)))",
        "TERMINAL_EVENT_ARTIFACT_V1",
        "MAX_TERMINAL_EVENT_ARTIFACT_BYTES_V1=92274688",
        "admission_max_terminal_event_count_i<=180224",
        "pre-reserved at exactly checked_mul(admission_max_terminal_event_count_i,512) bytes",
        "TerminalEventArtifactLeafEntryV1=U16LE(record_byte_count)||record is at most 354 bytes",
        "NodeEntryV1=U16LE(depth)||raw32(left_child)||raw32(right_child) is exactly 66 bytes",
        "RecoveryIndexEntryV1=U64LE(leaf_ordinal)||U64LE(chunk_offset) is exactly 16 bytes per leaf",
        "Each nonempty canonical chunk adds exactly one 120-byte header",
        "Define checked_ceil_div(n,1024)=checked_add(n,1023)/1024",
        "<= checked_mul(n,512) for every n>=1",
        "JSONL rows may index chunks, while every typed event remains a distinct authenticated leaf and payload",
        "At both maxima with p_i=1, one run reserves exactly 260046848 bytes and leaves 8388608 bytes governor slack",
        "global_reserved=checked_sum_i(run_reserve_i) must be <=268435456",
        "event-count/artifact/segment/governor cap+1, multiplication/addition/sum overflow, asymmetric/missing artifacts or segments",
        "primary_finalization_path_ns=checked_add(primary_prepare_and_submit_ns,checked_add(primary_queue_delay_ns,primary_terminal_segment_durable_service_ns))",
        "emergency_finalization_path_ns=checked_add(primary_prepare_and_submit_ns",
        "certified_terminal_persistence_bound_ns=checked_add(max(primary_finalization_path_ns,emergency_finalization_path_ns),clock_calibration_uncertainty_ns)",
        "elapsed primary prepare/queue/service/failure/fence time remains charged by the absolute tickets and admission equation, never restarted for free",
        "A receipt-qualified primary failure before finalization proves the primary grant had no canonical finalization effect",
        "forces a fresh FinalizationCoreV1 and Emergency-lane CommitBeforeGrantV1",
        "Before the emergency all-or-none terminal transaction commits, the attempt remains TerminalPersistencePending",
        "after its slot, both heads and execution.finalized become durable together",
        "no canonical reconstruction/recommit phase exists",
        "Later primary copyback is non-authoritative replication",
        "LaneEpochHeadV1=raw32(lane_epoch_head_digest)||U64LE(lane_epoch)",
        "lane_epoch_head_digest_0=BLAKE3::derive_key('org.frankensim.i13.lane-epoch-head.v1'",
        "only admitted pre-settlement transition is Primary/0 -> Emergency/1",
        "28 exact ASCII bytes I13_LANE_EPOCH_TRANSITION_V1 plus NUL",
        "exactly 295 bytes",
        "36 exact ASCII bytes I13_LANE_EPOCH_TRANSITION_RECEIPT_V1 plus NUL",
        "exactly 230 bytes",
        "successor head and this receipt are the authority pair",
        "rehashing transition intent alone never proves a committed failover",
        "recovered primary can never regain canonical append authority for that attempt",
        "TERMINAL_EVENT_CHAIN_V1: The lane-agnostic logical terminal authority is an AUTHENTICATED_SET_FRONTIER_V1 event set",
        "terminal generation=count=n only—it is never a logical clock",
        "CausalParentSetRoot=BLAKE3::derive_key('org.frankensim.i13.terminal-causal-parent-set.v1'",
        "TerminalCausalSchemaBytesV1 starts with the 29 exact ASCII bytes I13_TERMINAL_CAUSAL_SCHEMA_V1",
        "ParentExtractorAstV1 is the closed prefix grammar",
        "TerminalCausalSchemaRoot=BLAKE3::derive_key('org.frankensim.i13.terminal-causal-schema.v1'",
        "the AST root is depth 0, each framed child increments depth exactly once, depth 64 is admitted and 65 refused",
        "Admission statically evaluates every allowed extractor against the frozen registration maxima",
        "unless every event has at most 64 parents and its proof/work reservation fits",
        "Parent count/root must byte-equal the complete extractor result",
        "CanonicalParentMembershipMultiProofBytesV1=U64LE(entry_count)",
        "entry_count<=16384",
        "checked_add(8,checked_mul(entry_count,66))<=1081352 bytes",
        "Each entry depth is exactly 0..=255",
        "depth 0 requires an all-zero prefix",
        "encoded depth 256 is invalid",
        "duplicate/nonminimal entries, inconsistent shared nodes and nonzero masked suffix bits are refused",
        "self, future, duplicate, omitted, extra and cross-attempt parents are forbidden",
        "Each event logical_sequence=checked_add(max(parent logical sequences, with empty-parent maximum 0),1)",
        "rejects cycles, sorts the currently ready zero-indegree records by full bytes",
        "omit/add/swap/self/future/cycle parent",
        "incomparable events may share a sequence and stable full record bytes break every canonical sort tie",
        "ONTIME_RESOLUTION_RECEIPT_V1",
        "exactly 142 bytes",
        "TERMINAL_CAUSAL_BARRIER_V1",
        "30 exact ASCII bytes I13_TERMINAL_CAUSAL_BARRIER_V1 plus NUL",
        "exactly checked_add(151,checked_mul(parent_count,32)) bytes in 151..=2199",
        "typed parent_record_key list is the field selected by the terminal.causal_barrier ParentExtractorAstV1",
        "These bytes are the complete canonical_event_bytes",
        "While count(level_l)>64, partition level_l into consecutive chunks of at most 64 keys",
        "barrier_level_count is the number of generated levels",
        "barrier_idempotency_key_digest=BLAKE3::derive_key('org.frankensim.i13.terminal-barrier-idempotency.v1',TerminalCausalBarrierBytesV1)",
        "is the only permitted key for that record; no caller/lane key exists",
        "TerminalBarrierPlanDigest=BLAKE3::derive_key('org.frankensim.i13.terminal-barrier-plan.v1'",
        "final_parent_keys, which are embedded byte-for-byte in FinalizationEnvelopeBytesV1",
        "For snapshot_count=0..=64, level/event counts are zero",
        "for snapshot_count=65 they are 1/2",
        "for the maximum admitted nonbarrier count 177407 they are 2/2816",
        "checked_add(177407,checked_add(2816,1))=180224 exactly",
        "snapshot_count+barrier_event_count+1<=admission_max_terminal_event_count_i<=180224",
        "every barrier level plus execution.finalized inserts atomically in deterministic topological order",
        "maximum 2816 barrier records plus final record, their proofs, artifact writes and root CAS are explicitly included",
        "prepare-and-submit includes core/grant/barrier-plan validation, every sparse proof",
        "durable service includes every planned barrier/artifact insert plus the final envelope/root CAS",
        "FINALIZATION_ENVELOPE_V1",
        "U64LE(final_parent_count)||raw32(final_parent_set_root)||concat(raw32(final_parent_key) in lexical order)",
        "final_parent_count is in 0..=64",
        "PostSelectionCapacityBreachUnionV1 is the exact core-bound Absent one-byte or Present 33-byte union",
        "The fixed overhead excluding parent keys and the union payload is exactly 763 bytes",
        "breach_union_byte_count is exactly 1 or 33",
        "exact_finalized_record_byte_count=checked_add(checked_add(checked_add(763,breach_union_byte_count),checked_mul(final_parent_count,32)),immutable_evidence_payload_byte_count)",
        "DOWNSTREAM_RECEIPT_RECOVERY_V1",
    ] {
        assert!(
            policy.contains(exact),
            "acyclic finalization contract lost '{exact}'"
        );
    }
}

#[test]
fn i13_industrial_leaf_pins_typed_discharge_and_cancellation_safe_authority() {
    let draft = i13_draft();
    let industrial = draft
        .obligations
        .iter()
        .find(|row| row.leaf == "i13-industrial-validation-max")
        .expect("industrial validation leaf");
    let claim = draft
        .claims
        .iter()
        .find(|claim| claim.id == "i13-governed-industrial-machine-validation")
        .expect("industrial validation claim");
    let hypotheses = claim.hypotheses.join("\n");
    for required in [
        "predecessor manifest and transaction intent",
        "verified AmendmentRecord separately binds the final successor",
        "atomic transaction advances the authority head",
    ] {
        assert!(
            hypotheses.contains(required),
            "industrial claim hypothesis lost anti-cycle boundary '{required}'"
        );
    }
    for required in [
        "separately admitted i13-industrial-validation-max/governance attempt",
        "carrying no science or promotion authority while the waiver is live",
        "has completed GovernanceCommitted and CandidateFrozen",
        "one atomic RealizationCommitted fs-vvreg authority transaction",
        "commits before any candidate-side or adjudication-side protected access or reveal",
        "pre-discharge custodian/reviewer access remains restricted to the frozen GovernanceCommitted protocol",
        "After that governance attempt drains and finalizes with reconciled transition/access frontiers",
        "distinct i13-industrial-validation-max/science attempt may be admitted",
        "one-shot RevealedForAdjudication",
        "independent adjudication then precede Closed",
    ] {
        assert!(
            claim.activation.contains(required),
            "industrial claim activation lost lifecycle boundary '{required}'"
        );
    }
    assert!(!claim.activation.contains("closes before protected access"));
    for event in REQUIRED_GOVERNANCE_EVENTS {
        assert_eq!(
            industrial
                .obs_events
                .iter()
                .filter(|observed| *observed == event)
                .count(),
            1,
            "industrial leaf must declare {event} exactly once"
        );
    }

    for required in [
        "arbitrary-digest",
        "same-ID-fixture-only",
        "receipt-only",
        "waiver-only",
        "split-transaction",
        "exact raw/normalized roots",
        "censoring/missingness/dependence",
        "frozen expected QoIs/bands/simultaneous decision/adjudication code",
        "typed discharge envelope and receipt",
        "atomic amendment/authority-head transition",
        "generic manifest syntax or partial discharge grants no authority",
    ] {
        assert!(
            industrial.g0.contains(required),
            "industrial G0 lost discharge invariant '{required}'"
        );
    }

    for required in [
        "arbitrary same-ID External digest",
        "removing its waiver",
        "typed envelope",
        "verified receipt",
        "one atomic authority-head transaction",
        "cannot preserve discharge authority",
        "governance.transition or governance.protected_access event",
        "access-event-chain root is IntegrityFailed",
    ] {
        assert!(
            industrial
                .g3_relations
                .iter()
                .any(|relation| relation.contains(required)),
            "industrial G3 lost negative discharge relation '{required}'"
        );
    }

    for required in [
        "request-drain-finalize independently for the governance and science attempts",
        "cancel governance during discharge-envelope/receipt verification",
        "cancel science during one-shot capability issue/consume",
        "atomic amendment/authority-head commit",
        "Governance drain publishes one complete authenticated RealizationCommitted transaction",
        "durable commit point preceded the selected cancellation.proposed trigger",
        "science drain publishes one complete authenticated industrial-validation evidence transaction or none",
        "After the selected terminal trigger (the cancellation.proposed event for CancellationSeal) no new capability issue/use, protected read, reveal, adjudication, publication or authority-head advance may begin",
        "every registered transition reconciles to pre-proposal Committed or no-head-advance Refused/IntegrityFailed",
        "every in-flight Prepared access reconciles to exactly one terminal row without a post-proposal effect",
        "Deny capability issue and protected access unless the write-ahead Prepared",
        "cancellation or sink failure grants no new protected read, reveal, adjudication or authority-head advance",
        "governance cancellation before durable RealizationCommitted leaves science unadmitted",
        "pre-proposal durable commit remains but never auto-starts science",
        "Science cancellation cannot reopen governance, restore the waiver or mint another reveal attempt",
        "Resumed/forked protected-root, authority/access-head and adjudication state equals uninterrupted execution",
    ] {
        assert!(
            industrial.g4_schedule.contains(required),
            "industrial G4 lost cancellation-safe authority rule '{required}'"
        );
    }
}

#[test]
fn i13_slash_delimited_contract_tokens_survive_rust_line_continuations_exactly() {
    let draft = i13_draft();
    let mut contract = draft
        .fixtures
        .iter()
        .map(authored_spec)
        .collect::<Vec<_>>()
        .join("\n");
    for row in &draft.obligations {
        contract.push('\n');
        contract.push_str(row.g0);
        contract.push('\n');
        contract.push_str(row.g4_schedule);
    }
    for required in [
        "phase/branch/parallel-path/tap/terminal",
        "linkage/MMF/loss/end-turn",
        "domain/eigenvector/derivative",
        "residual/filter/source/QoI",
        "custody/link/domain/calibration/discharge/defect",
        "class/netlist/boundary/period",
        "grazing/simultaneous/cusp/hysteretic",
        "finite section/capacity/bend/separation",
        "global cover/bound/feasibility",
        "admission/implementation/proof-kernel-or-TCB",
        "FailureBundle/countermodel",
    ] {
        assert!(
            contract.contains(required),
            "line continuation corrupted exact token '{required}'"
        );
    }
}

#[test]
fn i13_theorem_ratchets_preserve_bold_targets_without_false_authority() {
    let draft = i13_draft();
    let policy = authored_spec(fixture(&draft.fixtures, CAMPAIGN_POLICY_FIXTURE));
    let diffchar = authored_spec(fixture(
        &draft.fixtures,
        "i13-differential-character-theorem-card",
    ));
    let route = authored_spec(fixture(&draft.fixtures, "i13-topology-route-theorem-card"));
    let grammar = authored_spec(fixture(&draft.fixtures, "i13-maximal-adversary-grammar"));

    for required in [
        "complete typed proposition/grammar AST",
        "all definitions and transitive dependencies",
        "runtime-premise map",
        "deterministic total AST-to-Lean/AST-to-enumerator translation",
        "sorryAx",
        "complete transitive closure {propext,Quot.sound,Classical.choice}",
    ] {
        assert!(
            policy.contains(required),
            "formal projection lost '{required}'"
        );
    }
    assert!(diffchar.contains("ordinary macroscopic real lane is R"));
    assert!(diffchar.contains("R/Lambda lane exists only"));
    assert!(diffchar.contains("physically specified charge/flux lattice Lambda"));
    assert!(diffchar.contains("force, loss and performance require separate"));
    assert!(diffchar.contains("nontrivial physical R3 winding example"));
    assert!(diffchar.contains("tagged abstract torsion example"));

    assert!(route.contains("topological class feasibility"));
    assert!(route.contains("graph/capacity embedding"));
    assert!(route.contains("knot/link/isotopy/geometric realization"));
    assert!(route.contains("process realization"));
    assert!(route.contains("constructive witness"));
    assert!(route.contains("arbitrary 3-D domains and every integral class are excluded"));

    assert!(grammar.contains("16 topology instances x 16 design/material instances x 16 event/route instances x 16 theorem-boundary instances = 65536 raw tuples"));
    assert!(grammar.contains("This prose target has no exhaustive authority"));
    assert!(
        grammar.contains("exactly 4096 trajectories x 4096 evaluations = 16777216 evaluations")
    );
    assert!(grammar.contains("explicitly non-exhaustive"));
    assert!(grammar.contains("Search survival never proves a theorem"));

    let continuation = draft
        .claims
        .iter()
        .find(|claim| claim.id == "i13-certified-topology-change-continuation")
        .expect("continuation claim");
    let continuation_text = format!(
        "{}\n{}\n{}",
        continuation.statement,
        continuation.hypotheses.join("\n"),
        continuation.no_claim
    );
    for required in [
        "classical saltation",
        "transverse",
        "Bouligand",
        "Clarke",
        "symmetric finite difference",
        "universal canonical continuation",
    ] {
        assert!(continuation_text.contains(required));
    }
}

#[test]
fn i13_core_route_and_event_authority_does_not_smuggle_in_maximal_synthesis() {
    let draft = i13_draft();
    let route = draft
        .claims
        .iter()
        .find(|claim| claim.id == "i13-route-netlist-realization-audit")
        .expect("core route claim");
    let event = draft
        .claims
        .iter()
        .find(|claim| claim.id == "i13-topology-event-lineage-invalidation")
        .expect("core event claim");
    let synthesis = draft
        .claims
        .iter()
        .find(|claim| claim.id == "i13-integral-winding-class-synthesis")
        .expect("maximal synthesis claim");
    let continuation = draft
        .claims
        .iter()
        .find(|claim| claim.id == "i13-certified-topology-change-continuation")
        .expect("maximal continuation claim");

    assert_eq!(route.ambition, Ambition::Solid);
    assert!(route.statement.contains("finite-cross-section coil route"));
    assert!(route.no_claim.contains("does not synthesize"));
    assert!(!route.statement.contains("necessary-and-sufficient"));
    assert_eq!(synthesis.ambition, Ambition::Moonshot);
    assert!(synthesis.statement.contains("mixed-integer synthesizer"));
    assert!(
        synthesis
            .no_claim
            .contains("not a connected finite-radius coil route")
    );

    assert_eq!(event.ambition, Ambition::Solid);
    assert!(event.statement.contains("deterministically invalidated"));
    assert!(
        event
            .no_claim
            .contains("do not prove a canonical class map")
    );
    assert_eq!(continuation.ambition, Ambition::Moonshot);
    assert!(
        continuation
            .statement
            .contains("conservative physical transfer")
    );
}

#[test]
fn i13_metamorphic_laws_preserve_semantic_order_and_require_proved_monotonicity() {
    let frozen = i13_draft().freeze().expect("freeze");
    let row = |leaf: &str| {
        frozen
            .obligations()
            .iter()
            .find(|row| row.leaf() == leaf)
            .unwrap_or_else(|| panic!("missing leaf '{leaf}'"))
    };

    let adjoint = row("i13-adjoint-optimizer-core");
    assert!(
        adjoint
            .g3_relations()
            .iter()
            .any(|relation| relation.contains("separately proved subset enclosure"))
    );
    assert!(
        adjoint
            .g3_relations()
            .iter()
            .any(|relation| relation.contains("merely smaller solver tolerance"))
    );

    let robust = row("i13-robust-qoi-core");
    assert!(robust.g0().contains("immutable logical observation order"));
    assert!(
        robust
            .g3_relations()
            .iter()
            .any(|relation| relation.contains("worker assignment or result-arrival order"))
    );
    assert!(
        robust
            .g3_relations()
            .iter()
            .any(|relation| relation.contains("explicitly exchangeable law"))
    );
    assert!(
        robust
            .g5_matrix()
            .contains("immutable logical observation order")
    );

    let industrial = row("i13-industrial-validation-max");
    assert!(
        industrial
            .g3_relations()
            .iter()
            .any(|relation| relation.contains("immutable logical test/time ids"))
    );
    assert!(
        industrial
            .g5_matrix()
            .contains("immutable logical test/time order")
    );

    let global = row("i13-global-robust-max-leaf");
    let global_claim = frozen
        .claim("i13-global-robust-manufacturable-optimum")
        .expect("global robust claim");
    let global_hypotheses = global_claim.hypotheses.join("\n");
    for required in [
        "ObjectiveSense={Minimize,Maximize}",
        "exact maximize-to-canonical-minimize sign transform",
        "bound orientation",
        "incumbent comparison",
        "globality gap formula",
    ] {
        assert!(
            global_hypotheses.contains(required),
            "global claim lost objective-direction contract '{required}'"
        );
    }
    let global_fixture = authored_spec(fixture(frozen.fixtures(), "i13-global-robust-max"));
    for required in [
        "OBJECTIVE CARD: freeze ObjectiveSense={Minimize,Maximize}",
        "Minimize mapping f_canonical=f_original",
        "Maximize mapping f_canonical=-f_original",
        "Lower bounds, feasible-incumbent upper bounds, original-sense comparison",
        "absolute_gap=upper_canonical-lower_canonical",
        "objective-sense/sign swap",
        "reversed bound",
    ] {
        assert!(
            global_fixture.contains(required),
            "global fixture lost objective-direction clause '{required}'"
        );
    }
    for required in [
        "objective-sense twins",
        "exact canonical minimize transform",
        "bound orientation",
        "objective-gap arithmetic",
    ] {
        assert!(
            global.g0().contains(required),
            "global G0 lost '{required}'"
        );
    }
    assert!(
        global
            .g3_relations()
            .iter()
            .any(|relation| relation.contains("intersected with the inherited parent lower bound"))
    );
    assert!(
        global
            .g3_relations()
            .iter()
            .any(|relation| relation.contains("cannot weaken"))
    );
    assert!(global.g3_relations().iter().any(|relation| {
        relation.contains("after the exact maximize-to-minimize transform it cannot decrease")
    }));

    let integer = row("i13-integer-synthesis-max-leaf");
    assert!(
        !integer.decks().contains(&"i13-topology-route-theorem-card"),
        "route-theorem revisions must not revoke route-independent integral synthesis"
    );
}

#[test]
fn i13_identity_is_stable_order_invariant_and_chunk_assembly_equivalent() {
    let one_shot = i13_draft();
    let frozen = one_shot.clone().freeze().expect("freeze");
    assert_eq!(
        frozen.digest(),
        i13_draft().freeze().expect("refreeze").digest()
    );

    let mut reordered = i13_draft();
    reordered.claims.reverse();
    reordered.fixtures.reverse();
    reordered.obligations.reverse();
    reordered.waivers.reverse();
    assert_eq!(
        reordered.freeze().expect("reordered freeze").digest(),
        frozen.digest()
    );

    // G4 precursor only: this clones partial in-memory state at deterministic
    // boundaries. It does not claim process restart, persistent checkpoint,
    // cancellation, corruption detection, or request-drain-finalize evidence.
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
        staged = staged.clone();
    }
    for chunk in one_shot.fixtures.chunks(5) {
        staged.fixtures.extend_from_slice(chunk);
        staged = staged.clone();
    }
    for chunk in one_shot.obligations.chunks(2) {
        staged.obligations.extend_from_slice(chunk);
        staged = staged.clone();
    }
    staged.waivers.extend_from_slice(&one_shot.waivers);
    assert_eq!(
        staged.freeze().expect("chunked freeze").digest(),
        frozen.digest()
    );
}

#[test]
fn i13_g3_mutations_refuse_or_move_frozen_authority() {
    let baseline = i13_draft().freeze().expect("freeze").digest();

    let mut weakened = i13_draft();
    weakened
        .claims
        .iter_mut()
        .find(|claim| claim.id == "i13-relative-differential-character-invariants")
        .expect("differential-character theorem")
        .hypotheses = &["some finite complex is supplied"];
    assert_ne!(
        weakened.freeze().expect("weakened authority").digest(),
        baseline
    );

    let mut repartitioned = i13_draft();
    repartitioned
        .fixtures
        .iter_mut()
        .find(|fixture| fixture.id == "i13-global-robust-max-holdout")
        .expect("global robust heldout")
        .partition = Partition::Development;
    assert_ne!(
        repartitioned
            .freeze()
            .expect("partition amendment remains structurally representable")
            .digest(),
        baseline,
        "heldout-to-development movement must change frozen authority identity"
    );

    let mut missing_hypotheses = i13_draft();
    missing_hypotheses
        .claims
        .iter_mut()
        .find(|claim| claim.id == "i13-relative-differential-character-invariants")
        .expect("differential-character theorem")
        .hypotheses = &[];
    assert!(matches!(
        missing_hypotheses.freeze(),
        Err(FreezeRefusal::BlankField {
            field: "claim.hypotheses",
            ..
        })
    ));

    let mut correlated_oracle = i13_draft();
    correlated_oracle
        .claims
        .iter_mut()
        .find(|claim| claim.id == "i13-terminal-relative-class-audit")
        .expect("relative topology claim")
        .oracle
        .independent = false;
    assert!(matches!(
        correlated_oracle.freeze(),
        Err(FreezeRefusal::ProductionOracleReuse { .. })
    ));

    let mut relaxed = i13_draft();
    relaxed
        .claims
        .iter_mut()
        .find(|claim| claim.id == "i13-energy-consistent-material-interpolation")
        .expect("material interpolation")
        .tolerance = ToleranceSemantics::Interval { lo: 0.0, hi: 2.0 };
    assert_ne!(
        relaxed.freeze().expect("relaxed authority").digest(),
        baseline
    );

    let mut swapped_holdout = i13_draft();
    swapped_holdout
        .fixtures
        .iter_mut()
        .find(|pin| pin.id == "i13-maximal-adversary-max-holdout")
        .expect("maximal adversary holdout")
        .source = FixtureSource::AuthoredSpec {
        spec: "unauthorized post-result maximal holdout replacement",
    };
    assert_ne!(
        swapped_holdout
            .freeze()
            .expect("replacement authority")
            .digest(),
        baseline
    );

    let mut missing_policy = i13_draft();
    missing_policy
        .fixtures
        .retain(|pin| pin.id != CAMPAIGN_POLICY_FIXTURE);
    assert!(matches!(
        missing_policy.freeze(),
        Err(FreezeRefusal::OrphanDeck { deck, .. }) if deck == CAMPAIGN_POLICY_FIXTURE
    ));
}

#[test]
#[allow(clippy::too_many_lines)]
fn i13_amendments_invalidate_exact_reverse_dependency_authority() {
    let frozen = i13_draft().freeze().expect("freeze");

    let mut version_only = i13_draft();
    version_only.version = 2;
    let (version_successor, version_record) =
        frozen.amend(version_only).expect("version-only amendment");
    assert!(version_record.invalidated.is_empty());
    assert_eq!(
        (version_record.from_version, version_record.to_version),
        (1, 2)
    );
    assert_eq!(version_record.from_digest, frozen.digest());
    assert_eq!(version_record.to_digest, version_successor.digest());
    assert_ne!(version_record.from_digest, version_record.to_digest);

    let mut changed_claim = i13_draft();
    changed_claim.version = 2;
    changed_claim
        .claims
        .iter_mut()
        .find(|claim| claim.id == "i13-fixed-topology-coupled-adjoint")
        .expect("adjoint claim")
        .tolerance = ToleranceSemantics::Interval { lo: 0.0, hi: 2.0 };
    let (_, claim_record) = frozen.amend(changed_claim).expect("claim amendment");
    assert_eq!(
        claim_record
            .invalidated
            .into_iter()
            .collect::<BTreeSet<_>>(),
        BTreeSet::from([
            "i13-adjoint-optimizer-core".to_string(),
            "i13-fixed-topology-coupled-adjoint".to_string(),
        ])
    );

    for (fixture_id, expected_claim, expected_leaf) in [
        (
            "i13-manufacturing-insulation-core-holdout",
            "i13-manufacturability-insulation-audit",
            "i13-manufacturability-core",
        ),
        (
            "i13-route-realization-core-holdout",
            "i13-route-netlist-realization-audit",
            "i13-route-realization-core",
        ),
        (
            "i13-robust-qoi-core-holdout",
            "i13-robust-guarded-machine-qois",
            "i13-robust-qoi-core",
        ),
        (
            "i13-multiphysics-closure-core-holdout",
            "i13-power-force-thermal-stress-closure",
            "i13-multiphysics-closure-core",
        ),
    ] {
        let mut changed_split_holdout = i13_draft();
        changed_split_holdout.version = 2;
        changed_split_holdout
            .fixtures
            .iter_mut()
            .find(|pin| pin.id == fixture_id)
            .unwrap_or_else(|| panic!("missing split holdout '{fixture_id}'"))
            .source = FixtureSource::AuthoredSpec {
            spec: "replacement independently governed split holdout",
        };
        let (_, split_record) = frozen
            .amend(changed_split_holdout)
            .expect("split holdout amendment");
        assert_eq!(
            split_record
                .invalidated
                .into_iter()
                .collect::<BTreeSet<_>>(),
            BTreeSet::from([expected_claim.to_string(), expected_leaf.to_string()]),
            "holdout '{fixture_id}' crossed an evidence-owner boundary"
        );
    }

    let mut changed_holdout = i13_draft();
    changed_holdout.version = 2;
    changed_holdout
        .fixtures
        .iter_mut()
        .find(|pin| pin.id == "i13-topology-continuation-max-holdout")
        .expect("continuation holdout")
        .source = FixtureSource::AuthoredSpec {
        spec: "replacement topology-continuation holdout with a new semantic identity",
    };
    let (_, holdout_record) = frozen.amend(changed_holdout).expect("holdout amendment");
    assert_eq!(
        holdout_record
            .invalidated
            .into_iter()
            .collect::<BTreeSet<_>>(),
        BTreeSet::from([
            "i13-certified-topology-change-continuation".to_string(),
            "i13-topology-continuation-max-leaf".to_string(),
        ])
    );

    let mut changed_theorem_card = i13_draft();
    changed_theorem_card.version = 2;
    changed_theorem_card
        .fixtures
        .iter_mut()
        .find(|pin| pin.id == "i13-differential-character-theorem-card")
        .expect("differential-character theorem card")
        .source = FixtureSource::AuthoredSpec {
        spec: "replacement differential-character proposition target",
    };
    let (_, theorem_record) = frozen
        .amend(changed_theorem_card)
        .expect("theorem-card amendment");
    assert_eq!(
        theorem_record
            .invalidated
            .into_iter()
            .collect::<BTreeSet<_>>(),
        BTreeSet::from([
            "i13-differential-character-theorem-max".to_string(),
            "i13-maximal-falsifier-max-leaf".to_string(),
            "i13-maximal-theorem-counterexample-search".to_string(),
            "i13-relative-differential-character-invariants".to_string(),
        ])
    );

    let mut changed_route_theorem_card = i13_draft();
    changed_route_theorem_card.version = 2;
    changed_route_theorem_card
        .fixtures
        .iter_mut()
        .find(|pin| pin.id == "i13-topology-route-theorem-card")
        .expect("topology-route theorem card")
        .source = FixtureSource::AuthoredSpec {
        spec: "replacement topology-to-route proposition target",
    };
    let (_, route_theorem_record) = frozen
        .amend(changed_route_theorem_card)
        .expect("route-theorem-card amendment");
    assert_eq!(
        route_theorem_record
            .invalidated
            .into_iter()
            .collect::<BTreeSet<_>>(),
        BTreeSet::from([
            "i13-maximal-falsifier-max-leaf".to_string(),
            "i13-maximal-theorem-counterexample-search".to_string(),
            "i13-topology-route-theorem-max".to_string(),
            "i13-topology-to-route-performance-theorem".to_string(),
        ])
    );

    let mut changed_falsifier_holdout = i13_draft();
    changed_falsifier_holdout.version = 2;
    changed_falsifier_holdout
        .fixtures
        .iter_mut()
        .find(|pin| pin.id == "i13-maximal-adversary-max-holdout")
        .expect("maximal falsifier holdout")
        .source = FixtureSource::AuthoredSpec {
        spec: "replacement maximal falsifier holdout",
    };
    let (_, falsifier_record) = frozen
        .amend(changed_falsifier_holdout)
        .expect("falsifier-deck amendment");
    assert_eq!(
        falsifier_record
            .invalidated
            .into_iter()
            .collect::<BTreeSet<_>>(),
        BTreeSet::from([
            "i13-maximal-falsifier-max-leaf".to_string(),
            "i13-maximal-theorem-counterexample-search".to_string(),
        ])
    );

    let mut changed_waiver = i13_draft();
    changed_waiver.version = 2;
    changed_waiver.waivers[0].promotion_effect =
        "test-only changed industrial-validation waiver scope";
    let (_, waiver_record) = frozen.amend(changed_waiver).expect("waiver amendment");
    assert_eq!(
        waiver_record
            .invalidated
            .into_iter()
            .collect::<BTreeSet<_>>(),
        BTreeSet::from([
            "i13-governed-industrial-machine-validation".to_string(),
            "i13-industrial-validation-max".to_string(),
        ])
    );

    // Generic manifest structure deliberately accepts this same-ID external
    // replacement so reverse-dependency invalidation remains testable. The
    // campaign policy is the authority boundary: an all-zero root plus waiver
    // removal is not a verified industrial-data discharge.
    let mut structural_only_discharge = i13_draft();
    structural_only_discharge.version = 2;
    structural_only_discharge.waivers.clear();
    structural_only_discharge.fixtures.push(FixturePin {
        id: "i13-external-industrial-ap242-pack",
        source: FixtureSource::External {
            digest_hex: "0000000000000000000000000000000000000000000000000000000000000000",
        },
        partition: Partition::HeldOut,
    });
    let (structural_successor, structural_record) = frozen
        .amend(structural_only_discharge)
        .expect("same-ID all-zero external replacement remains structurally representable");
    assert!(structural_successor.waivers().is_empty());
    let structural_pack = fixture(
        structural_successor.fixtures(),
        "i13-external-industrial-ap242-pack",
    );
    assert_eq!(structural_pack.partition, Partition::HeldOut);
    assert!(matches!(
        structural_pack.source,
        FixtureSource::External {
            digest_hex: "0000000000000000000000000000000000000000000000000000000000000000"
        }
    ));
    assert_eq!(
        structural_record
            .invalidated
            .into_iter()
            .collect::<BTreeSet<_>>(),
        BTreeSet::from([
            "i13-governed-industrial-machine-validation".to_string(),
            "i13-industrial-validation-max".to_string(),
        ])
    );
    let successor_policy = authored_spec(fixture(
        structural_successor.fixtures(),
        CAMPAIGN_POLICY_FIXTURE,
    ));
    for required in [
        "raw or all-zero digest",
        "arbitrary same-ID External fixture",
        "waiver removal alone",
        "structurally frozen successor",
        "structural successor lacking the verified fs-vvreg transaction grants no discharge or promotion authority",
    ] {
        assert!(
            successor_policy.contains(required),
            "structural successor lost negative-authority guard '{required}'"
        );
    }

    let mut changed_policy = i13_draft();
    changed_policy.version = 2;
    changed_policy
        .fixtures
        .iter_mut()
        .find(|pin| pin.id == CAMPAIGN_POLICY_FIXTURE)
        .expect("policy")
        .source = FixtureSource::AuthoredSpec {
        spec: "I13_CAMPAIGN_POLICY_V2 unauthorized no-collapse edit",
    };
    let (_, policy_record) = frozen.amend(changed_policy).expect("policy amendment");
    assert_eq!(
        policy_record
            .invalidated
            .iter()
            .map(String::as_str)
            .collect::<BTreeSet<_>>(),
        all_authority_ids()
    );
    assert_eq!(policy_record.invalidated.len(), 31);

    let mut changed_explicits = i13_draft();
    changed_explicits.version = 2;
    changed_explicits.explicits.capabilities =
        "test-only changed I13 capability authority; every descendant invalidates";
    let (_, explicits_record) = frozen.amend(changed_explicits).expect("explicit amendment");
    assert_eq!(
        explicits_record
            .invalidated
            .iter()
            .map(String::as_str)
            .collect::<BTreeSet<_>>(),
        all_authority_ids()
    );
    assert_eq!(explicits_record.invalidated.len(), 31);

    let mut changed_title = i13_draft();
    changed_title.version = 2;
    changed_title.title = "test-only changed I13 campaign semantics";
    let (_, title_record) = frozen.amend(changed_title).expect("title amendment");
    assert_eq!(
        title_record
            .invalidated
            .iter()
            .map(String::as_str)
            .collect::<BTreeSet<_>>(),
        all_authority_ids()
    );
}

#[test]
fn i13_replay_and_evidence_policy_is_agent_complete_and_bounded() {
    let policy = policy_spec();
    for required in [
        "typed canonical JSONL observability mirror",
        "manifest/claim/leaf/fixture/candidate/oracle/toolchain/ISA/seed/checkpoint ids",
        "timing/memory/cancellation latency",
        "event kind is at most 256 bytes",
        "complete encoded record including envelope and LF is at most 1048576 bytes",
        "durable flush is at most 1024 rows and 4194304 encoded bytes",
        "run/scope is at most 65536 rows and 67108864 encoded bytes",
        "governor retains at most 268435456 bytes",
        "non-borrowably reserves the base worst-case terminal rows and bytes",
        "each later observer, descendant, governance-transition or protected-access registration must atomically acquire its incremental reservation",
        "one terminal transition row for each permitted governance transition",
        "Prepared-plus-terminal pair for each permitted protected-access operation",
        "ordinary deterministic event schedule plus that exact row/byte reserve",
        "ordinary in-memory buffer is exactly one flush",
        "independently bounded priority writer",
        "stops at its next bounded tile boundary and durably flushes the prefix",
        "No created event is silently dropped, overwritten, truncated, arrival-time sampled",
        "lossless and non-coalescible",
        "deterministically aggregated before event creation by frozen logical tile/shard/window",
        "16777216-evaluation supergrammar retains full replay evidence",
        "receipt-qualified sink refusal contributes an InfrastructureFailed/bit4 witness",
        "slow sink consumes the frozen wall budget",
        "durably persisted wall-deadline receipt contributes a TimedOut/bit3 witness",
        "Either blocks promotion and, on CanonicalTerminalPath, yields PartialEvidence only when a durable authenticated prefix/FailureBundle exists and NoEvidence otherwise",
        "only FINAL_EXECUTION_DISPOSITION_JOIN_V1 derives final disposition",
        "Protected access is denied unless its Prepared row is durable",
        "If the reserved primary terminal transaction cannot persist, the independent emergency terminal journal and authenticated durable-prefix recovery routes are mandatory",
        "while recovery remains possible the attempt is TerminalPersistencePending",
        "Exhaustion of all such routes is TerminalPersistenceCatastrophe, not a canonical NoEvidence result",
        "out-of-band incident diagnostic is explicitly non-authoritative",
        "Observer, descendant, governance-transition and protected-access registration each atomically consumes its predeclared row/byte slots in both lane reservations",
        "no registration can expand p_i",
        "governance.transition consumes one terminal event/receipt plus refusal/failure capacity per permitted transition",
        "governance.protected_access consumes one Prepared-plus-terminal pair per permitted operation",
        "Terminal selection freezes all four membership roots/counts and their exact mirrored reservations",
        "primary/emergency segments cannot borrow from each other",
        "One emergency segment reserves the exact terminal slot plus both downstream receipts",
        "other emergency segments carry typed prefix events under EMERGENCY_PREFIX_FAILOVER_V1",
        "p_i disjoint 4194304-byte segments in the primary priority writer and p_i equally sized segments in the emergency priority writer",
        "global_reserved=checked_sum_i(run_reserve_i) must be <=268435456",
        "p_i and the artifact reservation never grow after admission",
        "priority_observation_latency_ns=checked_add(checked_add(worst_case_request_to_next_priority_poll_ns,worst_case_observer_reaction_and_enqueue_ns)",
        "priority_segment_count covers the full serial path from primary-request creation through the slowest sealed observer's durable receipt",
        "cap+1, asymmetric capacity, arithmetic failure, borrowing, priority starvation, stale service receipt or non-independent emergency service",
        "Only failure of both mirrored lanes and every authenticated replay route permits TerminalPersistenceCatastrophe",
        "No path permits silent loss or favorable authority",
        "G5 compares canonical logical sequences, typed identity bytes, event-chain/Merkle roots",
        "calibrated timestamps remain telemetry only",
        "ranked next actions",
        "complete authenticated evidence transaction or none",
        "every minimized FailureBundle or countermodel",
        "quiet-host receipt",
        "cross-ISA bit identity is not required",
    ] {
        assert!(
            policy.contains(required),
            "policy lost operational clause '{required}'"
        );
    }
    let frozen = i13_draft().freeze().expect("freeze");
    for row in frozen.obligations() {
        assert!(row.replay_command().contains("--replay <artifact-id>"));
        assert!(row.dsr_lane().ends_with("dsr quality --tool frankensim"));
        assert!(row.g4_schedule().contains("publishes"));
        assert!(row.g4_schedule().contains("or none"));
    }
}

#[test]
fn i13_terminal_authority_prior_ratchets_are_exact_and_explicitly_superseded_where_needed() {
    let policy = policy_spec();
    assert_eq!(b"I13_TERMINAL_EVENT_SEMANTIC_AUTHORITY_V1\0".len(), 41);
    assert_eq!(b"I13_RECEIPT_TO_TERMINAL_BINDING_V1\0".len(), 35);
    assert_eq!(b"I13_CAPACITY_BREACH_FACT_V1\0".len(), 28);
    assert_eq!(b"I13_CAPACITY_BREACH_SELECTION_V2\0".len(), 33);
    assert_eq!(b"I13_POST_SELECTION_BREACH_ANCHOR_V1\0".len(), 36);
    assert_eq!(b"I13_LANE_COMMON_PREFIX_CERTIFICATE_V1\0".len(), 38);
    assert_eq!(b"I13_LANE_APPEND_ACK_V1\0".len(), 23);
    assert_eq!(28 + 1 + 2 + 1 + 2 * 32 + 3 * 8 + 32 + 1, 153);
    assert_eq!(153 + 76, 229);
    assert_eq!(16 + 8 + 229 + 256 * 32, 8_445);
    assert_eq!(35 + 3 * 32 + 2, 133);
    assert_eq!(
        33 + 2 * 32 + 2 + 2 * 32 + 2 * 8 + 32 + 8 + 76 + 4 * 40 + 32 + 2 * 8 + 32,
        535
    );
    assert_eq!(36 + 2 * 8 + 32 + 8 + 2 * 32 + 2 + 2 * 32 + 8 + 32, 262);
    assert_eq!(1 + 32 + 2, 35);
    assert_eq!(177_408 + 2_816 + 1, 180_225);
    assert_eq!(180_225 * 512, 92_275_200);
    assert_eq!(67_108_864 + 2 * 92_275_200 + 2 * 4_194_304, 260_047_872);
    assert_eq!(268_435_456 - 260_047_872, 8_387_584);
    assert_eq!(1_047_813 - 35, 1_047_778);
    assert_eq!(295 + 32, 327);
    assert_eq!(38 + 2 * 8 + 3 * 32 + 2 * 8 + 2 * 32 + 32 + 2 * 8 + 32, 310);
    assert_eq!(8 + 16_384 * 66, 1_081_352);
    assert_eq!(
        23 + 2 * 8 + 1 + 2 * 32 + 2 * 8 + 32 + 8 + 32 + 2 * 8 + 32,
        240
    );

    for exact in [
        "TERMINAL_CAUSAL_SCHEMA_SEMANTIC_AUTHORITY_V1",
        "40 exact ASCII bytes I13_TERMINAL_EVENT_SEMANTIC_AUTHORITY_V1 followed by one byte 0x00",
        "raw32(CanonicalEventSchemaDigest)||raw32(event_semantic_conformance_digest)",
        "permissive all-None profile for a causal kind",
        "ReceiptToTerminalRecordBindingBytesV1 is the 34 exact ASCII bytes I13_RECEIPT_TO_TERMINAL_BINDING_V1",
        "exactly 133 bytes",
        "sparse-Merkle map keyed exactly by authority_receipt_digest",
        "AuthenticatedBreachFactFrontier=4||U16LE(field_ordinal)",
        "max(terminal-parent sequences, consumed external dependency sequences, with empty maximum 0)",
        "CANONICAL_PARENT_MULTIPROOF_RECONSTRUCTION_V1",
        "Build the union trie of their 256-bit sparse-Merkle paths",
        "if T occupies both children, the proof contains no entry and both child hashes are computed",
        "opposite child is canonical empty_{d+1}",
        "consumes every and only required nonempty sibling entry",
        "checked_mul(64,256)=16384 internal nodes",
        "CAPACITY_BREACH_FACT_FRONTIER_V1",
        "27 exact ASCII bytes I13_CAPACITY_BREACH_FACT_V1 followed by one byte 0x00",
        "exactly 153 or 229 bytes",
        "state is Open, BreachCollecting or Closed",
        "count is in 1..=257",
        "maximum 229-byte update proof is exactly 8445 bytes",
        "CapacityBreachSelectionReceiptBytesV2 is the 32 exact ASCII bytes I13_CAPACITY_BREACH_SELECTION_V2",
        "exactly 535 bytes",
        "older singular CandidateCapacityBreachBytesV1 and 501-byte receipt describe only the opening fact/legacy draft",
        "POST_SELECTION_BREACH_ANCHOR_AND_CAP_V2",
        "PostSelectionCapacityBreachUnionV2 is exactly Absent=0 or Present=1||raw32(closed_post_selection_breach_fact_root)||U16LE(breach_fact_count)",
        "one or 35 bytes",
        "35 exact ASCII bytes I13_POST_SELECTION_BREACH_ANCHOR_V1 followed by one byte 0x00",
        "exactly 262 bytes",
        "Ordinary nonbarrier admission retains its full 177407-event capacity",
        "MAX_TERMINAL_EVENT_COUNT_V2=checked_add(177408,checked_add(2816,1))=180225",
        "MAX_TERMINAL_EVENT_ARTIFACT_BYTES_V2=checked_mul(180225,512)=92275200",
        "260047872",
        "8387584 bytes slack",
        "zero-parent Present admits at most 1047778 payload bytes",
        "CAPACITY_BREACH_SELECTION_SUBBATCH_V1",
        "terminal_selection_transaction_max_logical_sequence=M_sel is the maximum of exactly SelectionBatchV1",
        "PostSelectionWitnessBatchV1",
        "explicitly excluded from SelectionBatchV1",
        "LANE_COMMON_PREFIX_FAILOVER_V2",
        "37 exact ASCII bytes I13_LANE_COMMON_PREFIX_CERTIFICATE_V1 followed by one byte 0x00",
        "exactly 310 bytes",
        "LaneEpochTransitionBytesV2 is the existing 295-byte V1 projection",
        "exactly 327 bytes",
        "compare-and-check the exact committed Emergency/1 head unchanged as an equality fence",
        "LANE_APPEND_ACK_AND_V2_DIGEST_AUTHORITY_V1",
        "22 exact ASCII bytes I13_LANE_APPEND_ACK_V1 followed by one byte 0x00",
        "exactly 240 bytes",
        "lane_append_ack_digest=BLAKE3::derive_key('org.frankensim.i13.lane-append-ack.v1',LaneAppendAckBytesV1)",
        "terminal_generation=terminal_count=artifact_event_count",
        "durable_interval_lower=min(primary_lower,emergency_lower)",
        "lane_common_prefix_certificate_digest=BLAKE3::derive_key('org.frankensim.i13.lane-common-prefix-certificate.v1',LaneCommonPrefixCertificateBytesV1)",
        "lane_epoch_transition_digest_v2=BLAKE3::derive_key('org.frankensim.i13.lane-epoch-transition.v2',LaneEpochTransitionBytesV2)",
        "successor_lane_epoch_head_digest_v2=BLAKE3::derive_key('org.frankensim.i13.lane-epoch-head.v1'",
        "no V1 digest alias has V2 authority",
        "239/240/241 acknowledgement size",
    ] {
        assert!(policy.contains(exact), "terminal ratchet lost '{exact}'");
    }
    assert!(policy.contains("I13_FRESH_V2_AUTHORITY_PRECEDENCE"));
    assert!(policy.contains("92275200-byte structure-only lane reserve"));
    assert!(policy.contains("1346073904 equation"));
}

#[test]
fn i13_terminal_v2_authority_is_acyclic_closed_and_capacity_complete() {
    let policy = policy_spec();

    for (bytes, expected) in [
        (
            b"I13_TERMINAL_CAPACITY_RESERVATION_INTENT_V2\0".as_slice(),
            44,
        ),
        (
            b"I13_TERMINAL_CAPACITY_RESERVATION_RECEIPT_V2\0".as_slice(),
            45,
        ),
        (b"I13_TERMINAL_AUTHORITY_ADMISSION_CORE_V2\0".as_slice(), 41),
        (b"I13_TERMINAL_AUTHORITY_ACTIVATION_V2\0".as_slice(), 37),
        (b"I13_ARBITRATED_SELECTION_V2\0".as_slice(), 28),
        (b"I13_TERMINAL_SELECTION_RECEIPT_V2\0".as_slice(), 34),
        (b"I13_BREACH_FACT_ARTIFACT_V1\0".as_slice(), 28),
        (b"I13_TERMINAL_EVENT_ARTIFACT_CHUNK_V2\0".as_slice(), 37),
        (b"I13_DURABLE_APPEND_OPERATION_RECEIPT_V1\0".as_slice(), 40),
        (b"I13_PREFIX_VERIFICATION_RECEIPT_V1\0".as_slice(), 35),
        (b"I13_SEMANTIC_INVENTORY_V2\0".as_slice(), 26),
        (b"I13_SEMANTIC_COMPLETENESS_MATRIX_V2\0".as_slice(), 36),
        (b"I13_EVENT_SEMANTIC_CONFORMANCE_V2\0".as_slice(), 34),
        (b"I13_TERMINAL_EVENT_SEMANTIC_AUTHORITY_V2\0".as_slice(), 41),
        (b"I13_TERMINAL_AUTHORITY_HEAD_V2\0".as_slice(), 31),
        (b"I13_HOLDOUT_GENERATOR_CONTRACT_V2\0".as_slice(), 34),
        (b"I13_HOLDOUT_PAYLOAD_SCHEMA_V2\0".as_slice(), 30),
        (b"I13_TERMINAL_CAUSAL_SCHEMA_V2\0".as_slice(), 30),
        (b"I13_FINALIZATION_CORE_V2\0".as_slice(), 25),
        (b"I13_FINALIZATION_ENVELOPE_V2\0".as_slice(), 29),
        (b"I13_EMERGENCY_TERMINAL_SLOT_V2\0".as_slice(), 31),
        (b"I13_CONFORMANCE_EVALUATOR_SET_V2\0".as_slice(), 33),
        (b"I13_GOVERNED_KAT_SET_V2\0".as_slice(), 24),
        (b"I13_SEMANTIC_SOURCE_CATALOG_V2\0".as_slice(), 31),
        (b"I13_SEMANTIC_SOURCE_OCCURRENCE_SET_V2\0".as_slice(), 38),
        (
            b"I13_SEMANTIC_OCCURRENCE_CONFORMANCE_RECEIPT_V2\0".as_slice(),
            47,
        ),
        (
            b"I13_SEMANTIC_OCCURRENCE_CONFORMANCE_SET_V2\0".as_slice(),
            43,
        ),
        (b"I13_TERMINAL_EVENT_SEMANTIC_CATALOG_V2\0".as_slice(), 39),
        (
            b"I13_SEMANTIC_CATALOG_GENERATION_RECEIPT_V2\0".as_slice(),
            43,
        ),
        (b"I13_SEMANTIC_CATALOG_RECEIPT_SET_V2\0".as_slice(), 36),
        (b"I13_RECEIPT_INDEX_UPDATE_BATCH_V2\0".as_slice(), 34),
        (b"I13_GLOBAL_TERMINAL_GOVERNOR_HEAD_V2\0".as_slice(), 37),
        (b"I13_GLOBAL_GOVERNOR_AUTHORITY_V2\0".as_slice(), 33),
        (b"I13_ACTIVE_RESERVATION_SET_V2\0".as_slice(), 30),
        (b"I13_LANE_RESERVATION_RECEIPT_V2\0".as_slice(), 32),
        (b"I13_LANE_SERVICE_RESERVATION_HEAD_V2\0".as_slice(), 37),
        (b"I13_LANE_RESERVATION_ISSUANCE_RECEIPT_V2\0".as_slice(), 41),
        (b"I13_LANE_RESERVATION_ATTEMPT_SET_V2\0".as_slice(), 36),
        (b"I13_LANE_RESERVATION_AUTHORITY_SET_V2\0".as_slice(), 38),
        (b"I13_TERMINAL_EVENT_PAYLOAD_ARTIFACT_V2\0".as_slice(), 39),
        (b"I13_RECEIPT_INDEX_ARTIFACT_V2\0".as_slice(), 30),
        (b"I13_OPERATION_PAYLOAD_AUTHORITY_V2\0".as_slice(), 35),
        (b"I13_OPERATION_PAYLOAD_SCHEMA_MATRIX_V2\0".as_slice(), 39),
        (b"I13_TERMINAL_APPEND_DELTA_V2\0".as_slice(), 29),
        (b"I13_TERMINAL_APPEND_REQUEST_V2\0".as_slice(), 31),
        (b"I13_TERMINAL_DURABILITY_COMMIT_V2\0".as_slice(), 34),
        (b"I13_BREACH_PERSISTENCE_HEAD_V2\0".as_slice(), 31),
        (
            b"I13_BREACH_ARTIFACT_PERSISTENCE_RECEIPT_V2\0".as_slice(),
            43,
        ),
        (
            b"I13_BREACH_ARTIFACT_MIRROR_CERTIFICATE_V2\0".as_slice(),
            42,
        ),
        (b"I13_BREACH_FACT_DURABILITY_BUNDLE_V2\0".as_slice(), 37),
        (
            b"I13_TERMINAL_CAPACITY_SETTLEMENT_RECEIPT_V2\0".as_slice(),
            44,
        ),
        (b"I13_HOLDOUT_REALIZATION_COMMITMENT_V2\0".as_slice(), 38),
        (b"I13_TERMINAL_PUBLICATION_AUTHORITY_V2\0".as_slice(), 38),
        (b"I13_TERMINAL_DURABILITY_HEAD_V2\0".as_slice(), 32),
        (b"I13_PUBLICATION_JOURNAL_RECEIPT_V2\0".as_slice(), 35),
        (b"I13_TERMINAL_PUBLICATION_CERTIFICATE_V2\0".as_slice(), 40),
        (b"I13_CAPABILITY_QUIESCENCE_RECEIPT_V2\0".as_slice(), 37),
        (b"I13_LANE_RESERVATION_RELEASE_RECEIPT_V2\0".as_slice(), 40),
        (b"I13_LANE_CLOSED_PLAN_SET_V2\0".as_slice(), 28),
        (b"I13_CAPABILITY_STATE_HEAD_V2\0".as_slice(), 29),
        (b"I13_CAPABILITY_REVOCATION_RECEIPT_V2\0".as_slice(), 37),
        (b"I13_PREACTIVATION_ABORT_EVIDENCE_V2\0".as_slice(), 36),
        (b"I13_RETAINED_ARTIFACT_SET_V2\0".as_slice(), 29),
        (b"I13_RETENTION_STORE_RECEIPT_SET_V2\0".as_slice(), 35),
        (b"I13_CANONICAL_FINALIZATION_EVIDENCE_V2\0".as_slice(), 39),
        (
            b"I13_TERMINAL_PERSISTENCE_CATASTROPHE_EVIDENCE_V2\0".as_slice(),
            49,
        ),
        (
            b"I13_DURABLE_RETENTION_TRANSFER_RECEIPT_V2\0".as_slice(),
            42,
        ),
        (b"I13_REVOKED_CAPABILITY_SET_V2\0".as_slice(), 30),
        (b"I13_FAILURE_INVENTORY_V2\0".as_slice(), 25),
        (b"I13_REPLAY_ROUTE_LEDGER_V2\0".as_slice(), 27),
        (b"I13_RETENTION_STORE_HEAD_V2\0".as_slice(), 28),
        (b"I13_RETENTION_STORE_RECEIPT_V2\0".as_slice(), 31),
        (b"I13_RETENTION_SET_INVENTORY_V2\0".as_slice(), 31),
        (b"I13_LANE_ISSUANCE_LEASE_V2\0".as_slice(), 27),
        (b"I13_PREACTIVATION_SERVICE_CATASTROPHE_V2\0".as_slice(), 41),
        (
            b"I13_PROTOCOL_FIXED_TERMINAL_EVENT_REGISTRY_V2\0".as_slice(),
            46,
        ),
        (b"I13_OPERATION_PAYLOAD_ARTIFACT_V2\0".as_slice(), 34),
        (b"I13_TERMINAL_APPEND_DELTA_ARTIFACT_V2\0".as_slice(), 38),
        (b"I13_APPEND_METADATA_ARTIFACT_V2\0".as_slice(), 32),
        (b"I13_APPEND_PHASE_HEAD_V2\0".as_slice(), 25),
        (b"I13_APPEND_PHASE_LEASE_RECEIPT_V2\0".as_slice(), 34),
        (b"I13_APPEND_PHASE_RECEIPT_ARTIFACT_V2\0".as_slice(), 37),
        (b"I13_APPEND_PHASE_LEASE_V3\0".as_slice(), 26),
        (b"I13_APPEND_COMPONENT_BYTE_PLAN_V3\0".as_slice(), 34),
        (b"I13_APPEND_PHASE_HEAD_V3\0".as_slice(), 25),
        (b"I13_APPEND_PHASE_LEASE_RECEIPT_V3\0".as_slice(), 34),
        (b"I13_APPEND_PHASE_RECEIPT_ARTIFACT_V3\0".as_slice(), 37),
        (b"I13_APPEND_PHASE_HEAD_ARTIFACT_V3\0".as_slice(), 34),
        (b"I13_APPEND_PHASE_LEASE_ARTIFACT_V3\0".as_slice(), 35),
        (b"I13_APPEND_BYTE_PLANS_ARTIFACT_V3\0".as_slice(), 34),
        (b"I13_APPEND_PHASE_AUTHORITY_SET_V3\0".as_slice(), 34),
        (b"I13_APPEND_SLOT_LEDGER_V3\0".as_slice(), 26),
        (
            b"I13_LOGICAL_APPEND_SLOT_ALLOCATION_INTENT_V3\0".as_slice(),
            45,
        ),
        (
            b"I13_LOGICAL_APPEND_SLOT_AUTHORITY_HEAD_V3\0".as_slice(),
            42,
        ),
        (
            b"I13_LOGICAL_APPEND_SLOT_ALLOCATION_RECEIPT_V3\0".as_slice(),
            46,
        ),
        (
            b"I13_LOGICAL_APPEND_SLOT_HISTORY_ARTIFACT_V3\0".as_slice(),
            44,
        ),
        (b"I13_APPEND_PHASE_AUTHORITY_HEAD_V3\0".as_slice(), 35),
        (b"I13_LANE_CONTROL_CLOSURE_ARTIFACT_V3\0".as_slice(), 37),
        (b"I13_LANE_BREACH_SLOT_SETTLEMENT_V3\0".as_slice(), 35),
        (b"I13_LANE_BREACH_SLOT_RESERVATION_V3\0".as_slice(), 36),
        (b"I13_LANE_RELEASE_AUTHORIZATION_V3\0".as_slice(), 34),
        (b"I13_INFRASTRUCTURE_CAPACITY_SET_V3\0".as_slice(), 35),
        (
            b"I13_INFRASTRUCTURE_SERVICE_CAPACITY_INTENT_V3\0".as_slice(),
            46,
        ),
        (
            b"I13_INFRASTRUCTURE_SERVICE_CAPACITY_HEAD_V3\0".as_slice(),
            44,
        ),
        (
            b"I13_INFRASTRUCTURE_SERVICE_CAPACITY_UPDATE_RECEIPT_V3\0".as_slice(),
            54,
        ),
        (
            b"I13_TERMINAL_CAPACITY_PROFILE_EXTENSION_V3\0".as_slice(),
            43,
        ),
        (
            b"I13_TERMINAL_CAPACITY_RESERVATION_INTENT_V3\0".as_slice(),
            44,
        ),
        (
            b"I13_TERMINAL_CAPACITY_RESERVATION_RECEIPT_V3\0".as_slice(),
            45,
        ),
        (b"I13_LANE_RESERVATION_RECEIPT_V3\0".as_slice(), 32),
        (b"I13_ORDINARY_EVIDENCE_ARTIFACT_V3\0".as_slice(), 34),
        (b"I13_LANE_SERVICE_CAPACITY_HEAD_V3\0".as_slice(), 34),
        (
            b"I13_LANE_SERVICE_CAPACITY_UPDATE_RECEIPT_V3\0".as_slice(),
            44,
        ),
        (b"I13_TERMINAL_ROUTE_FAILURE_MATRIX_V2\0".as_slice(), 37),
        (b"I13_TERMINAL_REPLAY_ROUTE_REGISTRY_V2\0".as_slice(), 38),
        (b"I13_TERMINAL_ROUTE_SPEC_V2\0".as_slice(), 27),
        (b"I13_ADMISSIBLE_REPLAY_ROUTE_UNIVERSE_V2\0".as_slice(), 40),
        (b"I13_TERMINAL_FAILURE_RECEIPT_V2\0".as_slice(), 32),
        (b"I13_ROUTE_INAPPLICABILITY_PROOF_V2\0".as_slice(), 35),
        (b"I13_TERMINAL_ROUTE_OUTCOME_RECEIPT_V2\0".as_slice(), 38),
        (b"I13_ARM_RETENTION_REQUIREMENT_SET_V2\0".as_slice(), 37),
        (b"I13_STORE_CAPACITY_HEAD_V2\0".as_slice(), 27),
        (b"I13_RETAINED_ARTIFACT_SET_V3\0".as_slice(), 29),
        (b"I13_ARM_RETENTION_REQUIREMENT_SET_V3\0".as_slice(), 37),
        (b"I13_RETENTION_CONFORMANCE_RECEIPT_V3\0".as_slice(), 37),
        (b"I13_STORE_CAPACITY_HEAD_V3\0".as_slice(), 27),
        (b"I13_STORE_CAPACITY_RESERVATION_INTENT_V3\0".as_slice(), 41),
        (b"I13_STORE_CAPACITY_UPDATE_RECEIPT_V3\0".as_slice(), 37),
        (
            b"I13_STORE_CAPACITY_UPDATE_PROOF_BUNDLE_V3\0".as_slice(),
            42,
        ),
        (b"I13_CONTROL_PLANE_LEDGER_ENTRY_V3\0".as_slice(), 34),
        (b"I13_CONTROL_PLANE_LEDGER_HEAD_V3\0".as_slice(), 33),
        (
            b"I13_TERMINAL_CAPACITY_SETTLEMENT_RECEIPT_V3\0".as_slice(),
            44,
        ),
        (b"I13_CONTROL_PLANE_SETTLEMENT_MARKER_V3\0".as_slice(), 39),
        (b"I13_BREACH_EVIDENCE_BUNDLE_V3\0".as_slice(), 30),
        (b"I13_LANE_RESERVATION_RELEASE_RECEIPT_V3\0".as_slice(), 40),
        (b"I13_RETENTION_SOURCE_SHARD_RECEIPT_V3\0".as_slice(), 38),
        (
            b"I13_RETENTION_SOURCE_SHARD_RECEIPT_SET_V3\0".as_slice(),
            42,
        ),
        (
            b"I13_RETENTION_PRODUCER_INSERTION_RECEIPT_V3\0".as_slice(),
            44,
        ),
        (b"I13_CANONICAL_FINALIZATION_EVIDENCE_V3\0".as_slice(), 39),
        (b"I13_DESTINATION_RETENTION_RECEIPT_SET_V3\0".as_slice(), 41),
        (
            b"I13_DURABLE_RETENTION_TRANSFER_RECEIPT_V3\0".as_slice(),
            42,
        ),
        (b"I13_APPEND_METADATA_ARTIFACT_V3\0".as_slice(), 32),
        (
            b"I13_APPEND_PHASE_CATASTROPHE_SURVIVOR_SET_V3\0".as_slice(),
            45,
        ),
        (b"I13_APPEND_PHASE_CLOSURE_SUM_V3\0".as_slice(), 32),
        (b"I13_APPEND_PHASE_LATEST_PAIR_SNAPSHOT_V3\0".as_slice(), 41),
        (b"I13_BOOTSTRAP_EVIDENCE_BUNDLE_V3\0".as_slice(), 33),
        (
            b"I13_BOOTSTRAP_RESERVATION_FAILURE_PAYLOAD_V3\0".as_slice(),
            45,
        ),
        (b"I13_BREACH_SLOT_AUTHORITY_MEMBER_V3\0".as_slice(), 36),
        (b"I13_CAPACITY_INDEX_HEAD_V3\0".as_slice(), 27),
        (
            b"I13_CAPACITY_INDEX_RESERVATION_SUBJECT_V3\0".as_slice(),
            42,
        ),
        (b"I13_CAPACITY_INDEX_ALLOCATION_REQUEST_V3\0".as_slice(), 41),
        (b"I13_CAPACITY_INDEX_ALLOCATION_RECEIPT_V3\0".as_slice(), 41),
        (
            b"I13_CAPACITY_INDEX_ALLOCATION_EVIDENCE_BUNDLE_V3\0".as_slice(),
            49,
        ),
        (
            b"I13_CAPACITY_INDEX_ALLOCATION_RECEIPT_PAIR_V3\0".as_slice(),
            46,
        ),
        (
            b"I13_CAPACITY_INDEX_ALLOCATION_AUTHORITY_SET_V3\0".as_slice(),
            47,
        ),
        (b"I13_CAPACITY_RECEIPT_JOURNAL_HEAD_V3\0".as_slice(), 37),
        (b"I13_CAPACITY_RECEIPT_SLOT_RECORD_V3\0".as_slice(), 36),
        (
            b"I13_CAPACITY_RECEIPT_JOURNAL_TRANSITION_V3\0".as_slice(),
            43,
        ),
        (b"I13_CAPACITY_RECEIPT_JOURNAL_WRITE_V3\0".as_slice(), 38),
        (
            b"I13_CAPACITY_MUTATION_PREPARED_PAYLOAD_V3\0".as_slice(),
            42,
        ),
        (
            b"I13_CAPACITY_MUTATION_DURABILITY_AUTHORITY_SET_V3\0".as_slice(),
            50,
        ),
        (
            b"I13_CAPACITY_OPERATION_EVIDENCE_BUNDLE_V3\0".as_slice(),
            42,
        ),
        (
            b"I13_CAPACITY_STATE_UPDATE_PROOF_BUNDLE_V3\0".as_slice(),
            42,
        ),
        (b"I13_COMBINED_RETENTION_ARCHIVE_SET_V3\0".as_slice(), 38),
        (b"I13_CONTROL_ARCHIVE_ANCHOR_PAYLOAD_V3\0".as_slice(), 38),
        (b"I13_CONTROL_ARCHIVE_QUORUM_MEMBER_SET_V3\0".as_slice(), 41),
        (b"I13_CONTROL_ARCHIVE_QUORUM_RECEIPT_V3\0".as_slice(), 38),
        (b"I13_CONTROL_ARCHIVE_SET_V3\0".as_slice(), 27),
        (b"I13_CONTROL_EVIDENCE_CUTOFF_RECEIPT_V3\0".as_slice(), 39),
        (b"I13_CONTROL_EVIDENCE_PLAN_V3\0".as_slice(), 29),
        (b"I13_CONTROL_PHYSICAL_LAYOUT_V3\0".as_slice(), 31),
        (b"I13_CONTROL_PLANE_ARTIFACT_V3\0".as_slice(), 30),
        (
            b"I13_CONTROL_PLANE_LEDGER_APPEND_RECEIPT_V3\0".as_slice(),
            43,
        ),
        (
            b"I13_CONTROL_PLANE_POSTSETTLEMENT_BUNDLE_V3\0".as_slice(),
            43,
        ),
        (
            b"I13_CONTROL_PLANE_PRESETTLEMENT_EVIDENCE_V3\0".as_slice(),
            44,
        ),
        (b"I13_FAILURE_INVENTORY_V3\0".as_slice(), 25),
        (b"I13_GLOBAL_ACTIVE_RESERVATION_SET_V3\0".as_slice(), 37),
        (
            b"I13_GLOBAL_CONTROL_CLOSURE_APPEND_RECEIPT_V3\0".as_slice(),
            45,
        ),
        (
            b"I13_GLOBAL_CONTROL_CLOSURE_DURABILITY_BUNDLE_V3\0".as_slice(),
            48,
        ),
        (b"I13_GLOBAL_CONTROL_CLOSURE_ENTRY_V3\0".as_slice(), 36),
        (
            b"I13_GLOBAL_CONTROL_CLOSURE_JOURNAL_HEAD_V3\0".as_slice(),
            43,
        ),
        (b"I13_GLOBAL_GOVERNOR_CAPACITY_HEAD_V3\0".as_slice(), 37),
        (
            b"I13_GLOBAL_GOVERNOR_TRANSITION_WITNESS_V3\0".as_slice(),
            42,
        ),
        (b"I13_GLOBAL_JOURNAL_EXPIRY_RECEIPT_V3\0".as_slice(), 37),
        (b"I13_LOGICAL_SLOT_SURVIVOR_SUM_V3\0".as_slice(), 33),
        (b"I13_NORMAL_ARCHIVE_CLOSURE_PAYLOAD_V3\0".as_slice(), 38),
        (b"I13_OPERATION_PAYLOAD_ARTIFACT_V3\0".as_slice(), 34),
        (b"I13_OPERATION_PAYLOAD_AUTHORITY_V3\0".as_slice(), 35),
        (b"I13_OPERATION_PAYLOAD_SCHEMA_MATRIX_V3\0".as_slice(), 39),
        (b"I13_PUBLICATION_CERTIFICATE_ARTIFACT_V3\0".as_slice(), 40),
        (b"I13_REPLAY_ARCHIVE_RECEIPT_V3\0".as_slice(), 30),
        (b"I13_REPLAY_METADATA_CAPACITY_HEAD_V3\0".as_slice(), 37),
        (b"I13_REPLAY_METADATA_LEASE_V3\0".as_slice(), 29),
        (b"I13_REPLAY_METADATA_LOCATOR_V3\0".as_slice(), 31),
        (
            b"I13_REPLAY_METADATA_RESERVATION_RECEIPT_V3\0".as_slice(),
            43,
        ),
        (b"I13_RETAINED_STORE_ARTIFACT_SET_V3\0".as_slice(), 35),
        (b"I13_RETENTION_INVENTORY_CLOSURE_V3\0".as_slice(), 35),
        (
            b"I13_RETENTION_RECONSTRUCTION_CERTIFICATE_SET_V3\0".as_slice(),
            48,
        ),
        (
            b"I13_RETENTION_RECONSTRUCTION_CERTIFICATE_V3\0".as_slice(),
            44,
        ),
        (b"I13_RETENTION_REGISTRY_ARTIFACT_V3\0".as_slice(), 35),
        (b"I13_RETENTION_REGISTRY_HEAD_V3\0".as_slice(), 31),
        (b"I13_RETENTION_ROLE_POLICY_MATRIX_V3\0".as_slice(), 36),
        (b"I13_RETENTION_SOURCE_INVENTORY_HEAD_V3\0".as_slice(), 39),
        (b"I13_RETENTION_WRITE_RECEIPT_SET_V3\0".as_slice(), 35),
        (b"I13_TERMINAL_APPEND_DELTA_ARTIFACT_V3\0".as_slice(), 38),
        (b"I13_TERMINAL_DURABILITY_HEAD_V3\0".as_slice(), 32),
        (b"I13_TERMINAL_ROUTE_REQUEST_V2\0".as_slice(), 30),
        (b"I13_CAPACITY_AUTHORITY_BINDING_MATRIX_V4\0".as_slice(), 41),
        (b"I13_CAPACITY_AUTHORITY_ROSTER_MATRIX_V4\0".as_slice(), 40),
        (b"I13_CAPACITY_INDEX_HEAD_V4\0".as_slice(), 27),
        (
            b"I13_CAPACITY_INDEX_RESERVATION_SUBJECT_V4\0".as_slice(),
            42,
        ),
        (b"I13_CAPACITY_INDEX_ALLOCATION_REQUEST_V4\0".as_slice(), 41),
        (b"I13_CAPACITY_INDEX_ALLOCATION_RECEIPT_V4\0".as_slice(), 41),
        (
            b"I13_CAPACITY_INDEX_ALLOCATION_EVIDENCE_BUNDLE_V4\0".as_slice(),
            49,
        ),
        (b"I13_CAPACITY_INDEX_DECISION_HEAD_V4\0".as_slice(), 36),
        (b"I13_CAPACITY_INDEX_DECISION_ENTRY_V4\0".as_slice(), 37),
        (b"I13_CAPACITY_INDEX_LOCAL_DECISION_V4\0".as_slice(), 37),
        (b"I13_CAPACITY_INDEX_COPY_COORDINATOR_V4\0".as_slice(), 39),
        (b"I13_CAPACITY_INDEX_DECISION_APPEND_V4\0".as_slice(), 38),
        (b"I13_CAPACITY_INDEX_ALLOCATION_OUTCOME_V4\0".as_slice(), 41),
        (
            b"I13_CAPACITY_INDEX_ALLOCATION_RECEIPT_PAIR_V4\0".as_slice(),
            46,
        ),
        (
            b"I13_CAPACITY_INDEX_ALLOCATION_AUTHORITY_SET_V4\0".as_slice(),
            47,
        ),
        (b"I13_CAPACITY_INDEX_ORPHAN_HEAD_V4\0".as_slice(), 34),
        (b"I13_CAPACITY_INDEX_ORPHAN_PROOF_V4\0".as_slice(), 35),
        (b"I13_CAPACITY_RESPONSE_LOG_HEAD_V4\0".as_slice(), 34),
        (b"I13_CAPACITY_RESPONSE_LOG_ENTRY_V4\0".as_slice(), 35),
        (
            b"I13_CAPACITY_RESPONSE_LOG_APPEND_RECEIPT_V4\0".as_slice(),
            44,
        ),
        (
            b"I13_CAPACITY_OPERATION_PROTOCOL_MATRIX_V4\0".as_slice(),
            42,
        ),
        (
            b"I13_CAPACITY_TRANSACTION_COORDINATOR_HEAD_V4\0".as_slice(),
            45,
        ),
        (b"I13_CAPACITY_TRANSACTION_WAL_HEAD_V4\0".as_slice(), 37),
        (b"I13_CAPACITY_TRANSACTION_WAL_ENTRY_V4\0".as_slice(), 38),
        (b"I13_CAPACITY_TRANSACTION_CHAIN_STEP_V4\0".as_slice(), 39),
        (b"I13_CAPACITY_TRANSACTION_WAL_RECORD_V4\0".as_slice(), 39),
        (b"I13_CAPACITY_TRANSACTION_UPDATE_PROOF_V4\0".as_slice(), 41),
        (b"I13_CAPACITY_MUTATION_COMMIT_INTENT_V4\0".as_slice(), 39),
        (
            b"I13_CAPACITY_MUTATION_PREPARED_PAYLOAD_V4\0".as_slice(),
            42,
        ),
        (b"I13_CAPACITY_TRANSACTION_RECEIPT_V4\0".as_slice(), 36),
        (
            b"I13_CAPACITY_TRANSACTION_WAL_APPEND_RECEIPT_V4\0".as_slice(),
            47,
        ),
        (
            b"I13_CAPACITY_TRANSACTION_DURABILITY_SET_V4\0".as_slice(),
            43,
        ),
        (
            b"I13_CAPACITY_MUTATION_PRECOMMIT_MANIFEST_V4\0".as_slice(),
            44,
        ),
        (
            b"I13_CAPACITY_MUTATION_EVIDENCE_MANIFEST_V4\0".as_slice(),
            43,
        ),
        (
            b"I13_CAPACITY_MUTATION_EVIDENCE_ARCHIVE_V4\0".as_slice(),
            42,
        ),
        (
            b"I13_CONTROL_EVIDENCE_PRODUCER_INVENTORY_V4\0".as_slice(),
            43,
        ),
        (b"I13_CONTROL_EVIDENCE_STAGING_HEAD_V4\0".as_slice(), 37),
        (
            b"I13_CONTROL_EVIDENCE_STAGING_PROOF_BUNDLE_V4\0".as_slice(),
            45,
        ),
        (
            b"I13_CONTROL_EVIDENCE_STAGING_INSERTION_V4\0".as_slice(),
            42,
        ),
        (
            b"I13_CONTROL_EVIDENCE_STAGING_SEAL_PREPARED_V4\0".as_slice(),
            46,
        ),
        (b"I13_CONTROL_EVIDENCE_STAGING_ARTIFACT_V4\0".as_slice(), 41),
        (b"I13_CONTROL_EVIDENCE_STAGING_LOG_PAGE_V4\0".as_slice(), 41),
        (b"I13_CONTROL_EVIDENCE_CUTOFF_V4\0".as_slice(), 31),
        (b"I13_GLOBAL_CLOSED_MEMBER_V4\0".as_slice(), 28),
        (b"I13_GLOBAL_CLOSED_EVIDENCE_MAP_HEAD_V4\0".as_slice(), 39),
        (b"I13_GLOBAL_JOURNAL_EXPIRY_REQUEST_V4\0".as_slice(), 37),
        (
            b"I13_GLOBAL_JOURNAL_EXPIRY_PROOF_BUNDLE_V4\0".as_slice(),
            42,
        ),
        (b"I13_GLOBAL_CLOSURE_EVIDENCE_ARTIFACT_V4\0".as_slice(), 40),
    ] {
        assert_eq!(bytes.len(), expected);
    }

    assert_eq!(44 + 2 * 8 + 32 + 1 + 8 * 8 + 4 * 32 + 8 + 2 * 32, 357);
    assert_eq!(
        45 + 2 * 32 + 11 * 8 + 2 * 32 + 2 * 8 + 1 + 2 * 8 + 32 + 5 * 32 + 2 * 8,
        502
    );
    assert_eq!(41 + 2 * 8 + 1 + 3 * 32 + 5 * 32 + 12 * 8, 410);
    assert_eq!(37 + 2 * 8 + 3 * 32 + 2 * 8 + 1 + 2 * 8 + 32, 214);
    assert_eq!(
        67_108_864u64 + 2 * 201_326_592 + 4 * 269_409_356,
        1_547_399_472
    );
    assert_eq!(1_547_399_472u64 + 4 * 4_194_304, 1_564_176_688);
    assert_eq!(2_147_483_648u64 - 1_564_176_688, 583_306_960);
    assert_eq!(1_547_399_472u64 + 2 * 71 * 4_194_304, 2_142_990_640);
    assert_eq!(1_547_399_472u64 + 2 * 72 * 4_194_304, 2_151_379_248);
    assert_eq!(111 + 257 * (8 + 229 + 32 + 8 + 8 + 1_048_000), 269_409_356);
    assert_eq!(
        28 + 2 * 8 + 5 * 32 + 4 * 8 + 76 + 32 + 8 + 4 * 40 + 32 + 2 * 8 + 32,
        592
    );
    assert_eq!(34 + 1 + 8 + 592, 635);
    assert_eq!(535 + 32, 567);
    assert_eq!(34 + 1 + 8 + 567, 610);
    assert_eq!(37 + 3 * 8 + 3 * 32, 157);
    assert_eq!(157 + 1_024 * 354, 362_653);
    assert!(354 * 180_225 + 157 * 177 <= 92_275_200);
    assert_eq!(40 + 1 + 8 * 32 + 9 * 8 + 1, 370);
    assert_eq!(
        35 + 2 * 8 + 5 * 32 + 2 * 8 + 3 * 32 + 1 + 2 * 8 + 2 * 32,
        404
    );
    assert_eq!(30 + 7 * 32 + 3 * 8, 278);
    assert_eq!(16 + 8 + 133 + 256 * 32, 8_349);
    assert_eq!(31 + 32 + 8 + 8 + 32 + 8 + 8, 127);
    assert_eq!(33 + 2 * 32 + 2 + 2 * 32, 163);
    assert_eq!(24 + 1 + 3 * 32 + 8, 129);
    assert_eq!(34 + 2 * 8 + 6 * 32, 242);
    assert_eq!(41 + 2 * 32 + 2, 107);
    assert_eq!(2 * 8 + 2 * 32, 80);
    assert_eq!(7 * 32 + 8 + 3 * 8, 256);
    assert_eq!(43 + 9 * 32 + 2 + 2 + 1, 336);
    assert_eq!(36 + 4 * 32 + 2 + 2 * 32, 230);
    assert_eq!(34 + 2 * 32 + 3 * 8 + 1_024 * (8 + 8_349), 8_557_690);
    assert_eq!(37 + 2 * 32 + 3 * 8, 125);
    assert_eq!(32 + 8 * 8 + 1 + 3 * 32, 193);
    assert_eq!(126 + 180_225 * 173, 31_179_051);
    assert_eq!(63_827_439 + 67_108_864 + 31_179_051, 162_115_354);
    assert_eq!(201_326_592 - 162_115_354, 39_211_238);
    assert_eq!(31 + 5 * 8 + 1 + 11 * 32 + 2 * 127, 678);
    assert_eq!(34 + 7 * 8 + 13 * 32 + 127 + 1, 634);
    assert_eq!(43 + 7 * 8 + 2 + 2 + 8 * 32 + 1, 360);
    assert_eq!(42 + 5 * 8 + 2 + 2 + 8 * 32, 342);
    assert_eq!(37 + 2 * 8 + 2 + 2 + 4 * 32, 185);
    assert_eq!(44 + 11 * 8 + 2 + 12 * 32, 518);
    assert_eq!(518 + 2 * 32, 582);
    assert_eq!(38 + 4 * 8 + 7 * 32, 294);
    assert_eq!(34 + 1 + 8 + 567, 610);
    assert_eq!(34, 32 + 2);
    assert_eq!(98, 32 + 2 + 32 + 32);
    assert_eq!(90, 32 + 3 * 8 + 2 + 32);
    assert_eq!(522, 8 + 2 + 16 * 32);
    assert_eq!(33 + 7 * 32 + 8, 265);
    assert_eq!(38 + 65_536 * 48, 3_145_766);
    assert_eq!(37 + 2 * 32 + 3 * 8, 125);
    assert_eq!(41 + 8 * 32 + 9 * 8 + 2, 371);
    assert_eq!(36 + 2 * 8 + 32 + 2 + 2 * (1 + 1 + 32), 154);
    assert_eq!(38 + 2 * 8 + 32 + 2 + 2 * (1 + 32), 154);
    assert_eq!(410 + 32 + 8, 450);
    assert_eq!(35 + 1 + 32 + 8 + 8, 84);
    assert_eq!(31 + 2 * 8 + 2 + 8 + 1 + 3 * 32, 154);
    assert_eq!(38 + 8 * 32, 294);
    assert_eq!(32 + 3 * 8 + 5 * 32, 216);
    assert_eq!(35 + 3 + 6 * 8 + 8 * 32, 342);
    assert_eq!(40 + 2 + 6 * 8 + 10 * 32, 410);
    assert_eq!(
        2 * 678 + 2 * 370 + 2 * 240 + 404 + 634 + 216 + 2 * 342 + 410,
        4_924
    );
    assert_eq!(4_924 * 4_096, 20_168_704);
    assert_eq!(162_115_354 + 20_168_704, 182_284_058);
    assert_eq!(201_326_592 - 182_284_058, 19_042_534);
    assert_eq!(37 + 6 * 8 + 8 * 32 + 1, 342);
    assert_eq!(40 + 9 * 8 + 9 * 32 + 3, 403);
    assert_eq!(4 * 32 + 8, 136);
    assert_eq!(36 + 4 * 8 + 10 * 32 + 1, 389);
    assert_eq!(1 + 2 + 32 + 8, 43);
    assert_eq!(29 + 4 * 8, 61);
    assert_eq!(1 + 2 * 32 + 8 + 32, 105);
    assert_eq!(35 + 2 * 8 + 2, 53);
    assert_eq!(39 + 6 * 8 + 9 * 32, 375);
    assert_eq!(49 + 5 * 8 + 13 * 32, 505);
    assert_eq!(42 + 7 * 8 + 12 * 32 + 2, 484);
    assert_eq!(30 + 2 * 8 + 2, 48);
    assert_eq!(2 + 2 * 32, 66);
    assert_eq!(25 + 2 * 8 + 1 + 2, 44);
    assert_eq!(2 + 2 + 32 + 1, 37);
    assert_eq!(27 + 2 * 8 + 32 + 2, 77);
    assert_eq!(28 + 2 * 32 + 3 * 8, 116);
    assert_eq!(31 + 7 * 8 + 8 * 32 + 2, 345);
    assert_eq!(345 + 32, 377);
    assert_eq!(47 + 1 + 2 * (2 * 8) + 7 * 32 + 1, 305);
    assert_eq!(43 + 32 + 8, 83);
    assert_eq!(39 + 32 + 2 + 6 * (1 + 32), 271);
    assert_eq!(450 + 32, 482);
    assert_eq!(32 + 2 * 8 + 32 + 3 * 8, 104);
    assert_eq!(29 + 5 * 8, 69);
    assert_eq!(69 * 4_096 + 32 * (180_225 + 180_225), 11_817_024);
    assert_eq!(162_115_354 + 20_168_704 + 11_817_024, 194_101_082);
    assert_eq!(201_326_592 - 194_101_082, 7_225_510);
    assert_eq!(28 + 32 + 8, 68);
    assert_eq!(125 + 8 + 32, 165);
    assert_eq!(29 + 3 * 8 + 32 + 8, 93);
    assert_eq!(37 + 6 * 8 + 7 * 32 + 2, 311);
    assert_eq!(80 + 65_534 * 64, 4_194_256);
    assert_eq!(32 + 8, 40);
    assert_eq!(31 + 32 + 8, 71);
    assert_eq!(116 + 2 * 8, 132);
    assert_eq!(71 + 8 + 40, 119);
    assert_eq!(27 + 4 * 8 + 1 + 5 * 32, 220);
    assert_eq!(41 + 4 * 8 + 1 + 10 * 32, 394);
    assert_eq!(8_557_690 + 32, 8_557_722);
    assert_eq!(48 + 21 * 26 + 460, 1_054);
    assert_eq!(39 + 32 + 2 + 7 * (1 + 32), 304);
    assert_eq!(482 + 2 * 32, 546);
    assert_eq!(34 + 2 * 8 + 32 + 2 * 8, 98);
    assert_eq!(38 + 2 * 8 + 32 + 2 * 8, 102);
    assert_eq!(102u64 + 4_096 * 8 + 11_817_024, 11_849_894);
    assert_eq!(32 + 2 * 8 + 32 + 2 * 8, 96);
    assert_eq!(8 + 32 + 8 + 4_926, 4_974);
    assert_eq!(96u64 + 4_096 * 4_974, 20_373_600);
    assert_eq!(216 + 2, 218);
    assert_eq!(4_924 + 2, 4_926);
    assert_eq!(
        162_115_354u64 + 20_373_600 + 11_849_894 + 67_108_864,
        261_447_712
    );
    assert_eq!(268_435_456u64 - 261_447_712, 6_987_744);
    assert_eq!(6u64 * (4_194_304 + 48), 25_166_112);
    assert_eq!(5u64 * (1_048_576 + 1_024 * 80) + 7_464_912, 13_117_392);
    assert_eq!(5u64 * 362_653 + 997_689, 2_810_954);
    assert_eq!(5u64 * 65_613 + 122_989, 451_054);
    assert_eq!(6u64 * 1_024 * 173, 1_062_912);
    assert_eq!(6u64 * 4_974, 29_844);
    assert_eq!(4_194_304u64 + 48, 4_194_352);
    assert_eq!(7_239_552u64 + 2_817 * 80, 7_464_912);
    assert_eq!(122_981u64 + 8, 122_989);
    assert_eq!(26 + 2 * 8 + 1 + 2 * 2 + 1 + 6 * 32 + 8 + 4 * 32, 376);
    assert_eq!(34 + 2 * 8 + 1 + 2 * 2 + 1 + 6 * 8 + 8 + 32 + 8 + 1, 153);
    assert_eq!(25 + 2 * 8 + 1 + 8 + 1 + 2 * 32 + 2 * 8 + 32, 163);
    assert_eq!(
        34 + 2 * 8 + 2 + 2 * 32 + 2 * 8 + 2 * 32 + 8 + 3 * 32 + 1 + 2 * 8,
        317
    );
    assert_eq!(37 + 2 * 8 + 1 + 32 + 2 * 8, 102);
    assert_eq!(102u64 + 16_384 * (8 + 317), 5_324_902);
    assert_eq!(34 + 2 * 8 + 1 + 32 + 2 * 8 + 2 + 1, 102);
    assert_eq!(102u64 + 16_385 * (8 + 163), 2_801_937);
    assert_eq!(35 + 2 * 8 + 1 + 32 + 2 * 8, 100);
    assert_eq!(100u64 + 4_096 * (8 + 376), 1_572_964);
    assert_eq!(34 + 2 * 8 + 1 + 32 + 2 * 8, 99);
    assert_eq!(99u64 + 4_096 * (8 + 153), 659_555);
    assert_eq!(34 + 2 * 8 + 4 * 32, 178);
    assert_eq!(26 + 2 * 8 + 32 + 8 + 2 * 2 + 3 + 512 + 32, 633);
    assert_eq!(45 + 2 * 8 + 32 + 1 + 2 * 2 + 5 * 32, 258);
    assert_eq!(42 + 32 + 2 * 8 + 8 + 32 + 2 + 32, 164);
    assert_eq!(46 + 32 + 2 * 32 + 2 * 8 + 2 * 32 + 3 * 32 + 1 + 2 * 8, 335);
    assert_eq!(35 + 2 * 8 + 1 + 8 + 8 * 32 + 4 * 8, 348);
    assert_eq!(
        5_324_902u64 + 2_801_937 + 1_572_964 + 659_555 + 126_402,
        10_485_760
    );
    assert_eq!(35 + 2 * 8 + 1 + 32 + 2 + 2 * (8 + 138), 378);
    assert_eq!(36 + 2 * 8 + 1 + 2 + 2 * (8 + 73), 217);
    assert_eq!(34 + 2 * 8 + 2 + 9 * 32 + 2 * 8, 356);
    assert_eq!(403 + 2 * 32, 467);
    assert_eq!(261_447_712u64 + 6_987_744, 268_435_456);
    assert_eq!(35 + 32 + 2 * 8 + 2 + 8 * (8 + 265), 2_269);
    assert_eq!(43 + 32 + 2 * 8 + 3 * 8 + 10 * 32, 435);
    assert_eq!(357 + 32, 389);
    assert_eq!(502 + 32, 534);
    assert_eq!(193 + 32, 225);
    assert_eq!(34 + 2 * 8 + 32 + 2 * 8, 98);
    assert_eq!(34 + 32 + 6 * 8 + 3 * 32 + 8, 218);
    assert_eq!(2 * 32 + 1 + 7 * 8, 121);
    assert_eq!(3 * 32 + 8 * 8, 160);
    assert_eq!(2 * 8 + 2 * 8 + 8 + 160 + 256 * 32, 8_392);
    assert_eq!(2 * 8 + 2 * 8 + 8 + 72 + 256 * 32, 8_304);
    assert_eq!(
        44 + 2 * 8 + 1 + 32 + 1 + 7 * 32 + 2 * 8 + 3 * 32 + 1 + 2 * 8,
        447
    );
    assert_eq!(
        4u64 * 269_409_356
            + 2 * 268_435_456
            + 2 * 134_217_728
            + 2 * 10_485_760
            + 67_108_864
            + 67_108_864
            + 8_388_608,
        2_046_521_648
    );
    assert_eq!(
        268_435_456u64 + 134_217_728 + 2 * 269_409_356 + 10_485_760,
        951_957_656
    );
    assert_eq!(2_046_521_648u64 + 2 * 8_388_608, 2_063_298_864);
    assert_eq!(2_147_483_648u64 - 2_063_298_864, 84_184_784);
    assert_eq!(2_046_521_648u64 + 12 * 8_388_608, 2_147_184_944);
    assert_eq!(2_147_483_648u64 - 2_147_184_944, 298_704);
    assert_eq!(2_046_521_648u64 + 13 * 8_388_608, 2_155_573_552);
    assert_eq!(37 + 32 + 2 + 11 * (8 + 1 + 2 + 32) + 2 * 99, 742);
    assert_eq!(38 + 2 * 32 + 2 + 4 * (8 + 2 + 4 * 32), 656);
    assert_eq!(27 + 2 * 8 + 2 * 2 + 11 * 32, 399);
    assert_eq!(40 + 2 * 8 + 4 * 32 + 2 + 4 * (8 + 399), 1_814);
    assert_eq!(
        32 + 2 * 8 + 2 * 32 + 2 * 2 + 1 + 1 + 2 + 1 + 12 * 32 + 2 * 8,
        521
    );
    assert_eq!(35 + 2 * 8 + 2 * 32 + 2 * 2 + 5 * 32 + 1 + 2 * 8, 296);
    assert_eq!(38 + 2 * 8 + 2 * 32 + 2 * 2 + 2 + 11 * 32 + 2 * 8, 492);
    assert_eq!(77 + 4 * (8 + 69), 385);
    assert_eq!(505 + 3 * 32, 601);
    assert_eq!(27 + 32 + 4 * 8 + 32, 123);
    assert_eq!(32 + 2 * 8 + 1 + 3 * 32 + 2 * 8 + 2 * 32, 225);
    assert_eq!(225 + 2 * 32, 289);
    assert_eq!(27 + 32 + 7 * 8 + 3 * 32, 211);
    assert_eq!(4 * 8 + 2 * (8 + 289) + 256 * 32, 8_818);
    assert_eq!(42 + 2 * 8 + 1 + 1 + 8, 68);
    assert_eq!(68 + (8 + 8_248) + (8 + 8_818), 17_150);
    assert_eq!(37 + 2 * 8 + 1 + 11 * 32 + 1 + 2 * 8, 423);
    assert_eq!(34 + 2 * 8 + 1 + 2 * 32 + 8 + 8, 131);
    assert_eq!(33 + 32 + 3 * 8 + 32, 121);
    assert_eq!(39 + 2 * 8 + 4 * 32 + 2 * 8 + 3 * 32 + 1 + 2 * 8, 312);
    assert_eq!(35 + 32 + 4 * 8 + 32, 131);
    assert_eq!(30 + 2 * 8 + 1, 47);
    assert_eq!(1 + 32 + 8 + 3 * 32, 137);
    assert_eq!(47 + 2 * (8 + 137), 337);
    assert_eq!(37 + 2 * 8 + 1 + 2, 56);
    assert_eq!(56 + 13 * 57, 797);
    assert_eq!(
        38 + 2 * 8 + 1 + 32 + 1 + 3 * 32 + 2 * 8 + 4 * 32 + 1 + 2 * 8,
        345
    );
    assert_eq!(42 + 2 * 8 + 1 + 32 + 2 + 2 * 8, 109);
    assert_eq!(109u64 + 65_535 * (8 + 32), 2_621_509);
    assert_eq!(37 + 2 * 8 + 1 + 4 * 32 + 2 + 8 + 3 * 32 + 1 + 2 * 8, 305);
    assert_eq!(29 + 1 + 4 * 8, 62);
    assert_eq!(62u64 + 65_535 * 51, 3_342_347);

    for exact in [
        "TERMINAL_CAPACITY_RESERVATION_AND_ACTIVATION_DAG_V2",
        "request->reservation intent->reservation receipt->admission core->activation receipt->context DAG",
        "TERMINAL_CAPACITY_GOVERNOR_IDENTITY_AND_GENERAL_EQUATION_V2",
        "GlobalTerminalGovernorHeadBytesV2 is the 36 exact ASCII bytes I13_GLOBAL_TERMINAL_GOVERNOR_HEAD_V2",
        "single global ledger",
        "p_i is in 2..=71",
        "global_governor_bytes=2147483648",
        "exactly 357 bytes",
        "exactly 502 bytes",
        "exactly 410 bytes",
        "exactly 214 bytes",
        "1547399472+8388608*p_i",
        "p_i=2 reserves 1564176688 and leaves 583306960",
        "p_i=71 reserves 2142990640 and leaves 4493008",
        "p_i=72 requires 2151379248 and is refused",
        "earlier 325-byte TerminalAuthorityAdmissionBytesV2",
        "non-authoritative superseded draft",
        "V2_FRONTIER_CONTEXT_AND_SELECTION_AUTHORITY",
        "CandidateFrontierContextV2=raw32(terminal_authority_context_v2)",
        "WitnessFrontierContextV2=raw32(terminal_authority_context_v2)",
        "BreachFactFrontierContextV2=raw32(terminal_authority_context_v2)",
        "raw32(terminal_authority_context_v2) is an exact prefix",
        "ArbitratedSelectionReceiptBytesV2 is the 27 exact ASCII bytes I13_ARBITRATED_SELECTION_V2",
        "exactly 592 bytes",
        "TerminalSelectionReceiptBytesV2 is the 33 exact ASCII bytes I13_TERMINAL_SELECTION_RECEIPT_V2",
        "ArbitratedSelectionReceiptBytesV2 remains 592 bytes and its union arm remains 635",
        "CapacityBreach arm is therefore exactly 610 bytes",
        "501-byte receipt, TerminalSelectionReceiptV1",
        "decode/replay-only",
        "BREACH_FACT_ARTIFACT_AUTHORITY_V1",
        "complete artifact is at most 111+257*(8+229+32+8+8+1048000)=269409356 bytes",
        "exactly four 269409356-byte reservations",
        "TERMINAL_EVENT_ARTIFACT_AND_APPEND_RECEIPT_V2",
        "TerminalEventArtifactRootV2=BLAKE3::derive_key('org.frankensim.i13.terminal-event-artifact.v2'",
        "TERMINAL_DURABILITY_ARTIFACT_BUNDLE_V2",
        "TerminalEventPayloadArtifactBytesV2 is the 38 exact ASCII bytes I13_TERMINAL_EVENT_PAYLOAD_ARTIFACT_V2",
        "ReceiptIndexArtifactBytesV2 is the 29 exact ASCII bytes I13_RECEIPT_INDEX_ARTIFACT_V2",
        "totaling 162115354 and leaving 39211238 bytes",
        "DurableAppendOperationReceiptBytesV1 is the 39 exact ASCII bytes I13_DURABLE_APPEND_OPERATION_RECEIPT_V1",
        "exactly 370 bytes",
        "PrefixVerificationReceiptBytesV1 is the 34 exact ASCII bytes I13_PREFIX_VERIFICATION_RECEIPT_V1",
        "exactly 404 bytes",
        "historical certificate field named independent_prefix_verifier_capability_digest byte-equals this receipt digest",
        "TERMINAL_OPERATIONAL_STORAGE_AND_RECOVERY_V2",
        "expected_witness_frontier_head_digest normatively means expected_witness_frontier_root",
        "GOVERNED_SUBARTIFACT_AND_COMPLETENESS_AUTHORITY_V2",
        "SEMANTIC_INVENTORY_EVALUATOR_AND_KAT_CLOSURE_V2",
        "SemanticInventoryDigestV2=BLAKE3::derive_key('org.frankensim.i13.semantic-inventory.v2'",
        "whose normalized root is True is refused",
        "TerminalRecordKeyList=2,AuthorityReceiptDigest=3 or BreachFactFrontierReference=4",
        "tag2 resolution contributes the exact producer_terminal_record_key to CausalParentSetRoot",
        "IndependentConformanceEvaluatorSetBytesV2 is the 32 exact ASCII bytes I13_CONFORMANCE_EVALUATOR_SET_V2",
        "GovernedKatSetBytesV2 is the 23 exact ASCII bytes I13_GOVERNED_KAT_SET_V2",
        "GeneratorExecutableIdentity payload is those 96 bytes plus raw32(sandbox_capability_profile_digest), exactly 128 bytes",
        "CaseOrderRule payload is raw32(SemanticInventoryDigestV2)",
        "exactly 81 bytes",
        "LabelSealingPolicy payload is raw32(SemanticInventoryDigestV2)",
        "exactly 163 bytes",
        "SemanticCompletenessMatrixRootV2=BLAKE3::derive_key('org.frankensim.i13.semantic-completeness-matrix.v2'",
        "SEMANTIC_COMPLETENESS_ALL_AND_ONLY_GENERATOR_AUTHORITY_V2",
        "scope_kind is the closed enum ExternalSemanticGroup=1,TerminalEventKind=2,HoldoutTargetSlot=3",
        "SemanticCatalogGenerationReceiptSetBytesV2 is the 35 exact ASCII bytes I13_SEMANTIC_CATALOG_RECEIPT_SET_V2",
        "HOLDOUT_COMMITMENT_SCHEMA_RECEIPT_CLOSURE_V2",
        "generator_digest in the historical HOLDOUT_REALIZATION commitment normatively byte-equals holdout_generator_contract_digest",
        "case_order_digest byte-equals case_order_rule_digest",
        "HoldoutPayloadSchemaDigestV2=BLAKE3::derive_key('org.frankensim.i13.holdout-payload-schema.v2'",
        "HOLDOUT_COMMITMENT_SCHEMA_AND_REALIZED_MERKLE_AUTHORITY_V2",
        "HoldoutRealizationCommitmentBytesV2 is the 37 exact ASCII bytes I13_HOLDOUT_REALIZATION_COMMITMENT_V2",
        "RealizedMerkleRootV2=BLAKE3::derive_key('org.frankensim.i13.holdout-realized-merkle.root.v2'",
        "TERMINAL_EVENT_SEMANTIC_AUTHORITY_AND_RECEIPT_INDEX_V2",
        "TerminalEventSemanticAuthorityRootV2=BLAKE3::derive_key('org.frankensim.i13.terminal-event-semantic-authority.v2'",
        "ReceiptIndexUpdateProofBytesV2",
        "at most 8349 bytes",
        "fixed and excludes mutable receipt-index root/count",
        "TERMINAL_AUTHORITY_JOINT_HEAD_SIZE_CORRECTION_V2",
        "exact 127-byte projection",
        "TERMINAL_CAUSAL_SCHEMA_AND_EVENT_AUTHORITY_V2",
        "TerminalCausalSchemaRootV2=BLAKE3::derive_key('org.frankensim.i13.terminal-causal-schema.v2'",
        "historical TERMINAL_EVENT_CHAIN_V1 definition",
        "TERMINAL_EVENT_DIGEST_CAUSAL_CONTEXT_AND_FRONTIER_CORRECTION_V2",
        "TerminalEventDigestProjectionV2 is exactly raw32(terminal_authority_context_v2)",
        "TerminalFrontierContextV2 is exactly raw32(terminal_authority_context_v2)",
        "TERMINAL_APPEND_PREPARE_MIRROR_AND_PUBLICATION_AUTHORITY_V2",
        "TerminalAppendRequestBytesV2 is the 30 exact ASCII bytes I13_TERMINAL_APPEND_REQUEST_V2",
        "exactly 678 bytes",
        "TerminalDurabilityCommitBytesV2 is the 33 exact ASCII bytes I13_TERMINAL_DURABILITY_COMMIT_V2",
        "exactly 634 bytes",
        "independently valid certificate is mandatory before visibility",
        "TERMINAL_SELECTION_AND_BREACH_DURABILITY_CLOSURE_V2",
        "BreachFactDurabilityBundleBytesV2 is the 36 exact ASCII bytes I13_BREACH_FACT_DURABILITY_BUNDLE_V2",
        "TERMINAL_CAPACITY_SETTLEMENT_AUTHORITY_V2",
        "TerminalCapacitySettlementReceiptBytesV2 is the 43 exact ASCII bytes I13_TERMINAL_CAPACITY_SETTLEMENT_RECEIPT_V2",
        "exactly 518 bytes",
        "I13_FRESH_V2_AUTHORITY_PRECEDENCE",
        "TERMINAL_DURABILITY_CAPACITY_FIELD_EQUALITIES_V2",
        "I13_FRESH_V2_DOUBLE_UTF8_FRAME_SIZE_CORRECTION",
        "SEMANTIC_KAT_NORMALIZATION_AND_MATRIX_CLOSURE_V2",
        "PositiveKatCasePayloadV2 is exactly raw32(SemanticInventoryDigestV2)",
        "SemanticSourceCatalogBytesV2 is the 30 exact ASCII bytes I13_SEMANTIC_SOURCE_CATALOG_V2",
        "exactly 336 bytes",
        "exactly 230 bytes",
        "SEMANTIC_SOURCE_OCCURRENCE_AND_ROLE_CLOSURE_V2",
        "SemanticSourceOccurrenceSetBytesV2 is the 37 exact ASCII bytes I13_SEMANTIC_SOURCE_OCCURRENCE_SET_V2",
        "all repeated lifecycle-event occurrences",
        "An Ordinary non-scalar atom has no scalar-kind field",
        "GLOBAL_GOVERNOR_AND_LANE_ISSUANCE_CLOSURE_V2",
        "GlobalGovernorAuthorityBytesV2 is the 32 exact ASCII bytes I13_GLOBAL_GOVERNOR_AUTHORITY_V2",
        "ActiveReservationSetBytesV2 is the 29 exact ASCII bytes I13_ACTIVE_RESERVATION_SET_V2",
        "checked_sum(member.run_reserve_bytes)=GlobalTerminalGovernorHeadBytesV2.global_reserved_bytes",
        "LaneReservationIssuanceReceiptBytesV2 is the 40 exact ASCII bytes I13_LANE_RESERVATION_ISSUANCE_RECEIPT_V2",
        "LaneReservationAttemptSetBytesV2 is the 35 exact ASCII bytes I13_LANE_RESERVATION_ATTEMPT_SET_V2",
        "Corrected TerminalAuthorityAdmissionCoreBytesV2 appends raw32(LaneReservationAuthoritySetRootV2)",
        "exactly 450 bytes",
        "ARTIFACT_CONTEXT_OPERATION_AND_BATCH_CLOSURE_V2",
        "ReceiptIndexContextDigestV2=BLAKE3::derive_key('org.frankensim.i13.receipt-index-context.v2'",
        "ReceiptIndexUpdateBatchDigestV2=BLAKE3::derive_key('org.frankensim.i13.receipt-index-update-batch.v2'",
        "Each canonical event may contain at most 1048576 bytes",
        "TerminalEventArtifactChunkBytesV2 carries no lane field",
        "OperationPayloadAuthorityBytesV2 is the 34 exact ASCII bytes I13_OPERATION_PAYLOAD_AUTHORITY_V2",
        "BREACH_PERSISTENCE_HEAD_AND_SLOT_CLOSURE_V2",
        "BreachPersistenceHeadBytesV2 is the 30 exact ASCII bytes I13_BREACH_PERSISTENCE_HEAD_V2",
        "Each attempt/lane/phase has a distinct genesis",
        "TERMINAL_PUBLICATION_AND_METADATA_CLOSURE_V2",
        "TerminalPublicationAuthorityBytesV2 is the 37 exact ASCII bytes I13_TERMINAL_PUBLICATION_AUTHORITY_V2",
        "TerminalDurabilityHeadBytesV2 is the 31 exact ASCII bytes I13_TERMINAL_DURABILITY_HEAD_V2",
        "PublicationJournalReceiptBytesV2 is the 34 exact ASCII bytes I13_PUBLICATION_JOURNAL_RECEIPT_V2",
        "TerminalPublicationCertificateBytesV2 is the 39 exact ASCII bytes I13_TERMINAL_PUBLICATION_CERTIFICATE_V2",
        "MAX_TERMINAL_APPEND_COUNT_V2=4096",
        "exactly 4924 bytes per append",
        "leaving 7225510 bytes",
        "SETTLEMENT_EVIDENCE_RETENTION_AND_RELEASE_CLOSURE_V2",
        "CapabilityQuiescenceReceiptBytesV2 is the 36 exact ASCII bytes I13_CAPABILITY_QUIESCENCE_RECEIPT_V2",
        "LaneReservationReleaseReceiptBytesV2 is the 39 exact ASCII bytes I13_LANE_RESERVATION_RELEASE_RECEIPT_V2",
        "PreActivationAbortEvidenceBytesV2 is the 35 exact ASCII bytes I13_PREACTIVATION_ABORT_EVIDENCE_V2",
        "CanonicalFinalizationEvidenceBytesV2 is the 38 exact ASCII bytes I13_CANONICAL_FINALIZATION_EVIDENCE_V2",
        "TerminalPersistenceCatastropheEvidenceBytesV2 is the 48 exact ASCII bytes I13_TERMINAL_PERSISTENCE_CATASTROPHE_EVIDENCE_V2",
        "DurableRetentionTransferReceiptBytesV2 is the 41 exact ASCII bytes I13_DURABLE_RETENTION_TRANSFER_RECEIPT_V2",
        "FINAL_STATE_MACHINE_AND_ENUMERATION_CORRECTIONS_V2",
        "OperationPayloadAuthorityBytesV2.operation_kind is normatively CandidateSelection=1,CapacityBreachSelection=2,BarrierBatch=3,Finalization=4,Failover=5 or Recovery=6",
        "RevokedCapabilitySetBytesV2 is the 29 exact ASCII bytes I13_REVOKED_CAPABILITY_SET_V2",
        "the historically named expected_old_published_authority_digest_or_zero field normatively byte-equals the fetched predecessor TerminalDurabilityHeadDigestV2",
        "FailureInventoryBytesV2 is the 24 exact ASCII bytes I13_FAILURE_INVENTORY_V2",
        "ReplayRouteLedgerBytesV2 is the 26 exact ASCII bytes I13_REPLAY_ROUTE_LEDGER_V2",
        "RetentionStoreHeadBytesV2 is the 27 exact ASCII bytes I13_RETENTION_STORE_HEAD_V2",
        "RetentionStoreReceiptBytesV2 is the 30 exact ASCII bytes I13_RETENTION_STORE_RECEIPT_V2",
        "FINAL_SEMANTIC_CONSTRAINT_AND_SOURCE_MEET_CLOSURE_V2",
        "source_authority_digest_v2=BLAKE3::derive_key('org.frankensim.i13.semantic-source-authority.v2'",
        "Scope2 is exactly the set union of every structured ManifestDraft.obligations[*].obs_events UTF-8 value",
        "SemanticOccurrenceConformanceReceiptBytesV2 is the 46 exact ASCII bytes I13_SEMANTIC_OCCURRENCE_CONFORMANCE_RECEIPT_V2",
        "receipt_count=checked_mul(2,occurrence_count)",
        "target_atom_ordinal in 1..=atom_count",
        "TerminalRecordKeyList uses role byte 9 plus exactly Sequence(min=0,max=64,order=LexicalUnique,child=Scalar(kind=8))",
        "CanonicalScalarKindV2 extends V1 only with kinds 8 and 10",
        "FINAL_STRUCTURED_PROTOCOL_EVENT_REGISTRY_CLOSURE_V3",
        "`I13_PROTOCOL_FIXED_TERMINAL_EVENTS_V2` is the only structured protocol-owned event-kind registry",
        "294+21=315 occurrences and 64 unique exact UTF-8 event keys",
        "cancellation.proposed, terminal.causal_barrier and terminal.post_selection_capacity_breach",
        "FINAL_DURABLE_ARTIFACT_PHASE_AND_CAPACITY_CLOSURE_V3",
        "OrdinaryEventBatch=1,CandidateSelection=2,CapacityBreachSelection=3,BarrierBatch=4,Finalization=5,Failover=6 or Recovery=7",
        "OperationPayloadSchemaMatrixBytesV3",
        "exactly 304 bytes",
        "OperationPayloadArtifactBytesV3",
        "TerminalAppendDeltaArtifactBytesV3",
        "AppendMetadataArtifactBytesV3",
        "exactly 218 bytes",
        "leaving exactly 6987744 nonborrowable safety-slack bytes",
        "Open admits at most 4090 nonclosing OrdinaryEventBatch appends",
        "preserves at least 25166112 operation-payload-artifact bytes",
        "AppendPhaseLeaseBytesV3",
        "AppendPhaseHeadBytesV3",
        "exactly 317 bytes",
        "InfrastructureCapacityReservationSetBytesV3",
        "InfrastructureCapacityReservationSetBytesV3 is the 34 exact ASCII bytes I13_INFRASTRUCTURE_CAPACITY_SET_V3 followed by one byte 0x00, raw32(predecessor_manifest_digest), U64LE(attempt_id), U64LE(epoch), U16LE(member_count=9), and nine FRAME_BYTES(one byte service_kind||raw32(service_identity)||raw32(owner_identity)||raw32(implementation_identity)||raw32(durable_medium_identity)||raw32(failure_domain_identity)||U64LE(reserved_bytes)||raw32(capacity_reservation_receipt_digest)||raw32(retention_policy_digest_or_zero)||raw32(destination_roster_digest_or_zero)) rows in service-kind order, exactly 2542 bytes",
        "exactly 546 bytes",
        "Each lane service reserves 951957656+4194304*p_i bytes",
        "run_reserve_i=2046521648+8388608*p_i",
        "p_i in 2..=12",
        "p_i=2 reserves 2063298864 and leaves 84184784",
        "p_i=12 reserves 2147184944 and leaves 298704",
        "p_i=13 requires 2155573552 and is refused",
        "FINAL_SERVICE_CAPACITY_AND_QUARANTINE_CLOSURE_V3",
        "LaneServiceCapacityHeadBytesV3",
        "checked_add(active_reserved_bytes,quarantined_reserved_bytes)<=service_capacity_bytes",
        "exactly 8392 bytes",
        "exactly 8304 bytes",
        "LaneServiceCapacityUpdateReceiptBytesV3",
        "Each lane's preactivation cleanup is a tagged sum",
        "FINAL_REPLAY_ROUTE_FAILURE_AND_CATASTROPHE_CLOSURE_V3",
        "TerminalRouteFailureMatrixBytesV2",
        "exactly 742 bytes",
        "TerminalReplayRouteRegistryBytesV2",
        "PrimaryExactReplay=1",
        "Corrected ReplayRouteLedgerBytesV3",
        "exactly 385 bytes",
        "exactly 633 bytes",
        "FINAL_BRANCH_RETENTION_STORE_AND_CONTROL_LEDGER_CLOSURE_V3",
        "conjunctive at every multi-parent node except the names in `I13_FRESH_V2_AUTHORITY_TAGGED_SUM_NODES`",
        "ArmRetentionRequirementSetBytesV3",
        "Retention source authority is shard-aware",
        "StoreCapacityHeadBytesV3",
        "exactly 377 bytes",
        "ControlPlaneSettlementLedgerV3",
        "FINAL_CONDITIONAL_BREACH_SLOT_CLOSURE_V3",
        "phase2 bundle is a prerequisite of Finalization only when Present",
        "zero, one or two distinct lane-neutral breach artifact digests",
        "zero, two or four lane-specific persistence receipts",
        "I13_FRESH_V2_PRECEDENCE_ADDENDUM",
        "546-byte admission-core",
    ] {
        assert!(policy.contains(exact), "V2 authority lost '{exact}'");
    }

    assert_eq!(I13_PROTOCOL_FIXED_TERMINAL_EVENTS_V2.len(), 21);
    assert_eq!(
        I13_PROTOCOL_FIXED_TERMINAL_EVENTS_V2
            .iter()
            .map(|kind| kind.len())
            .sum::<usize>(),
        460
    );
    assert!(
        I13_PROTOCOL_FIXED_TERMINAL_EVENTS_V2
            .windows(2)
            .all(|pair| pair[0].as_bytes() < pair[1].as_bytes())
    );
    let draft = i13_draft();
    let obs_occurrence_count = draft
        .obligations
        .iter()
        .map(|obligation| obligation.obs_events.len())
        .sum::<usize>();
    let obs_kinds = draft
        .obligations
        .iter()
        .flat_map(|obligation| obligation.obs_events.iter().copied())
        .collect::<BTreeSet<_>>();
    let complete_scope2 = obs_kinds
        .iter()
        .copied()
        .chain(I13_PROTOCOL_FIXED_TERMINAL_EVENTS_V2.iter().copied())
        .collect::<BTreeSet<_>>();
    assert_eq!(obs_occurrence_count, 294);
    assert_eq!(obs_kinds.len(), 61);
    assert_eq!(
        obs_occurrence_count + I13_PROTOCOL_FIXED_TERMINAL_EVENTS_V2.len(),
        315
    );
    assert_eq!(complete_scope2.len(), 64);
    for registry_only in [
        "cancellation.proposed",
        "terminal.causal_barrier",
        "terminal.post_selection_capacity_breach",
    ] {
        assert!(!obs_kinds.contains(registry_only));
        assert!(complete_scope2.contains(registry_only));
    }

    let structured_registry_clause = policy
        .split_once("FINAL_STRUCTURED_PROTOCOL_EVENT_REGISTRY_CLOSURE_V3:")
        .expect("final structured protocol registry clause")
        .1
        .split_once("FINAL_DURABLE_ARTIFACT_PHASE_AND_CAPACITY_CLOSURE_V3:")
        .expect("structured registry clause boundary")
        .0;
    assert!(structured_registry_clause.contains("exactly 1054 bytes"));
    assert!(structured_registry_clause.contains("/65535/i"));
    assert!(structured_registry_clause.contains("294+21=315 occurrences"));

    let final_typed_append_clause = policy
        .split_once("FINAL_V3_TYPED_APPEND_CAPACITY_AND_PHASE_AUTHORITY_CLOSURE:")
        .expect("final V3 typed append/capacity/phase clause")
        .1
        .split_once("FINAL_V3_EXACT_ROUTE_RETENTION_STORE_AND_SETTLEMENT_SCHEMA_CLOSURE:")
        .expect("final V3 typed append clause boundary")
        .0;
    let final_exact_schema_clause = policy
        .split_once("FINAL_V3_EXACT_ROUTE_RETENTION_STORE_AND_SETTLEMENT_SCHEMA_CLOSURE:")
        .expect("final V3 route/retention/store/settlement clause")
        .1
        .split_once("I13_FRESH_V2_AUTHORITY_PRECEDENCE:")
        .expect("final V3 schema clause boundary")
        .0;

    let durable_capacity_clause = policy
        .split_once("FINAL_DURABLE_ARTIFACT_PHASE_AND_CAPACITY_CLOSURE_V3:")
        .expect("final durable capacity clause")
        .1
        .split_once("FINAL_SERVICE_CAPACITY_AND_QUARANTINE_CLOSURE_V3:")
        .expect("durable capacity clause boundary")
        .0;
    for exact in [
        "exactly 304 bytes",
        "at most 67108864 bytes",
        "at most 11849894 bytes",
        "at most 20373600 bytes",
        "261447712 bytes",
        "exactly 6987744 bytes",
        "exactly 546 bytes",
        "InfrastructureCapacityReservationSetBytesV3 is the 34 exact ASCII bytes I13_INFRASTRUCTURE_CAPACITY_SET_V3 followed by one byte 0x00, raw32(predecessor_manifest_digest), U64LE(attempt_id), U64LE(epoch), U16LE(member_count=9), and nine FRAME_BYTES(one byte service_kind||raw32(service_identity)||raw32(owner_identity)||raw32(implementation_identity)||raw32(durable_medium_identity)||raw32(failure_domain_identity)||U64LE(reserved_bytes)||raw32(capacity_reservation_receipt_digest)||raw32(retention_policy_digest_or_zero)||raw32(destination_roster_digest_or_zero)) rows in service-kind order, exactly 2542 bytes",
        "run_reserve_i=2046521648+8388608*p_i",
        "25166112 operation-payload-artifact bytes including six 48-byte rows",
        "13117392 terminal-payload-artifact bytes including five at-most-1024-record entry-overhead blocks",
        "2810954 structural-artifact bytes",
        "451054 framed delta-artifact bytes",
        "1062912 receipt-index-artifact bytes",
        "29844 append-metadata bytes",
        "4194352 operation-payload-artifact bytes including its row",
        "7464912 terminal-payload-artifact bytes including 80 bytes of entry overhead for each of 2817 records",
        "997689 structural bytes",
        "122989 framed delta bytes",
        "state Idle=0,Acquired=1,Prepared=2 or Quarantined=3",
        "transition_kind Acquire=1,SealPrepared=2,CommitRelease=3,AbortReclaim=4,Quarantine=5 or RecoveryRelease=6",
    ] {
        assert!(durable_capacity_clause.contains(exact));
    }

    for exact in [
        "BarrierBatch is legal only as a matrix-validated nested subprojection of Finalization and never as a top-level append, artifact row, phase lease or durability generation",
        "the six top-level tagged-sum arms exclude it",
        "AppendPhaseLeaseBytesV3 is the 25 exact ASCII bytes I13_APPEND_PHASE_LEASE_V3 followed by one byte 0x00",
        "U16LE(primary_append_slot_ordinal), U16LE(secondary_append_slot_ordinal_or_zero)",
        "raw32(LogicalAppendSlotAllocationReceiptDigestV3)",
        "exactly 376 bytes",
        "org.frankensim.i13.append-phase-lease.v3",
        "ComponentBytePlanBytesV3 is the 33-ASCII-plus-NUL header I13_APPEND_COMPONENT_BYTE_PLAN_V3",
        "six U64 component maxima in structural,event-payload,receipt-index,operation-payload,append-delta,append-metadata order",
        "exactly 153 bytes",
        "AppendPhaseHeadBytesV3 uses header I13_APPEND_PHASE_HEAD_V3 plus NUL, the exact 163-byte four-state projection already specified",
        "AppendPhaseTransitionReceiptBytesV3 uses header I13_APPEND_PHASE_LEASE_RECEIPT_V3 plus NUL",
        "exactly 317 bytes under 'org.frankensim.i13.append-phase-transition-receipt.v3'",
        "AppendPhaseReceiptArtifactBytesV3 uses header I13_APPEND_PHASE_RECEIPT_ARTIFACT_V3 and its 102-byte fixed projection plus at most 16384 FRAME_BYTES(the exact 317-byte receipts), at most 5324902 bytes",
        "AppendPhaseHeadArtifactBytesV3 uses the 33 exact ASCII bytes I13_APPEND_PHASE_HEAD_ARTIFACT_V3 followed by NUL",
        "U16LE(head_byte_count=163),one byte completeness=Closed=1",
        "exactly checked_add(102,checked_mul(head_count,171)) and at most 2801937 bytes",
        "head_count=receipt_count+1",
        "AppendSlotLedgerBytesV3 uses header I13_APPEND_SLOT_LEDGER_V3 plus NUL",
        "the exact 512-byte 4096-slot bitmap",
        "exactly 633 bytes",
        "LogicalAppendSlotAllocationIntentBytesV3 uses header I13_LOGICAL_APPEND_SLOT_ALLOCATION_INTENT_V3 plus NUL",
        "exactly 258 bytes",
        "LogicalAppendSlotAuthorityHeadBytesV3 uses header I13_LOGICAL_APPEND_SLOT_AUTHORITY_HEAD_V3 plus NUL",
        "exactly 164 bytes",
        "LogicalAppendSlotAllocationReceiptBytesV3 uses header I13_LOGICAL_APPEND_SLOT_ALLOCATION_RECEIPT_V3 plus NUL",
        "exactly 335 bytes",
        "One replicated verifier-owned CAS allocates each logical slot pair exactly once before either lane Acquire",
        "both lane leases byte-equal the same allocation receipt",
        "A one-lane Acquire failure leaves the global slot burned and forces deterministic paired abort/quarantine reconciliation",
        "AppendPhaseLeaseArtifactBytesV3 uses the 34-ASCII-plus-NUL header I13_APPEND_PHASE_LEASE_ARTIFACT_V3",
        "at most 1572964 bytes",
        "AppendComponentPlanArtifactBytesV3 uses the 33-ASCII-plus-NUL header I13_APPEND_BYTE_PLANS_ARTIFACT_V3",
        "at most 659555 bytes",
        "A successful append has exactly Acquire,SealPrepared,CommitRelease",
        "any burned slot has at most four non-replay transitions",
        "LaneControlClosureArtifactBytesV3 uses header I13_LANE_CONTROL_CLOSURE_ARTIFACT_V3 plus NUL",
        "at most 126402 bytes",
        "AppendPhaseAuthoritySetBytesV3 is header I13_APPEND_PHASE_AUTHORITY_SET_V3 plus NUL",
        "exactly 178 bytes",
        "5324902+2801937+1572964+659555+126402=10485760 bytes",
        "none is charged to the separate 268435456-byte durability bundle",
    ] {
        assert!(
            final_typed_append_clause.contains(exact),
            "final typed append/phase closure lost '{exact}'"
        );
    }

    for exact in [
        "TerminalCapacityProfileExtensionBytesV3 uses the 42 exact ASCII bytes I13_TERMINAL_CAPACITY_PROFILE_EXTENSION_V3 plus NUL",
        "exactly 435 bytes",
        "It is a pre-intent proposal containing no committed receipt, infrastructure-set root or final plan digest",
        "TerminalCapacityReservationIntentBytesV3 appends raw32(TerminalCapacityProfileExtensionRootV3) to the final 357-byte V2 form",
        "exactly 389 bytes",
        "TerminalCapacityReservationReceiptBytesV3 appends that same root to the final 502-byte V2 form",
        "exactly 534 bytes",
        "Only after the global receipt commits may the two plan identities and every physical service/store/control reservation be derived",
        "LaneReservationReceiptBytesV3 appends raw32(the lane's LaneReservationPlanIdentityDigestV3) to the 193-byte V2 form",
        "exactly 225 bytes",
        "OrdinaryEvidenceArtifactBytesV3 uses the 33 exact ASCII bytes I13_ORDINARY_EVIDENCE_ARTIFACT_V3 plus NUL",
        "at most 67108864 bytes",
        "LaneReservationPlanIdentityDigestV3 uses domain 'org.frankensim.i13.lane-reservation-plan-identity.v3' over the exact 121-byte projection",
        "U64LE(lane_control_bytes=10485760)",
        "LaneBreachSlotReservationBytesV3 uses the 35 exact ASCII bytes I13_LANE_BREACH_SLOT_RESERVATION_V3 plus NUL",
        "two phase-ordered FRAME_BYTES(one byte phase||raw32(slot_reservation_receipt)||raw32(genesis_BreachPersistenceHeadDigestV2)||U64LE(max_artifact_bytes)) rows, exactly 217 bytes",
        "LaneActiveCapacityMemberBytesV3 is raw32(plan identity)||raw32(LaneReservationReceiptDigestV3)||raw32(LaneBreachSlotReservationRootV3)||the seven U64 component values||U64LE(total_reserved_bytes), exactly 160 bytes",
        "total_reserved_bytes=951957656+4194304*p_i",
        "LaneActiveCapacityMapRootV3 is a 256-level authenticated sparse map keyed by plan identity and valued by the complete member",
        "LaneActiveCapacityUpdateProofBytesV3 binds old/new counts, old/new reserved totals, FRAME_BYTES(the 160-byte member) and 256 siblings, exactly 8392 bytes",
        "QuarantinedCapacityMemberBytesV3 is raw32(LaneActiveCapacityMemberDigestV3)||U64LE(total_reserved_bytes)||raw32(PreActivationServiceCatastropheDigestV2), exactly 72 bytes",
        "its proof remains 8304 bytes",
        "LaneServiceCapacityUpdateReceiptBytesV3 is exactly 447 bytes with the seven raw32 fields",
        "plan identity,active-member-or-zero,operation-evidence-or-zero,expected head,committed head,closed-plan proof,active-or-quarantine proof",
        "LaneReleaseAuthorizationBytesV3 uses the 33 exact ASCII bytes I13_LANE_RELEASE_AUTHORIZATION_V3 plus NUL",
        "exactly 356 bytes",
        "ActiveRelease operation-evidence is the pre-CAS LaneReleaseAuthorizationDigestV3",
        "The 218-byte LaneServiceCapacityHeadBytesV3 is an admission-control-plane replicated authority",
        "LaneBreachSlotSettlementMemberBytesV3 is one byte phase PreSelection=1 or PostSelection=2",
        "exactly 138 bytes",
        "LaneBreachSlotSettlementBytesV3 uses the 34 exact ASCII bytes I13_LANE_BREACH_SLOT_SETTLEMENT_V3 plus NUL",
        "exactly 378 bytes",
        "LaneReservationReleaseReceiptBytesV3 uses header I13_LANE_RESERVATION_RELEASE_RECEIPT_V3 and appends raw32(LaneBreachSlotSettlementRootV3)||raw32(LaneServiceCapacityUpdateReceiptDigestV3) to the V2-form fields, exactly 467 bytes",
    ] {
        assert!(
            final_typed_append_clause.contains(exact),
            "service-capacity closure lost exact schema '{exact}'"
        );
    }

    let route_clause = policy
        .split_once("FINAL_REPLAY_ROUTE_FAILURE_AND_CATASTROPHE_CLOSURE_V3:")
        .expect("final replay-route clause")
        .1
        .split_once("FINAL_BRANCH_RETENTION_STORE_AND_CONTROL_LEDGER_CLOSURE_V3:")
        .expect("replay-route clause boundary")
        .0;
    for exact in [
        "exactly 742 bytes",
        "exactly 656 bytes",
        "exact 1814-byte four-spec universe",
        "PrimaryExactReplay=1 with stages [1,2,3,4,5,6]",
        "MirroredPublicationReplay=2 with [1,5,6]",
        "EmergencyFailover=3 with [1,7,8,9,4,5,6]",
        "DurablePrefixRecovery=4 with [1,10,11,4,5,6]",
        "outcome Exhausted=1,Inapplicable=2 or Succeeded=3",
        "exactly 385 bytes",
        "exactly 633 bytes",
    ] {
        assert!(route_clause.contains(exact));
    }
    for exact in [
        "TerminalRouteSpecBytesV2 is the 26 exact ASCII bytes I13_TERMINAL_ROUTE_SPEC_V2 followed by one byte 0x00, U64LE(attempt_id), U64LE(epoch), U16LE(route_kind), U16LE(route_ordinal=0), raw32(terminal_authority_context_v2), raw32(activation_receipt_digest_v2), raw32(starting_published_authority_digest), raw32(target_operation_payload_digest), raw32(precondition_predicate_digest), raw32(program_schema_digest), raw32(success_predicate_digest), raw32(capability_set_root), raw32(clock_digest), raw32(deadline_ticket_digest), and raw32(idempotency_key_digest), exactly 399 bytes",
        "AdmissibleReplayRouteUniverseBytesV2 is the 39 exact ASCII bytes I13_ADMISSIBLE_REPLAY_ROUTE_UNIVERSE_V2 followed by one byte 0x00, raw32(predecessor_manifest_digest), U64LE(attempt_id), U64LE(epoch), raw32(activation_receipt_digest_v2), raw32(terminal_authority_context_v2), raw32(TerminalReplayRouteRegistryRootV2), U16LE(route_count=4), and four FRAME_BYTES(the exact 399-byte TerminalRouteSpecBytesV2) rows in registry order, exactly 1814 bytes",
        "TerminalRouteRequestBytesV2 is the 29 exact ASCII bytes I13_TERMINAL_ROUTE_REQUEST_V2 followed by one byte 0x00",
        "route_kind 1..=4 and route_ordinal=0 are never confused with physical append slots 4092..=4095",
        "org.frankensim.i13.admissible-replay-route-universe.v2",
    ] {
        assert!(
            final_exact_schema_clause.contains(exact),
            "final route universe lost '{exact}'"
        );
    }
    for exact in [
        "TerminalFailureReceiptBytesV2 uses exact order header/NUL,attempt,epoch,universe,spec,U16 kind,U16 ordinal,scope,stage,U16 failure_kind,effect,affected artifact,route request,evidence schema,evidence,expected pre-state,observed post-state,service,capability,clock,idempotency,independent verifier,signature set,interval, exactly 521 bytes",
        "RouteInapplicabilityProofBytesV2 uses exact order header/NUL,attempt,epoch,universe,spec,U16 kind,U16 ordinal,precondition,observed state,evaluator,one byte PredicateFalse=1,idempotency,signature set,interval, exactly 296 bytes",
        "TerminalRouteOutcomeReceiptBytesV2 uses exact order header/NUL,attempt,epoch,universe,spec,U16 kind,U16 ordinal,outcome,effect,starting authority,operation payload,request,failure-or-zero,inapplicability-or-zero,resulting authority-or-zero,observed route state,capability,idempotency,verifier,signature set,interval, exactly 492 bytes",
        "FailureInventoryBytesV3 uses header I13_FAILURE_INVENTORY_V3",
        "failure_scope Global=0,Primary=1 or Emergency=2",
        "A receipt occurs in exactly one inventory",
    ] {
        assert!(
            final_exact_schema_clause.contains(exact),
            "final route wire schema lost '{exact}'"
        );
    }

    let allowed_route_failures: [(&str, u8, &[u8]); 11] = [
        ("AuthorityAndFence", 1, &[1, 2, 3]),
        ("PrimaryAppend", 2, &[1, 2, 3, 4, 5, 6, 7, 8, 9]),
        (
            "MirrorAndPrefixVerification",
            3,
            &[1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        ),
        (
            "DurabilityCommit",
            4,
            &[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
        ),
        (
            "DualJournalPublication",
            5,
            &[1, 2, 3, 4, 5, 6, 7, 11, 12, 14],
        ),
        (
            "IndependentPublicationVerification",
            6,
            &[1, 2, 8, 9, 11, 12, 15],
        ),
        ("FailoverTransition", 7, &[1, 2, 3, 4, 5, 7, 11, 12]),
        ("EmergencyPrefix", 8, &[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]),
        (
            "EmergencyFinalization",
            9,
            &[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
        ),
        (
            "DownstreamReceiptRecovery",
            10,
            &[1, 2, 3, 4, 5, 8, 9, 11, 12],
        ),
        (
            "ArtifactReconstruction",
            11,
            &[1, 2, 3, 8, 9, 11, 12, 13, 15],
        ),
    ];
    assert_eq!(
        allowed_route_failures
            .iter()
            .map(|(_, _, codes)| codes.len())
            .sum::<usize>(),
        99
    );
    for (stage_name, stage_ordinal, codes) in allowed_route_failures {
        assert!(
            codes.windows(2).all(|pair| pair[0] < pair[1]),
            "stage {stage_name} failure codes must be unique and sorted"
        );
        assert!(codes.iter().all(|code| (1..=15).contains(code)));
        let encoded_codes = codes
            .iter()
            .map(u8::to_string)
            .collect::<Vec<_>>()
            .join(",");
        let exact_set = format!("S{stage_ordinal}={{{encoded_codes}}}");
        assert!(
            final_exact_schema_clause.contains(&exact_set),
            "route-failure matrix lost exact set {stage_name}/{exact_set}"
        );
    }

    for exact in [
        "RetainedArtifactMemberBytesV3 keeps the exact 43-byte lane/role/digest/byte-count projection but uses RetentionArtifactRoleV3",
        "RetainedArtifactSetBytesV3 uses header I13_RETAINED_ARTIFACT_SET_V3, set_kind PreEvidenceSurvivors=1,TransferSet=2,RequirementSource=3 or SourceShard=4",
        "exactly checked_add(62,checked_mul(artifact_count,51))",
        "ArmRetentionRequirementSetBytesV3 uses header I13_ARM_RETENTION_REQUIREMENT_SET_V3",
        "every required_source_artifact_set_root fetches a RequirementSource set",
        "role5 is present iff either logical breach union is Present",
        "role6 is present iff any route outcome is Exhausted or Succeeded or any Failover/Recovery lease was acquired",
        "Role10's source set is the canonical singleton containing CanonicalFinalizationEvidenceBytesV3",
        "RetentionConformanceReceiptBytesV3 uses header I13_RETENTION_CONFORMANCE_RECEIPT_V3",
        "U16LE(requirement_count),U64LE(artifact_count)",
        "exactly 305 bytes",
        "RetentionSourceShardReceiptBytesV3 uses the 37 exact ASCII bytes I13_RETENTION_SOURCE_SHARD_RECEIPT_V3 plus NUL",
        "raw32(source_store_identity),raw32(TransferSetRootV3),raw32(SourceShardSetRootV3),raw32(expected_source_inventory_head),raw32(committed_source_inventory_head),raw32(source_capability),raw32(idempotency),raw32(signature_set)",
        "exactly 345 bytes",
        "RetentionSourceShardReceiptSetBytesV3 uses the 41 exact ASCII bytes I13_RETENTION_SOURCE_SHARD_RECEIPT_SET_V3 plus NUL",
        "exactly checked_add(109,checked_mul(receipt_count,40))",
        "Its SourceShard sets form an exact duplicate-free partition",
        "StoreCapacityReservationMemberBytesV3 has two exact tagged arms",
        "Held=1 is exactly 225 bytes with used_bytes=0",
        "Retained=2 appends raw32(RetainedArtifactSetRootV3)||raw32(RetentionStoreReceiptDigestV3), exactly 289 bytes",
        "StoreCapacityHeadBytesV3 uses header I13_STORE_CAPACITY_HEAD_V3 plus NUL",
        "exactly 211 bytes",
        "StoreCapacityUpdateProofBundleBytesV3 uses the 41 exact ASCII bytes I13_STORE_CAPACITY_UPDATE_PROOF_BUNDLE_V3 plus NUL",
        "at most 17150 bytes",
        "ClosedIntentProofBytesV3 is exactly 8248 bytes",
        "ActiveCapacityProofBytesV3 binds old/new counts and reserved totals plus framed old/new Held-or-Retained members and one 256-sibling path, at most 8818 bytes",
        "QuarantineProofBytesV3 is exactly 8304 bytes",
        "StoreCapacityUpdateReceiptBytesV3 uses header I13_STORE_CAPACITY_UPDATE_RECEIPT_V3 plus NUL",
        "eleven raw32 fields store,member-or-zero,expected head,committed head,proof-bundle,evidence-or-zero,policy,roster,capability,idempotency,signature",
        "exactly 423 bytes",
        "RetentionStoreReceiptBytesV3 appends only raw32(ReserveCommittedUpdateReceiptDigestV3) to the prior 345-byte write receipt, exactly 377 bytes",
        "The acyclic order is ReserveCommitted -> exact storage CAS/write receipt -> CommitRetained",
        "ControlPlaneLedgerEntryBytesV3 uses header I13_CONTROL_PLANE_LEDGER_ENTRY_V3 plus NUL",
        "exactly checked_add(131,payload_byte_count)",
        "ControlPlaneLedgerHeadBytesV3 uses header I13_CONTROL_PLANE_LEDGER_HEAD_V3 plus NUL",
        "exactly 121 bytes",
        "TerminalCapacitySettlementReceiptBytesV3 uses header I13_TERMINAL_CAPACITY_SETTLEMENT_RECEIPT_V3",
        "exactly 582 bytes",
        "ControlPlaneSettlementMarkerBytesV3 uses header I13_CONTROL_PLANE_SETTLEMENT_MARKER_V3 plus NUL",
        "exactly 312 bytes",
        "The only order is append complete pre-settlement evidence E, commit the global settlement receipt binding E and its successor head, then append the marker binding that receipt",
        "A replicated cross-attempt ControlPlaneCapacityHeadBytesV3",
        "exactly 131 bytes",
        "pre-reserves the full 67108864-byte cleanup slice before activation",
        "BreachEvidenceBundleBytesV3 uses header I13_BREACH_EVIDENCE_BUNDLE_V3 plus NUL",
        "exactly checked_add(47,checked_mul(phase_count,145))",
        "exactly 0,1 or 2 distinct lane-neutral logical breach-artifact digests, 0,2 or 4 lane-specific persistence-receipt preimages and 0,1 or 2 mirror certificates",
        "Any 0/2/4-artifact interpretation",
    ] {
        assert!(
            final_exact_schema_clause.contains(exact),
            "retention/control closure lost exact schema '{exact}'"
        );
    }
    assert!(policy.contains(
        "Roles are TerminalStructural=1,TerminalPayload=2,ReceiptIndex=3,AppendMetadataAndDelta=4,BreachEvidence=5,FailureRecoveryEvidence=6,FinalEnvelope=7,PublicationEvidence=8,QuiescenceEvidence=9,SettlementEvidence=10,OperationPayload=11,OrdinaryObservability=12 and CleanupAuthority=13"
    ));

    let corrected_event = policy
        .split_once("TERMINAL_EVENT_DIGEST_CAUSAL_CONTEXT_AND_FRONTIER_CORRECTION_V2:")
        .expect("fresh V2 event correction")
        .1;
    let digest_projection = corrected_event
        .split_once("terminal_event_digest_v2=")
        .expect("corrected event digest")
        .0;
    assert!(digest_projection.contains("raw32(terminal_authority_context_v2)"));
    assert!(!digest_projection.contains("raw32(terminal_authority_admission_digest_v2)"));
    let fixed_context = corrected_event
        .split_once("TerminalFrontierContextV2 is exactly ")
        .expect("fixed terminal context")
        .1
        .split_once(", exactly 112 bytes")
        .expect("fixed context boundary")
        .0;
    assert!(!fixed_context.contains("receipt_index_root"));
    assert!(!fixed_context.contains("terminal_selection_receipt_digest"));

    let mut nodes = BTreeSet::new();
    let mut indegree = BTreeMap::new();
    let mut outgoing: BTreeMap<&str, Vec<&str>> = BTreeMap::new();
    for &(prerequisite, consumer) in I13_FRESH_V2_AUTHORITY_DAG_EDGES {
        assert_ne!(prerequisite, consumer, "self edge at {prerequisite}");
        assert!(
            nodes.insert((prerequisite, consumer)),
            "duplicate authority edge {prerequisite}->{consumer}"
        );
        indegree.entry(prerequisite).or_insert(0usize);
        *indegree.entry(consumer).or_insert(0usize) += 1;
        outgoing.entry(prerequisite).or_default().push(consumer);
    }
    let mut ready = indegree
        .iter()
        .filter_map(|(&node, &degree)| (degree == 0).then_some(node))
        .collect::<BTreeSet<_>>();
    let mut visited = 0usize;
    while let Some(node) = ready.iter().next().copied() {
        ready.remove(node);
        visited += 1;
        for &consumer in outgoing.get(node).into_iter().flatten() {
            let degree = indegree.get_mut(consumer).expect("consumer node");
            *degree = degree.checked_sub(1).expect("nonzero indegree");
            if *degree == 0 {
                ready.insert(consumer);
            }
        }
    }
    assert_eq!(
        visited,
        indegree.len(),
        "fresh V2 authority graph contains a backward/cyclic edge"
    );

    assert_eq!(
        I13_FRESH_V2_AUTHORITY_DAG_EDGES.len(),
        1_373,
        "fresh authority graph must be the frozen V4 prerequisite graph"
    );
    assert_eq!(indegree.len(), 630, "fresh V4 authority node count drifted");
    for required_edge in [
        (
            "protocol_fixed_terminal_event_registry",
            "semantic_source_occurrence_sets",
        ),
        ("semantic_source_occurrence_sets", "semantic_source_catalog"),
        ("semantic_occurrence_conformance_set", "semantic_member_set"),
        ("infrastructure_capacity_reservation_set", "admission_core"),
        ("lane_reservation_authority_set", "admission_core"),
        (
            "published_terminal_authority_n",
            "primary_append_request_n_plus_1",
        ),
        (
            "operation_payload_arm_sum_n_plus_1",
            "event_prepare_n_plus_1",
        ),
        (
            "operation_payload_arm_sum_n_plus_1",
            "durability_commit_n_plus_1",
        ),
        (
            "operation_payload_artifact_n",
            "operation_payload_artifact_n_plus_1",
        ),
        (
            "terminal_append_delta_artifact_n",
            "terminal_append_delta_artifact_n_plus_1",
        ),
        (
            "append_metadata_artifact_n",
            "append_metadata_artifact_n_plus_1",
        ),
        (
            "primary_lane_reservation",
            "primary_breach_artifact_receipt",
        ),
        (
            "postselection_breach_durability_bundle",
            "postselection_breach_union_sum",
        ),
        (
            "for_each_capacity_transaction_durability_set_v4",
            "for_each_capacity_effect_visibility_gate_v4",
        ),
        (
            "global_admission_capacity_transaction_v4",
            "global_admission_transaction_receipt",
        ),
        (
            "bootstrap_capacity_copy_a_pre_effect_absent",
            "bootstrap_capacity_copy_a_decision_sum",
        ),
        (
            "bootstrap_capacity_copy_b_pre_effect_absent",
            "bootstrap_capacity_copy_b_decision_sum",
        ),
        (
            "control_evidence_staging_seal_durability_set_v4",
            "control_evidence_cutoff_v4",
        ),
        (
            "control_evidence_cutoff_v4",
            "control_plane_pre_settlement_evidence",
        ),
        ("replay_route_ledger", "catastrophe_evidence"),
        ("route_outcome_receipt_set", "replay_route_ledger"),
        (
            "normal_global_closure_transaction",
            "normal_global_closure_authority_set",
        ),
        (
            "normal_global_closure_authority_set",
            "normal_closure_evidence_archive_v4",
        ),
        (
            "normal_closure_evidence_artifact_v4",
            "global_closure_evidence_artifact_sum",
        ),
        (
            "global_closure_evidence_artifact_sum",
            "global_journal_expiry_request",
        ),
        (
            "global_journal_expiry_request",
            "global_journal_expiry_proof_bundle",
        ),
        (
            "global_journal_expiry_proof_bundle",
            "global_journal_expiry_prepared_payload",
        ),
        (
            "bootstrap_capacity_allocation_request_a",
            "bootstrap_capacity_prepare_decision_a",
        ),
        (
            "bootstrap_capacity_allocation_request_b",
            "bootstrap_capacity_prepare_decision_b",
        ),
        (
            "bootstrap_capacity_prepare_decision_a",
            "bootstrap_capacity_prepare_entry_a",
        ),
        (
            "bootstrap_capacity_prepare_decision_b",
            "bootstrap_capacity_prepare_entry_b",
        ),
        (
            "bootstrap_capacity_prepare_entry_a",
            "bootstrap_capacity_decision_head_a_prepared",
        ),
        (
            "bootstrap_capacity_prepare_entry_b",
            "bootstrap_capacity_decision_head_b_prepared",
        ),
        (
            "bootstrap_capacity_decision_head_a_prepared",
            "bootstrap_capacity_copy_coordinator_head_a_prepared",
        ),
        (
            "bootstrap_capacity_decision_head_b_prepared",
            "bootstrap_capacity_copy_coordinator_head_b_prepared",
        ),
        (
            "bootstrap_capacity_copy_coordinator_head_a_prepared",
            "bootstrap_capacity_prepare_append_receipt_a",
        ),
        (
            "bootstrap_capacity_copy_coordinator_head_b_prepared",
            "bootstrap_capacity_prepare_append_receipt_b",
        ),
        (
            "bootstrap_capacity_prepare_append_receipt_a",
            "bootstrap_capacity_local_decision_a",
        ),
        (
            "bootstrap_capacity_prepare_append_receipt_b",
            "bootstrap_capacity_local_decision_b",
        ),
        (
            "bootstrap_capacity_local_decision_a",
            "bootstrap_capacity_local_terminal_entry_a",
        ),
        (
            "bootstrap_capacity_local_decision_b",
            "bootstrap_capacity_local_terminal_entry_b",
        ),
        (
            "bootstrap_capacity_local_terminal_entry_a",
            "bootstrap_capacity_decision_head_a_terminal",
        ),
        (
            "bootstrap_capacity_local_terminal_entry_b",
            "bootstrap_capacity_decision_head_b_terminal",
        ),
        (
            "bootstrap_capacity_decision_head_a_terminal",
            "bootstrap_capacity_copy_coordinator_head_a_terminal",
        ),
        (
            "bootstrap_capacity_decision_head_b_terminal",
            "bootstrap_capacity_copy_coordinator_head_b_terminal",
        ),
        (
            "bootstrap_capacity_copy_coordinator_head_a_terminal",
            "bootstrap_capacity_local_terminal_append_receipt_a",
        ),
        (
            "bootstrap_capacity_copy_coordinator_head_b_terminal",
            "bootstrap_capacity_local_terminal_append_receipt_b",
        ),
        (
            "bootstrap_capacity_local_terminal_append_receipt_a",
            "bootstrap_capacity_allocation_receipt_a",
        ),
        (
            "bootstrap_capacity_local_terminal_append_receipt_b",
            "bootstrap_capacity_allocation_receipt_b",
        ),
        (
            "bootstrap_capacity_allocation_outcome_v4",
            "bootstrap_capacity_aggregate_outcome_entry_a",
        ),
        (
            "bootstrap_capacity_allocation_outcome_v4",
            "bootstrap_capacity_aggregate_outcome_entry_b",
        ),
        (
            "bootstrap_capacity_aggregate_outcome_entry_a",
            "bootstrap_capacity_decision_head_a_outcome",
        ),
        (
            "bootstrap_capacity_aggregate_outcome_entry_b",
            "bootstrap_capacity_decision_head_b_outcome",
        ),
        (
            "bootstrap_capacity_decision_head_a_outcome",
            "bootstrap_capacity_copy_coordinator_head_a_outcome",
        ),
        (
            "bootstrap_capacity_decision_head_b_outcome",
            "bootstrap_capacity_copy_coordinator_head_b_outcome",
        ),
        (
            "bootstrap_capacity_copy_coordinator_head_a_outcome",
            "bootstrap_capacity_aggregate_outcome_append_receipt_a",
        ),
        (
            "bootstrap_capacity_copy_coordinator_head_b_outcome",
            "bootstrap_capacity_aggregate_outcome_append_receipt_b",
        ),
        (
            "bootstrap_capacity_aggregate_outcome_append_receipt_a",
            "bootstrap_capacity_allocation_evidence_bundle",
        ),
        (
            "bootstrap_capacity_aggregate_outcome_append_receipt_b",
            "bootstrap_capacity_allocation_evidence_bundle",
        ),
        (
            "capacity_allocation_initial_pending_outcome_v4",
            "capacity_orphan_reconciliation_prepare_decision_v4",
        ),
        (
            "capacity_orphan_reconciliation_prepare_decision_v4",
            "capacity_orphan_ordinal_member_v4",
        ),
        (
            "capacity_allocation_initial_pending_outcome_v4",
            "capacity_orphan_ordinal_member_v4",
        ),
        (
            "capacity_orphan_ordinal_member_v4",
            "capacity_orphan_insertion_proof_v4",
        ),
        (
            "capacity_orphan_insertion_proof_v4",
            "capacity_orphan_candidate_root_v4",
        ),
        (
            "capacity_orphan_candidate_root_v4",
            "capacity_orphan_reconciliation_local_terminal_v4",
        ),
        (
            "capacity_orphan_reconciliation_local_terminal_v4",
            "capacity_orphan_terminal_append_receipt_v4",
        ),
        (
            "capacity_orphan_terminal_append_receipt_v4",
            "capacity_allocation_final_reconciled_no_effect_outcome_v4",
        ),
    ] {
        assert!(
            nodes.contains(&required_edge),
            "fresh V2 graph lost required edge {} -> {}",
            required_edge.0,
            required_edge.1
        );
    }

    for forbidden_edge in [
        (
            "control_plane_settlement_entry_sum",
            "control_plane_settlement_ledger_successor",
        ),
        (
            "preactivation_abort_evidence",
            "preactivation_global_active_reservation_set_n_plus_1",
        ),
        (
            "normal_global_closure_authority_set",
            "normal_global_closure_entry",
        ),
        (
            "bootstrap_global_closure_authority_set",
            "bootstrap_global_closure_entry",
        ),
        (
            "global_expiry_durability_authority_set",
            "global_journal_expiry_request",
        ),
        (
            "for_each_capacity_prepared_payload_v4",
            "for_each_capacity_effect_visibility_gate_v4",
        ),
        (
            "capacity_orphan_reconciliation_local_terminal_v4",
            "capacity_orphan_ordinal_member_v4",
        ),
        (
            "capacity_allocation_final_reconciled_no_effect_outcome_v4",
            "capacity_orphan_ordinal_member_v4",
        ),
        (
            "capacity_allocation_final_reconciled_no_effect_outcome_v4",
            "capacity_orphan_insertion_proof_v4",
        ),
    ] {
        assert!(
            !nodes.contains(&forbidden_edge),
            "fresh V4 graph restored forbidden/backward edge {} -> {}",
            forbidden_edge.0,
            forbidden_edge.1
        );
    }

    let tagged_sum_nodes = I13_FRESH_V2_AUTHORITY_TAGGED_SUM_NODES
        .iter()
        .copied()
        .collect::<BTreeSet<_>>();
    let expected_tagged_sum_arms = BTreeMap::from([
        ("primary_lane_attempt_outcome_sum", 2usize),
        ("emergency_lane_attempt_outcome_sum", 2),
        ("primary_preactivation_cleanup_sum", 3),
        ("emergency_preactivation_cleanup_sum", 3),
        ("operation_payload_arm_sum_n_plus_1", 6),
        ("postselection_breach_union_sum", 2),
        ("catastrophe_latest_published_authority_sum", 2),
        ("control_plane_settlement_entry_sum", 2),
        ("control_plane_postsettlement_receipt_sum", 2),
        ("final_breach_requirement_source_sum", 2),
        ("final_failure_recovery_requirement_source_sum", 2),
        ("catastrophe_breach_survivor_source_sum", 2),
        ("cleanup_store_failure_sum", 3),
        ("logical_slot_survivor_sum", 2),
        ("append_phase_closure_sum", 2),
        ("control_settlement_transport_sum", 2),
        ("retained_store_artifact_set_sum", 2),
        ("normal_h2_quiescence_source_sum", 2),
        ("global_closure_evidence_artifact_sum", 2),
        ("global_expiry_ticket_preimage_sum", 2),
        ("global_expiry_status_source_sum", 2),
        ("bootstrap_capacity_copy_a_decision_sum", 2),
        ("bootstrap_capacity_copy_b_decision_sum", 2),
        ("bootstrap_capacity_legal_copy_combination_sum", 3),
        ("bootstrap_capacity_initial_outcome_durability_sum", 3),
        ("for_each_capacity_legal_copy_combination_sum", 3),
        ("for_each_capacity_initial_outcome_durability_sum", 3),
        ("for_each_capacity_trace_completion_sum", 3),
        ("bootstrap_capacity_phase3_missing_a_evidence_sum", 4),
        ("bootstrap_capacity_phase3_missing_b_evidence_sum", 4),
        ("for_each_capacity_phase3_missing_a_evidence_sum", 4),
        ("for_each_capacity_phase3_missing_b_evidence_sum", 4),
        ("bootstrap_capacity_reserved_trace_completion_sum", 2),
        ("for_each_capacity_no_effect_completion_sum", 3),
        ("for_each_capacity_terminal_authority_sum", 2),
        ("for_each_capacity_repair_entry_sum", 2),
        ("for_each_capacity_phase3_append_repair_completion_sum", 2),
        ("capacity_allocation_affected_copy_trace_sum", 2),
        ("for_each_capacity_phase1_prepare_entry_a_sum", 2),
        ("for_each_capacity_phase1_prepare_entry_b_sum", 2),
        ("for_each_capacity_recovered_phase3_completion_a_sum", 2),
        ("for_each_capacity_recovered_phase3_completion_b_sum", 2),
        ("for_each_capacity_phase6_completion_sum", 4),
    ]);
    assert_eq!(
        tagged_sum_nodes,
        expected_tagged_sum_arms.keys().copied().collect()
    );
    for (&sum_node, &expected_arm_count) in &expected_tagged_sum_arms {
        let actual_arm_count = I13_FRESH_V2_AUTHORITY_DAG_EDGES
            .iter()
            .filter(|(_, consumer)| *consumer == sum_node)
            .count();
        assert_eq!(
            actual_arm_count, expected_arm_count,
            "tagged sum {sum_node} lost or gained an arm"
        );
    }
    for &(_, consumer) in I13_FRESH_V2_AUTHORITY_DAG_EDGES {
        if consumer.ends_with("_sum") || consumer.contains("_sum_") {
            assert!(
                tagged_sum_nodes.contains(consumer),
                "unregistered OR node {consumer} would be treated as a conjunction"
            );
        }
    }

    for (start, target) in [
        (
            "semantic_source_occurrence_sets",
            "published_terminal_authority_n_plus_1",
        ),
        ("execution_request", "published_terminal_authority_n_plus_1"),
        (
            "published_terminal_authority_n",
            "published_terminal_authority_n_plus_1",
        ),
        ("capacity_authority_roster_matrix_v4", "admission_core"),
        (
            "for_each_capacity_transaction_durability_set_v4",
            "global_admission_transaction_receipt",
        ),
        (
            "control_evidence_staging_member_set_v4",
            "control_plane_pre_settlement_evidence",
        ),
        (
            "global_closure_evidence_artifact_sum",
            "global_expiry_durability_authority_set",
        ),
    ] {
        let mut frontier = BTreeSet::from([start]);
        let mut reached = BTreeSet::new();
        while let Some(node) = frontier.iter().next().copied() {
            frontier.remove(node);
            if !reached.insert(node) {
                continue;
            }
            for &consumer in outgoing.get(node).into_iter().flatten() {
                frontier.insert(consumer);
            }
        }
        assert!(
            reached.contains(target),
            "fresh V2 graph lost required path {start} -> {target}"
        );
    }
}

#[test]
fn i13_v4_finite_capacity_wal_and_closed_evidence_authority_is_exact() {
    let policy = policy_spec();
    let final_v4 = policy
        .split_once("FINAL_V4_FINITE_CAPACITY_WAL_AND_CLOSED_EVIDENCE_AUTHORITY:")
        .expect("terminal fresh V4 authority clause")
        .1
        .split_once("I13_FRESH_V2_AUTHORITY_PRECEDENCE:")
        .expect("terminal fresh V4 authority boundary")
        .0;
    let (base_v4, corrected_allocation_and_after) = final_v4
        .split_once("FINAL_V4_DURABLE_ALLOCATION_JOURNAL_CORRECTION:")
        .expect("terminal V4 allocation-journal correction");
    let (allocation_correction, corrected_coordinator_and_after) = corrected_allocation_and_after
        .split_once("FINAL_V4_COORDINATOR_AND_STAGING_CAPACITY_CORRECTION:")
        .expect("terminal V4 coordinator correction");
    let (coordinator_correction, wire_and_after) = corrected_coordinator_and_after
        .split_once("FINAL_V4_WIRE_AND_ACCOUNTING_DISAMBIGUATION:")
        .expect("terminal V4 wire correction");
    let (wire_correction, gauntlet_obligations) = wire_and_after
        .split_once("The V4 decision journals,response logs,WALs")
        .expect("terminal V4 Gauntlet boundary");

    for header in [
        "I13_CAPACITY_AUTHORITY_BINDING_MATRIX_V4",
        "I13_CAPACITY_AUTHORITY_ROSTER_MATRIX_V4",
        "I13_CAPACITY_INDEX_ALLOCATION_EVIDENCE_BUNDLE_V4",
        "I13_CAPACITY_INDEX_ALLOCATION_RECEIPT_V4",
        "I13_CAPACITY_INDEX_ALLOCATION_REQUEST_V4",
        "I13_CAPACITY_INDEX_ALLOCATION_AUTHORITY_SET_V4",
        "I13_CAPACITY_INDEX_ALLOCATION_OUTCOME_V4",
        "I13_CAPACITY_INDEX_ALLOCATION_RECEIPT_PAIR_V4",
        "I13_CAPACITY_INDEX_COPY_COORDINATOR_V4",
        "I13_CAPACITY_INDEX_DECISION_APPEND_V4",
        "I13_CAPACITY_INDEX_DECISION_ENTRY_V4",
        "I13_CAPACITY_INDEX_DECISION_HEAD_V4",
        "I13_CAPACITY_INDEX_HEAD_V4",
        "I13_CAPACITY_INDEX_LOCAL_DECISION_V4",
        "I13_CAPACITY_INDEX_ORPHAN_HEAD_V4",
        "I13_CAPACITY_INDEX_ORPHAN_PROOF_V4",
        "I13_CAPACITY_INDEX_RESERVATION_SUBJECT_V4",
        "I13_CAPACITY_MUTATION_COMMIT_INTENT_V4",
        "I13_CAPACITY_MUTATION_PREPARED_PAYLOAD_V4",
        "I13_CAPACITY_MUTATION_EVIDENCE_ARCHIVE_V4",
        "I13_CAPACITY_MUTATION_EVIDENCE_MANIFEST_V4",
        "I13_CAPACITY_MUTATION_PRECOMMIT_MANIFEST_V4",
        "I13_CAPACITY_OPERATION_PROTOCOL_MATRIX_V4",
        "I13_CAPACITY_RESPONSE_LOG_APPEND_RECEIPT_V4",
        "I13_CAPACITY_RESPONSE_LOG_ENTRY_V4",
        "I13_CAPACITY_RESPONSE_LOG_HEAD_V4",
        "I13_CAPACITY_TRANSACTION_CHAIN_STEP_V4",
        "I13_CAPACITY_TRANSACTION_COORDINATOR_HEAD_V4",
        "I13_CAPACITY_TRANSACTION_DURABILITY_SET_V4",
        "I13_CAPACITY_TRANSACTION_RECEIPT_V4",
        "I13_CAPACITY_TRANSACTION_UPDATE_PROOF_V4",
        "I13_CAPACITY_TRANSACTION_WAL_APPEND_RECEIPT_V4",
        "I13_CAPACITY_TRANSACTION_WAL_ENTRY_V4",
        "I13_CAPACITY_TRANSACTION_WAL_HEAD_V4",
        "I13_CAPACITY_TRANSACTION_WAL_RECORD_V4",
        "I13_CONTROL_EVIDENCE_CUTOFF_V4",
        "I13_CONTROL_EVIDENCE_PRODUCER_INVENTORY_V4",
        "I13_CONTROL_EVIDENCE_STAGING_HEAD_V4",
        "I13_CONTROL_EVIDENCE_STAGING_INSERTION_V4",
        "I13_CONTROL_EVIDENCE_STAGING_PROOF_BUNDLE_V4",
        "I13_CONTROL_EVIDENCE_STAGING_SEAL_PREPARED_V4",
        "I13_CONTROL_EVIDENCE_STAGING_ARTIFACT_V4",
        "I13_CONTROL_EVIDENCE_STAGING_LOG_PAGE_V4",
        "I13_GLOBAL_CLOSED_EVIDENCE_MAP_HEAD_V4",
        "I13_GLOBAL_CLOSED_MEMBER_V4",
        "I13_GLOBAL_CLOSURE_EVIDENCE_ARTIFACT_V4",
        "I13_GLOBAL_JOURNAL_EXPIRY_PROOF_BUNDLE_V4",
        "I13_GLOBAL_JOURNAL_EXPIRY_REQUEST_V4",
    ] {
        assert!(
            final_v4.contains(header),
            "terminal V4 authority lost wire header {header}"
        );
    }

    for exact in [
        "terminal fresh authority and supersedes every conflicting V3 capacity-index quota,shared-slot journal,transaction,bootstrap,expiry,staging or downstream-binding definition above",
        "exactly fifteen logical authorities in fixed order",
        "one Global lifecycle authority; seven ReplayMetadata authorities",
        "four physical Generic authorities",
        "three physical Store authorities",
        "authority count=15",
        "exactly 2609 bytes under 'org.frankensim.i13.capacity-authority-roster-matrix.v4'",
        "response-log quota=262144",
        "transaction-WAL quota=131072",
        "exactly266 bytes under 'org.frankensim.i13.capacity-index-reservation-subject.v4'",
        "Rows are strictly (tag,authority)-ordered and duplicate-free",
        "Its closed disposition state machine is Unused=0,Prepared=1,IntentDurable=2,Committed=3,AbortedNoEffect=4,ConflictLost=5,ReconciledCommitted=6,ReconciledAborted=7 or PendingRepair=8",
        "Thus one coordinator-head CAS is the sole atomic visibility point for all one- or two-component operations",
        "An immutable Prepare is never discarded: it reaches exactly one Committed,AbortedNoEffect,ConflictLost,ReconciledCommitted or ReconciledAborted terminal transaction",
        "The postcommit manifest excludes its own durability set",
        "neither manifest hashes a future object",
        "If K is the deterministic attempted roster prefix,1<=K<=15",
        "A<=2K<=30",
        "M is at most 1+2*(K-1)<=29",
        "item_count=A+M+4<=63,not a guessed12/24/384 cap",
        "Bootstrap payload is at most19321 and its global closure entry at most19454 bytes",
        "Bootstrap never enters the normal ControlPlane E/H2/H3 path",
        "AdmissionCore and lane/service issuance bind GlobalAdmission durability,not the566-byte PreparedOnly primary",
        "H3 cannot precede H2 durability",
        "A separate owner-index map keyed only by extent digest permits exactly one PhysicalOwner",
        "No insertion may follow that durability set",
        "The ControlPlane E layout contains the exact staging artifact,exact seal durability-set projection and cutoff as three special PhysicalOwner extents",
        "Zero dispositions and references never enter the positive-length gap-free physical layout",
        "the H1/marker/H2/postreceipt/H3/footer tail is2489,and slack is63047",
        "crash at every local-index/response/WAL/coordinator boundary",
        "Prepared->every terminal disposition",
        "no backward digest among precommit manifest,Intent,chain step,WAL,candidate head,postcommit manifest and durability set",
        "every stated capacity equality/cap+1",
    ] {
        assert!(
            final_v4.contains(exact),
            "terminal V4 authority lost exact requirement '{exact}'"
        );
    }

    for exact in [
        "24576 nonborrowable 32768-byte decision slots, exactly805306368 bytes",
        "total fixed allocation baseline of2415919104 bytes",
        "U16LE(reserved_subject_count<=4096)",
        "used_bytes=checked_mul(record_count,32768)",
        "record_count<=checked_mul(6,reserved_subject_count)",
        "first durable Prepare for a new subject atomically increments reserved_subject_count,sets only bit1,appends exactly one decision entry and assigns physical_ordinal=old_reserved_subject_count",
        "new_decision_generation=old+1,new_record_count=old+1,new_used=old_used+32768",
        "Exact-key replay returns the byte-identical receipt and advances neither head/count",
        "exact signed/interval-bearing CapacityIndexAllocationOutcomeBytesV4 is appended as phase3 to every locally present preprovisioned decision journal before any newly allocated response log is trusted",
        "phase1 plus phase3 -> orphan member/proof/candidate root -> phase4 append receipt -> phase5 append receipt -> phase6 final-outcome append receipt",
    ] {
        assert!(
            allocation_correction.contains(exact),
            "terminal allocation correction lost proof obligation '{exact}'"
        );
    }

    for exact in [
        "new_generation=checked_add(old_generation,1),new_committed_transaction_count=checked_add(old_committed_transaction_count,1)",
        "Exact-key replay changes none; abort before CAS has no head/count authority",
        "CapacityTransactionWalRecordBytesV4 is therefore checked_add(676,proof_bytes)",
        "CapacityTransactionDurabilitySetBytesV4 is checked_add(3299,checked_add(primary_bytes,checked_add(checked_mul(component_count,138),proof_bytes)))",
        "Corrected normal maxima are20582 ReplayMetadata,20633 Generic,20648 Store,20624 ReplayArchive,20759 GlobalAdmission,20807 H2,20648 Expiry,20653 StagingInsertion and20547 StagingSeal; two-component GlobalClosure is29211",
        "at most114657",
        "at most114822",
        "at most123062",
        "at most131669",
        "checked_add(Smax,checked_add(checked_mul(21319,Nmax),checked_add(Lmax,checked_add(20547,66505))))<=67042864",
        "Equality passes; any byte,count,page,decision-slot or index-cap +1",
        "checked_add(fixed_cleanup_subject_count,page_count)<=4096",
    ] {
        assert!(
            coordinator_correction.contains(exact),
            "terminal coordinator/staging correction lost proof obligation '{exact}'"
        );
    }

    for exact in [
        "their checked cardinality sum equals DecisionHead.reserved_subject_count",
        "RefusedNoEffect therefore burns its ordinal",
        "exact length is checked_add(106,checked_add(L,checked_mul(member_count,21261)))",
        "L=byte_len(FRAME_BYTES(the exact page table))",
        "double-framing or omitting the table frame is invalid",
        "mutually biject under the generated cross-product",
    ] {
        assert!(
            wire_correction.contains(exact),
            "terminal wire/accounting correction lost proof obligation '{exact}'"
        );
    }

    for exact in [
        "CapacityAuthorityRosterMatrixBytesV4 uses the 39 exact ASCII bytes I13_CAPACITY_AUTHORITY_ROSTER_MATRIX_V4 plus NUL",
        "exactly 2609 bytes under 'org.frankensim.i13.capacity-authority-roster-matrix.v4'",
        "CapacityIndexReservationSubjectBytesV4 uses the distinct 41 exact ASCII bytes I13_CAPACITY_INDEX_RESERVATION_SUBJECT_V4 plus NUL",
        "exactly266 bytes under 'org.frankensim.i13.capacity-index-reservation-subject.v4'",
        "CapacityResponseLogHeadBytesV4 uses the 33 exact ASCII bytes I13_CAPACITY_RESPONSE_LOG_HEAD_V4 plus NUL",
        "CapacityResponseLogEntryBytesV4 uses the 34 exact ASCII bytes I13_CAPACITY_RESPONSE_LOG_ENTRY_V4 plus NUL",
        "exactly checked_add(272,payload_bytes)",
        "CapacityResponseLogAppendReceiptBytesV4 uses the 43 exact ASCII bytes I13_CAPACITY_RESPONSE_LOG_APPEND_RECEIPT_V4 plus NUL",
        "exactly257 bytes",
        "CapacityTransactionWalEntryBytesV4 uses the distinct 37 exact ASCII bytes I13_CAPACITY_TRANSACTION_WAL_ENTRY_V4 plus NUL",
        "exactly checked_add(128,record_bytes) under 'org.frankensim.i13.capacity-transaction-wal.entry.v4'",
        "CapacityTransactionWalAppendReceiptBytesV4 uses the 46 exact ASCII bytes I13_CAPACITY_TRANSACTION_WAL_APPEND_RECEIPT_V4 plus NUL",
        "exactly260 bytes",
        "CapacityMutationCommitIntentBytesV4 uses the 38 exact ASCII bytes I13_CAPACITY_MUTATION_COMMIT_INTENT_V4 plus NUL",
        "exactly checked_add(370,checked_mul(component_count,138))",
        "CapacityTransactionReceiptBytesV4 uses the 35 exact ASCII bytes I13_CAPACITY_TRANSACTION_RECEIPT_V4 plus NUL",
        "exactly359 bytes",
        "CapacityMutationPrecommitEvidenceManifestBytesV4",
        "exactly checked_add(137,checked_mul(object_count,94)) under 'org.frankensim.i13.capacity-mutation-precommit-manifest.v4'",
        "CapacityMutationEvidenceManifestBytesV4",
        "exactly checked_add(136,checked_mul(object_count,94)) under 'org.frankensim.i13.capacity-mutation-evidence-manifest.v4'",
        "CapacityIndexCopyCoordinatorHeadBytesV4 uses the 38 exact ASCII bytes I13_CAPACITY_INDEX_COPY_COORDINATOR_V4 plus NUL",
        "exactly272 bytes",
        "CapacityAuthorityBindingMatrixBytesV4",
        "exactly checked_add(75,checked_mul(row_count,141)) under 'org.frankensim.i13.capacity-authority-binding-matrix.v4'",
        "ControlEvidenceProducerInventoryBytesV4",
        "exactly checked_add(77,checked_mul(row_count,163)) under 'org.frankensim.i13.control-evidence-producer-inventory.v4'",
        "ControlEvidenceStagingHeadBytesV4",
        "exactly215 bytes",
        "ControlEvidenceStagingInsertionPreparedBytesV4",
        "exactly460 bytes",
        "ControlEvidenceStagingSealPreparedBytesV4",
        "exactly354 bytes under 'org.frankensim.i13.control-evidence-staging-seal-prepared.v4'",
        "ControlEvidenceCutoffBytesV4",
        "exactly225 bytes",
    ] {
        assert!(
            base_v4.contains(exact),
            "base V4 authority lost nonsuperseded schema lock '{exact}'"
        );
    }

    for (authority_clause, schema, boundary, exact_size_or_formula, domain) in [
        (
            allocation_correction,
            "CapacityIndexDecisionHeadBytesV4",
            "It requires used_bytes",
            "exactly228 bytes",
            "org.frankensim.i13.capacity-index-decision-head.v4",
        ),
        (
            allocation_correction,
            "CapacityIndexAllocationRequestBytesV4",
            "CapacityIndexLocalDecisionBytesV4",
            "exactly426 bytes",
            "org.frankensim.i13.capacity-index-allocation-request.v4",
        ),
        (
            allocation_correction,
            "CapacityIndexLocalDecisionBytesV4",
            "CapacityIndexAllocationReceiptBytesV4",
            "exactly395 bytes",
            "org.frankensim.i13.capacity-index-local-decision.v4",
        ),
        (
            allocation_correction,
            "CapacityIndexAllocationReceiptBytesV4",
            "Every field is normative",
            "exactly636 bytes",
            "org.frankensim.i13.capacity-index-allocation-receipt.v4",
        ),
        (
            allocation_correction,
            "CapacityIndexDecisionEntryBytesV4",
            "CapacityIndexDecisionAppendReceiptBytesV4",
            "checked_add(328,payload_bytes)",
            "org.frankensim.i13.capacity-index-decision-entry.v4",
        ),
        (
            allocation_correction,
            "CapacityIndexDecisionAppendReceiptBytesV4",
            "Every append checks",
            "exactly335 bytes",
            "org.frankensim.i13.capacity-index-decision-append-receipt.v4",
        ),
        (
            allocation_correction,
            "CapacityIndexAllocationEvidenceBundleBytesV4",
            "A clean Reserved trace",
            "U64LE(total_framed_evidence_bytes)",
            "org.frankensim.i13.capacity-index-allocation-evidence-bundle.v4",
        ),
        (
            coordinator_correction,
            "CapacityTransactionChainStepBytesV4",
            "A fresh coordinator CAS",
            "exactly346 bytes",
            "org.frankensim.i13.capacity-transaction-chain-step.v4",
        ),
        (
            coordinator_correction,
            "CapacityMutationPreparedPayloadBytesV4",
            "ControlEvidenceStagingProofBundleV4",
            "checked_add(140,primary_bytes)",
            "org.frankensim.i13.capacity-mutation-prepared-payload.v4",
        ),
        (
            coordinator_correction,
            "ControlEvidenceStagingProofBundleV4",
            "ControlEvidenceStagingArtifactBytesV4",
            "at most16716 bytes",
            "org.frankensim.i13.control-evidence-staging-proof-bundle.v4",
        ),
        (
            coordinator_correction,
            "ControlEvidenceStagingArtifactBytesV4",
            "Each generation row",
            "declared attempt,epoch,arm,final Open head,member count and total physical bytes fields",
            "org.frankensim.i13.control-evidence-staging-artifact.v4",
        ),
        (
            coordinator_correction,
            "ControlEvidenceStagingLogPageRowBytesV4",
            "The page table is",
            "checked_add(159,allocation_authority_set_bytes)",
            "org.frankensim.i13.control-evidence-staging-log-page.v4",
        ),
    ] {
        let schema_scope = authority_clause
            .split_once(schema)
            .unwrap_or_else(|| panic!("terminal V4 authority lost schema {schema}"))
            .1
            .split_once(boundary)
            .unwrap_or_else(|| panic!("terminal V4 schema {schema} lost boundary {boundary}"))
            .0;
        assert!(
            schema_scope.contains(exact_size_or_formula),
            "terminal V4 schema {schema} lost size/formula {exact_size_or_formula}"
        );
        assert!(
            schema_scope.contains(domain),
            "terminal V4 schema {schema} lost domain {domain}"
        );
    }

    for (authority_clause, schema, semantic_lock) in [
        (
            base_v4,
            "CapacityAuthorityBindingRowBytesV4",
            "prepared-only forbidden=1",
        ),
        (
            base_v4,
            "ControlEvidenceProducerInventoryBytesV4",
            "kind-ordered unique FRAME_BYTES rows",
        ),
        (
            base_v4,
            "ControlEvidenceStagingInsertionPreparedBytesV4",
            "the 460-byte projection is never called a receipt or treated as committed",
        ),
        (
            coordinator_correction,
            "GlobalClosedInsertionProofBytesV4",
            "at most123062",
        ),
    ] {
        let schema_start = authority_clause
            .find(schema)
            .unwrap_or_else(|| panic!("terminal V4 authority lost schema {schema}"));
        let schema_scope = authority_clause[schema_start..]
            .chars()
            .take(1_500)
            .collect::<String>();
        assert!(
            schema_scope.contains(semantic_lock),
            "terminal V4 schema {schema} lost semantic lock {semantic_lock}"
        );
    }

    for (header, expected_len) in [
        (b"I13_CAPACITY_INDEX_ALLOCATION_REQUEST_V4\0".as_slice(), 41),
        (b"I13_CAPACITY_INDEX_ALLOCATION_RECEIPT_V4\0".as_slice(), 41),
        (
            b"I13_CAPACITY_INDEX_ALLOCATION_EVIDENCE_BUNDLE_V4\0".as_slice(),
            49,
        ),
        (
            b"I13_CAPACITY_MUTATION_PREPARED_PAYLOAD_V4\0".as_slice(),
            42,
        ),
        (
            b"I13_CONTROL_EVIDENCE_STAGING_PROOF_BUNDLE_V4\0".as_slice(),
            45,
        ),
        (b"I13_CONTROL_EVIDENCE_STAGING_ARTIFACT_V4\0".as_slice(), 41),
        (b"I13_CONTROL_EVIDENCE_STAGING_LOG_PAGE_V4\0".as_slice(), 41),
    ] {
        assert_eq!(header.len(), expected_len);
    }

    // Exact wire sums are independent witnesses to the prose totals. A field
    // can no longer disappear while a later sentence happens to retain the
    // same decimal number.
    assert_eq!(36 + 2 * 32 + 1 + 8 + 3 * 2 + 2 * 8 + 3 * 32 + 1, 228);
    assert_eq!(41 + 2 * 8 + 11 * 32 + 1 + 2 * 8, 426);
    assert_eq!(
        37 + 2 * 8 + 2 * 32 + 1 + 2 + 5 * 32 + 1 + 2 + 3 * 32 + 2 * 8,
        395
    );
    assert_eq!(
        41 + 2 * 8 + 16 * 32 + 2 + 2 * 8 + 2 * 8 + 2 * 8 + 1 + 2 * 8,
        636
    );
    assert_eq!(
        38 + 2 * 8 + 3 * 32 + 4 * 32 + 4 * 8 + 4 * 2 + 2 * 8 + 1,
        335
    );
    assert_eq!(49 + 2 * 8 + 2 + 2 + 8, 77);
    assert_eq!(32 + 2 + 32 + 32 + 32 + 2, 132);
    assert_eq!(42 + 2 * 8 + 32 + 2 + 32 + 8 + 8, 140);
    assert_eq!(45 + 2 * 8 + 1 + 1 + 8 + 2 * 8 + 8_348 + 8_280 + 1, 16_716);
    assert_eq!(41 + 2 * 8 + 1 + 2 + 8 + 2 + 8 + 8 + 8 + 2 * 32 + 1, 159);
    assert_eq!(41 + 2 * 8 + 1 + 32 + 2 * 8, 106);

    assert_eq!(24_576u64, 6 * 4_096);
    assert_eq!(24_576u64 * 32_768, 805_306_368);
    assert_eq!(4_096u64 * (262_144 + 131_072), 1_610_612_736);
    assert_eq!(805_306_368u64 + 1_610_612_736, 2_415_919_104);
    assert_eq!(328 + 395, 723);
    assert_eq!(328 + 428, 756);

    let clean_raw_objects = 266
        + 2 * 426
        + 2 * 636
        + 4 * 180
        + 8 * 228
        + 8 * 272
        + 2 * 8_330
        + 4 * 723
        + 2 * 756
        + 6 * 335;
    assert_eq!(clean_raw_objects, 30_184);
    assert_eq!(clean_raw_objects + 39 * 8, 30_496);
    assert_eq!(30_496 + 77, 30_573);
    assert_eq!(30_573 + 215, 30_788);
    assert_eq!(30_573 + 849, 31_422);

    let reconciled_raw_objects = 266
        + 2 * 426
        + 2 * 636
        + 4 * 180
        + 14 * 228
        + 14 * 272
        + 2 * 8_330
        + 8 * 723
        + 4 * 756
        + 12 * 335;
    assert_eq!(reconciled_raw_objects, 39_598);
    assert_eq!(reconciled_raw_objects + 63 * 8, 40_102);
    assert_eq!(40_102 + 77, 40_179);
    assert_eq!(40_179 + 215, 40_394);
    assert_eq!(40_179 + 849, 41_028);
    assert_eq!(40_394 + 272, 40_666);

    assert_eq!(61 + 63 * 299, 18_898);
    assert_eq!(18_898 + 423, 19_321);
    assert_eq!(19_321 + 133, 19_454);
    assert_eq!(39 + 2 * 8 + 8 * 32 + 2 + 2 * 8 + 2 * 8 + 1, 346);
    assert_eq!(660 + 16, 676);
    assert_eq!(3_283 + 16, 3_299);
    assert_eq!(456 + 19_454 + 29_211 + 65_536, 114_657);
    assert_eq!(165 + 114_657, 114_822);
    assert_eq!(8_240 + 114_822, 123_062);
    assert_eq!(8_607 + 123_062, 131_669);
    assert_eq!(106 + 225 + 3 * 58 + 131 + 333 + 65_536, 66_505);
    assert_eq!(8 + (8 + 116) + (8 + 460) + (8 + 20_653), 21_261);
    assert_eq!(21_261 + 58, 21_319);
    assert_eq!(1_298 + 131, 1_429);
    assert_eq!(121 + 402 + 121 + 1_429 + 121 + 295, 2_489);

    let nmax = 1u64;
    let lmax = 8u64;
    let fixed_and_indexed = 21_319u64
        .checked_mul(nmax)
        .and_then(|bytes| bytes.checked_add(lmax))
        .and_then(|bytes| bytes.checked_add(20_547))
        .and_then(|bytes| bytes.checked_add(66_505))
        .expect("representable control bound");
    let smax = 67_042_864u64
        .checked_sub(fixed_and_indexed)
        .expect("nonnegative equality witness");
    let equality = smax
        .checked_add(fixed_and_indexed)
        .expect("representable equality witness");
    assert_eq!(equality, 67_042_864);
    assert_eq!(equality.checked_add(1), Some(67_042_865));

    for obligation in [
        "exact-byte/domain KAT twins",
        "empty/first/max/max+1 maps,logs,decision slots and closed-evidence slots",
        "crash at every local-index/response/WAL/coordinator boundary",
        "one-sided repair",
        "two-PreEffectAbsent rejection",
        "exact-key replay/conflict",
        "partial allocation/orphan tombstone",
        "full durability-projection generation replay and sealed insertion refusal",
        "Until those proofs pass,the feature stays [M]",
        "no transaction,closure,expiry or all-and-only certificate is promoted",
    ] {
        assert!(
            gauntlet_obligations.contains(obligation),
            "terminal V4 Gauntlet lost no-claim/test obligation '{obligation}'"
        );
    }

    assert!(policy.contains(
        "FINAL_V4_FINITE_CAPACITY_WAL_AND_CLOSED_EVIDENCE_AUTHORITY are part of the final fresh authority set"
    ));
    assert!(policy.contains(
        "FINAL_V4_DURABLE_ALLOCATION_JOURNAL_CORRECTION,FINAL_V4_COORDINATOR_AND_STAGING_CAPACITY_CORRECTION and FINAL_V4_WIRE_AND_ACCOUNTING_DISAMBIGUATION are terminal coequal members of the fresh authority set"
    ));
}

fn sha256_for_i13_source_lock(input: &[u8]) -> [u8; 32] {
    const INITIAL: [u32; 8] = [
        0x6a09_e667,
        0xbb67_ae85,
        0x3c6e_f372,
        0xa54f_f53a,
        0x510e_527f,
        0x9b05_688c,
        0x1f83_d9ab,
        0x5be0_cd19,
    ];
    const ROUND: [u32; 64] = [
        0x428a_2f98,
        0x7137_4491,
        0xb5c0_fbcf,
        0xe9b5_dba5,
        0x3956_c25b,
        0x59f1_11f1,
        0x923f_82a4,
        0xab1c_5ed5,
        0xd807_aa98,
        0x1283_5b01,
        0x2431_85be,
        0x550c_7dc3,
        0x72be_5d74,
        0x80de_b1fe,
        0x9bdc_06a7,
        0xc19b_f174,
        0xe49b_69c1,
        0xefbe_4786,
        0x0fc1_9dc6,
        0x240c_a1cc,
        0x2de9_2c6f,
        0x4a74_84aa,
        0x5cb0_a9dc,
        0x76f9_88da,
        0x983e_5152,
        0xa831_c66d,
        0xb003_27c8,
        0xbf59_7fc7,
        0xc6e0_0bf3,
        0xd5a7_9147,
        0x06ca_6351,
        0x1429_2967,
        0x27b7_0a85,
        0x2e1b_2138,
        0x4d2c_6dfc,
        0x5338_0d13,
        0x650a_7354,
        0x766a_0abb,
        0x81c2_c92e,
        0x9272_2c85,
        0xa2bf_e8a1,
        0xa81a_664b,
        0xc24b_8b70,
        0xc76c_51a3,
        0xd192_e819,
        0xd699_0624,
        0xf40e_3585,
        0x106a_a070,
        0x19a4_c116,
        0x1e37_6c08,
        0x2748_774c,
        0x34b0_bcb5,
        0x391c_0cb3,
        0x4ed8_aa4a,
        0x5b9c_ca4f,
        0x682e_6ff3,
        0x748f_82ee,
        0x78a5_636f,
        0x84c8_7814,
        0x8cc7_0208,
        0x90be_fffa,
        0xa450_6ceb,
        0xbef9_a3f7,
        0xc671_78f2,
    ];

    let bit_len = u64::try_from(input.len())
        .expect("I13 source length fits u64")
        .checked_mul(8)
        .expect("I13 source bit length fits u64");
    let mut padded = input.to_vec();
    padded.push(0x80);
    while padded.len() % 64 != 56 {
        padded.push(0);
    }
    padded.extend_from_slice(&bit_len.to_be_bytes());

    let mut state = INITIAL;
    for chunk in padded.chunks_exact(64) {
        let mut words = [0u32; 64];
        for (index, word) in words[..16].iter_mut().enumerate() {
            let offset = index * 4;
            *word = u32::from_be_bytes(
                chunk[offset..offset + 4]
                    .try_into()
                    .expect("SHA-256 word is four bytes"),
            );
        }
        for index in 16..64 {
            let s0 = words[index - 15].rotate_right(7)
                ^ words[index - 15].rotate_right(18)
                ^ (words[index - 15] >> 3);
            let s1 = words[index - 2].rotate_right(17)
                ^ words[index - 2].rotate_right(19)
                ^ (words[index - 2] >> 10);
            words[index] = words[index - 16]
                .wrapping_add(s0)
                .wrapping_add(words[index - 7])
                .wrapping_add(s1);
        }

        let [mut a, mut b, mut c, mut d, mut e, mut f, mut g, mut h] = state;
        for index in 0..64 {
            let choice = (e & f) ^ ((!e) & g);
            let majority = (a & b) ^ (a & c) ^ (b & c);
            let sigma0 = a.rotate_right(2) ^ a.rotate_right(13) ^ a.rotate_right(22);
            let sigma1 = e.rotate_right(6) ^ e.rotate_right(11) ^ e.rotate_right(25);
            let temporary1 = h
                .wrapping_add(sigma1)
                .wrapping_add(choice)
                .wrapping_add(ROUND[index])
                .wrapping_add(words[index]);
            let temporary2 = sigma0.wrapping_add(majority);
            h = g;
            g = f;
            f = e;
            e = d.wrapping_add(temporary1);
            d = c;
            c = b;
            b = a;
            a = temporary1.wrapping_add(temporary2);
        }
        state[0] = state[0].wrapping_add(a);
        state[1] = state[1].wrapping_add(b);
        state[2] = state[2].wrapping_add(c);
        state[3] = state[3].wrapping_add(d);
        state[4] = state[4].wrapping_add(e);
        state[5] = state[5].wrapping_add(f);
        state[6] = state[6].wrapping_add(g);
        state[7] = state[7].wrapping_add(h);
    }

    let mut digest = [0u8; 32];
    for (index, word) in state.into_iter().enumerate() {
        digest[index * 4..index * 4 + 4].copy_from_slice(&word.to_be_bytes());
    }
    digest
}

fn lower_hex(bytes: &[u8]) -> String {
    const DIGITS: &[u8; 16] = b"0123456789abcdef";
    let mut encoded = String::with_capacity(bytes.len() * 2);
    for byte in bytes {
        encoded.push(DIGITS[usize::from(*byte >> 4)] as char);
        encoded.push(DIGITS[usize::from(*byte & 0x0f)] as char);
    }
    encoded
}

fn source_clause<'a>(source: &'a str, start: &str, end: &str) -> &'a str {
    source
        .split_once(start)
        .unwrap_or_else(|| panic!("I13 source lost terminal clause {start}"))
        .1
        .split_once(end)
        .unwrap_or_else(|| panic!("I13 terminal clause {start} lost boundary {end}"))
        .0
}

#[test]
fn i13_terminal_source_snapshot_and_authority_clauses_are_exact() {
    const SOURCE: &[u8] = include_bytes!("../src/i13.rs");
    const EXPECTED_SHA256: &str =
        "ee92bd9742a47bc2e7bb5c2d87267036a91e7c3e5e7e58a5c3303f6af48e0bd6";

    assert_eq!(
        lower_hex(&sha256_for_i13_source_lock(b"abc")),
        "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad",
        "source-lock SHA-256 implementation lost its independent KAT"
    );
    assert_eq!(SOURCE.len(), 1_025_512);
    assert_eq!(
        lower_hex(&sha256_for_i13_source_lock(SOURCE)),
        EXPECTED_SHA256,
        "I13 authority source changed; review every locked byte and update deliberately"
    );

    let source = std::str::from_utf8(SOURCE).expect("I13 source is canonical UTF-8");
    for terminal_clause in [
        "FINAL_V4_CAPABILITY_SIGNATURE_AND_STORAGE_CLOSURE",
        "FINAL_V4_EXCEPTION_ATOMIC_REPLAY_AND_CAPACITY_CLOSURE",
        "FINAL_V4_RESOLVER_HISTORY_AND_MAXIMUM_CLOSURE",
        "FINAL_V4_DIGEST_DOMAIN_REGISTRY_CLOSURE",
        "CAPABILITY_AND_SIGNATURE_DOMAINS_V4",
        "QUIESCENCE_STORAGE_AND_AVAILABILITY_DOMAINS_V4",
        "EXCEPTION_AND_RECOVERY_DOMAINS_V4",
        "SCHEMA_TO_DOMAIN_ALIAS_V4",
        "CANONICAL_ZERO_AND_RESOLUTION_RULES_V4",
    ] {
        assert!(
            source.contains(terminal_clause),
            "I13 source lost terminal authority clause {terminal_clause}"
        );
    }

    let capability = source_clause(
        source,
        "FINAL_V4_CAPABILITY_SIGNATURE_AND_STORAGE_CLOSURE:",
        "FINAL_V4_EXCEPTION_ATOMIC_REPLAY_AND_CAPACITY_CLOSURE:",
    );
    for lock in [
        "only the272-byte CopyAtomicHead mutates copy allocation state",
        "stable capacity_subject_slot_ordinal assigned from the frozen admitted subject set before execution,never by a speculative held ordinal",
        "Each subject owns a dense depth-4 vector of exactly16 such rows",
        "S exactly16 and therefore6571 bytes",
        "S<16 is forbidden",
        "digest then derives the transaction-chain/latest-record fields and final service head",
        "The record never contains that final head and is not signed by the same intent layer",
        "Later proofs,records and committed heads bind prior full signed objects but are never signed backward by the same layer",
        "Empty padding layers are forbidden",
        "Authority resolves byte-identically through the roster and every receipt",
        "Its unsigned projection zeroes each current-layer signature root,then the one quiescence layer signs the aggregate and the finalized set installs that root",
        "Record digest then derives chain/latest and the final head; the record contains no final head",
        "The final signature layer signs the unsigned508-byte AvailabilityFenceEvidence,not this enclosing bundle",
        "no SignatureObjectSet is embedded inside it",
        "Bmax=601+434+2981+10474+7998+2500+1576=26564 bytes",
    ] {
        assert!(
            capability.contains(lock),
            "capability/storage authority lost semantic lock '{lock}'"
        );
    }

    let exception = source_clause(
        source,
        "FINAL_V4_EXCEPTION_ATOMIC_REPLAY_AND_CAPACITY_CLOSURE:",
        "FINAL_V4_RESOLVER_HISTORY_AND_MAXIMUM_CLOSURE:",
    );
    for lock in [
        "receipt_count=absence_count<=4096 and chain_count<=3*4096",
        "Slot and PreEffect-receipt vectors are distinct depth-12 trees with distinct nonzero empty/leaf/depth-node/root domains",
        "both paths must reconstruct canonical empty leaves and the two old inner roots in the retained old core",
        "The receipt leaf does not yet exist",
        "Aggregate counts/bytes are not asserted by this proof",
        "FRAME_BYTES(the exact old core) opens them against an exact272-byte expected CopyAtomicHead",
        "It never contains the committed CopyAtomicHead",
        "excludes the receipt-dependent new roots",
        "the receipt digest becomes the new receipt-vector leaf and the dual witness derives its successor root",
        "Phase3/6 preserve that vector",
        "the final core uses the receipt-vector/global-chain successors and only then derives the sole CopyAtomicHead successor",
        "one other-copy phase3 exception is Empty and both later phase6 exceptions are Present",
        "Every Direct phase3 and phase6 history arm consumes the corresponding artifact",
        "exception arms consume their append artifact",
        "At most one copy can be recovered in a Reserved attempt",
        "the zero-durable double-absence case issues no recovery capability and has no effect authority",
        "Staged losers are drained and reclaimed only after proving no committed copy/service/store head references them",
        "retries cannot grow unbounded histories",
        "The stage lease is recoverable and becomes sticky after the first external service commit",
    ] {
        assert!(
            exception.contains(lock),
            "exception atomic-replay authority lost semantic lock '{lock}'"
        );
    }

    let resolver = source_clause(
        source,
        "FINAL_V4_RESOLVER_HISTORY_AND_MAXIMUM_CLOSURE:",
        "FINAL_V4_DIGEST_DOMAIN_REGISTRY_CLOSURE:",
    );
    for lock in [
        "request/static registries -> recoverable stage lease -> optional expiry signature -> signed capability intent",
        "committed capability heads -> signed quiescence aggregate -> signed storage intent",
        "committed storage head -> signed508-byte availability evidence",
        "PreEffect candidate/envelope with old core and dual vacancy -> slot proof -> candidate nonrecursive copy state -> conditional receipt",
        "receipt-vector/global-chain successors -> final core/root -> the sole CopyAtomicHead CAS",
        "Only the winning final CAS grants copy-state authority",
        "PreEffect successor equations are exact",
        "An unavailable phase3/6 append preserves decision,index,orphan,reservation,absence and receipt-vector projections",
        "MAX_RESERVED_ALLOCATION_EVIDENCE_BUNDLE_BYTES_V4=198055 with at most74 items",
        "Until those proofs are green this complete authority remains [M]",
    ] {
        assert!(
            resolver.contains(lock),
            "resolver/maximum authority lost semantic lock '{lock}'"
        );
    }
}

#[test]
fn i13_terminal_capability_storage_and_availability_wire_sums_are_locked() {
    let source = include_str!("../src/i13.rs");
    let capability = source_clause(
        source,
        "FINAL_V4_CAPABILITY_SIGNATURE_AND_STORAGE_CLOSURE:",
        "FINAL_V4_EXCEPTION_ATOMIC_REPLAY_AND_CAPACITY_CLOSURE:",
    );

    for (header, expected_len) in [
        (
            b"I13_CAPACITY_INDEX_CAPABILITY_SUBJECT_HEAD_V4\0".as_slice(),
            46,
        ),
        (
            b"I13_CAPACITY_INDEX_CAPABILITY_SERVICE_HEAD_V4\0".as_slice(),
            46,
        ),
        (
            b"I13_CAPACITY_INDEX_CAPABILITY_FENCE_INTENT_V4\0".as_slice(),
            46,
        ),
        (b"I13_CAPACITY_INDEX_EXPIRY_QUORUM_V4\0".as_slice(), 36),
        (
            b"I13_CAPACITY_INDEX_CAPABILITY_FENCE_COMMITV4\0".as_slice(),
            45,
        ),
        (b"I13_CAPACITY_INDEX_WRITER_ROSTER_V4\0".as_slice(), 36),
        (
            b"I13_CAPACITY_INDEX_SIGNATURE_OBJECT_SET_V4\0".as_slice(),
            43,
        ),
        (
            b"I13_CAPACITY_INDEX_CANONICAL_SIGNATURE_SET_V4\0".as_slice(),
            46,
        ),
        (
            b"I13_CAPACITY_INDEX_SIGNATURE_LAYER_SET_V4\0".as_slice(),
            42,
        ),
        (
            b"I13_CAPACITY_INDEX_WRITER_QUIESCENCE_RECEIPT_V4\0".as_slice(),
            48,
        ),
        (b"I13_CAPACITY_INDEX_QUIESCENCE_SET_V4\0".as_slice(), 37),
        (b"I13_CAPACITY_INDEX_STORAGE_FENCE_HEAD_V4\0".as_slice(), 41),
        (
            b"I13_CAPACITY_INDEX_STORAGE_FENCE_INTENT_V4\0".as_slice(),
            43,
        ),
        (
            b"I13_CAPACITY_INDEX_STORAGE_FENCE_COMMIT_V4\0".as_slice(),
            43,
        ),
        (
            b"I13_CAPACITY_INDEX_AVAILABILITY_BUNDLE_V4\0".as_slice(),
            42,
        ),
    ] {
        assert_eq!(header.len(), expected_len);
    }

    // Independent field arithmetic prevents a stable decimal total from
    // masking a silently added, removed, or widened field.
    assert_eq!(2 + 2 * 32 + 1 + 4 * 8 + 3 * 32, 195);
    assert_eq!(
        46 + 2 * 8 + 3 * 32 + 1 + 2 * 8 + 1 + 2 * 32 + 2 + 32 + 2 * 2 + 2 * 32 + 1,
        343
    );
    assert_eq!(46 + 32 + 8 + 32 + 2 + 3 * 32 + 1, 217);
    assert_eq!(2 + 4 * 32 + 12 * 32 + 1, 515);
    assert_eq!(
        46 + 2 * 8 + 3 * 32 + 1 + 2 + 5 * 32 + 1 + 2 + 2 * 32 + 2 * 8 + 32 + 1,
        437
    );
    assert_eq!(75 + 406 * 16, 6_571);
    assert_eq!(2 + 32 + 2 + 2 + 1 + 32 + 1, 72);
    assert_eq!(3 * 32 + 2 * 8, 112);
    assert_eq!(
        36 + 2 * 8 + 3 * 32 + 1 + 3 * 32 + 2 + 6 * 8 + 1 + 32 + 1,
        329
    );
    assert_eq!(329 + 120 * 7, 1_169);
    assert_eq!(
        45 + 2 * 8 + 3 * 32 + 1 + 6 * 32 + 6 * 8 + 2 + 1 + 2 + 5 * 32 + 2 * 8 + 2,
        581
    );
    assert_eq!(445 + 2 * 351 + 523 + 589 + 2 * 225, 2_709);
    assert_eq!(1 + 2_709 + 80, 2_790);
    assert_eq!(1 + 2_709 + 6_579, 9_289);
    assert_eq!(1 + 2_709 + 6_579 + 1_177, 10_466);

    assert_eq!(2 + 2 + 1 + 1 + 3 * 32 + 3 * 8 + 32, 158);
    assert_eq!(
        36 + 2 * 8 + 3 * 32 + 1 + 4 * 32 + 2 + 1 + 2 * 2 + 32 + 1,
        317
    );
    assert_eq!(317 + 166 * 16, 2_973);
    assert_eq!(43 + 2 * 8 + 3 * 32 + 1 + 1 + 2 + 32 + 1, 192);
    assert_eq!(192 + 43 * 32, 1_568);
    assert_eq!(1 + 4 + 64, 69);
    assert_eq!(46 + 2 * 8 + 2 * 32 + 1 + 1 + 2 * 32 + 3, 195);
    assert_eq!(195 + 77 * 8, 811);
    assert_eq!(1 + 1 + 1 + 2 * 32, 67);
    assert_eq!(42 + 2 * 8 + 3 * 32 + 1 + 1 + 1 + 32, 189);
    assert_eq!(189 + 75 * 4, 489);
    assert_eq!(189 + 75 * 5, 564);
    assert_eq!(4 * 819 + (200 * 4 + 43 * 32) + 497, 5_949);
    assert_eq!(5 * 819 + (200 * 5 + 43 * 32) + 572, 7_043);

    assert_eq!(
        48 + 2 * 8 + 3 * 32 + 1 + 2 + 2 * 32 + 3 * 8 + 2 * 32 + 2 * 8 + 4 * 32 + 2 * 8 + 2,
        477
    );
    assert_eq!(37 + 2 * 8 + 2 * 32 + 3 * 32 + 1 + 2 * 2 + 8 + 4, 230);
    assert_eq!(230 + (8 + 477) * 16, 7_990);

    assert_eq!(2 * 32 + 3 + 3 * 8 + 3 * 32, 187);
    assert_eq!(41 + 32 + 8 + 32 + 2 + 3 * 32 + 1, 212);
    assert_eq!(2 + 2 * (8 + 187) + 2 * 32 + 12 * 32 + 1, 841);
    assert_eq!(43 + 2 * 8 + 3 * 32 + 1 + 2 + 10 * 32 + 2 * 8 + 1, 495);
    assert_eq!(2 * (8 + 187) + 32, 422);
    assert_eq!(43 + 7 * 32 + 2 * 8 + 2 * (8 + 187) + 2 * 8 + 2, 691);
    assert_eq!(1 + (8 + 495) + 2 * (8 + 212) + (8 + 841) + (8 + 691), 2_492);

    assert_eq!(
        42 + 16 * 32 + 2 * 8 + 3 + 2 + 2 + 2 + 1 + 4 + 2 * 8 + 1,
        601
    );
    assert_eq!(601 + 434 + 2_981 + 10_474 + 7_998 + 2_500 + 1_576, 26_564);

    for (schema, size, domain) in [
        (
            "CapacityIndexStorageStateMemberBytesV4",
            "exactly187 bytes",
            "org.frankensim.i13.capacity-index-storage-state-member.v4",
        ),
        (
            "CapacityIndexStorageFenceHeadBytesV4",
            "exactly212 bytes",
            "org.frankensim.i13.capacity-index-storage-fence-head.v4",
        ),
        (
            "Present->Present update proof",
            "exactly841 bytes",
            "org.frankensim.i13.capacity-index-storage-fence-update-proof.v4",
        ),
        (
            "CapacityIndexStorageFenceIntentBytesV4",
            "exactly495 bytes",
            "org.frankensim.i13.capacity-index-storage-fence-intent.v4",
        ),
        (
            "CapacityIndexStableExtentProjectionBytesV4",
            "exactly422 bytes",
            "org.frankensim.i13.capacity-index-stable-extent-projection.v4",
        ),
        (
            "CapacityIndexStorageFenceCommitRecordBytesV4",
            "exactly691 bytes",
            "org.frankensim.i13.capacity-index-storage-fence-commit.v4",
        ),
        (
            "StorageDispositionBytesV4",
            "exactly2492 bytes",
            "org.frankensim.i13.capacity-index-storage-disposition.v4",
        ),
    ] {
        assert!(capability.contains(schema), "terminal source lost {schema}");
        assert!(
            capability.contains(size),
            "terminal source lost {schema} {size}"
        );
        assert!(
            source.contains(domain),
            "terminal source lost domain {domain}"
        );
    }
}

#[test]
fn i13_terminal_exception_graph_replay_capacity_and_maxima_are_locked() {
    let source = include_str!("../src/i13.rs");
    let exception = source_clause(
        source,
        "FINAL_V4_EXCEPTION_ATOMIC_REPLAY_AND_CAPACITY_CLOSURE:",
        "FINAL_V4_RESOLVER_HISTORY_AND_MAXIMUM_CLOSURE:",
    );
    let resolver = source_clause(
        source,
        "FINAL_V4_RESOLVER_HISTORY_AND_MAXIMUM_CLOSURE:",
        "FINAL_V4_DIGEST_DOMAIN_REGISTRY_CLOSURE:",
    );

    for (header, expected_len) in [
        (
            b"I13_CAPACITY_INDEX_EXCEPTION_COMMIT_RECEIPT_V4\0".as_slice(),
            47,
        ),
        (b"I13_CAPACITY_EVIDENCE_CHARGE_RECEIPT_V4\0".as_slice(), 40),
    ] {
        assert_eq!(header.len(), expected_len);
    }

    assert_eq!(5 * 32 + 5 * 2 + 8, 178);
    assert_eq!(32 + 32 + 1 + 8 + 178 + 32, 283);
    assert_eq!(2 + 12 * 32 + 12 * 32 + 1, 771);
    assert_eq!(2 + 1 + (8 + 180) + 2 * 32 + 12 * 32 + 1, 640);
    assert_eq!(640 + 8 + 180, 828);
    assert_eq!(559 + 26_564, 27_123);
    assert_eq!(9_130 + 27_123 + (8 + 178) + (8 + 771), 37_218);
    assert_eq!(
        47 + 1 + 2 * 8 + 2 * 32 + 1 + 2 + 4 * 32 + 2 * 8 + 2 * 2 + 2 * 4 + 2 * 32 + 32 + 1,
        384
    );
    assert_eq!(2 * 32 + 1 + 2 * 2 + 2 + 3 * 32 + 32, 199);
    assert_eq!((8 + 37_218) + (8 + 640) + (8 + 384), 38_266);
    assert_eq!(1_728 + 26_564, 28_292);
    assert_eq!((8 + 28_292) + (8 + 640) + (8 + 384) + (8 + 178), 29_526);
    assert_eq!((8 + 28_292) + (8 + 828) + (8 + 384) + (8 + 178), 29_714);
    assert_eq!(29_526 + 8, 29_534);
    assert_eq!(29_714 + 8, 29_722);
    assert_eq!(29_534 + 2 * 29_722, 88_978);
    assert_eq!(2 + 1 + 1 + 2 * 32 + 12 * 32 + 1, 453);
    assert_eq!(453 + (8 + 180), 641);
    assert_eq!((8 + 641) + (8 + 178), 835);
    assert_eq!(835 + 8, 843);

    assert_eq!((8 + 37_218) + 2 * (8 + 28_292), 93_826);
    assert_eq!(131_072 - 93_826, 37_246);
    assert_eq!(4_096u64 * 131_072, 536_870_912);
    assert_eq!(
        1_040 + 3 * 1_414 + 4 * 1_115 + 3 * 7_051 + 3 * 395 + 64,
        32_144
    );
    assert_eq!(32_768 - 32_144, 624);
    assert_eq!(40 + 2 * 8 + 3 * 32 + 1 + 2 + 2 * 32 + 2 * 4 + 5 * 32, 387);
    assert_eq!(4_096u64 * 32_768, 134_217_728);
    assert_eq!(2_415_919_104u64 + 536_870_912 + 134_217_728, 3_087_007_744);

    assert_eq!(
        40_179 + (1_784 - 723) + (8 + 228) + (8 + 272) + (8 + 38_266),
        80_030
    );
    assert_eq!(80_030 + 88_978 + 4 * (8 + 7_043) + (8 + 835), 198_055);
    assert_eq!(63 + 3 + 3 + 4 + 1, 74);
    assert_eq!(198_055 + 215, 198_270);
    assert_eq!(198_055 + 849, 198_904);
    assert_eq!(198_270 + 272, 198_542);
    assert_eq!(198_904 + 159, 199_063);
    assert_eq!(199_063 + 8, 199_071);
    assert_eq!(262_144 - 198_055, 64_089);

    for lock in [
        "exactly771 bytes",
        "exactly640 bytes for Empty and828 for Present",
        "at most9130+27123+186+779=37218 bytes",
        "exactly384 bytes",
        "exactly199 bytes",
        "at most38266 bytes",
        "at most29526 bytes for Empty and29714 for Present",
        "its outer evidence item is29534 or29722",
        "retained sum is29534+2*29722=88978",
        "at most835 bytes and its outer item843",
        "93826,leaving37246 bytes",
        "32144<=32768,leaving624 bytes",
        "exactly387 bytes and FRAME_BYTES size395",
        "complete per-copy substrate is2415919104+536870912+4096*32768=3087007744 bytes",
    ] {
        assert!(
            exception.contains(lock),
            "exception/capacity authority lost arithmetic lock '{lock}'"
        );
    }
    for lock in [
        "for80030 bytes/66 items",
        "Add one Empty and two Present exception items88978 bytes/3 items",
        "four FRAME(SignatureMaterialArtifact7043)=28204 bytes/4 items",
        "one FRAME(SnapshotArtifact835)=843 bytes/1 item",
        "MAX_RESERVED_ALLOCATION_EVIDENCE_BUNDLE_BYTES_V4=198055 with at most74 items",
        "pair=198270,allocation authority set=198904,pair response-entry wrapper=198542,raw page row=199063 and framed page row=199071",
        "MAX_TERMINAL_NO_EFFECT_ALLOCATION_EVIDENCE_BUNDLE_BYTES_V4 remains262144 and is effect-ineligible",
        "A double-recovery construction is invalid",
        "no unpriced second recovered chain exists",
    ] {
        assert!(
            resolver.contains(lock),
            "resolver maximum authority lost arithmetic lock '{lock}'"
        );
    }

    for edge in [
        (
            "for_each_capacity_phase3_coordinator_head_a_v4",
            "for_each_capacity_phase3_slot_snapshot_artifact_a_v4",
        ),
        (
            "capacity_index_exception_evidence_store_v4",
            "for_each_capacity_phase3_slot_snapshot_artifact_a_v4",
        ),
        (
            "for_each_capacity_phase3_slot_snapshot_artifact_a_v4",
            "for_each_capacity_recovered_phase3_direct_a_v4",
        ),
        (
            "for_each_capacity_phase3_coordinator_head_b_v4",
            "for_each_capacity_phase3_slot_snapshot_artifact_b_v4",
        ),
        (
            "capacity_index_exception_evidence_store_v4",
            "for_each_capacity_phase3_slot_snapshot_artifact_b_v4",
        ),
        (
            "for_each_capacity_phase3_slot_snapshot_artifact_b_v4",
            "for_each_capacity_recovered_phase3_direct_b_v4",
        ),
        (
            "for_each_capacity_phase6_required_local_mask_v4",
            "for_each_capacity_phase6_slot_snapshot_proof_set_v4",
        ),
        (
            "for_each_capacity_phase6_head_set_v4",
            "for_each_capacity_phase6_slot_snapshot_proof_set_v4",
        ),
        (
            "capacity_index_exception_evidence_store_v4",
            "for_each_capacity_phase6_slot_snapshot_proof_set_v4",
        ),
        (
            "for_each_capacity_phase6_slot_snapshot_proof_set_v4",
            "for_each_capacity_phase6_exception_history_partition_v4",
        ),
    ] {
        assert!(
            I13_FRESH_V2_AUTHORITY_DAG_EDGES.contains(&edge),
            "I13 machine DAG lost direct-slot snapshot edge {edge:?}"
        );
        assert!(
            !I13_FRESH_V2_AUTHORITY_DAG_EDGES.contains(&(edge.1, edge.0)),
            "I13 machine DAG introduced snapshot backedge {edge:?}"
        );
    }

    for domain in [
        "org.frankensim.i13.capacity-index-availability-preimage-bundle.v4",
        "org.frankensim.i13.capacity-index-capability-fence-commit.v4",
        "org.frankensim.i13.capacity-index-storage-fence-commit.v4",
        "org.frankensim.i13.capacity-index-durable-absence.v4",
        "org.frankensim.i13.capacity-index-signature-material.v4",
        "org.frankensim.i13.capacity-index-pre-effect-dual-vacancy-witness.v4",
        "org.frankensim.i13.capacity-index-exception-slot-update-proof.v4",
        "org.frankensim.i13.capacity-index-canonical-signature-member.v4",
        "org.frankensim.i13.capacity-index-coordinator-core.v4",
        "org.frankensim.i13.capacity-index-candidate-coordinator-state.v4",
        "org.frankensim.i13.capacity-index-exception-commit-receipt.v4",
        "org.frankensim.i13.capacity-index-exception-chain-step.v4",
        "org.frankensim.i13.capacity-index-direct-slot-snapshot-proof.v4",
        "org.frankensim.i13.capacity-index-direct-slot-snapshot-artifact.v4",
        "org.frankensim.i13.capacity-index-evidence-store-charge-receipt.v4",
    ] {
        assert!(
            source.contains(domain),
            "I13 source lost terminal digest domain {domain}"
        );
    }
    assert!(source.contains(
        "Every other named XBytesV4 maps by its exact schema stem to the same literal listed above. No two semantically distinct schema digests share a domain"
    ));
}
