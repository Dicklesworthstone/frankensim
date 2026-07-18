//! Focused I14 multirung EMC/harness VerificationManifest conformance.
//!
//! These tests pin authored authority, partitions, execution leaves, theorem
//! ratchets, mutation sensitivity, canonical identity, and amendment blast
//! radius. They execute no physics solver, laboratory campaign, standard,
//! proof kernel, or exhaustive search and therefore mint no such evidence.

use fs_vmanifest::{
    Ambition, CampaignTier, ClaimPolarity, ClaimSpec, ContentHash, FixturePin, FixtureSource,
    FreezeRefusal, GauntletTier, I14_CANCELLATION_CARD_V2_KAT_HEX,
    I14_CANONICAL_TERMINAL_RESULT_V1_KAT_HEX, I14_CANONICAL_TERMINAL_RESULT_V2_KAT_HEX,
    I14_DRAIN_TRIGGER_ENCODING_V2_KAT_HEX, I14_INFRASTRUCTURE_FAILURE_ONSET_ENCODING_V2_KAT_HEX,
    I14_MAX_CANCELLATION_REQUESTS_V1, I14_MAX_OBSERVER_TILES_V1, I14_MAX_SCOPE_ANCESTRY_V1,
    I14_MAX_TERMINAL_ARBITRATION_PAIRS_V2, I14_MAX_TERMINAL_BOUNDARIES_V2,
    I14_MAX_WATCHDOG_OBSERVATIONS_V1, I14_TELEMETRY_ENVELOPE_V1_KAT_HEX,
    I14_TELEMETRY_ENVELOPE_V2_KAT_HEX, I14_TERMINAL_PREFIX_V2_KAT_HEX,
    I14_TERMINAL_STATUS_TABLE_V1_DIGEST_HEX, I14_TERMINAL_STATUS_TABLE_V1_TUPLES,
    I14_WATCHDOG_RAW_TRACE_V2_KAT_HEX, I14ArtifactCategoryV1, I14CancellationCardInputV2,
    I14CancellationCardRefusalV2, I14CancellationObservationV1, I14CancellationRequestStateV1,
    I14CancellationRequestV1, I14CancellationTierV2, I14CanonicalResultRefusalV1,
    I14CanonicalResultRefusalV2, I14CanonicalTerminalResultInputV1,
    I14CanonicalTerminalResultInputV2, I14ClaimAdjudication, I14DomainApplicability,
    I14DrainTriggerV2, I14EvidenceCompleteness, I14EvidenceIntegrity, I14ExecutionDisposition,
    I14ExternalHeartbeatCoverageV2, I14InfrastructureFailureOnsetV2,
    I14InfrastructureFailureSourceV2, I14InputValidity, I14LateEventTailV2,
    I14LifecycleCauseClassV2, I14LifecycleFailureV2, I14LifecycleRefusalV2, I14OperationalSupport,
    I14ReceiptValidity, I14RetentionClassV1, I14SpawnFrontierEvidenceV2,
    I14TelemetryEnvelopeInputV1, I14TelemetryEnvelopeInputV2, I14TelemetryEnvelopeRefusalV1,
    I14TelemetryEnvelopeRefusalV2, I14TerminalBoundaryDecisionV1, I14TerminalBoundaryRecordV2,
    I14TerminalBoundaryTraceV2, I14TerminalBoundaryV1, I14TerminalCauseRefusalV1,
    I14TerminalLifecycleTraceV2, I14TerminalStatusV1, I14TerminalTraceOutcomeV2,
    I14TerminalTraceRefusalV2, I14TilePollCoverageV2, I14TimedLogicalEventV2,
    I14TotalResourceUnitV2, I14WatchdogCoverageV2, I14WatchdogObservationKindV1,
    I14WatchdogObservationV1, ManifestDraft, Partition, ToleranceSemantics,
    i14_admit_cancellation_card_v2, i14_canonical_terminal_result_digest_v1,
    i14_canonical_terminal_result_digest_v2, i14_canonical_terminal_result_v1,
    i14_canonical_terminal_result_v2, i14_draft, i14_drain_trigger_encoding_v2,
    i14_evaluate_terminal_status_v1, i14_infrastructure_failure_onset_encoding_v2,
    i14_retention_rule_v1, i14_select_first_terminal_boundary_v2, i14_select_terminal_boundary_v1,
    i14_telemetry_envelope_digest_v1, i14_telemetry_envelope_digest_v2,
    i14_terminal_status_table_digest_v1, i14_watchdog_raw_trace_digest_v2,
};
use std::collections::{BTreeMap, BTreeSet};

const POLICY: &str = "i14-campaign-policy-v1";
const ACCEPTANCE_POLICY: &str = "i14-acceptance-arithmetic-policy-v1";
const THEOREM_POLICY: &str = "i14-theorem-formalization-policy-v1";
const EM_CONVENTION_CARD: &str = "i14-em-convention-card-v1";

const SOLID_CLAIMS: &[&str] = &[
    "i14-fullwave-problem-convention-admission",
    "i14-harnessgraph-identity-connectivity",
    "i14-synthetic-ap242-adapter-mechanics",
];

const FRONTIER_CLAIMS: &[&str] = &[
    "i14-bearing-current-hybrid-path",
    "i14-core-fidelity-crosswalk-routing",
    "i14-emc-uq-inference-mechanics",
    "i14-fixed-regime-adjoint-closure",
    "i14-governed-laboratory-emc-validation",
    "i14-governed-standards-crosswalk",
    "i14-ground-bond-shield-current-closure",
    "i14-immunity-victim-mode-ledger",
    "i14-mtl-passive-causal-propagation",
    "i14-mtl-rlgc-operator-admission",
    "i14-peec-extraction-power-mor",
    "i14-switching-source-probe-semantics",
];

const MOONSHOT_CLAIMS: &[&str] = &[
    "i14-certified-fidelity-descent-theorem",
    "i14-cover-refinement-naturality-theorem",
    "i14-emc-safety-case-integration",
    "i14-exterior-bem-formulation-correctness",
    "i14-fullwave-feec-stability-energy",
    "i14-governed-emc-reliability-validation",
    "i14-hypercohomology-obstruction-localization-theorem",
    "i14-kyp-sheaf-passivity-bridge-theorem",
    "i14-maximal-counterexample-search",
    "i14-maxwell-fmm-acceleration-envelope",
    "i14-passive-causal-sheaf-composition-theorem",
    "i14-production-bearing-population-reliability",
    "i14-robust-mitigation-heldout",
];

type ClaimAuthority = (
    &'static str,
    Ambition,
    GauntletTier,
    &'static str,
    &'static str,
    ToleranceSemantics,
    &'static str,
);

const CLAIM_AUTHORITY: &[ClaimAuthority] = &[
    (
        "i14-harnessgraph-identity-connectivity",
        Ambition::Solid,
        GauntletTier::G0,
        "native_harness_identity_connectivity_and_orientation_verdict",
        "bit",
        ToleranceSemantics::Exact,
        "i14.oracle.harnessgraph.v1 at fs-vmanifest-oracles/i14/harnessgraph.rs::reconstruct_and_compare_incidence",
    ),
    (
        "i14-synthetic-ap242-adapter-mechanics",
        Ambition::Solid,
        GauntletTier::G0,
        "synthetic_ap242_subset_adapter_identity_transform_loss_and_refusal_verdict",
        "bit",
        ToleranceSemantics::Exact,
        "i14.oracle.ap242_adapter.v1 at fs-vmanifest-oracles/i14/ap242_adapter.rs::reconstruct_occurrences_and_losses",
    ),
    (
        "i14-fullwave-problem-convention-admission",
        Ambition::Solid,
        GauntletTier::G0,
        "fullwave_semantic_admission_and_explicit_crosswalk_verdict",
        "bit",
        ToleranceSemantics::Exact,
        "i14.oracle.fullwave_schema.v1 at fs-vmanifest-oracles/i14/fullwave_schema.rs::check_conventions",
    ),
    (
        "i14-mtl-rlgc-operator-admission",
        Ambition::Frontier,
        GauntletTier::G2,
        "maximum_preregistered_normalized_rlgc_structure_source_and_validity_defect",
        "1",
        ToleranceSemantics::Interval { lo: 0.0, hi: 1.0 },
        "i14.oracle.rlgc.v1 at fs-vmanifest-oracles/i14/rlgc.rs::check_operator_and_quotient",
    ),
    (
        "i14-mtl-passive-causal-propagation",
        Ambition::Frontier,
        GauntletTier::G2,
        "maximum_preregistered_normalized_mtl_wave_port_power_and_time_frequency_discrepancy",
        "1",
        ToleranceSemantics::Interval { lo: 0.0, hi: 1.0 },
        "i14.oracle.mtl_wave.v1 at fs-vmanifest-oracles/i14/mtl_wave.rs::solve_analytic_and_balance",
    ),
    (
        "i14-peec-extraction-power-mor",
        Ambition::Frontier,
        GauntletTier::G2,
        "maximum_preregistered_normalized_peec_charge_gauge_power_passivity_and_reduction_defect",
        "1",
        ToleranceSemantics::Interval { lo: 0.0, hi: 1.0 },
        "i14.oracle.peec.v1 at fs-vmanifest-oracles/i14/peec.rs::extract_stamp_and_adjudicate",
    ),
    (
        "i14-ground-bond-shield-current-closure",
        Ambition::Frontier,
        GauntletTier::G2,
        "maximum_preregistered_normalized_return_current_shield_transfer_and_loss_ownership_defect",
        "1",
        ToleranceSemantics::Interval { lo: 0.0, hi: 1.0 },
        "i14.oracle.ground_shield.v1 at fs-vmanifest-oracles/i14/ground_shield.rs::reconcile_paths",
    ),
    (
        "i14-bearing-current-hybrid-path",
        Ambition::Frontier,
        GauntletTier::G2,
        "maximum_preregistered_normalized_shaft_voltage_bearing_current_event_and_energy_defect",
        "1",
        ToleranceSemantics::Interval { lo: 0.0, hi: 1.0 },
        "i14.oracle.bearing_current.v1 at fs-vmanifest-oracles/i14/bearing_current.rs::check_hybrid_path",
    ),
    (
        "i14-core-fidelity-crosswalk-routing",
        Ambition::Frontier,
        GauntletTier::G3,
        "core_route_admit_escalate_unknown_and_crosswalk_power_verdict",
        "bit",
        ToleranceSemantics::Exact,
        "i14.oracle.core_router.v1 at fs-vmanifest-oracles/i14/core_router.rs::adjudicate_crosswalks",
    ),
    (
        "i14-switching-source-probe-semantics",
        Ambition::Frontier,
        GauntletTier::G2,
        "maximum_preregistered_normalized_source_spectrum_probe_and_power_chain_defect",
        "1",
        ToleranceSemantics::Interval { lo: 0.0, hi: 1.0 },
        "i14.oracle.source_probe.v1 at fs-vmanifest-oracles/i14/source_probe.rs::reconstruct_and_balance",
    ),
    (
        "i14-immunity-victim-mode-ledger",
        Ambition::Frontier,
        GauntletTier::G3,
        "source_to_victim_power_event_upset_and_recovery_verdict",
        "bit",
        ToleranceSemantics::Exact,
        "i14.oracle.victim_chain.v1 at fs-vmanifest-oracles/i14/victim_chain.rs::replay_interventions",
    ),
    (
        "i14-fixed-regime-adjoint-closure",
        Ambition::Frontier,
        GauntletTier::G3,
        "maximum_preregistered_normalized_tangent_adjoint_taylor_and_independent_derivative_defect",
        "1",
        ToleranceSemantics::Interval { lo: 0.0, hi: 1.0 },
        "i14.oracle.adjoint.v1 at fs-vmanifest-oracles/i14/adjoint.rs::check_fixed_regime_derivative",
    ),
    (
        "i14-emc-uq-inference-mechanics",
        Ambition::Frontier,
        GauntletTier::G3,
        "uq_coverage_optional_stopping_tail_and_escalation_verdict",
        "bit",
        ToleranceSemantics::Exact,
        "i14.oracle.uq.v1 at fs-vmanifest-oracles/i14/uq.rs::audit_sampling_and_coverage",
    ),
    (
        "i14-fullwave-feec-stability-energy",
        Ambition::Moonshot,
        GauntletTier::G2,
        "maximum_preregistered_normalized_fullwave_stability_gauss_spectral_power_and_dispersion_defect",
        "1",
        ToleranceSemantics::Interval { lo: 0.0, hi: 1.0 },
        "i14.oracle.fullwave_feec.v1 at fs-vmanifest-oracles/i14/fullwave_feec.rs::adjudicate_stability_and_balance",
    ),
    (
        "i14-exterior-bem-formulation-correctness",
        Ambition::Moonshot,
        GauntletTier::G2,
        "maximum_preregistered_normalized_bem_formulation_trace_power_and_dense_qoi_defect",
        "1",
        ToleranceSemantics::Interval { lo: 0.0, hi: 1.0 },
        "i14.oracle.maxwell_bem.v1 at fs-vmanifest-oracles/i14/maxwell_bem.rs::dense_trace_and_farfield_check",
    ),
    (
        "i14-maxwell-fmm-acceleration-envelope",
        Ambition::Moonshot,
        GauntletTier::G2,
        "maximum_preregistered_normalized_fmm_dense_matvec_solved_qoi_and_performance_envelope_defect",
        "1",
        ToleranceSemantics::Interval { lo: 0.0, hi: 1.0 },
        "i14.oracle.maxwell_fmm.v1 at fs-vmanifest-oracles/i14/maxwell_fmm.rs::dense_envelope_and_crossover_check",
    ),
    (
        "i14-certified-fidelity-descent-theorem",
        Ambition::Moonshot,
        GauntletTier::G3,
        "independent_fidelity_descent_theorem_and_runtime_majorant_verdict",
        "bit",
        ToleranceSemantics::Exact,
        "i14.oracle.fidelity_descent.v1 composite proofs/i14/FidelityDescent.lean::emcFidelityDescent plus fs-vmanifest-oracles/i14/fidelity_descent.rs::check_runtime_majorant",
    ),
    (
        "i14-robust-mitigation-heldout",
        Ambition::Moonshot,
        GauntletTier::G3,
        "blind_robust_mitigation_improvement_all_guards_and_independent_reproduction_verdict",
        "bit",
        ToleranceSemantics::Exact,
        "i14.oracle.robust_mitigation.v1 at fs-vmanifest-oracles/i14/robust_mitigation.rs::blind_reconstruct_and_compare",
    ),
    (
        "i14-emc-safety-case-integration",
        Ambition::Moonshot,
        GauntletTier::G3,
        "scoped_emc_evidence_to_hazard_traceability_and_no_laundering_verdict",
        "bit",
        ToleranceSemantics::Exact,
        "i14.oracle.assurance.v1 at fs-vmanifest-oracles/i14/assurance.rs::audit_hazard_edges",
    ),
    (
        "i14-governed-standards-crosswalk",
        Ambition::Frontier,
        GauntletTier::G3,
        "scoped_exact_edition_clause_crosswalk_and_loss_accounting_verdict",
        "bit",
        ToleranceSemantics::Exact,
        "i14.oracle.standards_crosswalk.v1 at fs-vmanifest-oracles/i14/standards_crosswalk.rs::reconstruct_clause_and_adapter_edges",
    ),
    (
        "i14-governed-laboratory-emc-validation",
        Ambition::Frontier,
        GauntletTier::G3,
        "scoped_asbuilt_laboratory_emc_validation_vector_and_integrity_verdict",
        "bit",
        ToleranceSemantics::Exact,
        "i14.oracle.lab_validation.v1 at fs-vmanifest-oracles/i14/lab_validation.rs::reconstruct_calibration_and_qoi_vector",
    ),
    (
        "i14-production-bearing-population-reliability",
        Ambition::Moonshot,
        GauntletTier::G3,
        "production_bearing_population_reliability_coverage_and_integrity_verdict",
        "bit",
        ToleranceSemantics::Exact,
        "i14.oracle.bearing_population.v1 at fs-vmanifest-oracles/i14/bearing_population.rs::audit_frame_events_and_coverage",
    ),
    (
        "i14-governed-emc-reliability-validation",
        Ambition::Moonshot,
        GauntletTier::G3,
        "governed_emc_population_reliability_coverage_and_integrity_verdict",
        "bit",
        ToleranceSemantics::Exact,
        "i14.oracle.emc_reliability.v1 at fs-vmanifest-oracles/i14/emc_reliability.rs::audit_frame_events_and_anytime_coverage",
    ),
    (
        "i14-passive-causal-sheaf-composition-theorem",
        Ambition::Moonshot,
        GauntletTier::G3,
        "independent_passive_causal_sheaf_cosheaf_composition_checker_verdict",
        "bit",
        ToleranceSemantics::Exact,
        "i14.oracle.passive_composition.v1 composite proofs/i14/PassiveCausalComposition.lean::assemblyPassivity plus fs-vmanifest-oracles/i14/passive_composition.rs::check_runtime_diagrams",
    ),
    (
        "i14-hypercohomology-obstruction-localization-theorem",
        Ambition::Moonshot,
        GauntletTier::G3,
        "hypercohomology_obstruction_gluing_witness_and_minimal_support_checker_verdict",
        "bit",
        ToleranceSemantics::Exact,
        "i14.oracle.hypercohomology_obstruction.v1 composite proofs/i14/HypercohomologyObstruction.lean::relativeGluingObstruction plus fs-vmanifest-oracles/i14/hypercohomology.rs::check_total_complex_and_support",
    ),
    (
        "i14-cover-refinement-naturality-theorem",
        Ambition::Moonshot,
        GauntletTier::G3,
        "cover_refinement_section_obstruction_power_and_passivity_naturality_verdict",
        "bit",
        ToleranceSemantics::Exact,
        "i14.oracle.cover_refinement.v1 composite proofs/i14/CoverRefinement.lean::passivityNaturality plus fs-vmanifest-oracles/i14/cover_refinement.rs::check_comparison_diagrams",
    ),
    (
        "i14-kyp-sheaf-passivity-bridge-theorem",
        Ambition::Moonshot,
        GauntletTier::G3,
        "generalized_pr_descriptor_kyp_storage_supply_and_local_relation_bridge_verdict",
        "bit",
        ToleranceSemantics::Exact,
        "i14.oracle.kyp_sheaf_bridge.v1 composite proofs/i14/KypSheafBridge.lean::generalizedPrToLocalDissipativity plus fs-vmanifest-oracles/i14/kyp_bridge.rs::check_descriptor_witness",
    ),
    (
        "i14-maximal-counterexample-search",
        Ambition::Moonshot,
        GauntletTier::G3,
        "exact_nonvacuity_coverage_zero_genuine_countermodel_and_integrity_verdict",
        "bit",
        ToleranceSemantics::Exact,
        "i14.oracle.maximal_falsifier.v1 at fs-vmanifest-oracles/i14/maximal_falsifier.rs::verify_enumeration_membership_and_minimize",
    ),
];

const UNIT_CASES: [&str; 9] = [
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
        "i14-harnessgraph-core",
        CampaignTier::Core,
        &["i14-harnessgraph-identity-connectivity"],
    ),
    (
        "i14-ap242-adapter-core",
        CampaignTier::Core,
        &["i14-synthetic-ap242-adapter-mechanics"],
    ),
    (
        "i14-rlgc-operator-core",
        CampaignTier::Core,
        &["i14-mtl-rlgc-operator-admission"],
    ),
    (
        "i14-mtl-propagation-core",
        CampaignTier::Core,
        &["i14-mtl-passive-causal-propagation"],
    ),
    (
        "i14-peec-network",
        CampaignTier::Core,
        &["i14-peec-extraction-power-mor"],
    ),
    (
        "i14-ground-shield-core",
        CampaignTier::Core,
        &["i14-ground-bond-shield-current-closure"],
    ),
    (
        "i14-bearing-current-core",
        CampaignTier::Core,
        &["i14-bearing-current-hybrid-path"],
    ),
    (
        "i14-source-probe-core",
        CampaignTier::Core,
        &["i14-switching-source-probe-semantics"],
    ),
    (
        "i14-victim-mode-core",
        CampaignTier::Core,
        &["i14-immunity-victim-mode-ledger"],
    ),
    (
        "i14-fidelity-routing-core",
        CampaignTier::Core,
        &["i14-core-fidelity-crosswalk-routing"],
    ),
    (
        "i14-fixed-regime-adjoint-core",
        CampaignTier::Core,
        &["i14-fixed-regime-adjoint-closure"],
    ),
    (
        "i14-uq-inference-core",
        CampaignTier::Core,
        &["i14-emc-uq-inference-mechanics"],
    ),
    (
        "i14-fullwave-schema-core",
        CampaignTier::Core,
        &["i14-fullwave-problem-convention-admission"],
    ),
    (
        "i14-fullwave-feec-max",
        CampaignTier::Max,
        &["i14-fullwave-feec-stability-energy"],
    ),
    (
        "i14-bem-formulation-max",
        CampaignTier::Max,
        &["i14-exterior-bem-formulation-correctness"],
    ),
    (
        "i14-fmm-acceleration-max",
        CampaignTier::Max,
        &["i14-maxwell-fmm-acceleration-envelope"],
    ),
    (
        "i14-fidelity-descent-theorem-max",
        CampaignTier::Max,
        &["i14-certified-fidelity-descent-theorem"],
    ),
    (
        "i14-robust-mitigation-max",
        CampaignTier::Max,
        &["i14-robust-mitigation-heldout"],
    ),
    (
        "i14-safety-case-integration-max",
        CampaignTier::Max,
        &["i14-emc-safety-case-integration"],
    ),
    (
        "i14-passive-composition-theorem-max",
        CampaignTier::Max,
        &["i14-passive-causal-sheaf-composition-theorem"],
    ),
    (
        "i14-hypercohomology-obstruction-max",
        CampaignTier::Max,
        &["i14-hypercohomology-obstruction-localization-theorem"],
    ),
    (
        "i14-cover-refinement-naturality-max",
        CampaignTier::Max,
        &["i14-cover-refinement-naturality-theorem"],
    ),
    (
        "i14-kyp-sheaf-bridge-max",
        CampaignTier::Max,
        &["i14-kyp-sheaf-passivity-bridge-theorem"],
    ),
    (
        "i14-maximal-falsifier-max",
        CampaignTier::Max,
        &["i14-maximal-counterexample-search"],
    ),
    (
        "i14-standards-crosswalk-max",
        CampaignTier::Max,
        &["i14-governed-standards-crosswalk"],
    ),
    (
        "i14-laboratory-validation-max",
        CampaignTier::Max,
        &["i14-governed-laboratory-emc-validation"],
    ),
    (
        "i14-bearing-population-validation-max",
        CampaignTier::Max,
        &["i14-production-bearing-population-reliability"],
    ),
    (
        "i14-emc-reliability-validation-max",
        CampaignTier::Max,
        &["i14-governed-emc-reliability-validation"],
    ),
];

const HOLDOUTS: &[(&str, CampaignTier, &str)] = &[
    (
        "i14-harnessgraph-core-holdout",
        CampaignTier::Core,
        "i14-harnessgraph-core",
    ),
    (
        "i14-ap242-adapter-core-holdout",
        CampaignTier::Core,
        "i14-ap242-adapter-core",
    ),
    (
        "i14-rlgc-operator-core-holdout",
        CampaignTier::Core,
        "i14-rlgc-operator-core",
    ),
    (
        "i14-mtl-bundle-core-holdout",
        CampaignTier::Core,
        "i14-mtl-propagation-core",
    ),
    (
        "i14-peec-mtl-crossover-core-holdout",
        CampaignTier::Core,
        "i14-fidelity-routing-core",
    ),
    (
        "i14-peec-network-core-holdout",
        CampaignTier::Core,
        "i14-peec-network",
    ),
    (
        "i14-shield-ground-core-holdout",
        CampaignTier::Core,
        "i14-ground-shield-core",
    ),
    (
        "i14-bearing-current-core-holdout",
        CampaignTier::Core,
        "i14-bearing-current-core",
    ),
    (
        "i14-source-probe-core-holdout",
        CampaignTier::Core,
        "i14-source-probe-core",
    ),
    (
        "i14-victim-upset-core-holdout",
        CampaignTier::Core,
        "i14-victim-mode-core",
    ),
    (
        "i14-adjoint-core-holdout",
        CampaignTier::Core,
        "i14-fixed-regime-adjoint-core",
    ),
    (
        "i14-uq-core-holdout",
        CampaignTier::Core,
        "i14-uq-inference-core",
    ),
    (
        "i14-fullwave-schema-core-holdout",
        CampaignTier::Core,
        "i14-fullwave-schema-core",
    ),
    (
        "i14-fullwave-max-holdout",
        CampaignTier::Max,
        "i14-fullwave-feec-max",
    ),
    (
        "i14-bem-formulation-max-holdout",
        CampaignTier::Max,
        "i14-bem-formulation-max",
    ),
    (
        "i14-fmm-acceleration-max-holdout",
        CampaignTier::Max,
        "i14-fmm-acceleration-max",
    ),
    (
        "i14-mitigation-max-holdout",
        CampaignTier::Max,
        "i14-robust-mitigation-max",
    ),
    (
        "i14-emc-reliability-max-holdout",
        CampaignTier::Max,
        "i14-emc-reliability-validation-max",
    ),
    (
        "i14-laboratory-validation-max-holdout",
        CampaignTier::Max,
        "i14-laboratory-validation-max",
    ),
    (
        "i14-bearing-population-max-holdout",
        CampaignTier::Max,
        "i14-bearing-population-validation-max",
    ),
];

const WAIVER_CONSUMERS: &[(&str, &[&str])] = &[
    (
        "i14-external-standards-edition-clause-pack",
        &["i14-standards-crosswalk-max"],
    ),
    (
        "i14-external-emc-laboratory-calibration-pack",
        &["i14-laboratory-validation-max"],
    ),
    (
        "i14-external-asbuilt-specimen-geometry-pack",
        &["i14-laboratory-validation-max"],
    ),
    (
        "i14-external-bearing-population-reliability-pack",
        &["i14-bearing-population-validation-max"],
    ),
    (
        "i14-external-bearing-population-metrology-pack",
        &["i14-bearing-population-validation-max"],
    ),
    (
        "i14-external-emc-reliability-population-pack",
        &["i14-emc-reliability-validation-max"],
    ),
    (
        "i14-external-blind-mitigation-custody-pack",
        &["i14-robust-mitigation-max"],
    ),
];

fn claim<'a>(draft: &'a ManifestDraft, id: &str) -> &'a ClaimSpec {
    draft
        .claims
        .iter()
        .find(|candidate| candidate.id == id)
        .unwrap_or_else(|| panic!("missing I14 claim '{id}'"))
}

fn authored_spec<'a>(draft: &'a ManifestDraft, id: &str) -> &'a str {
    let fixture = draft
        .fixtures
        .iter()
        .find(|candidate| candidate.id == id)
        .unwrap_or_else(|| panic!("missing I14 fixture '{id}'"));
    match fixture.source {
        FixtureSource::AuthoredSpec { spec } => spec,
        FixtureSource::External { .. } => panic!("I14 fixture '{id}' must be authored"),
    }
}

fn authority_ids(draft: &ManifestDraft) -> BTreeSet<&'static str> {
    draft
        .claims
        .iter()
        .map(|claim| claim.id)
        .chain(draft.obligations.iter().map(|row| row.leaf))
        .collect()
}

fn assert_exact_lifecycle_event_membership(draft: &ManifestDraft, leaf: &str, expected: &[&str]) {
    let row = draft
        .obligations
        .iter()
        .find(|row| row.leaf == leaf)
        .unwrap_or_else(|| panic!("missing I14 obligation '{leaf}'"));
    let domains = expected
        .iter()
        .map(|event| {
            event
                .split_once('.')
                .map(|(domain, _)| domain)
                .unwrap_or_else(|| panic!("{leaf} lifecycle event {event} has no domain prefix"))
        })
        .collect::<BTreeSet<_>>();
    let is_relevant = |event: &str| {
        domains.iter().any(|domain| {
            event.starts_with(domain) && event.as_bytes().get(domain.len()) == Some(&b'.')
        })
    };
    let actual = row
        .obs_events
        .iter()
        .copied()
        .filter(|event| is_relevant(event))
        .collect::<BTreeSet<_>>();
    let expected = expected.iter().copied().collect::<BTreeSet<_>>();
    assert_eq!(
        actual, expected,
        "{leaf} must carry exactly the frozen {domains:?} lifecycle vocabulary"
    );
    assert_eq!(
        row.obs_events
            .iter()
            .filter(|event| is_relevant(event))
            .count(),
        actual.len(),
        "{leaf} must not duplicate a {domains:?} lifecycle event"
    );
}

fn assert_ordered_contract_tokens(context: &str, text: &str, expected: &[&str]) {
    let mut cursor = 0;
    for token in expected {
        let relative = text[cursor..]
            .find(token)
            .unwrap_or_else(|| panic!("{context} omits ordered token {token}"));
        cursor += relative + token.len();
    }
}

fn i14_core_cancellation_card_input() -> I14CancellationCardInputV2 {
    I14CancellationCardInputV2 {
        tier: I14CancellationTierV2::Core,
        semantic_work_unit_digest: ContentHash([0x41; 32]),
        campaign_wall_budget_ns: 10_000_000_000,
        logical_memory_ceiling_bytes: 32 * 1_073_741_824,
        total_resource_unit: I14TotalResourceUnitV2::GraphTraceRecords,
        total_resource_ceiling: 4_096,
        maximum_resource_tile_quantum: 64,
        resource_budget_authority_digest: ContentHash([0x45; 32]),
        logical_partition_spec_digest: ContentHash([0x46; 32]),
        execution_environment_digest: ContentHash([0x47; 32]),
        logical_tile_response_bound_ns: 200_000_000,
        indivisible_item_response_bound_ns: 20_000_000,
        external_heartbeat_bound_ns: 25_000_000,
        external_child_catalog_digest: Some(ContentHash([0x48; 32])),
    }
}

fn i14_core_cancellation_card() -> fs_vmanifest::I14CancellationCardV2 {
    i14_admit_cancellation_card_v2(i14_core_cancellation_card_input())
        .expect("valid Core cancellation card")
}

fn i14_terminal_status(execution: I14ExecutionDisposition) -> I14TerminalStatusV1 {
    I14TerminalStatusV1 {
        execution,
        claim: I14ClaimAdjudication::Pending,
        completeness: I14EvidenceCompleteness::CompleteEvidence,
        integrity: I14EvidenceIntegrity::IntegrityVerified,
        input: I14InputValidity::WellFormedInput,
        domain: I14DomainApplicability::Admitted,
        support: I14OperationalSupport::SupportedOperation,
        receipt: I14ReceiptValidity::WellFormedReceipt,
    }
}

fn i14_infrastructure_onset(
    source: I14InfrastructureFailureSourceV2,
    logical_sequence: u64,
    monotonic_ns: u64,
    receipt_byte: u8,
) -> I14InfrastructureFailureOnsetV2 {
    I14InfrastructureFailureOnsetV2 {
        source,
        event: I14TimedLogicalEventV2 {
            logical_sequence,
            monotonic_ns,
        },
        verification_receipt_digest: ContentHash([receipt_byte; 32]),
    }
}

#[test]
#[allow(clippy::too_many_lines)]
fn i14_seed_freezes_with_exact_claim_and_execution_authority() {
    const EXACT_EXIT_PRECEDENCE: &str = "ReceiptValidity=MalformedReceipt or InputValidity=MalformedInput or EvidenceIntegrity=IntegrityFailed -> 60; ExecutionDisposition=InfrastructureFailed -> 70; ExecutionDisposition=Cancelled -> 20; ExecutionDisposition=TimedOut -> 21; ExecutionDisposition=BudgetExhausted -> 22; ClaimAdjudication=Failed -> 10; ClaimAdjudication=Refuted -> 40; DomainApplicability in {OutOfDomain,Indeterminate} or OperationalSupport=UnsupportedOperation -> 30; EvidenceCompleteness in {PartialEvidence,NoEvidence} -> 50; ClaimAdjudication=Unknown -> 30; 0 only for ExecutionDisposition=Completed with every requested ClaimAdjudication=Supported, EvidenceCompleteness=CompleteEvidence, EvidenceIntegrity=IntegrityVerified, InputValidity=WellFormedInput, DomainApplicability=Admitted, OperationalSupport=SupportedOperation and ReceiptValidity=WellFormedReceipt";
    let draft = i14_draft();
    assert_eq!(draft.initiative, "I14");
    assert_eq!(draft.version, 1);
    assert!(draft.title.contains("Multirung EMC/harness gate"));
    assert!(draft.explicits.versions.contains("fs-vmanifest schema v2"));
    assert!(draft.explicits.versions.contains(
        "I14 TerminalStatusTruthTable schema v1 with code-pinned domain-separated digest"
    ));
    assert!(
        draft
            .explicits
            .versions
            .contains("I14 legacy local scoped terminal-cause selector schema v1")
    );
    assert!(draft.explicits.versions.contains(
        "I14 legacy single-boundary canonical terminal-result schema v1 under org.frankensim.i14.terminal-result.v1"
    ));
    assert!(draft.explicits.versions.contains(
        "I14 legacy noncanonical telemetry-envelope schema v1 under org.frankensim.i14.telemetry-envelope.v1"
    ));
    for schema in [
        "I14 validated cancellation-card schema v2 under org.frankensim.i14.cancellation-card.v2",
        "I14 recomputed clock-free logical execution trace schema v2 under org.frankensim.i14.logical-execution-trace.v2",
        "I14 authoritative genesis-to-first-terminal trace schema v2 under org.frankensim.i14.terminal-trace.v2",
        "I14 authoritative canonical terminal-result schema v2 under org.frankensim.i14.terminal-result.v2",
        "I14 receipt-bound noncanonical telemetry-envelope schema v2 under org.frankensim.i14.telemetry-envelope.v2",
        "I14 complete raw-watchdog trace schema v2 under org.frankensim.i14.watchdog-raw-trace.v2",
        "I14 typed four-variant drain-trigger/lifecycle schema v2",
    ] {
        assert!(
            draft.explicits.versions.contains(schema),
            "missing {schema}"
        );
    }
    assert!(
        draft
            .explicits
            .versions
            .contains("I14 typed artifact-retention policy schema v1")
    );
    assert_eq!(draft.claims.len(), CLAIM_AUTHORITY.len());
    assert_eq!(draft.fixtures.len(), 53);
    assert_eq!(draft.obligations.len(), LEAF_MAP.len());
    assert_eq!(draft.waivers.len(), 7);

    let expected_lattice = [
        (Ambition::Solid, SOLID_CLAIMS),
        (Ambition::Frontier, FRONTIER_CLAIMS),
        (Ambition::Moonshot, MOONSHOT_CLAIMS),
    ];
    for (ambition, ids) in expected_lattice {
        let actual = draft
            .claims
            .iter()
            .filter(|claim| claim.ambition == ambition)
            .map(|claim| claim.id)
            .collect::<BTreeSet<_>>();
        assert_eq!(actual, ids.iter().copied().collect(), "{ambition:?}");
    }
    let refutations = draft
        .claims
        .iter()
        .filter(|claim| claim.polarity == ClaimPolarity::Refutation)
        .map(|claim| claim.id)
        .collect::<Vec<_>>();
    assert_eq!(refutations, ["i14-maximal-counterexample-search"]);

    for (id, ambition, tier, qoi, unit, tolerance, oracle) in CLAIM_AUTHORITY {
        let actual = claim(&draft, id);
        assert_eq!(actual.ambition, *ambition, "{id} ambition");
        assert_eq!(actual.evidence_tier, *tier, "{id} Gauntlet tier");
        assert_eq!(actual.qoi, *qoi, "{id} QoI");
        assert_eq!(actual.unit, *unit, "{id} unit");
        assert_eq!(actual.tolerance, *tolerance, "{id} tolerance");
        assert_eq!(actual.oracle.identity, *oracle, "{id} oracle");
        assert!(actual.oracle.independent, "{id} production-oracle reuse");
        assert!(!actual.oracle.tcb_overlap.trim().is_empty(), "{id} TCB");
        assert!(!actual.hypotheses.is_empty(), "{id} hypotheses");
        for (field, text) in [
            ("statement", actual.statement),
            ("activation", actual.activation),
            ("kill", actual.kill),
            ("fallback", actual.fallback),
            ("no_claim", actual.no_claim),
        ] {
            assert!(!text.trim().is_empty(), "{id} omits {field}");
        }
    }

    let fixture_ids = draft
        .fixtures
        .iter()
        .map(|pin| pin.id)
        .collect::<BTreeSet<_>>();
    let waiver_ids = draft
        .waivers
        .iter()
        .map(|waiver| waiver.subject)
        .collect::<BTreeSet<_>>();
    let mut coverage = BTreeMap::<&str, usize>::new();
    for (expected_leaf, expected_tier, expected_claims) in LEAF_MAP {
        let row = draft
            .obligations
            .iter()
            .find(|row| row.leaf == *expected_leaf)
            .unwrap_or_else(|| panic!("missing I14 leaf '{expected_leaf}'"));
        assert_eq!(row.tier, *expected_tier, "{} tier", row.leaf);
        assert_eq!(
            row.claims_covered.iter().copied().collect::<BTreeSet<_>>(),
            expected_claims.iter().copied().collect(),
            "{} mapping",
            row.leaf
        );
        for id in row.claims_covered {
            *coverage.entry(id).or_default() += 1;
        }
        assert_eq!(
            row.unit_cases.iter().copied().collect::<BTreeSet<_>>(),
            UNIT_CASES.into_iter().collect(),
            "{} must own all nine unit-case classes",
            row.leaf
        );
        assert!(row.decks.contains(&POLICY), "{} omits policy", row.leaf);
        assert!(
            row.decks.contains(&ACCEPTANCE_POLICY),
            "{} omits acceptance arithmetic",
            row.leaf
        );
        if !matches!(row.leaf, "i14-harnessgraph-core" | "i14-ap242-adapter-core") {
            assert!(
                row.decks.contains(&EM_CONVENTION_CARD),
                "{} omits the shared EM convention card",
                row.leaf
            );
        }
        for deck in row.decks {
            assert!(
                fixture_ids.contains(deck) || waiver_ids.contains(deck),
                "{} has orphan deck {deck}",
                row.leaf
            );
        }
        let suffix = row
            .leaf
            .strip_prefix("i14-")
            .expect("I14 leaf")
            .replace('-', "_");
        let entry = format!("scripts/e2e/leapfrog/i14_{suffix}.sh");
        assert_eq!(row.entry_point, entry, "{} entry", row.leaf);
        assert_eq!(
            row.dsr_lane,
            format!(
                "env FRANKENSIM_VMANIFEST_LEAF={} dsr quality --tool frankensim",
                row.leaf
            )
        );
        assert_eq!(
            row.replay_command,
            format!("{} --replay <artifact-id>", row.entry_point)
        );
        assert!(row.g0.contains("generators:"), "{} G0 generators", row.leaf);
        assert!(
            row.g0.contains("validity predicates:"),
            "{} G0 predicates",
            row.leaf
        );
        assert!(row.g0.contains("laws:"), "{} G0 laws", row.leaf);
        assert!(row.g0.contains("shrinkers:"), "{} G0 shrinkers", row.leaf);
        assert!(row.g0.contains("replay"), "{} G0 replay", row.leaf);
        assert!(row.g3_relations.len() >= 3, "{} G3 depth", row.leaf);
        for token in [
            "request cancellation",
            "drain",
            "finalize",
            "checkpoint",
            "resume",
            "fork",
            "Core total ceilings are 4096 graph/trace records",
            "16384 field/quadrature unknowns",
            "1024 search/tree nodes",
            "256 formal declarations",
            "Core logical poll tile contains at most 64 graph/trace records",
            "256 field/quadrature unknowns",
            "16 search/tree nodes",
            "or 4 formal declarations",
            "Core asynchronous watchdog polls at intervals <=25 ms",
            "without changing logical tile membership or order",
            "Core request-to-observation <=250 ms",
            "drain-trigger-to-drained <=2000 ms",
            "drained-to-finalized <=2000 ms",
            "at most 32 in-flight children",
            "Max total ceilings are 16384 graph/trace records",
            "65536 field/quadrature unknowns",
            "4096 search/tree nodes",
            "1024 formal declarations",
            "Max logical poll tile contains at most 256 graph/trace records",
            "1024 field/quadrature unknowns",
            "64 search/tree nodes",
            "or 16 formal declarations",
            "Max asynchronous watchdog polls at intervals <=100 ms",
            "Max request-to-observation <=1000 ms",
            "drain-trigger-to-drained <=8000 ms",
            "drained-to-finalized <=8000 ms",
            "at most 128 in-flight children",
            "Admission refuses any indivisible item or external heartbeat",
            "any logical tile whose admitted worst-case request-to-boundary bound exceeds its tier request-to-observation SLO",
            "execution.requested means job admission, not cancellation",
            "cancellation.requested binds request id",
            "cancellation.observed binds the same request id",
            "optional latest-completed-boundary ordinal",
            "whose logical sequence precedes the candidate boundary",
            "scope root occurs in that candidate's acyclic root-to-leaf scope ancestry",
            "request slices and observer-tile catalogs are canonicalized before validation",
            "malformed refusals are input-order invariant",
            "multiple relevant requests are totally ordered by their unique logical sequence",
            "I14_MAX_CANCELLATION_REQUESTS_V1=16384",
            "I14_MAX_SCOPE_ANCESTRY_V1=256",
            "I14_MAX_OBSERVER_TILES_V1=128",
            "i14_select_terminal_boundary_v1 validates those predicates",
            "V1 proves only local arbitration at a caller-supplied boundary and has zero promotion authority",
            "I14CancellationCardV2 admits Core, ordinary Max, or the explicit MaxTheoremFalsifier subtype",
            "hard logical-memory ceiling, exact count resource kind/ceiling/tile quantum",
            "resource authority, deterministic partition, execution-environment fingerprint and external-child catalog",
            "preserves the frozen 90-minute/18-hour/24-hour envelopes",
            "refuses logical-tile, indivisible-item or external-heartbeat response bounds wider than the frozen tier SLO",
            "exact 250 ms/1000 ms request-deadline delta",
            "tier-specific 32/128 observer and in-flight-child caps",
            "unique scope identities",
            "strict request-before-observation causality with strict-pre-boundary participation in the frozen cut",
            "admitted observing-tile identity",
            "observation's latest-completed boundary strictly precedes the candidate",
            "a retained observation later than the candidate must name a latest-completed boundary ordinal at least the candidate ordinal",
            "Calibrated monotonic timestamps must be nondecreasing in globally unique logical-sequence order across requests, observations and the candidate boundary",
            "request-to-observation is their timestamp delta",
            "drain-trigger-to-drained starts at the referenced on-time observation, first nanosecond after the missed inclusive request deadline, structurally admitted receipt-bound infrastructure-failure onset whose receipt HELM/ledger authenticates, or drain-start timestamp",
            "CancellationObserved, ObservationTimeoutDrain, InfrastructureFailure, or NonCancellationDrain respectively",
            "drained-to-finalized begins there and ends at execution.finalized",
            "An actual observation closes the spawn frontier at its logical event",
            "without one, request/deadline alone does not synthesize a logical cut and timeout/failure drain start closes it",
            "first terminal-eligible logical tile boundary in deterministic logical tile/event order",
            "i14_select_first_terminal_boundary_v2 requires a nonempty trace beginning at genesis ordinal 0",
            "I14_MAX_TERMINAL_BOUNDARIES_V2=4096 contiguous ordinals",
            "I14_MAX_TERMINAL_ARBITRATION_PAIRS_V2=1048576 boundary/request pairs",
            "no supplied record after the first Selected decision",
            "opaque nonterminal frontier certificate cannot be converted into a terminal result",
            "Temporal precedence across boundaries is distinct from same-boundary cause precedence",
            "makes that candidate ineligible until the same request is observed",
            "missed observation deadline selects TimedOut",
            "InfrastructureFailed > TimedOut > Cancelled > BudgetExhausted > Completed",
            "overall campaign wall-clock expiry",
            "only normal completion of all scheduled work",
            "completion exactly at the ceiling remains Completed",
            "raw timing fields never change logical tile membership or boundary ordinal and are excluded from canonical-result bytes",
            "when their recorded semantics select a different trigger or terminal cause, the canonical disposition and digest honestly change",
            "Late triggers are retained but never rewrite an already finalized boundary",
            "telemetry-only verification-receipt-bound late-event tail whose authentication is deferred to HELM/ledger",
            "I14TerminalLifecycleTraceV2 makes execution-started, drain-started, drained and finalized events mandatory",
            "I14DrainTriggerV2={CancellationObserved,ObservationTimeoutDrain,InfrastructureFailure,NonCancellationDrain}",
            "trigger is the earliest effective in-scope on-time observation, first nanosecond after a missed inclusive observation deadline or structurally admitted receipt-bound infrastructure-failure onset",
            "HELM/ledger authenticates the onset receipt",
            "Canonical bytes bind the derived variant plus request id or infrastructure-onset logical sequence",
            "Every terminal path requires I14SpawnFrontierEvidenceV2",
            "audit binds the unconditional child semantic root",
            "separate child raw root is telemetry-only",
            "independent child semantic-verification receipt identity is canonical",
            "Watchdog, tile-poll and heartbeat evidence likewise split canonical semantic roots from telemetry-only raw roots and independent semantic-verification receipt identities",
            "scheduler polls before item 0 and immediately before and after every logical tile",
            "Local code validates only structural/internal consistency, caps, trigger arithmetic and failure reflection",
            "HELM/ledger authenticates the corresponding receipts; local code only validates structural binding and consistency",
            "typed infrastructure-onset witness binds source={WatchdogCoverage,TilePollCoverage,ExternalHeartbeatCoverage,DescendantDrain,SpawnAfterFrontierClosure,Supervisor,Authentication,DrainProtocol,PublicationProtocol}",
            "causal tie rank second (Infrastructure=0, Observation=1, Timeout=2; distinct from wire tags)",
            "Timeout onset is derived, never caller-selected",
            "earliest first nanosecond outside the inclusive trigger-to-drained or drained-to-finalized cap",
            "timeout logical field is not a synthetic onset event",
            "exact first latch-boundary logical sequence",
            "first nanosecond after a missed observation deadline selects TimedOut",
            "effective time first, causal tie rank second",
            "causal logical sequence third and stable identity fourth",
            "I14_DRAIN_TRIGGER_ENCODING_V2_KAT_HEX pins all four exact trigger tags and U64LE payloads",
            "I14_INFRASTRUCTURE_FAILURE_ONSET_ENCODING_V2_KAT_HEX pins the present canonical presence/sequence/source/receipt form",
            "schema-valid projection itself has no promotion authority until the consuming HELM/ledger verifier authenticates",
            "Real watchdog, descendant-drain or spawn-after-frontier failure remains hashable only when selected InfrastructureFailed reports it",
            "real trigger-to-drained or drained-to-finalized deadline failure remains hashable only when selected TimedOut or higher-priority InfrastructureFailed reports it",
            "eight orthogonal axes exactly",
            "ExecutionDisposition={Completed,Cancelled,TimedOut,BudgetExhausted,InfrastructureFailed}",
            "ClaimAdjudication={Pending,Supported,Failed,Refuted,Unknown}",
            "EvidenceCompleteness={CompleteEvidence,PartialEvidence,NoEvidence}",
            "EvidenceIntegrity={IntegrityVerified,IntegrityFailed}",
            "InputValidity={WellFormedInput,MalformedInput}",
            "DomainApplicability={Admitted,OutOfDomain,Indeterminate}",
            "OperationalSupport={SupportedOperation,UnsupportedOperation}",
            "ReceiptValidity={WellFormedReceipt,MalformedReceipt}",
            "EvidenceIntegrity covers evidence bytes, custody and checker integrity",
            "InputValidity covers structural input validity",
            "ReceiptValidity covers schema, causal-event and cross-axis consistency",
            "Combination validity is fail-closed",
            "MalformedInput forces DomainApplicability=Indeterminate before any other combination check",
            "ClaimAdjudication in {Supported,Failed,Refuted} requires ExecutionDisposition=Completed",
            "InputValidity=WellFormedInput",
            "ClaimAdjudication=Pending requires a non-Completed disposition or non-complete evidence",
            "MalformedInput",
            "permits only ClaimAdjudication in {Unknown,Pending}",
            "every forbidden combination is ReceiptValidity=MalformedReceipt",
            "I14TerminalStatusV1 and i14_evaluate_terminal_status_v1 implement TerminalStatusTruthTableV1",
            "exhaustively enumerates all 3600 Cartesian tuples",
            "retains the raw producer tuple",
            "records every normalization action",
            "raw non-Indeterminate domain is a cross-axis contradiction",
            "i14_terminal_status_table_digest_v1 against I14_TERMINAL_STATUS_TABLE_V1_DIGEST_HEX",
            "ClaimAdjudication=Unknown -> 30",
            "Pending has no standalone primary exit",
            "Every declared exit code has a well-formed-receipt witness",
            "every precedence disjunct has an explicit witness including the necessarily malformed-receipt branch",
            "never scientific refutation",
            "Logical event identity, causal order and semantic payload are deterministic",
            "I14CanonicalTerminalResultInputV2 plus i14_canonical_terminal_result_v2 and i14_canonical_terminal_result_digest_v2",
            "revalidate the complete genesis-to-first-terminal prefix, validated cancellation card and terminal lifecycle",
            "recompute the clock-free logical execution root from canonical boundary/resource/work/request semantics",
            "request-inclusive terminal-prefix digest",
            "raw and normalized terminal axes, normalization actions, exit projection, semantic payload digest",
            "immutable strict-pre-boundary request cut",
            "A nonterminal frontier has no result constructor, and post-terminal boundaries are refused",
            "canonical terminal-result digest binds the canonicalized logical event/cause trace",
            "selected disposition/cause, request ids, logical sequences, boundary ordinals and semantic payload",
            "excluding calibrated monotonic timestamps, deadlines, live watchdog arrival and clock calibration",
            "i14_telemetry_envelope_digest_v2 revalidates the authoritative result through the same local implementation",
            "every raw request/observation timestamp and deadline",
            "I14_MAX_WATCHDOG_OBSERVATIONS_V1=4096",
            "I14_CANCELLATION_CARD_V2_KAT_HEX, I14_TERMINAL_PREFIX_V2_KAT_HEX, I14_CANONICAL_TERMINAL_RESULT_V2_KAT_HEX and I14_TELEMETRY_ENVELOPE_V2_KAT_HEX pin the four V2 byte layouts",
            "I14_WATCHDOG_RAW_TRACE_V2_KAT_HEX pins the complete multi-kind raw-watchdog byte layout",
            "a checked-in KAT alone is not evidence that an independent encoder reproduced it",
            "they do not by themselves prove an independent encoder exists or agrees",
            "I14_CANONICAL_TERMINAL_RESULT_V1_KAT_HEX/I14_TELEMETRY_ENVELOPE_V1_KAT_HEX remain only for legacy ledger readability",
            "identical recorded logical event/cause trace",
            "identical bound verification-receipt identities",
            "different live watchdog/expiry trace may honestly select a different disposition",
            "bounded to 64 KiB after schema-aware redaction",
            "governed digest/slot plus typed explicit-truncation metadata",
            "Licensed text, secrets, specimen identities, governed holdout/validation/population bytes before or after controlled reveal, and derived sensitive slices never enter a public artifact",
            "Sanitized manifest, adjudication receipt, logs, oracle output and replay capsule map exactly to I14_EVIDENCE_DURABLE",
            "sanitized minimized counterexamples, refutations and FailureBundles map exactly to I14_FAILURE_PERMANENT",
            "raw licensed/secret/specimen/governed-holdout bytes, derived sensitive slices and unredacted diagnostic slices map exactly to encrypted capability-controlled I14_GOVERNED_RESTRICTED",
            "I14ArtifactCategoryV1 and i14_retention_rule_v1 exhaustively map every named category",
            "complete-access-ledger requirement",
            "Schema-aware sanitization applies before retention to events, manifests, adjudication receipts, logs, oracle output, replay capsules, minimized counterexamples, refutations and FailureBundles",
            "retention tail must preserve the complete access ledger and class-specific retention/erasure decision",
        ] {
            assert!(
                row.g4_schedule.contains(token),
                "{} G4 omits {token}",
                row.leaf
            );
        }
        assert!(
            row.g4_schedule.contains(EXACT_EXIT_PRECEDENCE),
            "{} drifts the total ordered terminal-exit mapping",
            row.leaf
        );
        for token in [
            "threads {1,2,7}",
            "shards {1,4}",
            "Apple-aarch64",
            "x86_64",
            "bitwise comparison",
            "canonical terminal-result digest",
            "binds the canonicalized logical event/cause trace",
            "selected disposition/cause, request ids, logical sequences, boundary ordinals, lifecycle projection and semantic payload",
            "identical recorded logical event/cause trace",
            "identical topology/ISA/toolchain fingerprint",
            "excluding calibrated monotonic timestamps, deadlines, live watchdog arrival and clock calibration",
            "telemetry-envelope digest adds those noncanonical timing/calibration fields",
            "telemetry-envelope digest adds those noncanonical timing/calibration fields and is intentionally not bitwise compared",
        ] {
            assert!(
                row.g5_matrix.contains(token),
                "{} G5 omits {token}",
                row.leaf
            );
        }
        for event in [
            "execution.requested",
            "cancellation.requested",
            "cancellation.observed",
            "execution.cancelled",
            "execution.drained",
            "execution.finalized",
            "checkpoint.saved",
            "checkpoint.resumed",
            "checkpoint.forked",
            "evidence.adjudication_receipt",
            "evidence.failure_bundle",
        ] {
            assert_eq!(
                row.obs_events
                    .iter()
                    .filter(|candidate| **candidate == event)
                    .count(),
                1,
                "{} must emit {event} exactly once in its event schema",
                row.leaf,
            );
        }
    }
    assert_eq!(coverage.len(), draft.claims.len());
    for claim in &draft.claims {
        assert_eq!(coverage.get(claim.id), Some(&1), "{} coverage", claim.id);
    }

    assert_eq!(
        draft
            .waivers
            .iter()
            .map(|waiver| waiver.subject)
            .collect::<BTreeSet<_>>(),
        WAIVER_CONSUMERS
            .iter()
            .map(|(subject, _)| *subject)
            .collect()
    );
    for (subject, expected_consumers) in WAIVER_CONSUMERS {
        let waiver = draft
            .waivers
            .iter()
            .find(|waiver| waiver.subject == *subject)
            .unwrap_or_else(|| panic!("missing waiver {subject}"));
        assert!(
            waiver.predicate.contains(subject),
            "{subject} predicate omits its literal subject identity"
        );
        for token in [
            "fs-vvreg",
            "typed DischargeReceipt",
            "same-ID typed External discharge-envelope root",
            "FrozenManifest::amend",
            "AmendmentRecord",
            "removes this Waiver row",
            "IntegrityFailed",
            "atomic",
        ] {
            assert!(
                waiver.predicate.contains(token),
                "{subject} predicate omits {token}"
            );
        }
        assert!(
            waiver.promotion_effect.contains("only")
                && waiver.promotion_effect.contains("NoPromotionAuthority"),
            "{subject} scope"
        );
        let consumers = draft
            .obligations
            .iter()
            .filter(|row| row.decks.contains(subject))
            .map(|row| row.leaf)
            .collect::<BTreeSet<_>>();
        assert_eq!(
            consumers,
            expected_consumers.iter().copied().collect(),
            "{subject} blast radius"
        );
    }
    for row in draft
        .obligations
        .iter()
        .filter(|row| row.tier == CampaignTier::Core)
    {
        assert!(
            row.decks.iter().all(|deck| !waiver_ids.contains(deck)),
            "synthetic Core leaf {} consumes external waiver",
            row.leaf
        );
    }

    let frozen = draft.freeze().expect("the I14 seed must freeze");
    assert_eq!(frozen.initiative(), "I14");
    assert_eq!(frozen.version(), 1);
    assert_eq!(frozen.claims().len(), 28);
    assert_eq!(frozen.fixtures().len(), 53);
    assert_eq!(frozen.obligations().len(), 28);
    assert_eq!(frozen.waivers().len(), 7);
}

#[test]
#[allow(clippy::too_many_lines)]
fn i14_terminal_status_truth_table_is_exhaustive_typed_and_digest_pinned() {
    let expected_exit_codes = BTreeSet::from([0_u8, 10, 20, 21, 22, 30, 40, 50, 60, 70]);
    let mut observed_exit_codes = BTreeSet::new();
    let mut well_formed_receipt_exit_witnesses = BTreeSet::new();
    let mut well_formed_pending_exit_codes = BTreeSet::new();
    let mut tuples = 0usize;

    for execution in I14ExecutionDisposition::ALL {
        for claim in I14ClaimAdjudication::ALL {
            for completeness in I14EvidenceCompleteness::ALL {
                for integrity in I14EvidenceIntegrity::ALL {
                    for input in I14InputValidity::ALL {
                        for domain in I14DomainApplicability::ALL {
                            for support in I14OperationalSupport::ALL {
                                for receipt in I14ReceiptValidity::ALL {
                                    let raw = I14TerminalStatusV1 {
                                        execution,
                                        claim,
                                        completeness,
                                        integrity,
                                        input,
                                        domain,
                                        support,
                                        receipt,
                                    };
                                    let evaluated = i14_evaluate_terminal_status_v1(raw);
                                    let normalized = evaluated.normalized;
                                    assert_eq!(evaluated.raw, raw);
                                    let domain_was_forced = input
                                        == I14InputValidity::MalformedInput
                                        && domain != I14DomainApplicability::Indeterminate;
                                    assert_eq!(
                                        evaluated.normalization.domain_forced_indeterminate,
                                        domain_was_forced
                                    );
                                    assert_eq!(
                                        evaluated.normalization.receipt_marked_malformed,
                                        receipt == I14ReceiptValidity::WellFormedReceipt
                                            && normalized.receipt
                                                == I14ReceiptValidity::MalformedReceipt
                                    );
                                    tuples += 1;
                                    observed_exit_codes.insert(evaluated.exit_code);

                                    if normalized.receipt == I14ReceiptValidity::WellFormedReceipt {
                                        well_formed_receipt_exit_witnesses
                                            .insert(evaluated.exit_code);
                                    }
                                    if input == I14InputValidity::MalformedInput {
                                        assert_eq!(
                                            normalized.domain,
                                            I14DomainApplicability::Indeterminate,
                                            "malformed input must make applicability indeterminate: {raw:?}"
                                        );
                                        if domain_was_forced {
                                            assert_eq!(
                                                normalized.receipt,
                                                I14ReceiptValidity::MalformedReceipt,
                                                "raw malformed-input/domain contradiction was hidden: {raw:?}"
                                            );
                                        }
                                    } else {
                                        assert_eq!(normalized.domain, domain);
                                    }
                                    if matches!(
                                        normalized.claim,
                                        I14ClaimAdjudication::Supported
                                            | I14ClaimAdjudication::Failed
                                            | I14ClaimAdjudication::Refuted
                                    ) && normalized.receipt
                                        == I14ReceiptValidity::WellFormedReceipt
                                    {
                                        assert_eq!(
                                            normalized.execution,
                                            I14ExecutionDisposition::Completed
                                        );
                                        assert_eq!(
                                            normalized.completeness,
                                            I14EvidenceCompleteness::CompleteEvidence
                                        );
                                        assert_eq!(
                                            normalized.integrity,
                                            I14EvidenceIntegrity::IntegrityVerified
                                        );
                                        assert_eq!(
                                            normalized.input,
                                            I14InputValidity::WellFormedInput
                                        );
                                        assert_eq!(
                                            normalized.domain,
                                            I14DomainApplicability::Admitted
                                        );
                                        assert_eq!(
                                            normalized.support,
                                            I14OperationalSupport::SupportedOperation
                                        );
                                    }
                                    if normalized.claim == I14ClaimAdjudication::Pending
                                        && normalized.receipt
                                            == I14ReceiptValidity::WellFormedReceipt
                                    {
                                        well_formed_pending_exit_codes.insert(evaluated.exit_code);
                                        assert!(
                                            normalized.execution
                                                != I14ExecutionDisposition::Completed
                                                || normalized.completeness
                                                    != I14EvidenceCompleteness::CompleteEvidence,
                                            "completed complete-evidence Pending tuple was accepted: {raw:?}"
                                        );
                                    }
                                    if normalized.receipt == I14ReceiptValidity::WellFormedReceipt
                                        && (normalized.execution
                                            != I14ExecutionDisposition::Completed
                                            || normalized.completeness
                                                != I14EvidenceCompleteness::CompleteEvidence
                                            || normalized.integrity
                                                != I14EvidenceIntegrity::IntegrityVerified
                                            || normalized.input
                                                != I14InputValidity::WellFormedInput
                                            || normalized.domain
                                                != I14DomainApplicability::Admitted
                                            || normalized.support
                                                != I14OperationalSupport::SupportedOperation)
                                    {
                                        assert!(matches!(
                                            normalized.claim,
                                            I14ClaimAdjudication::Unknown
                                                | I14ClaimAdjudication::Pending
                                        ));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    assert_eq!(tuples, I14_TERMINAL_STATUS_TABLE_V1_TUPLES);
    assert_eq!(observed_exit_codes, expected_exit_codes);
    assert_eq!(well_formed_receipt_exit_witnesses, expected_exit_codes);
    assert_eq!(
        well_formed_pending_exit_codes,
        BTreeSet::from([20_u8, 21, 22, 30, 50, 60, 70])
    );
    assert_eq!(
        i14_terminal_status_table_digest_v1().to_hex(),
        I14_TERMINAL_STATUS_TABLE_V1_DIGEST_HEX
    );
}

#[test]
fn i14_drain_trigger_and_infrastructure_witness_wire_tables_are_exhaustively_pinned() {
    let trigger_vectors = [
        (I14DrainTriggerV2::NonCancellationDrain, vec![0x00]),
        (
            I14DrainTriggerV2::CancellationObserved {
                request_id: 0x0807_0605_0403_0201,
            },
            vec![0x01, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08],
        ),
        (
            I14DrainTriggerV2::ObservationTimeoutDrain {
                request_id: 0x1817_1615_1413_1211,
            },
            vec![0x02, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17, 0x18],
        ),
        (
            I14DrainTriggerV2::InfrastructureFailure {
                onset_logical_sequence: 0x2827_2625_2423_2221,
            },
            vec![0x03, 0x21, 0x22, 0x23, 0x24, 0x25, 0x26, 0x27, 0x28],
        ),
    ];
    let mut concatenated = Vec::new();
    for (trigger, expected) in trigger_vectors {
        let encoded = i14_drain_trigger_encoding_v2(trigger);
        assert_eq!(encoded, expected, "drain-trigger exact bytes drifted");
        concatenated.extend_from_slice(&encoded);
    }
    let to_hex = |bytes: &[u8]| {
        bytes
            .iter()
            .map(|byte| format!("{byte:02x}"))
            .collect::<String>()
    };
    assert_eq!(to_hex(&concatenated), I14_DRAIN_TRIGGER_ENCODING_V2_KAT_HEX);

    let sources = [
        (I14InfrastructureFailureSourceV2::WatchdogCoverage, 0),
        (I14InfrastructureFailureSourceV2::TilePollCoverage, 1),
        (
            I14InfrastructureFailureSourceV2::ExternalHeartbeatCoverage,
            2,
        ),
        (I14InfrastructureFailureSourceV2::DescendantDrain, 3),
        (
            I14InfrastructureFailureSourceV2::SpawnAfterFrontierClosure,
            4,
        ),
        (I14InfrastructureFailureSourceV2::Supervisor, 5),
        (I14InfrastructureFailureSourceV2::Authentication, 6),
        (I14InfrastructureFailureSourceV2::DrainProtocol, 7),
        (I14InfrastructureFailureSourceV2::PublicationProtocol, 8),
    ];
    for (source, expected_tag) in sources {
        assert_eq!(source as u8, expected_tag);
        let encoded = i14_infrastructure_failure_onset_encoding_v2(Some(i14_infrastructure_onset(
            source,
            0x0807_0605_0403_0201,
            0xdead_beef,
            0xa7,
        )));
        assert_eq!(encoded.len(), 42);
        assert_eq!(encoded[0], 1, "present witness tag drifted");
        assert_eq!(
            encoded[9], expected_tag,
            "infrastructure source wire tag drifted for {source:?}"
        );
    }
    assert_eq!(
        i14_infrastructure_failure_onset_encoding_v2(None),
        vec![0],
        "absent infrastructure onset must preserve its one-byte form"
    );
    let witness = i14_infrastructure_onset(
        I14InfrastructureFailureSourceV2::Supervisor,
        0x0807_0605_0403_0201,
        0xdead_beef,
        0xa7,
    );
    assert_eq!(
        to_hex(&i14_infrastructure_failure_onset_encoding_v2(Some(witness))),
        I14_INFRASTRUCTURE_FAILURE_ONSET_ENCODING_V2_KAT_HEX,
        "canonical witness bytes must exclude calibrated onset time"
    );
}

#[test]
#[allow(clippy::too_many_lines)]
fn i14_terminal_exit_precedence_has_explicit_well_formed_witnesses() {
    let supported = I14TerminalStatusV1 {
        execution: I14ExecutionDisposition::Completed,
        claim: I14ClaimAdjudication::Supported,
        completeness: I14EvidenceCompleteness::CompleteEvidence,
        integrity: I14EvidenceIntegrity::IntegrityVerified,
        input: I14InputValidity::WellFormedInput,
        domain: I14DomainApplicability::Admitted,
        support: I14OperationalSupport::SupportedOperation,
        receipt: I14ReceiptValidity::WellFormedReceipt,
    };
    let pending = I14TerminalStatusV1 {
        claim: I14ClaimAdjudication::Pending,
        ..supported
    };
    let unknown = I14TerminalStatusV1 {
        claim: I14ClaimAdjudication::Unknown,
        ..supported
    };
    let witnesses = [
        (supported, 0),
        (
            I14TerminalStatusV1 {
                claim: I14ClaimAdjudication::Failed,
                ..supported
            },
            10,
        ),
        (
            I14TerminalStatusV1 {
                claim: I14ClaimAdjudication::Refuted,
                ..supported
            },
            40,
        ),
        (unknown, 30),
        (
            I14TerminalStatusV1 {
                completeness: I14EvidenceCompleteness::PartialEvidence,
                ..pending
            },
            50,
        ),
        (
            I14TerminalStatusV1 {
                integrity: I14EvidenceIntegrity::IntegrityFailed,
                ..unknown
            },
            60,
        ),
        (
            I14TerminalStatusV1 {
                execution: I14ExecutionDisposition::InfrastructureFailed,
                ..unknown
            },
            70,
        ),
        (
            I14TerminalStatusV1 {
                execution: I14ExecutionDisposition::Cancelled,
                ..pending
            },
            20,
        ),
        (
            I14TerminalStatusV1 {
                execution: I14ExecutionDisposition::TimedOut,
                ..pending
            },
            21,
        ),
        (
            I14TerminalStatusV1 {
                execution: I14ExecutionDisposition::BudgetExhausted,
                ..pending
            },
            22,
        ),
    ];
    for (status, expected_exit) in witnesses {
        let evaluated = i14_evaluate_terminal_status_v1(status);
        assert_eq!(evaluated.exit_code, expected_exit, "{status:?}");
        assert_eq!(
            evaluated.normalized.receipt,
            I14ReceiptValidity::WellFormedReceipt,
            "witness must not rely on a malformed receipt: {status:?}"
        );
    }

    let precedence_witnesses = [
        (
            I14TerminalStatusV1 {
                receipt: I14ReceiptValidity::MalformedReceipt,
                execution: I14ExecutionDisposition::InfrastructureFailed,
                ..unknown
            },
            60,
        ),
        (
            I14TerminalStatusV1 {
                execution: I14ExecutionDisposition::InfrastructureFailed,
                completeness: I14EvidenceCompleteness::NoEvidence,
                ..unknown
            },
            70,
        ),
        (
            I14TerminalStatusV1 {
                domain: I14DomainApplicability::OutOfDomain,
                completeness: I14EvidenceCompleteness::NoEvidence,
                ..unknown
            },
            30,
        ),
        (
            I14TerminalStatusV1 {
                domain: I14DomainApplicability::OutOfDomain,
                completeness: I14EvidenceCompleteness::NoEvidence,
                ..pending
            },
            30,
        ),
        (
            I14TerminalStatusV1 {
                completeness: I14EvidenceCompleteness::NoEvidence,
                ..unknown
            },
            50,
        ),
    ];
    for (status, expected_exit) in precedence_witnesses {
        assert_eq!(
            i14_evaluate_terminal_status_v1(status).exit_code,
            expected_exit,
            "precedence drift for {status:?}"
        );
    }
}

#[test]
#[allow(clippy::too_many_lines)]
fn i14_terminal_cause_selector_is_scoped_bounded_causal_and_total() {
    let scope = [1_u64, 2, 3];
    let observer_tiles = [41_u64, 42, 43];
    let completed = I14TerminalBoundaryV1 {
        boundary_ordinal: 10,
        logical_sequence: 100,
        monotonic_ns: 1_000,
        scope_ancestry: &scope,
        admitted_observer_tile_ids: &observer_tiles,
        infrastructure_failed: false,
        timed_out: false,
        budget_exhausted: false,
        completed: true,
    };
    let observed_request = I14CancellationRequestV1 {
        request_id: 7,
        scope_root: 2,
        logical_sequence: 10,
        requested_monotonic_ns: 100,
        observation_deadline_ns: 500,
        observation: Some(I14CancellationObservationV1 {
            logical_sequence: 20,
            monotonic_ns: 200,
            observing_tile_id: 41,
            latest_completed_boundary_ordinal: Some(5),
        }),
    };
    let pending_request = I14CancellationRequestV1 {
        request_id: 8,
        logical_sequence: 30,
        observation_deadline_ns: 1_100,
        observation: None,
        ..observed_request
    };

    assert_eq!(
        i14_select_terminal_boundary_v1(completed, &[]),
        Ok(I14TerminalBoundaryDecisionV1::Selected {
            disposition: I14ExecutionDisposition::Completed,
            request_id: None,
        })
    );
    let not_terminal = I14TerminalBoundaryV1 {
        completed: false,
        ..completed
    };
    assert_eq!(
        i14_select_terminal_boundary_v1(not_terminal, &[]),
        Ok(I14TerminalBoundaryDecisionV1::NotTerminal)
    );
    assert_eq!(
        i14_select_terminal_boundary_v1(not_terminal, &[pending_request]),
        Ok(I14TerminalBoundaryDecisionV1::NotTerminal),
        "a pending cancellation request defers only a normal terminal candidate"
    );
    assert_eq!(
        i14_select_terminal_boundary_v1(
            I14TerminalBoundaryV1 {
                timed_out: true,
                ..not_terminal
            },
            &[],
        ),
        Ok(I14TerminalBoundaryDecisionV1::Selected {
            disposition: I14ExecutionDisposition::TimedOut,
            request_id: None,
        })
    );
    assert_eq!(
        i14_select_terminal_boundary_v1(
            I14TerminalBoundaryV1 {
                budget_exhausted: true,
                ..completed
            },
            &[],
        ),
        Ok(I14TerminalBoundaryDecisionV1::Selected {
            disposition: I14ExecutionDisposition::BudgetExhausted,
            request_id: None,
        })
    );
    assert_eq!(
        i14_select_terminal_boundary_v1(completed, &[observed_request]),
        Ok(I14TerminalBoundaryDecisionV1::Selected {
            disposition: I14ExecutionDisposition::Cancelled,
            request_id: Some(7),
        })
    );
    assert_eq!(
        i14_select_terminal_boundary_v1(
            I14TerminalBoundaryV1 {
                boundary_ordinal: 0,
                ..completed
            },
            &[I14CancellationRequestV1 {
                observation: Some(I14CancellationObservationV1 {
                    latest_completed_boundary_ordinal: None,
                    ..observed_request.observation.expect("observation")
                }),
                ..observed_request
            }],
        ),
        Ok(I14TerminalBoundaryDecisionV1::Selected {
            disposition: I14ExecutionDisposition::Cancelled,
            request_id: Some(7),
        }),
        "None represents observation before the first completed boundary"
    );
    assert_eq!(
        i14_select_terminal_boundary_v1(completed, &[pending_request]),
        Ok(I14TerminalBoundaryDecisionV1::DeferredByCancellation { request_id: 8 })
    );
    assert_eq!(
        i14_select_terminal_boundary_v1(
            completed,
            &[I14CancellationRequestV1 {
                scope_root: 99,
                ..observed_request
            }],
        ),
        Ok(I14TerminalBoundaryDecisionV1::Selected {
            disposition: I14ExecutionDisposition::Completed,
            request_id: None,
        }),
        "an unrelated scope must not cancel this boundary"
    );
    assert_eq!(
        i14_select_terminal_boundary_v1(
            completed,
            &[I14CancellationRequestV1 {
                scope_root: 99,
                observation: Some(I14CancellationObservationV1 {
                    observing_tile_id: 99,
                    ..observed_request.observation.expect("observation")
                }),
                ..observed_request
            }],
        ),
        Ok(I14TerminalBoundaryDecisionV1::Selected {
            disposition: I14ExecutionDisposition::Completed,
            request_id: None,
        }),
        "an out-of-scope observation uses its own scope's observer catalog"
    );
    assert_eq!(
        i14_select_terminal_boundary_v1(
            completed,
            &[I14CancellationRequestV1 {
                logical_sequence: 110,
                requested_monotonic_ns: 1_100,
                observation_deadline_ns: 1_500,
                observation: Some(I14CancellationObservationV1 {
                    logical_sequence: 120,
                    monotonic_ns: 1_200,
                    observing_tile_id: 42,
                    latest_completed_boundary_ordinal: Some(11),
                }),
                ..observed_request
            }],
        ),
        Ok(I14TerminalBoundaryDecisionV1::Selected {
            disposition: I14ExecutionDisposition::Completed,
            request_id: None,
        }),
        "a logically later request must not rewrite this boundary"
    );

    let missed_request = I14CancellationRequestV1 {
        request_id: 6,
        logical_sequence: 5,
        observation_deadline_ns: 150,
        observation: None,
        ..observed_request
    };
    assert_eq!(
        i14_select_terminal_boundary_v1(completed, &[observed_request, missed_request]),
        Ok(I14TerminalBoundaryDecisionV1::Selected {
            disposition: I14ExecutionDisposition::TimedOut,
            request_id: Some(6),
        }),
        "one missed request must outrank another observed cancellation"
    );
    let later_observed = I14CancellationRequestV1 {
        request_id: 9,
        logical_sequence: 40,
        requested_monotonic_ns: 250,
        observation: Some(I14CancellationObservationV1 {
            logical_sequence: 50,
            monotonic_ns: 300,
            observing_tile_id: 43,
            latest_completed_boundary_ordinal: Some(6),
        }),
        ..observed_request
    };
    assert_eq!(
        i14_select_terminal_boundary_v1(completed, &[later_observed, observed_request]),
        Ok(I14TerminalBoundaryDecisionV1::Selected {
            disposition: I14ExecutionDisposition::Cancelled,
            request_id: Some(7),
        }),
        "multiple observed requests use deterministic logical-sequence order"
    );
    assert_eq!(
        i14_select_terminal_boundary_v1(
            I14TerminalBoundaryV1 {
                infrastructure_failed: true,
                timed_out: true,
                budget_exhausted: true,
                ..completed
            },
            &[observed_request, missed_request],
        ),
        Ok(I14TerminalBoundaryDecisionV1::Selected {
            disposition: I14ExecutionDisposition::InfrastructureFailed,
            request_id: None,
        })
    );

    assert_eq!(
        i14_select_terminal_boundary_v1(
            I14TerminalBoundaryV1 {
                scope_ancestry: &[],
                ..completed
            },
            &[],
        ),
        Err(I14TerminalCauseRefusalV1::EmptyScopeAncestry)
    );
    let deep_scope = vec![0_u64; I14_MAX_SCOPE_ANCESTRY_V1 + 1];
    assert_eq!(
        i14_select_terminal_boundary_v1(
            I14TerminalBoundaryV1 {
                scope_ancestry: &deep_scope,
                ..completed
            },
            &[],
        ),
        Err(I14TerminalCauseRefusalV1::ScopeAncestryTooDeep {
            depth: I14_MAX_SCOPE_ANCESTRY_V1 + 1,
            cap: I14_MAX_SCOPE_ANCESTRY_V1,
        })
    );
    assert_eq!(
        i14_select_terminal_boundary_v1(
            I14TerminalBoundaryV1 {
                scope_ancestry: &[1, 2, 1],
                ..completed
            },
            &[],
        ),
        Err(I14TerminalCauseRefusalV1::DuplicateScopeId { scope_id: 1 })
    );
    let too_many_observer_tiles = vec![0_u64; I14_MAX_OBSERVER_TILES_V1 + 1];
    assert_eq!(
        i14_select_terminal_boundary_v1(
            I14TerminalBoundaryV1 {
                admitted_observer_tile_ids: &too_many_observer_tiles,
                ..completed
            },
            &[],
        ),
        Err(I14TerminalCauseRefusalV1::TooManyObserverTiles {
            count: I14_MAX_OBSERVER_TILES_V1 + 1,
            cap: I14_MAX_OBSERVER_TILES_V1,
        })
    );
    assert_eq!(
        i14_select_terminal_boundary_v1(
            I14TerminalBoundaryV1 {
                admitted_observer_tile_ids: &[41, 41],
                ..completed
            },
            &[],
        ),
        Err(I14TerminalCauseRefusalV1::DuplicateObserverTileId {
            observing_tile_id: 41,
        })
    );
    let max_scope = (0..I14_MAX_SCOPE_ANCESTRY_V1)
        .map(|index| u64::try_from(index).expect("scope index fits u64"))
        .collect::<Vec<_>>();
    let max_observer_tiles = (0..I14_MAX_OBSERVER_TILES_V1)
        .map(|index| u64::try_from(index).expect("observer index fits u64"))
        .collect::<Vec<_>>();
    let max_requests = (0..I14_MAX_CANCELLATION_REQUESTS_V1)
        .map(|index| {
            let logical_index = u64::try_from(index).expect("request index fits u64");
            I14CancellationRequestV1 {
                request_id: logical_index,
                scope_root: logical_index
                    % u64::try_from(I14_MAX_SCOPE_ANCESTRY_V1).expect("scope cap fits u64"),
                logical_sequence: logical_index * 2 + 1,
                requested_monotonic_ns: logical_index,
                observation_deadline_ns: 200_000,
                observation: None,
            }
        })
        .collect::<Vec<_>>();
    assert_eq!(
        i14_select_terminal_boundary_v1(
            I14TerminalBoundaryV1 {
                logical_sequence: 100_000,
                monotonic_ns: 100_000,
                scope_ancestry: &max_scope,
                admitted_observer_tile_ids: &max_observer_tiles,
                ..completed
            },
            &max_requests,
        ),
        Ok(I14TerminalBoundaryDecisionV1::DeferredByCancellation { request_id: 0 }),
        "all three public selector caps are inclusive"
    );
    let duplicate_observers_forward = i14_select_terminal_boundary_v1(
        I14TerminalBoundaryV1 {
            admitted_observer_tile_ids: &[43, 43, 41, 41, 42],
            ..completed
        },
        &[],
    );
    let duplicate_observers_reverse = i14_select_terminal_boundary_v1(
        I14TerminalBoundaryV1 {
            admitted_observer_tile_ids: &[42, 41, 41, 43, 43],
            ..completed
        },
        &[],
    );
    assert_eq!(
        duplicate_observers_forward, duplicate_observers_reverse,
        "malformed observer-catalog refusals must be order invariant"
    );
    assert_eq!(
        duplicate_observers_forward,
        Err(I14TerminalCauseRefusalV1::DuplicateObserverTileId {
            observing_tile_id: 41,
        })
    );
    let too_many = vec![pending_request; I14_MAX_CANCELLATION_REQUESTS_V1 + 1];
    assert_eq!(
        i14_select_terminal_boundary_v1(completed, &too_many),
        Err(I14TerminalCauseRefusalV1::TooManyRequests {
            count: I14_MAX_CANCELLATION_REQUESTS_V1 + 1,
            cap: I14_MAX_CANCELLATION_REQUESTS_V1,
        })
    );
    assert_eq!(
        i14_select_terminal_boundary_v1(
            completed,
            &[
                observed_request,
                I14CancellationRequestV1 {
                    logical_sequence: 30,
                    observation: None,
                    ..observed_request
                },
            ],
        ),
        Err(I14TerminalCauseRefusalV1::DuplicateRequestId { request_id: 7 })
    );
    assert_eq!(
        i14_select_terminal_boundary_v1(
            completed,
            &[
                observed_request,
                I14CancellationRequestV1 {
                    request_id: 9,
                    observation: None,
                    ..observed_request
                },
            ],
        ),
        Err(I14TerminalCauseRefusalV1::DuplicateLogicalSequence {
            logical_sequence: 10,
        })
    );
    let malformed_a = I14CancellationRequestV1 {
        request_id: 1,
        logical_sequence: 10,
        observation: None,
        ..observed_request
    };
    let malformed_b = I14CancellationRequestV1 {
        request_id: 2,
        logical_sequence: 10,
        observation: None,
        ..observed_request
    };
    let malformed_c = I14CancellationRequestV1 {
        request_id: 1,
        logical_sequence: 30,
        observation: None,
        ..observed_request
    };
    let forward =
        i14_select_terminal_boundary_v1(completed, &[malformed_a, malformed_b, malformed_c]);
    let reverse =
        i14_select_terminal_boundary_v1(completed, &[malformed_c, malformed_b, malformed_a]);
    assert_eq!(
        forward, reverse,
        "malformed refusals must be order invariant"
    );
    assert_eq!(
        forward,
        Err(I14TerminalCauseRefusalV1::DuplicateLogicalSequence {
            logical_sequence: 10,
        })
    );
    let tied_valid = I14CancellationRequestV1 {
        observation: None,
        ..observed_request
    };
    let tied_malformed = I14CancellationRequestV1 {
        observation_deadline_ns: 99,
        ..tied_valid
    };
    let tied_forward = i14_select_terminal_boundary_v1(completed, &[tied_valid, tied_malformed]);
    let tied_reverse = i14_select_terminal_boundary_v1(completed, &[tied_malformed, tied_valid]);
    assert_eq!(
        tied_forward, tied_reverse,
        "equal request-id/sequence payload ties must have an order-invariant refusal"
    );
    assert_eq!(
        tied_forward,
        Err(I14TerminalCauseRefusalV1::DeadlineBeforeRequest { request_id: 7 })
    );
    assert_eq!(
        i14_select_terminal_boundary_v1(
            completed,
            &[I14CancellationRequestV1 {
                observation_deadline_ns: 99,
                observation: None,
                ..observed_request
            }],
        ),
        Err(I14TerminalCauseRefusalV1::DeadlineBeforeRequest { request_id: 7 })
    );
    assert_eq!(
        i14_select_terminal_boundary_v1(
            completed,
            &[I14CancellationRequestV1 {
                observation: Some(I14CancellationObservationV1 {
                    logical_sequence: 10,
                    monotonic_ns: 100,
                    observing_tile_id: 41,
                    latest_completed_boundary_ordinal: Some(5),
                }),
                ..observed_request
            }],
        ),
        Err(I14TerminalCauseRefusalV1::ObservationBeforeRequest { request_id: 7 })
    );
    assert_eq!(
        i14_select_terminal_boundary_v1(
            completed,
            &[I14CancellationRequestV1 {
                observation_deadline_ns: 2_000,
                observation: Some(I14CancellationObservationV1 {
                    logical_sequence: 20,
                    monotonic_ns: 1_001,
                    observing_tile_id: 41,
                    latest_completed_boundary_ordinal: Some(5),
                }),
                ..observed_request
            }],
        ),
        Err(I14TerminalCauseRefusalV1::NonMonotonicLogicalTimestamp {
            earlier_logical_sequence: 20,
            earlier_monotonic_ns: 1_001,
            later_logical_sequence: 100,
            later_monotonic_ns: 1_000,
        })
    );
    assert_eq!(
        i14_select_terminal_boundary_v1(
            completed,
            &[I14CancellationRequestV1 {
                observation: Some(I14CancellationObservationV1 {
                    latest_completed_boundary_ordinal: Some(10),
                    ..observed_request.observation.expect("observation")
                }),
                ..observed_request
            }],
        ),
        Err(I14TerminalCauseRefusalV1::ObservationBoundaryNotBeforeCandidate { request_id: 7 })
    );
    assert_eq!(
        i14_select_terminal_boundary_v1(
            completed,
            &[I14CancellationRequestV1 {
                observation: Some(I14CancellationObservationV1 {
                    observing_tile_id: 99,
                    ..observed_request.observation.expect("observation")
                }),
                ..observed_request
            }],
        ),
        Err(I14TerminalCauseRefusalV1::UnknownObserverTile {
            request_id: 7,
            observing_tile_id: 99,
        })
    );
    assert_eq!(
        i14_select_terminal_boundary_v1(
            completed,
            &[I14CancellationRequestV1 {
                requested_monotonic_ns: 1_001,
                observation_deadline_ns: 1_500,
                observation: None,
                ..observed_request
            }],
        ),
        Err(I14TerminalCauseRefusalV1::NonMonotonicLogicalTimestamp {
            earlier_logical_sequence: 10,
            earlier_monotonic_ns: 1_001,
            later_logical_sequence: 100,
            later_monotonic_ns: 1_000,
        })
    );
    assert_eq!(
        i14_select_terminal_boundary_v1(
            completed,
            &[I14CancellationRequestV1 {
                logical_sequence: 110,
                requested_monotonic_ns: 999,
                observation_deadline_ns: 1_500,
                observation: None,
                ..observed_request
            }],
        ),
        Err(I14TerminalCauseRefusalV1::NonMonotonicLogicalTimestamp {
            earlier_logical_sequence: 100,
            earlier_monotonic_ns: 1_000,
            later_logical_sequence: 110,
            later_monotonic_ns: 999,
        })
    );
    assert_eq!(
        i14_select_terminal_boundary_v1(
            completed,
            &[I14CancellationRequestV1 {
                observation_deadline_ns: 2_000,
                observation: Some(I14CancellationObservationV1 {
                    logical_sequence: 120,
                    monotonic_ns: 999,
                    observing_tile_id: 41,
                    latest_completed_boundary_ordinal: Some(10),
                }),
                ..observed_request
            }],
        ),
        Err(I14TerminalCauseRefusalV1::NonMonotonicLogicalTimestamp {
            earlier_logical_sequence: 100,
            earlier_monotonic_ns: 1_000,
            later_logical_sequence: 120,
            later_monotonic_ns: 999,
        })
    );
    let timestamp_inversion_a = I14CancellationRequestV1 {
        request_id: 1,
        logical_sequence: 10,
        requested_monotonic_ns: 100,
        observation_deadline_ns: 1_500,
        observation: None,
        ..observed_request
    };
    let timestamp_inversion_b = I14CancellationRequestV1 {
        request_id: 2,
        logical_sequence: 30,
        requested_monotonic_ns: 90,
        ..timestamp_inversion_a
    };
    let timestamp_inversion_forward =
        i14_select_terminal_boundary_v1(completed, &[timestamp_inversion_a, timestamp_inversion_b]);
    let timestamp_inversion_reverse =
        i14_select_terminal_boundary_v1(completed, &[timestamp_inversion_b, timestamp_inversion_a]);
    assert_eq!(
        timestamp_inversion_forward, timestamp_inversion_reverse,
        "cross-request timestamp refusals must be presentation-order invariant"
    );
    assert_eq!(
        timestamp_inversion_forward,
        Err(I14TerminalCauseRefusalV1::NonMonotonicLogicalTimestamp {
            earlier_logical_sequence: 10,
            earlier_monotonic_ns: 100,
            later_logical_sequence: 30,
            later_monotonic_ns: 90,
        })
    );
}

#[test]
#[allow(clippy::too_many_lines)]
fn i14_canonical_result_and_telemetry_separate_semantics_from_raw_timing() {
    let scope = [1_u64, 2];
    let observers = [42_u64, 41];
    let sorted_observers = [41_u64, 42];
    let boundary = I14TerminalBoundaryV1 {
        boundary_ordinal: 10,
        logical_sequence: 100,
        monotonic_ns: 1_000,
        scope_ancestry: &scope,
        admitted_observer_tile_ids: &observers,
        infrastructure_failed: false,
        timed_out: false,
        budget_exhausted: false,
        completed: true,
    };
    let observed = I14CancellationRequestV1 {
        request_id: 7,
        scope_root: 1,
        logical_sequence: 10,
        requested_monotonic_ns: 100,
        observation_deadline_ns: 500,
        observation: Some(I14CancellationObservationV1 {
            logical_sequence: 20,
            monotonic_ns: 200,
            observing_tile_id: 41,
            latest_completed_boundary_ordinal: Some(5),
        }),
    };
    let out_of_scope = I14CancellationRequestV1 {
        request_id: 8,
        scope_root: 99,
        logical_sequence: 30,
        requested_monotonic_ns: 300,
        observation_deadline_ns: 900,
        observation: None,
    };
    let late = I14CancellationRequestV1 {
        request_id: 9,
        scope_root: 1,
        logical_sequence: 110,
        requested_monotonic_ns: 1_100,
        observation_deadline_ns: 1_500,
        observation: None,
    };
    let status = I14TerminalStatusV1 {
        execution: I14ExecutionDisposition::Cancelled,
        claim: I14ClaimAdjudication::Pending,
        completeness: I14EvidenceCompleteness::PartialEvidence,
        integrity: I14EvidenceIntegrity::IntegrityVerified,
        input: I14InputValidity::WellFormedInput,
        domain: I14DomainApplicability::Admitted,
        support: I14OperationalSupport::SupportedOperation,
        receipt: I14ReceiptValidity::WellFormedReceipt,
    };
    let semantic_payload_digest = ContentHash([0x5a; 32]);
    let requests = [late, out_of_scope, observed];
    let canonical_input = I14CanonicalTerminalResultInputV1 {
        boundary,
        requests: &requests,
        terminal_status: status,
        semantic_payload_digest,
    };
    let canonical = i14_canonical_terminal_result_v1(canonical_input)
        .expect("selected terminal result must canonicalize");
    assert_eq!(canonical.boundary_ordinal(), 10);
    assert_eq!(canonical.boundary_logical_sequence(), 100);
    assert_eq!(canonical.scope_ancestry(), &scope);
    assert_eq!(canonical.admitted_observer_tile_ids(), &sorted_observers);
    assert_eq!(
        canonical.decision(),
        I14TerminalBoundaryDecisionV1::Selected {
            disposition: I14ExecutionDisposition::Cancelled,
            request_id: Some(7),
        }
    );
    assert_eq!(
        canonical.requests().len(),
        2,
        "late requests stay outside the frozen cut"
    );
    assert_eq!(
        canonical.requests()[0].state,
        I14CancellationRequestStateV1::ObservedWithinDeadline
    );
    assert_eq!(
        canonical.requests()[1].state,
        I14CancellationRequestStateV1::OutOfScope
    );
    assert_eq!(canonical.terminal_evaluation().exit_code, 20);
    assert_eq!(canonical.semantic_payload_digest(), semantic_payload_digest);
    assert_eq!(
        canonical.digest(),
        i14_canonical_terminal_result_digest_v1(canonical_input).expect("digest convenience path")
    );
    assert_eq!(
        canonical.digest().to_hex(),
        I14_CANONICAL_TERMINAL_RESULT_V1_KAT_HEX,
        "canonical terminal-result v1 byte layout drifted"
    );

    let shifted_boundary = I14TerminalBoundaryV1 {
        monotonic_ns: 2_000,
        admitted_observer_tile_ids: &sorted_observers,
        ..boundary
    };
    let shifted_observed = I14CancellationRequestV1 {
        requested_monotonic_ns: 150,
        observation_deadline_ns: 600,
        observation: Some(I14CancellationObservationV1 {
            monotonic_ns: 250,
            ..observed.observation.expect("observation")
        }),
        ..observed
    };
    let shifted_out_of_scope = I14CancellationRequestV1 {
        requested_monotonic_ns: 350,
        observation_deadline_ns: 1_900,
        ..out_of_scope
    };
    let shifted_late = I14CancellationRequestV1 {
        requested_monotonic_ns: 2_100,
        observation_deadline_ns: 2_500,
        ..late
    };
    let shifted_requests = [shifted_observed, shifted_late, shifted_out_of_scope];
    let shifted_canonical = i14_canonical_terminal_result_v1(I14CanonicalTerminalResultInputV1 {
        boundary: shifted_boundary,
        requests: &shifted_requests,
        terminal_status: status,
        semantic_payload_digest,
    })
    .expect("timing-equivalent result");
    assert_eq!(
        canonical, shifted_canonical,
        "raw timing and presentation shifts preserving semantic relations must not change canonical identity"
    );

    let late_malformed = I14CancellationRequestV1 {
        request_id: 7,
        logical_sequence: 120,
        requested_monotonic_ns: 1_200,
        observation_deadline_ns: 1_100,
        observation: None,
        ..late
    };
    let requests_with_late_malformed = [observed, out_of_scope, late_malformed];
    let canonical_with_late_malformed =
        i14_canonical_terminal_result_v1(I14CanonicalTerminalResultInputV1 {
            boundary,
            requests: &requests_with_late_malformed,
            terminal_status: status,
            semantic_payload_digest,
        })
        .expect("post-boundary corruption cannot rewrite an already finalized canonical cut");
    assert_eq!(canonical, canonical_with_late_malformed);

    for boundary_sequence_collision in [
        I14CancellationRequestV1 {
            logical_sequence: boundary.logical_sequence,
            requested_monotonic_ns: boundary.monotonic_ns,
            ..late
        },
        I14CancellationRequestV1 {
            observation: Some(I14CancellationObservationV1 {
                logical_sequence: boundary.logical_sequence,
                monotonic_ns: boundary.monotonic_ns,
                latest_completed_boundary_ordinal: Some(9),
                ..observed.observation.expect("observation")
            }),
            ..observed
        },
    ] {
        assert_eq!(
            i14_canonical_terminal_result_v1(I14CanonicalTerminalResultInputV1 {
                boundary,
                requests: &[boundary_sequence_collision],
                terminal_status: status,
                semantic_payload_digest,
            }),
            Err(I14CanonicalResultRefusalV1::TerminalCause(
                I14TerminalCauseRefusalV1::DuplicateLogicalSequence {
                    logical_sequence: boundary.logical_sequence,
                }
            )),
            "an event equal to the boundary sequence is a collision, not a late record"
        );
    }

    let late_out_of_scope_observation = I14CancellationRequestV1 {
        observation: Some(I14CancellationObservationV1 {
            logical_sequence: 120,
            monotonic_ns: 1_200,
            observing_tile_id: 42,
            latest_completed_boundary_ordinal: Some(10),
        }),
        ..out_of_scope
    };
    let requests_with_late_observation = [observed, late_out_of_scope_observation, late];
    assert_eq!(
        canonical,
        i14_canonical_terminal_result_v1(I14CanonicalTerminalResultInputV1 {
            boundary,
            requests: &requests_with_late_observation,
            terminal_status: status,
            semantic_payload_digest,
        })
        .expect("a post-boundary observation is outside the canonical cut")
    );

    let changed_payload =
        i14_canonical_terminal_result_digest_v1(I14CanonicalTerminalResultInputV1 {
            semantic_payload_digest: ContentHash([0x5b; 32]),
            ..canonical_input
        })
        .expect("changed semantic payload remains valid");
    assert_ne!(canonical.digest(), changed_payload);
    let changed_cause =
        i14_canonical_terminal_result_digest_v1(I14CanonicalTerminalResultInputV1 {
            boundary: I14TerminalBoundaryV1 {
                timed_out: true,
                ..boundary
            },
            terminal_status: I14TerminalStatusV1 {
                execution: I14ExecutionDisposition::TimedOut,
                ..status
            },
            ..canonical_input
        })
        .expect("changed terminal cause remains valid");
    assert_ne!(canonical.digest(), changed_cause);

    let pending = I14CancellationRequestV1 {
        request_id: 11,
        scope_root: 1,
        logical_sequence: 40,
        requested_monotonic_ns: 400,
        observation_deadline_ns: 1_500,
        observation: None,
    };
    assert_eq!(
        i14_canonical_terminal_result_v1(I14CanonicalTerminalResultInputV1 {
            boundary,
            requests: &[pending],
            terminal_status: status,
            semantic_payload_digest,
        }),
        Err(I14CanonicalResultRefusalV1::BoundaryNotTerminal {
            decision: I14TerminalBoundaryDecisionV1::DeferredByCancellation { request_id: 11 },
        })
    );
    assert_eq!(
        i14_canonical_terminal_result_v1(I14CanonicalTerminalResultInputV1 {
            boundary: I14TerminalBoundaryV1 {
                completed: false,
                ..boundary
            },
            requests: &[],
            terminal_status: status,
            semantic_payload_digest,
        }),
        Err(I14CanonicalResultRefusalV1::BoundaryNotTerminal {
            decision: I14TerminalBoundaryDecisionV1::NotTerminal,
        })
    );
    assert_eq!(
        i14_canonical_terminal_result_v1(I14CanonicalTerminalResultInputV1 {
            terminal_status: I14TerminalStatusV1 {
                execution: I14ExecutionDisposition::Completed,
                claim: I14ClaimAdjudication::Supported,
                completeness: I14EvidenceCompleteness::CompleteEvidence,
                ..status
            },
            ..canonical_input
        }),
        Err(I14CanonicalResultRefusalV1::ExecutionDispositionMismatch {
            selected: I14ExecutionDisposition::Cancelled,
            receipt: I14ExecutionDisposition::Completed,
        })
    );

    let watchdogs = [
        I14WatchdogObservationV1 {
            observation_id: 2,
            kind: I14WatchdogObservationKindV1::ExternalHeartbeat,
            monotonic_ns: 950,
        },
        I14WatchdogObservationV1 {
            observation_id: 1,
            kind: I14WatchdogObservationKindV1::Poll,
            monotonic_ns: 900,
        },
    ];
    let reversed_watchdogs = [watchdogs[1], watchdogs[0]];
    let reordered_requests = [observed, out_of_scope, late];
    let telemetry_input = I14TelemetryEnvelopeInputV1 {
        boundary,
        requests: &requests,
        terminal_status: status,
        semantic_payload_digest,
        watchdog_observations: &watchdogs,
        clock_calibration_artifact: ContentHash([0x33; 32]),
    };
    let telemetry =
        i14_telemetry_envelope_digest_v1(telemetry_input).expect("valid telemetry envelope");
    assert_eq!(
        telemetry.to_hex(),
        I14_TELEMETRY_ENVELOPE_V1_KAT_HEX,
        "telemetry-envelope v1 byte layout drifted"
    );
    assert_eq!(
        telemetry,
        i14_telemetry_envelope_digest_v1(I14TelemetryEnvelopeInputV1 {
            boundary: I14TerminalBoundaryV1 {
                admitted_observer_tile_ids: &sorted_observers,
                ..boundary
            },
            requests: &reordered_requests,
            watchdog_observations: &reversed_watchdogs,
            ..telemetry_input
        })
        .expect("presentation-order-equivalent telemetry"),
    );
    assert_ne!(
        telemetry,
        i14_telemetry_envelope_digest_v1(I14TelemetryEnvelopeInputV1 {
            boundary: shifted_boundary,
            requests: &shifted_requests,
            ..telemetry_input
        })
        .expect("timing-shifted telemetry"),
        "raw timing changes remain visible in telemetry"
    );
    assert_ne!(
        telemetry,
        i14_telemetry_envelope_digest_v1(I14TelemetryEnvelopeInputV1 {
            clock_calibration_artifact: ContentHash([0x34; 32]),
            ..telemetry_input
        })
        .expect("changed calibration artifact"),
    );
    for (field, changed_watchdogs) in [
        (
            "watchdog kind",
            [
                I14WatchdogObservationV1 {
                    kind: I14WatchdogObservationKindV1::DeadlineSample,
                    ..watchdogs[0]
                },
                watchdogs[1],
            ],
        ),
        (
            "watchdog monotonic timestamp",
            [
                I14WatchdogObservationV1 {
                    monotonic_ns: watchdogs[0].monotonic_ns + 1,
                    ..watchdogs[0]
                },
                watchdogs[1],
            ],
        ),
    ] {
        assert_ne!(
            telemetry,
            i14_telemetry_envelope_digest_v1(I14TelemetryEnvelopeInputV1 {
                watchdog_observations: &changed_watchdogs,
                ..telemetry_input
            })
            .unwrap_or_else(|error| panic!("valid {field} mutation: {error:?}")),
            "telemetry identity omits {field}"
        );
    }

    let late_observed = I14CancellationRequestV1 {
        observation: Some(I14CancellationObservationV1 {
            logical_sequence: 120,
            monotonic_ns: 1_200,
            observing_tile_id: 41,
            latest_completed_boundary_ordinal: Some(10),
        }),
        ..late
    };
    let late_observed_changed_fields = I14CancellationRequestV1 {
        scope_root: 2,
        observation: Some(I14CancellationObservationV1 {
            observing_tile_id: 42,
            latest_completed_boundary_ordinal: Some(11),
            ..late_observed.observation.expect("late observation")
        }),
        ..late_observed
    };
    let late_trace_a = [observed, out_of_scope, late_observed];
    let late_trace_b = [observed, out_of_scope, late_observed_changed_fields];
    let late_telemetry_a = i14_telemetry_envelope_digest_v1(I14TelemetryEnvelopeInputV1 {
        requests: &late_trace_a,
        ..telemetry_input
    })
    .expect("first late raw trace");
    let late_telemetry_b = i14_telemetry_envelope_digest_v1(I14TelemetryEnvelopeInputV1 {
        requests: &late_trace_b,
        ..telemetry_input
    })
    .expect("second late raw trace");
    assert_ne!(late_telemetry_a, late_telemetry_b);
    for (field, changed_late_request) in [
        (
            "late request scope root",
            I14CancellationRequestV1 {
                scope_root: 2,
                ..late_observed
            },
        ),
        (
            "late observation tile",
            I14CancellationRequestV1 {
                observation: Some(I14CancellationObservationV1 {
                    observing_tile_id: 42,
                    ..late_observed.observation.expect("late observation")
                }),
                ..late_observed
            },
        ),
        (
            "late latest-completed boundary",
            I14CancellationRequestV1 {
                observation: Some(I14CancellationObservationV1 {
                    latest_completed_boundary_ordinal: Some(11),
                    ..late_observed.observation.expect("late observation")
                }),
                ..late_observed
            },
        ),
    ] {
        let changed_trace = [observed, out_of_scope, changed_late_request];
        let changed_telemetry = i14_telemetry_envelope_digest_v1(I14TelemetryEnvelopeInputV1 {
            requests: &changed_trace,
            ..telemetry_input
        })
        .unwrap_or_else(|error| panic!("valid {field} mutation: {error:?}"));
        assert_ne!(
            late_telemetry_a, changed_telemetry,
            "telemetry identity omits {field}"
        );
        assert_eq!(
            canonical.digest(),
            i14_canonical_terminal_result_digest_v1(I14CanonicalTerminalResultInputV1 {
                requests: &changed_trace,
                ..canonical_input
            })
            .unwrap_or_else(|error| panic!("canonical late {field} mutation: {error:?}")),
            "late {field} must remain outside the canonical cut"
        );
    }
    for latest_completed_boundary_ordinal in [None, Some(9)] {
        let invalid_late_observation = I14CancellationRequestV1 {
            observation: Some(I14CancellationObservationV1 {
                latest_completed_boundary_ordinal,
                ..late_observed.observation.expect("late observation")
            }),
            ..late_observed
        };
        let invalid_late_trace = [observed, out_of_scope, invalid_late_observation];
        assert_eq!(
            i14_telemetry_envelope_digest_v1(I14TelemetryEnvelopeInputV1 {
                requests: &invalid_late_trace,
                ..telemetry_input
            }),
            Err(I14TelemetryEnvelopeRefusalV1::RecordedTrace(
                I14TerminalCauseRefusalV1::ObservationBoundaryBehindCandidate {
                    request_id: late_observed.request_id,
                    latest_completed_boundary_ordinal,
                    candidate_boundary_ordinal: boundary.boundary_ordinal,
                }
            )),
            "a post-boundary observation cannot claim the candidate was incomplete"
        );
    }
    assert_eq!(
        canonical.digest(),
        i14_canonical_terminal_result_digest_v1(I14CanonicalTerminalResultInputV1 {
            requests: &late_trace_b,
            ..canonical_input
        })
        .expect("late field changes remain outside canonical cut")
    );

    assert_eq!(
        i14_telemetry_envelope_digest_v1(I14TelemetryEnvelopeInputV1 {
            requests: &requests_with_late_malformed,
            ..telemetry_input
        }),
        Err(I14TelemetryEnvelopeRefusalV1::RecordedTrace(
            I14TerminalCauseRefusalV1::DuplicateRequestId { request_id: 7 }
        ))
    );
    assert_eq!(
        i14_telemetry_envelope_digest_v1(I14TelemetryEnvelopeInputV1 {
            boundary: I14TerminalBoundaryV1 {
                completed: false,
                ..boundary
            },
            requests: &[],
            ..telemetry_input
        }),
        Err(I14TelemetryEnvelopeRefusalV1::CanonicalResult(
            I14CanonicalResultRefusalV1::BoundaryNotTerminal {
                decision: I14TerminalBoundaryDecisionV1::NotTerminal,
            }
        ))
    );
    assert_eq!(
        i14_telemetry_envelope_digest_v1(I14TelemetryEnvelopeInputV1 {
            requests: &[],
            ..telemetry_input
        }),
        Err(I14TelemetryEnvelopeRefusalV1::CanonicalResult(
            I14CanonicalResultRefusalV1::ExecutionDispositionMismatch {
                selected: I14ExecutionDisposition::Completed,
                receipt: I14ExecutionDisposition::Cancelled,
            }
        ))
    );
    let duplicate_watchdogs = [watchdogs[0], watchdogs[0]];
    assert_eq!(
        i14_telemetry_envelope_digest_v1(I14TelemetryEnvelopeInputV1 {
            watchdog_observations: &duplicate_watchdogs,
            ..telemetry_input
        }),
        Err(I14TelemetryEnvelopeRefusalV1::DuplicateWatchdogObservationId { observation_id: 2 })
    );
    let max_watchdogs = (0..I14_MAX_WATCHDOG_OBSERVATIONS_V1)
        .map(|index| I14WatchdogObservationV1 {
            observation_id: u64::try_from(index).expect("watchdog index fits u64"),
            kind: I14WatchdogObservationKindV1::DeadlineSample,
            monotonic_ns: u64::try_from(index).expect("watchdog time fits u64"),
        })
        .collect::<Vec<_>>();
    assert!(
        i14_telemetry_envelope_digest_v1(I14TelemetryEnvelopeInputV1 {
            watchdog_observations: &max_watchdogs,
            ..telemetry_input
        })
        .is_ok(),
        "watchdog cap is inclusive"
    );
    let too_many_watchdogs = vec![
        I14WatchdogObservationV1 {
            observation_id: 0,
            kind: I14WatchdogObservationKindV1::Poll,
            monotonic_ns: 0,
        };
        I14_MAX_WATCHDOG_OBSERVATIONS_V1 + 1
    ];
    assert_eq!(
        i14_telemetry_envelope_digest_v1(I14TelemetryEnvelopeInputV1 {
            watchdog_observations: &too_many_watchdogs,
            ..telemetry_input
        }),
        Err(I14TelemetryEnvelopeRefusalV1::TooManyWatchdogObservations {
            count: I14_MAX_WATCHDOG_OBSERVATIONS_V1 + 1,
            cap: I14_MAX_WATCHDOG_OBSERVATIONS_V1,
        })
    );
    let too_many_late_requests = vec![late; I14_MAX_CANCELLATION_REQUESTS_V1 + 1];
    assert_eq!(
        i14_canonical_terminal_result_v1(I14CanonicalTerminalResultInputV1 {
            requests: &too_many_late_requests,
            ..canonical_input
        }),
        Err(I14CanonicalResultRefusalV1::TerminalCause(
            I14TerminalCauseRefusalV1::TooManyRequests {
                count: I14_MAX_CANCELLATION_REQUESTS_V1 + 1,
                cap: I14_MAX_CANCELLATION_REQUESTS_V1,
            }
        ))
    );
}

#[test]
#[allow(clippy::too_many_lines)]
fn i14_v2_terminal_authority_binds_first_boundary_card_and_lifecycle() {
    let card = i14_core_cancellation_card();
    assert_eq!(card.tier(), I14CancellationTierV2::Core);
    assert_eq!(card.tier().child_and_observer_cap(), 32);
    assert_eq!(card.tier().request_to_observation_ns(), 250_000_000);
    assert_eq!(card.tier().watchdog_quantum_ns(), 25_000_000);
    assert_eq!(card.tier().trigger_to_drained_ns(), 2_000_000_000);
    assert_eq!(card.tier().drained_to_finalized_ns(), 2_000_000_000);
    assert_eq!(
        card.digest().to_hex(),
        I14_CANCELLATION_CARD_V2_KAT_HEX,
        "cancellation-card V2 byte layout drifted"
    );
    assert_eq!(
        i14_admit_cancellation_card_v2(I14CancellationCardInputV2 {
            campaign_wall_budget_ns: 0,
            ..i14_core_cancellation_card_input()
        }),
        Err(I14CancellationCardRefusalV2::ZeroCampaignWallBudget)
    );
    assert_eq!(
        i14_admit_cancellation_card_v2(I14CancellationCardInputV2 {
            logical_tile_response_bound_ns: 250_000_001,
            ..i14_core_cancellation_card_input()
        }),
        Err(
            I14CancellationCardRefusalV2::LogicalTileResponseExceedsTier {
                declared_ns: 250_000_001,
                cap_ns: 250_000_000,
            }
        )
    );
    assert_eq!(
        i14_admit_cancellation_card_v2(I14CancellationCardInputV2 {
            logical_memory_ceiling_bytes: 32 * 1_073_741_824 + 1,
            ..i14_core_cancellation_card_input()
        }),
        Err(
            I14CancellationCardRefusalV2::LogicalMemoryCeilingExceedsTier {
                declared_bytes: 32 * 1_073_741_824 + 1,
                cap_bytes: 32 * 1_073_741_824,
            }
        )
    );
    assert_eq!(
        i14_admit_cancellation_card_v2(I14CancellationCardInputV2 {
            total_resource_ceiling: 4_097,
            ..i14_core_cancellation_card_input()
        }),
        Err(
            I14CancellationCardRefusalV2::TotalResourceCeilingExceedsTier {
                declared: 4_097,
                cap: 4_096,
            }
        )
    );
    assert_eq!(
        i14_admit_cancellation_card_v2(I14CancellationCardInputV2 {
            maximum_resource_tile_quantum: 65,
            ..i14_core_cancellation_card_input()
        }),
        Err(
            I14CancellationCardRefusalV2::ResourceTileQuantumExceedsTier {
                declared: 65,
                cap: 64,
            }
        )
    );
    for malformed in [
        I14CancellationCardInputV2 {
            logical_tile_response_bound_ns: 0,
            ..i14_core_cancellation_card_input()
        },
        I14CancellationCardInputV2 {
            indivisible_item_response_bound_ns: 0,
            ..i14_core_cancellation_card_input()
        },
    ] {
        assert!(
            i14_admit_cancellation_card_v2(malformed).is_err(),
            "a zero demonstrated response bound has no cancellation authority"
        );
    }
    assert_eq!(
        i14_admit_cancellation_card_v2(I14CancellationCardInputV2 {
            external_child_catalog_digest: None,
            ..i14_core_cancellation_card_input()
        }),
        Err(I14CancellationCardRefusalV2::ExternalChildPolicyMismatch)
    );

    let scope = [1_u64, 2];
    let observers = [10_u64, 11];
    let boundary_0 = I14TerminalBoundaryV1 {
        boundary_ordinal: 0,
        logical_sequence: 20,
        monotonic_ns: 100_000_000,
        scope_ancestry: &scope,
        admitted_observer_tile_ids: &observers,
        infrastructure_failed: false,
        timed_out: false,
        budget_exhausted: false,
        completed: false,
    };
    let boundary_1 = I14TerminalBoundaryV1 {
        boundary_ordinal: 1,
        logical_sequence: 40,
        monotonic_ns: 300_000_000,
        ..boundary_0
    };
    let boundary_2 = I14TerminalBoundaryV1 {
        boundary_ordinal: 2,
        logical_sequence: 100,
        monotonic_ns: 800_000_000,
        completed: true,
        ..boundary_0
    };
    let records = [
        I14TerminalBoundaryRecordV2 {
            boundary: boundary_0,
            in_flight_children: 2,
            last_watchdog_poll_monotonic_ns: 90_000_000,
            total_resource_consumed: 100,
            work_items_remaining: 2,
            rejected_next_work_resource_cost: None,
            resource_ledger_prefix_digest: ContentHash([0xa0; 32]),
            work_frontier_prefix_digest: ContentHash([0xb0; 32]),
        },
        I14TerminalBoundaryRecordV2 {
            boundary: boundary_1,
            in_flight_children: 2,
            last_watchdog_poll_monotonic_ns: 290_000_000,
            total_resource_consumed: 200,
            work_items_remaining: 1,
            rejected_next_work_resource_cost: None,
            resource_ledger_prefix_digest: ContentHash([0xa1; 32]),
            work_frontier_prefix_digest: ContentHash([0xb1; 32]),
        },
        I14TerminalBoundaryRecordV2 {
            boundary: boundary_2,
            in_flight_children: 0,
            last_watchdog_poll_monotonic_ns: 790_000_000,
            total_resource_consumed: 300,
            work_items_remaining: 0,
            rejected_next_work_resource_cost: None,
            resource_ledger_prefix_digest: ContentHash([0xa2; 32]),
            work_frontier_prefix_digest: ContentHash([0xb2; 32]),
        },
    ];
    let request = I14CancellationRequestV1 {
        request_id: 7,
        scope_root: 2,
        logical_sequence: 30,
        requested_monotonic_ns: 200_000_000,
        observation_deadline_ns: 450_000_000,
        observation: Some(I14CancellationObservationV1 {
            logical_sequence: 50,
            monotonic_ns: 400_000_000,
            observing_tile_id: 10,
            latest_completed_boundary_ordinal: Some(1),
        }),
    };
    let requests = [request];
    let trace = I14TerminalBoundaryTraceV2 {
        logical_execution_verification_receipt_digest: ContentHash([0x61; 32]),
        cancellation_card: card,
        campaign_started_monotonic_ns: 0,
        boundaries: &records,
        requests: &requests,
    };
    let maximum_scope = (0..I14_MAX_SCOPE_ANCESTRY_V1)
        .map(|index| u64::try_from(index).expect("scope index fits u64"))
        .collect::<Vec<_>>();
    let maximum_observers = (0..card.tier().child_and_observer_cap())
        .map(|index| u64::try_from(index).expect("observer index fits u64"))
        .collect::<Vec<_>>();
    let maximum_catalog_records = [
        I14TerminalBoundaryRecordV2 {
            boundary: I14TerminalBoundaryV1 {
                scope_ancestry: &maximum_scope,
                admitted_observer_tile_ids: &maximum_observers,
                ..boundary_0
            },
            ..records[0]
        },
        I14TerminalBoundaryRecordV2 {
            boundary: I14TerminalBoundaryV1 {
                scope_ancestry: &maximum_scope,
                admitted_observer_tile_ids: &maximum_observers,
                ..boundary_1
            },
            ..records[1]
        },
        I14TerminalBoundaryRecordV2 {
            boundary: I14TerminalBoundaryV1 {
                scope_ancestry: &maximum_scope,
                admitted_observer_tile_ids: &maximum_observers,
                ..boundary_2
            },
            ..records[2]
        },
    ];
    assert!(
        matches!(
            i14_select_first_terminal_boundary_v2(I14TerminalBoundaryTraceV2 {
                boundaries: &maximum_catalog_records,
                ..trace
            }),
            Ok(I14TerminalTraceOutcomeV2::Selected(_))
        ),
        "the exact ancestry and tier observer caps remain admissible"
    );
    let overdeep_scope = (0..=I14_MAX_SCOPE_ANCESTRY_V1)
        .map(|index| u64::try_from(index).expect("scope index fits u64"))
        .collect::<Vec<_>>();
    let overdeep_records = [
        records[0],
        records[1],
        I14TerminalBoundaryRecordV2 {
            boundary: I14TerminalBoundaryV1 {
                scope_ancestry: &overdeep_scope,
                ..boundary_2
            },
            ..records[2]
        },
    ];
    assert_eq!(
        i14_select_first_terminal_boundary_v2(I14TerminalBoundaryTraceV2 {
            boundaries: &overdeep_records,
            ..trace
        }),
        Err(I14TerminalTraceRefusalV2::ScopeAncestryTooDeep {
            boundary_ordinal: 2,
            depth: I14_MAX_SCOPE_ANCESTRY_V1 + 1,
            cap: I14_MAX_SCOPE_ANCESTRY_V1,
        }),
        "V2 refuses ancestry cap+1 anywhere in the trace before cloning or comparing the path"
    );
    let observer_cap = card.tier().child_and_observer_cap();
    let excess_observers = (0..=observer_cap)
        .map(|index| u64::try_from(index).expect("observer index fits u64"))
        .collect::<Vec<_>>();
    let excess_observer_records = [
        records[0],
        records[1],
        I14TerminalBoundaryRecordV2 {
            boundary: I14TerminalBoundaryV1 {
                admitted_observer_tile_ids: &excess_observers,
                ..boundary_2
            },
            ..records[2]
        },
    ];
    assert_eq!(
        i14_select_first_terminal_boundary_v2(I14TerminalBoundaryTraceV2 {
            boundaries: &excess_observer_records,
            ..trace
        }),
        Err(I14TerminalTraceRefusalV2::TooManyObserverTilesForTier {
            boundary_ordinal: 2,
            count: observer_cap + 1,
            cap: observer_cap,
        }),
        "V2 refuses observer cap+1 anywhere in the trace before cloning or sorting the catalog"
    );
    let selection = match i14_select_first_terminal_boundary_v2(trace)
        .expect("valid genesis-to-first-terminal trace")
    {
        I14TerminalTraceOutcomeV2::Selected(selection) => selection,
        I14TerminalTraceOutcomeV2::Frontier(_) => panic!("terminal trace returned frontier"),
    };
    assert_eq!(selection.selected_index(), 2);
    assert_eq!(selection.boundary_count(), 3);
    assert_eq!(
        selection.prefix_digest().to_hex(),
        I14_TERMINAL_PREFIX_V2_KAT_HEX,
        "selected terminal-prefix V2 byte layout drifted"
    );
    assert_eq!(
        selection.decision(),
        I14TerminalBoundaryDecisionV1::Selected {
            disposition: I14ExecutionDisposition::Cancelled,
            request_id: Some(7),
        }
    );
    let out_of_scope_request = I14CancellationRequestV1 {
        request_id: 8,
        scope_root: 999,
        logical_sequence: 35,
        requested_monotonic_ns: 250_000_000,
        observation_deadline_ns: 500_000_000,
        observation: None,
    };
    let request_sensitive_trace = [request, out_of_scope_request];
    let request_sensitive_selection =
        match i14_select_first_terminal_boundary_v2(I14TerminalBoundaryTraceV2 {
            requests: &request_sensitive_trace,
            ..trace
        })
        .expect("semantically distinct request trace")
        {
            I14TerminalTraceOutcomeV2::Selected(selection) => selection,
            I14TerminalTraceOutcomeV2::Frontier(_) => panic!("terminal trace returned frontier"),
        };
    assert_eq!(request_sensitive_selection.decision(), selection.decision());
    assert_ne!(
        request_sensitive_selection.prefix_digest(),
        selection.prefix_digest(),
        "a decision-neutral request mutation must still change complete-prefix identity"
    );

    let lifecycle = I14TerminalLifecycleTraceV2 {
        execution_started: I14TimedLogicalEventV2 {
            logical_sequence: 1,
            monotonic_ns: 10_000_000,
        },
        drain_trigger: I14DrainTriggerV2::CancellationObserved { request_id: 7 },
        drain_started: I14TimedLogicalEventV2 {
            logical_sequence: 60,
            monotonic_ns: 500_000_000,
        },
        drained: I14TimedLogicalEventV2 {
            logical_sequence: 70,
            monotonic_ns: 600_000_000,
        },
        finalized: I14TimedLogicalEventV2 {
            logical_sequence: 80,
            monotonic_ns: 700_000_000,
        },
        active_children_at_drain_start: 2,
        drained_children: 2,
        child_lifecycle_semantic_trace_digest: ContentHash([0x77; 32]),
        child_lifecycle_raw_trace_digest: ContentHash([0x93; 32]),
        child_lifecycle_verification_receipt_digest: ContentHash([0x76; 32]),
        spawn_frontier_audit: Some(I14SpawnFrontierEvidenceV2 {
            request_id: Some(7),
            scheduler_semantic_trace_digest: ContentHash([0x77; 32]),
            last_child_spawn: Some(I14TimedLogicalEventV2 {
                logical_sequence: 25,
                monotonic_ns: 150_000_000,
            }),
            post_frontier_spawn_count: 0,
        }),
        watchdog_coverage: I14WatchdogCoverageV2 {
            poll_count: 33,
            first_poll_monotonic_ns: 10_000_000,
            last_poll_monotonic_ns: 790_000_000,
            maximum_poll_gap_ns: 25_000_000,
            watchdog_semantic_trace_digest: ContentHash([0x88; 32]),
            watchdog_raw_trace_digest: ContentHash([0x8e; 32]),
            watchdog_verification_receipt_digest: ContentHash([0x89; 32]),
        },
        tile_poll_coverage: I14TilePollCoverageV2 {
            admitted_tile_count: 3,
            fully_bracketed_tile_count: 3,
            before_item_zero_poll_observed: true,
            tile_poll_semantic_trace_digest: ContentHash([0x8c; 32]),
            tile_poll_raw_trace_digest: ContentHash([0x8f; 32]),
            tile_poll_verification_receipt_digest: ContentHash([0x8d; 32]),
        },
        external_heartbeat_coverage: I14ExternalHeartbeatCoverageV2 {
            admitted_external_children: 1,
            fully_covered_external_children: 1,
            maximum_concurrent_external_children: 1,
            heartbeat_count: 1,
            maximum_heartbeat_gap_ns: 25_000_000,
            heartbeat_semantic_trace_digest: ContentHash([0x8a; 32]),
            heartbeat_raw_trace_digest: ContentHash([0x92; 32]),
            heartbeat_verification_receipt_digest: ContentHash([0x8b; 32]),
            all_termination_acks_observed: true,
            atomic_publication_verified: true,
        },
        first_infrastructure_failure_onset: None,
        first_timeout_failure_onset_logical_sequence: None,
        first_timeout_failure_onset_monotonic_ns: None,
    };
    let canonical_input = I14CanonicalTerminalResultInputV2 {
        trace,
        lifecycle,
        terminal_status: i14_terminal_status(I14ExecutionDisposition::Cancelled),
        semantic_payload_digest: ContentHash([0x5a; 32]),
    };
    let canonical = i14_canonical_terminal_result_v2(canonical_input)
        .expect("valid authoritative canonical result");
    assert_eq!(canonical.selected_index(), 2);
    assert_eq!(canonical.boundary_count(), 3);
    assert_eq!(
        canonical.terminal_prefix_digest(),
        selection.prefix_digest()
    );
    assert_eq!(canonical.lifecycle().drained_children(), 2);
    assert_eq!(
        canonical
            .lifecycle()
            .child_lifecycle_semantic_trace_digest(),
        ContentHash([0x77; 32])
    );
    assert_eq!(
        canonical
            .lifecycle()
            .child_lifecycle_verification_receipt_digest(),
        ContentHash([0x76; 32])
    );
    assert_eq!(
        canonical.lifecycle().drain_trigger(),
        I14DrainTriggerV2::CancellationObserved { request_id: 7 }
    );
    assert!(!canonical.lifecycle().watchdog_slo_breached());
    assert!(!canonical.lifecycle().external_heartbeat_slo_breached());
    assert_eq!(
        canonical
            .lifecycle()
            .first_timeout_failure_onset_logical_sequence(),
        None,
        "an on-time cancellation has no lifecycle-timeout latch"
    );
    assert_eq!(
        i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
            lifecycle: I14TerminalLifecycleTraceV2 {
                child_lifecycle_semantic_trace_digest: ContentHash([0x78; 32]),
                ..lifecycle
            },
            ..canonical_input
        }),
        Err(I14CanonicalResultRefusalV2::Lifecycle(
            I14LifecycleRefusalV2::SpawnAuditTraceDigestMismatch
        )),
        "the mandatory spawn-frontier audit and child semantic trace must share identity"
    );
    assert_eq!(
        i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
            lifecycle: I14TerminalLifecycleTraceV2 {
                drained_children: lifecycle.active_children_at_drain_start + 1,
                ..lifecycle
            },
            ..canonical_input
        }),
        Err(I14CanonicalResultRefusalV2::Lifecycle(
            I14LifecycleRefusalV2::DrainedChildrenExceedActive {
                active: 2,
                drained: 3,
            }
        )),
        "a drained cut cannot account for children absent at drain start"
    );
    assert_eq!(
        i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
            lifecycle: I14TerminalLifecycleTraceV2 {
                external_heartbeat_coverage: I14ExternalHeartbeatCoverageV2 {
                    maximum_concurrent_external_children: 2,
                    ..lifecycle.external_heartbeat_coverage
                },
                ..lifecycle
            },
            ..canonical_input
        }),
        Err(I14CanonicalResultRefusalV2::Lifecycle(
            I14LifecycleRefusalV2::ExternalHeartbeatCoverageInconsistent
        )),
        "external peak concurrency cannot exceed the admitted population"
    );
    assert_eq!(
        i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
            lifecycle: I14TerminalLifecycleTraceV2 {
                watchdog_coverage: I14WatchdogCoverageV2 {
                    poll_count: 2,
                    first_poll_monotonic_ns: 780_000_000,
                    last_poll_monotonic_ns: 790_000_000,
                    maximum_poll_gap_ns: 11_000_000,
                    ..lifecycle.watchdog_coverage
                },
                ..lifecycle
            },
            ..canonical_input
        }),
        Err(I14CanonicalResultRefusalV2::Lifecycle(
            I14LifecycleRefusalV2::WatchdogCoverageInconsistent {
                poll_count: 2,
                span_ns: 10_000_000,
                maximum_poll_gap_ns: 11_000_000,
            }
        )),
        "the claimed maximum consecutive gap cannot exceed the complete trace span"
    );
    assert_eq!(
        canonical.digest(),
        i14_canonical_terminal_result_digest_v2(canonical_input)
            .expect("V2 digest convenience path")
    );
    assert_eq!(
        canonical.digest().to_hex(),
        I14_CANONICAL_TERMINAL_RESULT_V2_KAT_HEX,
        "canonical terminal-result V2 byte layout drifted"
    );

    let watchdog_samples = [
        I14WatchdogObservationV1 {
            observation_id: 2,
            kind: I14WatchdogObservationKindV1::ExternalHeartbeat,
            monotonic_ns: 750_000_000,
        },
        I14WatchdogObservationV1 {
            observation_id: 1,
            kind: I14WatchdogObservationKindV1::Poll,
            monotonic_ns: 500_000_000,
        },
    ];
    let telemetry_input = I14TelemetryEnvelopeInputV2 {
        terminal: canonical_input,
        watchdog_observations: &watchdog_samples,
        watchdog_samples_complete: false,
        late_event_tail: I14LateEventTailV2 {
            event_count: 0,
            first_logical_sequence: None,
            last_logical_sequence: None,
            semantic_trace_digest: ContentHash([0x9a; 32]),
            verification_receipt_digest: ContentHash([0x9b; 32]),
        },
        clock_calibration_artifact: ContentHash([0x99; 32]),
    };
    let telemetry =
        i14_telemetry_envelope_digest_v2(telemetry_input).expect("valid noncanonical telemetry");
    assert_eq!(
        telemetry.to_hex(),
        I14_TELEMETRY_ENVELOPE_V2_KAT_HEX,
        "telemetry-envelope V2 byte layout drifted"
    );
    let raw_root_only_lifecycle = I14TerminalLifecycleTraceV2 {
        child_lifecycle_raw_trace_digest: ContentHash([0xa3; 32]),
        watchdog_coverage: I14WatchdogCoverageV2 {
            watchdog_raw_trace_digest: ContentHash([0xa4; 32]),
            ..lifecycle.watchdog_coverage
        },
        tile_poll_coverage: I14TilePollCoverageV2 {
            tile_poll_raw_trace_digest: ContentHash([0xa5; 32]),
            ..lifecycle.tile_poll_coverage
        },
        external_heartbeat_coverage: I14ExternalHeartbeatCoverageV2 {
            heartbeat_raw_trace_digest: ContentHash([0xa6; 32]),
            ..lifecycle.external_heartbeat_coverage
        },
        ..lifecycle
    };
    let raw_root_only_terminal = I14CanonicalTerminalResultInputV2 {
        lifecycle: raw_root_only_lifecycle,
        ..canonical_input
    };
    assert_eq!(
        canonical.digest(),
        i14_canonical_terminal_result_digest_v2(raw_root_only_terminal)
            .expect("raw-root-only canonical mutation"),
        "raw child/watchdog/tile/heartbeat roots stay outside canonical identity"
    );
    assert_ne!(
        telemetry,
        i14_telemetry_envelope_digest_v2(I14TelemetryEnvelopeInputV2 {
            terminal: raw_root_only_terminal,
            ..telemetry_input
        })
        .expect("raw-root-only telemetry mutation"),
        "every raw evidence root must remain visible in telemetry identity"
    );
    assert_eq!(
        i14_telemetry_envelope_digest_v2(I14TelemetryEnvelopeInputV2 {
            watchdog_samples_complete: true,
            ..telemetry_input
        }),
        Err(I14TelemetryEnvelopeRefusalV2::CompleteWatchdogSampleMismatch),
        "a diagnostic subset cannot masquerade as the complete watchdog trace"
    );
    let mut complete_watchdog_samples = (0_u64..=32)
        .map(|index| I14WatchdogObservationV1 {
            observation_id: 100 + index,
            kind: I14WatchdogObservationKindV1::Poll,
            monotonic_ns: if index <= 28 {
                10_000_000 + 25_000_000 * index
            } else {
                710_000_000 + 20_000_000 * (index - 28)
            },
        })
        .collect::<Vec<_>>();
    complete_watchdog_samples.extend([
        I14WatchdogObservationV1 {
            observation_id: 200,
            kind: I14WatchdogObservationKindV1::ExternalHeartbeat,
            monotonic_ns: 750_000_000,
        },
        I14WatchdogObservationV1 {
            observation_id: 201,
            kind: I14WatchdogObservationKindV1::DeadlineSample,
            monotonic_ns: 760_000_000,
        },
    ]);
    let complete_watchdog_root = i14_watchdog_raw_trace_digest_v2(&complete_watchdog_samples)
        .expect("complete raw watchdog trace");
    assert_eq!(
        complete_watchdog_root.to_hex(),
        I14_WATCHDOG_RAW_TRACE_V2_KAT_HEX,
        "complete multi-kind raw-watchdog V2 byte layout drifted"
    );
    let mut reversed_watchdog_samples = complete_watchdog_samples.clone();
    reversed_watchdog_samples.reverse();
    assert_eq!(
        i14_watchdog_raw_trace_digest_v2(&reversed_watchdog_samples)
            .expect("presentation-reordered watchdog trace"),
        complete_watchdog_root,
        "watchdog presentation order carries no raw-trace identity"
    );
    assert_eq!(
        i14_watchdog_raw_trace_digest_v2(&[
            complete_watchdog_samples[0],
            complete_watchdog_samples[0],
        ]),
        Err(
            I14TelemetryEnvelopeRefusalV2::DuplicateWatchdogObservationId {
                observation_id: 100,
            }
        )
    );
    let over_cap_watchdog_samples =
        vec![complete_watchdog_samples[0]; I14_MAX_WATCHDOG_OBSERVATIONS_V1 + 1];
    assert_eq!(
        i14_watchdog_raw_trace_digest_v2(&over_cap_watchdog_samples),
        Err(I14TelemetryEnvelopeRefusalV2::TooManyWatchdogObservations {
            count: I14_MAX_WATCHDOG_OBSERVATIONS_V1 + 1,
            cap: I14_MAX_WATCHDOG_OBSERVATIONS_V1,
        })
    );
    let complete_watchdog_lifecycle = I14TerminalLifecycleTraceV2 {
        watchdog_coverage: I14WatchdogCoverageV2 {
            watchdog_raw_trace_digest: complete_watchdog_root,
            ..lifecycle.watchdog_coverage
        },
        ..lifecycle
    };
    let complete_watchdog_records = [
        I14TerminalBoundaryRecordV2 {
            last_watchdog_poll_monotonic_ns: 85_000_000,
            ..records[0]
        },
        I14TerminalBoundaryRecordV2 {
            last_watchdog_poll_monotonic_ns: 285_000_000,
            ..records[1]
        },
        records[2],
    ];
    let complete_watchdog_terminal = I14CanonicalTerminalResultInputV2 {
        trace: I14TerminalBoundaryTraceV2 {
            boundaries: &complete_watchdog_records,
            ..trace
        },
        lifecycle: complete_watchdog_lifecycle,
        ..canonical_input
    };
    let complete_watchdog_telemetry = I14TelemetryEnvelopeInputV2 {
        terminal: complete_watchdog_terminal,
        watchdog_observations: &complete_watchdog_samples,
        watchdog_samples_complete: true,
        ..telemetry_input
    };
    assert!(
        i14_telemetry_envelope_digest_v2(complete_watchdog_telemetry).is_ok(),
        "a complete summary-consistent trace with its derived raw root is admitted"
    );
    assert_eq!(
        i14_telemetry_envelope_digest_v2(I14TelemetryEnvelopeInputV2 {
            terminal: I14CanonicalTerminalResultInputV2 {
                lifecycle: complete_watchdog_lifecycle,
                ..canonical_input
            },
            watchdog_observations: &complete_watchdog_samples,
            watchdog_samples_complete: true,
            ..telemetry_input
        }),
        Err(
            I14TelemetryEnvelopeRefusalV2::CompleteWatchdogBoundarySnapshotMismatch {
                boundary_ordinal: 0,
                expected_last_poll_monotonic_ns: Some(85_000_000),
                found_last_poll_monotonic_ns: 90_000_000,
            }
        ),
        "a boundary cannot fabricate a newer poll snapshot than the complete raw trace"
    );
    let later_fabricated_watchdog_records = [
        complete_watchdog_records[0],
        records[1],
        complete_watchdog_records[2],
    ];
    assert_eq!(
        i14_telemetry_envelope_digest_v2(I14TelemetryEnvelopeInputV2 {
            terminal: I14CanonicalTerminalResultInputV2 {
                trace: I14TerminalBoundaryTraceV2 {
                    boundaries: &later_fabricated_watchdog_records,
                    ..trace
                },
                lifecycle: complete_watchdog_lifecycle,
                ..canonical_input
            },
            watchdog_observations: &complete_watchdog_samples,
            watchdog_samples_complete: true,
            ..telemetry_input
        }),
        Err(
            I14TelemetryEnvelopeRefusalV2::CompleteWatchdogBoundarySnapshotMismatch {
                boundary_ordinal: 1,
                expected_last_poll_monotonic_ns: Some(285_000_000),
                found_last_poll_monotonic_ns: 290_000_000,
            }
        ),
        "complete-stream reconciliation must inspect every boundary snapshot"
    );
    let wrong_root = ContentHash([0xfe; 32]);
    assert_eq!(
        i14_telemetry_envelope_digest_v2(I14TelemetryEnvelopeInputV2 {
            terminal: I14CanonicalTerminalResultInputV2 {
                lifecycle: I14TerminalLifecycleTraceV2 {
                    watchdog_coverage: I14WatchdogCoverageV2 {
                        watchdog_raw_trace_digest: wrong_root,
                        ..complete_watchdog_lifecycle.watchdog_coverage
                    },
                    ..complete_watchdog_lifecycle
                },
                ..complete_watchdog_terminal
            },
            ..complete_watchdog_telemetry
        }),
        Err(
            I14TelemetryEnvelopeRefusalV2::CompleteWatchdogRawTraceDigestMismatch {
                expected: complete_watchdog_root,
                found: wrong_root,
            }
        )
    );
    let mut interior_shifted_watchdog_samples = complete_watchdog_samples.clone();
    interior_shifted_watchdog_samples[29].monotonic_ns += 1_000_000;
    let interior_shifted_root =
        i14_watchdog_raw_trace_digest_v2(&interior_shifted_watchdog_samples)
            .expect("interior-shifted raw watchdog trace");
    assert_ne!(interior_shifted_root, complete_watchdog_root);
    assert_eq!(
        i14_telemetry_envelope_digest_v2(I14TelemetryEnvelopeInputV2 {
            watchdog_observations: &interior_shifted_watchdog_samples,
            ..complete_watchdog_telemetry
        }),
        Err(
            I14TelemetryEnvelopeRefusalV2::CompleteWatchdogRawTraceDigestMismatch {
                expected: interior_shifted_root,
                found: complete_watchdog_root,
            }
        ),
        "an interior raw-time mutation cannot hide behind unchanged count/endpoints/max-gap"
    );
    let late_tail = I14LateEventTailV2 {
        event_count: 1,
        first_logical_sequence: Some(101),
        last_logical_sequence: Some(101),
        semantic_trace_digest: ContentHash([0x9c; 32]),
        verification_receipt_digest: ContentHash([0x9d; 32]),
    };
    let late_telemetry = i14_telemetry_envelope_digest_v2(I14TelemetryEnvelopeInputV2 {
        late_event_tail: late_tail,
        ..telemetry_input
    })
    .expect("structurally valid verification-receipt-bound post-terminal telemetry tail");
    assert_ne!(late_telemetry, telemetry);
    assert_eq!(
        i14_telemetry_envelope_digest_v2(I14TelemetryEnvelopeInputV2 {
            late_event_tail: I14LateEventTailV2 {
                last_logical_sequence: Some(999),
                ..late_tail
            },
            ..telemetry_input
        }),
        Err(I14TelemetryEnvelopeRefusalV2::LateEventTailInconsistent),
        "a singleton tail must identify the same first and last event"
    );
    assert_eq!(
        i14_telemetry_envelope_digest_v2(I14TelemetryEnvelopeInputV2 {
            late_event_tail: I14LateEventTailV2 {
                first_logical_sequence: Some(100),
                last_logical_sequence: Some(100),
                ..late_tail
            },
            ..telemetry_input
        }),
        Err(
            I14TelemetryEnvelopeRefusalV2::LateEventTailNotAfterTerminal {
                first_logical_sequence: 100,
                terminal_logical_sequence: 100,
            }
        )
    );

    let shifted_boundaries = [
        I14TerminalBoundaryRecordV2 {
            boundary: I14TerminalBoundaryV1 {
                monotonic_ns: boundary_0.monotonic_ns + 1_000_000,
                ..boundary_0
            },
            last_watchdog_poll_monotonic_ns: 91_000_000,
            ..records[0]
        },
        I14TerminalBoundaryRecordV2 {
            boundary: I14TerminalBoundaryV1 {
                monotonic_ns: boundary_1.monotonic_ns + 1_000_000,
                ..boundary_1
            },
            last_watchdog_poll_monotonic_ns: 291_000_000,
            ..records[1]
        },
        I14TerminalBoundaryRecordV2 {
            boundary: I14TerminalBoundaryV1 {
                monotonic_ns: boundary_2.monotonic_ns + 1_000_000,
                ..boundary_2
            },
            last_watchdog_poll_monotonic_ns: 791_000_000,
            ..records[2]
        },
    ];
    let shifted_requests = [I14CancellationRequestV1 {
        requested_monotonic_ns: request.requested_monotonic_ns + 1_000_000,
        observation_deadline_ns: request.observation_deadline_ns + 1_000_000,
        observation: request
            .observation
            .map(|observation| I14CancellationObservationV1 {
                monotonic_ns: observation.monotonic_ns + 1_000_000,
                ..observation
            }),
        ..request
    }];
    let shifted_lifecycle = I14TerminalLifecycleTraceV2 {
        execution_started: I14TimedLogicalEventV2 {
            monotonic_ns: lifecycle.execution_started.monotonic_ns + 1_000_000,
            ..lifecycle.execution_started
        },
        drain_started: I14TimedLogicalEventV2 {
            monotonic_ns: lifecycle.drain_started.monotonic_ns + 1_000_000,
            ..lifecycle.drain_started
        },
        drained: I14TimedLogicalEventV2 {
            monotonic_ns: lifecycle.drained.monotonic_ns + 1_000_000,
            ..lifecycle.drained
        },
        finalized: I14TimedLogicalEventV2 {
            monotonic_ns: lifecycle.finalized.monotonic_ns + 1_000_000,
            ..lifecycle.finalized
        },
        spawn_frontier_audit: lifecycle.spawn_frontier_audit.map(|audit| {
            I14SpawnFrontierEvidenceV2 {
                last_child_spawn: audit.last_child_spawn.map(|event| I14TimedLogicalEventV2 {
                    monotonic_ns: event.monotonic_ns + 1_000_000,
                    ..event
                }),
                ..audit
            }
        }),
        child_lifecycle_raw_trace_digest: ContentHash([0x94; 32]),
        watchdog_coverage: I14WatchdogCoverageV2 {
            first_poll_monotonic_ns: 11_000_000,
            last_poll_monotonic_ns: 791_000_000,
            watchdog_raw_trace_digest: ContentHash([0x95; 32]),
            ..lifecycle.watchdog_coverage
        },
        tile_poll_coverage: I14TilePollCoverageV2 {
            tile_poll_raw_trace_digest: ContentHash([0x96; 32]),
            ..lifecycle.tile_poll_coverage
        },
        external_heartbeat_coverage: I14ExternalHeartbeatCoverageV2 {
            heartbeat_raw_trace_digest: ContentHash([0x97; 32]),
            ..lifecycle.external_heartbeat_coverage
        },
        ..lifecycle
    };
    let shifted_input = I14CanonicalTerminalResultInputV2 {
        trace: I14TerminalBoundaryTraceV2 {
            campaign_started_monotonic_ns: 1_000_000,
            boundaries: &shifted_boundaries,
            requests: &shifted_requests,
            ..trace
        },
        lifecycle: shifted_lifecycle,
        ..canonical_input
    };
    assert_eq!(
        canonical.digest(),
        i14_canonical_terminal_result_digest_v2(shifted_input).expect("timing-equivalent V2 trace"),
        "raw clock shifts preserving all semantic relations remain outside canonical identity"
    );
    assert_ne!(
        telemetry,
        i14_telemetry_envelope_digest_v2(I14TelemetryEnvelopeInputV2 {
            terminal: shifted_input,
            ..telemetry_input
        })
        .expect("timing-shifted telemetry"),
        "raw clock shifts remain visible in telemetry identity"
    );

    let terminal_children_records = [
        records[0],
        records[1],
        I14TerminalBoundaryRecordV2 {
            in_flight_children: 1,
            ..records[2]
        },
    ];
    assert_eq!(
        i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
            trace: I14TerminalBoundaryTraceV2 {
                boundaries: &terminal_children_records,
                ..trace
            },
            ..canonical_input
        }),
        Err(I14CanonicalResultRefusalV2::TerminalTrace(
            I14TerminalTraceRefusalV2::TerminalBoundaryHasInFlightChildren {
                boundary_ordinal: 2,
                count: 1,
            }
        )),
        "a terminal receipt cannot coexist with an undrained child"
    );

    let post_drained_live_children_lifecycle = I14TerminalLifecycleTraceV2 {
        drain_trigger: I14DrainTriggerV2::NonCancellationDrain,
        drain_started: I14TimedLogicalEventV2 {
            logical_sequence: 31,
            monotonic_ns: 200_000_000,
        },
        drained: I14TimedLogicalEventV2 {
            logical_sequence: 35,
            monotonic_ns: 250_000_000,
        },
        spawn_frontier_audit: Some(I14SpawnFrontierEvidenceV2 {
            request_id: None,
            ..lifecycle.spawn_frontier_audit.expect("spawn audit")
        }),
        ..lifecycle
    };
    assert_eq!(
        i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
            trace: I14TerminalBoundaryTraceV2 {
                requests: &[],
                ..trace
            },
            lifecycle: post_drained_live_children_lifecycle,
            terminal_status: i14_terminal_status(I14ExecutionDisposition::Completed),
            ..canonical_input
        }),
        Err(I14CanonicalResultRefusalV2::Lifecycle(
            I14LifecycleRefusalV2::BoundaryAfterDrainedHasInFlightChildren {
                boundary_ordinal: 1,
                count: 2,
            }
        )),
        "every boundary after the descendant-drained cut must remain child-free"
    );

    let exact_ceiling_completion_records = [
        records[0],
        records[1],
        I14TerminalBoundaryRecordV2 {
            total_resource_consumed: 4_096,
            work_items_remaining: 0,
            rejected_next_work_resource_cost: None,
            ..records[2]
        },
    ];
    let exact_ceiling_trace = I14TerminalBoundaryTraceV2 {
        boundaries: &exact_ceiling_completion_records,
        requests: &[],
        ..trace
    };
    let exact_ceiling_lifecycle = I14TerminalLifecycleTraceV2 {
        drain_trigger: I14DrainTriggerV2::NonCancellationDrain,
        child_lifecycle_semantic_trace_digest: ContentHash([0x90; 32]),
        child_lifecycle_raw_trace_digest: ContentHash([0x98; 32]),
        child_lifecycle_verification_receipt_digest: ContentHash([0x91; 32]),
        spawn_frontier_audit: Some(I14SpawnFrontierEvidenceV2 {
            request_id: None,
            scheduler_semantic_trace_digest: ContentHash([0x90; 32]),
            last_child_spawn: None,
            post_frontier_spawn_count: 0,
        }),
        ..lifecycle
    };
    let exact_ceiling_input = I14CanonicalTerminalResultInputV2 {
        trace: exact_ceiling_trace,
        lifecycle: exact_ceiling_lifecycle,
        terminal_status: i14_terminal_status(I14ExecutionDisposition::Completed),
        ..canonical_input
    };
    let exact_ceiling_result = i14_canonical_terminal_result_v2(exact_ceiling_input)
        .expect("completion exactly at the resource ceiling");
    assert_eq!(
        exact_ceiling_result.local_result().decision(),
        I14TerminalBoundaryDecisionV1::Selected {
            disposition: I14ExecutionDisposition::Completed,
            request_id: None,
        },
        "completion at the exact ceiling is not budget exhaustion"
    );
    let substituted_child_trace_result =
        i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
            lifecycle: I14TerminalLifecycleTraceV2 {
                child_lifecycle_semantic_trace_digest: ContentHash([0x92; 32]),
                child_lifecycle_raw_trace_digest: ContentHash([0x99; 32]),
                child_lifecycle_verification_receipt_digest: ContentHash([0x9a; 32]),
                spawn_frontier_audit: Some(I14SpawnFrontierEvidenceV2 {
                    scheduler_semantic_trace_digest: ContentHash([0x92; 32]),
                    ..exact_ceiling_lifecycle
                        .spawn_frontier_audit
                        .expect("spawn audit")
                }),
                ..exact_ceiling_lifecycle
            },
            ..exact_ceiling_input
        })
        .expect("identity-distinct complete child trace");
    assert_ne!(
        exact_ceiling_result.digest(),
        substituted_child_trace_result.digest(),
        "equal child counts cannot hide substitution of the child-lifecycle trace"
    );

    let budget_terminal = I14TerminalBoundaryRecordV2 {
        boundary: I14TerminalBoundaryV1 {
            budget_exhausted: true,
            completed: false,
            ..boundary_2
        },
        total_resource_consumed: 4_090,
        work_items_remaining: 1,
        rejected_next_work_resource_cost: Some(64),
        ..records[2]
    };
    let budget_records = [records[0], records[1], budget_terminal];
    let budget_trace = I14TerminalBoundaryTraceV2 {
        boundaries: &budget_records,
        requests: &[],
        ..trace
    };
    assert!(matches!(
        i14_select_first_terminal_boundary_v2(budget_trace),
        Ok(I14TerminalTraceOutcomeV2::Selected(selection))
            if selection.decision()
                == I14TerminalBoundaryDecisionV1::Selected {
                    disposition: I14ExecutionDisposition::BudgetExhausted,
                    request_id: None,
                }
    ));
    let unsupported_budget_records = [
        records[0],
        records[1],
        I14TerminalBoundaryRecordV2 {
            rejected_next_work_resource_cost: None,
            ..budget_terminal
        },
    ];
    assert_eq!(
        i14_select_first_terminal_boundary_v2(I14TerminalBoundaryTraceV2 {
            boundaries: &unsupported_budget_records,
            requests: &[],
            ..trace
        }),
        Err(I14TerminalTraceRefusalV2::BudgetExhaustionMismatch {
            boundary_ordinal: 2,
            expected: false,
            found: true,
        }),
        "a caller boolean without a rejected-next-work witness is not budget evidence"
    );

    let malformed_deadline_requests = [I14CancellationRequestV1 {
        observation_deadline_ns: request.observation_deadline_ns - 1,
        ..request
    }];
    assert_eq!(
        i14_select_first_terminal_boundary_v2(I14TerminalBoundaryTraceV2 {
            requests: &malformed_deadline_requests,
            ..trace
        }),
        Err(I14TerminalTraceRefusalV2::ObservationDeadlineSloMismatch {
            request_id: 7,
            declared_delta_ns: 249_999_999,
            required_delta_ns: 250_000_000,
        })
    );
    let pre_campaign_requests = [I14CancellationRequestV1 {
        requested_monotonic_ns: 0,
        observation_deadline_ns: 250_000_000,
        ..request
    }];
    assert_eq!(
        i14_select_first_terminal_boundary_v2(I14TerminalBoundaryTraceV2 {
            campaign_started_monotonic_ns: 1,
            requests: &pre_campaign_requests,
            ..trace
        }),
        Err(I14TerminalTraceRefusalV2::RequestBeforeCampaign {
            request_id: 7,
            requested_ns: 0,
            campaign_started_ns: 1,
        })
    );

    let trailing_after_completion = [
        records[0],
        I14TerminalBoundaryRecordV2 {
            boundary: I14TerminalBoundaryV1 {
                completed: true,
                ..boundary_1
            },
            in_flight_children: 0,
            work_items_remaining: 0,
            ..records[1]
        },
        I14TerminalBoundaryRecordV2 {
            boundary: I14TerminalBoundaryV1 {
                infrastructure_failed: true,
                ..boundary_2
            },
            ..records[2]
        },
    ];
    assert_eq!(
        i14_select_first_terminal_boundary_v2(I14TerminalBoundaryTraceV2 {
            boundaries: &trailing_after_completion,
            requests: &[],
            ..trace
        }),
        Err(I14TerminalTraceRefusalV2::BoundaryAfterFirstTerminal {
            first_terminal_ordinal: 1,
            trailing_ordinal: 2,
        }),
        "a later higher-priority cause cannot erase an earlier terminal boundary"
    );
    let early_observed_request = I14CancellationRequestV1 {
        observation: request
            .observation
            .map(|observation| I14CancellationObservationV1 {
                latest_completed_boundary_ordinal: Some(0),
                ..observation
            }),
        ..request
    };
    let early_observed_requests = [early_observed_request];
    let trailing_after_cancellation = [
        records[0],
        I14TerminalBoundaryRecordV2 {
            boundary: I14TerminalBoundaryV1 {
                logical_sequence: 60,
                monotonic_ns: 500_000_000,
                completed: true,
                ..boundary_1
            },
            in_flight_children: 0,
            work_items_remaining: 0,
            last_watchdog_poll_monotonic_ns: 490_000_000,
            ..records[1]
        },
        I14TerminalBoundaryRecordV2 {
            boundary: I14TerminalBoundaryV1 {
                infrastructure_failed: true,
                ..boundary_2
            },
            ..records[2]
        },
    ];
    assert_eq!(
        i14_select_first_terminal_boundary_v2(I14TerminalBoundaryTraceV2 {
            boundaries: &trailing_after_cancellation,
            requests: &early_observed_requests,
            ..trace
        }),
        Err(I14TerminalTraceRefusalV2::BoundaryAfterFirstTerminal {
            first_terminal_ordinal: 1,
            trailing_ordinal: 2,
        }),
        "earlier observed cancellation wins before a later infrastructure failure"
    );
    match i14_select_first_terminal_boundary_v2(I14TerminalBoundaryTraceV2 {
        boundaries: &records[..2],
        requests: &[],
        ..trace
    })
    .expect("valid nonterminal frontier")
    {
        I14TerminalTraceOutcomeV2::Frontier(frontier) => {
            assert_eq!(frontier.boundary_count(), 2);
            assert_eq!(frontier.last_boundary_ordinal(), 1);
        }
        I14TerminalTraceOutcomeV2::Selected(_) => panic!("frontier selected a terminal"),
    }

    let slow_terminal = I14TerminalBoundaryRecordV2 {
        boundary: I14TerminalBoundaryV1 {
            monotonic_ns: 3_000_000_000,
            ..boundary_2
        },
        last_watchdog_poll_monotonic_ns: 2_990_000_000,
        ..records[2]
    };
    let slow_records = [records[0], records[1], slow_terminal];
    let slow_lifecycle = I14TerminalLifecycleTraceV2 {
        drained: I14TimedLogicalEventV2 {
            monotonic_ns: 2_500_000_001,
            ..lifecycle.drained
        },
        finalized: I14TimedLogicalEventV2 {
            monotonic_ns: 2_700_000_001,
            ..lifecycle.finalized
        },
        watchdog_coverage: I14WatchdogCoverageV2 {
            poll_count: 121,
            last_poll_monotonic_ns: 2_990_000_000,
            ..lifecycle.watchdog_coverage
        },
        first_timeout_failure_onset_logical_sequence: Some(100),
        first_timeout_failure_onset_monotonic_ns: Some(2_400_000_001),
        ..lifecycle
    };
    let exact_observed_trigger_lifecycle = I14TerminalLifecycleTraceV2 {
        drained: I14TimedLogicalEventV2 {
            monotonic_ns: 2_400_000_000,
            ..lifecycle.drained
        },
        finalized: I14TimedLogicalEventV2 {
            monotonic_ns: 2_600_000_000,
            ..lifecycle.finalized
        },
        watchdog_coverage: slow_lifecycle.watchdog_coverage,
        ..lifecycle
    };
    let exact_observed_trigger =
        i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
            trace: I14TerminalBoundaryTraceV2 {
                boundaries: &slow_records,
                ..trace
            },
            lifecycle: exact_observed_trigger_lifecycle,
            ..canonical_input
        })
        .expect("the observed-trigger drain SLO accepts its exact inclusive cap");
    assert!(
        !exact_observed_trigger
            .lifecycle()
            .trigger_to_drained_slo_breached()
    );
    let cap_plus_one_observed_trigger_lifecycle = I14TerminalLifecycleTraceV2 {
        drained: I14TimedLogicalEventV2 {
            monotonic_ns: 2_400_000_001,
            ..lifecycle.drained
        },
        finalized: I14TimedLogicalEventV2 {
            monotonic_ns: 2_600_000_001,
            ..lifecycle.finalized
        },
        watchdog_coverage: slow_lifecycle.watchdog_coverage,
        first_timeout_failure_onset_logical_sequence: Some(100),
        first_timeout_failure_onset_monotonic_ns: Some(2_400_000_001),
        ..lifecycle
    };
    let cap_plus_one_observed_records = [
        records[0],
        records[1],
        I14TerminalBoundaryRecordV2 {
            boundary: I14TerminalBoundaryV1 {
                timed_out: true,
                ..slow_terminal.boundary
            },
            ..slow_terminal
        },
    ];
    let cap_plus_one_observed_trigger =
        i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
            trace: I14TerminalBoundaryTraceV2 {
                boundaries: &cap_plus_one_observed_records,
                ..trace
            },
            lifecycle: cap_plus_one_observed_trigger_lifecycle,
            terminal_status: i14_terminal_status(I14ExecutionDisposition::TimedOut),
            ..canonical_input
        })
        .expect("observed-trigger cap-plus-one is honestly reflected as timed out");
    assert!(
        cap_plus_one_observed_trigger
            .lifecycle()
            .trigger_to_drained_slo_breached()
    );
    assert_eq!(
        cap_plus_one_observed_trigger
            .lifecycle()
            .first_timeout_failure_onset_logical_sequence(),
        Some(100),
        "the first selected boundary at or after the derived cap-plus-one onset is canonical"
    );
    for (onset_logical_sequence, onset_monotonic_ns) in
        [(Some(100_u64), None), (None, Some(2_400_000_001_u64))]
    {
        assert_eq!(
            i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
                lifecycle: I14TerminalLifecycleTraceV2 {
                    first_timeout_failure_onset_logical_sequence: onset_logical_sequence,
                    first_timeout_failure_onset_monotonic_ns: onset_monotonic_ns,
                    ..lifecycle
                },
                ..canonical_input
            }),
            Err(I14CanonicalResultRefusalV2::Lifecycle(
                I14LifecycleRefusalV2::FailureOnsetAxisPresenceMismatch {
                    cause: I14LifecycleCauseClassV2::TimedOut,
                    onset_logical_sequence,
                    onset_monotonic_ns,
                }
            )),
            "timeout latch identity and calibrated onset must be present or absent together"
        );
    }
    assert_eq!(
        i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
            trace: I14TerminalBoundaryTraceV2 {
                boundaries: &slow_records,
                ..trace
            },
            lifecycle: I14TerminalLifecycleTraceV2 {
                first_timeout_failure_onset_monotonic_ns: Some(2_500_000_001),
                ..slow_lifecycle
            },
            ..canonical_input
        }),
        Err(I14CanonicalResultRefusalV2::Lifecycle(
            I14LifecycleRefusalV2::FailureOnsetTimestampMismatch {
                cause: I14LifecycleCauseClassV2::TimedOut,
                expected_monotonic_ns: Some(2_400_000_001),
                found_monotonic_ns: Some(2_500_000_001),
            }
        )),
        "a caller cannot delay the derived timeout onset past the first failing nanosecond"
    );
    assert_eq!(
        i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
            trace: I14TerminalBoundaryTraceV2 {
                boundaries: &slow_records,
                ..trace
            },
            lifecycle: I14TerminalLifecycleTraceV2 {
                first_timeout_failure_onset_logical_sequence: Some(65),
                ..slow_lifecycle
            },
            ..canonical_input
        }),
        Err(I14CanonicalResultRefusalV2::Lifecycle(
            I14LifecycleRefusalV2::FailureOnsetLogicalSequenceMismatch {
                cause: I14LifecycleCauseClassV2::TimedOut,
                expected_logical_sequence: Some(100),
                found_logical_sequence: Some(65),
            }
        )),
        "a caller cannot invent a canonical logical event between timeout onset and its first latch boundary"
    );
    assert_eq!(
        i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
            trace: I14TerminalBoundaryTraceV2 {
                boundaries: &slow_records,
                ..trace
            },
            lifecycle: slow_lifecycle,
            ..canonical_input
        }),
        Err(I14CanonicalResultRefusalV2::Lifecycle(
            I14LifecycleRefusalV2::LifecycleFailureNotReflected {
                failure: I14LifecycleFailureV2::TriggerToDrainedDeadline,
                selected: I14ExecutionDisposition::Cancelled,
            }
        ))
    );
    let reflected_timeout_records = [
        records[0],
        records[1],
        I14TerminalBoundaryRecordV2 {
            boundary: I14TerminalBoundaryV1 {
                timed_out: true,
                ..slow_terminal.boundary
            },
            ..slow_terminal
        },
    ];
    assert!(
        i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
            trace: I14TerminalBoundaryTraceV2 {
                boundaries: &reflected_timeout_records,
                ..trace
            },
            lifecycle: slow_lifecycle,
            terminal_status: i14_terminal_status(I14ExecutionDisposition::TimedOut),
            ..canonical_input
        })
        .is_ok(),
        "honestly reported lifecycle deadline failure remains hashable evidence"
    );
    let exact_drained_to_finalized_lifecycle = I14TerminalLifecycleTraceV2 {
        finalized: I14TimedLogicalEventV2 {
            monotonic_ns: 2_600_000_000,
            ..lifecycle.finalized
        },
        watchdog_coverage: slow_lifecycle.watchdog_coverage,
        ..lifecycle
    };
    let exact_drained_to_finalized =
        i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
            trace: I14TerminalBoundaryTraceV2 {
                boundaries: &slow_records,
                ..trace
            },
            lifecycle: exact_drained_to_finalized_lifecycle,
            ..canonical_input
        })
        .expect("drained-to-finalized accepts its exact inclusive cap");
    assert!(
        !exact_drained_to_finalized
            .lifecycle()
            .drained_to_finalized_slo_breached()
    );
    let over_cap_drained_to_finalized_lifecycle = I14TerminalLifecycleTraceV2 {
        finalized: I14TimedLogicalEventV2 {
            monotonic_ns: 2_600_000_001,
            ..exact_drained_to_finalized_lifecycle.finalized
        },
        first_timeout_failure_onset_logical_sequence: Some(100),
        first_timeout_failure_onset_monotonic_ns: Some(2_600_000_001),
        ..exact_drained_to_finalized_lifecycle
    };
    assert_eq!(
        i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
            trace: I14TerminalBoundaryTraceV2 {
                boundaries: &slow_records,
                ..trace
            },
            lifecycle: over_cap_drained_to_finalized_lifecycle,
            ..canonical_input
        }),
        Err(I14CanonicalResultRefusalV2::Lifecycle(
            I14LifecycleRefusalV2::LifecycleFailureNotReflected {
                failure: I14LifecycleFailureV2::DrainedToFinalizedDeadline,
                selected: I14ExecutionDisposition::Cancelled,
            }
        ))
    );
    let reflected_drained_to_finalized =
        i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
            trace: I14TerminalBoundaryTraceV2 {
                boundaries: &reflected_timeout_records,
                ..trace
            },
            lifecycle: over_cap_drained_to_finalized_lifecycle,
            terminal_status: i14_terminal_status(I14ExecutionDisposition::TimedOut),
            ..canonical_input
        })
        .expect("drained-to-finalized cap-plus-one is honestly reflected");
    assert!(
        reflected_drained_to_finalized
            .lifecycle()
            .drained_to_finalized_slo_breached()
    );

    let post_spawn_lifecycle = I14TerminalLifecycleTraceV2 {
        spawn_frontier_audit: Some(I14SpawnFrontierEvidenceV2 {
            last_child_spawn: Some(I14TimedLogicalEventV2 {
                logical_sequence: 55,
                monotonic_ns: 450_000_000,
            }),
            post_frontier_spawn_count: 1,
            ..lifecycle.spawn_frontier_audit.expect("spawn audit")
        }),
        first_infrastructure_failure_onset: Some(i14_infrastructure_onset(
            I14InfrastructureFailureSourceV2::SpawnAfterFrontierClosure,
            56,
            450_000_000,
            0xa1,
        )),
        ..lifecycle
    };
    let out_of_order_spawn_lifecycle = I14TerminalLifecycleTraceV2 {
        spawn_frontier_audit: Some(I14SpawnFrontierEvidenceV2 {
            last_child_spawn: Some(lifecycle.finalized),
            ..lifecycle.spawn_frontier_audit.expect("spawn audit")
        }),
        ..lifecycle
    };
    assert_eq!(
        i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
            lifecycle: out_of_order_spawn_lifecycle,
            ..canonical_input
        }),
        Err(I14CanonicalResultRefusalV2::Lifecycle(
            I14LifecycleRefusalV2::SpawnEventOrder
        ))
    );
    assert_eq!(
        i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
            lifecycle: I14TerminalLifecycleTraceV2 {
                spawn_frontier_audit: Some(I14SpawnFrontierEvidenceV2 {
                    last_child_spawn: Some(I14TimedLogicalEventV2 {
                        logical_sequence: 75,
                        monotonic_ns: 650_000_000,
                    }),
                    ..lifecycle.spawn_frontier_audit.expect("spawn audit")
                }),
                ..lifecycle
            },
            ..canonical_input
        }),
        Err(I14CanonicalResultRefusalV2::Lifecycle(
            I14LifecycleRefusalV2::SpawnEventOrder
        )),
        "no child may spawn after the final descendant-drained cut"
    );
    assert_eq!(
        i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
            lifecycle: post_spawn_lifecycle,
            ..canonical_input
        }),
        Err(I14CanonicalResultRefusalV2::Lifecycle(
            I14LifecycleRefusalV2::LifecycleFailureNotReflected {
                failure: I14LifecycleFailureV2::SpawnAfterFrontierClosure,
                selected: I14ExecutionDisposition::Cancelled,
            }
        ))
    );
    let incomplete_external_lifecycle = I14TerminalLifecycleTraceV2 {
        external_heartbeat_coverage: I14ExternalHeartbeatCoverageV2 {
            fully_covered_external_children: 0,
            heartbeat_count: 0,
            maximum_heartbeat_gap_ns: 0,
            ..lifecycle.external_heartbeat_coverage
        },
        first_infrastructure_failure_onset: Some(i14_infrastructure_onset(
            I14InfrastructureFailureSourceV2::ExternalHeartbeatCoverage,
            60,
            500_000_000,
            0xa2,
        )),
        ..lifecycle
    };
    assert_eq!(
        i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
            lifecycle: incomplete_external_lifecycle,
            ..canonical_input
        }),
        Err(I14CanonicalResultRefusalV2::Lifecycle(
            I14LifecycleRefusalV2::LifecycleFailureNotReflected {
                failure: I14LifecycleFailureV2::ExternalHeartbeatCoverage,
                selected: I14ExecutionDisposition::Cancelled,
            }
        )),
        "missing external-child coverage cannot promote as a clean cancellation"
    );
    let incomplete_tile_poll_lifecycle = I14TerminalLifecycleTraceV2 {
        tile_poll_coverage: I14TilePollCoverageV2 {
            fully_bracketed_tile_count: 2,
            ..lifecycle.tile_poll_coverage
        },
        first_infrastructure_failure_onset: Some(i14_infrastructure_onset(
            I14InfrastructureFailureSourceV2::TilePollCoverage,
            60,
            500_000_000,
            0xa3,
        )),
        ..lifecycle
    };
    assert_eq!(
        i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
            lifecycle: incomplete_tile_poll_lifecycle,
            ..canonical_input
        }),
        Err(I14CanonicalResultRefusalV2::Lifecycle(
            I14LifecycleRefusalV2::LifecycleFailureNotReflected {
                failure: I14LifecycleFailureV2::TilePollCoverage,
                selected: I14ExecutionDisposition::Cancelled,
            }
        )),
        "missing before/after-tile polls cannot promote as a clean cancellation"
    );
    assert_eq!(
        i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
            lifecycle: I14TerminalLifecycleTraceV2 {
                tile_poll_coverage: I14TilePollCoverageV2 {
                    admitted_tile_count: 0,
                    fully_bracketed_tile_count: 0,
                    before_item_zero_poll_observed: true,
                    ..lifecycle.tile_poll_coverage
                },
                ..lifecycle
            },
            ..canonical_input
        }),
        Err(I14CanonicalResultRefusalV2::Lifecycle(
            I14LifecycleRefusalV2::TilePollCoverageInconsistent
        )),
        "the empty tile population has one canonical false before-item-zero encoding"
    );
    let infrastructure_records = [
        records[0],
        records[1],
        I14TerminalBoundaryRecordV2 {
            boundary: I14TerminalBoundaryV1 {
                infrastructure_failed: true,
                ..boundary_2
            },
            ..records[2]
        },
    ];
    assert_eq!(
        i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
            trace: I14TerminalBoundaryTraceV2 {
                boundaries: &infrastructure_records,
                ..trace
            },
            lifecycle,
            terminal_status: i14_terminal_status(I14ExecutionDisposition::InfrastructureFailed),
            ..canonical_input
        }),
        Err(I14CanonicalResultRefusalV2::Lifecycle(
            I14LifecycleRefusalV2::FailureOnsetPresenceMismatch {
                cause: I14LifecycleCauseClassV2::InfrastructureFailed,
                failure_present: true,
                onset_logical_sequence: None,
                onset_monotonic_ns: None,
            }
        )),
        "an infrastructure-selected boundary requires a receipt-bound onset witness"
    );
    let spurious_infrastructure_lifecycle = I14TerminalLifecycleTraceV2 {
        drain_trigger: I14DrainTriggerV2::InfrastructureFailure {
            onset_logical_sequence: 45,
        },
        first_infrastructure_failure_onset: Some(i14_infrastructure_onset(
            I14InfrastructureFailureSourceV2::Supervisor,
            45,
            350_000_000,
            0xa6,
        )),
        ..lifecycle
    };
    assert_eq!(
        i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
            lifecycle: spurious_infrastructure_lifecycle,
            ..canonical_input
        }),
        Err(I14CanonicalResultRefusalV2::Lifecycle(
            I14LifecycleRefusalV2::FailureOnsetPresenceMismatch {
                cause: I14LifecycleCauseClassV2::InfrastructureFailed,
                failure_present: false,
                onset_logical_sequence: Some(45),
                onset_monotonic_ns: Some(350_000_000),
            }
        )),
        "a generic infrastructure witness cannot fabricate an unselected failure"
    );
    let pre_execution_infrastructure_lifecycle = I14TerminalLifecycleTraceV2 {
        drain_trigger: I14DrainTriggerV2::InfrastructureFailure {
            onset_logical_sequence: 0,
        },
        first_infrastructure_failure_onset: Some(i14_infrastructure_onset(
            I14InfrastructureFailureSourceV2::Supervisor,
            0,
            lifecycle.execution_started.monotonic_ns,
            0xa6,
        )),
        ..lifecycle
    };
    assert_eq!(
        i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
            trace: I14TerminalBoundaryTraceV2 {
                boundaries: &infrastructure_records,
                ..trace
            },
            lifecycle: pre_execution_infrastructure_lifecycle,
            terminal_status: i14_terminal_status(I14ExecutionDisposition::InfrastructureFailed),
            ..canonical_input
        }),
        Err(I14CanonicalResultRefusalV2::Lifecycle(
            I14LifecycleRefusalV2::FailureOnsetOrder {
                cause: I14LifecycleCauseClassV2::InfrastructureFailed,
                onset_logical_sequence: 0,
                onset_monotonic_ns: 10_000_000,
            }
        )),
        "an infrastructure onset cannot precede the execution it purports to terminate"
    );
    let generic_infrastructure_lifecycle = I14TerminalLifecycleTraceV2 {
        drain_trigger: I14DrainTriggerV2::InfrastructureFailure {
            onset_logical_sequence: 45,
        },
        first_infrastructure_failure_onset: Some(i14_infrastructure_onset(
            I14InfrastructureFailureSourceV2::Supervisor,
            45,
            350_000_000,
            0xa7,
        )),
        ..lifecycle
    };
    let generic_infrastructure =
        i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
            trace: I14TerminalBoundaryTraceV2 {
                boundaries: &infrastructure_records,
                ..trace
            },
            lifecycle: generic_infrastructure_lifecycle,
            terminal_status: i14_terminal_status(I14ExecutionDisposition::InfrastructureFailed),
            ..canonical_input
        })
        .expect("receipt-bound supervisor failure needs no fabricated coverage defect");
    assert_eq!(
        generic_infrastructure
            .lifecycle()
            .first_infrastructure_failure_onset_logical_sequence(),
        Some(45)
    );
    assert_eq!(
        generic_infrastructure
            .lifecycle()
            .first_infrastructure_failure_source(),
        Some(I14InfrastructureFailureSourceV2::Supervisor)
    );
    assert_eq!(
        generic_infrastructure
            .lifecycle()
            .infrastructure_failure_verification_receipt_digest(),
        Some(ContentHash([0xa7; 32]))
    );
    let no_requests: [I14CancellationRequestV1; 0] = [];
    let coalesced_infrastructure_lifecycle = I14TerminalLifecycleTraceV2 {
        drain_trigger: I14DrainTriggerV2::InfrastructureFailure {
            onset_logical_sequence: lifecycle.drain_started.logical_sequence,
        },
        spawn_frontier_audit: Some(I14SpawnFrontierEvidenceV2 {
            request_id: None,
            ..lifecycle.spawn_frontier_audit.expect("spawn audit")
        }),
        first_infrastructure_failure_onset: Some(i14_infrastructure_onset(
            I14InfrastructureFailureSourceV2::Supervisor,
            lifecycle.drain_started.logical_sequence,
            lifecycle.drain_started.monotonic_ns,
            0xa7,
        )),
        ..lifecycle
    };
    let coalesced_infrastructure =
        i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
            trace: I14TerminalBoundaryTraceV2 {
                boundaries: &infrastructure_records,
                requests: &no_requests,
                ..trace
            },
            lifecycle: coalesced_infrastructure_lifecycle,
            terminal_status: i14_terminal_status(I14ExecutionDisposition::InfrastructureFailed),
            ..canonical_input
        })
        .expect("coalesced infrastructure-onset/drain-start event remains causal");
    assert_eq!(
        coalesced_infrastructure.lifecycle().drain_trigger(),
        I14DrainTriggerV2::InfrastructureFailure {
            onset_logical_sequence: lifecycle.drain_started.logical_sequence,
        }
    );
    assert_eq!(
        i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
            trace: I14TerminalBoundaryTraceV2 {
                boundaries: &infrastructure_records,
                requests: &no_requests,
                ..trace
            },
            lifecycle: I14TerminalLifecycleTraceV2 {
                first_infrastructure_failure_onset: Some(i14_infrastructure_onset(
                    I14InfrastructureFailureSourceV2::Supervisor,
                    lifecycle.drain_started.logical_sequence,
                    lifecycle.drain_started.monotonic_ns - 1,
                    0xa7,
                )),
                ..coalesced_infrastructure_lifecycle
            },
            terminal_status: i14_terminal_status(I14ExecutionDisposition::InfrastructureFailed),
            ..canonical_input
        }),
        Err(I14CanonicalResultRefusalV2::Lifecycle(
            I14LifecycleRefusalV2::LifecycleTimestampOrder
        )),
        "drain-start aliasing requires exact sequence-and-time coalescence"
    );
    for (event_kind, aliased_event) in [
        ("execution.started", lifecycle.execution_started),
        (
            "cancellation.requested",
            I14TimedLogicalEventV2 {
                logical_sequence: request.logical_sequence,
                monotonic_ns: request.requested_monotonic_ns,
            },
        ),
        (
            "child.spawn",
            lifecycle
                .spawn_frontier_audit
                .expect("spawn audit")
                .last_child_spawn
                .expect("last child spawn"),
        ),
        (
            "cancellation.observed",
            I14TimedLogicalEventV2 {
                logical_sequence: request.observation.expect("observation").logical_sequence,
                monotonic_ns: request.observation.expect("observation").monotonic_ns,
            },
        ),
        ("execution.drained", lifecycle.drained),
        ("execution.finalized", lifecycle.finalized),
        (
            "terminal.boundary",
            I14TimedLogicalEventV2 {
                logical_sequence: boundary_2.logical_sequence,
                monotonic_ns: boundary_2.monotonic_ns,
            },
        ),
    ] {
        assert_eq!(
            i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
                trace: I14TerminalBoundaryTraceV2 {
                    boundaries: &infrastructure_records,
                    ..trace
                },
                lifecycle: I14TerminalLifecycleTraceV2 {
                    first_infrastructure_failure_onset: Some(i14_infrastructure_onset(
                        I14InfrastructureFailureSourceV2::Supervisor,
                        aliased_event.logical_sequence,
                        aliased_event.monotonic_ns,
                        0xa7,
                    )),
                    ..generic_infrastructure_lifecycle
                },
                terminal_status: i14_terminal_status(I14ExecutionDisposition::InfrastructureFailed),
                ..canonical_input
            }),
            Err(I14CanonicalResultRefusalV2::Lifecycle(
                I14LifecycleRefusalV2::DuplicateTraceLogicalSequence {
                    logical_sequence: aliased_event.logical_sequence,
                }
            )),
            "infrastructure onset must not alias {event_kind}, even at the same timestamp"
        );
    }
    assert_eq!(
        i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
            trace: I14TerminalBoundaryTraceV2 {
                boundaries: &infrastructure_records,
                ..trace
            },
            lifecycle: I14TerminalLifecycleTraceV2 {
                first_infrastructure_failure_onset: Some(i14_infrastructure_onset(
                    I14InfrastructureFailureSourceV2::Supervisor,
                    45,
                    350_000_000,
                    0x00,
                )),
                ..generic_infrastructure_lifecycle
            },
            terminal_status: i14_terminal_status(I14ExecutionDisposition::InfrastructureFailed),
            ..canonical_input
        }),
        Err(I14CanonicalResultRefusalV2::Lifecycle(
            I14LifecycleRefusalV2::InfrastructureFailureVerificationReceiptDigestZero {
                source: I14InfrastructureFailureSourceV2::Supervisor,
            }
        ))
    );
    assert_eq!(
        i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
            trace: I14TerminalBoundaryTraceV2 {
                boundaries: &infrastructure_records,
                ..trace
            },
            lifecycle: I14TerminalLifecycleTraceV2 {
                first_infrastructure_failure_onset: Some(i14_infrastructure_onset(
                    I14InfrastructureFailureSourceV2::WatchdogCoverage,
                    45,
                    350_000_000,
                    0xa8,
                )),
                ..generic_infrastructure_lifecycle
            },
            terminal_status: i14_terminal_status(I14ExecutionDisposition::InfrastructureFailed),
            ..canonical_input
        }),
        Err(I14CanonicalResultRefusalV2::Lifecycle(
            I14LifecycleRefusalV2::InfrastructureFailureSourceMismatch {
                source: I14InfrastructureFailureSourceV2::WatchdogCoverage,
            }
        ))
    );
    let slow_infrastructure_terminal = I14TerminalBoundaryRecordV2 {
        boundary: I14TerminalBoundaryV1 {
            infrastructure_failed: true,
            ..slow_terminal.boundary
        },
        ..slow_terminal
    };
    let slow_infrastructure_records = [records[0], records[1], slow_infrastructure_terminal];
    let exact_infrastructure_trigger_lifecycle = I14TerminalLifecycleTraceV2 {
        drained: I14TimedLogicalEventV2 {
            monotonic_ns: 2_350_000_000,
            ..lifecycle.drained
        },
        finalized: I14TimedLogicalEventV2 {
            monotonic_ns: 2_550_000_000,
            ..lifecycle.finalized
        },
        watchdog_coverage: slow_lifecycle.watchdog_coverage,
        ..generic_infrastructure_lifecycle
    };
    let exact_infrastructure_trigger =
        i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
            trace: I14TerminalBoundaryTraceV2 {
                boundaries: &slow_infrastructure_records,
                ..trace
            },
            lifecycle: exact_infrastructure_trigger_lifecycle,
            terminal_status: i14_terminal_status(I14ExecutionDisposition::InfrastructureFailed),
            ..canonical_input
        })
        .expect("infrastructure-trigger drain accepts its exact inclusive cap");
    assert!(
        !exact_infrastructure_trigger
            .lifecycle()
            .trigger_to_drained_slo_breached()
    );
    let over_cap_infrastructure_trigger_lifecycle = I14TerminalLifecycleTraceV2 {
        drained: I14TimedLogicalEventV2 {
            monotonic_ns: 2_350_000_001,
            ..exact_infrastructure_trigger_lifecycle.drained
        },
        first_timeout_failure_onset_logical_sequence: Some(100),
        first_timeout_failure_onset_monotonic_ns: Some(2_350_000_001),
        ..exact_infrastructure_trigger_lifecycle
    };
    assert_eq!(
        i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
            trace: I14TerminalBoundaryTraceV2 {
                boundaries: &slow_infrastructure_records,
                ..trace
            },
            lifecycle: over_cap_infrastructure_trigger_lifecycle,
            terminal_status: i14_terminal_status(I14ExecutionDisposition::InfrastructureFailed),
            ..canonical_input
        }),
        Err(I14CanonicalResultRefusalV2::Lifecycle(
            I14LifecycleRefusalV2::LifecycleFailureNotReflected {
                failure: I14LifecycleFailureV2::TriggerToDrainedDeadline,
                selected: I14ExecutionDisposition::InfrastructureFailed,
            }
        )),
        "an infrastructure cause alone cannot hide its separate drain timeout"
    );
    let reflected_slow_infrastructure_records = [
        records[0],
        records[1],
        I14TerminalBoundaryRecordV2 {
            boundary: I14TerminalBoundaryV1 {
                timed_out: true,
                ..slow_infrastructure_terminal.boundary
            },
            ..slow_infrastructure_terminal
        },
    ];
    let reflected_over_cap_infrastructure =
        i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
            trace: I14TerminalBoundaryTraceV2 {
                boundaries: &reflected_slow_infrastructure_records,
                ..trace
            },
            lifecycle: over_cap_infrastructure_trigger_lifecycle,
            terminal_status: i14_terminal_status(I14ExecutionDisposition::InfrastructureFailed),
            ..canonical_input
        })
        .expect("infrastructure cap-plus-one timeout is separately reflected");
    assert!(
        reflected_over_cap_infrastructure
            .lifecycle()
            .trigger_to_drained_slo_breached()
    );
    assert!(
        i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
            trace: I14TerminalBoundaryTraceV2 {
                boundaries: &infrastructure_records,
                ..trace
            },
            lifecycle: post_spawn_lifecycle,
            terminal_status: i14_terminal_status(I14ExecutionDisposition::InfrastructureFailed),
            ..canonical_input
        })
        .is_ok(),
        "honestly reported spawn-after-frontier failure remains hashable evidence"
    );
    assert!(
        i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
            trace: I14TerminalBoundaryTraceV2 {
                boundaries: &infrastructure_records,
                ..trace
            },
            lifecycle: incomplete_external_lifecycle,
            terminal_status: i14_terminal_status(I14ExecutionDisposition::InfrastructureFailed),
            ..canonical_input
        })
        .is_ok(),
        "reflected external-heartbeat failure remains retainable evidence"
    );
    let infrastructure_before_observation_lifecycle = I14TerminalLifecycleTraceV2 {
        drain_trigger: I14DrainTriggerV2::InfrastructureFailure {
            onset_logical_sequence: 45,
        },
        first_infrastructure_failure_onset: Some(i14_infrastructure_onset(
            I14InfrastructureFailureSourceV2::ExternalHeartbeatCoverage,
            45,
            350_000_000,
            0xa2,
        )),
        ..incomplete_external_lifecycle
    };
    let infrastructure_before_observation =
        i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
            trace: I14TerminalBoundaryTraceV2 {
                boundaries: &infrastructure_records,
                ..trace
            },
            lifecycle: infrastructure_before_observation_lifecycle,
            terminal_status: i14_terminal_status(I14ExecutionDisposition::InfrastructureFailed),
            ..canonical_input
        })
        .expect("infrastructure onset before cancellation observation triggers drain");
    assert_eq!(
        infrastructure_before_observation
            .lifecycle()
            .drain_trigger(),
        I14DrainTriggerV2::InfrastructureFailure {
            onset_logical_sequence: 45,
        }
    );
    let observation_before_infrastructure_lifecycle = I14TerminalLifecycleTraceV2 {
        first_infrastructure_failure_onset: Some(i14_infrastructure_onset(
            I14InfrastructureFailureSourceV2::ExternalHeartbeatCoverage,
            55,
            450_000_000,
            0xa2,
        )),
        ..incomplete_external_lifecycle
    };
    let observation_before_infrastructure =
        i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
            trace: I14TerminalBoundaryTraceV2 {
                boundaries: &infrastructure_records,
                ..trace
            },
            lifecycle: observation_before_infrastructure_lifecycle,
            terminal_status: i14_terminal_status(I14ExecutionDisposition::InfrastructureFailed),
            ..canonical_input
        })
        .expect("cancellation observation before infrastructure onset remains the drain trigger");
    assert_eq!(
        observation_before_infrastructure
            .lifecycle()
            .drain_trigger(),
        I14DrainTriggerV2::CancellationObserved { request_id: 7 }
    );
    let infrastructure_before_observation_input = I14CanonicalTerminalResultInputV2 {
        trace: I14TerminalBoundaryTraceV2 {
            boundaries: &infrastructure_records,
            ..trace
        },
        lifecycle: infrastructure_before_observation_lifecycle,
        terminal_status: i14_terminal_status(I14ExecutionDisposition::InfrastructureFailed),
        ..canonical_input
    };
    let shifted_infrastructure_onset_input = I14CanonicalTerminalResultInputV2 {
        lifecycle: I14TerminalLifecycleTraceV2 {
            first_infrastructure_failure_onset: Some(i14_infrastructure_onset(
                I14InfrastructureFailureSourceV2::ExternalHeartbeatCoverage,
                45,
                360_000_000,
                0xa2,
            )),
            ..infrastructure_before_observation_lifecycle
        },
        ..infrastructure_before_observation_input
    };
    assert_eq!(
        i14_canonical_terminal_result_digest_v2(infrastructure_before_observation_input)
            .expect("infrastructure-trigger canonical baseline"),
        i14_canonical_terminal_result_digest_v2(shifted_infrastructure_onset_input)
            .expect("raw infrastructure-onset shift"),
        "raw onset timing stays outside canonical identity when causality is unchanged"
    );
    assert_ne!(
        i14_telemetry_envelope_digest_v2(I14TelemetryEnvelopeInputV2 {
            terminal: infrastructure_before_observation_input,
            ..telemetry_input
        })
        .expect("infrastructure-trigger telemetry baseline"),
        i14_telemetry_envelope_digest_v2(I14TelemetryEnvelopeInputV2 {
            terminal: shifted_infrastructure_onset_input,
            ..telemetry_input
        })
        .expect("shifted infrastructure-onset telemetry"),
        "raw infrastructure-onset timing must remain visible in telemetry identity"
    );
    assert_eq!(
        i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
            trace: I14TerminalBoundaryTraceV2 {
                boundaries: &infrastructure_records,
                ..trace
            },
            lifecycle: I14TerminalLifecycleTraceV2 {
                drain_trigger: I14DrainTriggerV2::InfrastructureFailure {
                    onset_logical_sequence: 31,
                },
                first_infrastructure_failure_onset: Some(i14_infrastructure_onset(
                    I14InfrastructureFailureSourceV2::ExternalHeartbeatCoverage,
                    31,
                    210_000_000,
                    0xa2,
                )),
                ..incomplete_external_lifecycle
            },
            terminal_status: i14_terminal_status(I14ExecutionDisposition::InfrastructureFailed),
            ..canonical_input
        }),
        Err(I14CanonicalResultRefusalV2::Lifecycle(
            I14LifecycleRefusalV2::FailureNotLatchedAtFirstBoundary {
                cause: I14LifecycleCauseClassV2::InfrastructureFailed,
                onset_logical_sequence: 31,
                required_boundary_ordinal: 1,
                selected_boundary_ordinal: 2,
            }
        )),
        "an asynchronous failure must terminate at the first boundary after its onset"
    );

    let earlier_request_later_observation = I14CancellationRequestV1 {
        request_id: 7,
        logical_sequence: 25,
        requested_monotonic_ns: 150_000_000,
        observation_deadline_ns: 400_000_000,
        observation: Some(I14CancellationObservationV1 {
            logical_sequence: 55,
            monotonic_ns: 390_000_000,
            observing_tile_id: 10,
            latest_completed_boundary_ordinal: Some(1),
        }),
        ..request
    };
    let later_request_earlier_observation = I14CancellationRequestV1 {
        request_id: 8,
        logical_sequence: 35,
        requested_monotonic_ns: 250_000_000,
        observation_deadline_ns: 500_000_000,
        observation: Some(I14CancellationObservationV1 {
            logical_sequence: 45,
            monotonic_ns: 350_000_000,
            observing_tile_id: 11,
            latest_completed_boundary_ordinal: Some(1),
        }),
        ..request
    };
    let observation_race_requests = [
        earlier_request_later_observation,
        later_request_earlier_observation,
    ];
    let earliest_observation_lifecycle = I14TerminalLifecycleTraceV2 {
        drain_trigger: I14DrainTriggerV2::CancellationObserved { request_id: 8 },
        spawn_frontier_audit: Some(I14SpawnFrontierEvidenceV2 {
            request_id: Some(8),
            last_child_spawn: Some(I14TimedLogicalEventV2 {
                logical_sequence: 34,
                monotonic_ns: 200_000_000,
            }),
            ..lifecycle.spawn_frontier_audit.expect("spawn audit")
        }),
        ..lifecycle
    };
    let observation_race_input = I14CanonicalTerminalResultInputV2 {
        trace: I14TerminalBoundaryTraceV2 {
            requests: &observation_race_requests,
            ..trace
        },
        lifecycle: earliest_observation_lifecycle,
        ..canonical_input
    };
    let observation_race_result = i14_canonical_terminal_result_v2(observation_race_input)
        .expect("earliest effective observation triggers drain");
    assert_eq!(
        observation_race_result.local_result().decision(),
        I14TerminalBoundaryDecisionV1::Selected {
            disposition: I14ExecutionDisposition::Cancelled,
            request_id: Some(7),
        },
        "cause tie-break remains earliest request even when drain starts at an earlier observation"
    );
    assert_eq!(
        i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
            lifecycle: I14TerminalLifecycleTraceV2 {
                drain_trigger: I14DrainTriggerV2::CancellationObserved { request_id: 7 },
                spawn_frontier_audit: Some(I14SpawnFrontierEvidenceV2 {
                    request_id: Some(7),
                    ..earliest_observation_lifecycle
                        .spawn_frontier_audit
                        .expect("spawn audit")
                }),
                ..earliest_observation_lifecycle
            },
            ..observation_race_input
        }),
        Err(I14CanonicalResultRefusalV2::Lifecycle(
            I14LifecycleRefusalV2::DrainTriggerMismatch {
                expected: I14DrainTriggerV2::CancellationObserved { request_id: 8 },
                found: I14DrainTriggerV2::CancellationObserved { request_id: 7 },
            }
        )),
        "a later request cannot hide the first effective cancellation cut"
    );

    let missed_before_observed_requests = [
        I14CancellationRequestV1 {
            request_id: 7,
            logical_sequence: 25,
            requested_monotonic_ns: 150_000_000,
            observation_deadline_ns: 400_000_000,
            observation: None,
            ..request
        },
        I14CancellationRequestV1 {
            request_id: 8,
            logical_sequence: 35,
            requested_monotonic_ns: 250_000_000,
            observation_deadline_ns: 500_000_000,
            observation: Some(I14CancellationObservationV1 {
                logical_sequence: 45,
                monotonic_ns: 425_000_000,
                observing_tile_id: 11,
                latest_completed_boundary_ordinal: Some(1),
            }),
            ..request
        },
    ];
    let missed_before_observed_lifecycle = I14TerminalLifecycleTraceV2 {
        drain_trigger: I14DrainTriggerV2::ObservationTimeoutDrain { request_id: 7 },
        spawn_frontier_audit: Some(I14SpawnFrontierEvidenceV2 {
            request_id: Some(8),
            last_child_spawn: Some(I14TimedLogicalEventV2 {
                logical_sequence: 34,
                monotonic_ns: 200_000_000,
            }),
            ..lifecycle.spawn_frontier_audit.expect("spawn audit")
        }),
        ..lifecycle
    };
    let missed_before_observed_input = I14CanonicalTerminalResultInputV2 {
        trace: I14TerminalBoundaryTraceV2 {
            requests: &missed_before_observed_requests,
            ..trace
        },
        lifecycle: missed_before_observed_lifecycle,
        terminal_status: i14_terminal_status(I14ExecutionDisposition::TimedOut),
        ..canonical_input
    };
    let missed_before_observed = i14_canonical_terminal_result_v2(missed_before_observed_input)
        .expect("earlier missed deadline triggers timeout drain");
    assert_eq!(
        missed_before_observed.local_result().decision(),
        I14TerminalBoundaryDecisionV1::Selected {
            disposition: I14ExecutionDisposition::TimedOut,
            request_id: Some(7),
        }
    );
    assert_eq!(
        missed_before_observed.lifecycle().drain_trigger(),
        I14DrainTriggerV2::ObservationTimeoutDrain { request_id: 7 },
        "a later on-time observation cannot outrank an earlier effective deadline"
    );
    assert_eq!(
        i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
            lifecycle: I14TerminalLifecycleTraceV2 {
                drain_trigger: I14DrainTriggerV2::CancellationObserved { request_id: 8 },
                ..missed_before_observed_lifecycle
            },
            ..missed_before_observed_input
        }),
        Err(I14CanonicalResultRefusalV2::Lifecycle(
            I14LifecycleRefusalV2::DrainTriggerMismatch {
                expected: I14DrainTriggerV2::ObservationTimeoutDrain { request_id: 7 },
                found: I14DrainTriggerV2::CancellationObserved { request_id: 8 },
            }
        ))
    );

    let simultaneous_timeout_infrastructure_request = I14CancellationRequestV1 {
        observation: None,
        ..request
    };
    let simultaneous_timeout_infrastructure_requests =
        [simultaneous_timeout_infrastructure_request];
    let simultaneous_timeout_infrastructure_records = [
        records[0],
        records[1],
        I14TerminalBoundaryRecordV2 {
            boundary: I14TerminalBoundaryV1 {
                infrastructure_failed: true,
                timed_out: true,
                ..boundary_2
            },
            ..records[2]
        },
    ];
    let simultaneous_timeout_infrastructure_lifecycle = I14TerminalLifecycleTraceV2 {
        drain_trigger: I14DrainTriggerV2::InfrastructureFailure {
            onset_logical_sequence: 55,
        },
        spawn_frontier_audit: Some(I14SpawnFrontierEvidenceV2 {
            request_id: None,
            ..lifecycle.spawn_frontier_audit.expect("spawn audit")
        }),
        first_infrastructure_failure_onset: Some(i14_infrastructure_onset(
            I14InfrastructureFailureSourceV2::Supervisor,
            55,
            450_000_001,
            0xab,
        )),
        ..lifecycle
    };
    let simultaneous_timeout_infrastructure_input = I14CanonicalTerminalResultInputV2 {
        trace: I14TerminalBoundaryTraceV2 {
            boundaries: &simultaneous_timeout_infrastructure_records,
            requests: &simultaneous_timeout_infrastructure_requests,
            ..trace
        },
        lifecycle: simultaneous_timeout_infrastructure_lifecycle,
        terminal_status: i14_terminal_status(I14ExecutionDisposition::InfrastructureFailed),
        ..canonical_input
    };
    let simultaneous_timeout_infrastructure =
        i14_canonical_terminal_result_v2(simultaneous_timeout_infrastructure_input)
            .expect("infrastructure causal rank wins at the same effective nanosecond as timeout");
    assert_eq!(
        simultaneous_timeout_infrastructure
            .lifecycle()
            .drain_trigger(),
        I14DrainTriggerV2::InfrastructureFailure {
            onset_logical_sequence: 55,
        }
    );
    assert_eq!(
        simultaneous_timeout_infrastructure
            .lifecycle()
            .first_infrastructure_failure_onset_logical_sequence(),
        Some(55)
    );
    assert_eq!(
        simultaneous_timeout_infrastructure
            .lifecycle()
            .first_infrastructure_failure_source(),
        Some(I14InfrastructureFailureSourceV2::Supervisor)
    );
    assert_eq!(
        simultaneous_timeout_infrastructure
            .lifecycle()
            .infrastructure_failure_verification_receipt_digest(),
        Some(ContentHash([0xab; 32]))
    );
    let source_changed = i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
        lifecycle: I14TerminalLifecycleTraceV2 {
            first_infrastructure_failure_onset: Some(i14_infrastructure_onset(
                I14InfrastructureFailureSourceV2::Authentication,
                55,
                450_000_001,
                0xab,
            )),
            ..simultaneous_timeout_infrastructure_lifecycle
        },
        ..simultaneous_timeout_infrastructure_input
    })
    .expect("a structurally valid alternate generic failure source remains hashable");
    assert_ne!(
        source_changed.digest(),
        simultaneous_timeout_infrastructure.digest(),
        "the closed infrastructure source is canonical identity"
    );
    let verification_receipt_changed =
        i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
            lifecycle: I14TerminalLifecycleTraceV2 {
                first_infrastructure_failure_onset: Some(i14_infrastructure_onset(
                    I14InfrastructureFailureSourceV2::Supervisor,
                    55,
                    450_000_001,
                    0xac,
                )),
                ..simultaneous_timeout_infrastructure_lifecycle
            },
            ..simultaneous_timeout_infrastructure_input
        })
        .expect("a distinct nonzero verification-receipt identity remains structurally hashable");
    assert_ne!(
        verification_receipt_changed.digest(),
        simultaneous_timeout_infrastructure.digest(),
        "the independent infrastructure verification-receipt identity is canonical"
    );

    let unobserved_request = I14CancellationRequestV1 {
        observation: None,
        ..request
    };
    let unobserved_requests = [unobserved_request];
    let timeout_lifecycle = I14TerminalLifecycleTraceV2 {
        drain_trigger: I14DrainTriggerV2::ObservationTimeoutDrain { request_id: 7 },
        spawn_frontier_audit: Some(I14SpawnFrontierEvidenceV2 {
            request_id: None,
            last_child_spawn: Some(I14TimedLogicalEventV2 {
                logical_sequence: 55,
                monotonic_ns: 400_000_000,
            }),
            ..lifecycle.spawn_frontier_audit.expect("spawn audit")
        }),
        ..lifecycle
    };
    let timeout_input = I14CanonicalTerminalResultInputV2 {
        trace: I14TerminalBoundaryTraceV2 {
            requests: &unobserved_requests,
            ..trace
        },
        lifecycle: timeout_lifecycle,
        terminal_status: i14_terminal_status(I14ExecutionDisposition::TimedOut),
        ..canonical_input
    };
    assert!(
        i14_canonical_terminal_result_v2(timeout_input).is_ok(),
        "an unobserved request does not close the spawn frontier before timeout drain starts"
    );
    let exact_timeout_trigger_lifecycle = I14TerminalLifecycleTraceV2 {
        drain_trigger: I14DrainTriggerV2::ObservationTimeoutDrain { request_id: 7 },
        drained: I14TimedLogicalEventV2 {
            monotonic_ns: 2_450_000_001,
            ..lifecycle.drained
        },
        finalized: I14TimedLogicalEventV2 {
            monotonic_ns: 2_650_000_000,
            ..lifecycle.finalized
        },
        spawn_frontier_audit: timeout_lifecycle.spawn_frontier_audit,
        watchdog_coverage: slow_lifecycle.watchdog_coverage,
        ..lifecycle
    };
    let exact_timeout_trigger =
        i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
            trace: I14TerminalBoundaryTraceV2 {
                boundaries: &slow_records,
                requests: &unobserved_requests,
                ..trace
            },
            lifecycle: exact_timeout_trigger_lifecycle,
            terminal_status: i14_terminal_status(I14ExecutionDisposition::TimedOut),
            ..canonical_input
        })
        .expect("the timeout-trigger drain SLO accepts its exact inclusive cap");
    assert!(
        !exact_timeout_trigger
            .lifecycle()
            .trigger_to_drained_slo_breached()
    );
    assert_eq!(
        exact_timeout_trigger
            .lifecycle()
            .first_timeout_failure_onset_logical_sequence(),
        None,
        "the request-observation timeout is the drain cause, not a fabricated lifecycle-SLO onset"
    );
    let cap_plus_one_timeout_trigger_lifecycle = I14TerminalLifecycleTraceV2 {
        drained: I14TimedLogicalEventV2 {
            monotonic_ns: 2_450_000_002,
            ..exact_timeout_trigger_lifecycle.drained
        },
        finalized: I14TimedLogicalEventV2 {
            monotonic_ns: 2_650_000_002,
            ..exact_timeout_trigger_lifecycle.finalized
        },
        first_timeout_failure_onset_logical_sequence: Some(100),
        first_timeout_failure_onset_monotonic_ns: Some(2_450_000_002),
        ..exact_timeout_trigger_lifecycle
    };
    let cap_plus_one_timeout_trigger =
        i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
            trace: I14TerminalBoundaryTraceV2 {
                boundaries: &reflected_timeout_records,
                requests: &unobserved_requests,
                ..trace
            },
            lifecycle: cap_plus_one_timeout_trigger_lifecycle,
            terminal_status: i14_terminal_status(I14ExecutionDisposition::TimedOut),
            ..canonical_input
        })
        .expect("timeout-trigger cap-plus-one is honestly reflected as timed out");
    assert!(
        cap_plus_one_timeout_trigger
            .lifecycle()
            .trigger_to_drained_slo_breached()
    );
    assert_eq!(
        cap_plus_one_timeout_trigger
            .lifecycle()
            .first_timeout_failure_onset_logical_sequence(),
        Some(100),
        "the cap-plus-one timeout-trigger breach latches at the first later boundary"
    );

    let non_cancellation_lifecycle = I14TerminalLifecycleTraceV2 {
        drain_trigger: I14DrainTriggerV2::NonCancellationDrain,
        drained: I14TimedLogicalEventV2 {
            monotonic_ns: 2_500_000_000,
            ..lifecycle.drained
        },
        finalized: I14TimedLogicalEventV2 {
            monotonic_ns: 2_700_000_000,
            ..lifecycle.finalized
        },
        spawn_frontier_audit: Some(I14SpawnFrontierEvidenceV2 {
            request_id: None,
            ..lifecycle.spawn_frontier_audit.expect("spawn audit")
        }),
        watchdog_coverage: slow_lifecycle.watchdog_coverage,
        ..lifecycle
    };
    let non_cancellation_exact =
        i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
            trace: I14TerminalBoundaryTraceV2 {
                boundaries: &slow_records,
                requests: &[],
                ..trace
            },
            lifecycle: non_cancellation_lifecycle,
            terminal_status: i14_terminal_status(I14ExecutionDisposition::Completed),
            ..canonical_input
        })
        .expect("the non-cancellation drain SLO accepts its exact inclusive cap");
    assert!(
        !non_cancellation_exact
            .lifecycle()
            .trigger_to_drained_slo_breached()
    );

    let non_cancellation_timeout_records = [
        records[0],
        records[1],
        I14TerminalBoundaryRecordV2 {
            boundary: I14TerminalBoundaryV1 {
                timed_out: true,
                ..slow_terminal.boundary
            },
            ..slow_terminal
        },
    ];
    let non_cancellation_over_cap_lifecycle = I14TerminalLifecycleTraceV2 {
        drained: I14TimedLogicalEventV2 {
            monotonic_ns: 2_500_000_001,
            ..non_cancellation_lifecycle.drained
        },
        first_timeout_failure_onset_logical_sequence: Some(100),
        first_timeout_failure_onset_monotonic_ns: Some(2_500_000_001),
        ..non_cancellation_lifecycle
    };
    let non_cancellation_over_cap =
        i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
            trace: I14TerminalBoundaryTraceV2 {
                boundaries: &non_cancellation_timeout_records,
                requests: &[],
                ..trace
            },
            lifecycle: non_cancellation_over_cap_lifecycle,
            terminal_status: i14_terminal_status(I14ExecutionDisposition::TimedOut),
            ..canonical_input
        })
        .expect("non-cancellation cap-plus-one drain is honestly reported as timed out");
    assert!(
        non_cancellation_over_cap
            .lifecycle()
            .trigger_to_drained_slo_breached()
    );
    let timeout_pre_onset_boundary = I14TerminalBoundaryRecordV2 {
        boundary: I14TerminalBoundaryV1 {
            boundary_ordinal: 2,
            logical_sequence: 65,
            monotonic_ns: 2_500_000_000,
            completed: false,
            ..boundary_0
        },
        in_flight_children: 0,
        last_watchdog_poll_monotonic_ns: 2_490_000_000,
        total_resource_consumed: 225,
        work_items_remaining: 1,
        rejected_next_work_resource_cost: None,
        resource_ledger_prefix_digest: ContentHash([0xa8; 32]),
        work_frontier_prefix_digest: ContentHash([0xb8; 32]),
    };
    let timeout_at_onset_boundary = I14TerminalBoundaryRecordV2 {
        boundary: I14TerminalBoundaryV1 {
            boundary_ordinal: 3,
            logical_sequence: 90,
            monotonic_ns: 2_500_000_001,
            completed: false,
            ..boundary_0
        },
        in_flight_children: 0,
        last_watchdog_poll_monotonic_ns: 2_500_000_001,
        total_resource_consumed: 250,
        work_items_remaining: 1,
        rejected_next_work_resource_cost: None,
        resource_ledger_prefix_digest: ContentHash([0xa9; 32]),
        work_frontier_prefix_digest: ContentHash([0xb9; 32]),
    };
    let timeout_after_unlatched_boundary_records = [
        records[0],
        records[1],
        timeout_pre_onset_boundary,
        timeout_at_onset_boundary,
        I14TerminalBoundaryRecordV2 {
            boundary: I14TerminalBoundaryV1 {
                boundary_ordinal: 4,
                timed_out: true,
                ..slow_terminal.boundary
            },
            ..slow_terminal
        },
    ];
    assert_eq!(
        i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
            trace: I14TerminalBoundaryTraceV2 {
                boundaries: &timeout_after_unlatched_boundary_records,
                requests: &[],
                ..trace
            },
            lifecycle: I14TerminalLifecycleTraceV2 {
                finalized: I14TimedLogicalEventV2 {
                    logical_sequence: 95,
                    ..non_cancellation_over_cap_lifecycle.finalized
                },
                first_timeout_failure_onset_logical_sequence: Some(90),
                ..non_cancellation_over_cap_lifecycle
            },
            terminal_status: i14_terminal_status(I14ExecutionDisposition::TimedOut),
            ..canonical_input
        }),
        Err(I14CanonicalResultRefusalV2::Lifecycle(
            I14LifecycleRefusalV2::FailureNotLatchedAtFirstBoundary {
                cause: I14LifecycleCauseClassV2::TimedOut,
                onset_logical_sequence: 90,
                required_boundary_ordinal: 3,
                selected_boundary_ordinal: 4,
            }
        )),
        "a boundary one nanosecond before timeout onset is admissible, but a boundary stamped exactly at onset must latch it"
    );
    let slow_timeout_lifecycle = I14TerminalLifecycleTraceV2 {
        drain_trigger: I14DrainTriggerV2::ObservationTimeoutDrain { request_id: 7 },
        spawn_frontier_audit: Some(I14SpawnFrontierEvidenceV2 {
            request_id: None,
            ..slow_lifecycle.spawn_frontier_audit.expect("spawn audit")
        }),
        first_timeout_failure_onset_logical_sequence: Some(100),
        first_timeout_failure_onset_monotonic_ns: Some(2_450_000_002),
        ..slow_lifecycle
    };
    let slow_timeout = i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
        trace: I14TerminalBoundaryTraceV2 {
            boundaries: &reflected_timeout_records,
            requests: &unobserved_requests,
            ..trace
        },
        lifecycle: slow_timeout_lifecycle,
        terminal_status: i14_terminal_status(I14ExecutionDisposition::TimedOut),
        ..canonical_input
    })
    .expect("timeout trigger whose drain SLO is honestly reflected");
    assert!(
        slow_timeout.lifecycle().trigger_to_drained_slo_breached(),
        "unobserved timeout drain is measured from the first nanosecond after its inclusive deadline"
    );
    assert_eq!(
        i14_canonical_terminal_result_v2(I14CanonicalTerminalResultInputV2 {
            lifecycle: I14TerminalLifecycleTraceV2 {
                spawn_frontier_audit: Some(I14SpawnFrontierEvidenceV2 {
                    last_child_spawn: Some(I14TimedLogicalEventV2 {
                        logical_sequence: 65,
                        monotonic_ns: 550_000_000,
                    }),
                    post_frontier_spawn_count: 1,
                    ..timeout_lifecycle.spawn_frontier_audit.expect("spawn audit")
                }),
                first_infrastructure_failure_onset: Some(i14_infrastructure_onset(
                    I14InfrastructureFailureSourceV2::SpawnAfterFrontierClosure,
                    65,
                    550_000_000,
                    0xa4,
                )),
                ..timeout_lifecycle
            },
            ..timeout_input
        }),
        Err(I14CanonicalResultRefusalV2::Lifecycle(
            I14LifecycleRefusalV2::LifecycleFailureNotReflected {
                failure: I14LifecycleFailureV2::SpawnAfterFrontierClosure,
                selected: I14ExecutionDisposition::TimedOut,
            }
        )),
        "a spawn at or after timeout-triggered drain is an infrastructure failure"
    );

    let too_many_core_observers = (0_u64..33).collect::<Vec<_>>();
    let core_observer_record = [I14TerminalBoundaryRecordV2 {
        boundary: I14TerminalBoundaryV1 {
            admitted_observer_tile_ids: &too_many_core_observers,
            ..boundary_0
        },
        ..records[0]
    }];
    assert_eq!(
        i14_select_first_terminal_boundary_v2(I14TerminalBoundaryTraceV2 {
            boundaries: &core_observer_record,
            requests: &[],
            ..trace
        }),
        Err(I14TerminalTraceRefusalV2::TooManyObserverTilesForTier {
            boundary_ordinal: 0,
            count: 33,
            cap: 32,
        })
    );

    let max_card = i14_admit_cancellation_card_v2(I14CancellationCardInputV2 {
        tier: I14CancellationTierV2::Max,
        semantic_work_unit_digest: ContentHash([0x42; 32]),
        campaign_wall_budget_ns: 60_000_000_000,
        logical_tile_response_bound_ns: 1_000_000_000,
        indivisible_item_response_bound_ns: 100_000_000,
        external_heartbeat_bound_ns: 100_000_000,
        ..i14_core_cancellation_card_input()
    })
    .expect("valid Max cancellation card");
    assert_eq!(
        i14_admit_cancellation_card_v2(I14CancellationCardInputV2 {
            tier: I14CancellationTierV2::Max,
            campaign_wall_budget_ns: 86_400_000_000_000,
            semantic_work_unit_digest: ContentHash([0x43; 32]),
            logical_tile_response_bound_ns: 1_000_000_000,
            indivisible_item_response_bound_ns: 100_000_000,
            external_heartbeat_bound_ns: 100_000_000,
            ..i14_core_cancellation_card_input()
        }),
        Err(
            I14CancellationCardRefusalV2::CampaignWallBudgetExceedsTier {
                declared_ns: 86_400_000_000_000,
                cap_ns: 64_800_000_000_000,
            }
        )
    );
    assert!(
        i14_admit_cancellation_card_v2(I14CancellationCardInputV2 {
            tier: I14CancellationTierV2::MaxTheoremFalsifier,
            campaign_wall_budget_ns: 86_400_000_000_000,
            semantic_work_unit_digest: ContentHash([0x44; 32]),
            logical_tile_response_bound_ns: 1_000_000_000,
            indivisible_item_response_bound_ns: 100_000_000,
            external_heartbeat_bound_ns: 100_000_000,
            ..i14_core_cancellation_card_input()
        })
        .is_ok(),
        "the explicit theorem/falsifier subtype preserves the 24-hour Max envelope"
    );
    let max_observers = (0_u64..128).collect::<Vec<_>>();
    let max_record = [I14TerminalBoundaryRecordV2 {
        boundary: I14TerminalBoundaryV1 {
            monotonic_ns: 1,
            admitted_observer_tile_ids: &max_observers,
            ..boundary_0
        },
        in_flight_children: 128,
        last_watchdog_poll_monotonic_ns: 1,
        total_resource_consumed: 0,
        work_items_remaining: 1,
        rejected_next_work_resource_cost: None,
        resource_ledger_prefix_digest: ContentHash([0xc0; 32]),
        work_frontier_prefix_digest: ContentHash([0xd0; 32]),
    }];
    assert!(matches!(
        i14_select_first_terminal_boundary_v2(I14TerminalBoundaryTraceV2 {
            cancellation_card: max_card,
            boundaries: &max_record,
            requests: &[],
            ..trace
        }),
        Ok(I14TerminalTraceOutcomeV2::Frontier(_))
    ));

    let cap_records = (0..I14_MAX_TERMINAL_BOUNDARIES_V2)
        .map(|index| {
            let ordinal = u64::try_from(index).expect("boundary index fits u64");
            I14TerminalBoundaryRecordV2 {
                boundary: I14TerminalBoundaryV1 {
                    boundary_ordinal: ordinal,
                    logical_sequence: ordinal + 1,
                    monotonic_ns: ordinal * 1_000_000,
                    ..boundary_0
                },
                in_flight_children: 0,
                last_watchdog_poll_monotonic_ns: ordinal * 1_000_000,
                total_resource_consumed: ordinal.min(4_095),
                work_items_remaining: 1,
                rejected_next_work_resource_cost: None,
                resource_ledger_prefix_digest: ContentHash([0xc1; 32]),
                work_frontier_prefix_digest: ContentHash([0xd1; 32]),
            }
        })
        .collect::<Vec<_>>();
    assert!(matches!(
        i14_select_first_terminal_boundary_v2(I14TerminalBoundaryTraceV2 {
            boundaries: &cap_records,
            requests: &[],
            ..trace
        }),
        Ok(I14TerminalTraceOutcomeV2::Frontier(_))
    ));
    let arbitration_boundaries = (0_u64..65)
        .map(|ordinal| I14TerminalBoundaryRecordV2 {
            boundary: I14TerminalBoundaryV1 {
                boundary_ordinal: ordinal,
                logical_sequence: 20_000 + ordinal,
                monotonic_ns: 1_000_000_000 + ordinal,
                ..boundary_0
            },
            in_flight_children: 0,
            last_watchdog_poll_monotonic_ns: 1_000_000_000 + ordinal,
            total_resource_consumed: ordinal,
            work_items_remaining: 1,
            rejected_next_work_resource_cost: None,
            resource_ledger_prefix_digest: ContentHash([0xc2; 32]),
            work_frontier_prefix_digest: ContentHash([0xd2; 32]),
        })
        .collect::<Vec<_>>();
    let arbitration_requests = (0_u64..16_384)
        .map(|index| I14CancellationRequestV1 {
            request_id: index,
            scope_root: 999,
            logical_sequence: index + 1,
            requested_monotonic_ns: index,
            observation_deadline_ns: index + 250_000_000,
            observation: None,
        })
        .collect::<Vec<_>>();
    let arbitration_pairs = arbitration_boundaries.len() * arbitration_requests.len();
    assert!(arbitration_pairs > I14_MAX_TERMINAL_ARBITRATION_PAIRS_V2);
    assert_eq!(
        i14_select_first_terminal_boundary_v2(I14TerminalBoundaryTraceV2 {
            boundaries: &arbitration_boundaries,
            requests: &arbitration_requests,
            ..trace
        }),
        Err(I14TerminalTraceRefusalV2::ArbitrationWorkBudgetExceeded {
            boundary_count: 65,
            request_count: 16_384,
            pair_count: arbitration_pairs,
            cap: I14_MAX_TERMINAL_ARBITRATION_PAIRS_V2,
        }),
        "the reference selector refuses Cartesian sort work before entering its repeated V1 loop"
    );
    let mut over_cap_records = cap_records;
    over_cap_records.push(*over_cap_records.last().expect("cap trace is nonempty"));
    assert_eq!(
        i14_select_first_terminal_boundary_v2(I14TerminalBoundaryTraceV2 {
            boundaries: &over_cap_records,
            requests: &[],
            ..trace
        }),
        Err(I14TerminalTraceRefusalV2::TooManyBoundaries {
            count: I14_MAX_TERMINAL_BOUNDARIES_V2 + 1,
            cap: I14_MAX_TERMINAL_BOUNDARIES_V2,
        })
    );
}

#[test]
fn i14_retention_policy_is_typed_exhaustive_and_fail_closed() {
    let expected = [
        (
            I14ArtifactCategoryV1::Event,
            I14RetentionClassV1::EvidenceDurable,
        ),
        (
            I14ArtifactCategoryV1::Manifest,
            I14RetentionClassV1::EvidenceDurable,
        ),
        (
            I14ArtifactCategoryV1::AdjudicationReceipt,
            I14RetentionClassV1::EvidenceDurable,
        ),
        (
            I14ArtifactCategoryV1::Log,
            I14RetentionClassV1::EvidenceDurable,
        ),
        (
            I14ArtifactCategoryV1::OracleOutput,
            I14RetentionClassV1::EvidenceDurable,
        ),
        (
            I14ArtifactCategoryV1::ReplayCapsule,
            I14RetentionClassV1::EvidenceDurable,
        ),
        (
            I14ArtifactCategoryV1::MinimizedCounterexample,
            I14RetentionClassV1::FailurePermanent,
        ),
        (
            I14ArtifactCategoryV1::Refutation,
            I14RetentionClassV1::FailurePermanent,
        ),
        (
            I14ArtifactCategoryV1::FailureBundle,
            I14RetentionClassV1::FailurePermanent,
        ),
        (
            I14ArtifactCategoryV1::LicensedBytes,
            I14RetentionClassV1::GovernedRestricted,
        ),
        (
            I14ArtifactCategoryV1::SecretBytes,
            I14RetentionClassV1::GovernedRestricted,
        ),
        (
            I14ArtifactCategoryV1::SpecimenIdentity,
            I14RetentionClassV1::GovernedRestricted,
        ),
        (
            I14ArtifactCategoryV1::GovernedHoldoutBytes,
            I14RetentionClassV1::GovernedRestricted,
        ),
        (
            I14ArtifactCategoryV1::DerivedSensitiveSlice,
            I14RetentionClassV1::GovernedRestricted,
        ),
        (
            I14ArtifactCategoryV1::UnredactedDiagnosticSlice,
            I14RetentionClassV1::GovernedRestricted,
        ),
    ];
    assert_eq!(
        I14ArtifactCategoryV1::ALL,
        expected.map(|(category, _)| category)
    );
    for (category, expected_class) in expected {
        let rule = i14_retention_rule_v1(category);
        assert_eq!(rule.class, expected_class, "{category:?}");
        assert!(rule.complete_access_ledger, "{category:?}");
        assert!(rule.retain_retention_or_erasure_decision, "{category:?}");
        match expected_class {
            I14RetentionClassV1::EvidenceDurable | I14RetentionClassV1::FailurePermanent => {
                assert!(rule.sanitize_before_retention, "{category:?}");
                assert!(!rule.encrypted_capability_controlled, "{category:?}");
            }
            I14RetentionClassV1::GovernedRestricted => {
                assert!(!rule.sanitize_before_retention, "{category:?}");
                assert!(rule.encrypted_capability_controlled, "{category:?}");
            }
        }
    }
}

#[test]
fn i14_seed_algebra_is_fully_framed_endian_stable_and_lane_distinct() {
    let draft = i14_draft();
    let seeds = draft.explicits.seeds;
    for token in [
        "Philox 4x32-10",
        "d=BLAKE3::derive_key('org.frankensim.i14.fixture-stream.v1', alias_utf8)",
        "k0=LE32(d[0..4])",
        "k1=LE32(d[4..8])",
        "c0=low32(semantic_case_index)",
        "c1=high32(semantic_case_index)",
        "c2=low32(output_block_ordinal)",
        "c3=high32(output_block_ordinal)",
        "lane selection is the Philox output-word index",
        "never folded into the counter",
        "native-endian casts are forbidden",
        "Development indices 0..=4095",
        "PublicReplayCore 65536..=69631",
        "PublicReplayMax 131072..=135167",
        "falsifier indices 196608..=212991",
        "candidate-before-protected-access slot/envelope/discharge transaction",
        "GovernedPhysicalSlot permits logged least-privilege AuthorizedCalibration and AsBuiltModelInstantiation",
        "only to frozen calibration/model-input strata",
        "requires CandidateFrozen before candidate-side validation access or commitment opening",
        "joined-root/envelope/discharge transaction before adjudication",
        "GovernedStandards alone permits least-privilege licensed-input AuthorizedConstruction",
        "after GovernanceCommitted and before CandidateFrozen",
        "post-construction, pre-adjudication same-ID envelope/discharge transaction",
        "remains NoPromotionAuthority",
    ] {
        assert!(
            seeds.contains(token),
            "Five Explicits seed algebra omits {token}"
        );
    }
    for token in [
        "Core total resource ceilings",
        "Core logical poll tiles contain at most 64 graph/trace records",
        "asynchronous cancellation watchdog polls at intervals <=25 ms",
        "without ending or repartitioning a logical tile",
        "drain-trigger-to-drained <=2 s",
        "drained-to-finalized <=2 s",
        "Max total resource ceilings are four times",
        "Max logical poll tiles contain at most four times the Core item counts",
        "asynchronous watchdog polls at intervals <=100 ms",
        "without changing logical tile membership or order",
        "drain-trigger-to-drained <=8 s",
        "drained-to-finalized <=8 s",
        "structurally admitted receipt-bound infrastructure-failure onset whose receipt HELM/ledger authenticates",
        "CancellationObserved, ObservationTimeoutDrain, InfrastructureFailure, or NonCancellationDrain respectively",
        "Overall campaign wall-budget expiry is TimedOut",
        "heartbeat within the tier watchdog quantum",
        "logical tile lacking a demonstrated tier response bound",
    ] {
        assert!(
            draft.explicits.budgets.contains(token),
            "Five Explicits budgets omit {token}"
        );
    }
}

#[test]
fn i14_acceptance_arithmetic_is_preregistered_componentwise_and_fail_closed() {
    let draft = i14_draft();
    let policy = authored_spec(&draft, ACCEPTANCE_POLICY);
    for token in [
        "ACCEPTANCE_ARITHMETIC_V2",
        "pre-candidate FrozenManifest::amend successor",
        "one content-addressed AcceptanceCard per numeric claim and per independently scored component",
        "X=[x_lo,x_hi]",
        "HardUpper(u) is Supported iff x_hi<=u",
        "HardLower(l) iff x_lo>=l",
        "HardInterval(l,u)",
        "tau=max(a,rho*r)",
        "SoftEquality",
        "SoftUpper uses",
        "SoftLower uses",
        "SoftInterval uses",
        "Every hard predicate and every preregistered applicable soft component are conjoined",
        "max_i e_i only for the soft-score vector",
        "directed outward rounding",
        "minimal polynomial plus rational isolating interval",
        "Missing, undefined, nonfinite, invalid enclosure",
        "Exact-bit claims remain exact",
        "full soft component/enclosure vector rather than only a maximum",
    ] {
        assert!(policy.contains(token), "acceptance policy omits {token}");
    }
    for forbidden_hard_tolerance_laundering in [
        "max(0,x-u)/max(a,r)",
        "max(0,l-x)/max(a,r)",
        "max(l-x,0,x-u)/max(a,r)",
    ] {
        assert!(
            !policy.contains(forbidden_hard_tolerance_laundering),
            "hard comparator is incorrectly softened by {forbidden_hard_tolerance_laundering}"
        );
    }
    for row in &draft.obligations {
        assert!(row.decks.contains(&ACCEPTANCE_POLICY), "{}", row.leaf);
    }
    for claim in &draft.claims {
        match claim.tolerance {
            ToleranceSemantics::Exact => assert_eq!(claim.unit, "bit", "{}", claim.id),
            ToleranceSemantics::Interval { lo, hi } => {
                assert_eq!((lo, hi), (0.0, 1.0), "{}", claim.id);
                assert_eq!(claim.unit, "1", "{}", claim.id);
            }
            _ => panic!(
                "I14 claim {} uses an unauthorized tolerance class",
                claim.id
            ),
        }
    }
}

#[test]
#[allow(clippy::too_many_lines)]
fn i14_load_bearing_emc_physics_and_no_claim_clauses_are_pinned() {
    let draft = i14_draft();
    for token in [
        "canonical positive algebraic root alpha^2=2",
        "minimal polynomial x^2-2",
        "outward-rounded enclosures",
        "complex power is V I*",
        "one-half factor",
        "radiation",
        "LossOwnershipReceipt",
        "X(s)=integral_0^infinity x(t) exp(-s t) dt",
        "Re(s)<0",
        "Re(s)>0",
        "non-tangential",
        "boundary-axis poles",
    ] {
        assert!(
            draft.explicits.units.contains(token) || draft.explicits.versions.contains(token),
            "Five Explicits omit {token}"
        );
    }

    let native_harness = claim(&draft, "i14-harnessgraph-identity-connectivity");
    assert!(native_harness.statement.contains("native HarnessGraph"));
    assert!(native_harness.no_claim.contains("adapter correctness"));
    let adapter = claim(&draft, "i14-synthetic-ap242-adapter-mechanics");
    assert!(
        adapter
            .statement
            .contains("legally authored AP242-shaped supported subset")
    );
    assert!(adapter.statement.contains("every information loss"));
    assert!(
        adapter
            .no_claim
            .contains("not licensed-edition AP242 conformance")
    );

    let schema = claim(&draft, "i14-fullwave-problem-convention-admission");
    for token in [
        "FullWaveProblem",
        "EmConventionCard",
        "phasor",
        "RMS",
        "source",
        "relative phases",
        "global reference-phase",
        "outgoing-wave branch",
        "boundary",
        "dissipation/radiation power ownership",
    ] {
        assert!(schema.statement.contains(token), "schema omits {token}");
    }
    for token in [
        "exp(-i k r)/(4 pi r)",
        "simultaneously conjugates",
        "Hermitian dissipative inequality",
    ] {
        assert!(
            schema
                .hypotheses
                .iter()
                .any(|hypothesis| hypothesis.contains(token)),
            "schema hypotheses omit {token}"
        );
    }
    assert!(schema.no_claim.contains("solver accuracy"));
    let convention = authored_spec(&draft, EM_CONVENTION_CARD);
    for token in [
        "unique positive algebraic root alpha^2=2",
        "x(t)=Re{alpha X_+ exp(+i omega t)}",
        "global time-origin/reference-phase shift",
        "relative source/probe phase is physical frozen data",
        "X_-=conj(X_+)",
        "curl E=-i omega B",
        "curl H=J_impressed+(sigma+i omega epsilon)E",
        "RMS terminal complex power is V I*",
        "port current is positive into",
        "G_+(r)=exp(-i k r)/(4 pi r)",
        "Re(k)>=0 and Im(k)<=0",
        "conductivity has PSD Hermitian part",
        "epsilon and mu are negative semidefinite",
        "PML stretch signs",
        "complex-adjoint conjugate transpose",
        "certified outward enclosure",
    ] {
        assert!(
            convention.contains(token),
            "EM convention card omits {token}"
        );
    }

    let rlgc = claim(&draft, "i14-mtl-rlgc-operator-admission");
    for token in [
        "R+i omega L",
        "G+i omega C",
        "Z(s)",
        "Y(s)",
        "reference",
        "quotient",
        "band dissipativity",
    ] {
        assert!(
            rlgc.statement.contains(token) || rlgc.hypotheses.iter().any(|h| h.contains(token)),
            "RLGC authority omits {token}"
        );
    }
    for token in [
        "quasi-static or nondispersive",
        "passive internal-state/convolution realization",
        "frequency-derivative energy formulas",
        "Sampled nonnegative Hermitian parts",
        "boundary-axis pole",
        "PSD Hermitian residue",
    ] {
        assert!(
            rlgc.hypotheses.iter().any(|h| h.contains(token)),
            "RLGC hypotheses omit {token}"
        );
    }
    let mtl = claim(&draft, "i14-mtl-passive-causal-propagation");
    assert!(mtl.statement.contains("dV/dz=-(R+i omega L)I"));
    assert!(mtl.statement.contains("dI/dz=-(G+i omega C)V"));
    assert!(mtl.no_claim.contains("global passivity/causality theorem"));

    let peec = claim(&draft, "i14-peec-extraction-power-mor");
    for token in [
        "potential coefficient P",
        "capacitance C",
        "Lorenz",
        "retardation",
    ] {
        assert!(
            peec.statement.contains(token) || peec.hypotheses.iter().any(|h| h.contains(token)),
            "PEEC authority omits {token}"
        );
    }
    assert!(peec.kill.contains("P/C conflation"));

    let shield = claim(&draft, "i14-ground-bond-shield-current-closure");
    assert!(shield.statement.contains("exactly one owner"));
    assert!(shield.statement.contains("radiation boundary flux"));
    assert!(shield.no_claim.contains("legal EMC compliance"));

    let fullwave = claim(&draft, "i14-fullwave-feec-stability-energy");
    assert!(fullwave.statement.contains("discrete de Rham"));
    assert!(fullwave.hypotheses.iter().any(|h| h.contains("inf-sup")));
    assert!(
        fullwave
            .hypotheses
            .iter()
            .any(|h| h.contains("exact incidence is not used as a no-spurious-mode theorem"))
    );
    assert!(fullwave.no_claim.contains("PML"));
    assert!(authored_spec(&draft, "i14-fullwave-pml-dispersion").contains("late-time"));

    let bem = claim(&draft, "i14-exterior-bem-formulation-correctness");
    assert!(bem.statement.contains("dense"));
    for token in [
        "closed PEC",
        "PMCHWT",
        "Mueller",
        "JMCFIE",
        "open PEC screens",
        "screen EFIE",
        "edge-singularity",
    ] {
        assert!(bem.statement.contains(token), "BEM claim omits {token}");
    }
    assert!(bem.no_claim.contains("FMM acceleration"));
    assert!(authored_spec(&draft, "i14-fullwave-bem-scattering").contains("dense high-precision"));
    let fmm = claim(&draft, "i14-maxwell-fmm-acceleration-envelope");
    for token in [
        "dense",
        "certified outward error envelope",
        "matvec",
        "solved QoI",
        "crossover",
    ] {
        assert!(
            fmm.statement.contains(token)
                || fmm
                    .hypotheses
                    .iter()
                    .any(|hypothesis| hypothesis.contains(token)),
            "FMM acceleration claim omits {token}"
        );
    }
    assert!(fmm.no_claim.contains("discretization"));
    let fmm_leaf = draft
        .obligations
        .iter()
        .find(|row| row.leaf == "i14-fmm-acceleration-max")
        .expect("FMM acceleration leaf");
    assert!(fmm_leaf.g3_relations.iter().any(|relation| {
        relation.contains("cannot widen the declared certified envelope")
            && relation.contains("point error may be nonmonotone")
    }));

    let adjoint = claim(&draft, "i14-fixed-regime-adjoint-closure");
    assert!(adjoint.statement.contains("Complex states"));
    assert!(adjoint.statement.contains("Wirtinger convention"));
    assert!(adjoint.statement.contains("fixed EMC rung"));
    assert!(adjoint.hypotheses.iter().any(|h| h.contains("event")));
    assert!(adjoint.no_claim.contains("topology derivative"));

    let uq = claim(&draft, "i14-emc-uq-inference-mechanics");
    assert!(uq.hypotheses.iter().any(|h| h.contains("aleatory")));
    assert!(uq.hypotheses.iter().any(|h| h.contains("epistemic")));
    assert!(
        uq.hypotheses
            .iter()
            .any(|h| h.contains("public deterministic fixtures"))
    );
    assert!(uq.no_claim.contains("physical-population reliability"));
    let reliability = claim(&draft, "i14-governed-emc-reliability-validation");
    assert!(reliability.statement.contains("governed"));
    assert!(
        reliability
            .hypotheses
            .iter()
            .any(|h| h.contains("sampling frame"))
    );
    assert!(reliability.no_claim.contains("universal"));

    let safety = claim(&draft, "i14-emc-safety-case-integration");
    assert!(safety.hypotheses.iter().any(|h| h.contains("Unknown")));
    assert!(
        safety
            .activation
            .contains("complete atomic discharge transaction")
    );
    assert!(
        safety
            .activation
            .contains("independently adjudicated governed producer receipt")
    );
    assert!(safety.activation.contains("Unknown/NoPromotionAuthority"));
    for token in [
        "legal compliance",
        "product certification",
        "regulatory approval",
    ] {
        assert!(safety.no_claim.contains(token));
    }

    let standards = claim(&draft, "i14-governed-standards-crosswalk");
    assert!(standards.statement.contains("exact governed editions"));
    assert!(standards.statement.contains("loss-accounting"));
    assert!(standards.no_claim.contains("legal advice"));
    let laboratory = claim(&draft, "i14-governed-laboratory-emc-validation");
    assert!(
        laboratory
            .statement
            .contains("scoped physical-validation vector")
    );
    for token in [
        "AuthorizedCalibration",
        "AsBuiltModelInstantiation",
        "salted or equivalently hiding calibration/model-input content commitments",
        "hiding validation source-universe/frame commitment",
        "disjoint-membership commitment",
        "exact validation-selection algorithm",
        "pre-candidate secret-seed/VRF",
        "receive only opaque commitment identities",
        "never validation bytes, membership witnesses, labels, aggregates, derived statistics, selection outputs or commitment-opening material",
        "before any candidate-side protected validation access",
        "separately addressable calibration, model-input and validation roots",
        "membership proofs",
        "mutual disjoint-membership proof",
        "non-adaptive selection proof",
        "contamination receipt naming every audited principal and transitive capability",
    ] {
        assert!(
            laboratory
                .hypotheses
                .iter()
                .any(|hypothesis| hypothesis.contains(token)),
            "laboratory claim omits {token}"
        );
    }
    assert!(
        laboratory
            .hypotheses
            .iter()
            .any(|h| h.contains("AcceptanceCard"))
    );
    assert!(laboratory.no_claim.contains("population reliability"));
    let laboratory_leaf = draft
        .obligations
        .iter()
        .find(|row| row.leaf == "i14-laboratory-validation-max")
        .expect("laboratory validation leaf");
    assert!(laboratory_leaf.g3_relations.iter().any(|relation| {
        relation.contains("immutable timestamps")
            && relation.contains("sequence numbers")
            && relation.contains("acquisition-clock metadata")
            && relation.contains("canonical reconstruction")
    }));
    for token in [
        "root-only/no-capability pack-commitment registration",
        "GovernanceCommitted hiding calibration/model-input commitments",
        "AuthorizedCalibration least-privilege ledger",
        "AsBuiltModelInstantiation derivations",
        "CandidateFrozen candidate/checker/AcceptanceCard roots",
        "independent contamination receipt naming audited principals and transitive capabilities",
        "membership/disjointness/non-adaptive-selection proofs",
        "atomic RealizationCommitted authority-head advance",
        "validation-only reveal and Closed receipt",
    ] {
        assert!(
            laboratory_leaf.g0.contains(token),
            "laboratory G0 omits {token}"
        );
    }
    for token in [
        "giving any candidate builder/fitter/checker/threshold principal or transitive capability validation bytes",
        "commitment-opening material before CandidateFrozen",
        "validation source universe",
        "pre-candidate secret-seed/VRF commitment",
        "partially committing the two discharge envelopes and receipts",
    ] {
        assert!(
            laboratory_leaf
                .g3_relations
                .iter()
                .any(|relation| relation.contains(token)),
            "laboratory G3 omits {token}"
        );
    }
    for token in [
        "root-only/no-capability opaque pack-commitment registration",
        "GovernanceCommitted commitment freeze",
        "AuthorizedCalibration access",
        "calibration checking",
        "AsBuiltModelInstantiation",
        "CandidateFrozen sealing",
        "contamination checking",
        "validation realization",
        "atomic RealizationCommitted authority-head commit",
        "validation-only reveal",
        "Closed finalization",
        "forbid post-cancel access, reveal and authority-head advancement",
    ] {
        assert!(
            laboratory_leaf.g4_schedule.contains(token),
            "laboratory G4 omits {token}"
        );
    }
    assert_ordered_contract_tokens(
        "laboratory G4 lifecycle",
        laboratory_leaf.g4_schedule,
        &[
            "root-only/no-capability opaque pack-commitment registration",
            "GovernanceCommitted commitment freeze",
            "AuthorizedCalibration access",
            "calibration checking",
            "AsBuiltModelInstantiation",
            "CandidateFrozen sealing",
            "contamination checking",
            "validation realization",
            "separately addressable root/proof verification",
            "atomic RealizationCommitted authority-head commit",
            "validation-only reveal",
            "trace reconstruction",
            "component comparison",
            "adjudication",
            "Closed finalization",
        ],
    );
    let expected_lifecycle_events = [
        "laboratory.pack_commitment_registered",
        "laboratory.governance_committed",
        "laboratory.calibration_authorized",
        "laboratory.calibration_checked",
        "laboratory.model_instantiated",
        "laboratory.candidate_frozen",
        "laboratory.contamination_checked",
        "laboratory.validation_realized",
        "laboratory.authority_committed",
        "laboratory.validation_revealed",
        "laboratory.qoi_adjudicated",
        "laboratory.closed",
    ];
    assert_exact_lifecycle_event_membership(
        &draft,
        "i14-laboratory-validation-max",
        &expected_lifecycle_events,
    );
    let population = claim(&draft, "i14-production-bearing-population-reliability");
    assert!(population.statement.contains("production bearing"));
    assert!(
        population
            .hypotheses
            .iter()
            .any(|h| h.contains("censoring"))
    );
    assert!(population.no_claim.contains("warranty"));
}

#[test]
#[allow(clippy::too_many_lines)]
fn i14_governed_blind_and_standards_leaves_execute_their_authority_lifecycles() {
    let draft = i14_draft();
    let standards = draft
        .obligations
        .iter()
        .find(|row| row.leaf == "i14-standards-crosswalk-max")
        .expect("standards leaf");
    for token in [
        "GovernanceCommitted publisher/edition/corrigendum/license/source-root/scope identities",
        "least-privilege AuthorizedConstruction read/output ledger",
        "CandidateFrozen crosswalk/toolchain/checker/AcceptanceCard roots",
        "atomic StandardsAuthorityCommitted authority-head transaction",
        "independent StandardsAdjudicated reconstruction and Closed receipt",
        "no licensed byte or derived output is accessed before GovernanceCommitted",
        "partial/split retirement grants no authority",
    ] {
        assert!(standards.g0.contains(token), "standards G0 omits {token}");
    }
    for token in [
        "reordering GovernanceCommitted, AuthorizedConstruction, CandidateFrozen, StandardsAuthorityCommitted, StandardsAdjudicated and Closed",
        "changing publisher, edition, corrigendum, scope, principal, sandbox, procedure or disclosure filter",
        "same-ID-envelope-only installation, waiver-only retirement or a split authority-head transition",
    ] {
        assert!(
            standards
                .g3_relations
                .iter()
                .any(|relation| relation.contains(token)),
            "standards G3 omits {token}"
        );
    }
    for token in [
        "GovernanceCommitted freeze",
        "AuthorizedConstruction access",
        "CandidateFrozen sealing",
        "atomic StandardsAuthorityCommitted authority-head commit",
        "independent StandardsAdjudicated reconstruction",
        "Closed finalization",
        "forbid post-cancel licensed access, publication, adjudication and authority-head advancement",
    ] {
        assert!(
            standards.g4_schedule.contains(token),
            "standards G4 omits {token}"
        );
    }
    assert_ordered_contract_tokens(
        "standards G4 lifecycle",
        standards.g4_schedule,
        &[
            "GovernanceCommitted freeze",
            "AuthorizedConstruction access",
            "exact-edition admission",
            "clause/occurrence reconstruction",
            "loss/disclosure audit",
            "CandidateFrozen sealing",
            "same-ID envelope/receipt verification",
            "atomic StandardsAuthorityCommitted authority-head commit",
            "independent StandardsAdjudicated reconstruction",
            "Closed finalization",
        ],
    );
    assert_exact_lifecycle_event_membership(
        &draft,
        "i14-standards-crosswalk-max",
        &[
            "standards.governance_committed",
            "standards.construction_authorized",
            "standards.edition_admitted",
            "standards.clause_mapped",
            "adapter.loss_audited",
            "standards.candidate_frozen",
            "standards.authority_committed",
            "standards.adjudicated",
            "standards.closed",
        ],
    );

    let mitigation = draft
        .obligations
        .iter()
        .find(|row| row.leaf == "i14-robust-mitigation-max")
        .expect("robust mitigation leaf");
    for token in [
        "GovernanceCommitted generator/acquisition/access protocol and AcceptanceCards",
        "immutable CandidateFrozen candidate/model/toolchain/checker roots",
        "independent custodian realization",
        "atomic RealizationCommitted authority-head transaction",
        "RevealedForAdjudication one-shot access",
        "complete adjudication/access/attempt ledger and Closed receipt",
        "every candidate-side transitive capability never read a protected root",
        "partial/split retirement grants no authority",
    ] {
        assert!(mitigation.g0.contains(token), "mitigation G0 omits {token}");
    }
    for token in [
        "reordering GovernanceCommitted, CandidateFrozen, custodian realization, RealizationCommitted, RevealedForAdjudication and Closed",
        "changing the frozen generator, selection mechanism, candidate, checker, AcceptanceCard",
        "slot-only replacement, same-ID-envelope-only installation, waiver-only retirement",
    ] {
        assert!(
            mitigation
                .g3_relations
                .iter()
                .any(|relation| relation.contains(token)),
            "mitigation G3 omits {token}"
        );
    }
    for token in [
        "GovernanceCommitted freeze",
        "CandidateFrozen sealing",
        "independent custodian realization",
        "atomic RealizationCommitted authority-head commit",
        "RevealedForAdjudication access",
        "Closed finalization",
        "forbid post-cancel protected access, reveal and authority-head advancement",
    ] {
        assert!(
            mitigation.g4_schedule.contains(token),
            "mitigation G4 omits {token}"
        );
    }
    assert_ordered_contract_tokens(
        "mitigation G4 lifecycle",
        mitigation.g4_schedule,
        &[
            "GovernanceCommitted freeze",
            "design search",
            "adjoint/UQ evaluation",
            "fidelity escalation",
            "CandidateFrozen sealing",
            "independent custodian realization",
            "slot/envelope/receipt verification",
            "atomic RealizationCommitted authority-head commit",
            "RevealedForAdjudication access",
            "heldout adjudication",
            "guard checking",
            "Closed finalization",
        ],
    );
    assert_exact_lifecycle_event_membership(
        &draft,
        "i14-robust-mitigation-max",
        &[
            "mitigation.governance_committed",
            "mitigation.candidate_frozen",
            "mitigation.custodian_realized",
            "mitigation.authority_committed",
            "mitigation.revealed_for_adjudication",
            "mitigation.holdout_adjudicated",
            "mitigation.guard_checked",
            "mitigation.closed",
        ],
    );
}

#[test]
#[allow(clippy::too_many_lines)]
fn i14_governed_population_leaves_execute_atomic_candidate_first_lifecycles() {
    let draft = i14_draft();
    let emc = draft
        .obligations
        .iter()
        .find(|row| row.leaf == "i14-emc-reliability-validation-max")
        .expect("EMC reliability population leaf");
    for token in [
        "GovernanceCommitted population frame",
        "immutable CandidateFrozen model/toolchain/checker/AcceptanceCard roots",
        "independent custodian realization",
        "one atomic RealizationCommitted authority-head transaction",
        "one-shot RevealedForAdjudication",
        "complete access/attempt/adjudication ledger in Closed",
        "no candidate-side principal or transitive capability reads a protected population root",
        "partial/split retirement grants no authority",
    ] {
        assert!(emc.g0.contains(token), "EMC population G0 omits {token}");
    }
    for token in [
        "reordering GovernanceCommitted, CandidateFrozen, custodian realization, RealizationCommitted, RevealedForAdjudication and Closed",
        "changing the frozen population frame, acquisition/selection mechanism, candidate, checker, AcceptanceCard",
        "slot-only replacement, same-ID-envelope-only installation, waiver-only retirement, split authority-head transition, sequential reveal or post-reveal amendment",
    ] {
        assert!(
            emc.g3_relations
                .iter()
                .any(|relation| relation.contains(token)),
            "EMC population G3 omits {token}"
        );
    }
    assert_ordered_contract_tokens(
        "EMC population G4 lifecycle",
        emc.g4_schedule,
        &[
            "GovernanceCommitted freeze",
            "CandidateFrozen sealing",
            "independent custodian realization",
            "atomic RealizationCommitted authority-head commit",
            "RevealedForAdjudication access",
            "adjudication",
            "Closed finalization",
        ],
    );
    for token in [
        "forbid post-cancel protected access, reveal, adjudication and authority-head advancement",
        "unless the exact atomic stage had already committed",
        "complete frame/event/censoring/access/attempt/stopping ledgers",
    ] {
        assert!(
            emc.g4_schedule.contains(token),
            "EMC population G4 omits {token}"
        );
    }
    assert_exact_lifecycle_event_membership(
        &draft,
        "i14-emc-reliability-validation-max",
        &[
            "emc.governance_committed",
            "emc.candidate_frozen",
            "emc.custodian_realized",
            "emc.authority_committed",
            "emc.revealed_for_adjudication",
            "emc.reliability_event_audited",
            "emc.reliability_coverage_checked",
            "emc.reliability_adjudicated",
            "emc.closed",
        ],
    );

    let bearing = draft
        .obligations
        .iter()
        .find(|row| row.leaf == "i14-bearing-population-validation-max")
        .expect("bearing population leaf");
    for token in [
        "GovernanceCommitted sampling frame/estimand/configuration/event/inspection definitions",
        "immutable CandidateFrozen model/toolchain/checker/AcceptanceCard roots",
        "independent custodian joint realization",
        "two distinct same-ID discharge envelopes and verified receipts",
        "both waiver retirements and authority-head advance in one atomic RealizationCommitted coupled transaction",
        "one-shot joint RevealedForAdjudication",
        "complete access/attempt/adjudication ledger in Closed",
        "both packs commit or neither grants authority",
    ] {
        assert!(
            bearing.g0.contains(token),
            "bearing population G0 omits {token}"
        );
    }
    for token in [
        "reordering GovernanceCommitted, CandidateFrozen, custodian joint realization, coupled RealizationCommitted, joint RevealedForAdjudication and Closed",
        "changing the frozen frame, metrology, acquisition/selection mechanism, candidate, checker, AcceptanceCard",
        "slot-only replacement, either same-ID envelope alone, either waiver retirement alone, split authority-head transitions, sequential pack reveal or post-reveal amendment",
    ] {
        assert!(
            bearing
                .g3_relations
                .iter()
                .any(|relation| relation.contains(token)),
            "bearing population G3 omits {token}"
        );
    }
    assert_ordered_contract_tokens(
        "bearing population G4 lifecycle",
        bearing.g4_schedule,
        &[
            "GovernanceCommitted freeze",
            "CandidateFrozen sealing",
            "independent custodian joint realization",
            "atomic coupled RealizationCommitted authority-head commit",
            "joint RevealedForAdjudication access",
            "metrology audit",
            "adjudication",
            "Closed finalization",
        ],
    );
    for token in [
        "forbid post-cancel protected access, either-pack reveal, adjudication and authority-head advancement",
        "unless the exact coupled atomic stage had already committed",
        "complete frame/event/metrology/censoring/access/attempt/stopping ledgers",
    ] {
        assert!(
            bearing.g4_schedule.contains(token),
            "bearing population G4 omits {token}"
        );
    }
    assert_exact_lifecycle_event_membership(
        &draft,
        "i14-bearing-population-validation-max",
        &[
            "bearing.governance_committed",
            "bearing.candidate_frozen",
            "bearing.custodian_realized",
            "bearing.coupled_authority_committed",
            "bearing.revealed_for_adjudication",
            "bearing.metrology_audited",
            "bearing.event_audited",
            "bearing.coverage_checked",
            "bearing.reliability_adjudicated",
            "bearing.closed",
        ],
    );
}

#[test]
#[allow(clippy::too_many_lines)]
fn i14_theorem_ratchets_are_ambitious_formal_and_falsifiable() {
    let draft = i14_draft();
    let policy = authored_spec(&draft, THEOREM_POLICY);
    for token in [
        "pre-candidate successor",
        "canonical proposition-and-definition AST",
        "deterministic total AST-to-formal translation",
        "runtime-premise schemas",
        "complete transitive axiom closure",
        "{propext, Quot.sound, Classical.choice}",
        "sorryAx",
        "strength-matched nonvacuity",
        "AUTHORITY_CONTRADICTION",
        "GenuineCountermodel",
        "quarantines proof and runtime authority",
    ] {
        assert!(policy.contains(token), "theorem policy omits {token}");
    }
    for forbidden_shortcut in [
        "Bare cohomology equivalence",
        "component passivity",
        "naming a Dirac structure",
    ] {
        assert!(policy.contains(forbidden_shortcut));
    }

    let composition = claim(&draft, "i14-passive-causal-sheaf-composition-theorem");
    for token in [
        "compatibility sheaf",
        "balance cosheaf",
        "relative-hypercohomology",
        "orientation/dualizing cosheaf",
        "discrete Green-Stokes identity",
        "maximally isotropic lossless Dirac relation",
        "nondegenerate ambient split-signature power pairing",
        "Causal dissipative",
        "higher-overlap",
        "Zero unaccounted defect",
        "non-strict dominance",
        "robust strict passivity",
        "relative-boundary compatibility",
        "chain homotopies",
        "KYP/passive-realization bridge",
        "cover refinement",
        "hypercohomology or crosswalk-holonomy obstruction",
        "S(x(T))-S(x(0))",
    ] {
        assert!(
            composition.statement.contains(token)
                || composition.hypotheses.iter().any(|h| h.contains(token)),
            "composition theorem omits {token}"
        );
    }
    assert!(composition.activation.contains("pre-candidate successor"));
    assert!(composition.no_claim.contains("bare sheaf/cohomology"));
    assert!(composition.no_claim.contains("component passivity"));
    let composition_card = authored_spec(&draft, "i14-passive-composition-theorem-card");
    for token in [
        "higher-overlap cocycle/coherence",
        "relative-boundary compatibility",
        "chain-homotopy representative invariance",
        "maximally isotropic lossless Dirac relation",
        "nondegenerate ambient split-signature power pairing",
        "zero unaccounted defect",
        "integral d>=Delta_accounted",
        "integral d-Delta_accounted>=mu*norm_signal^2",
        "dependency-aware outward bounds",
        "KYP/passive-realization bridge",
        "cover-refinement comparison",
        "relative-hypercohomology",
    ] {
        assert!(
            composition_card.contains(token),
            "theorem card omits {token}"
        );
    }

    let theorem_ratchets = [
        (
            "i14-hypercohomology-obstruction-localization-theorem",
            "i14-hypercohomology-obstruction-card",
            &["relative-hypercohomology", "obstruction", "minimal"] as &[&str],
        ),
        (
            "i14-cover-refinement-naturality-theorem",
            "i14-cover-refinement-naturality-card",
            &[
                "refinement",
                "cofinal",
                "comparison natural transformations",
            ],
        ),
        (
            "i14-kyp-sheaf-passivity-bridge-theorem",
            "i14-kyp-sheaf-bridge-card",
            &["generalized-positive-real", "KYP", "storage"],
        ),
    ];
    for (claim_id, card_id, tokens) in theorem_ratchets {
        let theorem = claim(&draft, claim_id);
        let card = authored_spec(&draft, card_id);
        for token in tokens {
            assert!(
                theorem.statement.contains(token)
                    || theorem
                        .hypotheses
                        .iter()
                        .any(|hypothesis| hypothesis.contains(token))
                    || card.contains(token),
                "{claim_id}/{card_id} omits {token}"
            );
        }
        assert!(
            theorem.activation.contains("successor"),
            "{claim_id} activation"
        );
        assert!(
            theorem.kill.contains("AuthorityContradiction"),
            "{claim_id} contradiction policy"
        );
        assert!(
            card.contains("not a proposition AST or proof"),
            "{card_id} no-authority boundary"
        );
    }

    let descent = claim(&draft, "i14-certified-fidelity-descent-theorem");
    assert!(descent.activation.contains("pre-candidate"));
    assert!(descent.no_claim.contains("adjacent-rung agreement"));
    assert!(descent.kill.contains("countermodel"));

    let falsifier = claim(&draft, "i14-maximal-counterexample-search");
    assert_eq!(falsifier.polarity, ClaimPolarity::Refutation);
    assert!(
        falsifier
            .statement
            .contains("cardinality-proved finite rational microgrammar")
    );
    assert!(
        falsifier
            .hypotheses
            .iter()
            .any(|h| h.contains("rank/unrank"))
    );
    assert!(
        falsifier
            .hypotheses
            .iter()
            .any(|h| h.contains("GenuineCountermodel"))
    );
    assert!(
        falsifier
            .hypotheses
            .iter()
            .any(|h| h.contains("AuthorityContradiction"))
    );
    assert!(falsifier.kill.contains("exact unproved claim revision"));
    assert!(falsifier.kill.contains("proof/countermodel pair"));
    assert!(
        falsifier
            .no_claim
            .contains("version-1 prose has no exhaustive-search authority")
    );
    let grammar = authored_spec(&draft, "i14-theorem-falsifier-grammar");
    assert!(grammar.contains("VERSION-1 TARGET GRAMMAR"));
    assert!(grammar.contains("no exhaustive-search or theorem-survival authority"));
}

#[test]
#[allow(clippy::too_many_lines)]
fn i14_holdouts_are_exact_stage_local_and_do_not_launder_blind_authority() {
    let draft = i14_draft();
    let actual = draft
        .fixtures
        .iter()
        .filter(|fixture| fixture.partition == Partition::HeldOut)
        .map(|fixture| fixture.id)
        .collect::<BTreeSet<_>>();
    assert_eq!(
        actual,
        HOLDOUTS.iter().map(|(fixture, _, _)| *fixture).collect()
    );
    assert!(
        draft.fixtures.iter().all(|fixture| match fixture.source {
            FixtureSource::AuthoredSpec { spec } => {
                !spec.contains("HOLDOUT_KIND: GovernedStandards")
            }
            FixtureSource::External { .. } => true,
        }),
        "GovernedStandards is an authorized-input path, not an AuthoredSpec holdout"
    );

    for (fixture, expected_tier, expected_leaf) in HOLDOUTS {
        let consumers = draft
            .obligations
            .iter()
            .filter(|row| row.decks.contains(fixture))
            .collect::<Vec<_>>();
        assert_eq!(consumers.len(), 1, "{fixture} stage locality");
        assert_eq!(consumers[0].leaf, *expected_leaf, "{fixture} consumer");
        assert_eq!(consumers[0].tier, *expected_tier, "{fixture} tier");
        let spec = authored_spec(&draft, fixture);
        assert!(
            spec.contains("holdout") || spec.contains("HOLDOUT") || spec.contains("custodian"),
            "{fixture} does not name holdout governance"
        );
        match *expected_tier {
            CampaignTier::Core => assert!(
                spec.contains("core holdout range"),
                "{fixture} omits Core seed range"
            ),
            CampaignTier::Max if spec.contains("HOLDOUT_KIND: Governed") => {
                assert!(spec.contains("No public seed"));
                assert!(spec.contains("RealizationCommitted"));
                assert!(spec.contains("DischargeReceipt"));
                assert!(spec.contains("same-ID typed External discharge-envelope root"));
                assert!(spec.contains("AmendmentRecord"));
                assert!(spec.contains("authority head"));
                assert!(spec.contains("IntegrityFailed"));
            }
            CampaignTier::Max => assert!(
                spec.contains("maximal holdout range"),
                "{fixture} omits Max seed range"
            ),
            CampaignTier::Smoke => panic!("I14 has no smoke holdout"),
        }
    }

    let policy = authored_spec(&draft, POLICY);
    for token in [
        "EVIDENCE_AUTHORITY_CLASS={PublicReplayCore,PublicReplayMax,GovernedBlindSlot,GovernedPhysicalSlot,GovernedPopulationSlot,GovernedStandards}",
        "Partition and FixtureSource syntax alone carry no custody, discharge or promotion authority",
        "GOVERNED_EXTERNAL_BINDING",
        "typed DischargeReceipt",
        "same-ID typed External discharge-envelope root",
        "predecessor manifest digest",
        "domain-separated transaction-intent digest",
        "TRANSACTION_INTENT_V1",
        "org.frankensim.i14.governed-transaction-intent.v1",
        "canonical length-framed successor-intent projection",
        "the 25 exact ASCII bytes I14_TRANSACTION_INTENT_V1 followed by one byte 0x00",
        "U16LE(23)",
        "operation_kind U16LE where BlindPopulationRealization=1,PhysicalValidationRealization=2,StandardsAuthorityCommit=3,CoupledExternalRealization=4",
        "governance_stage U16LE where RealizationCommitted=1,StandardsAuthorityCommitted=2",
        "0x0006 expected_authority_head AuthorityHead",
        "0x0009 mutation_fence MutationFence",
        "0x0013 discharge_receipt_schema_digest Digest",
        "0x0014 authority_commit_receipt_schema_digest Digest",
        "0x0015 governance_protocol_schema_digest Digest",
        "MutationFence binding nonzero Digest(idempotency_key)||U64LE(attempt_id)||U64LE(capability_epoch)",
        "FRAME_BYTES(x)=U64LE(byte_len)||x",
        "ProtectedBindingList=U64LE(count)||concat(FRAME_BYTES(ProtectedBinding))",
        "FutureArtifact=FRAME_BYTES(Utf8(target_slot_id))||FRAME_BYTES(Utf8(artifact_role))||Digest(artifact_schema)||Utf8List(related_waiver_subjects)||digest_union",
        "FutureArtifactList=U64LE(count)||concat(FRAME_BYTES(FutureArtifact))",
        "Each governed-slot replacement's related-waiver set is exactly the sorted set of protected bindings naming its target slot/role/schema",
        "Every output, final-successor and AmendmentRecord slot id is pairwise distinct",
        "P is at most 16777216 bytes",
        "I14_ENVELOPE_DIGEST_PENDING_V1",
        "I14_AMENDMENT_RECORD_PENDING_V1",
        "Tag 0x0007 binds the already-existing predecessor-stage receipt",
        "P never serializes a newly created discharge-receipt digest, authority-commit-receipt digest, realized envelope digest, final successor digest, realized AmendmentRecord digest or new authority head",
        "TRANSACTION_OPERATION_MATRIX_V1",
        "BlindPopulationRealization requires governance_stage=RealizationCommitted, authority_scope='i14.blind-population-realization', at least one governed-slot replacement and at least one retired waiver",
        "PhysicalValidationRealization requires governance_stage=RealizationCommitted, authority_scope='i14.physical-validation-realization', at least one governed-slot replacement and at least one retired waiver",
        "StandardsAuthorityCommit requires governance_stage=StandardsAuthorityCommitted, authority_scope='i14.standards-authority-commit', exactly zero governed-slot replacements and at least one retired waiver",
        "CoupledExternalRealization requires governance_stage=RealizationCommitted, authority_scope='i14.coupled-external-realization', at least one governed-slot replacement and at least two retired waivers",
        "discharge-envelope count equals retired-waiver count",
        "No other operation/stage/scope/cardinality tuple is valid",
        "TRANSACTION_SCHEMA_REGISTRY_V1",
        "org.frankensim.i14.governed-schema-artifact.v1",
        "GovernanceSchemaSetRootV1=BLAKE3::derive_key('org.frankensim.i14.governance-schema-set.v1'",
        "The three exact case-sensitive roles are DischargeReceipt, AuthorityCommitReceipt and GovernanceProtocol",
        "U64LE(3)||FRAME_BYTES(Utf8('DischargeReceipt')||raw32(discharge_receipt_schema_digest))||FRAME_BYTES(Utf8('AuthorityCommitReceipt')||raw32(authority_commit_receipt_schema_digest))||FRAME_BYTES(Utf8('GovernanceProtocol')||raw32(governance_protocol_schema_digest))",
        "in exactly that role order",
        "The authenticated predecessor-stage receipt binds this exact root, each role/digest, the operation-matrix row, group identity and a closed role-addressed artifact-schema table",
        "transaction-time callers cannot substitute another schema or scope",
        "org.frankensim.i14.governed-artifact-schema.v1",
        "every ProtectedBinding/FutureArtifact artifact_schema must equal the table's exact ArtifactSchemaDigestV1 for that role",
        "Every fetched schema is at most 1048576 bytes, has nesting depth at most 64 and at most 4096 field/constraint records",
        "every DischargeReceipt is at most 1048576 bytes with nesting depth at most 64",
        "cap+1, depth+1, truncation, trailing bytes, role alias and schema swap are IntegrityFailed before mutation",
        "TRANSACTION_DECODER_DEPTH_V1",
        "The top-level schema or receipt value begins at nesting depth 0",
        "Depth 64 is admitted, any attempted entry to depth 65 is refused before allocating or decoding that child",
        "TRANSACTION_SCHEMA_ROLE_CONFORMANCE_V1",
        "the decoded AuthorityCommitReceipt schema must describe exactly the closed 361-byte AuthorityCommitReceiptBytesV1 layout with no optional or extension field",
        "The AuthorityCommitReceipt proves schema membership by binding transaction_intent",
        "A schema that omits, relaxes, renames, retypes or extends any required field/constraint is a schema mismatch",
        "org.frankensim.i14.discharge-receipt-set.v1",
        "receipt_count equals the retired-waiver count, is at most 4096",
        "org.frankensim.i14.realized-output-set.v1",
        "org.frankensim.i14.authority-head.v1",
        "RealizedArtifactBytesV1=FRAME_BYTES(Utf8(target_slot_id))||FRAME_BYTES(Utf8(artifact_role))||raw32(artifact_schema_digest)||Utf8List(related_waiver_subjects)||raw32(realized_artifact_digest)",
        "it is the exact realized image of every Pending governed-slot and discharge-envelope record and no other output",
        "The exact bytes of every content-addressed transaction, receipt and artifact schema must be available",
        "used to decode the corresponding receipt, governed-slot wrapper, discharge-envelope wrapper and realized output before mutation",
        "the 31 exact ASCII bytes I14_AUTHORITY_COMMIT_RECEIPT_V1 plus NUL",
        "exactly 361 bytes",
        "Every repeated receipt value must byte-equal its P/MutationFence value or the independently recomputed root, digest or head",
        "authority_capability_digest must byte-equal the authenticated capability that performs the commit",
        "compare-and-swaps the exact expected head/generation to the exact successor head/generation",
        "durably records the idempotency-key-to-intent-and-receipt mapping",
        "A crash before the single durable commit point leaves no effect",
        "a crash after it preserves the exact committed head, successor, waiver/output changes, idempotency mapping and receipt",
        "Exact-key byte-identical replay performs no second mutation or capability consumption",
        "To avoid a content-hash cycle",
        "two independently implemented encoders must reproduce published exact P, schema-set, receipt-set, output-set, successor-head and commit-receipt KATs",
        "per-field receipt mismatch",
        "crash-before-commit, crash-after-commit and partial-state mutation twins",
        "TRANSACTION_PARTIAL_STATE_RECOVERY_V1",
        "If recovery ever observes only a strict subset of the promised head/successor/waiver/output/idempotency/receipt state, that observation is itself durable corruption",
        "it must never be reported as a successful zero-effect refusal or silently completed from unauthenticated evidence",
        "Version 1 does not claim those independent encoders already exist",
        "this transaction path remains NoPromotionAuthority",
        "The AmendmentRecord binds predecessor and final successor digests",
        "a structurally frozen successor lacking the verified fs-vvreg transaction never proves discharge",
        "GOVERNED_BLIND_POPULATION_PATH",
        "GovernanceCommitted",
        "CandidateFrozen",
        "RealizationCommitted",
        "RevealedForAdjudication",
        "Closed",
        "independent custodian",
        "before protected outcome, label, aggregate or adjudication access",
        "inaccessible cases",
        "replace every governed AuthoredSpec commitment slot",
        "every waiver subject's same-ID typed External discharge-envelope root",
        "partial multi-pack discharge",
        "slot-only replacement",
        "waiver-only retirement is IntegrityFailed",
        "GOVERNED_PHYSICAL_VALIDATION_PATH",
        "protected calibration/model-input strata from an untouched validation stratum",
        "salted or equivalently hiding content/Merkle commitments",
        "hiding validation source-universe/frame commitment",
        "disjoint-membership commitment",
        "exact validation-selection algorithm",
        "pre-candidate secret-seed/VRF",
        "AuthorizedCalibration",
        "AsBuiltModelInstantiation",
        "complete read/output/derivation ledger",
        "receive only opaque commitment identities",
        "no validation bytes, membership witnesses, labels, aggregates, derived statistics, selection outputs or commitment-opening material",
        "independent contamination receipt names and audits every candidate-side principal and transitive capability",
        "separately addressable",
        "membership proofs to the pre-access commitments",
        "mutual disjoint-membership proof",
        "non-adaptive selection proof",
        "Candidate-side validation access or opening before freeze",
        "GOVERNED_STANDARDS_PATH",
        "GovernedStandards has no HeldOut AuthoredSpec commitment slot in manifest version 1",
        "AuthorizedConstruction",
        "CandidateFrozen occurs after that authorized construction",
        "StandardsAuthorityCommitted FrozenManifest::amend",
        "StandardsAdjudicated",
        "never blind, untouched-data, IID",
        "HOLDOUT_REALIZATION",
        "NoPromotionAuthority",
        "JOINT_EXTERNAL_TRANSACTION",
        "one atomic signed transaction",
        "partial retirement",
        "Sequential adaptation",
    ] {
        assert!(policy.contains(token), "campaign policy omits {token}");
    }
    let intent_fields = [
        ("0x0001", "initiative_id", "Utf8"),
        ("0x0002", "schema_identity", "Utf8"),
        ("0x0003", "operation_kind", "U16LE"),
        ("0x0004", "successor_version", "U64LE"),
        ("0x0005", "predecessor_manifest_digest", "Digest"),
        ("0x0006", "expected_authority_head", "AuthorityHead"),
        ("0x0007", "predecessor_stage_receipt_digest", "Digest"),
        ("0x0008", "candidate_freeze_commitment_digest", "Digest"),
        ("0x0009", "mutation_fence", "MutationFence"),
        ("0x000a", "governance_stage", "U16LE"),
        ("0x000b", "authority_scope", "Utf8"),
        ("0x000c", "coupled_transaction_group", "Utf8"),
        ("0x000d", "retired_waiver_subjects", "Utf8List"),
        ("0x000e", "protected_bindings", "ProtectedBindingList"),
        ("0x000f", "governed_slot_replacements", "FutureArtifactList"),
        ("0x0010", "discharge_envelope_slots", "FutureArtifactList"),
        ("0x0011", "final_successor_slot", "FutureDigest"),
        ("0x0012", "amendment_record_slot", "FutureDigest"),
        ("0x0013", "discharge_receipt_schema_digest", "Digest"),
        ("0x0014", "authority_commit_receipt_schema_digest", "Digest"),
        ("0x0015", "governance_protocol_schema_digest", "Digest"),
        ("0x0016", "access_ledger_prefix_digest", "Digest"),
        ("0x0017", "redaction_policy_digest", "Digest"),
    ];
    assert_eq!(intent_fields.len(), 23);
    let intent_header = b"I14_TRANSACTION_INTENT_V1\0";
    assert_eq!(intent_header.len(), 26);
    assert_eq!(intent_header.last(), Some(&0));
    assert!(!intent_header.ends_with(br"\0"));
    assert!(!policy.contains(r"I14_TRANSACTION_INTENT_V1\0"));
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
    let authority_receipt_header = b"I14_AUTHORITY_COMMIT_RECEIPT_V1\0";
    assert_eq!(authority_receipt_header.len(), 32);
    assert_eq!(authority_receipt_header.last(), Some(&0));
    assert_eq!(
        authority_receipt_header.len() + 3 * 8 + 9 * 32 + 2 * 8 + 1,
        361,
        "authority commit-receipt fixed-width arithmetic drifted"
    );
    let blind = authored_spec(&draft, "i14-mitigation-max-holdout");
    assert!(blind.contains("No public seed"));
    assert!(blind.contains("independent custodian"));
    assert!(blind.contains("GovernanceCommitted"));
    assert!(blind.contains("CandidateFrozen"));
    assert!(blind.contains("without protected holdout root, byte, label, aggregate"));
    assert!(blind.contains("public development fixtures"));
    assert!(blind.contains("frozen candidate-input permissions"));
    assert!(blind.contains("RealizationCommitted"));
    assert!(blind.contains("same-ID typed External discharge-envelope root"));
    assert!(blind.contains("complete authority transaction"));

    let physical = authored_spec(&draft, "i14-laboratory-validation-max-holdout");
    for token in [
        "disjoint calibration, as-built model-input and untouched validation strata",
        "salted or equivalently hiding calibration/model-input content commitments",
        "hiding validation source-universe/frame commitment",
        "disjoint-membership commitment",
        "exact validation-selection algorithm",
        "pre-candidate secret-seed/VRF",
        "AuthorizedCalibration",
        "AsBuiltModelInstantiation",
        "CandidateFrozen",
        "receive only opaque commitment identities",
        "no validation bytes, membership witnesses, labels, aggregates, derived statistics, selection outputs or commitment-opening material",
        "independent contamination receipt names and audits every candidate-side principal and transitive capability",
        "separately addressable calibration, model-input and validation roots",
        "membership proofs to the pre-access commitments",
        "mutual disjoint-membership proof",
        "non-adaptive selection proof",
        "before untouched validation reveal",
        "Candidate-side validation access or opening before freeze",
    ] {
        assert!(physical.contains(token), "physical slot omits {token}");
    }

    let class_specific_discharge: &[(&str, &[&str])] = &[
        (
            "i14-external-standards-edition-clause-pack",
            &[
                "GovernanceCommitted",
                "AuthorizedConstruction",
                "CandidateFrozen",
                "StandardsAuthorityCommitted FrozenManifest::amend",
                "no AuthoredSpec slot is replaced",
                "no blind, untouched, IID, crosswalk, conformance",
            ],
        ),
        (
            "i14-external-emc-laboratory-calibration-pack",
            &[
                "GovernanceCommitted, AuthorizedCalibration, AsBuiltModelInstantiation and CandidateFrozen",
                "i14-laboratory-validation-max-holdout",
                "joined typed External realized root",
                "i14-external-asbuilt-specimen-geometry-pack",
                "distinct receipts",
                "separately addressable calibration/model-input/untouched-validation roots",
                "immutable pre-access calibration/model-input content commitments",
                "validation source-universe root",
                "disjoint-membership commitment",
                "exact validation-selection algorithm",
                "pre-candidate secret-seed/VRF",
                "membership proofs to those commitments",
                "mutual disjoint-membership proof",
                "non-adaptive selection proof",
                "contamination receipt naming every audited candidate-side principal and transitive capability",
                "no forbidden validation payload/opening access",
                "candidate-side validation access or opening before freeze",
                "before untouched validation reveal",
            ],
        ),
        (
            "i14-external-asbuilt-specimen-geometry-pack",
            &[
                "GovernanceCommitted, AuthorizedCalibration, AsBuiltModelInstantiation and CandidateFrozen",
                "i14-laboratory-validation-max-holdout",
                "joined typed External realized root",
                "i14-external-emc-laboratory-calibration-pack",
                "distinct receipts",
                "separately addressable calibration/model-input/untouched-validation roots",
                "immutable pre-access calibration/model-input content commitments",
                "validation source-universe root",
                "disjoint-membership commitment",
                "exact validation-selection algorithm",
                "pre-candidate secret-seed/VRF",
                "membership proofs to those commitments",
                "mutual disjoint-membership proof",
                "non-adaptive selection proof",
                "contamination receipt naming every audited candidate-side principal and transitive capability",
                "no forbidden validation payload/opening access",
                "candidate-side validation access or opening before freeze",
                "before untouched validation reveal",
            ],
        ),
        (
            "i14-external-bearing-population-reliability-pack",
            &[
                "GovernanceCommitted and CandidateFrozen",
                "i14-bearing-population-max-holdout",
                "joined typed External realized root",
                "i14-external-bearing-population-metrology-pack",
                "distinct receipts",
                "before joint reveal",
            ],
        ),
        (
            "i14-external-bearing-population-metrology-pack",
            &[
                "GovernanceCommitted and CandidateFrozen",
                "i14-bearing-population-max-holdout",
                "joined typed External realized root",
                "i14-external-bearing-population-reliability-pack",
                "distinct receipts",
                "before joint reveal",
            ],
        ),
        (
            "i14-external-blind-mitigation-custody-pack",
            &[
                "GovernanceCommitted and CandidateFrozen",
                "i14-mitigation-max-holdout",
                "atomic RealizationCommitted FrozenManifest::amend transaction",
                "one-shot joint reveal",
                "envelope-only replacement",
            ],
        ),
        (
            "i14-external-emc-reliability-population-pack",
            &[
                "GovernanceCommitted and CandidateFrozen",
                "i14-emc-reliability-max-holdout",
                "typed External realized root",
                "before joint reveal",
                "envelope-only replacement",
            ],
        ),
    ];
    for (subject, tokens) in class_specific_discharge {
        let waiver = draft
            .waivers
            .iter()
            .find(|waiver| waiver.subject == *subject)
            .unwrap_or_else(|| panic!("missing waiver {subject}"));
        for token in *tokens {
            assert!(
                waiver.predicate.contains(token),
                "{subject} predicate omits {token}"
            );
        }
    }
    let standards_waiver = draft
        .waivers
        .iter()
        .find(|waiver| waiver.subject == "i14-external-standards-edition-clause-pack")
        .expect("standards waiver");
    assert!(
        standards_waiver
            .expiry
            .contains("before any ungoverned licensed-input access")
    );
    assert!(
        standards_waiver
            .expiry
            .contains("before i14-standards-crosswalk-max independent adjudication or promotion")
    );
    assert!(standards_waiver.expiry.contains(
        "AuthorizedConstruction access is permitted only through the frozen governed capability"
    ));
    for trigger in [
        "candidate",
        "builder/toolchain/checker/AcceptanceCard",
        "envelope/receipt",
        "authority-head",
    ] {
        assert!(
            standards_waiver.expiry.contains(trigger),
            "standards expiry omits review trigger {trigger}"
        );
    }
    for subject in [
        "i14-external-emc-laboratory-calibration-pack",
        "i14-external-asbuilt-specimen-geometry-pack",
    ] {
        let waiver = draft
            .waivers
            .iter()
            .find(|waiver| waiver.subject == subject)
            .unwrap_or_else(|| panic!("missing physical waiver {subject}"));
        assert!(waiver.expiry.contains("before any ungoverned"));
        assert!(
            waiver
                .expiry
                .contains("before i14-laboratory-validation-max untouched-validation reveal, independent adjudication or promotion")
        );
        assert!(waiver.expiry.contains(
            "AuthorizedCalibration and AsBuiltModelInstantiation access is permitted only through the frozen governed capability"
        ));
        for trigger in [
            "pre-access content commitment",
            "validation source-universe/frame commitment",
            "disjoint-membership commitment",
            "validation-selection algorithm",
            "seed/VRF or equivalent non-adaptive commitment",
            "membership/disjointness/non-adaptive proof scheme",
            "audited principal/transitive capability",
            "candidate/model/toolchain/checker/AcceptanceCard identity",
            "discharge envelope or receipt",
            "authority head",
            "custodian signature",
            "independent-adjudicator signature",
        ] {
            assert!(
                waiver.expiry.contains(trigger),
                "{subject} expiry omits review trigger {trigger}"
            );
        }
    }
    for (subject, reveal_scope) in [
        (
            "i14-external-bearing-population-metrology-pack",
            "joint RevealedForAdjudication",
        ),
        (
            "i14-external-bearing-population-reliability-pack",
            "joint RevealedForAdjudication",
        ),
        (
            "i14-external-blind-mitigation-custody-pack",
            "one-shot RevealedForAdjudication",
        ),
        (
            "i14-external-emc-reliability-population-pack",
            "one-shot RevealedForAdjudication",
        ),
    ] {
        let waiver = draft
            .waivers
            .iter()
            .find(|waiver| waiver.subject == subject)
            .unwrap_or_else(|| panic!("missing governed waiver {subject}"));
        for trigger in [
            reveal_scope,
            "acquisition/selection mechanism",
            "candidate/model/toolchain/checker/AcceptanceCard",
            "envelope/receipt",
            "authority-head",
            "custody change",
        ] {
            assert!(
                waiver.expiry.contains(trigger),
                "{subject} expiry omits review trigger {trigger}"
            );
        }
    }
}

#[test]
fn i14_canonical_identity_ignores_all_set_and_top_level_order() {
    let expected = i14_draft().freeze().expect("baseline freeze");
    let mut permuted = i14_draft();
    permuted.claims.reverse();
    permuted.fixtures.reverse();
    permuted.obligations.reverse();
    permuted.waivers.reverse();
    for row in &mut permuted.obligations {
        let mut claims = row.claims_covered.to_vec();
        claims.reverse();
        row.claims_covered = Box::leak(claims.into_boxed_slice());
        let mut unit_cases = row.unit_cases.to_vec();
        unit_cases.reverse();
        row.unit_cases = Box::leak(unit_cases.into_boxed_slice());
        let mut decks = row.decks.to_vec();
        decks.reverse();
        row.decks = Box::leak(decks.into_boxed_slice());
        let mut relations = row.g3_relations.to_vec();
        relations.reverse();
        row.g3_relations = Box::leak(relations.into_boxed_slice());
        let mut events = row.obs_events.to_vec();
        events.reverse();
        row.obs_events = Box::leak(events.into_boxed_slice());
    }
    let actual = permuted.freeze().expect("permuted freeze");
    assert_eq!(actual.digest(), expected.digest());
    assert_eq!(actual, expected);
}

#[test]
fn i14_mutations_refuse_or_change_exact_authored_identity() {
    let baseline = i14_draft().freeze().expect("baseline freeze").digest();

    let mut no_hypotheses = i14_draft();
    no_hypotheses.claims[0].hypotheses = &[];
    assert!(matches!(
        no_hypotheses.freeze(),
        Err(FreezeRefusal::BlankField {
            field: "claim.hypotheses",
            ..
        })
    ));

    let mut correlated = i14_draft();
    correlated
        .claims
        .iter_mut()
        .find(|claim| claim.id == "i14-fullwave-feec-stability-energy")
        .expect("fullwave claim")
        .oracle
        .independent = false;
    assert!(matches!(
        correlated.freeze(),
        Err(FreezeRefusal::ProductionOracleReuse { claim })
            if claim == "i14-fullwave-feec-stability-energy"
    ));

    let mut tolerance = i14_draft();
    tolerance
        .claims
        .iter_mut()
        .find(|claim| claim.id == "i14-mtl-rlgc-operator-admission")
        .expect("RLGC claim")
        .tolerance = ToleranceSemantics::Interval { lo: 0.0, hi: 2.0 };
    assert_ne!(
        tolerance.freeze().expect("tolerance freeze").digest(),
        baseline
    );

    let mut repartitioned = i14_draft();
    repartitioned
        .fixtures
        .iter_mut()
        .find(|fixture| fixture.id == "i14-source-probe-core-holdout")
        .expect("heldout")
        .partition = Partition::Development;
    assert_ne!(
        repartitioned.freeze().expect("partition freeze").digest(),
        baseline
    );

    let mut policy = i14_draft();
    policy
        .fixtures
        .iter_mut()
        .find(|fixture| fixture.id == POLICY)
        .expect("policy")
        .source = FixtureSource::AuthoredSpec {
        spec: "unauthorized weaker I14 campaign policy",
    };
    assert_ne!(policy.freeze().expect("policy freeze").digest(), baseline);

    let mut relation = i14_draft();
    relation
        .obligations
        .iter_mut()
        .find(|row| row.leaf == "i14-passive-composition-theorem-max")
        .expect("theorem leaf")
        .g3_relations = &["unauthorized weakened relation"];
    assert_ne!(
        relation.freeze().expect("relation freeze").digest(),
        baseline
    );
}

#[test]
#[allow(clippy::too_many_lines)]
fn i14_amendments_invalidate_targeted_consumers_and_global_policy() {
    let predecessor = i14_draft();
    let all = authority_ids(&predecessor);
    let frozen = predecessor.freeze().expect("baseline freeze");

    let mut version_only = i14_draft();
    version_only.version = 2;
    let (_, record) = frozen.amend(version_only).expect("version amendment");
    assert!(record.invalidated.is_empty());

    let mut claim_change = i14_draft();
    claim_change.version = 2;
    claim_change
        .claims
        .iter_mut()
        .find(|claim| claim.id == "i14-mtl-rlgc-operator-admission")
        .expect("RLGC claim")
        .statement = "successor RLGC claim with deliberately changed semantic authority";
    let (_, record) = frozen.amend(claim_change).expect("claim amendment");
    assert_eq!(
        record.invalidated,
        vec!["i14-mtl-rlgc-operator-admission", "i14-rlgc-operator-core",]
    );

    let mut heldout_change = i14_draft();
    heldout_change.version = 2;
    heldout_change
        .fixtures
        .iter_mut()
        .find(|fixture| fixture.id == "i14-source-probe-core-holdout")
        .expect("source holdout")
        .source = FixtureSource::AuthoredSpec {
        spec: "successor source/probe holdout corpus",
    };
    let (_, record) = frozen.amend(heldout_change).expect("holdout amendment");
    assert_eq!(
        record.invalidated,
        vec![
            "i14-source-probe-core",
            "i14-switching-source-probe-semantics",
        ]
    );

    let mut theorem_card_change = i14_draft();
    theorem_card_change.version = 2;
    theorem_card_change
        .fixtures
        .iter_mut()
        .find(|fixture| fixture.id == "i14-passive-composition-theorem-card")
        .expect("passive theorem card")
        .source = FixtureSource::AuthoredSpec {
        spec: "successor passive-composition theorem target",
    };
    let (_, record) = frozen
        .amend(theorem_card_change)
        .expect("theorem-card amendment");
    assert_eq!(
        record.invalidated.into_iter().collect::<BTreeSet<_>>(),
        BTreeSet::from([
            "i14-cover-refinement-naturality-max".to_string(),
            "i14-cover-refinement-naturality-theorem".to_string(),
            "i14-hypercohomology-obstruction-localization-theorem".to_string(),
            "i14-hypercohomology-obstruction-max".to_string(),
            "i14-maximal-counterexample-search".to_string(),
            "i14-maximal-falsifier-max".to_string(),
            "i14-passive-causal-sheaf-composition-theorem".to_string(),
            "i14-passive-composition-theorem-max".to_string(),
        ])
    );

    for (card_id, theorem_id, leaf_id) in [
        (
            "i14-hypercohomology-obstruction-card",
            "i14-hypercohomology-obstruction-localization-theorem",
            "i14-hypercohomology-obstruction-max",
        ),
        (
            "i14-cover-refinement-naturality-card",
            "i14-cover-refinement-naturality-theorem",
            "i14-cover-refinement-naturality-max",
        ),
        (
            "i14-kyp-sheaf-bridge-card",
            "i14-kyp-sheaf-passivity-bridge-theorem",
            "i14-kyp-sheaf-bridge-max",
        ),
    ] {
        let mut theorem_card_change = i14_draft();
        theorem_card_change.version = 2;
        theorem_card_change
            .fixtures
            .iter_mut()
            .find(|fixture| fixture.id == card_id)
            .unwrap_or_else(|| panic!("missing theorem card {card_id}"))
            .source = FixtureSource::AuthoredSpec {
            spec: "successor independently ratcheted theorem target",
        };
        let (_, record) = frozen
            .amend(theorem_card_change)
            .unwrap_or_else(|error| panic!("{card_id} amendment: {error}"));
        assert_eq!(
            record.invalidated.into_iter().collect::<BTreeSet<_>>(),
            BTreeSet::from([
                theorem_id.to_string(),
                leaf_id.to_string(),
                "i14-maximal-counterexample-search".to_string(),
                "i14-maximal-falsifier-max".to_string(),
            ]),
            "{card_id} amendment blast radius"
        );
    }

    let mut falsifier_change = i14_draft();
    falsifier_change.version = 2;
    falsifier_change
        .fixtures
        .iter_mut()
        .find(|fixture| fixture.id == "i14-theorem-falsifier-grammar")
        .expect("falsifier grammar")
        .source = FixtureSource::AuthoredSpec {
        spec: "successor executable falsifier grammar target",
    };
    let (_, record) = frozen.amend(falsifier_change).expect("falsifier amendment");
    assert_eq!(
        record.invalidated.into_iter().collect::<BTreeSet<_>>(),
        BTreeSet::from([
            "i14-maximal-counterexample-search".to_string(),
            "i14-maximal-falsifier-max".to_string(),
        ])
    );

    let waiver_blast_radii: &[(&str, &[&str])] = &[
        (
            "i14-external-standards-edition-clause-pack",
            &[
                "i14-governed-standards-crosswalk",
                "i14-standards-crosswalk-max",
            ],
        ),
        (
            "i14-external-emc-laboratory-calibration-pack",
            &[
                "i14-governed-laboratory-emc-validation",
                "i14-laboratory-validation-max",
            ],
        ),
        (
            "i14-external-asbuilt-specimen-geometry-pack",
            &[
                "i14-governed-laboratory-emc-validation",
                "i14-laboratory-validation-max",
            ],
        ),
        (
            "i14-external-bearing-population-reliability-pack",
            &[
                "i14-bearing-population-validation-max",
                "i14-production-bearing-population-reliability",
            ],
        ),
        (
            "i14-external-bearing-population-metrology-pack",
            &[
                "i14-bearing-population-validation-max",
                "i14-production-bearing-population-reliability",
            ],
        ),
        (
            "i14-external-emc-reliability-population-pack",
            &[
                "i14-emc-reliability-validation-max",
                "i14-governed-emc-reliability-validation",
            ],
        ),
        (
            "i14-external-blind-mitigation-custody-pack",
            &["i14-robust-mitigation-heldout", "i14-robust-mitigation-max"],
        ),
    ];
    for (subject, expected) in waiver_blast_radii {
        let mut changed = i14_draft();
        changed.version = 2;
        changed
            .waivers
            .iter_mut()
            .find(|waiver| waiver.subject == *subject)
            .unwrap_or_else(|| panic!("missing waiver {subject}"))
            .promotion_effect = "test-only changed scoped waiver authority";
        let (_, record) = frozen.amend(changed).expect("waiver amendment");
        assert_eq!(
            record.invalidated.into_iter().collect::<BTreeSet<_>>(),
            expected.iter().map(|id| (*id).to_string()).collect(),
            "{subject} amendment blast radius"
        );
    }

    let governed_subject = "i14-external-blind-mitigation-custody-pack";
    let mut syntactic_replacement = i14_draft();
    syntactic_replacement.version = 2;
    syntactic_replacement
        .waivers
        .retain(|waiver| waiver.subject != governed_subject);
    syntactic_replacement.fixtures.push(FixturePin {
        id: governed_subject,
        source: FixtureSource::External {
            digest_hex: "1111111111111111111111111111111111111111111111111111111111111111",
        },
        partition: Partition::HeldOut,
    });
    let (_, record) = frozen
        .amend(syntactic_replacement)
        .expect("same-id waiver-to-external-fixture syntax change");
    assert_eq!(
        record.invalidated.into_iter().collect::<BTreeSet<_>>(),
        BTreeSet::from([
            "i14-robust-mitigation-heldout".to_string(),
            "i14-robust-mitigation-max".to_string(),
        ])
    );

    let mut slot_only_replacement = i14_draft();
    slot_only_replacement.version = 2;
    slot_only_replacement
        .fixtures
        .iter_mut()
        .find(|fixture| fixture.id == "i14-mitigation-max-holdout")
        .expect("blind commitment slot")
        .source = FixtureSource::External {
        digest_hex: "3333333333333333333333333333333333333333333333333333333333333333",
    };
    let (_, record) = frozen
        .amend(slot_only_replacement)
        .expect("slot-only syntax change");
    assert_eq!(
        record.invalidated.into_iter().collect::<BTreeSet<_>>(),
        BTreeSet::from([
            "i14-robust-mitigation-heldout".to_string(),
            "i14-robust-mitigation-max".to_string(),
        ])
    );

    let standards_subject = "i14-external-standards-edition-clause-pack";
    let mut standards_envelope_replacement = i14_draft();
    standards_envelope_replacement.version = 2;
    standards_envelope_replacement
        .waivers
        .retain(|waiver| waiver.subject != standards_subject);
    // `Partition::HeldOut` is the only protected-fixture syntax in schema v2;
    // the campaign policy explicitly denies it custody, access, or holdout authority.
    standards_envelope_replacement.fixtures.push(FixturePin {
        id: standards_subject,
        source: FixtureSource::External {
            digest_hex: "2222222222222222222222222222222222222222222222222222222222222222",
        },
        partition: Partition::HeldOut,
    });
    let (_, record) = frozen
        .amend(standards_envelope_replacement)
        .expect("standards same-id waiver-to-external-fixture syntax change");
    assert_eq!(
        record.invalidated.into_iter().collect::<BTreeSet<_>>(),
        BTreeSet::from([
            "i14-governed-standards-crosswalk".to_string(),
            "i14-standards-crosswalk-max".to_string(),
        ])
    );

    let laboratory_subject = "i14-external-emc-laboratory-calibration-pack";
    let mut partial_laboratory_retirement = i14_draft();
    partial_laboratory_retirement.version = 2;
    partial_laboratory_retirement
        .waivers
        .retain(|waiver| waiver.subject != laboratory_subject);
    partial_laboratory_retirement.fixtures.push(FixturePin {
        id: laboratory_subject,
        source: FixtureSource::External {
            digest_hex: "4444444444444444444444444444444444444444444444444444444444444444",
        },
        partition: Partition::HeldOut,
    });
    partial_laboratory_retirement
        .fixtures
        .iter_mut()
        .find(|fixture| fixture.id == "i14-laboratory-validation-max-holdout")
        .expect("laboratory commitment slot")
        .source = FixtureSource::External {
        digest_hex: "5555555555555555555555555555555555555555555555555555555555555555",
    };
    let (_, record) = frozen
        .amend(partial_laboratory_retirement)
        .expect("partial laboratory retirement is structurally representable");
    assert_eq!(
        record.invalidated.into_iter().collect::<BTreeSet<_>>(),
        BTreeSet::from([
            "i14-governed-laboratory-emc-validation".to_string(),
            "i14-laboratory-validation-max".to_string(),
        ])
    );

    let bearing_subject = "i14-external-bearing-population-reliability-pack";
    let mut partial_bearing_retirement = i14_draft();
    partial_bearing_retirement.version = 2;
    partial_bearing_retirement
        .waivers
        .retain(|waiver| waiver.subject != bearing_subject);
    partial_bearing_retirement.fixtures.push(FixturePin {
        id: bearing_subject,
        source: FixtureSource::External {
            digest_hex: "6666666666666666666666666666666666666666666666666666666666666666",
        },
        partition: Partition::HeldOut,
    });
    partial_bearing_retirement
        .fixtures
        .iter_mut()
        .find(|fixture| fixture.id == "i14-bearing-population-max-holdout")
        .expect("bearing commitment slot")
        .source = FixtureSource::External {
        digest_hex: "7777777777777777777777777777777777777777777777777777777777777777",
    };
    let (_, record) = frozen
        .amend(partial_bearing_retirement)
        .expect("partial bearing retirement is structurally representable");
    assert_eq!(
        record.invalidated.into_iter().collect::<BTreeSet<_>>(),
        BTreeSet::from([
            "i14-bearing-population-validation-max".to_string(),
            "i14-production-bearing-population-reliability".to_string(),
        ])
    );
    let policy_draft = i14_draft();
    let campaign_policy = authored_spec(&policy_draft, POLICY);
    assert!(campaign_policy.contains("same-ID fixture alone"));
    assert!(campaign_policy.contains("slot replacement alone"));
    assert!(campaign_policy.contains("typed DischargeReceipt"));
    assert!(campaign_policy.contains("verified fs-vvreg transaction never proves discharge"));
    assert!(campaign_policy.contains("partial multi-pack discharge"));
    assert!(campaign_policy.contains("partial retirement"));
    assert!(campaign_policy.contains("waiver-only retirement is IntegrityFailed"));

    let mut convention_change = i14_draft();
    convention_change.version = 2;
    convention_change
        .fixtures
        .iter_mut()
        .find(|fixture| fixture.id == EM_CONVENTION_CARD)
        .expect("EM convention card")
        .source = FixtureSource::AuthoredSpec {
        spec: "successor EM convention card with deliberately changed sign authority",
    };
    let (_, record) = frozen
        .amend(convention_change)
        .expect("EM-convention amendment");
    let mut expected_convention_blast_radius = all.clone();
    for graph_only_authority in [
        "i14-harnessgraph-identity-connectivity",
        "i14-harnessgraph-core",
        "i14-synthetic-ap242-adapter-mechanics",
        "i14-ap242-adapter-core",
    ] {
        assert!(expected_convention_blast_radius.remove(graph_only_authority));
    }
    assert_eq!(
        record
            .invalidated
            .iter()
            .map(String::as_str)
            .collect::<BTreeSet<_>>(),
        expected_convention_blast_radius
    );
    assert_eq!(record.invalidated.len(), 52);

    let mut acceptance_change = i14_draft();
    acceptance_change.version = 2;
    acceptance_change
        .fixtures
        .iter_mut()
        .find(|fixture| fixture.id == ACCEPTANCE_POLICY)
        .expect("acceptance policy")
        .source = FixtureSource::AuthoredSpec {
        spec: "successor acceptance arithmetic with deliberately changed authority",
    };
    let (_, record) = frozen
        .amend(acceptance_change)
        .expect("acceptance-policy amendment");
    assert_eq!(
        record
            .invalidated
            .iter()
            .map(String::as_str)
            .collect::<BTreeSet<_>>(),
        all
    );
    assert_eq!(record.invalidated.len(), 56);

    let mut policy_change = i14_draft();
    policy_change.version = 2;
    policy_change
        .fixtures
        .iter_mut()
        .find(|fixture| fixture.id == POLICY)
        .expect("campaign policy")
        .source = FixtureSource::AuthoredSpec {
        spec: "I14 campaign policy successor with intentionally changed global authority",
    };
    let (_, record) = frozen.amend(policy_change).expect("policy amendment");
    assert_eq!(
        record
            .invalidated
            .iter()
            .map(String::as_str)
            .collect::<BTreeSet<_>>(),
        all
    );
    assert_eq!(record.invalidated.len(), 56);

    let mut title_change = i14_draft();
    title_change.version = 2;
    title_change.title = "successor global I14 campaign authority";
    let (_, record) = frozen.amend(title_change).expect("title amendment");
    assert_eq!(
        record
            .invalidated
            .iter()
            .map(String::as_str)
            .collect::<BTreeSet<_>>(),
        all
    );
}
