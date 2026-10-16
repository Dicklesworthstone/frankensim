//! G0/G3 fidelity-campaign extraction, fit, freshness, and ledger battery.

use std::collections::BTreeMap;
use std::time::Instant;

use fs_blake3::{ContentHash, hash_bytes};
use fs_conduction::{
    EMISSIVITY_DIMS, LinearizedSurfaceRadiation, SURFACE_EMISSIVITY_PROPERTY, SurfaceEmissivity,
};
use fs_convection::{CorrelationId, CorrelationInputs, ThermalDirection, evaluate};
use fs_evidence::ValidityDomain as EvidenceValidityDomain;
use fs_ladder::{
    ClosedInterval, FidelityGraph, FidelityNode, ModelCardRef, ModelId, QoiId, RegimeAxis,
    TransferRef,
};
use fs_ledger::{EdgeRole, Ledger};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, SelectionPolicy, UncertaintyModel,
};
use fs_plan::{
    CAMPAIGN_ARTIFACT_KIND, CampaignAuthority, CampaignError, CampaignGap, CampaignRun,
    EdgeProbeCampaign, FIDELITY_COST_ARTIFACT_KIND, FIDELITY_DISCREPANCY_ARTIFACT_KIND,
    FIDELITY_GRAPH_ARTIFACT_KIND, FreshnessReason, RunPartition, fit_fidelity_campaign,
    record_fidelity_campaign,
};

#[derive(Clone, Copy)]
struct EdgeModels {
    source: ModelId,
    target: ModelId,
    transfer: TransferRef,
}

fn id(label: &str) -> ContentHash {
    hash_bytes(label.as_bytes())
}

fn models(source: &str, target: &str) -> EdgeModels {
    EdgeModels {
        source: ModelId::new(id(source)),
        target: ModelId::new(id(target)),
        transfer: TransferRef::new(id(&format!("transfer:{source}->{target}"))),
    }
}

fn graph(name: &str, edges: &[EdgeModels]) -> FidelityGraph {
    let mut graph = FidelityGraph::new(name).expect("graph");
    for edge in edges {
        for (model, label) in [(edge.source, "source"), (edge.target, "target")] {
            if graph.node(model).is_none() {
                graph
                    .add_node(
                        FidelityNode::new(
                            model,
                            ModelCardRef::new(id(&format!("card:{model}"))),
                            format!("{label}-{}", &model.to_string()[..12]),
                        )
                        .expect("node"),
                    )
                    .expect("unique node");
            }
        }
    }
    graph
}

fn authority(edges: &[EdgeModels]) -> CampaignAuthority {
    let mut model_builds = BTreeMap::new();
    for edge in edges {
        model_builds.insert(edge.source, id(&format!("build:{}", edge.source)));
        model_builds.insert(edge.target, id(&format!("build:{}", edge.target)));
    }
    CampaignAuthority {
        corpus: id("cooling-corpus-level-a-execution-manifest-v1"),
        corpus_version: 1,
        machine_fingerprint: b"test-machine/aarch64/same-process".to_vec(),
        model_builds,
    }
}

fn regime(entries: &[(&str, f64, f64)]) -> BTreeMap<RegimeAxis, ClosedInterval> {
    entries
        .iter()
        .map(|(axis, lower, upper)| {
            (
                RegimeAxis::new(*axis).expect("axis"),
                ClosedInterval::new(*lower, *upper).expect("interval"),
            )
        })
        .collect()
}

fn params(entries: &[(&str, f64)]) -> BTreeMap<String, f64> {
    entries
        .iter()
        .map(|(name, value)| ((*name).to_string(), *value))
        .collect()
}

fn measured<T>(f: impl FnOnce() -> T) -> (T, f64) {
    let start = Instant::now();
    let value = f();
    (value, start.elapsed().as_secs_f64().max(1.0e-9))
}

fn run(
    label: &str,
    partition: RunPartition,
    params: BTreeMap<String, f64>,
    problem_size: f64,
    qois: (f64, f64, Option<f64>),
    costs: (f64, f64),
) -> CampaignRun {
    CampaignRun {
        run_id: id(&format!("paired-execution-receipt:{label}")),
        case_id: label.to_string(),
        partition,
        params,
        problem_size,
        source_qoi: qois.0,
        target_qoi: qois.1,
        reference_qoi: qois.2,
        source_cost_s: costs.0,
        target_cost_s: costs.1,
    }
}

fn emissivity_card() -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
            value: PropertyValue::Scalar {
                value: 0.8,
                dims: EMISSIVITY_DIMS,
            },
            validity: EvidenceValidityDomain::unconstrained().with("T", 250.0, 500.0),
            uncertainty: UncertaintyModel::HalfWidth {
                half_width: 0.01,
                confidence: 0.95,
            },
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: Vec::new(),
            provenance: Provenance {
                source: "f85xj.10.2 executing radiation campaign fixture".to_string(),
                license: "internal-test-use".to_string(),
                artifact: None,
            },
        })
        .expect("emissivity claim");
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "fixture-alloy".to_string(),
            phase: "solid".to_string(),
            process: "black-anodized".to_string(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .expect("material card")
}

fn radiation_campaign(edge: EdgeModels) -> EdgeProbeCampaign {
    let emissivity = SurfaceEmissivity::from_card(
        "radiator",
        &emissivity_card(),
        320.0,
        SelectionPolicy::SingleClaimOnly,
    )
    .expect("emissivity");
    let model = LinearizedSurfaceRadiation::new("radiator", emissivity, 320.0, 310.0, 25.0)
        .expect("linearized model");
    let temperatures = [300.0, 305.0, 315.0, 325.0, 335.0, 340.0];
    let runs = temperatures
        .into_iter()
        .enumerate()
        .map(|(index, temperature)| {
            let (source_point, source_cost) =
                measured(|| model.evaluate(temperature).expect("source radiation"));
            let (target_point, target_cost) =
                measured(|| model.evaluate(temperature).expect("target radiation"));
            run(
                &format!("radiation-temperature-{temperature:.0}K"),
                if index < 4 {
                    RunPartition::Fit
                } else {
                    RunPartition::HeldOut
                },
                params(&[
                    ("surface_temperature_k", temperature),
                    ("ambient_temperature_k", 310.0),
                    ("emissivity", 0.8),
                ]),
                1.0,
                (
                    source_point.linearized_outward_flux_w_m2,
                    target_point.nonlinear_outward_flux_w_m2,
                    Some(target_point.nonlinear_outward_flux_w_m2),
                ),
                (source_cost, target_cost),
            )
        })
        .collect();
    EdgeProbeCampaign {
        source: edge.source,
        target: edge.target,
        qoi: QoiId::new("outward-radiative-flux").expect("qoi"),
        qoi_unit: "W/m2".to_string(),
        transfer: edge.transfer,
        regime_bin: regime(&[
            ("surface_temperature_k", 295.0, 345.0),
            ("ambient_temperature_k", 310.0, 310.0),
            ("emissivity", 0.8, 0.8),
        ]),
        runs,
    }
}

fn correlation_campaign(
    edge: EdgeModels,
    source_id: CorrelationId,
    target_id: CorrelationId,
    cases: &[(f64, f64, Option<f64>)],
    length_ratio: Option<f64>,
) -> EdgeProbeCampaign {
    let runs = cases
        .iter()
        .enumerate()
        .map(|(index, (reynolds, prandtl, length))| {
            let inputs = CorrelationInputs::forced(*reynolds, *prandtl)
                .with_direction(ThermalDirection::CoolingFluid);
            let inputs = length.map_or(inputs, |ratio| inputs.with_length_ratio(ratio));
            let (source, source_cost) =
                measured(|| evaluate(source_id, inputs).expect("source correlation"));
            let (target, target_cost) =
                measured(|| evaluate(target_id, inputs).expect("target correlation"));
            let mut point = params(&[
                ("Re", *reynolds),
                ("Pr", *prandtl),
                ("Pe", reynolds * prandtl),
            ]);
            if let Some(ratio) = length {
                point.insert("L_over_Dh".to_string(), *ratio);
            }
            run(
                &format!("{}-vs-{}-{index}", source_id.name(), target_id.name()),
                if index < 4 {
                    RunPartition::Fit
                } else {
                    RunPartition::HeldOut
                },
                point,
                *reynolds,
                (source.evidence().value, target.evidence().value, None),
                (source_cost, target_cost),
            )
        })
        .collect();
    let mut axes = vec![
        ("Re", cases[0].0, cases[cases.len() - 1].0),
        ("Pr", 0.6, 10.0),
    ];
    if let Some(ratio) = length_ratio {
        axes.push(("L_over_Dh", ratio, ratio));
    }
    EdgeProbeCampaign {
        source: edge.source,
        target: edge.target,
        qoi: QoiId::new("nusselt-number").expect("qoi"),
        qoi_unit: "1".to_string(),
        transfer: edge.transfer,
        regime_bin: regime(&axes),
        runs,
    }
}

#[test]
fn full_campaign_logs_three_executing_edges_fits_predicates_and_graph_diff() {
    let radiation = models("radiation-linearized", "radiation-full-t4");
    let developing = models("duct-laminar-cwt", "duct-laminar-hausen");
    let turbulent = models("duct-dittus-boelter", "duct-gnielinski");
    let edges = [radiation, developing, turbulent];
    let graph = graph("electronics-cooling-fidelity", &edges);
    let before = graph.identity();
    let campaign = fit_fidelity_campaign(
        "cooling-adjacent-probes-2026-07",
        graph,
        authority(&edges),
        vec![
            radiation_campaign(radiation),
            correlation_campaign(
                developing,
                CorrelationId::CircularDuctLaminarCwt,
                CorrelationId::CircularDuctHausen,
                &[
                    (500.0, 0.7, Some(100.0)),
                    (800.0, 1.0, Some(100.0)),
                    (1_100.0, 2.0, Some(100.0)),
                    (1_400.0, 4.0, Some(100.0)),
                    (1_700.0, 7.0, Some(100.0)),
                    (2_000.0, 10.0, Some(100.0)),
                ],
                Some(100.0),
            ),
            correlation_campaign(
                turbulent,
                CorrelationId::DittusBoelter,
                CorrelationId::Gnielinski,
                &[
                    (10_000.0, 0.7, Some(100.0)),
                    (20_000.0, 1.0, Some(100.0)),
                    (40_000.0, 2.0, Some(100.0)),
                    (60_000.0, 4.0, Some(100.0)),
                    (80_000.0, 7.0, Some(100.0)),
                    (100_000.0, 10.0, Some(100.0)),
                ],
                Some(100.0),
            ),
        ],
        vec![
            CampaignGap {
                source: ModelId::new(id("correlation-nu")),
                target: ModelId::new(id("steady-rans")),
                qoi: QoiId::new("junction-temperature").expect("qoi"),
                reason: "no retained paired correlation/RANS execution set".to_string(),
            },
            CampaignGap {
                source: ModelId::new(id("churchill-chu-vertical-plate")),
                target: ModelId::new(id("thermal-lbm-rayleigh-benard")),
                qoi: QoiId::new("nusselt-number").expect("qoi"),
                reason: "current correlation is a vertical external plate while thermal LBM is a horizontal periodic Rayleigh-Benard slab; no shared-validity geometry exists"
                    .to_string(),
            },
            CampaignGap {
                source: ModelId::new(id("homogenized-pcb")),
                target: ModelId::new(id("resolved-pcb")),
                qoi: QoiId::new("junction-temperature").expect("qoi"),
                reason: "resolved PCB rung is not implemented".to_string(),
            },
        ],
    )
    .expect("fit complete cooling campaign");

    assert_eq!(campaign.graph_before, before);
    assert_ne!(campaign.graph_after(), before);
    assert_eq!(campaign.edges.len(), 3);
    assert_eq!(campaign.graph.edges().len(), 3);
    assert_eq!(campaign.gaps.len(), 3);
    assert_eq!(
        campaign
            .edges
            .iter()
            .filter(|edge| edge.informativeness_supported)
            .count(),
        1,
        "only full T^4 has an independent held-out reference"
    );
    assert!(
        campaign
            .edges
            .iter()
            .filter(|edge| !edge.informativeness_supported)
            .all(|edge| edge.edge.informativeness().predicates().is_unknown())
    );

    let ledger = Ledger::open(":memory:").expect("campaign ledger");
    let receipt =
        record_fidelity_campaign(&ledger, &campaign, 1_000, 2_000).expect("atomic retention");
    assert_eq!(receipt.edge_artifacts.len(), 3);
    assert_eq!(
        ledger
            .artifact_info(&receipt.graph_before.hash)
            .expect("graph-before info")
            .expect("graph-before exists")
            .kind,
        FIDELITY_GRAPH_ARTIFACT_KIND
    );
    assert_eq!(
        ledger
            .artifact_info(&receipt.graph_after.hash)
            .expect("graph-after info")
            .expect("graph-after exists")
            .kind,
        FIDELITY_GRAPH_ARTIFACT_KIND
    );
    assert_eq!(
        ledger
            .artifact_info(&receipt.campaign.hash)
            .expect("campaign info")
            .expect("campaign exists")
            .kind,
        CAMPAIGN_ARTIFACT_KIND
    );
    for (fitted, (discrepancy, cost)) in campaign.edges.iter().zip(&receipt.edge_artifacts) {
        assert_eq!(discrepancy.hash, fitted.discrepancy_artifact_id());
        assert_eq!(cost.hash, fitted.cost_artifact_id());
        assert_eq!(
            ledger
                .artifact_info(&discrepancy.hash)
                .expect("discrepancy info")
                .expect("discrepancy exists")
                .kind,
            FIDELITY_DISCREPANCY_ARTIFACT_KIND
        );
        assert_eq!(
            ledger
                .artifact_info(&cost.hash)
                .expect("cost info")
                .expect("cost exists")
                .kind,
            FIDELITY_COST_ARTIFACT_KIND
        );
        assert!(
            ledger
                .edge_exists(receipt.op, &discrepancy.hash, EdgeRole::Out)
                .expect("discrepancy lineage")
        );
    }
    assert!(ledger.lint().expect("lint").is_clean());
    println!(
        "{{\"suite\":\"fs-plan-fidelity-campaign\",\"edges\":3,\"gaps\":3,\
         \"informativeness_supported\":1,\"graph_before\":\"{}\",\
         \"graph_after\":\"{}\",\"ledger_op\":{},\"verdict\":\"pass\"}}",
        campaign.graph_before,
        campaign.graph_after(),
        receipt.op
    );
}

fn synthetic_campaign(runs: Vec<CampaignRun>) -> (CampaignAuthority, EdgeProbeCampaign) {
    let edge = models("synthetic-cheap", "synthetic-reference");
    (
        authority(&[edge]),
        EdgeProbeCampaign {
            source: edge.source,
            target: edge.target,
            qoi: QoiId::new("synthetic-qoi").expect("qoi"),
            qoi_unit: "K".to_string(),
            transfer: edge.transfer,
            regime_bin: regime(&[("x", 1.0, 6.0)]),
            runs,
        },
    )
}

fn synthetic_runs() -> Vec<CampaignRun> {
    (1..=6)
        .map(|index| {
            let size = f64::from(index);
            run(
                &format!("synthetic-{index}"),
                if index <= 4 {
                    RunPartition::Fit
                } else {
                    RunPartition::HeldOut
                },
                params(&[("x", size)]),
                size,
                (9.0, 10.0, Some(10.0)),
                (0.01 * size, 0.10 * size),
            )
        })
        .collect()
}

#[test]
fn synthetic_truth_recovery_is_partitioned_and_order_deterministic() {
    let edge = models("synthetic-cheap", "synthetic-reference");
    let runs = synthetic_runs();
    let (authority_a, input_a) = synthetic_campaign(runs.clone());
    let first = fit_fidelity_campaign(
        "synthetic-recovery",
        graph("synthetic-graph", &[edge]),
        authority_a,
        vec![input_a],
        Vec::new(),
    )
    .expect("first fit");

    let mut reversed = runs;
    reversed.reverse();
    let (authority_b, input_b) = synthetic_campaign(reversed);
    let second = fit_fidelity_campaign(
        "synthetic-recovery",
        graph("synthetic-graph", &[edge]),
        authority_b,
        vec![input_b],
        Vec::new(),
    )
    .expect("reordered fit");

    let fit = first.edges.first().expect("one fitted edge");
    assert!((fit.discrepancy_band.mean_observed_rel - 0.1).abs() < 1.0e-15);
    assert!((fit.discrepancy_band.max_observed_rel - 0.1).abs() < 1.0e-15);
    assert!((fit.held_out_discrepancy_coverage - 1.0).abs() <= f64::EPSILON);
    assert!((fit.source_cost_calibration - 1.0).abs() <= f64::EPSILON);
    assert!((fit.target_cost_calibration - 1.0).abs() <= f64::EPSILON);
    assert!(fit.informativeness_supported);
    assert_eq!(first.graph_after(), second.graph_after());
    assert_eq!(
        fit.discrepancy_artifact_id(),
        second
            .edges
            .first()
            .expect("one reordered fitted edge")
            .discrepancy_artifact_id()
    );
    assert_eq!(
        fit.cost_artifact_id(),
        second
            .edges
            .first()
            .expect("one reordered fitted edge")
            .cost_artifact_id()
    );
}

#[test]
fn extraction_refuses_missing_holdout_and_out_of_bin_runs() {
    let edge = models("synthetic-cheap", "synthetic-reference");
    let mut no_holdout = synthetic_runs();
    for run in &mut no_holdout {
        run.partition = RunPartition::Fit;
    }
    let (authority, input) = synthetic_campaign(no_holdout);
    let error = fit_fidelity_campaign(
        "no-holdout",
        graph("synthetic-graph", &[edge]),
        authority,
        vec![input],
        Vec::new(),
    )
    .expect_err("missing holdout must refuse");
    assert!(matches!(
        error,
        CampaignError::Invalid {
            field: "held-out partition",
            ..
        }
    ));

    let mut out_of_bin = synthetic_runs();
    out_of_bin
        .first_mut()
        .expect("synthetic runs")
        .params
        .insert("x".to_string(), 7.0);
    let (authority, input) = synthetic_campaign(out_of_bin);
    let error = fit_fidelity_campaign(
        "out-of-bin",
        graph("synthetic-graph", &[edge]),
        authority,
        vec![input],
        Vec::new(),
    )
    .expect_err("regime escape must refuse");
    assert!(matches!(
        error,
        CampaignError::Invalid {
            field: "run regime",
            ..
        }
    ));
}

#[test]
fn freshness_flags_corpus_machine_and_kernel_drift_without_mutating_fit() {
    let edge = models("synthetic-cheap", "synthetic-reference");
    let (authority, input) = synthetic_campaign(synthetic_runs());
    let campaign = fit_fidelity_campaign(
        "freshness",
        graph("synthetic-graph", &[edge]),
        authority.clone(),
        vec![input],
        Vec::new(),
    )
    .expect("fit");
    assert!(!campaign.assess_freshness(&authority).is_stale());

    let mut current = authority;
    current.corpus = id("new-corpus");
    current.corpus_version += 1;
    current.machine_fingerprint = b"different-machine".to_vec();
    current
        .model_builds
        .insert(edge.target, id("new-target-build"));
    let freshness = campaign.assess_freshness(&current);
    assert!(freshness.is_stale());
    assert!(
        freshness
            .reasons
            .contains(&FreshnessReason::CorpusIdentityChanged)
    );
    assert!(freshness.reasons.contains(&FreshnessReason::MachineChanged));
    assert!(freshness.reasons.iter().any(|reason| matches!(
        reason,
        FreshnessReason::CorpusVersionChanged {
            fitted: 1,
            current: 2
        }
    )));
    assert!(freshness.reasons.iter().any(|reason| matches!(
        reason,
        FreshnessReason::ModelBuildChanged { model } if *model == edge.target
    )));
}
