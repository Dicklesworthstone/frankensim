//! Empirical population of contextual [`fs_ladder::FidelityGraph`] edges.
//!
//! The graph crate owns declarations and selection. This L6 module owns the
//! orchestration boundary that turns exact paired executions into retained
//! discrepancy/cost artifacts and a replayable graph diff. A populated edge is
//! still not an informativeness claim: the target ordering remains `UNKNOWN`
//! unless every held-out pair carries an independent reference and the target
//! is no worse on every one (and strictly better on at least one).

use std::collections::{BTreeMap, BTreeSet};
use std::error::Error;
use std::fmt;

use fs_blake3::{ContentHash, hash_bytes};
use fs_evidence::{DiscrepancyBand, DiscrepancyModel, FidelityPair, FitError};
use fs_ladder::{
    ClosedInterval, ContextClause, ContextPredicateSet, CostModelRef, DiscrepancyModelRef,
    FidelityEdge, FidelityGraph, FidelityGraphId, GraphError, Informativeness, ModelId, QoiId,
    QoiSelector, RegimeAxis, TransferRef, ValidityDomain,
};
use fs_ledger::{EdgeRole, FiveExplicits, Ledger, LedgerError, OpOutcome, PutReceipt};

use crate::{CostModel, CostObservation, CostPrediction, CostRefusal, MIN_OBS};

/// Canonical campaign artifact schema.
pub const CAMPAIGN_SCHEMA_VERSION: u16 = 1;
/// Retained exact paired-run discrepancy fit.
pub const FIDELITY_DISCREPANCY_ARTIFACT_KIND: &str = "fs-plan-fidelity-discrepancy-model-v1";
/// Retained exact paired-run source/target cost fit.
pub const FIDELITY_COST_ARTIFACT_KIND: &str = "fs-plan-fidelity-cost-model-v1";
/// Canonical [`FidelityGraph`] bytes before and after population.
pub const FIDELITY_GRAPH_ARTIFACT_KIND: &str = "fs-ladder-fidelity-graph-v1";
/// Retained campaign summary and graph diff.
pub const CAMPAIGN_ARTIFACT_KIND: &str = "fs-plan-fidelity-campaign-v1";

/// Maximum paired executions in one edge fit.
pub const MAX_CAMPAIGN_RUNS_PER_EDGE: usize = 4_096;
/// Maximum edges fitted by one atomic campaign.
pub const MAX_CAMPAIGN_EDGES: usize = 256;
/// Maximum explicit acquisition gaps in one campaign.
pub const MAX_CAMPAIGN_GAPS: usize = 1_024;
/// Maximum model/build identities in one campaign authority.
pub const MAX_CAMPAIGN_MODELS: usize = 512;
/// Maximum named parameters or regime axes in one paired execution.
pub const MAX_CAMPAIGN_DIMENSIONS: usize = 64;
/// Maximum visible-ASCII bytes in a campaign/case/unit/reason identity.
pub const MAX_CAMPAIGN_NAME_BYTES: usize = 256;
/// Maximum opaque machine-fingerprint bytes.
pub const MAX_MACHINE_FINGERPRINT_BYTES: usize = 4_096;

/// Partition role of one exact paired execution.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum RunPartition {
    /// Used to fit discrepancy and source/target cost models.
    Fit,
    /// Withheld from all fits and used only for audits.
    HeldOut,
}

/// One exact adjacent-fidelity probe execution.
#[derive(Debug, Clone, PartialEq)]
pub struct CampaignRun {
    /// Content identity of the complete execution receipt.
    pub run_id: ContentHash,
    /// Stable corpus case identity.
    pub case_id: String,
    /// Fit or held-out role, fixed before fitting.
    pub partition: RunPartition,
    /// Exact shared parameter point.
    pub params: BTreeMap<String, f64>,
    /// Positive scalar feature consumed by [`CostModel`].
    pub problem_size: f64,
    /// Source/cheap-model QoI.
    pub source_qoi: f64,
    /// Target/finer-model QoI.
    pub target_qoi: f64,
    /// Independent reference QoI, when one exists.
    pub reference_qoi: Option<f64>,
    /// Measured source wall seconds.
    pub source_cost_s: f64,
    /// Measured target wall seconds.
    pub target_cost_s: f64,
}

/// Exact campaign-wide freshness authority.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CampaignAuthority {
    /// Content identity of the complete corpus manifest.
    pub corpus: ContentHash,
    /// Monotone corpus schema/content generation.
    pub corpus_version: u64,
    /// Machine fingerprint under which cost observations were measured.
    pub machine_fingerprint: Vec<u8>,
    /// Exact model implementation/build identities.
    pub model_builds: BTreeMap<ModelId, ContentHash>,
}

/// One contextual candidate edge and its paired executions.
#[derive(Debug, Clone, PartialEq)]
pub struct EdgeProbeCampaign {
    /// Less-informative candidate endpoint.
    pub source: ModelId,
    /// More-informative candidate endpoint.
    pub target: ModelId,
    /// Quantity compared by every run.
    pub qoi: QoiId,
    /// Exact visible unit identity for the QoI.
    pub qoi_unit: String,
    /// Transfer implementation/configuration used between endpoints.
    pub transfer: TransferRef,
    /// One explicit regime bin. All runs must lie inside it.
    pub regime_bin: BTreeMap<RegimeAxis, ClosedInterval>,
    /// Exact paired executions.
    pub runs: Vec<CampaignRun>,
}

/// A candidate comparison that remains an explicit acquisition gap.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CampaignGap {
    /// Candidate source endpoint.
    pub source: ModelId,
    /// Candidate target endpoint.
    pub target: ModelId,
    /// QoI whose evidence is missing.
    pub qoi: QoiId,
    /// Actionable reason the edge was not populated.
    pub reason: String,
}

/// One fitted edge plus exact artifact payloads.
#[derive(Debug, Clone)]
pub struct FittedCampaignEdge {
    /// Graph edge added by this fit.
    pub edge: FidelityEdge,
    /// Observed discrepancy statistics over fit pairs.
    pub discrepancy_band: DiscrepancyBand,
    /// Whether held-out discrepancy stayed within the fit maximum.
    pub held_out_discrepancy_coverage: f64,
    /// Source cost-band held-out coverage.
    pub source_cost_calibration: f64,
    /// Target cost-band held-out coverage.
    pub target_cost_calibration: f64,
    /// Whether independent held-out references established target ordering.
    pub informativeness_supported: bool,
    /// Canonical retained discrepancy artifact bytes.
    discrepancy_artifact: Vec<u8>,
    /// Canonical retained two-endpoint cost artifact bytes.
    cost_artifact: Vec<u8>,
}

impl FittedCampaignEdge {
    /// Exact discrepancy artifact identity referenced by the graph edge.
    #[must_use]
    pub fn discrepancy_artifact_id(&self) -> ContentHash {
        hash_bytes(&self.discrepancy_artifact)
    }

    /// Exact cost artifact identity referenced by the graph edge.
    #[must_use]
    pub fn cost_artifact_id(&self) -> ContentHash {
        hash_bytes(&self.cost_artifact)
    }
}

/// Pure fit result, ready for one atomic ledger publication.
#[derive(Debug, Clone)]
pub struct FittedCampaign {
    /// Bounded campaign identity.
    pub name: String,
    /// Freshness authority captured by the fit.
    pub authority: CampaignAuthority,
    /// Graph identity before population.
    pub graph_before: FidelityGraphId,
    /// Populated graph.
    pub graph: FidelityGraph,
    /// Fitted edges in deterministic edge-id order.
    pub edges: Vec<FittedCampaignEdge>,
    /// Unpopulated candidates retained explicitly.
    pub gaps: Vec<CampaignGap>,
    /// Exact canonical graph bytes before population.
    graph_before_bytes: Vec<u8>,
    /// Exact canonical campaign summary bytes.
    campaign_artifact: Vec<u8>,
}

impl FittedCampaign {
    /// Graph identity after population.
    #[must_use]
    pub fn graph_after(&self) -> FidelityGraphId {
        self.graph.identity()
    }

    /// Assess exact corpus, machine, and model-build freshness.
    #[must_use]
    pub fn assess_freshness(&self, current: &CampaignAuthority) -> CampaignFreshness {
        let mut reasons = Vec::new();
        if current.corpus != self.authority.corpus {
            reasons.push(FreshnessReason::CorpusIdentityChanged);
        }
        if current.corpus_version != self.authority.corpus_version {
            reasons.push(FreshnessReason::CorpusVersionChanged {
                fitted: self.authority.corpus_version,
                current: current.corpus_version,
            });
        }
        if current.machine_fingerprint != self.authority.machine_fingerprint {
            reasons.push(FreshnessReason::MachineChanged);
        }
        for (model, fitted) in &self.authority.model_builds {
            match current.model_builds.get(model) {
                Some(now) if now == fitted => {}
                Some(_) => reasons.push(FreshnessReason::ModelBuildChanged { model: *model }),
                None => reasons.push(FreshnessReason::ModelBuildMissing { model: *model }),
            }
        }
        CampaignFreshness { reasons }
    }
}

/// Exact reason a retained campaign must be treated as stale.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FreshnessReason {
    /// The complete corpus manifest identity changed.
    CorpusIdentityChanged,
    /// The corpus generation changed.
    CorpusVersionChanged {
        /// Version captured by the fit.
        fitted: u64,
        /// Version supplied by the current authority.
        current: u64,
    },
    /// Cost observations came from a different machine.
    MachineChanged,
    /// A model implementation/configuration changed.
    ModelBuildChanged {
        /// Changed model.
        model: ModelId,
    },
    /// The current authority omitted a fitted model.
    ModelBuildMissing {
        /// Missing model.
        model: ModelId,
    },
}

/// Freshness verdict with every exact cause retained.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CampaignFreshness {
    /// Empty means exact authority agreement.
    pub reasons: Vec<FreshnessReason>,
}

impl CampaignFreshness {
    /// Whether any authority component changed.
    #[must_use]
    pub fn is_stale(&self) -> bool {
        !self.reasons.is_empty()
    }
}

/// Atomic ledger publication receipt.
#[derive(Debug, Clone)]
pub struct CampaignLedgerReceipt {
    /// Finished ledger operation.
    pub op: i64,
    /// Graph artifact before the campaign.
    pub graph_before: PutReceipt,
    /// Graph artifact after the campaign.
    pub graph_after: PutReceipt,
    /// Per-edge `(discrepancy, cost)` artifact receipts.
    pub edge_artifacts: Vec<(PutReceipt, PutReceipt)>,
    /// Campaign summary artifact.
    pub campaign: PutReceipt,
}

/// Structured refusal from campaign extraction, fit, graph, or persistence.
#[derive(Debug)]
pub enum CampaignError {
    /// A bounded/canonical input invariant failed.
    Invalid {
        /// Offending field.
        field: &'static str,
        /// Actionable problem.
        problem: String,
    },
    /// Discrepancy fit refused.
    Discrepancy(FitError),
    /// Cost fit or audit refused.
    Cost(CostRefusal),
    /// Graph construction refused.
    Graph(GraphError),
    /// Ledger storage refused.
    Ledger(Box<LedgerError>),
    /// A primary persistence failure was followed by rollback failure.
    LedgerCleanup {
        /// Primary failure.
        primary: String,
        /// Cleanup failure.
        rollback: Box<LedgerError>,
    },
}

impl fmt::Display for CampaignError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid { field, problem } => {
                write!(f, "fidelity campaign refused `{field}`: {problem}")
            }
            Self::Discrepancy(error) => error.fmt(f),
            Self::Cost(error) => error.fmt(f),
            Self::Graph(error) => error.fmt(f),
            Self::Ledger(error) => error.fmt(f),
            Self::LedgerCleanup { primary, rollback } => write!(
                f,
                "fidelity campaign ledger write failed ({primary}); rollback also failed ({rollback})"
            ),
        }
    }
}

impl Error for CampaignError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::Discrepancy(error) => Some(error),
            Self::Cost(error) => Some(error),
            Self::Graph(error) => Some(error),
            Self::Ledger(error) => Some(error.as_ref()),
            Self::LedgerCleanup { rollback, .. } => Some(rollback.as_ref()),
            Self::Invalid { .. } => None,
        }
    }
}

impl From<FitError> for CampaignError {
    fn from(value: FitError) -> Self {
        Self::Discrepancy(value)
    }
}

impl From<CostRefusal> for CampaignError {
    fn from(value: CostRefusal) -> Self {
        Self::Cost(value)
    }
}

impl From<GraphError> for CampaignError {
    fn from(value: GraphError) -> Self {
        Self::Graph(value)
    }
}

impl From<LedgerError> for CampaignError {
    fn from(value: LedgerError) -> Self {
        Self::Ledger(Box::new(value))
    }
}

/// Fit paired executions and add their exact artifact references to a graph.
///
/// # Errors
///
/// Refuses malformed authority, missing graph nodes/build identities, mixed
/// parameter schemas, unbalanced partitions, invalid costs/QoIs, runs outside
/// the declared regime bin, duplicate run identities, or any lower-layer fit
/// or graph error.
pub fn fit_fidelity_campaign(
    name: impl Into<String>,
    mut graph: FidelityGraph,
    authority: CampaignAuthority,
    edge_campaigns: Vec<EdgeProbeCampaign>,
    mut gaps: Vec<CampaignGap>,
) -> Result<FittedCampaign, CampaignError> {
    let name = name.into();
    validate_identity("campaign name", &name)?;
    validate_authority(&authority)?;
    if edge_campaigns.is_empty() || edge_campaigns.len() > MAX_CAMPAIGN_EDGES {
        return Err(invalid(
            "edge campaigns",
            format!(
                "need 1..={MAX_CAMPAIGN_EDGES} edge campaigns, got {}",
                edge_campaigns.len()
            ),
        ));
    }
    if gaps.len() > MAX_CAMPAIGN_GAPS {
        return Err(invalid(
            "campaign gaps",
            format!("{} exceeds bounded maximum {MAX_CAMPAIGN_GAPS}", gaps.len()),
        ));
    }
    validate_gaps(&gaps)?;

    let graph_before = graph.identity();
    let graph_before_bytes = graph.canonical_bytes();
    let mut fitted = Vec::with_capacity(edge_campaigns.len());
    for campaign in edge_campaigns {
        fitted.push(fit_edge(&graph, &authority, campaign)?);
    }
    fitted.sort_by_key(|edge| edge.edge.id());
    for edge in &fitted {
        graph.add_edge(edge.edge.clone())?;
    }
    gaps.sort_by(|left, right| {
        left.source
            .cmp(&right.source)
            .then(left.target.cmp(&right.target))
            .then(left.qoi.cmp(&right.qoi))
            .then(left.reason.cmp(&right.reason))
    });
    let campaign_artifact = campaign_bytes(&name, &authority, graph_before, &graph, &fitted, &gaps);
    Ok(FittedCampaign {
        name,
        authority,
        graph_before,
        graph,
        edges: fitted,
        gaps,
        graph_before_bytes,
        campaign_artifact,
    })
}

struct RunPartitions<'a> {
    fit: Vec<&'a CampaignRun>,
    held_out: Vec<&'a CampaignRun>,
}

struct CalibratedCosts {
    source: CostModel,
    target: CostModel,
    source_calibration: f64,
    target_calibration: f64,
}

struct EdgeContext {
    validity: ValidityDomain,
    informativeness: Informativeness,
    informativeness_supported: bool,
}

fn validate_edge_campaign(
    graph: &FidelityGraph,
    authority: &CampaignAuthority,
    campaign: &EdgeProbeCampaign,
) -> Result<(), CampaignError> {
    if campaign.source == campaign.target {
        return Err(invalid("edge endpoints", "source and target are equal"));
    }
    for (field, model) in [
        ("source model", campaign.source),
        ("target model", campaign.target),
    ] {
        if graph.node(model).is_none() {
            return Err(invalid(
                field,
                format!("model {model} is absent from graph"),
            ));
        }
        if !authority.model_builds.contains_key(&model) {
            return Err(invalid(
                "model builds",
                format!("model {model} has no exact build identity"),
            ));
        }
    }
    validate_identity("QoI unit", &campaign.qoi_unit)?;
    if campaign.regime_bin.is_empty() || campaign.regime_bin.len() > MAX_CAMPAIGN_DIMENSIONS {
        return Err(invalid(
            "regime bin",
            format!(
                "need 1..={MAX_CAMPAIGN_DIMENSIONS} exact regime axes, got {}",
                campaign.regime_bin.len()
            ),
        ));
    }
    if campaign.runs.len() > MAX_CAMPAIGN_RUNS_PER_EDGE {
        return Err(invalid(
            "campaign runs",
            format!(
                "{} exceeds bounded maximum {MAX_CAMPAIGN_RUNS_PER_EDGE}",
                campaign.runs.len()
            ),
        ));
    }
    let mut seen = BTreeSet::new();
    for run in &campaign.runs {
        if !seen.insert(run.run_id) {
            return Err(invalid(
                "run id",
                format!("duplicate paired execution {}", run.run_id),
            ));
        }
        validate_run(run, &campaign.regime_bin)?;
    }
    Ok(())
}

fn partition_runs(campaign: &EdgeProbeCampaign) -> Result<RunPartitions<'_>, CampaignError> {
    let fit: Vec<_> = campaign
        .runs
        .iter()
        .filter(|run| run.partition == RunPartition::Fit)
        .collect();
    let held_out: Vec<_> = campaign
        .runs
        .iter()
        .filter(|run| run.partition == RunPartition::HeldOut)
        .collect();
    if fit.len() < MIN_OBS {
        return Err(invalid(
            "fit partition",
            format!(
                "{} paired runs cannot fit both cost models; need at least {MIN_OBS}",
                fit.len()
            ),
        ));
    }
    if held_out.is_empty() {
        return Err(invalid(
            "held-out partition",
            "at least one run must be withheld from all fits",
        ));
    }
    Ok(RunPartitions { fit, held_out })
}

fn fit_discrepancy(
    partitions: &RunPartitions<'_>,
) -> Result<(DiscrepancyBand, f64), CampaignError> {
    let pairs: Vec<_> = partitions
        .fit
        .iter()
        .map(|run| FidelityPair {
            params: run.params.clone(),
            lo_fi: run.source_qoi,
            hi_fi: run.target_qoi,
        })
        .collect();
    let discrepancy = DiscrepancyModel::fit(&pairs)?;
    let first_pair = pairs
        .first()
        .ok_or_else(|| invalid("fit partition", "discrepancy fit had no paired runs"))?;
    let discrepancy_band = discrepancy
        .query(&first_pair.params)
        .map_err(FitError::QueryOutOfDomain)?;
    let held_out_coverage = discrepancy_coverage(discrepancy_band, &partitions.held_out)?;
    Ok((discrepancy_band, held_out_coverage))
}

fn fit_cost_models(partitions: &RunPartitions<'_>) -> Result<CalibratedCosts, CampaignError> {
    let source_fit: Vec<_> = partitions
        .fit
        .iter()
        .map(|run| CostObservation {
            size: run.problem_size,
            cost_s: run.source_cost_s,
        })
        .collect();
    let target_fit: Vec<_> = partitions
        .fit
        .iter()
        .map(|run| CostObservation {
            size: run.problem_size,
            cost_s: run.target_cost_s,
        })
        .collect();
    let source_held_out: Vec<_> = partitions
        .held_out
        .iter()
        .map(|run| CostObservation {
            size: run.problem_size,
            cost_s: run.source_cost_s,
        })
        .collect();
    let target_held_out: Vec<_> = partitions
        .held_out
        .iter()
        .map(|run| CostObservation {
            size: run.problem_size,
            cost_s: run.target_cost_s,
        })
        .collect();
    let source_cost = CostModel::fit(&source_fit)?;
    let target_cost = CostModel::fit(&target_fit)?;
    let source_cost_calibration = source_cost.calibration(&source_held_out)?;
    let target_cost_calibration = target_cost.calibration(&target_held_out)?;
    Ok(CalibratedCosts {
        source: source_cost,
        target: target_cost,
        source_calibration: source_cost_calibration,
        target_calibration: target_cost_calibration,
    })
}

fn edge_context(
    campaign: &EdgeProbeCampaign,
    held_out: &[&CampaignRun],
) -> Result<EdgeContext, CampaignError> {
    let clause = ContextClause::new(
        QoiSelector::Exact(campaign.qoi.clone()),
        campaign
            .regime_bin
            .iter()
            .map(|(axis, interval)| (axis.clone(), *interval)),
    )?;
    let predicates = ContextPredicateSet::new([clause])?;
    let validity = ValidityDomain::new(predicates.clone());
    let informativeness_supported = held_out_supports_informativeness(held_out);
    let informativeness = if informativeness_supported {
        Informativeness::new(predicates)
    } else {
        Informativeness::unknown()
    };
    Ok(EdgeContext {
        validity,
        informativeness,
        informativeness_supported,
    })
}

fn fit_edge(
    graph: &FidelityGraph,
    authority: &CampaignAuthority,
    mut campaign: EdgeProbeCampaign,
) -> Result<FittedCampaignEdge, CampaignError> {
    campaign.runs.sort_by_key(|run| run.run_id);
    validate_edge_campaign(graph, authority, &campaign)?;
    let partitions = partition_runs(&campaign)?;
    let (discrepancy_band, held_out_discrepancy_coverage) = fit_discrepancy(&partitions)?;
    let costs = fit_cost_models(&partitions)?;
    let context = edge_context(&campaign, &partitions.held_out)?;

    let discrepancy_artifact = discrepancy_bytes(
        authority,
        &campaign,
        discrepancy_band,
        held_out_discrepancy_coverage,
    );
    let cost_artifact = cost_bytes(
        authority,
        &campaign,
        &costs.source,
        &costs.target,
        costs.source_calibration,
        costs.target_calibration,
    )?;
    let edge = FidelityEdge::new(
        campaign.source,
        campaign.target,
        CostModelRef::new(hash_bytes(&cost_artifact)),
        DiscrepancyModelRef::new(hash_bytes(&discrepancy_artifact)),
        campaign.transfer,
        context.validity,
        context.informativeness,
    )?;
    Ok(FittedCampaignEdge {
        edge,
        discrepancy_band,
        held_out_discrepancy_coverage,
        source_cost_calibration: costs.source_calibration,
        target_cost_calibration: costs.target_calibration,
        informativeness_supported: context.informativeness_supported,
        discrepancy_artifact,
        cost_artifact,
    })
}

fn validate_authority(authority: &CampaignAuthority) -> Result<(), CampaignError> {
    if authority.corpus_version == 0 {
        return Err(invalid("corpus version", "zero is not a published version"));
    }
    if authority.machine_fingerprint.is_empty()
        || authority.machine_fingerprint.len() > MAX_MACHINE_FINGERPRINT_BYTES
    {
        return Err(invalid(
            "machine fingerprint",
            format!(
                "need 1..={MAX_MACHINE_FINGERPRINT_BYTES} bytes, got {}",
                authority.machine_fingerprint.len()
            ),
        ));
    }
    if authority.model_builds.is_empty() || authority.model_builds.len() > MAX_CAMPAIGN_MODELS {
        return Err(invalid(
            "model builds",
            format!(
                "need 1..={MAX_CAMPAIGN_MODELS} exact build identities, got {}",
                authority.model_builds.len()
            ),
        ));
    }
    Ok(())
}

fn validate_gaps(gaps: &[CampaignGap]) -> Result<(), CampaignError> {
    for gap in gaps {
        if gap.source == gap.target {
            return Err(invalid("gap endpoints", "source and target are equal"));
        }
        validate_text("gap reason", &gap.reason)?;
    }
    Ok(())
}

fn validate_run(
    run: &CampaignRun,
    regime_bin: &BTreeMap<RegimeAxis, ClosedInterval>,
) -> Result<(), CampaignError> {
    validate_identity("case id", &run.case_id)?;
    if run.params.is_empty() || run.params.len() > MAX_CAMPAIGN_DIMENSIONS {
        return Err(invalid(
            "run parameters",
            format!(
                "need 1..={MAX_CAMPAIGN_DIMENSIONS} named parameters, got {}",
                run.params.len()
            ),
        ));
    }
    for (name, value) in &run.params {
        validate_identity("run parameter name", name)?;
        if !value.is_finite() {
            return Err(invalid("run parameter", format!("{name} is not finite")));
        }
    }
    for (axis, interval) in regime_bin {
        let value = run.params.get(axis.as_str()).ok_or_else(|| {
            invalid(
                "run parameters",
                format!("case {} omitted regime axis {axis}", run.case_id),
            )
        })?;
        if *value < interval.lower() || *value > interval.upper() {
            return Err(invalid(
                "run regime",
                format!(
                    "case {} has {axis}={value}, outside [{}, {}]",
                    run.case_id,
                    interval.lower(),
                    interval.upper()
                ),
            ));
        }
    }
    for (field, value) in [
        ("problem size", run.problem_size),
        ("source QoI", run.source_qoi),
        ("target QoI", run.target_qoi),
        ("source cost", run.source_cost_s),
        ("target cost", run.target_cost_s),
    ] {
        if !value.is_finite()
            || matches!(field, "problem size" | "source cost" | "target cost") && value <= 0.0
        {
            return Err(invalid(
                field,
                format!("{value} is outside its finite domain"),
            ));
        }
    }
    if run.reference_qoi.is_some_and(|value| !value.is_finite()) {
        return Err(invalid("reference QoI", "value is not finite"));
    }
    Ok(())
}

fn discrepancy_coverage(
    band: DiscrepancyBand,
    held_out: &[&CampaignRun],
) -> Result<f64, CampaignError> {
    let mut covered = 0usize;
    for run in held_out {
        let relative =
            (run.target_qoi - run.source_qoi).abs() / run.target_qoi.abs().max(f64::MIN_POSITIVE);
        if !relative.is_finite() {
            return Err(invalid(
                "held-out discrepancy",
                format!("case {} produced a non-finite relative error", run.case_id),
            ));
        }
        if relative <= band.max_observed_rel {
            covered += 1;
        }
    }
    Ok(covered as f64 / held_out.len() as f64)
}

fn held_out_supports_informativeness(held_out: &[&CampaignRun]) -> bool {
    let mut strictly_better = false;
    for run in held_out {
        let Some(reference) = run.reference_qoi else {
            return false;
        };
        let source_error = (run.source_qoi - reference).abs();
        let target_error = (run.target_qoi - reference).abs();
        if target_error > source_error {
            return false;
        }
        strictly_better |= target_error < source_error;
    }
    strictly_better
}

fn discrepancy_bytes(
    authority: &CampaignAuthority,
    campaign: &EdgeProbeCampaign,
    band: DiscrepancyBand,
    held_out_coverage: f64,
) -> Vec<u8> {
    let mut out = artifact_header(b"FSDISC01");
    put_authority(&mut out, authority);
    put_edge_header(&mut out, campaign);
    put_f64(&mut out, band.mean_observed_rel);
    put_f64(&mut out, band.max_observed_rel);
    put_f64(&mut out, held_out_coverage);
    put_runs(&mut out, &campaign.runs);
    out
}

fn cost_bytes(
    authority: &CampaignAuthority,
    campaign: &EdgeProbeCampaign,
    source: &CostModel,
    target: &CostModel,
    source_calibration: f64,
    target_calibration: f64,
) -> Result<Vec<u8>, CampaignError> {
    let mut out = artifact_header(b"FSCOST01");
    put_authority(&mut out, authority);
    put_edge_header(&mut out, campaign);
    put_f64(&mut out, source_calibration);
    put_f64(&mut out, target_calibration);
    let mut sizes: Vec<_> = campaign.runs.iter().map(|run| run.problem_size).collect();
    sizes.sort_by(f64::total_cmp);
    sizes.dedup_by(|left, right| left.to_bits() == right.to_bits());
    put_u32(&mut out, sizes.len() as u32);
    for size in sizes {
        put_f64(&mut out, size);
        put_prediction(&mut out, source.predict(size)?);
        put_prediction(&mut out, target.predict(size)?);
    }
    put_runs(&mut out, &campaign.runs);
    Ok(out)
}

fn campaign_bytes(
    name: &str,
    authority: &CampaignAuthority,
    before: FidelityGraphId,
    after: &FidelityGraph,
    edges: &[FittedCampaignEdge],
    gaps: &[CampaignGap],
) -> Vec<u8> {
    let mut out = artifact_header(b"FSCAMP01");
    put_str(&mut out, name);
    put_authority(&mut out, authority);
    put_hash(&mut out, before.hash());
    put_hash(&mut out, after.identity().hash());
    put_u32(&mut out, edges.len() as u32);
    for fitted in edges {
        put_hash(&mut out, fitted.edge.id().hash());
        put_hash(&mut out, fitted.cost_artifact_id());
        put_hash(&mut out, fitted.discrepancy_artifact_id());
        put_u8(&mut out, u8::from(fitted.informativeness_supported));
        put_f64(&mut out, fitted.held_out_discrepancy_coverage);
        put_f64(&mut out, fitted.source_cost_calibration);
        put_f64(&mut out, fitted.target_cost_calibration);
    }
    put_u32(&mut out, gaps.len() as u32);
    for gap in gaps {
        put_hash(&mut out, gap.source.hash());
        put_hash(&mut out, gap.target.hash());
        put_str(&mut out, gap.qoi.as_str());
        put_str(&mut out, &gap.reason);
    }
    out
}

/// Atomically retain graph-before, all fits, graph-after, and campaign summary.
///
/// The finished deterministic operation links the graph-before artifact as an
/// input and every newly retained artifact as an output. Caller-owned ledger
/// transactions are refused so a successful receipt always means this method
/// committed the complete publication.
///
/// # Errors
///
/// Returns [`CampaignError::Invalid`] when a caller transaction is already
/// open, or preserves the exact ledger/rollback failure.
pub fn record_fidelity_campaign(
    ledger: &Ledger,
    campaign: &FittedCampaign,
    t_start_ns: i64,
    t_end_ns: i64,
) -> Result<CampaignLedgerReceipt, CampaignError> {
    if ledger.in_transaction() {
        return Err(invalid(
            "ledger transaction",
            "caller-owned transaction is already open",
        ));
    }
    if t_end_ns < t_start_ns {
        return Err(invalid(
            "ledger time",
            "operation end precedes operation start",
        ));
    }
    ledger.begin()?;
    let result = record_inside_transaction(ledger, campaign, t_start_ns, t_end_ns);
    match result {
        Ok(receipt) => match ledger.commit() {
            Ok(()) => Ok(receipt),
            Err(primary) => match ledger.rollback() {
                Ok(()) => Err(CampaignError::Ledger(Box::new(primary))),
                Err(rollback) => Err(CampaignError::LedgerCleanup {
                    primary: primary.to_string(),
                    rollback: Box::new(rollback),
                }),
            },
        },
        Err(primary) => match ledger.rollback() {
            Ok(()) => Err(primary),
            Err(rollback) => Err(CampaignError::LedgerCleanup {
                primary: primary.to_string(),
                rollback: Box::new(rollback),
            }),
        },
    }
}

fn record_inside_transaction(
    ledger: &Ledger,
    campaign: &FittedCampaign,
    t_start_ns: i64,
    t_end_ns: i64,
) -> Result<CampaignLedgerReceipt, CampaignError> {
    let graph_before = ledger.put_artifact(
        FIDELITY_GRAPH_ARTIFACT_KIND,
        &campaign.graph_before_bytes,
        None,
    )?;
    let mut edge_artifacts = Vec::with_capacity(campaign.edges.len());
    for edge in &campaign.edges {
        let discrepancy = ledger.put_artifact(
            FIDELITY_DISCREPANCY_ARTIFACT_KIND,
            &edge.discrepancy_artifact,
            None,
        )?;
        let cost = ledger.put_artifact(FIDELITY_COST_ARTIFACT_KIND, &edge.cost_artifact, None)?;
        edge_artifacts.push((discrepancy, cost));
    }
    let graph_after_bytes = campaign.graph.canonical_bytes();
    let graph_after =
        ledger.put_artifact(FIDELITY_GRAPH_ARTIFACT_KIND, &graph_after_bytes, None)?;
    let campaign_receipt =
        ledger.put_artifact(CAMPAIGN_ARTIFACT_KIND, &campaign.campaign_artifact, None)?;

    let ir = format!(
        "{{\"op\":\"fit-fidelity-campaign\",\"schema\":{},\"campaign\":\"{}\",\
         \"qoi_units_bound_in_edge_artifacts\":true,\"edge_count\":{},\"gap_count\":{}}}",
        CAMPAIGN_SCHEMA_VERSION,
        campaign.name,
        campaign.edges.len(),
        campaign.gaps.len()
    );
    let versions = format!(
        "{{\"fs-plan\":\"{}\",\"campaign_schema\":{}}}",
        crate::VERSION,
        CAMPAIGN_SCHEMA_VERSION
    );
    let budget = format!(
        "{{\"max_edges\":{MAX_CAMPAIGN_EDGES},\"max_runs_per_edge\":{MAX_CAMPAIGN_RUNS_PER_EDGE}}}"
    );
    let capability =
        "{\"ledger_write\":true,\"fit\":\"bounded-pure\",\"informativeness\":\"held-out-only\"}";
    let seed = hash_bytes(&campaign.campaign_artifact);
    let explicits = FiveExplicits {
        seed: seed.as_bytes(),
        versions: &versions,
        budget: &budget,
        capability,
    };
    let op = ledger.begin_op(Some(campaign.name.as_bytes()), &ir, &explicits, t_start_ns)?;
    ledger.link(op, &graph_before.hash, EdgeRole::In)?;
    for (discrepancy, cost) in &edge_artifacts {
        ledger.link(op, &discrepancy.hash, EdgeRole::Out)?;
        ledger.link(op, &cost.hash, EdgeRole::Out)?;
    }
    ledger.link(op, &graph_after.hash, EdgeRole::Out)?;
    ledger.link(op, &campaign_receipt.hash, EdgeRole::Out)?;
    ledger.finish_op(op, OpOutcome::Ok, None, t_end_ns)?;
    Ok(CampaignLedgerReceipt {
        op,
        graph_before,
        graph_after,
        edge_artifacts,
        campaign: campaign_receipt,
    })
}

fn artifact_header(magic: &[u8]) -> Vec<u8> {
    let mut out = Vec::new();
    out.extend_from_slice(magic);
    put_u16(&mut out, CAMPAIGN_SCHEMA_VERSION);
    out
}

fn put_authority(out: &mut Vec<u8>, authority: &CampaignAuthority) {
    put_hash(out, authority.corpus);
    put_u64(out, authority.corpus_version);
    put_bytes(out, &authority.machine_fingerprint);
    put_u32(out, authority.model_builds.len() as u32);
    for (model, build) in &authority.model_builds {
        put_hash(out, model.hash());
        put_hash(out, *build);
    }
}

fn put_edge_header(out: &mut Vec<u8>, campaign: &EdgeProbeCampaign) {
    put_hash(out, campaign.source.hash());
    put_hash(out, campaign.target.hash());
    put_str(out, campaign.qoi.as_str());
    put_str(out, &campaign.qoi_unit);
    put_hash(out, campaign.transfer.hash());
    put_u32(out, campaign.regime_bin.len() as u32);
    for (axis, interval) in &campaign.regime_bin {
        put_str(out, axis.as_str());
        put_f64(out, interval.lower());
        put_f64(out, interval.upper());
    }
}

fn put_runs(out: &mut Vec<u8>, runs: &[CampaignRun]) {
    put_u32(out, runs.len() as u32);
    for run in runs {
        put_hash(out, run.run_id);
        put_str(out, &run.case_id);
        put_u8(
            out,
            match run.partition {
                RunPartition::Fit => 0,
                RunPartition::HeldOut => 1,
            },
        );
        put_u32(out, run.params.len() as u32);
        for (name, value) in &run.params {
            put_str(out, name);
            put_f64(out, *value);
        }
        put_f64(out, run.problem_size);
        put_f64(out, run.source_qoi);
        put_f64(out, run.target_qoi);
        match run.reference_qoi {
            None => put_u8(out, 0),
            Some(value) => {
                put_u8(out, 1);
                put_f64(out, value);
            }
        }
        put_f64(out, run.source_cost_s);
        put_f64(out, run.target_cost_s);
    }
}

fn put_prediction(out: &mut Vec<u8>, prediction: CostPrediction) {
    put_f64(out, prediction.p10);
    put_f64(out, prediction.p50);
    put_f64(out, prediction.p90);
    put_u64(out, prediction.n_obs as u64);
    put_u8(out, u8::from(prediction.extrapolated));
}

fn put_hash(out: &mut Vec<u8>, hash: ContentHash) {
    out.extend_from_slice(hash.as_bytes());
}

fn put_bytes(out: &mut Vec<u8>, value: &[u8]) {
    put_u32(out, value.len() as u32);
    out.extend_from_slice(value);
}

fn put_str(out: &mut Vec<u8>, value: &str) {
    put_bytes(out, value.as_bytes());
}

fn put_f64(out: &mut Vec<u8>, value: f64) {
    let canonical = if value == 0.0 { 0.0 } else { value };
    put_u64(out, canonical.to_bits());
}

fn put_u8(out: &mut Vec<u8>, value: u8) {
    out.push(value);
}

fn put_u16(out: &mut Vec<u8>, value: u16) {
    out.extend_from_slice(&value.to_le_bytes());
}

fn put_u32(out: &mut Vec<u8>, value: u32) {
    out.extend_from_slice(&value.to_le_bytes());
}

fn put_u64(out: &mut Vec<u8>, value: u64) {
    out.extend_from_slice(&value.to_le_bytes());
}

fn validate_identity(field: &'static str, value: &str) -> Result<(), CampaignError> {
    if value.is_empty()
        || value.len() > MAX_CAMPAIGN_NAME_BYTES
        || !value
            .bytes()
            .all(|byte| byte.is_ascii_alphanumeric() || b"._-:/[]()*^".contains(&byte))
    {
        return Err(invalid(
            field,
            format!("must be 1..={MAX_CAMPAIGN_NAME_BYTES} visible canonical ASCII bytes"),
        ));
    }
    Ok(())
}

fn validate_text(field: &'static str, value: &str) -> Result<(), CampaignError> {
    if value.trim().is_empty()
        || value.len() > MAX_CAMPAIGN_NAME_BYTES
        || value.chars().any(char::is_control)
    {
        return Err(invalid(
            field,
            format!("must be nonblank, control-free, and at most {MAX_CAMPAIGN_NAME_BYTES} bytes"),
        ));
    }
    Ok(())
}

fn invalid(field: &'static str, problem: impl Into<String>) -> CampaignError {
    CampaignError::Invalid {
        field,
        problem: problem.into(),
    }
}
