//! fs-report — automatic lab notebooks + semantic design diffs. Layer: L6.
//!
//! Reproducibility should be a SIDE EFFECT of running a study, not a virtue you
//! remember to practice. A [`LabNotebook`] is the automatic lab notebook: every
//! study emits a deterministic, human-readable report — provenance, prose,
//! Qty-labelled metrics (units on every value, P10), AND THE EXACT IR TO
//! REPRODUCE IT ([`LabNotebook::repro_ir`]). Because the render is deterministic
//! it is CONTENT-ADDRESSED ([`LabNotebook::content_hash`]), so replaying the IR
//! and re-rendering yields the same hash — the reproducibility loop closes by
//! construction.
//!
//! [`semantic_diff`] is the other half: a diff between two designs that is a
//! GEOMETRIC attribution ("lip curvature −18%, wall thinned 0.4 mm"), ranked by
//! significance — not a file diff.
//!
//! [`decision_headline_markdown`] projects an already-admitted
//! [`fs_session::DecisionAssessment`] into a compact reviewer-facing headline.
//! It does not recompute compliance, uncertainty, attribution, or action value.
//! [`project_decision_gate_markdown`] explains whether the declared project
//! context may use an already-computed tri-state verdict.
//! [`regime_no_claims_markdown`] projects final operating-envelope demotions
//! into the report's no-claim section without changing their evidence state.
//! [`retain_regime_demotions_in_package`] preserves those same receipts in an
//! evidence package as explicitly unbounded Estimated declarations.
//! [`project_regime_audit_outputs`] couples both projections so a product path
//! cannot silently retain only one side of the same audit.
//!
//! [`nominal_vs_as_built_markdown`] renders the design-build-measure-update
//! comparison: per-QoI nominal-versus-as-built shifts with the geometry
//! uncertainty term attributed, whether each shift is resolved by the pose
//! measurement at all, and whether the terms actually agree on one propagation
//! record — the only thing that makes the cross-QoI correlation real.

use core::fmt::Write as _;
use std::collections::BTreeMap;

use fs_evidence::{
    NoUsefulBound,
    action::ActionKind,
    uncertainty::{
        BudgetContribution, ComplianceVerdict, EngineeringUncertaintyKind,
        EngineeringUncertaintyTerm, RequirementRelation, TermValue,
    },
};
use fs_package::{Claim, EvidencePackage, PackageError};
use fs_project::ProjectDecisionAuthority;
use fs_regime::{OutputClaimReceipt, ProductOutputAudit};
use fs_session::DecisionAssessment;

/// Estimator identity used by package declarations that retain demotion receipts.
pub const REGIME_DEMOTION_PACKAGE_ESTIMATOR: &str = "fs-regime/output-demotion-package-receipt-v1";

/// Coupled human-readable and machine-readable projections of one regime audit.
#[derive(Debug, Clone, PartialEq)]
pub struct RegimeAuditOutputs {
    /// Deterministic reviewer-facing no-claim section, absent only in-domain.
    pub no_claims_markdown: Option<String>,
    /// Evidence package with every demotion receipt retained.
    pub package: EvidencePackage,
}

/// Refusal while adding final-envelope demotion receipts to a package.
#[derive(Debug)]
pub enum RegimePackageError {
    /// An existing claim reused the deterministic receipt claim id for
    /// different declaration bytes.
    ClaimIdConflict {
        /// Deterministic receipt claim identity.
        claim_id: String,
    },
    /// The resulting package exceeded its normal structural or transport gate.
    Package(PackageError),
}

impl core::fmt::Display for RegimePackageError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::ClaimIdConflict { claim_id } => write!(
                f,
                "package already contains a different claim at regime receipt id {claim_id:?}"
            ),
            Self::Package(error) => write!(f, "regime receipt package refused: {error}"),
        }
    }
}

impl std::error::Error for RegimePackageError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::ClaimIdConflict { .. } => None,
            Self::Package(error) => Some(error),
        }
    }
}

/// Deterministic package claim id for one receipt in one audit collection.
#[must_use]
pub fn regime_demotion_package_claim_id(
    audit: &ProductOutputAudit,
    receipt: &OutputClaimReceipt,
) -> String {
    format!(
        "regime-output-demotion/{:016x}/{}",
        audit.provenance.0,
        receipt.content_id()
    )
}

fn regime_demotion_package_statement(
    audit: &ProductOutputAudit,
    receipt: &OutputClaimReceipt,
) -> String {
    format!(
        "{{\"schema\":\"fs-regime-output-demotion-package-v1\",\"audit_provenance\":\"{:016x}\",\"receipt_content_id\":\"{}\",\"receipt\":{}}}",
        audit.provenance.0,
        receipt.content_id(),
        receipt.to_canonical_json()
    )
}

/// Retain every demoted audit receipt as an explicitly unbounded package claim.
///
/// Each declaration is `Estimated` with infinite dispersion. Its statement
/// binds the audit collection provenance, the domain-separated receipt
/// identity, and the exact canonical receipt JSON. Fully in-domain receipts
/// add nothing. Exact retries are idempotent; a deterministic claim-id
/// collision with different bytes refuses.
///
/// The wrapper is structural retention only. It does not use the portable
/// semantic-witness path, authenticate model-card authorities, or restore any
/// evidence color.
///
/// # Errors
///
/// Returns [`RegimePackageError::ClaimIdConflict`] for a conflicting retained
/// id or [`RegimePackageError::Package`] when the resulting package exceeds its
/// ordinary bounded transport gate.
pub fn retain_regime_demotions_in_package(
    mut package: EvidencePackage,
    audit: &ProductOutputAudit,
) -> Result<EvidencePackage, RegimePackageError> {
    let mut receipts = audit
        .receipts
        .iter()
        .filter(|receipt| receipt.demoted())
        .collect::<Vec<_>>();
    receipts.sort_by(|left, right| {
        left.qoi
            .cmp(&right.qoi)
            .then_with(|| left.content_id().cmp(&right.content_id()))
    });

    for receipt in receipts {
        let claim_id = regime_demotion_package_claim_id(audit, receipt);
        let claim = Claim::estimated(
            claim_id.clone(),
            regime_demotion_package_statement(audit, receipt),
            REGIME_DEMOTION_PACKAGE_ESTIMATOR,
            f64::INFINITY,
        );
        if let Some(existing) = package
            .declared_claims_unverified()
            .iter()
            .find(|candidate| candidate.id() == claim_id)
        {
            if existing != &claim {
                return Err(RegimePackageError::ClaimIdConflict { claim_id });
            }
            continue;
        }
        package = package.with_claim(claim);
    }

    package
        .try_merkle_root()
        .map_err(RegimePackageError::Package)?;
    Ok(package)
}

/// Project one final-envelope audit into both mandatory product outputs.
///
/// The returned Markdown and package are derived from the same immutable audit.
/// This prevents orchestration from rendering a demotion without retaining its
/// receipt, or retaining a receipt without surfacing the no-claim boundary.
/// Fully in-domain audits return no Markdown and leave the package unchanged.
///
/// # Errors
///
/// Returns the same bounded package or claim-id refusal as
/// [`retain_regime_demotions_in_package`]. No partial projection is returned.
pub fn project_regime_audit_outputs(
    package: EvidencePackage,
    audit: &ProductOutputAudit,
) -> Result<RegimeAuditOutputs, RegimePackageError> {
    let no_claims_markdown = regime_no_claims_markdown(audit);
    let package = retain_regime_demotions_in_package(package, audit)?;
    Ok(RegimeAuditOutputs {
        no_claims_markdown,
        package,
    })
}

/// Render every demoted final-envelope receipt in a deterministic no-claim section.
///
/// Fully in-domain QoIs are omitted. Each demoted entry includes the human
/// diagnosis, the strong receipt identity, and the exact canonical receipt JSON
/// needed for ledger/package handoff. This is presentation only: it cannot
/// restore color or authenticate model-card authorities.
#[must_use]
pub fn regime_no_claims_markdown(audit: &ProductOutputAudit) -> Option<String> {
    let mut receipts = audit
        .receipts
        .iter()
        .filter(|receipt| receipt.demoted())
        .collect::<Vec<_>>();
    if receipts.is_empty() {
        return None;
    }
    receipts.sort_by(|left, right| {
        left.qoi
            .cmp(&right.qoi)
            .then_with(|| left.content_id().cmp(&right.content_id()))
    });

    let mut output = String::from("## Operating-envelope no-claim boundaries\n\n");
    let _ = writeln!(
        output,
        "- **Audit collection provenance:** `{:016x}`\n",
        audit.provenance.0
    );
    for receipt in receipts {
        let Some(summary) = receipt.no_claim_markdown() else {
            continue;
        };
        let _ = writeln!(output, "{summary}");
        let _ = writeln!(
            output,
            "  - **Receipt identity:** `{}`",
            receipt.content_id()
        );
        output.push_str("  - **Exact canonical receipt:**\n\n");
        let _ = writeln!(
            output,
            "    ```json\n    {}\n    ```\n",
            receipt.to_canonical_json()
        );
    }
    output.push_str(
        "_Projection only: this section cannot authenticate model-card or calibration authorities, and an acknowledged override cannot restore evidence color._\n",
    );
    Some(output)
}

/// Render project requirement lineage and context gating for one lower-layer
/// verdict without assembling or recomputing a decision assessment.
#[must_use]
pub fn project_decision_gate_markdown(
    authority: &ProjectDecisionAuthority,
    compliance: &ComplianceVerdict,
) -> String {
    let requirement = authority.requirement();
    let scalar = requirement.scalar();
    let context = authority.context();
    let verdict = match compliance {
        ComplianceVerdict::Compliant { .. } => "compliant",
        ComplianceVerdict::NonCompliant { .. } => "non-compliant",
        ComplianceVerdict::Indeterminate { .. } => "indeterminate",
    };
    let gate_outcome = if matches!(compliance, ComplianceVerdict::Indeterminate { .. })
        && !context.permits_indeterminate()
    {
        "refused: this context requires a determinate assessment"
    } else {
        "admitted"
    };
    let relation = match scalar.relation() {
        RequirementRelation::AtMost => "at-most",
        RequirementRelation::AtLeast => "at-least",
    };

    format!(
        "## Project decision gate: `{qoi}`\n\n- **Project:** `{project}` (created `{created}`)\n- **Context of use:** {context_of_use}\n- **Intended decision:** {intended_decision}\n- **Gate:** `{gate}`\n- **Consequence:** `{consequence}`\n- **Lower-layer verdict:** `{verdict}`\n- **Gate outcome:** **{gate_outcome}**\n- **Effective requirement:** `{relation}` `{limit}` `{unit}`\n- **Requirement source:** `{requirement_kind}` document `{requirement_document}` version `{requirement_version}` locator `{requirement_locator}`; artifact `{requirement_artifact}`\n- **Safety-factor policy:** factor `{factor}` from `{factor_kind}` document `{factor_document}` version `{factor_version}` locator `{factor_locator}`; artifact `{factor_artifact}`\n- **Context artifact:** `{context_id}` at `{context_hash}`\n\n_Projection only: this block reports declared project intent and an existing verdict. It does not authenticate either source, recompute compliance, or turn an admitted scoping result into sign-off authority._\n",
        qoi = scalar.qoi(),
        project = context.project_name(),
        created = context.created(),
        context_of_use = context.context_of_use(),
        intended_decision = context.intended_decision(),
        gate = context.decision_gate().slug(),
        consequence = context.consequence().slug(),
        unit = scalar.unit(),
        limit = scalar.limit(),
        requirement_kind = requirement.source().kind().slug(),
        requirement_document = requirement.source().document(),
        requirement_version = requirement.source().version(),
        requirement_locator = requirement.source().locator(),
        requirement_artifact = scalar.provenance().digest(),
        factor = requirement.safety_factor().value(),
        factor_kind = requirement.safety_factor_source().kind().slug(),
        factor_document = requirement.safety_factor_source().document(),
        factor_version = requirement.safety_factor_source().version(),
        factor_locator = requirement.safety_factor_source().locator(),
        factor_artifact = requirement.safety_factor().policy().digest(),
        context_id = context.artifact().id(),
        context_hash = context.artifact().hash(),
    )
}

/// Render one already-validated decision assessment as deterministic Markdown.
///
/// The headline keeps the tri-state verdict, effective sourced requirement,
/// safety-factor policy, context, evidence identities, flip conditions, paired
/// attribution views, and replay root together. The indented audit projection
/// is the exact lower-layer explanation; this function is presentation only.
#[must_use]
pub fn decision_headline_markdown<Q>(assessment: &DecisionAssessment<Q>) -> String {
    let quantity = assessment.quantity();
    let unit = quantity.unit();
    let mut output = format!("## Decision headline: `{}`\n\n", quantity.qoi());

    render_verdict(assessment.compliance(), unit, &mut output);
    render_decision_authorities(assessment, &mut output);
    render_flip_conditions(assessment, unit, &mut output);
    render_attribution_headline(assessment, unit, &mut output);
    render_exact_audit(assessment, &mut output);
    output
}

fn render_verdict(verdict: &ComplianceVerdict, unit: &str, output: &mut String) {
    match verdict {
        ComplianceVerdict::Compliant { margin, .. } => {
            let _ = writeln!(
                output,
                "- **Verdict:** `compliant` with residual margin `{margin} {unit}`"
            );
        }
        ComplianceVerdict::NonCompliant { shortfall, .. } => {
            let _ = writeln!(
                output,
                "- **Verdict:** `non-compliant` with residual shortfall `{shortfall} {unit}`"
            );
        }
        ComplianceVerdict::Indeterminate {
            known_lower,
            known_upper,
            no_useful_bound,
            ..
        } => {
            let _ = writeln!(
                output,
                "- **Verdict:** `indeterminate`; known band `[{known_lower}, {known_upper}] {unit}`"
            );
            if let Some(refusal) = no_useful_bound {
                output.push('\n');
                output.push_str(&no_useful_bound_markdown(refusal));
            }
        }
    }
}

/// Render a typed usefulness refusal in its own visual class.
///
/// This projection cannot produce a certificate claim or a binary requirement
/// verdict.
#[must_use]
pub fn no_useful_bound_markdown(refusal: &NoUsefulBound) -> String {
    let criterion = refusal.criterion();
    format!(
        "### NoUsefulBound\n\n> The retained enclosure is valid, but it is not useful for the declared engineering decision.\n\n- **Cause:** `{}`\n- **Achieved enclosure:** `[{lower}, {upper}] {unit}`\n- **Achieved width:** `{width} {unit}`\n- **Required maximum width:** `{threshold} {unit}`\n- **Decision context:** `{context}`\n- **Suggested E09 reformulation:** `{suggestion}` ({suggestion_title})\n- **No-claim boundary:** no compliance verdict, scientific color, or certificate is minted from this refusal.\n",
        refusal.cause().code(),
        lower = refusal.interval().lower(),
        upper = refusal.interval().upper(),
        unit = criterion.unit(),
        width = refusal.width_achieved(),
        threshold = criterion.max_width(),
        context = criterion.decision_context(),
        suggestion = refusal.suggested_reformulation().code(),
        suggestion_title = refusal.suggested_reformulation().title(),
    )
}

fn render_decision_authorities<Q>(assessment: &DecisionAssessment<Q>, output: &mut String) {
    let quantity = assessment.quantity();
    let requirement = assessment.requirement();
    let scalar = requirement.scalar();
    let unit = quantity.unit();
    let relation = match scalar.relation() {
        RequirementRelation::AtMost => "at most",
        RequirementRelation::AtLeast => "at least",
    };
    let _ = writeln!(
        output,
        "- **Effective requirement:** `{}` — {relation} `{}` `{unit}`",
        scalar.id(),
        scalar.limit()
    );
    let _ = writeln!(
        output,
        "- **Declared safety factor:** `{}` (already reflected in the effective limit)",
        requirement.safety_factor().value()
    );
    let _ = writeln!(
        output,
        "- **Requirement authority:** `{}` document `{}` version `{}` locator `{}`; artifact `{}` at `{}`",
        requirement.source().kind().slug(),
        requirement.source().document(),
        requirement.source().version(),
        requirement.source().locator(),
        scalar.provenance().role(),
        scalar.provenance().digest()
    );
    let _ = writeln!(
        output,
        "- **Safety-factor policy:** `{}` document `{}` version `{}` locator `{}`; artifact `{}` at `{}`",
        requirement.safety_factor_source().kind().slug(),
        requirement.safety_factor_source().document(),
        requirement.safety_factor_source().version(),
        requirement.safety_factor_source().locator(),
        requirement.safety_factor().policy().role(),
        requirement.safety_factor().policy().digest()
    );
    let _ = writeln!(
        output,
        "- **Quantity evidence:** schema `{}` at `{}`",
        quantity.schema(),
        quantity.artifact()
    );
    let _ = writeln!(
        output,
        "- **Context of use:** `{}` at `{}`",
        assessment.context().id(),
        assessment.context().hash()
    );
    let _ = writeln!(
        output,
        "- **Decision assessment:** `{}`",
        assessment.content_hash()
    );
    let _ = writeln!(
        output,
        "- **Replay package:** `{}`\n",
        assessment.replay_package()
    );
}

fn render_flip_conditions<Q>(assessment: &DecisionAssessment<Q>, unit: &str, output: &mut String) {
    output.push_str("### What could change this verdict\n\n");
    if assessment.flip_conditions().is_empty() {
        output.push_str("_No admitted unknown is reported as verdict-flipping._\n\n");
    } else {
        for unknown in assessment.flip_conditions() {
            let _ = writeln!(
                output,
                "- `{}`: adverse magnitude `{}` `{unit}`; suggested evidence `{}`",
                unknown.kind().name(),
                unknown.required_magnitude(),
                action_kind_name(unknown.suggested_action())
            );
        }
        let _ = writeln!(
            output,
            "\nThe assessment retains {} explicit evidence recommendation(s).\n",
            assessment.actions().len()
        );
    }
}

fn render_attribution_headline<Q>(
    assessment: &DecisionAssessment<Q>,
    unit: &str,
    output: &mut String,
) {
    output.push_str("### Attribution headline\n\n");
    match assessment.largest_known_budget_link() {
        Some(link) => match link.contribution() {
            BudgetContribution::Known {
                conservative_half_width,
                share_of_known,
            } => {
                let share = share_of_known
                    .map_or_else(|| "not-defined".to_string(), |value| format!("{value}"));
                let _ = writeln!(
                    output,
                    "- Largest finite budget group: `{}`; half-width `{conservative_half_width} {unit}`; share `{share}`",
                    link.group().label()
                );
            }
            BudgetContribution::Unknown { .. } => {
                output.push_str("- Largest finite budget group: unavailable\n");
            }
        },
        None => output.push_str("- Largest finite budget group: unavailable\n"),
    }
    match assessment.strongest_decision_link() {
        Some(link) => {
            let _ = writeln!(
                output,
                "- Strongest one-group-at-a-time decision influence: `{}`; signed-separation shift `{}` `{unit}`",
                link.group().label(),
                link.influence()
            );
        }
        None => output.push_str("- Strongest decision influence: unavailable\n"),
    }
    if assessment.attribution().headline_disagrees() {
        output.push_str(
            "- **Paired-view warning:** the largest budget magnitude is not the strongest decision influence.\n",
        );
    }
}

fn render_exact_audit<Q>(assessment: &DecisionAssessment<Q>, output: &mut String) {
    output.push_str("\n### Exact audit projection\n\n");
    for line in assessment.render_explain().lines() {
        let _ = writeln!(output, "    {line}");
    }
    output.push_str(
        "\n_Projection only: this report does not certify evidence, resolve artifact hashes, recompute the verdict, or authenticate requirement authorities._\n",
    );
}

const fn action_kind_name(kind: ActionKind) -> &'static str {
    match kind {
        ActionKind::SolverTolerance => "solver-tolerance",
        ActionKind::MeshRefinement => "mesh-refinement",
        ActionKind::TimeRefinement => "time-refinement",
        ActionKind::RepresentationEscalation => "representation-escalation",
        ActionKind::UqSamples => "uq-samples",
        ActionKind::MaterialCouponTest => "material-coupon-test",
        ActionKind::SensorCampaign => "sensor-campaign",
        ActionKind::Falsification => "falsification",
        ActionKind::StandardsObligation => "standards-obligation",
        ActionKind::Refusal => "refusal",
        _ => "unsupported",
    }
}

/// A dimensioned quantity — a value with its unit (units on every value).
#[derive(Debug, Clone, PartialEq)]
pub struct Quantity {
    /// The numeric value.
    pub value: f64,
    /// The unit label (e.g. `"mm"`, `"kg"`, `"1/mm"`).
    pub unit: String,
}

impl Quantity {
    /// A quantity.
    #[must_use]
    pub fn new(value: f64, unit: impl Into<String>) -> Quantity {
        Quantity {
            value,
            unit: unit.into(),
        }
    }
}

/// One replayable operation of the reproducibility IR.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReproStep {
    /// The operation name.
    pub op: String,
    /// Its serialized arguments.
    pub args: Vec<String>,
}

/// A notebook block.
#[derive(Debug, Clone, PartialEq)]
pub enum Block {
    /// Free prose.
    Prose(String),
    /// A named, dimensioned metric.
    Metric {
        /// The metric name.
        name: String,
        /// The value + unit.
        quantity: Quantity,
    },
    /// A reproducibility step (part of the replay IR).
    Step(ReproStep),
}

/// An automatic lab notebook for a study.
#[derive(Debug, Clone, PartialEq)]
pub struct LabNotebook {
    /// The study title.
    pub title: String,
    /// The RNG seed (provenance).
    pub seed: u64,
    /// The toolchain / crate version (provenance).
    pub version: String,
    /// The report body.
    pub blocks: Vec<Block>,
}

impl LabNotebook {
    /// A new notebook with provenance.
    #[must_use]
    pub fn new(title: impl Into<String>, seed: u64, version: impl Into<String>) -> LabNotebook {
        LabNotebook {
            title: title.into(),
            seed,
            version: version.into(),
            blocks: Vec::new(),
        }
    }

    /// Append prose.
    pub fn prose(&mut self, text: impl Into<String>) -> &mut LabNotebook {
        self.blocks.push(Block::Prose(text.into()));
        self
    }

    /// Append a dimensioned metric.
    pub fn metric(
        &mut self,
        name: impl Into<String>,
        value: f64,
        unit: impl Into<String>,
    ) -> &mut LabNotebook {
        self.blocks.push(Block::Metric {
            name: name.into(),
            quantity: Quantity::new(value, unit),
        });
        self
    }

    /// Append a reproducibility step.
    pub fn step(&mut self, op: impl Into<String>, args: Vec<String>) -> &mut LabNotebook {
        self.blocks.push(Block::Step(ReproStep {
            op: op.into(),
            args,
        }));
        self
    }

    /// The metrics recorded (name + quantity).
    #[must_use]
    pub fn metrics(&self) -> Vec<(&str, &Quantity)> {
        self.blocks
            .iter()
            .filter_map(|b| match b {
                Block::Metric { name, quantity } => Some((name.as_str(), quantity)),
                _ => None,
            })
            .collect()
    }

    /// THE EXACT IR TO REPRODUCE the study — the ordered replay steps.
    #[must_use]
    pub fn repro_ir(&self) -> Vec<ReproStep> {
        self.blocks
            .iter()
            .filter_map(|b| match b {
                Block::Step(s) => Some(s.clone()),
                _ => None,
            })
            .collect()
    }

    /// The report rendered to Markdown (deterministic).
    #[must_use]
    pub fn render_markdown(&self) -> String {
        let mut s = String::new();
        let _ = writeln!(s, "# {}", self.title);
        let _ = writeln!(s);
        let _ = writeln!(s, "_seed: {} · version: {}_", self.seed, self.version);
        let _ = writeln!(s);
        for block in &self.blocks {
            match block {
                Block::Prose(t) => {
                    let _ = writeln!(s, "{t}");
                    let _ = writeln!(s);
                }
                Block::Metric { name, quantity } => {
                    let _ = writeln!(s, "- **{}**: {} {}", name, quantity.value, quantity.unit);
                }
                Block::Step(step) => {
                    let _ = writeln!(s, "- repro: `{}({})`", step.op, step.args.join(", "));
                }
            }
        }
        s
    }

    /// A content hash of the report STRUCTURE — a report is as
    /// content-addressed as any other ledger artifact. Canonical
    /// replay identity encoding (gp3.14): the former hash of the
    /// RENDERED Markdown was non-injective — a Prose block containing
    /// `- **name**: value unit` rendered byte-identically to a Metric
    /// block, so structurally different notebooks could share a
    /// content address (gated in the battery). The Markdown render
    /// remains the human artifact; the hash binds the typed fields.
    #[must_use]
    pub fn content_hash(&self) -> u64 {
        let mut b = fs_obs::ident::IdentityBuilder::new("lab-notebook")
            .str("title", &self.title)
            .u64("seed", self.seed)
            .str("version", &self.version);
        for block in &self.blocks {
            b = match block {
                Block::Prose(t) => b.str("prose", t),
                Block::Metric { name, quantity } => b
                    .str("metric", name)
                    .f64_bits("value", quantity.value)
                    .str("unit", &quantity.unit),
                Block::Step(step) => {
                    let mut sb = b.str("step_op", &step.op);
                    for arg in &step.args {
                        sb = sb.str("step_arg", arg);
                    }
                    sb
                }
            };
        }
        b.finish().root()
    }
}

/// A per-feature semantic difference between two designs.
#[derive(Debug, Clone, PartialEq)]
pub struct FeatureDelta {
    /// The feature name.
    pub name: String,
    /// The value before.
    pub before: f64,
    /// The value after.
    pub after: f64,
    /// The absolute change (`after − before`).
    pub abs_change: f64,
    /// The relative change (`abs_change / before`; `0` if `before == 0`).
    pub rel_change: f64,
    /// The unit.
    pub unit: String,
}

impl FeatureDelta {
    /// A human attribution string, e.g. `"wall_thickness: 2 mm → 1.6 mm (−20.0%)"`.
    #[must_use]
    pub fn describe(&self) -> String {
        let mut s = String::new();
        let _ = write!(
            s,
            "{}: {} {} → {} {} ({:+.1}%)",
            self.name,
            self.before,
            self.unit,
            self.after,
            self.unit,
            self.rel_change * 100.0
        );
        s
    }
}

/// A SEMANTIC (per-feature) diff between two designs described as
/// `feature → Quantity` maps: the changed features with absolute + relative
/// deltas, ranked by significance (largest relative change first). Not a file
/// diff — a geometric attribution.
#[must_use]
pub fn semantic_diff(
    before: &BTreeMap<String, Quantity>,
    after: &BTreeMap<String, Quantity>,
) -> Vec<FeatureDelta> {
    let mut deltas: Vec<FeatureDelta> = before
        .iter()
        .filter_map(|(name, b)| {
            after.get(name).map(|a| {
                let abs_change = a.value - b.value;
                let rel_change = if b.value == 0.0 {
                    0.0
                } else {
                    abs_change / b.value
                };
                FeatureDelta {
                    name: name.clone(),
                    before: b.value,
                    after: a.value,
                    abs_change,
                    rel_change,
                    unit: b.unit.clone(),
                }
            })
        })
        .collect();
    // rank by significance (largest |relative change| first); name as tiebreak.
    deltas.sort_by(|x, y| {
        y.rel_change
            .abs()
            .partial_cmp(&x.rel_change.abs())
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| x.name.cmp(&y.name))
    });
    deltas
}

// ---------------------------------------------------------------------------
// Nominal-versus-as-built QoI rendering (bead f85xj.12.2 item (d))
// ---------------------------------------------------------------------------

/// Refusal while building a nominal-versus-as-built comparison row.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AsBuiltDeltaError {
    /// The cited term is not the budget's geometry slot.
    NotGeometryTerm {
        /// The kind actually supplied.
        kind: EngineeringUncertaintyKind,
    },
    /// A rendered value was not finite.
    NonFiniteValue {
        /// Which field.
        field: &'static str,
    },
}

impl core::fmt::Display for AsBuiltDeltaError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            AsBuiltDeltaError::NotGeometryTerm { kind } => write!(
                f,
                "the as-built comparison must cite the geometry term, not the {} term",
                kind.name()
            ),
            AsBuiltDeltaError::NonFiniteValue { field } => {
                write!(f, "{field} is not finite")
            }
        }
    }
}

impl core::error::Error for AsBuiltDeltaError {}

/// One QoI compared between a nominal-placement solve and an as-built-placement
/// solve, carrying the geometry term that attributes the difference.
///
/// The term is borrowed rather than copied so the row cannot drift from the
/// budget it claims to describe.
#[derive(Debug, Clone, PartialEq)]
pub struct AsBuiltQoiDelta<'a> {
    qoi: String,
    unit: String,
    nominal: f64,
    as_built: f64,
    geometry: &'a EngineeringUncertaintyTerm,
}

impl<'a> AsBuiltQoiDelta<'a> {
    /// Build one comparison row.
    ///
    /// # Errors
    /// The cited term is not the geometry slot, or a value is not finite.
    pub fn try_new(
        qoi: impl Into<String>,
        unit: impl Into<String>,
        nominal: f64,
        as_built: f64,
        geometry: &'a EngineeringUncertaintyTerm,
    ) -> Result<Self, AsBuiltDeltaError> {
        if geometry.kind() != EngineeringUncertaintyKind::Geometry {
            return Err(AsBuiltDeltaError::NotGeometryTerm {
                kind: geometry.kind(),
            });
        }
        if !nominal.is_finite() {
            return Err(AsBuiltDeltaError::NonFiniteValue { field: "nominal" });
        }
        if !as_built.is_finite() {
            return Err(AsBuiltDeltaError::NonFiniteValue { field: "as_built" });
        }
        Ok(Self {
            qoi: qoi.into(),
            unit: unit.into(),
            nominal,
            as_built,
            geometry,
        })
    }

    /// The QoI name.
    #[must_use]
    pub fn qoi(&self) -> &str {
        &self.qoi
    }

    /// The budget unit.
    #[must_use]
    pub fn unit(&self) -> &str {
        &self.unit
    }

    /// The nominal-placement value.
    #[must_use]
    pub const fn nominal(&self) -> f64 {
        self.nominal
    }

    /// The as-built-placement value.
    #[must_use]
    pub const fn as_built(&self) -> f64 {
        self.as_built
    }

    /// The geometry term attributing the shift.
    #[must_use]
    pub const fn geometry(&self) -> &EngineeringUncertaintyTerm {
        self.geometry
    }

    /// The as-built minus nominal shift, in the budget unit.
    #[must_use]
    pub fn shift(&self) -> f64 {
        self.as_built - self.nominal
    }

    /// The geometry term's conservative half-width, when the term carries one.
    ///
    /// `None` whenever the term is `Unknown` — which is exactly what a rejected
    /// linearization spot-check produces upstream — or any representation that
    /// publishes no half-width.
    ///
    /// `IntervalBound` brackets the half-width itself rather than the value, so
    /// the conservative reading is its upper endpoint. That is the same choice
    /// `fs-evidence` makes when it aggregates a marginal, and taking the
    /// midpoint here instead would quietly report less uncertainty than the
    /// budget does.
    #[must_use]
    pub fn geometry_half_width(&self) -> Option<f64> {
        match self.geometry.value() {
            TermValue::Distribution(distribution) => Some(distribution.conservative_half_width),
            TermValue::IntervalBound { upper, .. } => Some(*upper),
            _ => None,
        }
    }

    /// Whether the shift is larger than the geometry uncertainty that the
    /// as-built measurement itself carries.
    ///
    /// `None` when the geometry term publishes no half-width, so the question
    /// is unanswerable rather than answered negatively. A `Some(false)` says
    /// the shift sits inside the pose-propagated band: the two solves differ,
    /// but this measurement cannot resolve the difference from its own noise.
    #[must_use]
    pub fn shift_resolved(&self) -> Option<bool> {
        self.geometry_half_width()
            .map(|half_width| self.shift().abs() > half_width)
    }

    /// A human attribution string.
    #[must_use]
    pub fn describe(&self) -> String {
        let mut s = String::new();
        let _ = write!(
            s,
            "{}: {} {} → {} {} (shift {:+} {})",
            self.qoi,
            self.nominal,
            self.unit,
            self.as_built,
            self.unit,
            self.shift(),
            self.unit
        );
        s
    }
}

/// Render a nominal-versus-as-built QoI table with the geometry term attributed
/// and the cross-QoI correlation structure stated.
///
/// Rows are rendered in the order supplied: the caller's declaration order is
/// the report order, so a given study renders identically on every run.
///
/// The correlation line is the load-bearing part. Upstream, one pose covariance
/// propagates into every QoI through a single content-addressed record, and each
/// per-QoI geometry term cites that record. This renderer therefore does not
/// assert correlation — it reports whether the terms in front of it actually
/// agree on one record, and says so plainly when they do not.
///
/// This is presentation only: it recomputes no solve, no propagation, and no
/// compliance verdict.
#[must_use]
pub fn nominal_vs_as_built_markdown(rows: &[AsBuiltQoiDelta<'_>]) -> String {
    let mut output = String::from("## Nominal versus as-built\n\n");
    if rows.is_empty() {
        output.push_str(
            "- **No rows:** no QoI was compared, so this section makes no as-built claim.\n",
        );
        return output;
    }

    for row in rows {
        let _ = writeln!(
            output,
            "- **`{}`:** `{} {}` → `{} {}`; shift `{:+} {}`",
            row.qoi(),
            row.nominal(),
            row.unit(),
            row.as_built(),
            row.unit(),
            row.shift(),
            row.unit()
        );
        let provenance = row.geometry().provenance();
        match row.geometry().value() {
            TermValue::Unknown { reason } => {
                let _ = writeln!(
                    output,
                    "    - geometry term: `unknown` ({reason}); the shift is NOT attributed"
                );
            }
            _ => {
                if let Some(half_width) = row.geometry_half_width() {
                    let _ = writeln!(
                        output,
                        "    - geometry term: conservative half-width `{half_width} {}`",
                        row.unit()
                    );
                }
            }
        }
        let resolution = match row.shift_resolved() {
            Some(true) => "resolved: the shift exceeds the pose-propagated half-width",
            Some(false) => "NOT resolved: the shift sits inside the pose-propagated half-width",
            None => "not decidable: the geometry term publishes no half-width",
        };
        let _ = writeln!(output, "    - {resolution}");
        let _ = writeln!(
            output,
            "    - attribution: `{}@{}`",
            provenance.role(),
            provenance.digest()
        );
    }

    output.push('\n');
    let mut records: Vec<String> = rows
        .iter()
        .map(|row| row.geometry().provenance().digest().to_hex())
        .collect();
    records.sort();
    records.dedup();
    match records.as_slice() {
        [single] if rows.len() > 1 => {
            let _ = writeln!(
                output,
                "- **Correlation:** all {} geometry terms cite one propagation record `{single}`, so the pose error is carried as a single correlated block rather than as independent per-QoI noise.",
                rows.len()
            );
        }
        [single] => {
            let _ = writeln!(
                output,
                "- **Correlation:** one QoI, one propagation record `{single}`; no cross-QoI correlation is exercised by this table."
            );
        }
        many => {
            let _ = writeln!(
                output,
                "- **Correlation:** the geometry terms cite {} distinct propagation records, so NO single correlated block covers this table and cross-QoI correlation is not established here.",
                many.len()
            );
        }
    }
    output.push_str(
        "- **No-claim boundary:** the shift is a difference between two computed values, not a measurement of the physical part. The geometry term bounds only the pose-propagated contribution and covers no discretization, parameter, boundary-condition, or model-form error. A shift inside the half-width is not evidence that nothing changed.\n",
    );
    output
}
