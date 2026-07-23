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

use core::fmt::Write as _;
use std::collections::BTreeMap;

use fs_evidence::{
    action::ActionKind,
    uncertainty::{BudgetContribution, ComplianceVerdict, RequirementRelation},
};
use fs_session::DecisionAssessment;

/// Render one already-validated decision assessment as deterministic Markdown.
///
/// The headline keeps the tri-state verdict, effective sourced requirement,
/// safety-factor policy, context, evidence identities, flip conditions, paired
/// attribution views, and replay root together. The indented audit projection
/// is the exact lower-layer explanation; this function is presentation only.
#[must_use]
pub fn decision_headline_markdown<Q>(assessment: &DecisionAssessment<Q>) -> String {
    let quantity = assessment.quantity();
    let requirement = assessment.requirement();
    let scalar = requirement.scalar();
    let unit = quantity.unit();
    let relation = match scalar.relation() {
        RequirementRelation::AtMost => "at most",
        RequirementRelation::AtLeast => "at least",
    };
    let mut output = format!("## Decision headline: `{}`\n\n", quantity.qoi());

    match assessment.compliance() {
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
            ..
        } => {
            let _ = writeln!(
                output,
                "- **Verdict:** `indeterminate`; known band `[{known_lower}, {known_upper}] {unit}`"
            );
        }
    }
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
        "- **Requirement authority:** `{}` at `{}`",
        scalar.provenance().role(),
        scalar.provenance().digest()
    );
    let _ = writeln!(
        output,
        "- **Safety-factor policy:** `{}` at `{}`",
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

    output.push_str("\n### Exact audit projection\n\n");
    for line in assessment.render_explain().lines() {
        let _ = writeln!(output, "    {line}");
    }
    output.push_str(
        "\n_Projection only: this report does not certify evidence, resolve artifact hashes, recompute the verdict, or authenticate requirement authorities._\n",
    );
    output
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
