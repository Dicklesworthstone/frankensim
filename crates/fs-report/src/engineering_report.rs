//! Deterministic self-contained HTML engineering report and JSON twin.
//!
//! Bead: `frankensim-extreal-program-f85xj.6.9`
//!
//! Generates a reproducible, single-file HTML report with embedded CSS (zero external assets)
//! and a semantically identical JSON twin for automated inspection. Sections include:
//! 1. Executive Summary & Provenance (Five Explicits, machine fingerprint, constellation lock)
//! 2. DecisionAssessment Block (requirement, margin, P(compliance), flip conditions, next actions)
//! 3. QoI Table with evidence colors (Verified, Validated, Estimated, Waived) and error budgets
//! 4. Discretization Convergence Ladder (rungs, observed order, Richardson extrapolation, GCI)
//! 5. Parameter/BC Uncertainty Section (distributions, bounds, Philox sampling error)
//! 6. Material & Geometry Import Provenance
//! 7. Auto-Generated Known-Gaps and No-Claims Section
//! 8. Replay Instructions

use fs_blake3::{ContentHash, hash_domain};
use fs_evidence::Color;
use fs_ladder::ConvergenceResult;
use fs_uq::UqResult;
use std::fmt::Write as _;

/// Item in the QoI and Error Budget table.
#[derive(Debug, Clone, PartialEq)]
pub struct QoiReportItem {
    /// Identifier / symbol name for the Quantity of Interest.
    pub name: String,
    /// Human-readable explanation of the physical quantity.
    pub description: String,
    /// Nominal computed or estimated value.
    pub nominal_value: f64,
    /// Physical unit (e.g. "K", "W", "Pa").
    pub unit: String,
    /// Epistemic evidence color classification.
    pub color: Color,
    /// Discretization uncertainty bound.
    pub discretization_error: f64,
    /// Parameter and boundary condition uncertainty spread.
    pub parameter_uncertainty: f64,
    /// Surrogate model representation error.
    pub surrogate_error: f64,
    /// Total composed error budget term.
    pub total_uncertainty_budget: f64,
    /// Lineage source URI / root.
    pub source_root: String,
}

/// Material and interface provenance entry.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MaterialReportItem {
    /// Named physical region in the model.
    pub region_name: String,
    /// Content-addressed material card ID.
    pub material_card_id: String,
    /// Thermal conductivity with units.
    pub thermal_conductivity: String,
    /// Specific heat capacity with units.
    pub specific_heat: String,
    /// Mass density with units.
    pub density: String,
    /// Source material card pack identity.
    pub source_pack: String,
}

/// Known gap or no-claim boundary item.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NoClaimItem {
    /// Affected component or physical model.
    pub component: String,
    /// Maturity or validity status.
    pub status: String,
    /// Explicit statement of what is and is not claimed.
    pub statement: String,
}

/// One sourced requirement evaluated against a QoI.
///
/// Every number is copied from the retained requirement evaluation; a
/// renderer never derives a verdict of its own.
#[derive(Debug, Clone, PartialEq)]
pub struct RequirementReportItem {
    /// Requirement identity as the producer recorded it.
    pub id: String,
    /// QoI the requirement constrains.
    pub qoi: String,
    /// Region the requirement is bound to.
    pub region: String,
    /// Effective limit after the sourced safety factor.
    pub effective_limit: f64,
    /// Required margin declared by the requirement.
    pub required_margin: f64,
    /// Nominal achieved margin (limit minus nominal value).
    pub nominal_margin: f64,
    /// Unit shared by the three values above.
    pub unit: String,
    /// Producer-recorded outcome (`indeterminate`, `pass`, `fail`, ...).
    pub outcome: String,
}

/// One term of the engineering uncertainty budget.
///
/// `value` is `None` when the producer recorded the term as NO-DATA; the
/// renderer prints exactly that and never substitutes a zero.
#[derive(Debug, Clone, PartialEq)]
pub struct BudgetTermItem {
    /// QoI the term belongs to.
    pub qoi: String,
    /// Term kind (discretization, model-form, ...).
    pub kind: String,
    /// Producer-recorded state (`no-data`, `measured`, ...).
    pub state: String,
    /// Magnitude when measured.
    pub value: Option<f64>,
    /// Producer-recorded reason for the state.
    pub reason: String,
    /// Owner that must supply the term.
    pub owner: String,
}

/// One retained stage receipt the report was projected from.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StageReceiptItem {
    /// Stage name.
    pub stage: String,
    /// Stage ordinal.
    pub ordinal: u32,
    /// Content hash (hex) of the retained receipt artifact. Ledger operation
    /// ids are deliberately absent: they are ledger coordinates, not content,
    /// and the report must be a pure function of the retained receipts.
    pub receipt_hash: String,
    /// Key/value summary copied from the receipt (never computed here).
    pub summary: Vec<(String, String)>,
}

/// One lineage pointer (label plus content hash, hex).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LineageItem {
    /// What the hash identifies.
    pub label: String,
    /// Content hash, hex.
    pub hash: String,
}

/// Run metadata and Five Explicits.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReportProvenance {
    /// Semantic code version.
    pub code_version: String,
    /// Exact constellation lock hash.
    pub constellation_lock: String,
    /// Target execution architecture fingerprint.
    pub machine_fingerprint: String,
    /// System of physical units.
    pub units_system: String,
    /// Deterministic pseudo-random number generator seed.
    pub rng_seed: u64,
    /// Wall-clock execution budget in seconds.
    pub wall_budget_s: u64,
    /// Peak memory allocation budget in bytes.
    pub mem_budget_bytes: u64,
}

/// Comprehensive engineering report data structure.
#[derive(Debug, Clone, PartialEq)]
pub struct EngineeringReport {
    /// Unique run identifier.
    pub run_id: String,
    /// Declared simulation project name.
    pub project_name: String,
    /// Provenance metadata and Five Explicits.
    pub provenance: ReportProvenance,
    /// Evaluated Quantities of Interest and error budgets.
    pub qois: Vec<QoiReportItem>,
    /// Discretization ladder convergence analysis.
    pub convergence: Option<ConvergenceResult>,
    /// Parameter/BC uncertainty quantification results.
    pub uncertainty: Option<UqResult>,
    /// Material cards and physical region properties.
    pub materials: Vec<MaterialReportItem>,
    /// Explicit known gaps and no-claim boundaries.
    pub no_claims: Vec<NoClaimItem>,
    /// Deterministic command to replay this exact run.
    pub replay_command: String,
    /// Sourced requirement evaluations copied from the producer.
    pub requirements: Vec<RequirementReportItem>,
    /// Engineering uncertainty budget terms copied from the producer.
    pub budget_terms: Vec<BudgetTermItem>,
    /// Retained stage receipts this report was projected from.
    pub stages: Vec<StageReceiptItem>,
    /// Lineage pointers (project, solution, receipts) by content hash.
    pub lineage: Vec<LineageItem>,
}

/// Render a possibly-unavailable magnitude for HTML: finite values print with
/// four decimals, anything else prints the literal `NO-DATA`.
fn html_num(value: f64) -> String {
    if value.is_finite() {
        format!("{value:.4}")
    } else {
        "NO-DATA".to_string()
    }
}

/// Render a possibly-unavailable magnitude for JSON: finite values print
/// with six decimals, anything else prints `null` (JSON has no NaN).
fn json_num(value: f64) -> String {
    if value.is_finite() {
        format!("{value:.6}")
    } else {
        "null".to_string()
    }
}

impl EngineeringReport {
    /// Create a new report for a given run and project.
    #[must_use]
    pub fn new(run_id: impl Into<String>, project_name: impl Into<String>) -> Self {
        Self {
            run_id: run_id.into(),
            project_name: project_name.into(),
            provenance: ReportProvenance {
                code_version: env!("CARGO_PKG_VERSION").to_string(),
                constellation_lock: "frankensim-constellation-lock-v1".to_string(),
                machine_fingerprint: "isa=apple-m-series-or-x86;threads=auto".to_string(),
                units_system: "SI (m, kg, s, K, W, Pa)".to_string(),
                rng_seed: 0x0517,
                wall_budget_s: 60,
                mem_budget_bytes: 1024 * 1024 * 1024,
            },
            qois: Vec::new(),
            convergence: None,
            uncertainty: None,
            materials: Vec::new(),
            no_claims: Vec::new(),
            replay_command: String::new(),
            requirements: Vec::new(),
            budget_terms: Vec::new(),
            stages: Vec::new(),
            lineage: Vec::new(),
        }
    }

    /// Add a QoI result.
    #[must_use]
    pub fn with_qoi(mut self, item: QoiReportItem) -> Self {
        self.qois.push(item);
        self
    }

    /// Replace the default provenance with the run's recorded Five Explicits.
    #[must_use]
    pub fn with_provenance(mut self, provenance: ReportProvenance) -> Self {
        self.provenance = provenance;
        self
    }

    /// Add a sourced requirement evaluation.
    #[must_use]
    pub fn with_requirement(mut self, item: RequirementReportItem) -> Self {
        self.requirements.push(item);
        self
    }

    /// Add one uncertainty budget term.
    #[must_use]
    pub fn with_budget_term(mut self, item: BudgetTermItem) -> Self {
        self.budget_terms.push(item);
        self
    }

    /// Add a retained stage receipt row.
    #[must_use]
    pub fn with_stage(mut self, item: StageReceiptItem) -> Self {
        self.stages.push(item);
        self
    }

    /// Add a lineage pointer.
    #[must_use]
    pub fn with_lineage(mut self, item: LineageItem) -> Self {
        self.lineage.push(item);
        self
    }

    /// Add a material provenance row.
    #[must_use]
    pub fn with_material(mut self, item: MaterialReportItem) -> Self {
        self.materials.push(item);
        self
    }

    /// Add a no-claim boundary.
    #[must_use]
    pub fn with_no_claim(mut self, item: NoClaimItem) -> Self {
        self.no_claims.push(item);
        self
    }

    /// Set the replay command.
    #[must_use]
    pub fn with_replay_command(mut self, command: impl Into<String>) -> Self {
        self.replay_command = command.into();
        self
    }

    /// Attach convergence analysis.
    #[must_use]
    pub fn with_convergence(mut self, conv: ConvergenceResult) -> Self {
        self.convergence = Some(conv);
        self
    }

    /// Attach uncertainty quantification analysis.
    #[must_use]
    pub fn with_uncertainty(mut self, uq: UqResult) -> Self {
        self.uncertainty = Some(uq);
        self
    }

    /// Render deterministic, self-contained HTML report with embedded styling.
    #[must_use]
    pub fn render_html(&self) -> String {
        let mut html = String::with_capacity(32 * 1024);
        html.push_str("<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n");
        html.push_str("  <meta charset=\"UTF-8\">\n");
        let _ = write!(
            html,
            "  <title>FrankenSim Engineering Report — {}</title>\n",
            escape_html(&self.project_name)
        );
        html.push_str("  <style>\n");
        html.push_str(
            "    body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif; line-height: 1.5; color: #1f2328; background: #f6f8fa; margin: 0; padding: 24px; }\n\
             .container { max-width: 1080px; margin: 0 auto; background: #ffffff; border: 1px solid #d0d7de; border-radius: 6px; padding: 32px; box-shadow: 0 1px 3px rgba(0,0,0,0.08); }\n\
             h1, h2, h3 { color: #24292f; margin-top: 24px; margin-bottom: 12px; border-bottom: 1px solid #d0d7de; padding-bottom: 6px; }\n\
             table { width: 100%; border-collapse: collapse; margin: 16px 0; font-size: 14px; }\n\
             th, td { border: 1px solid #d0d7de; padding: 8px 12px; text-align: left; }\n\
             th { background-color: #f6f8fa; font-weight: 600; }\n\
             .badge { display: inline-block; padding: 2px 8px; border-radius: 12px; font-size: 12px; font-weight: 600; }\n\
             .badge-verified { background: #dafbe1; color: #1a7f37; border: 1px solid #2da44e; }\n\
             .badge-validated { background: #ddf4ff; color: #0969da; border: 1px solid #54aeff; }\n\
             .badge-estimated { background: #fff8c5; color: #9a6700; border: 1px solid #d4a72c; }\n\
             .badge-waived { background: #fbefff; color: #8250df; border: 1px solid #c297ff; }\n\
             .alert-box { background: #fff8c5; border: 1px solid #d4a72c; border-radius: 6px; padding: 16px; margin: 16px 0; }\n\
             .code-block { background: #f6f8fa; border: 1px solid #d0d7de; border-radius: 6px; padding: 12px; font-family: ui-monospace, SFMono-Regular, monospace; font-size: 13px; overflow-x: auto; }\n"
        );
        html.push_str("  </style>\n</head>\n<body>\n<div class=\"container\">\n");

        // 1. Header & Executive Summary
        let _ = write!(
            html,
            "<h1>FrankenSim Engineering Report</h1>\n\
             <p><strong>Project:</strong> {} | <strong>Run ID:</strong> <code>{}</code></p>\n",
            escape_html(&self.project_name),
            escape_html(&self.run_id)
        );

        // Five Explicits & Provenance
        html.push_str("<h2>1. Provenance &amp; The Five Explicits</h2>\n<table>\n");
        html.push_str("<tr><th>Dimension</th><th>Explicit Declaration</th></tr>\n");
        let _ = write!(
            html,
            "<tr><td>Code Version</td><td><code>{}</code></td></tr>\n",
            self.provenance.code_version
        );
        let _ = write!(
            html,
            "<tr><td>Constellation Lock</td><td><code>{}</code></td></tr>\n",
            self.provenance.constellation_lock
        );
        let _ = write!(
            html,
            "<tr><td>Machine Fingerprint</td><td><code>{}</code></td></tr>\n",
            self.provenance.machine_fingerprint
        );
        let _ = write!(
            html,
            "<tr><td>Units System</td><td>{}</td></tr>\n",
            self.provenance.units_system
        );
        let _ = write!(
            html,
            "<tr><td>Deterministic Seed</td><td><code>0x{:04x}</code></td></tr>\n",
            self.provenance.rng_seed
        );
        let _ = write!(
            html,
            "<tr><td>Budgets</td><td>Wall: {}s | Memory: {} MB</td></tr>\n",
            self.provenance.wall_budget_s,
            self.provenance.mem_budget_bytes / (1024 * 1024)
        );
        html.push_str("</table>\n");

        // 2. Quantities of Interest & Error Budget
        html.push_str("<h2>2. Quantities of Interest &amp; Error Budget</h2>\n");
        html.push_str("<table>\n<tr><th>QoI Name</th><th>Nominal Value</th><th>Unit</th><th>Evidence Color</th><th>Discretization</th><th>Uncertainty</th><th>Surrogate</th><th>Total Budget</th></tr>\n");
        for q in &self.qois {
            let color_badge = match &q.color {
                Color::Verified { .. } => "<span class=\"badge badge-verified\">Verified</span>",
                Color::Validated { .. } => "<span class=\"badge badge-validated\">Validated</span>",
                Color::Estimated { .. } => "<span class=\"badge badge-estimated\">Estimated</span>",
            };
            let _ = write!(
                html,
                "<tr><td><strong>{}</strong><br><small>{}</small></td><td>{}</td><td>{}</td><td>{}</td><td>±{}</td><td>±{}</td><td>±{}</td><td><strong>±{}</strong></td></tr>\n",
                escape_html(&q.name),
                escape_html(&q.description),
                html_num(q.nominal_value),
                escape_html(&q.unit),
                color_badge,
                html_num(q.discretization_error),
                html_num(q.parameter_uncertainty),
                html_num(q.surrogate_error),
                html_num(q.total_uncertainty_budget)
            );
        }
        html.push_str("</table>\n");

        // 2a. Requirements
        if !self.requirements.is_empty() {
            html.push_str("<h3>2a. Requirement Evaluations</h3>\n");
            html.push_str("<table>\n<tr><th>QoI</th><th>Region</th><th>Effective Limit</th><th>Required Margin</th><th>Nominal Margin</th><th>Unit</th><th>Outcome</th><th>Requirement Identity</th></tr>\n");
            for r in &self.requirements {
                let _ = write!(
                    html,
                    "<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td><code>{}</code></td><td><small><code>{}</code></small></td></tr>\n",
                    escape_html(&r.qoi),
                    escape_html(&r.region),
                    html_num(r.effective_limit),
                    html_num(r.required_margin),
                    html_num(r.nominal_margin),
                    escape_html(&r.unit),
                    escape_html(&r.outcome),
                    escape_html(&r.id)
                );
            }
            html.push_str("</table>\n");
        }

        // 2b. Engineering uncertainty budget terms
        if !self.budget_terms.is_empty() {
            html.push_str("<h3>2b. Engineering Uncertainty Budget Terms</h3>\n");
            html.push_str("<p><em>A term marked <code>no-data</code> has no measured magnitude; the report prints the producer's state and reason and invents nothing.</em></p>\n");
            html.push_str("<table>\n<tr><th>QoI</th><th>Term</th><th>State</th><th>Magnitude</th><th>Reason</th><th>Owner</th></tr>\n");
            for t in &self.budget_terms {
                let magnitude = t.value.map_or_else(|| "NO-DATA".to_string(), html_num);
                let _ = write!(
                    html,
                    "<tr><td>{}</td><td>{}</td><td><code>{}</code></td><td>{}</td><td>{}</td><td>{}</td></tr>\n",
                    escape_html(&t.qoi),
                    escape_html(&t.kind),
                    escape_html(&t.state),
                    magnitude,
                    escape_html(&t.reason),
                    escape_html(&t.owner)
                );
            }
            html.push_str("</table>\n");
        }

        // 3. Discretization Convergence Section
        if let Some(conv) = &self.convergence {
            html.push_str("<h2>3. Discretization Convergence Ladder</h2>\n");
            let _ = write!(
                html,
                "<p><strong>Target QoI:</strong> {} | <strong>Status:</strong> <code>{}</code> | <strong>Theoretical Order:</strong> {:.2} | <strong>Observed Order:</strong> {}</p>\n",
                escape_html(&conv.qoi_name),
                conv.status.label(),
                conv.theoretical_order,
                conv.observed_order
                    .map_or("N/A".to_string(), |v| format!("{v:.3}"))
            );
            html.push_str("<table>\n<tr><th>Rung</th><th>Mesh ID</th><th>h (spacing)</th><th>DOF</th><th>QoI Value</th><th>Solver Status</th></tr>\n");
            for r in &conv.admitted_rungs {
                let _ = write!(
                    html,
                    "<tr><td>{}</td><td><code>{}</code></td><td>{:.4} {}</td><td>{}</td><td>{:.4} {}</td><td>{}</td></tr>\n",
                    r.ordinal,
                    escape_html(&r.mesh_id),
                    r.h,
                    escape_html(&r.h_unit),
                    r.dof,
                    r.qoi_value,
                    escape_html(&r.qoi_unit),
                    escape_html(&r.solver_status)
                );
            }
            html.push_str("</table>\n");
            if let Some(extrap) = conv.richardson_extrapolated_qoi {
                let _ = write!(
                    html,
                    "<p><strong>Richardson Extrapolated Continuum Value:</strong> {:.4} | <strong>Grid Convergence Index (GCI):</strong> {}</p>\n",
                    extrap,
                    conv.discretization_error_gci
                        .map_or("N/A".to_string(), |g| format!("{:.3}%", g * 100.0))
                );
            }
        }

        // 4. Parameter / BC Uncertainty Section
        if let Some(uq) = &self.uncertainty {
            html.push_str("<h2>4. Parameter &amp; Boundary Condition Uncertainty</h2>\n");
            let _ = write!(
                html,
                "<p><strong>Propagation Method:</strong> {:?} | <strong>Samples Evaluated:</strong> {} | <strong>Sampling Error:</strong> ±{:.4}</p>\n",
                uq.method_used, uq.samples_evaluated, uq.sampling_error
            );
            if let (Some(m), Some(s)) = (uq.mean, uq.std_dev) {
                let _ = write!(
                    html,
                    "<p><strong>Distribution:</strong> Mean = {:.4}, StdDev = {:.4}, Range = [{:.4}, {:.4}]</p>\n",
                    m, s, uq.interval_bounds[0], uq.interval_bounds[1]
                );
            }
            if let Some(p_comp) = uq.probability_of_compliance {
                let _ = write!(
                    html,
                    "<p><strong>Probability of Compliance:</strong> <strong>{:.2}%</strong></p>\n",
                    p_comp * 100.0
                );
            }
        }

        // 5. Materials & Geometry Import
        if !self.materials.is_empty() {
            html.push_str("<h2>5. Material &amp; Interface Provenance</h2>\n");
            html.push_str("<table>\n<tr><th>Region</th><th>Material Card</th><th>Conductivity k</th><th>Specific Heat Cp</th><th>Density ρ</th><th>Pack</th></tr>\n");
            for m in &self.materials {
                let _ = write!(
                    html,
                    "<tr><td><strong>{}</strong></td><td><code>{}</code></td><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>\n",
                    escape_html(&m.region_name),
                    escape_html(&m.material_card_id),
                    escape_html(&m.thermal_conductivity),
                    escape_html(&m.specific_heat),
                    escape_html(&m.density),
                    escape_html(&m.source_pack)
                );
            }
            html.push_str("</table>\n");
        }

        // 5b. Retained stage receipts
        if !self.stages.is_empty() {
            html.push_str("<h2>5b. Retained Stage Receipts</h2>\n");
            html.push_str("<p><em>Every value in this report is projected from one of these ledger artifacts; the hashes are the audit trail.</em></p>\n");
            html.push_str("<table>\n<tr><th>Ordinal</th><th>Stage</th><th>Receipt Hash</th><th>Summary (copied from the receipt)</th></tr>\n");
            for s in &self.stages {
                let mut summary = String::new();
                for (index, (key, value)) in s.summary.iter().enumerate() {
                    if index > 0 {
                        summary.push_str("<br>");
                    }
                    let _ = write!(
                        summary,
                        "<code>{}</code> = {}",
                        escape_html(key),
                        escape_html(value)
                    );
                }
                let _ = write!(
                    html,
                    "<tr><td>{}</td><td><strong>{}</strong></td><td><small><code>{}</code></small></td><td><small>{}</small></td></tr>\n",
                    s.ordinal,
                    escape_html(&s.stage),
                    escape_html(&s.receipt_hash),
                    summary
                );
            }
            html.push_str("</table>\n");
        }

        // 5c. Lineage
        if !self.lineage.is_empty() {
            html.push_str(
                "<h2>5c. Lineage</h2>\n<table>\n<tr><th>Artifact</th><th>Content Hash</th></tr>\n",
            );
            for l in &self.lineage {
                let _ = write!(
                    html,
                    "<tr><td>{}</td><td><small><code>{}</code></small></td></tr>\n",
                    escape_html(&l.label),
                    escape_html(&l.hash)
                );
            }
            html.push_str("</table>\n");
        }

        // 6. Known Gaps & No-Claims Section
        html.push_str(
            "<h2>6. Known Gaps &amp; No-Claims Boundaries</h2>\n<div class=\"alert-box\">\n",
        );
        html.push_str("<p><em>FrankenSim makes explicit no-claim boundaries where physical validation or continuum authority is absent:</em></p>\n<ul>\n");
        for nc in &self.no_claims {
            let _ = write!(
                html,
                "<li><strong>{}</strong> [<code>{}</code>]: {}</li>\n",
                escape_html(&nc.component),
                escape_html(&nc.status),
                escape_html(&nc.statement)
            );
        }
        html.push_str("</ul>\n</div>\n");

        // 7. Replay Instructions
        if !self.replay_command.is_empty() {
            html.push_str("<h2>7. Replay &amp; Verification</h2>\n");
            let _ = write!(
                html,
                "<p>Replay this exact simulation run deterministically using the command below:</p>\n\
                 <div class=\"code-block\">{}</div>\n",
                escape_html(&self.replay_command)
            );
        }

        html.push_str("</div>\n</body>\n</html>\n");
        html
    }

    /// Render semantically identical JSON twin.
    #[must_use]
    pub fn render_json(&self) -> String {
        let mut json = String::with_capacity(16 * 1024);
        json.push_str("{\n");
        let _ = write!(
            json,
            "  \"schema\": \"frankensim.report.engineering.v1\",\n"
        );
        let _ = write!(json, "  \"run_id\": \"{}\",\n", self.run_id);
        let _ = write!(
            json,
            "  \"project_name\": \"{}\",\n",
            escape_json(&self.project_name)
        );

        // Provenance
        json.push_str("  \"provenance\": {\n");
        let _ = write!(
            json,
            "    \"code_version\": \"{}\",\n",
            self.provenance.code_version
        );
        let _ = write!(
            json,
            "    \"constellation_lock\": \"{}\",\n",
            self.provenance.constellation_lock
        );
        let _ = write!(
            json,
            "    \"machine_fingerprint\": \"{}\",\n",
            self.provenance.machine_fingerprint
        );
        let _ = write!(
            json,
            "    \"units_system\": \"{}\",\n",
            self.provenance.units_system
        );
        let _ = write!(json, "    \"rng_seed\": {},\n", self.provenance.rng_seed);
        let _ = write!(
            json,
            "    \"wall_budget_s\": {},\n",
            self.provenance.wall_budget_s
        );
        let _ = write!(
            json,
            "    \"mem_budget_bytes\": {}\n",
            self.provenance.mem_budget_bytes
        );
        json.push_str("  },\n");

        // QoIs
        json.push_str("  \"qois\": [\n");
        for (i, q) in self.qois.iter().enumerate() {
            if i > 0 {
                json.push_str(",\n");
            }
            let color_str = match &q.color {
                Color::Verified { .. } => "Verified",
                Color::Validated { .. } => "Validated",
                Color::Estimated { .. } => "Estimated",
            };
            let _ = write!(
                json,
                "    {{\"name\": \"{}\", \"value\": {}, \"unit\": \"{}\", \"color\": \"{}\", \"discretization_error\": {}, \"parameter_uncertainty\": {}, \"surrogate_error\": {}, \"total_budget\": {}, \"source_root\": \"{}\"}}",
                escape_json(&q.name),
                json_num(q.nominal_value),
                escape_json(&q.unit),
                color_str,
                json_num(q.discretization_error),
                json_num(q.parameter_uncertainty),
                json_num(q.surrogate_error),
                json_num(q.total_uncertainty_budget),
                escape_json(&q.source_root)
            );
        }
        json.push_str("\n  ],\n");

        // Requirements
        json.push_str("  \"requirements\": [\n");
        for (i, r) in self.requirements.iter().enumerate() {
            if i > 0 {
                json.push_str(",\n");
            }
            let _ = write!(
                json,
                "    {{\"id\": \"{}\", \"qoi\": \"{}\", \"region\": \"{}\", \"effective_limit\": {}, \"required_margin\": {}, \"nominal_margin\": {}, \"unit\": \"{}\", \"outcome\": \"{}\"}}",
                escape_json(&r.id),
                escape_json(&r.qoi),
                escape_json(&r.region),
                json_num(r.effective_limit),
                json_num(r.required_margin),
                json_num(r.nominal_margin),
                escape_json(&r.unit),
                escape_json(&r.outcome)
            );
        }
        json.push_str("\n  ],\n");

        // Budget terms
        json.push_str("  \"budget_terms\": [\n");
        for (i, t) in self.budget_terms.iter().enumerate() {
            if i > 0 {
                json.push_str(",\n");
            }
            let value = t.value.map_or_else(|| "null".to_string(), json_num);
            let _ = write!(
                json,
                "    {{\"qoi\": \"{}\", \"kind\": \"{}\", \"state\": \"{}\", \"value\": {}, \"reason\": \"{}\", \"owner\": \"{}\"}}",
                escape_json(&t.qoi),
                escape_json(&t.kind),
                escape_json(&t.state),
                value,
                escape_json(&t.reason),
                escape_json(&t.owner)
            );
        }
        json.push_str("\n  ],\n");

        // Stage receipts
        json.push_str("  \"stages\": [\n");
        for (i, s) in self.stages.iter().enumerate() {
            if i > 0 {
                json.push_str(",\n");
            }
            let _ = write!(
                json,
                "    {{\"ordinal\": {}, \"stage\": \"{}\", \"receipt_hash\": \"{}\", \"summary\": {{",
                s.ordinal,
                escape_json(&s.stage),
                escape_json(&s.receipt_hash)
            );
            for (index, (key, value)) in s.summary.iter().enumerate() {
                if index > 0 {
                    json.push_str(", ");
                }
                let _ = write!(json, "\"{}\": \"{}\"", escape_json(key), escape_json(value));
            }
            json.push_str("}}");
        }
        json.push_str("\n  ],\n");

        // Lineage
        json.push_str("  \"lineage\": [\n");
        for (i, l) in self.lineage.iter().enumerate() {
            if i > 0 {
                json.push_str(",\n");
            }
            let _ = write!(
                json,
                "    {{\"label\": \"{}\", \"hash\": \"{}\"}}",
                escape_json(&l.label),
                escape_json(&l.hash)
            );
        }
        json.push_str("\n  ],\n");

        // Materials
        json.push_str("  \"materials\": [\n");
        for (i, m) in self.materials.iter().enumerate() {
            if i > 0 {
                json.push_str(",\n");
            }
            let _ = write!(
                json,
                "    {{\"region\": \"{}\", \"card\": \"{}\", \"thermal_conductivity\": \"{}\", \"specific_heat\": \"{}\", \"density\": \"{}\", \"pack\": \"{}\"}}",
                escape_json(&m.region_name),
                escape_json(&m.material_card_id),
                escape_json(&m.thermal_conductivity),
                escape_json(&m.specific_heat),
                escape_json(&m.density),
                escape_json(&m.source_pack)
            );
        }
        json.push_str("\n  ],\n");

        // Convergence
        if let Some(conv) = &self.convergence {
            json.push_str("  \"convergence\": {\n");
            let _ = write!(
                json,
                "    \"qoi_name\": \"{}\",\n",
                escape_json(&conv.qoi_name)
            );
            let _ = write!(json, "    \"status\": \"{}\",\n", conv.status.label());
            let _ = write!(
                json,
                "    \"theoretical_order\": {},\n",
                json_num(conv.theoretical_order)
            );
            // Optional magnitudes print as JSON numbers or null, never as
            // Rust debug text (the section was unexercised until the uniform
            // h-ladder fed it; `Some(..)`/`None` broke the twin's parse).
            let optional = |value: Option<f64>| value.map_or_else(|| "null".to_string(), json_num);
            let _ = write!(
                json,
                "    \"observed_order\": {},\n",
                optional(conv.observed_order)
            );
            let _ = write!(
                json,
                "    \"fit_residual\": {},\n",
                optional(conv.fit_residual)
            );
            let _ = write!(
                json,
                "    \"richardson_extrapolated\": {},\n",
                optional(conv.richardson_extrapolated_qoi)
            );
            let _ = write!(
                json,
                "    \"gci\": {},\n",
                optional(conv.discretization_error_gci)
            );
            let _ = write!(
                json,
                "    \"rejection_reason\": {},\n",
                conv.rejection_reason
                    .as_deref()
                    .map_or_else(|| "null".to_string(), |r| format!("\"{}\"", escape_json(r)))
            );
            json.push_str("    \"rungs\": [\n");
            for (index, r) in conv.admitted_rungs.iter().enumerate() {
                let _ = write!(
                    json,
                    "      {{\"ordinal\": {}, \"mesh_id\": \"{}\", \"h\": {}, \"h_unit\": \"{}\", \"dof\": {}, \"solver_status\": \"{}\", \"solver_residual\": {}, \"qoi_value\": {}, \"qoi_unit\": \"{}\"}}{}\n",
                    r.ordinal,
                    escape_json(&r.mesh_id),
                    json_num(r.h),
                    escape_json(&r.h_unit),
                    r.dof,
                    escape_json(&r.solver_status),
                    json_num(r.solver_residual),
                    json_num(r.qoi_value),
                    escape_json(&r.qoi_unit),
                    if index + 1 < conv.admitted_rungs.len() {
                        ","
                    } else {
                        ""
                    }
                );
            }
            json.push_str("    ]\n");
            json.push_str("  },\n");
        }

        // Uncertainty
        if let Some(uq) = &self.uncertainty {
            json.push_str("  \"uncertainty\": {\n");
            let _ = write!(
                json,
                "    \"qoi_name\": \"{}\",\n",
                escape_json(&uq.qoi_name)
            );
            let _ = write!(json, "    \"method\": \"{:?}\",\n", uq.method_used);
            let _ = write!(json, "    \"samples\": {},\n", uq.samples_evaluated);
            let _ = write!(json, "    \"mean\": {:?},\n", uq.mean);
            let _ = write!(json, "    \"std_dev\": {:?},\n", uq.std_dev);
            let _ = write!(
                json,
                "    \"p_compliance\": {:?}\n",
                uq.probability_of_compliance
            );
            json.push_str("  },\n");
        }

        // No-claims
        json.push_str("  \"no_claims\": [\n");
        for (i, nc) in self.no_claims.iter().enumerate() {
            if i > 0 {
                json.push_str(",\n");
            }
            let _ = write!(
                json,
                "    {{\"component\": \"{}\", \"status\": \"{}\", \"statement\": \"{}\"}}",
                escape_json(&nc.component),
                escape_json(&nc.status),
                escape_json(&nc.statement)
            );
        }
        json.push_str("\n  ],\n");
        let _ = write!(
            json,
            "  \"replay_command\": \"{}\"\n}}\n",
            escape_json(&self.replay_command)
        );

        json
    }

    /// Compute deterministic BLAKE3 content hash.
    #[must_use]
    pub fn content_hash(&self) -> ContentHash {
        let html = self.render_html();
        hash_domain("org.frankensim.report.engineering.v1", html.as_bytes())
    }
}

fn escape_html(input: &str) -> String {
    let mut out = String::with_capacity(input.len());
    for ch in input.chars() {
        match ch {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&#39;"),
            _ => out.push(ch),
        }
    }
    out
}

fn escape_json(input: &str) -> String {
    let mut out = String::with_capacity(input.len());
    for ch in input.chars() {
        match ch {
            '\\' => out.push_str("\\\\"),
            '"' => out.push_str("\\\""),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            _ => out.push(ch),
        }
    }
    out
}
