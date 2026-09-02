//! Seventh solve stage: project the retained stage receipts into a
//! deterministic engineering report, its JSON twin, and a checkable evidence
//! package, then retain all three as ledger artifacts.
//!
//! Authority doctrine (bead frankensim-rc-root-q61wp.12): this stage adds NO
//! physical, numerical, or validation authority. Every number in the report
//! is copied from a retained receipt that names its own producer; every claim
//! in the package carries the colour its producer recorded (today: Estimated
//! with an unbounded dispersion, because all eight engineering-uncertainty
//! terms are explicit NO-DATA). A receipt that does not parse, or that lacks a
//! field the report needs, refuses the stage instead of being papered over
//! with a literal — the 2026-08-25 fabricated report is the failure this
//! module exists to make impossible.
//!
//! The same module hosts the read path the `report` / `package` export verbs
//! use: they locate a completed run through the resume loader and hand back
//! the exact retained bytes, so an export can never differ from what the
//! solve stage sealed.

use fs_blake3::{ContentHash, hash_bytes};
use fs_checker::check as check_package;
use fs_evidence::Color;
use fs_exec::CancelGate;
use fs_ledger::{EdgeRole, Ledger, OpArtifactEdge};
use fs_package::{Claim, EvidencePackage, Provenance};
use fs_project::spec::ProjectSpec;
use fs_report::{
    BudgetTermItem, EngineeringReport, LineageItem, MaterialReportItem, NoClaimItem, QoiReportItem,
    ReportProvenance, RequirementReportItem, StageReceiptItem,
};

use super::{
    CompletedStage, EDGE_SCAN_CAP, EvidenceWork, InvocationWorkLedger, QOI_RECEIPT_SCHEMA,
    ResumeProof, RetainedSideArtifact, SolveDriverState, SolveRefusal, SolveRunId, SolveStage,
    VerifiedResume, invocation_work_refusal, load_latest_state, require_exact_edges,
};
use crate::import::json_string;
use crate::json_read::JsonValue;

/// Ledger artifact kind of the rendered HTML report.
pub(super) const REPORT_HTML_KIND: &str = "solve-report-html";
/// Ledger artifact kind of the JSON twin.
pub(super) const REPORT_JSON_KIND: &str = "solve-report-json";
/// Ledger artifact kind of the evidence package (format-9 JSON).
pub(super) const EVIDENCE_PACKAGE_KIND: &str = "solve-evidence-package";
/// Stage receipt schema.
pub(super) const REPORT_RECEIPT_SCHEMA: &str = "frankensim.cli.solve-report.v1";

/// Product of the report stage: the receipt plus the three retained artifacts.
#[derive(Debug)]
pub(super) struct ReportStageProduct {
    pub(super) receipt: String,
    pub(super) artifacts: Vec<RetainedSideArtifact>,
}

fn report_error(
    code: &'static str,
    what: impl Into<String>,
    fix: impl Into<String>,
) -> SolveRefusal {
    SolveRefusal::staged(code, SolveStage::Report, what, fix)
}

fn shape_error(stage: SolveStage, what: impl Into<String>) -> SolveRefusal {
    report_error(
        "cli-solve-report-receipt-shape",
        format!(
            "the retained {} receipt does not carry what the report needs: {}",
            stage.name(),
            what.into()
        ),
        "regenerate the run with the current producers; the report never substitutes a value for a missing field",
    )
}

/// A parsed stage receipt with its ledger coordinates.
struct LoadedReceipt {
    stage: SolveStage,
    completed: CompletedStage,
    value: JsonValue,
}

fn load_receipt(
    ledger: &Ledger,
    stage: SolveStage,
    completed: &CompletedStage,
    work: &EvidenceWork<'_>,
    run: SolveRunId,
) -> Result<LoadedReceipt, SolveRefusal> {
    let bytes = ledger
        .get_artifact(&completed.receipt)
        .map_err(|error| {
            report_error(
                "cli-solve-report-ledger",
                format!("reading the {} receipt failed: {error}", stage.name()),
                "repair or restore the ledger; the report reads only retained artifacts",
            )
        })?
        .ok_or_else(|| {
            report_error(
                "cli-solve-report-receipt-missing",
                format!(
                    "the {} receipt {} is not retained in this ledger",
                    stage.name(),
                    completed.receipt.to_hex()
                ),
                "resume the run in the ledger that retained its receipts",
            )
        })?;
    let charged = u64::try_from(bytes.len()).unwrap_or(u64::MAX);
    work.charge(charged)
        .map_err(|error| invocation_work_refusal(Some(run), Some(SolveStage::Report), error))?;
    let text =
        String::from_utf8(bytes).map_err(|_| shape_error(stage, "receipt bytes are not UTF-8"))?;
    let value = JsonValue::parse(&text).map_err(|error| shape_error(stage, error.to_string()))?;
    Ok(LoadedReceipt {
        stage,
        completed: *completed,
        value,
    })
}

fn required_str<'a>(
    stage: SolveStage,
    node: &'a JsonValue,
    key: &str,
) -> Result<&'a str, SolveRefusal> {
    node.str_field(key)
        .ok_or_else(|| shape_error(stage, format!("missing string field `{key}`")))
}

fn required_f64(stage: SolveStage, node: &JsonValue, key: &str) -> Result<f64, SolveRefusal> {
    node.f64_field(key)
        .ok_or_else(|| shape_error(stage, format!("missing numeric field `{key}`")))
}

fn required_first<'a>(
    stage: SolveStage,
    node: &'a JsonValue,
    key: &str,
) -> Result<&'a JsonValue, SolveRefusal> {
    node.get(key)
        .and_then(JsonValue::as_array)
        .and_then(<[JsonValue]>::first)
        .ok_or_else(|| shape_error(stage, format!("missing non-empty array `{key}`")))
}

/// Render one scalar receipt field for a summary row without inventing it:
/// numbers keep their exact spelling, strings are copied, anything else is
/// skipped by the caller.
fn scalar_text(node: &JsonValue) -> Option<String> {
    match node {
        JsonValue::Str(s) => Some(s.clone()),
        JsonValue::Number { raw, .. } => Some(raw.clone()),
        JsonValue::Bool(b) => Some(b.to_string()),
        JsonValue::Null | JsonValue::Array(_) | JsonValue::Object(_) => None,
    }
}

fn push_summary(
    summary: &mut Vec<(String, String)>,
    receipt: &JsonValue,
    label: &str,
    path: &[&str],
) {
    if let Some(text) = receipt.path(path).and_then(scalar_text) {
        summary.push((label.to_string(), text));
    }
}

fn stage_summary(stage: SolveStage, receipt: &JsonValue) -> Vec<(String, String)> {
    let mut summary = Vec::new();
    match stage {
        SolveStage::ImportVerify => {
            if let Some(entries) = receipt.get("entries").and_then(JsonValue::as_array) {
                summary.push(("verified_imports".to_string(), entries.len().to_string()));
            }
            push_summary(&mut summary, receipt, "authority", &["authority"]);
        }
        SolveStage::Assign => {
            if let Some(regions) = receipt.get("regions").and_then(JsonValue::as_array) {
                summary.push(("regions".to_string(), regions.len().to_string()));
            }
        }
        SolveStage::MaterialResolve => {
            if let Some(bindings) = receipt.get("bindings").and_then(JsonValue::as_array) {
                summary.push(("bindings".to_string(), bindings.len().to_string()));
            }
            push_summary(&mut summary, receipt, "pack_set_root", &["pack_set_root"]);
        }
        SolveStage::FlowNetwork => {
            for (label, path) in [
                ("flow_lo", &["operating_point", "flow_lo"][..]),
                ("flow_mid", &["operating_point", "flow_mid"][..]),
                ("flow_hi", &["operating_point", "flow_hi"][..]),
                ("pressure_lo", &["operating_point", "pressure_lo"][..]),
                ("pressure_hi", &["operating_point", "pressure_hi"][..]),
                ("leakage_fraction", &["leakage_fraction"][..]),
                ("vent_count", &["vent_count"][..]),
                ("authority", &["authority"][..]),
            ] {
                push_summary(&mut summary, receipt, label, path);
                if summary
                    .last()
                    .is_none_or(|(last, _)| last.as_str() != label)
                    && path.len() > 1
                {
                    push_summary(&mut summary, receipt, label, &path[1..]);
                }
            }
        }
        SolveStage::Conduction => {
            for (label, path) in [
                ("temperature_min", &["temperature", "min"][..]),
                ("temperature_max", &["temperature", "max"][..]),
                ("temperature_unit", &["temperature_unit"][..]),
                ("iterations", &["solver", "iterations"][..]),
                ("final_residual", &["solver", "final_residual"][..]),
                ("stop_reason", &["solver", "stop_reason"][..]),
                (
                    "energy_relative_closure",
                    &["energy", "relative_closure"][..],
                ),
                ("elements", &["mesh", "elements"][..]),
                ("vertices", &["mesh", "vertices"][..]),
                ("authority", &["authority"][..]),
            ] {
                push_summary(&mut summary, receipt, label, path);
                if summary
                    .last()
                    .is_none_or(|(last, _)| last.as_str() != label)
                    && path.len() > 1
                {
                    push_summary(&mut summary, receipt, label, &path[1..]);
                }
            }
        }
        SolveStage::Qoi => {
            push_summary(&mut summary, receipt, "authority", &["authority"]);
            push_summary(
                &mut summary,
                receipt,
                "composition_identity",
                &["composition_identity"],
            );
        }
        SolveStage::Report => {}
    }
    summary
}

fn material_rows(receipt: &JsonValue) -> Vec<MaterialReportItem> {
    let Some(bindings) = receipt.get("bindings").and_then(JsonValue::as_array) else {
        return Vec::new();
    };
    let pack = receipt
        .str_field("pack_set_root")
        .unwrap_or("not retained in receipt")
        .to_string();
    bindings
        .iter()
        .map(|binding| {
            let region = binding
                .str_field("target")
                .map(str::to_string)
                .or_else(|| binding.path(&["target", "name"]).and_then(scalar_text))
                .unwrap_or_else(|| "not retained in receipt".to_string());
            let card = binding
                .str_field("card")
                .map(str::to_string)
                .or_else(|| binding.str_field("card_identity").map(str::to_string))
                .or_else(|| binding.path(&["card", "identity"]).and_then(scalar_text))
                .unwrap_or_else(|| "not retained in receipt".to_string());
            let mut conductivity = "not retained in receipt".to_string();
            let mut specific_heat = "not retained in receipt".to_string();
            let mut density = "not retained in receipt".to_string();
            if let Some(properties) = binding
                .get("material_properties")
                .and_then(JsonValue::as_array)
            {
                for property in properties {
                    let Some(name) = property.str_field("property") else {
                        continue;
                    };
                    let lo = property.get("value_lo").and_then(scalar_text);
                    let hi = property.get("value_hi").and_then(scalar_text);
                    let dims = property.str_field("dims").unwrap_or("");
                    let text = match (lo, hi) {
                        (Some(lo), Some(hi)) if lo == hi => {
                            format!("{lo} {dims}").trim().to_string()
                        }
                        (Some(lo), Some(hi)) => format!("{lo}..{hi} {dims}").trim().to_string(),
                        (Some(one), None) | (None, Some(one)) => {
                            format!("{one} {dims}").trim().to_string()
                        }
                        (None, None) => continue,
                    };
                    if name.contains("conductivity") {
                        conductivity = text;
                    } else if name.contains("specific_heat") || name.contains("heat_capacity") {
                        specific_heat = text;
                    } else if name.contains("density") {
                        density = text;
                    }
                }
            }
            MaterialReportItem {
                region_name: region,
                material_card_id: card,
                thermal_conductivity: conductivity,
                specific_heat,
                density,
                source_pack: pack.clone(),
            }
        })
        .collect()
}

/// Spell an arbitrary declared name as a package identity token: every byte
/// outside fs-package's identity alphabet becomes `_`, and an empty result
/// becomes `undeclared`. The original spelling stays visible in the claim
/// statement, so nothing is lost — only the machine identity is normalised.
fn identity_token(text: &str) -> String {
    let token: String = text
        .chars()
        .map(|ch| {
            if ch.is_ascii_alphanumeric()
                || matches!(ch, '-' | '_' | '.' | '/' | ':' | '@' | '+' | '=')
            {
                ch
            } else {
                '_'
            }
        })
        .collect();
    if token.is_empty() {
        "undeclared".to_string()
    } else {
        token
    }
}

fn provenance_from_spec(spec: &ProjectSpec) -> ReportProvenance {
    let (constellation, workspace) = spec
        .versions
        .as_ref()
        .map_or(("undeclared", "undeclared"), |v| {
            (v.constellation.as_str(), v.workspace.as_str())
        });
    let wall_budget_s = spec
        .budgets
        .as_ref()
        .map_or(0, |b| b.solve_time.value.max(0.0).floor() as u64);
    let mem_budget_bytes = spec.budgets.as_ref().map_or(0, |b| b.memory_bytes);
    ReportProvenance {
        code_version: format!(
            "fs-cli {}; project workspace {workspace}",
            env!("CARGO_PKG_VERSION")
        ),
        constellation_lock: constellation.to_string(),
        machine_fingerprint: "not retained by the solve run (roofline lanes own machine identity)"
            .to_string(),
        units_system: "coherent SI base values (m, kg, s, K); temperatures in kelvin".to_string(),
        rng_seed: spec.seeds.as_ref().map_or(0, |s| s.root),
        wall_budget_s,
        mem_budget_bytes,
    }
}

/// Build the report, its JSON twin, and the evidence package from the retained
/// receipts of the completed six-stage prefix, then return them as side
/// artifacts plus the stage receipt.
#[allow(clippy::too_many_lines)]
pub(super) fn report_receipt(
    ledger: &Ledger,
    spec: &ProjectSpec,
    state: &SolveDriverState,
    run: SolveRunId,
    project_hash: ContentHash,
    work: EvidenceWork<'_>,
) -> Result<ReportStageProduct, SolveRefusal> {
    if work.is_requested() {
        return Err(SolveRefusal::plain(
            "cli-solve-cancelled",
            "cancellation was requested before the report stage started",
            "resume the run; no report artifact was sealed",
        ));
    }
    let prefix: Vec<SolveStage> = SolveStage::ALL
        .iter()
        .copied()
        .filter(|stage| *stage != SolveStage::Report)
        .collect();
    let mut loaded = Vec::with_capacity(prefix.len());
    for stage in prefix {
        let completed = state
            .completed
            .get(stage.ordinal() as usize)
            .filter(|completed| completed.ordinal == stage.ordinal())
            .ok_or_else(|| {
                report_error(
                    "cli-solve-report-prefix",
                    format!(
                        "the report stage requires the completed {} receipt",
                        stage.name()
                    ),
                    "rerun the ordered producer prefix; report projects only completed stages",
                )
            })?;
        loaded.push(load_receipt(ledger, stage, completed, &work, run)?);
    }
    let qoi = loaded
        .iter()
        .find(|receipt| receipt.stage == SolveStage::Qoi)
        .expect("prefix includes qoi");
    let conduction = loaded
        .iter()
        .find(|receipt| receipt.stage == SolveStage::Conduction)
        .expect("prefix includes conduction");
    let materials = loaded
        .iter()
        .find(|receipt| receipt.stage == SolveStage::MaterialResolve)
        .expect("prefix includes material-resolve");

    // ---- QoI receipt: the only source of every number the report states.
    let qoi_stage = SolveStage::Qoi;
    let qoi_schema = required_str(qoi_stage, &qoi.value, "schema")?;
    if qoi_schema != QOI_RECEIPT_SCHEMA {
        return Err(shape_error(
            qoi_stage,
            format!("schema `{qoi_schema}` is not `{QOI_RECEIPT_SCHEMA}`"),
        ));
    }
    let qoi_row = required_first(qoi_stage, &qoi.value, "qoi")?;
    let qoi_name = required_str(qoi_stage, qoi_row, "name")?;
    let qoi_region = required_str(qoi_stage, qoi_row, "region")?;
    let qoi_value = required_f64(qoi_stage, qoi_row, "value")?;
    let qoi_unit = required_str(qoi_stage, qoi_row, "unit")?;
    let qoi_color = required_str(qoi_stage, qoi_row, "color")?;
    let qoi_identity = required_str(qoi_stage, qoi_row, "identity")?;
    if qoi_color != "estimated" {
        return Err(shape_error(
            qoi_stage,
            format!("QoI colour `{qoi_color}` is not the estimate-only producer's colour"),
        ));
    }
    let requirement = required_first(qoi_stage, &qoi.value, "requirements")?;
    let requirement_id = required_str(qoi_stage, requirement, "id")?;
    let effective_limit = required_f64(qoi_stage, requirement, "effective_limit_kelvin")?;
    let required_margin = required_f64(qoi_stage, requirement, "required_margin_kelvin")?;
    let nominal_margin = required_f64(qoi_stage, requirement, "nominal_margin_kelvin")?;
    let outcome = required_str(qoi_stage, requirement, "outcome")?;
    let budget = required_first(qoi_stage, &qoi.value, "budget")?;
    let terms = budget
        .get("terms")
        .and_then(JsonValue::as_array)
        .ok_or_else(|| shape_error(qoi_stage, "missing `budget[0].terms`"))?;
    let mut budget_items = Vec::with_capacity(terms.len());
    let mut measured_terms = 0usize;
    for term in terms {
        let kind = required_str(qoi_stage, term, "kind")?;
        let term_state = required_str(qoi_stage, term, "state")?;
        let value = term.f64_field("value");
        if term_state != "no-data" && value.is_none() {
            return Err(shape_error(
                qoi_stage,
                format!("term `{kind}` is `{term_state}` but carries no magnitude"),
            ));
        }
        if value.is_some() {
            measured_terms += 1;
        }
        budget_items.push(BudgetTermItem {
            qoi: qoi_name.to_string(),
            kind: kind.to_string(),
            state: term_state.to_string(),
            value,
            reason: term.str_field("reason").unwrap_or("").to_string(),
            owner: term.str_field("owner").unwrap_or("").to_string(),
        });
    }
    let lineage = qoi
        .value
        .get("lineage")
        .ok_or_else(|| shape_error(qoi_stage, "missing `lineage`"))?;
    let lineage_project = required_str(qoi_stage, lineage, "project")?;
    if lineage_project != project_hash.to_hex() {
        return Err(shape_error(
            qoi_stage,
            format!(
                "lineage project `{lineage_project}` is not this run's project `{}`",
                project_hash.to_hex()
            ),
        ));
    }
    let lineage_conduction_receipt = required_str(qoi_stage, lineage, "conduction_receipt")?;
    if lineage_conduction_receipt != conduction.completed.receipt.to_hex() {
        return Err(shape_error(
            qoi_stage,
            "lineage conduction receipt does not match the completed conduction stage",
        ));
    }
    let lineage_solution = required_str(qoi_stage, lineage, "conduction_solution")?;
    let composition_identity = required_str(qoi_stage, &qoi.value, "composition_identity")?;

    // ---- Report structure.
    let project_name = spec
        .metadata
        .as_ref()
        .map_or("(unnamed project)", |m| m.name.as_str());
    let run_hex = run.to_hex();
    let mut report = EngineeringReport::new(run_hex.clone(), project_name)
        .with_provenance(provenance_from_spec(spec))
        .with_qoi(QoiReportItem {
            name: qoi_name.to_string(),
            description: format!(
                "region `{qoi_region}` maximum temperature from the request-selective conduction producer"
            ),
            nominal_value: qoi_value,
            unit: qoi_unit.to_string(),
            color: Color::Estimated {
                estimator: QOI_RECEIPT_SCHEMA.to_string(),
                dispersion: f64::NAN,
            },
            discretization_error: f64::NAN,
            parameter_uncertainty: f64::NAN,
            surrogate_error: f64::NAN,
            total_uncertainty_budget: f64::NAN,
            source_root: qoi_identity.to_string(),
        })
        .with_requirement(RequirementReportItem {
            id: requirement_id.to_string(),
            qoi: qoi_name.to_string(),
            region: qoi_region.to_string(),
            effective_limit,
            required_margin,
            nominal_margin,
            unit: "kelvin".to_string(),
            outcome: outcome.to_string(),
        });
    for item in budget_items {
        report = report.with_budget_term(item);
    }
    for receipt in &loaded {
        report = report.with_stage(StageReceiptItem {
            stage: receipt.stage.name().to_string(),
            ordinal: receipt.completed.ordinal,
            receipt_hash: receipt.completed.receipt.to_hex(),
            summary: stage_summary(receipt.stage, &receipt.value),
        });
        if let Some(no_claim) = receipt.value.str_field("no_claim") {
            report = report.with_no_claim(NoClaimItem {
                component: receipt.stage.name().to_string(),
                status: receipt
                    .value
                    .str_field("authority")
                    .unwrap_or("retained")
                    .to_string(),
                statement: no_claim.to_string(),
            });
        }
    }
    for material in material_rows(&materials.value) {
        report = report.with_material(material);
    }
    for (label, hash) in [
        ("project canonical source", project_hash.to_hex()),
        ("conduction solution field", lineage_solution.to_string()),
        ("conduction receipt", lineage_conduction_receipt.to_string()),
        ("QoI receipt", qoi.completed.receipt.to_hex()),
        ("requirement composition", composition_identity.to_string()),
    ] {
        report = report.with_lineage(LineageItem {
            label: label.to_string(),
            hash,
        });
    }
    report = report
        .with_no_claim(NoClaimItem {
            component: "report".to_string(),
            status: "projection".to_string(),
            statement: "this report projects retained stage receipts and adds no physical, numerical, or validation authority; the QoI is an Estimated candidate whose eight uncertainty terms are NO-DATA, so no binary compliance verdict, DWR bound, corpus validation, or L3/L4 maturity is claimed".to_string(),
        })
        .with_replay_command(format!(
            "frankensim solve <project.fsim> <ledger> --materials <pack>   # project canonical hash {}; run {run_hex}",
            project_hash.to_hex()
        ));
    let html = report.render_html();
    let json = report.render_json();
    let content_hash = report.content_hash();

    // ---- Evidence package: one Estimated claim per retained statement.
    // Package identities admit only `[A-Za-z0-9-_./:@+=]` (fs-package's
    // identity policy), so provenance and claim ids are spelled as tokens.
    let (constellation, workspace) = spec
        .versions
        .as_ref()
        .map_or(("undeclared", "undeclared"), |v| {
            (v.constellation.as_str(), v.workspace.as_str())
        });
    let package = EvidencePackage::new(Provenance::new(
        format!(
            "fs-cli@{}+workspace={}",
            env!("CARGO_PKG_VERSION"),
            identity_token(workspace)
        ),
        identity_token(constellation),
    ))
    .with_claim(Claim::estimated(
        format!("qoi.{}.{}", identity_token(qoi_name), identity_token(qoi_region)),
        format!(
            "{qoi_name} in region `{qoi_region}` = {qoi_value} {qoi_unit}; estimate-only candidate from {QOI_RECEIPT_SCHEMA} with all {} engineering-uncertainty terms NO-DATA (receipt {})",
            terms.len(),
            qoi.completed.receipt.to_hex()
        ),
        QOI_RECEIPT_SCHEMA,
        f64::INFINITY,
    ))
    .with_claim(Claim::estimated(
        format!(
            "requirement.{}.{}",
            identity_token(qoi_name),
            identity_token(qoi_region)
        ),
        format!(
            "requirement outcome `{outcome}`: effective limit {effective_limit} K, required margin {required_margin} K, nominal margin {nominal_margin} K (composition {composition_identity})"
        ),
        QOI_RECEIPT_SCHEMA,
        f64::INFINITY,
    ));
    let package_root = package.try_merkle_root().map_err(|error| {
        report_error(
            "cli-solve-report-package",
            format!("the evidence package did not seal: {error}"),
            "report the packaging defect; no partial package is retained",
        )
    })?;
    let package_json = package.to_json().map_err(|error| {
        report_error(
            "cli-solve-report-package",
            format!("the evidence package did not serialize: {error}"),
            "report the packaging defect; no partial package is retained",
        )
    })?;
    let checked = EvidencePackage::from_json(&package_json).map_err(|error| {
        report_error(
            "cli-solve-report-package-check",
            format!("the sealed package does not round-trip: {error}"),
            "report the packaging defect; the checker must accept the exact retained bytes",
        )
    })?;
    let check = check_package(&checked);
    if !check.passed() || check.merkle_root() != package_root {
        return Err(report_error(
            "cli-solve-report-package-check",
            "the solver-free checker refused the package this stage sealed",
            "report the packaging defect; a package the checker refuses is never retained",
        ));
    }

    // ---- Retain.
    let html_bytes = html.into_bytes();
    let json_bytes = json.into_bytes();
    let package_bytes = package_json.into_bytes();
    let total = html_bytes.len() + json_bytes.len() + package_bytes.len();
    work.charge(u64::try_from(total).unwrap_or(u64::MAX))
        .map_err(|error| invocation_work_refusal(Some(run), Some(SolveStage::Report), error))?;
    let html_hash = hash_bytes(&html_bytes);
    let json_hash = hash_bytes(&json_bytes);
    let package_hash = hash_bytes(&package_bytes);
    let receipt = format!(
        "{{\"schema\":{},\"run\":{},\"stage\":\"report\",\"project_hash\":{},\"report_html\":{},\"report_json\":{},\"report_content_hash\":{},\"package\":{},\"package_root\":{},\"checker\":{{\"passed\":true,\"protocol\":{}}},\"qoi_count\":1,\"verdict\":{},\"budget_terms_measured\":{measured_terms},\"budget_terms_total\":{},\"sources\":{{\"qoi_receipt\":{},\"conduction_receipt\":{},\"material_receipt\":{}}},\"authority\":\"projection-of-retained-receipts\",\"no_claim\":\"the report and package project retained stage receipts and add no physical, numerical, or validation authority; every claim keeps the colour its producer recorded (Estimated, unbounded dispersion)\"}}",
        json_string(REPORT_RECEIPT_SCHEMA),
        json_string(&run_hex),
        json_string(&project_hash.to_hex()),
        json_string(&html_hash.to_hex()),
        json_string(&json_hash.to_hex()),
        json_string(&content_hash.to_hex()),
        json_string(&package_hash.to_hex()),
        json_string(&package_root.to_hex()),
        fs_checker::CHECKER_PROTOCOL_VERSION,
        json_string(outcome),
        terms.len(),
        json_string(&qoi.completed.receipt.to_hex()),
        json_string(&conduction.completed.receipt.to_hex()),
        json_string(&materials.completed.receipt.to_hex()),
    );
    Ok(ReportStageProduct {
        receipt,
        artifacts: vec![
            RetainedSideArtifact {
                kind: REPORT_HTML_KIND,
                artifact: html_hash,
                bytes: html_bytes,
            },
            RetainedSideArtifact {
                kind: REPORT_JSON_KIND,
                artifact: json_hash,
                bytes: json_bytes,
            },
            RetainedSideArtifact {
                kind: EVIDENCE_PACKAGE_KIND,
                artifact: package_hash,
                bytes: package_bytes,
            },
        ],
    })
}

/// The exact retained bytes of a completed run, for the export verbs.
#[derive(Debug)]
pub(crate) struct CompletedRunExport {
    /// Run id, hex.
    pub(crate) run: String,
    /// Canonical project identity bound into the completed driver state.
    pub(crate) project_hash: String,
    /// `(stage, op id, receipt hash hex)` for every completed stage.
    pub(crate) stages: Vec<(&'static str, i64, String)>,
    /// The report stage receipt (JSON text).
    pub(crate) report_receipt: String,
    /// How the retained run was proven before export (`ResumeProof::as_str`).
    pub(crate) verification: &'static str,
    /// Retained HTML report bytes.
    pub(crate) report_html: Vec<u8>,
    /// Retained JSON twin bytes.
    pub(crate) report_json: Vec<u8>,
    /// Retained evidence package bytes (format-9 JSON).
    pub(crate) package_json: Vec<u8>,
}

/// Locate a completed run through the resume loader and return the exact
/// artifacts its report stage retained.
///
/// # Errors
/// [`SolveRefusal`] when the run id is malformed, unknown, incomplete, or when
/// the report stage's retained artifacts cannot be read back byte-for-byte.
pub(crate) fn load_completed_run(
    ledger: &Ledger,
    run_id_hex: &str,
) -> Result<CompletedRunExport, SolveRefusal> {
    let run = SolveRunId::parse_hex(run_id_hex).ok_or_else(|| {
        SolveRefusal::plain(
            "cli-solve-run-id",
            format!("`{run_id_hex}` is not a 64-hex run id"),
            "pass the run id printed by `frankensim solve`",
        )
    })?;
    let gate = CancelGate::new_clock_free();
    let invocation = InvocationWorkLedger::default();
    let work = EvidenceWork::new(&gate, None, &invocation);
    let VerifiedResume {
        state,
        last_expected_edges,
        ..
    } = {
        let t_state = std::time::Instant::now();
        let loaded = load_latest_state(ledger, run, work, ResumeProof::SealedEvidence)?;
        if std::env::var_os("FS_CLI_TRACE_EXPORT").is_some() {
            eprintln!(
                "TRACE export: load_latest_state {:.3}s",
                t_state.elapsed().as_secs_f64()
            );
        }
        loaded
    };
    if state.completed.len() < SolveStage::ALL.len() {
        let next =
            SolveStage::from_ordinal(u32::try_from(state.completed.len()).unwrap_or(u32::MAX))
                .map_or("unknown", SolveStage::name);
        return Err(SolveRefusal::plain(
            "cli-report-run-incomplete",
            format!(
                "run `{run_id_hex}` completed {} of {} stages; the report stage has not executed (next stage: {next})",
                state.completed.len(),
                SolveStage::ALL.len()
            ),
            format!(
                "resume with `frankensim solve --resume {run_id_hex} <ledger>` and rerun the export"
            ),
        ));
    }
    let report_stage = state
        .completed
        .get(SolveStage::Report.ordinal() as usize)
        .filter(|completed| completed.ordinal == SolveStage::Report.ordinal())
        .ok_or_else(|| {
            SolveRefusal::plain(
                "cli-report-run-incomplete",
                format!("run `{run_id_hex}` has no sealed report stage"),
                format!("resume with `frankensim solve --resume {run_id_hex} <ledger>`"),
            )
        })?;
    let ledger_refusal = |what: String| {
        SolveRefusal::plain(
            "cli-report-ledger",
            what,
            "repair or restore the ledger; exports read only retained artifacts",
        )
    };
    let receipt_bytes = ledger
        .get_artifact(&report_stage.receipt)
        .map_err(|error| ledger_refusal(format!("reading the report receipt failed: {error}")))?
        .ok_or_else(|| ledger_refusal("the report receipt is not retained".to_string()))?;
    let edges = ledger
        .op_artifact_edges_bounded(report_stage.op_id, EDGE_SCAN_CAP)
        .map_err(|error| ledger_refusal(format!("reading the report op edges failed: {error}")))?;
    if edges.truncated {
        return Err(ledger_refusal(format!(
            "the report op exceeds the {EDGE_SCAN_CAP}-edge verification cap"
        )));
    }
    require_authenticated_report_snapshot(
        run,
        report_stage.op_id,
        report_stage.receipt,
        &receipt_bytes,
        &edges.edges,
        &last_expected_edges,
        work,
        SolveStage::Report.ordinal() as usize,
    )?;
    let report_receipt = String::from_utf8(receipt_bytes)
        .map_err(|_| ledger_refusal("the report receipt is not UTF-8".to_string()))?;
    let mut report_html = None;
    let mut report_json = None;
    let mut package_json = None;
    for edge in &edges.edges {
        if edge.role != EdgeRole::Out {
            continue;
        }
        let Some(info) = ledger
            .artifact_info(&edge.artifact)
            .map_err(|error| ledger_refusal(format!("reading artifact info failed: {error}")))?
        else {
            continue;
        };
        let slot = match info.kind.as_str() {
            REPORT_HTML_KIND => &mut report_html,
            REPORT_JSON_KIND => &mut report_json,
            EVIDENCE_PACKAGE_KIND => &mut package_json,
            _ => continue,
        };
        let bytes = ledger
            .get_artifact(&edge.artifact)
            .map_err(|error| ledger_refusal(format!("reading a report artifact failed: {error}")))?
            .ok_or_else(|| ledger_refusal("a report artifact is not retained".to_string()))?;
        if hash_bytes(&bytes) != edge.artifact {
            return Err(ledger_refusal(format!(
                "retained artifact {} does not hash to its recorded identity",
                edge.artifact.to_hex()
            )));
        }
        if slot.replace(bytes).is_some() {
            return Err(ledger_refusal(format!(
                "the report op retained more than one `{}` artifact",
                info.kind
            )));
        }
    }
    let (Some(report_html), Some(report_json), Some(package_json)) =
        (report_html, report_json, package_json)
    else {
        return Err(ledger_refusal(
            "the report op did not retain all of html, json twin, and package".to_string(),
        ));
    };
    let stages = state
        .completed
        .iter()
        .map(|completed| {
            (
                SolveStage::from_ordinal(completed.ordinal).map_or("unknown", SolveStage::name),
                completed.op_id,
                completed.receipt.to_hex(),
            )
        })
        .collect();
    Ok(CompletedRunExport {
        run: run.to_hex(),
        project_hash: ContentHash(state.project).to_hex(),
        stages,
        report_receipt,
        report_html,
        report_json,
        package_json,
        verification: ResumeProof::SealedEvidence.as_str(),
    })
}

#[allow(clippy::too_many_arguments)]
fn require_authenticated_report_snapshot(
    run: SolveRunId,
    op: i64,
    expected_receipt: ContentHash,
    receipt_bytes: &[u8],
    edges: &[OpArtifactEdge],
    expected_edges: &[(EdgeRole, ContentHash)],
    work: EvidenceWork<'_>,
    candidate_index: usize,
) -> Result<(), SolveRefusal> {
    if hash_bytes(receipt_bytes) != expected_receipt {
        return Err(SolveRefusal::plain(
            "cli-report-ledger",
            "the reread report receipt does not hash to the identity re-attested by resume",
            "repair or restore the ledger; exports never follow a substituted report snapshot",
        ));
    }
    require_exact_edges(
        run,
        SolveStage::Report,
        op,
        edges,
        expected_edges,
        work,
        candidate_index,
    )
}

#[cfg(test)]
mod snapshot_tests {
    use super::*;

    #[test]
    fn reread_report_snapshot_refuses_receipt_and_edge_substitution() {
        let run = SolveRunId::parse_hex(&"11".repeat(32)).expect("run id");
        let receipt = hash_bytes(b"sealed report receipt");
        let output = hash_bytes(b"sealed report output");
        let expected_edges = [(EdgeRole::Out, receipt), (EdgeRole::Out, output)];
        let edges = [
            OpArtifactEdge {
                role: EdgeRole::Out,
                artifact: receipt,
            },
            OpArtifactEdge {
                role: EdgeRole::Out,
                artifact: output,
            },
        ];
        let gate = CancelGate::new_clock_free();
        let work = EvidenceWork::unmetered(&gate, None);

        require_authenticated_report_snapshot(
            run,
            7,
            receipt,
            b"sealed report receipt",
            &edges,
            &expected_edges,
            work,
            SolveStage::Report.ordinal() as usize,
        )
        .expect("the authenticated reread is accepted");

        let receipt_error = require_authenticated_report_snapshot(
            run,
            7,
            receipt,
            b"substituted report receipt",
            &edges,
            &expected_edges,
            work,
            SolveStage::Report.ordinal() as usize,
        )
        .expect_err("a substituted receipt is refused");
        assert_eq!(receipt_error.code, "cli-report-ledger");

        let substituted_edges = [
            edges[0],
            OpArtifactEdge {
                role: EdgeRole::Out,
                artifact: hash_bytes(b"substituted report output"),
            },
        ];
        require_authenticated_report_snapshot(
            run,
            7,
            receipt,
            b"sealed report receipt",
            &substituted_edges,
            &expected_edges,
            work,
            SolveStage::Report.ordinal() as usize,
        )
        .expect_err("a substituted edge set is refused");
    }
}
