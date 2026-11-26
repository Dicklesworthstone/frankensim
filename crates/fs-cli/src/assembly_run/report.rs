//! Probe evaluation, `EngineeringReport` synthesis, and format-9
//! `EvidencePackage` checking for the assembly pipeline (`.11.7` Stage 2).

use std::collections::BTreeMap;
use std::fmt::Write as _;

use fs_blake3::{ContentHash, hash_bytes};
use fs_checker::check as check_package;
use fs_evidence::Color;
use fs_package::{Claim, EvidencePackage, Provenance};
use fs_qty::Dims;
use fs_report::{
    EngineeringReport, LineageItem, NoClaimItem, QoiReportItem, ReportProvenance, StageReceiptItem,
};
use fs_scenario::assembly::draft::{ProbeDecl, ProbeKind, TimingEdge};
use fs_scenario::entity::EntityId;

use super::budget::{EvaluatedQoiBudget, QOI_BUDGET_SCHEMA};

use crate::import::json_string;
use crate::json_read::JsonValue;
use crate::solve::{SolveRefusal, SolveRunId};
use crate::stage_driver::{RetainedSideArtifact, staged_refusal};

pub(crate) const PROBE_SERIES_KIND: &str = "assembly-probe-series";
pub(crate) const ASSEMBLY_REPORT_HTML_KIND: &str = "assembly-report-html";
pub(crate) const ASSEMBLY_REPORT_JSON_KIND: &str = "assembly-report-json";
pub(crate) const ASSEMBLY_PACKAGE_KIND: &str = "assembly-evidence-package";
pub(crate) const ASSEMBLY_PROBE_SCHEMA: &str = "frankensim.cli.assembly-probe.v1";
pub(crate) const ASSEMBLY_REPORT_SCHEMA: &str = "frankensim.cli.assembly-report.v1";

const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);
const VELOCITY: Dims = Dims([1, 0, -1, 0, 0, 0]);
const TEMPERATURE: Dims = Dims([0, 0, 0, 1, 0, 0]);
const ANGULAR_VELOCITY: Dims = Dims([0, 0, -1, 0, 0, 0]);

#[derive(Debug, Clone)]
pub(crate) struct EntitySample {
    pub(crate) position_m: Option<f64>,
    pub(crate) velocity_m_s: Option<f64>,
    pub(crate) temperature_k: Option<f64>,
    /// A body moving in 3-D: velocity (m/s), then angular velocity (rad/s).
    pub(crate) motion: Option<[f64; 6]>,
}

#[derive(Debug, Clone)]
pub(crate) struct RunTrace {
    pub(crate) analysis: String,
    pub(crate) times_s: Vec<f64>,
    pub(crate) entities: BTreeMap<String, Vec<EntitySample>>,
    pub(crate) modal_modes: Vec<(f64, f64, f64)>,
    pub(crate) harmonic_points: Vec<(f64, f64, f64)>,
    pub(crate) stationary_norm: Option<f64>,
    pub(crate) stationary_state: Vec<(String, f64)>,
}

/// A trace sample's six-number `motion` list, when it has one.
fn motion_field(sample: &JsonValue) -> Option<[f64; 6]> {
    let list = sample.get("motion")?.as_array()?;
    let values: Option<Vec<f64>> = list.iter().map(JsonValue::as_f64).collect();
    values?.try_into().ok()
}

fn push_f64_list(out: &mut String, vals: &[f64]) {
    out.push('[');
    for (i, v) in vals.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        let _ = write!(out, "{v:.17e}");
    }
    out.push(']');
}

impl RunTrace {
    pub(crate) fn to_json(&self) -> String {
        let mut out = format!(
            "{{\"schema\":\"frankensim.cli.assembly-trace.v1\",\"analysis\":{},\"times_s\":",
            json_string(&self.analysis)
        );
        push_f64_list(&mut out, &self.times_s);
        out.push_str(",\"entities\":{");
        for (ei, (name, samples)) in self.entities.iter().enumerate() {
            if ei > 0 {
                out.push(',');
            }
            let _ = write!(out, "{}:[", json_string(name));
            for (si, s) in samples.iter().enumerate() {
                if si > 0 {
                    out.push(',');
                }
                let mut parts = Vec::new();
                if let Some(v) = s.position_m {
                    parts.push(format!("\"position_m\":{v:.17e}"));
                }
                if let Some(v) = s.velocity_m_s {
                    parts.push(format!("\"velocity_m_s\":{v:.17e}"));
                }
                if let Some(v) = s.temperature_k {
                    parts.push(format!("\"temperature_k\":{v:.17e}"));
                }
                if let Some(motion) = s.motion {
                    let mut list = String::new();
                    push_f64_list(&mut list, &motion);
                    parts.push(format!("\"motion\":{list}"));
                }
                let _ = write!(out, "{{{}}}", parts.join(","));
            }
            out.push(']');
        }
        out.push_str("},\"modal_modes\":[");
        for (i, (wn, hz, z)) in self.modal_modes.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            let _ = write!(out, "[{wn:.17e},{hz:.17e},{z:.17e}]");
        }
        out.push_str("],\"harmonic_points\":[");
        for (i, (hz, mag, ph)) in self.harmonic_points.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            let _ = write!(out, "[{hz:.17e},{mag:.17e},{ph:.17e}]");
        }
        out.push(']');
        if let Some(res) = self.stationary_norm {
            let _ = write!(out, ",\"stationary_residual\":{res:.17e}");
        }
        if !self.stationary_state.is_empty() {
            out.push_str(",\"stationary_state\":[");
            for (i, (label, val)) in self.stationary_state.iter().enumerate() {
                if i > 0 {
                    out.push(',');
                }
                let _ = write!(out, "[{},{val:.17e}]", json_string(label));
            }
            out.push(']');
        }
        out.push('}');
        out
    }

    pub(crate) fn from_json(text: &str) -> Result<Self, SolveRefusal> {
        let root = JsonValue::parse(text).map_err(|e| {
            staged_refusal(
                "cli-assembly-trace-parse",
                "probe",
                format!("trace JSON invalid: {e}"),
                "re-run assembly",
            )
        })?;
        let parse_triples = |key: &str| -> Vec<(f64, f64, f64)> {
            root.get(key)
                .and_then(JsonValue::as_array)
                .map(|a| {
                    a.iter()
                        .filter_map(|row| {
                            let r = row.as_array()?;
                            Some((
                                r.first()?.as_f64()?,
                                r.get(1)?.as_f64()?,
                                r.get(2)?.as_f64()?,
                            ))
                        })
                        .collect()
                })
                .unwrap_or_default()
        };
        let mut entities = BTreeMap::new();
        if let Some(JsonValue::Object(entries)) = root.get("entities") {
            for (k, v) in entries {
                if let Some(arr) = v.as_array() {
                    entities.insert(
                        k.clone(),
                        arr.iter()
                            .map(|o| EntitySample {
                                position_m: o.f64_field("position_m"),
                                velocity_m_s: o.f64_field("velocity_m_s"),
                                temperature_k: o.f64_field("temperature_k"),
                                motion: motion_field(o),
                            })
                            .collect(),
                    );
                }
            }
        }
        let stationary_state = root
            .get("stationary_state")
            .and_then(JsonValue::as_array)
            .map(|a| {
                a.iter()
                    .filter_map(|row| {
                        let r = row.as_array()?;
                        Some((r.first()?.as_str()?.to_string(), r.get(1)?.as_f64()?))
                    })
                    .collect()
            })
            .unwrap_or_default();
        Ok(Self {
            analysis: root
                .str_field("analysis")
                .unwrap_or("transient")
                .to_string(),
            times_s: root
                .get("times_s")
                .and_then(JsonValue::as_array)
                .map(|a| a.iter().filter_map(JsonValue::as_f64).collect())
                .unwrap_or_default(),
            entities,
            modal_modes: parse_triples("modal_modes"),
            harmonic_points: parse_triples("harmonic_points"),
            stationary_norm: root.f64_field("stationary_residual"),
            stationary_state,
        })
    }
}

#[derive(Debug, Clone)]
pub(crate) struct EvaluatedProbe {
    pub(crate) name: String,
    pub(crate) kind: String,
    pub(crate) unit: String,
    pub(crate) values: Vec<f64>,
    pub(crate) final_value: f64,
    pub(crate) peak_value: f64,
}

fn probe_refusal(probe: &str, why: impl Into<String>) -> SolveRefusal {
    staged_refusal(
        "cli-assembly-probe-unsupported",
        "probe",
        format!(
            "probe `{probe}` cannot be served by the lowered system: {}",
            why.into()
        ),
        "declare a lowering producing this observable or remove the probe; no synthetic values are fabricated",
    )
}

fn entity_series(
    trace: &RunTrace,
    names: &BTreeMap<EntityId, String>,
    entity: EntityId,
    probe: &str,
    get: impl Fn(&EntitySample) -> Option<f64>,
    label: &str,
) -> Result<Vec<f64>, SolveRefusal> {
    let ename = names
        .get(&entity)
        .ok_or_else(|| probe_refusal(probe, "unknown entity"))?;
    let samples = trace
        .entities
        .get(ename)
        .ok_or_else(|| probe_refusal(probe, format!("no trace for `{ename}`")))?;
    let vals: Vec<f64> = samples
        .iter()
        .map(|s| {
            get(s).ok_or_else(|| {
                let what = match (s.motion, label) {
                    (Some(_), "velocity_m_s") => {
                        "a 3-D velocity: probe `velocity-x`, `velocity-y` or `velocity-z`".into()
                    }
                    _ => format!("no lowered `{label}`"),
                };
                probe_refusal(probe, format!("`{ename}` has {what}"))
            })
        })
        .collect::<Result<_, _>>()?;
    if vals.is_empty() {
        return Err(probe_refusal(probe, format!("empty trace for `{ename}`")));
    }
    Ok(vals)
}

/// The index in [`EntitySample::motion`] and the unit of a body's component
/// `velocity-x|y|z` or `angular-velocity-x|y|z`, when `dims` are its own.
fn motion_component(name: &str, dims: Dims) -> Option<(usize, &'static str)> {
    let (quantity, axis) = name.rsplit_once('-')?;
    let k = ["x", "y", "z"].iter().position(|a| *a == axis)?;
    match quantity {
        "velocity" if dims == VELOCITY => Some((k, "m/s")),
        "angular-velocity" if dims == ANGULAR_VELOCITY => Some((k + 3, "rad/s")),
        _ => None,
    }
}

fn make_probe(name: impl Into<String>, kind: &str, unit: &str, values: Vec<f64>) -> EvaluatedProbe {
    let final_value = values.last().copied().unwrap_or(0.0);
    let peak = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    EvaluatedProbe {
        name: name.into(),
        kind: kind.to_string(),
        unit: unit.to_string(),
        values,
        final_value,
        peak_value: if peak.is_finite() { peak } else { 0.0 },
    }
}

#[allow(clippy::too_many_lines)]
pub(crate) fn evaluate_probes(
    probes: &[ProbeDecl],
    names: &BTreeMap<EntityId, String>,
    trace: &RunTrace,
) -> Result<(String, RetainedSideArtifact, Vec<EvaluatedProbe>), SolveRefusal> {
    let mut evaluated = Vec::new();
    if trace.analysis == "transient" {
        for p in probes {
            let (kind, unit, values) = match &p.kind {
                ProbeKind::Thermocouple {
                    location,
                    time_constant,
                    ..
                } => {
                    let raw = entity_series(
                        trace,
                        names,
                        *location,
                        &p.name,
                        |s| s.temperature_k,
                        "temperature_k",
                    )?;
                    let vals = if let Some(tau) = time_constant {
                        let mut out = vec![raw[0]];
                        for (i, &sample) in raw.iter().enumerate().skip(1) {
                            let dt = (trace.times_s[i] - trace.times_s[i - 1]).max(0.0);
                            let prev = *out.last().unwrap();
                            out.push(prev + (1.0 - (-dt / tau.value).exp()) * (sample - prev));
                        }
                        out
                    } else {
                        raw
                    };
                    ("thermocouple", "K", vals)
                }
                ProbeKind::Timing {
                    target,
                    threshold,
                    edge,
                } => {
                    let s = if threshold.dims == TEMPERATURE {
                        entity_series(
                            trace,
                            names,
                            *target,
                            &p.name,
                            |s| s.temperature_k,
                            "temperature_k",
                        )?
                    } else if threshold.dims == LENGTH {
                        entity_series(
                            trace,
                            names,
                            *target,
                            &p.name,
                            |s| s.position_m,
                            "position_m",
                        )?
                    } else if threshold.dims == VELOCITY {
                        entity_series(
                            trace,
                            names,
                            *target,
                            &p.name,
                            |s| s.velocity_m_s,
                            "velocity_m_s",
                        )?
                    } else {
                        return Err(probe_refusal(
                            &p.name,
                            "timing threshold dimensions match no lowered state",
                        ));
                    };
                    let mut crossings = Vec::new();
                    for i in 1..s.len() {
                        let (a, b) = (s[i - 1] - threshold.value, s[i] - threshold.value);
                        let hit = match edge {
                            TimingEdge::Rising => a < 0.0 && b >= 0.0,
                            TimingEdge::Falling => a > 0.0 && b <= 0.0,
                            TimingEdge::Both => (a < 0.0 && b >= 0.0) || (a > 0.0 && b <= 0.0),
                        };
                        if hit {
                            crossings.push(
                                trace.times_s[i - 1]
                                    + (-a / (b - a)) * (trace.times_s[i] - trace.times_s[i - 1]),
                            );
                        }
                    }
                    ("timing", "s", crossings)
                }
                ProbeKind::FieldQuantity {
                    target, quantity, ..
                } => match (quantity.name.as_str(), quantity.dims) {
                    ("position" | "displacement", d) if d == LENGTH => (
                        "field-quantity",
                        "m",
                        entity_series(
                            trace,
                            names,
                            *target,
                            &p.name,
                            |s| s.position_m,
                            "position_m",
                        )?,
                    ),
                    ("velocity", d) if d == VELOCITY => (
                        "field-quantity",
                        "m/s",
                        entity_series(
                            trace,
                            names,
                            *target,
                            &p.name,
                            |s| s.velocity_m_s,
                            "velocity_m_s",
                        )?,
                    ),
                    ("temperature", d) if d == TEMPERATURE => (
                        "field-quantity",
                        "K",
                        entity_series(
                            trace,
                            names,
                            *target,
                            &p.name,
                            |s| s.temperature_k,
                            "temperature_k",
                        )?,
                    ),
                    (name, dims) => {
                        let Some((k, unit)) = motion_component(name, dims) else {
                            return Err(probe_refusal(
                                &p.name,
                                format!("quantity `{}` is not lowered", quantity.name),
                            ));
                        };
                        let get = |s: &EntitySample| s.motion.map(|m| m[k]);
                        let series = entity_series(trace, names, *target, &p.name, get, name)?;
                        ("field-quantity", unit, series)
                    }
                },
                ProbeKind::Microphone { .. } => {
                    return Err(probe_refusal(&p.name, "no lowered acoustic field"));
                }
                ProbeKind::StrainGauge { .. } => {
                    return Err(probe_refusal(&p.name, "no lowered continuum strain field"));
                }
            };
            evaluated.push(make_probe(p.name.clone(), kind, unit, values));
        }
    } else if trace.analysis == "modal" {
        evaluated.push(make_probe(
            "modal.natural_rad_s",
            "modal-spectrum",
            "rad/s",
            trace.modal_modes.iter().map(|m| m.0).collect(),
        ));
        evaluated.push(make_probe(
            "modal.frequency_hz",
            "modal-spectrum",
            "Hz",
            trace.modal_modes.iter().map(|m| m.1).collect(),
        ));
        evaluated.push(make_probe(
            "modal.damping_ratio",
            "modal-spectrum",
            "1",
            trace.modal_modes.iter().map(|m| m.2).collect(),
        ));
    } else if trace.analysis == "harmonic" {
        evaluated.push(make_probe(
            "harmonic.transfer_magnitude",
            "harmonic-response",
            "SI",
            trace.harmonic_points.iter().map(|p| p.1).collect(),
        ));
    } else if let Some(res) = trace.stationary_norm {
        evaluated.push(make_probe(
            "stationary.residual_norm",
            "stationary-equilibrium",
            "SI",
            vec![res],
        ));
        for (label, val) in &trace.stationary_state {
            evaluated.push(make_probe(
                format!("stationary.state.{label}"),
                "stationary-equilibrium",
                "SI",
                vec![*val],
            ));
        }
    }

    let mut json = format!(
        "{{\"schema\":\"frankensim.cli.assembly-probe-series.v1\",\"analysis\":{},\"probes\":[",
        json_string(&trace.analysis)
    );
    for (i, ep) in evaluated.iter().enumerate() {
        if i > 0 {
            json.push(',');
        }
        let _ = write!(
            json,
            "{{\"name\":{},\"kind\":{},\"unit\":{},\"count\":{},\"final\":{:.17e},\"peak\":{:.17e},\"values\":",
            json_string(&ep.name),
            json_string(&ep.kind),
            json_string(&ep.unit),
            ep.values.len(),
            ep.final_value,
            ep.peak_value,
        );
        push_f64_list(&mut json, &ep.values);
        json.push('}');
    }
    json.push_str("]}");
    let bytes = json.into_bytes();
    let artifact = hash_bytes(&bytes);
    let receipt = format!(
        "{{\"schema\":{},\"analysis\":{},\"probe_count\":{},\"series_artifact\":{},\"authority\":\"Estimate\"}}",
        json_string(ASSEMBLY_PROBE_SCHEMA),
        json_string(&trace.analysis),
        evaluated.len(),
        json_string(&artifact.to_hex()),
    );
    Ok((
        receipt,
        RetainedSideArtifact {
            kind: PROBE_SERIES_KIND,
            artifact,
            bytes,
        },
        evaluated,
    ))
}

fn id_token(text: &str) -> String {
    let tok: String = text
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric()
                || matches!(c, '-' | '_' | '.' | '/' | ':' | '@' | '+' | '=')
            {
                c
            } else {
                '_'
            }
        })
        .collect();
    if tok.is_empty() {
        "undeclared".to_string()
    } else {
        tok
    }
}

pub(crate) struct AssemblyReportInputs<'a> {
    pub(crate) run: SolveRunId,
    pub(crate) assembly_name: &'a str,
    pub(crate) source_hash: ContentHash,
    pub(crate) constellation: &'a str,
    pub(crate) workspace: &'a str,
    pub(crate) seed: u64,
    pub(crate) wall_budget_s: f64,
    pub(crate) mem_budget_bytes: u64,
    pub(crate) analysis: &'a str,
    pub(crate) lowering_hash: ContentHash,
    pub(crate) plan_hash: &'a str,
    pub(crate) plan_artifact: ContentHash,
    pub(crate) ledger_artifact: Option<ContentHash>,
    pub(crate) energy_residual_j: Option<f64>,
    pub(crate) probe_series_artifact: ContentHash,
    pub(crate) probes: &'a [EvaluatedProbe],
    pub(crate) qoi_budget: Option<(ContentHash, &'a EvaluatedQoiBudget)>,
    pub(crate) stage_receipts: &'a [(&'static str, u32, ContentHash, String)],
}

#[allow(clippy::too_many_lines)]
pub(crate) fn build_assembly_report(
    inp: &AssemblyReportInputs<'_>,
) -> Result<(String, Vec<RetainedSideArtifact>), SolveRefusal> {
    let run_hex = inp.run.to_hex();
    let mut report = EngineeringReport::new(run_hex.clone(), inp.assembly_name.to_string())
        .with_provenance(ReportProvenance {
            code_version: format!(
                "fs-cli {}; workspace {}",
                env!("CARGO_PKG_VERSION"),
                inp.workspace
            ),
            constellation_lock: inp.constellation.to_string(),
            machine_fingerprint: "assembly-stage-driver (deterministic)".to_string(),
            units_system: "coherent SI (m, kg, s, K, A, mol)".to_string(),
            rng_seed: inp.seed,
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
            wall_budget_s: inp.wall_budget_s.max(0.0).floor() as u64,
            mem_budget_bytes: inp.mem_budget_bytes,
        });
    for ep in inp.probes {
        report = report.with_qoi(QoiReportItem {
            name: ep.name.clone(),
            description: format!(
                "{} ({}) from `{}` [{}]",
                ep.kind, ep.unit, inp.assembly_name, inp.analysis
            ),
            nominal_value: ep.final_value,
            unit: ep.unit.clone(),
            color: Color::Estimated {
                estimator: ASSEMBLY_PROBE_SCHEMA.to_string(),
                dispersion: f64::NAN,
            },
            discretization_error: f64::NAN,
            parameter_uncertainty: f64::NAN,
            surrogate_error: f64::NAN,
            total_uncertainty_budget: f64::NAN,
            source_root: inp.probe_series_artifact.to_hex(),
        });
    }
    for (name, ordinal, hash, summary) in inp.stage_receipts {
        report = report.with_stage(StageReceiptItem {
            stage: (*name).to_string(),
            ordinal: *ordinal,
            receipt_hash: hash.to_hex(),
            summary: vec![("detail".to_string(), summary.clone())],
        });
    }
    for (label, hash) in [
        ("canonical assembly source", inp.source_hash.to_hex()),
        ("lowering receipt", inp.lowering_hash.to_hex()),
        ("execution plan", inp.plan_artifact.to_hex()),
        ("probe series", inp.probe_series_artifact.to_hex()),
    ] {
        report = report.with_lineage(LineageItem {
            label: label.to_string(),
            hash,
        });
    }
    if let Some(lh) = inp.ledger_artifact {
        report = report.with_lineage(LineageItem {
            label: "window energy ledger".to_string(),
            hash: lh.to_hex(),
        });
    }
    if let Some((qh, qb)) = inp.qoi_budget {
        report = report.with_lineage(LineageItem {
            label: "QoI uncertainty budget".into(),
            hash: qh.to_hex(),
        });
        for r in &qb.requirements {
            report = report.with_requirement(r.clone());
        }
        for t in &qb.terms {
            report = report.with_budget_term(t.clone());
        }
    }
    let no_claim_stmt = if inp.qoi_budget.is_some() {
        "8-term QoI uncertainty budget evaluated; unmeasured terms remain NO-DATA and authority remains Estimate without external validation"
    } else {
        "projects retained stage receipts; Estimate-only authority, no 8-term QoI budget or experimental validation claimed"
    };
    report = report
        .with_no_claim(NoClaimItem {
            component: "assembly-report".to_string(),
            status: "Estimate".to_string(),
            statement: no_claim_stmt.to_string(),
        })
        .with_replay_command(format!("frankensim run <assembly.fsim> <ledger> --analysis {} # source {} run {run_hex}", inp.analysis, inp.source_hash.to_hex()));

    let (html_bytes, json_bytes, content_hash) = (
        report.render_html().into_bytes(),
        report.render_json().into_bytes(),
        report.content_hash(),
    );
    let mut package = EvidencePackage::new(Provenance::new(
        format!(
            "fs-cli@{}+workspace={}",
            env!("CARGO_PKG_VERSION"),
            id_token(inp.workspace)
        ),
        id_token(inp.constellation),
    ))
    .with_claim(Claim::estimated(
        format!("assembly.lowering.{}", id_token(inp.assembly_name)),
        format!(
            "lowering receipt {} plan hash {}",
            inp.lowering_hash.to_hex(),
            inp.plan_hash
        ),
        fs_system::LOWERING_RECEIPT_SCHEMA,
        f64::INFINITY,
    ));
    if let Some(defect) = inp.energy_residual_j {
        package = package.with_claim(Claim::estimated(
            format!("assembly.energy_ledger.{}", id_token(inp.assembly_name)),
            format!("window energy residual {defect:.3e} J"),
            "frankensim.window_ledger.v1",
            f64::INFINITY,
        ));
    }
    for ep in inp.probes {
        package = package.with_claim(Claim::estimated(
            format!("probe.{}", id_token(&ep.name)),
            format!(
                "{} ({}) final={:.12e} peak={:.12e}",
                ep.name, ep.unit, ep.final_value, ep.peak_value
            ),
            ASSEMBLY_PROBE_SCHEMA,
            f64::INFINITY,
        ));
    }
    if let Some((_, qb)) = inp.qoi_budget {
        for r in &qb.requirements {
            let desc = format!(
                "{} {} (margin={:.12e}, limit={:.12e} {})",
                r.qoi, r.outcome, r.nominal_margin, r.effective_limit, r.unit
            );
            package = package.with_claim(Claim::estimated(
                format!("requirement.{}", id_token(&r.id)),
                desc,
                QOI_BUDGET_SCHEMA,
                f64::INFINITY,
            ));
        }
    }
    let package_root = package.try_merkle_root().map_err(|e| {
        staged_refusal(
            "cli-assembly-report-package",
            "report",
            format!("package seal failed: {e}"),
            "report packaging defect",
        )
    })?;
    let package_json = package.to_json().map_err(|e| {
        staged_refusal(
            "cli-assembly-report-package",
            "report",
            format!("package serialize failed: {e}"),
            "report packaging defect",
        )
    })?;
    let parsed_pkg = EvidencePackage::from_json(&package_json).map_err(|e| {
        staged_refusal(
            "cli-assembly-report-package-check",
            "report",
            format!("package round-trip failed: {e}"),
            "report packaging defect",
        )
    })?;
    let check = check_package(&parsed_pkg);
    if !check.passed() || check.merkle_root() != package_root {
        return Err(staged_refusal(
            "cli-assembly-report-package-check",
            "report",
            "solver-free checker refused the assembly package",
            "a refused package is never retained",
        ));
    }

    let pkg_bytes = package_json.into_bytes();
    let (html_hash, json_hash, pkg_hash) = (
        hash_bytes(&html_bytes),
        hash_bytes(&json_bytes),
        hash_bytes(&pkg_bytes),
    );
    let receipt = format!(
        "{{\"schema\":{},\"run\":{},\"stage\":\"report\",\"assembly\":{},\"analysis\":{},\"source_hash\":{},\"lowering_receipt\":{},\"plan_hash\":{},\"report_html\":{},\"report_json\":{},\"report_content_hash\":{},\"package\":{},\"package_root\":{},\"checker\":{{\"passed\":true}},\"probe_count\":{},\"authority\":\"Estimate\"}}",
        json_string(ASSEMBLY_REPORT_SCHEMA),
        json_string(&run_hex),
        json_string(inp.assembly_name),
        json_string(inp.analysis),
        json_string(&inp.source_hash.to_hex()),
        json_string(&inp.lowering_hash.to_hex()),
        json_string(inp.plan_hash),
        json_string(&html_hash.to_hex()),
        json_string(&json_hash.to_hex()),
        json_string(&content_hash.to_hex()),
        json_string(&pkg_hash.to_hex()),
        json_string(&package_root.to_hex()),
        inp.probes.len(),
    );
    Ok((
        receipt,
        vec![
            RetainedSideArtifact {
                kind: ASSEMBLY_REPORT_HTML_KIND,
                artifact: html_hash,
                bytes: html_bytes,
            },
            RetainedSideArtifact {
                kind: ASSEMBLY_REPORT_JSON_KIND,
                artifact: json_hash,
                bytes: json_bytes,
            },
            RetainedSideArtifact {
                kind: ASSEMBLY_PACKAGE_KIND,
                artifact: pkg_hash,
                bytes: pkg_bytes,
            },
        ],
    ))
}

#[cfg(test)]
mod tests {
    use fs_qty::QtyAny;
    use fs_scenario::assembly::ParameterKey;
    use fs_scenario::{EntityCatalog, EntityDeclaration, FrameId};

    use super::*;

    #[test]
    fn a_body_in_3d_serves_its_motion_components_and_refuses_plain_velocity() {
        let mut catalog = EntityCatalog::new();
        let root = catalog.declare(EntityDeclaration::assembly("a")).unwrap();
        let part = EntityDeclaration::part(root, "disc");
        let disc = catalog.declare(part).unwrap();
        let names = BTreeMap::from([(disc, "disc".to_owned())]);
        let sample = EntitySample {
            position_m: None,
            velocity_m_s: None,
            temperature_k: None,
            motion: Some([1.0, 2.0, 3.0, 4.0, 5.0, 6.0]),
        };
        let trace = RunTrace {
            analysis: "transient".to_owned(),
            times_s: vec![0.0, 0.1],
            entities: BTreeMap::from([("disc".to_owned(), vec![sample; 2])]),
            modal_modes: Vec::new(),
            harmonic_points: Vec::new(),
            stationary_norm: None,
            stationary_state: Vec::new(),
        };
        // The motion survives the written trace.
        let trace = RunTrace::from_json(&trace.to_json()).unwrap();
        let probe = |quantity: &str, dims| {
            let kind = ProbeKind::FieldQuantity {
                target: disc,
                position: [QtyAny::new(0.0, LENGTH); 3],
                quantity: ParameterKey::new(quantity, dims),
            };
            let probe = ProbeDecl {
                name: quantity.to_owned(),
                frame: FrameId(0),
                sample_rate: QtyAny::new(10.0, Dims([0, 0, -1, 0, 0, 0])),
                kind,
            };
            evaluate_probes(&[probe], &names, &trace).map(|(_, _, probes)| probes[0].clone())
        };
        let spin = probe("angular-velocity-z", ANGULAR_VELOCITY).unwrap();
        assert_eq!((spin.values, spin.unit.as_str()), (vec![6.0; 2], "rad/s"));
        assert_eq!(probe("velocity-y", VELOCITY).unwrap().values, [2.0; 2]);
        let plain = probe("velocity", VELOCITY).unwrap_err().what;
        let want = "`disc` has a 3-D velocity: probe `velocity-x`, `velocity-y` or `velocity-z`";
        assert!(plain.ends_with(want), "{plain}");
        let wrong = probe("angular-velocity-z", VELOCITY).unwrap_err().what;
        let want = "quantity `angular-velocity-z` is not lowered";
        assert!(wrong.ends_with(want), "{wrong}");
    }
}
