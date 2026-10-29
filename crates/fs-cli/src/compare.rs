//! `compare`: project retained receipts without replaying physics
//! (bead frankensim-rc-root-q61wp.47).
//!
//! Runs need at least one common QoI name, not identical output/requirement
//! sets. Rows are joined by their declared keys; duplicate keys are malformed
//! receipts, never an invitation to select an arbitrary row. Missing values and
//! unavailable arithmetic are JSON null, not zero. A converted delta names its
//! coherent SI unit; the original nominal values and units remain untouched.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;
use std::path::Path;

use fs_blake3::ContentHash;
use fs_ledger::Ledger;

use crate::json_read::JsonValue;
use crate::report::{LoadedExport, load_export_from, open_export_ledger, refuse};
use crate::{CommandOutput, OutputMode, RESULT_SCHEMA, escape_text, exit, push_json_string};

const COMMAND: &str = "compare";
const NO_CLAIM: &str = "the comparison projects retained stage receipts of two completed runs and adds no physical, numerical, or validation authority; a changed value is a change between two Estimated candidates, not a verified effect";
const BINDING_KEYS: [&str; 3] = ["run", "project_hash", "import_op"];

struct QoiValue {
    unit: String,
    nominal: String,
    value: f64,
    color: String,
    identity: String,
}

struct QoiDiff {
    name: String,
    left: Option<QoiValue>,
    right: Option<QoiValue>,
    delta: Option<f64>,
    rel_delta: Option<f64>,
    delta_unit: Option<String>,
    reason: Option<&'static str>,
}

impl QoiDiff {
    fn classification(&self) -> &'static str {
        match (&self.left, &self.right) {
            (None, _) => "added",
            (_, None) => "removed",
            (Some(_), Some(_)) if self.reason.is_some() => "incomparable",
            (Some(left), Some(right)) => {
                if left.value == right.value && left.unit == right.unit && left.color == right.color {
                    "same"
                } else {
                    "changed"
                }
            }
        }
    }
    fn changed(&self) -> bool { self.classification() != "same" }
    fn identity_same(&self) -> bool {
        matches!((&self.left, &self.right), (Some(l), Some(r)) if l.identity == r.identity)
    }
}

struct RequirementValue {
    outcome: String,
    nominal_margin: String,
    effective_limit: String,
}

struct RequirementDiff {
    id: String,
    left: Option<RequirementValue>,
    right: Option<RequirementValue>,
}

impl RequirementDiff {
    fn verdict_changed(&self) -> bool {
        matches!((&self.left, &self.right), (Some(l), Some(r)) if l.outcome != r.outcome)
    }
    fn classification(&self) -> &'static str {
        match (&self.left, &self.right) {
            (None, _) => "added",
            (_, None) => "removed",
            (Some(l), Some(r)) if l.outcome == r.outcome
                && l.nominal_margin == r.nominal_margin && l.effective_limit == r.effective_limit => "same",
            _ => "changed",
        }
    }
    fn changed(&self) -> bool { self.classification() != "same" }
}

#[derive(PartialEq)]
struct TermValue {
    state: String,
    half_width: Option<String>,
}

struct TermDiff {
    kind: String,
    left: Option<TermValue>,
    right: Option<TermValue>,
}

impl TermDiff {
    fn changed(&self) -> bool { self.left != self.right }
    fn classification(&self) -> &'static str {
        match (&self.left, &self.right) {
            (None, _) => "added",
            (_, None) => "removed",
            _ if self.changed() => "changed",
            _ => "same",
        }
    }
}

struct PackDiff {
    kind: String,
    card_left: Option<String>,
    card_right: Option<String>,
    identity_left: Option<String>,
    identity_right: Option<String>,
}

impl PackDiff {
    fn changed(&self) -> bool {
        self.card_left != self.card_right || self.identity_left != self.identity_right
    }
}

struct StageDiff {
    stage: &'static str,
    left: String,
    right: String,
    differing_keys: Vec<String>,
}

impl StageDiff {
    fn status(&self) -> &'static str {
        if self.differing_keys.is_empty() {
            "unchanged (same receipt)"
        } else if self.differing_keys.iter().all(|key| BINDING_KEYS.contains(&key.as_str())) {
            "unchanged (same inputs; differs only by binding keys)"
        } else {
            "changed"
        }
    }
    fn changed(&self) -> bool { self.status() == "changed" }
}

struct Comparison {
    left_run: String,
    right_run: String,
    project_hash_left: String,
    project_hash_right: String,
    verification: &'static str,
    qoi: Vec<QoiDiff>,
    requirements: Vec<RequirementDiff>,
    terms: Vec<TermDiff>,
    pack_set_root_left: String,
    pack_set_root_right: String,
    packs: Vec<PackDiff>,
    stages: Vec<StageDiff>,
}

impl Comparison {
    fn same_project(&self) -> bool { self.project_hash_left == self.project_hash_right }
    fn changed(&self) -> bool {
        !self.same_project() || self.qoi.iter().any(QoiDiff::changed)
            || self.requirements.iter().any(RequirementDiff::changed)
            || self.terms.iter().any(TermDiff::changed)
            || self.packs.iter().any(PackDiff::changed)
            || self.pack_set_root_left != self.pack_set_root_right
            || self.stages.iter().any(StageDiff::changed)
    }
    fn summary(&self) -> String {
        if !self.changed() {
            return "identical runs: no differences in any retained receipt".to_string();
        }
        let mut parts = vec![format!("{} of {} stages changed",
            self.stages.iter().filter(|stage| stage.changed()).count(), self.stages.len())];
        if !self.same_project() {
            parts.push(format!("project hash {} -> {}", self.project_hash_left, self.project_hash_right));
        }
        for qoi in self.qoi.iter().filter(|qoi| qoi.changed()) {
            parts.push(format!("{} {} -> {} {} ({})", qoi.name,
                qoi.left.as_ref().map_or("absent", |v| v.nominal.as_str()),
                qoi.right.as_ref().map_or("absent", |v| v.nominal.as_str()),
                qoi.right.as_ref().map_or("", |v| v.unit.as_str()), qoi.classification()));
        }
        for requirement in self.requirements.iter().filter(|r| r.changed()) {
            match (&requirement.left, &requirement.right) {
                (Some(l), Some(r)) if !requirement.verdict_changed() => parts.push(format!(
                    "verdict {} stays {} (margin {} -> {})", requirement.id, r.outcome, l.nominal_margin, r.nominal_margin)),
                (Some(l), Some(r)) => parts.push(format!("verdict {} {} -> {}", requirement.id, l.outcome, r.outcome)),
                _ => parts.push(format!("requirement {} {}", requirement.id, requirement.classification())),
            }
        }
        for term in self.terms.iter().filter(|term| term.changed()) {
            parts.push(format!("budget term {} {} -> {}", term.kind,
                term.left.as_ref().map_or("absent", |v| v.state.as_str()),
                term.right.as_ref().map_or("absent", |v| v.state.as_str())));
        }
        for pack in self.packs.iter().filter(|pack| pack.changed()) {
            parts.push(format!("{} card {} -> {}", pack.kind,
                pack.card_left.as_deref().unwrap_or("absent"), pack.card_right.as_deref().unwrap_or("absent")));
        }
        parts.join("; ")
    }
}

fn shape_refusal(mode: OutputMode, subject: &str, why: impl Into<String>) -> CommandOutput {
    refuse(mode, COMMAND, exit::REFUSED, "cli-compare-receipt-shape", subject, why,
        "regenerate the run; compare reads retained receipts and never repairs one")
}

fn stage_receipt(ledger: &Ledger, loaded: &LoadedExport, stage: &str) -> Result<(String, JsonValue), String> {
    let mut hashes = loaded.export.stages.iter().filter(|(name, _, _)| *name == stage)
        .map(|(_, _, hash)| hash.as_str());
    let hash = hashes.next().ok_or_else(|| format!("run {} has no `{stage}` stage receipt", loaded.export.run))?;
    if hashes.next().is_some() {
        return Err(format!("run {} has multiple `{stage}` stage receipts", loaded.export.run));
    }
    let content = ContentHash::from_hex(hash).ok_or_else(|| format!("`{stage}` receipt hash `{hash}` is not a content hash"))?;
    let bytes = ledger.get_artifact(&content).map_err(|error| format!("reading the `{stage}` receipt failed: {error}"))?
        .ok_or_else(|| format!("the `{stage}` receipt {hash} is not retained in this ledger"))?;
    let text = String::from_utf8(bytes).map_err(|_| format!("the `{stage}` receipt bytes are not UTF-8"))?;
    let value = JsonValue::parse(&text).map_err(|error| format!("the `{stage}` receipt does not parse: {error}"))?;
    Ok((hash.to_string(), value))
}

fn required_str(value: &JsonValue, key: &str, context: &str) -> Result<String, String> {
    value.str_field(key).map(str::to_string).ok_or_else(|| format!("{context} is missing string field `{key}`"))
}

fn required_number(value: &JsonValue, key: &str, context: &str) -> Result<String, String> {
    let raw = value.get(key).and_then(JsonValue::number_raw)
        .ok_or_else(|| format!("{context} is missing numeric field `{key}`"))?;
    match raw.parse::<f64>() {
        Ok(number) if number.is_finite() => Ok(raw.to_string()),
        _ => Err(format!("{context} numeric field `{key}` is not a finite f64")),
    }
}

fn required_array<'a>(value: &'a JsonValue, key: &str, context: &str) -> Result<&'a [JsonValue], String> {
    value.get(key).and_then(JsonValue::as_array).ok_or_else(|| format!("{context} is missing array field `{key}`"))
}

/// Index each side independently so a duplicate cannot disappear in a union.
fn indexed<'a>(rows: &'a [JsonValue], key: &str, context: &str) -> Result<BTreeMap<String, &'a JsonValue>, String> {
    let mut index = BTreeMap::new();
    for row in rows {
        let name = required_str(row, key, context)?;
        if index.insert(name.clone(), row).is_some() {
            return Err(format!("{context} contains duplicate `{key}` `{name}`"));
        }
    }
    Ok(index)
}

fn union_keys<L, R>(left: &BTreeMap<String, L>, right: &BTreeMap<String, R>) -> BTreeSet<String> {
    left.keys().chain(right.keys()).cloned().collect()
}

fn qoi_value(row: &JsonValue) -> Result<QoiValue, String> {
    let nominal = required_number(row, "value", "QoI row")?;
    Ok(QoiValue {
        unit: required_str(row, "unit", "QoI row")?,
        value: nominal.parse::<f64>().map_err(|error| error.to_string())?,
        nominal,
        color: required_str(row, "color", "QoI row")?,
        identity: required_str(row, "identity", "QoI row")?,
    })
}

/// Differing unit spellings are converted by the canonical quantity parser,
/// not a second unit table. The common QoI name is the semantic join key;
/// dimensional equality alone never matches different named observables.
fn set_delta(row: &mut QoiDiff) {
    let (Some(left), Some(right)) = (&row.left, &row.right) else { return; };
    let affine = left.unit.trim() == "degC" || right.unit.trim() == "degC";
    let (l, r, unit) = if left.unit == right.unit {
        (left.value, right.value, if affine { "K".to_string() } else { left.unit.clone() })
    } else {
        let l = fs_qty::parse::parse_qty(&format!("{} {}", left.nominal, left.unit));
        let r = fs_qty::parse::parse_qty(&format!("{} {}", right.nominal, right.unit));
        match (l, r) {
            (Ok(l), Ok(r)) if l.dims == r.dims => (l.value, r.value, l.dims.unit_string()),
            (Ok(_), Ok(_)) => { row.reason = Some("incompatible-dimensions"); return; }
            _ => { row.reason = Some("unit-conversion-unavailable"); return; }
        }
    };
    row.delta_unit = Some(unit);
    let delta = r - l;
    if !delta.is_finite() {
        row.reason = Some("delta-overflow");
        return;
    }
    row.delta = Some(delta);
    // Zero baselines and affine origins do not define a relative change.
    if l != 0.0 && !affine {
        let relative = delta / l.abs();
        if relative.is_finite() { row.rel_delta = Some(relative); }
    }
}

fn diff_qoi(left: &JsonValue, right: &JsonValue) -> Result<Vec<QoiDiff>, String> {
    let left = indexed(required_array(left, "qoi", "left QoI receipt")?, "name", "left QoI row")?;
    let right = indexed(required_array(right, "qoi", "right QoI receipt")?, "name", "right QoI row")?;
    union_keys(&left, &right).into_iter().map(|name| {
        let mut row = QoiDiff {
            left: left.get(&name).map(|v| qoi_value(v)).transpose()?,
            right: right.get(&name).map(|v| qoi_value(v)).transpose()?,
            name, delta: None, rel_delta: None, delta_unit: None, reason: None,
        };
        set_delta(&mut row);
        Ok(row)
    }).collect()
}

fn requirement_value(row: &JsonValue) -> Result<RequirementValue, String> {
    Ok(RequirementValue {
        outcome: required_str(row, "outcome", "requirement row")?,
        nominal_margin: required_number(row, "nominal_margin_kelvin", "requirement row")?,
        effective_limit: required_number(row, "effective_limit_kelvin", "requirement row")?,
    })
}

fn diff_requirements(left: &JsonValue, right: &JsonValue) -> Result<Vec<RequirementDiff>, String> {
    let left = indexed(required_array(left, "requirements", "left QoI receipt")?, "id", "left requirement row")?;
    let right = indexed(required_array(right, "requirements", "right QoI receipt")?, "id", "right requirement row")?;
    union_keys(&left, &right).into_iter().map(|id| Ok(RequirementDiff {
        left: left.get(&id).map(|v| requirement_value(v)).transpose()?,
        right: right.get(&id).map(|v| requirement_value(v)).transpose()?, id,
    })).collect()
}

fn budget_terms<'a>(receipt: &'a JsonValue, side: &str) -> Result<&'a [JsonValue], String> {
    let budgets = required_array(receipt, "budget", &format!("{side} QoI receipt"))?;
    let budget = budgets.first().ok_or_else(|| format!("{side} QoI receipt carries no budget"))?;
    required_array(budget, "terms", &format!("{side} budget"))
}

fn term_value(row: &JsonValue) -> Result<TermValue, String> {
    let half_width = match row.get("upper_kelvin").and_then(JsonValue::number_raw) {
        Some(_) => Some(required_number(row, "upper_kelvin", "budget term")?),
        None => None,
    };
    Ok(TermValue { state: required_str(row, "state", "budget term")?, half_width })
}

fn diff_terms(left: &JsonValue, right: &JsonValue) -> Result<Vec<TermDiff>, String> {
    let left = indexed(budget_terms(left, "left")?, "kind", "left budget term")?;
    let right = indexed(budget_terms(right, "right")?, "kind", "right budget term")?;
    union_keys(&left, &right).into_iter().map(|kind| Ok(TermDiff {
        left: left.get(&kind).map(|v| term_value(v)).transpose()?,
        right: right.get(&kind).map(|v| term_value(v)).transpose()?, kind,
    })).collect()
}

fn diff_packs(left: &JsonValue, right: &JsonValue) -> Result<Vec<PackDiff>, String> {
    let packs_left = required_array(left, "packs", "left material receipt")?;
    let packs_right = required_array(right, "packs", "right material receipt")?;
    let mut kinds = Vec::new();
    for pack in packs_left.iter().chain(packs_right) {
        let kind = required_str(pack, "kind", "card pack row")?;
        if !kinds.contains(&kind) { kinds.push(kind); }
    }
    let find = |packs: &[JsonValue], kind: &str| -> Result<(Option<String>, Option<String>), String> {
        let Some(pack) = packs.iter().find(|pack| pack.str_field("kind") == Some(kind)) else {
            return Ok((None, None));
        };
        Ok((Some(required_str(pack, "card", "card pack row")?), Some(required_str(pack, "identity", "card pack row")?)))
    };
    let mut rows = Vec::with_capacity(kinds.len());
    for kind in kinds {
        let (card_left, identity_left) = find(packs_left, &kind)?;
        let (card_right, identity_right) = find(packs_right, &kind)?;
        rows.push(PackDiff { kind, card_left, card_right, identity_left, identity_right });
    }
    Ok(rows)
}

fn differing_keys(left: &JsonValue, right: &JsonValue) -> Vec<String> {
    if left == right { return Vec::new(); }
    let (Some(left), Some(right)) = (left.as_object(), right.as_object()) else {
        return vec!["<receipt>".to_string()];
    };
    let mut keys = Vec::new();
    for (key, value) in left {
        match right.iter().find(|(other, _)| other == key) {
            Some((_, other)) if other == value => {}
            _ => keys.push(key.clone()),
        }
    }
    for (key, _) in right {
        if !left.iter().any(|(other, _)| other == key) { keys.push(key.clone()); }
    }
    keys
}

enum BuildError { Incomparable(String), Shape(String) }
impl From<String> for BuildError {
    fn from(why: String) -> Self { Self::Shape(why) }
}

fn comparability(left: &JsonValue, right: &JsonValue) -> Result<(), BuildError> {
    let left = indexed(required_array(left, "qoi", "left QoI receipt")?, "name", "left QoI row")?;
    let right = indexed(required_array(right, "qoi", "right QoI receipt")?, "name", "right QoI row")?;
    if !left.keys().any(|name| right.contains_key(name)) {
        return Err(BuildError::Incomparable(format!("the runs measure different QoIs with no common name ({} vs {})",
            left.keys().cloned().collect::<Vec<_>>().join(","), right.keys().cloned().collect::<Vec<_>>().join(","))));
    }
    Ok(())
}

fn build(ledger: &Ledger, left: &LoadedExport, right: &LoadedExport) -> Result<Comparison, BuildError> {
    let (_, qoi_left) = stage_receipt(ledger, left, "qoi")?;
    let (_, qoi_right) = stage_receipt(ledger, right, "qoi")?;
    comparability(&qoi_left, &qoi_right)?;
    let (_, material_left) = stage_receipt(ledger, left, "material-resolve")?;
    let (_, material_right) = stage_receipt(ledger, right, "material-resolve")?;
    let mut stages = Vec::with_capacity(left.export.stages.len());
    for (stage, _, hash_left) in &left.export.stages {
        let (_, receipt_left) = stage_receipt(ledger, left, stage)?;
        let (hash_right, receipt_right) = stage_receipt(ledger, right, stage)?;
        stages.push(StageDiff { stage, left: hash_left.clone(), right: hash_right,
            differing_keys: differing_keys(&receipt_left, &receipt_right) });
    }
    if right.export.stages.len() != left.export.stages.len() {
        return Err(BuildError::Shape(format!("the runs completed different stage counts ({} vs {})",
            left.export.stages.len(), right.export.stages.len())));
    }
    Ok(Comparison {
        left_run: left.export.run.clone(), right_run: right.export.run.clone(),
        project_hash_left: left.export.project_hash.clone(), project_hash_right: right.export.project_hash.clone(),
        verification: left.export.verification, qoi: diff_qoi(&qoi_left, &qoi_right)?,
        requirements: diff_requirements(&qoi_left, &qoi_right)?, terms: diff_terms(&qoi_left, &qoi_right)?,
        pack_set_root_left: required_str(&material_left, "pack_set_root", "left material receipt")?,
        pack_set_root_right: required_str(&material_right, "pack_set_root", "right material receipt")?,
        packs: diff_packs(&material_left, &material_right)?, stages,
    })
}

fn push_opt(out: &mut String, value: Option<&str>) {
    match value { Some(v) => push_json_string(out, v), None => out.push_str("null") }
}
fn number(value: Option<f64>) -> String {
    value.map_or_else(|| "null".to_string(), |value| value.to_string())
}

fn render_json(comparison: &Comparison, subject: &str) -> String {
    let mut out = String::from("{\"schema\":");
    push_json_string(&mut out, RESULT_SCHEMA);
    out.push_str(",\"command\":\"compare\",\"status\":\"ok\",\"subject\":");
    push_json_string(&mut out, subject);
    for (key, value) in [("left_run", &comparison.left_run), ("right_run", &comparison.right_run),
        ("project_hash_left", &comparison.project_hash_left), ("project_hash_right", &comparison.project_hash_right)] {
        let _ = write!(out, ",\"{key}\":"); push_json_string(&mut out, value);
    }
    let _ = write!(out, ",\"same_project\":{},\"changed\":{},\"summary\":", comparison.same_project(), comparison.changed());
    push_json_string(&mut out, &comparison.summary());
    let _ = write!(out, ",\"qoi_count\":{},\"qoi_diffs\":[", comparison.qoi.len());
    for (index, row) in comparison.qoi.iter().enumerate() {
        if index > 0 { out.push(','); }
        out.push_str("{\"name\":"); push_json_string(&mut out, &row.name);
        let l = row.left.as_ref(); let r = row.right.as_ref();
        out.push_str(",\"unit_left\":"); push_opt(&mut out, l.map(|v| v.unit.as_str()));
        out.push_str(",\"unit_right\":"); push_opt(&mut out, r.map(|v| v.unit.as_str()));
        let _ = write!(out, ",\"nominal_left\":{},\"nominal_right\":{},\"delta\":{},\"rel_delta\":{},\"color_left\":",
            l.map_or("null", |v| v.nominal.as_str()), r.map_or("null", |v| v.nominal.as_str()), number(row.delta), number(row.rel_delta));
        push_opt(&mut out, l.map(|v| v.color.as_str()));
        out.push_str(",\"color_right\":"); push_opt(&mut out, r.map(|v| v.color.as_str()));
        out.push_str(",\"color_evolution\":");
        let evolution = match (l, r) {
            (Some(l), Some(r)) if l.color == r.color => "same".to_string(),
            (Some(l), Some(r)) => format!("{}->{}", l.color, r.color),
            _ => row.classification().to_string(),
        };
        push_json_string(&mut out, &evolution);
        let _ = write!(out, ",\"identity_same\":{},\"classification\":\"{}\",\"delta_unit\":", row.identity_same(), row.classification());
        push_opt(&mut out, row.delta_unit.as_deref());
        out.push_str(",\"comparison_reason\":"); push_opt(&mut out, row.reason);
        out.push('}');
    }
    out.push_str("],\"requirements\":[");
    for (index, row) in comparison.requirements.iter().enumerate() {
        if index > 0 { out.push(','); }
        out.push_str("{\"id\":"); push_json_string(&mut out, &row.id);
        let l = row.left.as_ref(); let r = row.right.as_ref();
        out.push_str(",\"outcome_left\":"); push_opt(&mut out, l.map(|v| v.outcome.as_str()));
        out.push_str(",\"outcome_right\":"); push_opt(&mut out, r.map(|v| v.outcome.as_str()));
        let _ = write!(out, ",\"verdict_changed\":{},\"nominal_margin_left\":{},\"nominal_margin_right\":{},\"effective_limit_left\":{},\"effective_limit_right\":{},\"classification\":\"{}\",\"changed\":{}}}",
            row.verdict_changed(), l.map_or("null", |v| v.nominal_margin.as_str()), r.map_or("null", |v| v.nominal_margin.as_str()),
            l.map_or("null", |v| v.effective_limit.as_str()), r.map_or("null", |v| v.effective_limit.as_str()), row.classification(), row.changed());
    }
    out.push_str("],\"budget_terms\":[");
    for (index, row) in comparison.terms.iter().enumerate() {
        if index > 0 { out.push(','); }
        out.push_str("{\"kind\":"); push_json_string(&mut out, &row.kind);
        out.push_str(",\"state_left\":"); push_opt(&mut out, row.left.as_ref().map(|v| v.state.as_str()));
        out.push_str(",\"state_right\":"); push_opt(&mut out, row.right.as_ref().map(|v| v.state.as_str()));
        let _ = write!(out, ",\"half_width_left_kelvin\":{},\"half_width_right_kelvin\":{},\"changed\":{},\"classification\":\"{}\"}}",
            row.left.as_ref().and_then(|v| v.half_width.as_deref()).unwrap_or("null"),
            row.right.as_ref().and_then(|v| v.half_width.as_deref()).unwrap_or("null"), row.changed(), row.classification());
    }
    out.push_str("],\"materials\":{\"pack_set_root_left\":"); push_json_string(&mut out, &comparison.pack_set_root_left);
    out.push_str(",\"pack_set_root_right\":"); push_json_string(&mut out, &comparison.pack_set_root_right);
    let _ = write!(out, ",\"changed\":{},\"packs\":[", comparison.pack_set_root_left != comparison.pack_set_root_right || comparison.packs.iter().any(PackDiff::changed));
    for (index, row) in comparison.packs.iter().enumerate() {
        if index > 0 { out.push(','); }
        out.push_str("{\"kind\":"); push_json_string(&mut out, &row.kind);
        for (key, value) in [("card_left", &row.card_left), ("card_right", &row.card_right), ("identity_left", &row.identity_left), ("identity_right", &row.identity_right)] {
            let _ = write!(out, ",\"{key}\":"); push_opt(&mut out, value.as_deref());
        }
        let _ = write!(out, ",\"changed\":{}}}", row.changed());
    }
    out.push_str("]},\"stages\":[");
    for (index, row) in comparison.stages.iter().enumerate() {
        if index > 0 { out.push(','); }
        out.push_str("{\"stage\":"); push_json_string(&mut out, row.stage);
        out.push_str(",\"left\":"); push_json_string(&mut out, &row.left);
        out.push_str(",\"right\":"); push_json_string(&mut out, &row.right);
        out.push_str(",\"status\":"); push_json_string(&mut out, row.status());
        out.push_str(",\"differing_keys\":[");
        for (i, key) in row.differing_keys.iter().enumerate() {
            if i > 0 { out.push(','); } push_json_string(&mut out, key);
        }
        out.push_str("]}");
    }
    out.push_str("],\"authority\":\"projection-of-retained-receipts\",\"verification\":");
    push_json_string(&mut out, comparison.verification);
    out.push_str(",\"no_claim\":"); push_json_string(&mut out, NO_CLAIM);
    out.push_str("}\n"); out
}

fn render_text(c: &Comparison, subject: &str) -> String {
    let mut out = format!("status=ok\ncommand=compare\nsubject={}\nleft_run={}\nright_run={}\nproject_hash_left={}\nproject_hash_right={}\nsame_project={}\nchanged={}\nsummary={}\nqoi_count={}\n",
        escape_text(subject), c.left_run, c.right_run, c.project_hash_left, c.project_hash_right, c.same_project(), c.changed(), escape_text(&c.summary()), c.qoi.len());
    for (i, row) in c.qoi.iter().enumerate() {
        let l = row.left.as_ref(); let r = row.right.as_ref();
        let _ = writeln!(out, "qoi.{i}={} {} {} -> {} {} delta={} delta_unit={} color={}->{} classification={} reason={}",
            escape_text(&row.name), l.map_or("null", |v| v.nominal.as_str()), escape_text(l.map_or("absent", |v| v.unit.as_str())),
            r.map_or("null", |v| v.nominal.as_str()), escape_text(r.map_or("absent", |v| v.unit.as_str())), number(row.delta),
            escape_text(row.delta_unit.as_deref().unwrap_or("unavailable")), escape_text(l.map_or("absent", |v| v.color.as_str())),
            escape_text(r.map_or("absent", |v| v.color.as_str())), row.classification(), row.reason.unwrap_or("none"));
    }
    for (i, row) in c.requirements.iter().enumerate() {
        let l = row.left.as_ref(); let r = row.right.as_ref();
        let _ = writeln!(out, "requirement.{i}={} outcome={}->{} nominal_margin_kelvin={}->{} classification={}", escape_text(&row.id),
            escape_text(l.map_or("absent", |v| v.outcome.as_str())), escape_text(r.map_or("absent", |v| v.outcome.as_str())),
            l.map_or("null", |v| v.nominal_margin.as_str()), r.map_or("null", |v| v.nominal_margin.as_str()), row.classification());
    }
    for (i, row) in c.terms.iter().enumerate() {
        let _ = writeln!(out, "budget_term.{i}={} state={}->{} half_width_kelvin={}->{} classification={}", escape_text(&row.kind),
            escape_text(row.left.as_ref().map_or("absent", |v| v.state.as_str())), escape_text(row.right.as_ref().map_or("absent", |v| v.state.as_str())),
            row.left.as_ref().and_then(|v| v.half_width.as_deref()).unwrap_or("null"), row.right.as_ref().and_then(|v| v.half_width.as_deref()).unwrap_or("null"), row.classification());
    }
    let _ = writeln!(out, "materials.pack_set_root={}->{}", c.pack_set_root_left, c.pack_set_root_right);
    for (i, row) in c.packs.iter().enumerate() {
        let _ = writeln!(out, "materials.pack.{i}={} card={}->{} changed={}", escape_text(&row.kind),
            escape_text(row.card_left.as_deref().unwrap_or("absent")), escape_text(row.card_right.as_deref().unwrap_or("absent")), row.changed());
    }
    for row in &c.stages { let _ = writeln!(out, "stage.{}={} {} -> {}", row.stage, row.status(), row.left, row.right); }
    let _ = writeln!(out, "authority=projection-of-retained-receipts\nverification={}\nno_claim={}", c.verification, escape_text(NO_CLAIM));
    out
}

/// Execute the `compare` verb. Both exports are verified before diffing.
#[must_use]
pub fn compare_path(left_run: &str, right_run: &str, ledger_path: Option<&Path>, mode: OutputMode) -> CommandOutput {
    let subject = format!("{left_run}..{right_run}");
    let (ledger, dir) = match open_export_ledger(COMMAND, &subject, ledger_path, mode) { Ok(v) => v, Err(output) => return output };
    let left = match load_export_from(COMMAND, left_run, &ledger, dir.clone(), mode) { Ok(v) => v, Err(output) => return output };
    let right = match load_export_from(COMMAND, right_run, &ledger, dir, mode) { Ok(v) => v, Err(output) => return output };
    let comparison = match build(&ledger, &left, &right) {
        Ok(v) => v,
        Err(BuildError::Shape(why)) => return shape_refusal(mode, &subject, why),
        Err(BuildError::Incomparable(why)) => return refuse(mode, COMMAND, exit::REFUSED,
            "cli-compare-project-mismatch", &subject,
            format!("the runs belong to different projects with no common QoI to diff: {why}; {} solved project {} and {} solved project {}",
                left.export.run, left.export.project_hash, right.export.run, right.export.project_hash),
            "compare runs sharing at least one named QoI; output additions, removals, and changed requirements are supported"),
    };
    let stdout = match mode { OutputMode::Json => render_json(&comparison, &subject), OutputMode::Text => render_text(&comparison, &subject) };
    CommandOutput { exit_code: exit::SUCCESS, stdout, stderr: String::new() }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn receipt(rows: &str) -> JsonValue { JsonValue::parse(&format!("{{\"qoi\":[{rows}],\"requirements\":[]}}" )).unwrap() }
    fn qoi(name: &str, value: &str, unit: &str) -> String {
        format!("{{\"name\":\"{name}\",\"value\":{value},\"unit\":\"{unit}\",\"color\":\"estimated\",\"identity\":\"i\"}}")
    }
    fn comparison(qoi: Vec<QoiDiff>) -> Comparison {
        Comparison { left_run: "a".into(), right_run: "b".into(), project_hash_left: "p".into(), project_hash_right: "p".into(),
            verification: "sealed-evidence", qoi, requirements: vec![], terms: vec![], pack_set_root_left: "m".into(),
            pack_set_root_right: "m".into(), packs: vec![], stages: vec![] }
    }

    #[test]
    fn g0_output_union_retains_additions_removals_and_shared_values() {
        let l = receipt(&format!("{},{}", qoi("common", "300", "K"), qoi("removed", "1", "W")));
        let r = receipt(&format!("{},{}", qoi("added", "2", "W"), qoi("common", "301", "K")));
        assert!(comparability(&l, &r).is_ok());
        let rows = diff_qoi(&l, &r).unwrap();
        assert_eq!(rows.iter().map(|r| (r.name.as_str(), r.classification())).collect::<Vec<_>>(),
            vec![("added", "added"), ("common", "changed"), ("removed", "removed")]);
        assert_eq!(rows[1].delta, Some(1.0));
        assert_eq!(rows[0].delta, None);
        let out = render_json(&comparison(rows), "test");
        assert!(JsonValue::parse(&out).is_ok(), "{out}");
        assert!(out.contains("\"nominal_left\":null"));
        assert!(out.contains("\"nominal_right\":null"));
    }

    #[test]
    fn g0_disjoint_projects_and_duplicate_keys_are_not_silently_joined() {
        let l = receipt(&qoi("temperature", "300", "K"));
        let r = receipt(&qoi("pressure", "300", "Pa"));
        assert!(matches!(comparability(&l, &r), Err(BuildError::Incomparable(_))));
        let duplicate = receipt(&format!("{},{}", qoi("temperature", "300", "K"), qoi("temperature", "301", "K")));
        assert!(matches!(comparability(&l, &duplicate), Err(BuildError::Shape(_))));
        assert!(diff_qoi(&duplicate, &l).is_err());
    }

    #[test]
    fn g3_unit_rescaling_and_affine_temperature_differences() {
        for (l, lu, r, ru, delta, du) in [("1", "kPa", "1500", "Pa", 500.0, "kg·m^-1·s^-2"),
            ("20", "degC", "303.15", "K", 10.0, "K"), ("20", "degC", "30", "degC", 10.0, "K")] {
            let rows = diff_qoi(&receipt(&qoi("q", l, lu)), &receipt(&qoi("q", r, ru))).unwrap();
            assert!((rows[0].delta.unwrap() - delta).abs() < 1e-10);
            assert_eq!(rows[0].delta_unit.as_deref(), Some(du));
            if lu == "degC" { assert_eq!(rows[0].rel_delta, None); }
        }
    }

    #[test]
    fn g0_unavailable_deltas_remain_null_without_losing_other_rows() {
        for (l, lu, r, ru, why) in [("1", "K", "2", "Pa", "incompatible-dimensions"),
            ("1", "unknown", "2", "Pa", "unit-conversion-unavailable"),
            ("-1e308", "K", "1e308", "K", "delta-overflow")] {
            let left = receipt(&format!("{},{}", qoi("bad", l, lu), qoi("good", "3", "W")));
            let right = receipt(&format!("{},{}", qoi("bad", r, ru), qoi("good", "4", "W")));
            let rows = diff_qoi(&left, &right).unwrap();
            assert_eq!(rows[0].classification(), "incomparable");
            assert_eq!(rows[0].reason, Some(why));
            assert_eq!(rows[0].delta, None);
            assert_eq!(rows[1].delta, Some(1.0));
            let out = render_json(&comparison(rows), "test");
            assert!(JsonValue::parse(&out).is_ok(), "{out}");
            assert!(out.contains("\"delta\":null,\"rel_delta\":null"));
        }
    }

    #[test]
    fn g0_zero_baseline_and_relative_overflow_are_not_reported_as_zero() {
        for (l, r) in [("0", "1"), ("1e-300", "1e300")] {
            let rows = diff_qoi(&receipt(&qoi("q", l, "W")), &receipt(&qoi("q", r, "W"))).unwrap();
            assert!(rows[0].delta.unwrap() > 0.0);
            assert_eq!(rows[0].rel_delta, None);
        }
    }

    #[test]
    fn g0_numeric_spelling_is_not_a_physical_change() {
        let rows = diff_qoi(&receipt(&qoi("q", "1", "W")), &receipt(&qoi("q", "1.0", "W"))).unwrap();
        assert_eq!(rows[0].classification(), "same");
        assert_eq!(rows[0].delta, Some(0.0));
        assert_eq!(rows[0].rel_delta, Some(0.0));
        // A non-finite spelling never reaches the diff: the receipt reader refuses it.
        assert!(JsonValue::parse(&format!("{{\"qoi\":[{}],\"requirements\":[]}}", qoi("q", "1e999", "W"))).is_err());
    }

    #[test]
    fn g0_changed_requirement_set_is_a_union_not_a_project_refusal() {
        let parse = |id: &str| JsonValue::parse(&format!("{{\"qoi\":[{}],\"requirements\":[{{\"id\":\"{id}\",\"outcome\":\"indeterminate\",\"nominal_margin_kelvin\":2,\"effective_limit_kelvin\":310}}]}}", qoi("q", "308", "K"))).unwrap();
        let l = parse("old"); let r = parse("new");
        assert!(comparability(&l, &r).is_ok());
        let rows = diff_requirements(&l, &r).unwrap();
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0].classification(), "added");
        assert_eq!(rows[1].classification(), "removed");
        assert!(rows.iter().all(|r| r.changed() && !r.verdict_changed()));
    }

    #[test]
    fn g0_budget_terms_retain_additions_and_removals() {
        let parse = |kind: &str| JsonValue::parse(&format!("{{\"budget\":[{{\"terms\":[{{\"kind\":\"{kind}\",\"state\":\"no-data\"}}]}}]}}" )).unwrap();
        let rows = diff_terms(&parse("old"), &parse("new")).unwrap();
        assert_eq!(rows[0].classification(), "added");
        assert_eq!(rows[1].classification(), "removed");
    }

    #[test]
    fn g3_swapping_operands_reverses_compatible_deltas() {
        let l = receipt(&qoi("q", "1", "kPa")); let r = receipt(&qoi("q", "1500", "Pa"));
        let forward = diff_qoi(&l, &r).unwrap(); let reverse = diff_qoi(&r, &l).unwrap();
        assert_eq!(forward[0].delta, reverse[0].delta.map(|v| -v));
        assert_eq!(forward[0].delta_unit, reverse[0].delta_unit);
    }

    #[test]
    fn differing_keys_names_exactly_the_keys_whose_values_differ() {
        let left = JsonValue::parse("{\"schema\":\"s\",\"run\":\"a\",\"verified\":[1,2]}").unwrap();
        let right = JsonValue::parse("{\"schema\":\"s\",\"run\":\"b\",\"verified\":[1,2]}").unwrap();
        assert!(differing_keys(&left, &left).is_empty());
        assert_eq!(differing_keys(&left, &right), vec!["run".to_string()]);
        let extra = JsonValue::parse("{\"schema\":\"s\",\"run\":\"a\",\"verified\":[1,3],\"x\":1}").unwrap();
        assert_eq!(differing_keys(&left, &extra), vec!["verified".to_string(), "x".to_string()]);
    }

    #[test]
    fn stage_status_treats_only_binding_keys_as_unchanged_inputs() {
        let mut stage = StageDiff { stage: "qoi", left: "a".into(), right: "b".into(), differing_keys: vec![] };
        assert_eq!(stage.status(), "unchanged (same receipt)");
        stage.differing_keys = BINDING_KEYS.iter().map(|s| (*s).to_string()).collect();
        assert!(!stage.changed());
        stage.differing_keys.push("value".into());
        assert!(stage.changed());
    }
}
