//! `compare` verb: diff two completed runs of the same project by their
//! retained receipts (bead frankensim-rc-root-q61wp.47).
//!
//! Like `report` and `package`, this verb never replays physics. Both runs
//! are located through the resume loader, their report receipts are bound to
//! the retained artifacts, and every row below is read from a retained stage
//! receipt: the QoI receipt (values, colours, requirement outcomes, budget
//! terms), the material-resolve receipt (card identities), and the seven
//! stage receipt hashes. A row is "changed" only when the two receipts
//! disagree; a stage whose receipt differs solely in binding keys (`run`,
//! `project_hash`, `import_op`: which run, project, and import op the receipt
//! belongs to, not what the stage computed) is reported as unchanged in its
//! inputs, because that is what the receipt bytes say, not because the verb
//! inferred it.
//!
//! A design change (another card pack, fidelity, or geometry import) is a
//! different canonical project, so two runs with different project hashes are
//! compared and the hash change is one of the rows. Two runs are refused as
//! different projects, naming both hashes, only when they have no common
//! requirement to diff against: their QoI names or requirement identities
//! differ.

use std::fmt::Write as _;
use std::path::Path;

use fs_blake3::ContentHash;
use fs_ledger::Ledger;

use crate::json_read::JsonValue;
use crate::report::{LoadedExport, load_export_from, open_export_ledger, refuse};
use crate::{CommandOutput, OutputMode, RESULT_SCHEMA, escape_text, exit, push_json_string};

const COMMAND: &str = "compare";
const NO_CLAIM: &str = "the comparison projects retained stage receipts of two completed runs and adds no physical, numerical, or validation authority; a changed value is a change between two Estimated candidates, not a verified effect";

/// One QoI row across the two runs.
struct QoiDiff {
    name: String,
    unit_left: String,
    unit_right: String,
    nominal_left: String,
    nominal_right: String,
    delta: f64,
    rel_delta: f64,
    color_left: String,
    color_right: String,
    identity_same: bool,
}

impl QoiDiff {
    /// The physical row changed: value, unit, or colour. Identity is reported
    /// separately because it carries lineage (project and run binding), and a
    /// bit-identical value under a renamed project is not a changed result.
    fn changed(&self) -> bool {
        self.nominal_left != self.nominal_right
            || self.color_left != self.color_right
            || self.unit_left != self.unit_right
    }
}

/// One requirement row across the two runs.
struct RequirementDiff {
    id: String,
    outcome_left: String,
    outcome_right: String,
    nominal_margin_left: String,
    nominal_margin_right: String,
    effective_limit_left: String,
    effective_limit_right: String,
}

impl RequirementDiff {
    fn verdict_changed(&self) -> bool {
        self.outcome_left != self.outcome_right
    }
    fn changed(&self) -> bool {
        self.verdict_changed()
            || self.nominal_margin_left != self.nominal_margin_right
            || self.effective_limit_left != self.effective_limit_right
    }
}

/// One engineering-uncertainty budget term across the two runs.
struct TermDiff {
    kind: String,
    state_left: String,
    state_right: String,
    /// Upper half-width in kelvin when the term is an interval; `None` when
    /// the producer recorded NO-DATA.
    half_width_left: Option<String>,
    half_width_right: Option<String>,
}

impl TermDiff {
    fn changed(&self) -> bool {
        self.state_left != self.state_right || self.half_width_left != self.half_width_right
    }
}

/// One admitted card pack across the two runs, matched by pack kind.
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

/// One solve stage across the two runs.
struct StageDiff {
    stage: &'static str,
    left: String,
    right: String,
    /// Top-level receipt keys whose values differ (empty when the receipts
    /// are byte-identical).
    differing_keys: Vec<String>,
}

/// Receipt keys that name which run, project, and import op a receipt
/// belongs to, not what the stage computed. Two receipts differing only in
/// these keys saw the same inputs.
const BINDING_KEYS: [&str; 3] = ["run", "project_hash", "import_op"];

impl StageDiff {
    fn status(&self) -> &'static str {
        if self.differing_keys.is_empty() {
            "unchanged (same receipt)"
        } else if self
            .differing_keys
            .iter()
            .all(|key| BINDING_KEYS.contains(&key.as_str()))
        {
            "unchanged (same inputs; differs only by binding keys)"
        } else {
            "changed"
        }
    }
    fn changed(&self) -> bool {
        self.status() == "changed"
    }
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
    fn same_project(&self) -> bool {
        self.project_hash_left == self.project_hash_right
    }

    fn changed(&self) -> bool {
        !self.same_project()
            || self.qoi.iter().any(QoiDiff::changed)
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
        let mut parts = Vec::new();
        let changed_stages = self.stages.iter().filter(|stage| stage.changed()).count();
        parts.push(format!(
            "{changed_stages} of {} stages changed",
            self.stages.len()
        ));
        if !self.same_project() {
            parts.push(format!(
                "project hash {} -> {}",
                self.project_hash_left, self.project_hash_right
            ));
        }
        for qoi in self.qoi.iter().filter(|qoi| qoi.changed()) {
            parts.push(format!(
                "{} {} -> {} {}",
                qoi.name, qoi.nominal_left, qoi.nominal_right, qoi.unit_right
            ));
        }
        for requirement in &self.requirements {
            if requirement.verdict_changed() {
                parts.push(format!(
                    "verdict {} {} -> {}",
                    requirement.id, requirement.outcome_left, requirement.outcome_right
                ));
            } else if requirement.changed() {
                parts.push(format!(
                    "verdict {} stays {} (margin {} -> {})",
                    requirement.id,
                    requirement.outcome_left,
                    requirement.nominal_margin_left,
                    requirement.nominal_margin_right
                ));
            }
        }
        for term in self.terms.iter().filter(|term| term.changed()) {
            parts.push(format!(
                "budget term {} {} -> {}",
                term.kind, term.state_left, term.state_right
            ));
        }
        for pack in self.packs.iter().filter(|pack| pack.changed()) {
            parts.push(format!(
                "{} card {} -> {}",
                pack.kind,
                pack.card_left.as_deref().unwrap_or("absent"),
                pack.card_right.as_deref().unwrap_or("absent")
            ));
        }
        parts.join("; ")
    }
}

fn shape_refusal(mode: OutputMode, subject: &str, why: impl Into<String>) -> CommandOutput {
    refuse(
        mode,
        COMMAND,
        exit::REFUSED,
        "cli-compare-receipt-shape",
        subject,
        why,
        "regenerate the run; compare reads retained receipts and never repairs one",
    )
}

/// Read the retained receipt of `stage` for a loaded run as parsed JSON.
fn stage_receipt(
    ledger: &Ledger,
    loaded: &LoadedExport,
    stage: &str,
) -> Result<(String, JsonValue), String> {
    let mut hashes = loaded
        .export
        .stages
        .iter()
        .filter(|(name, _, _)| *name == stage)
        .map(|(_, _, hash)| hash.as_str());
    let Some(hash) = hashes.next() else {
        return Err(format!(
            "run {} has no `{stage}` stage receipt",
            loaded.export.run
        ));
    };
    if hashes.next().is_some() {
        return Err(format!(
            "run {} has multiple `{stage}` stage receipts",
            loaded.export.run
        ));
    }
    let content = ContentHash::from_hex(hash)
        .ok_or_else(|| format!("`{stage}` receipt hash `{hash}` is not a content hash"))?;
    let bytes = ledger
        .get_artifact(&content)
        .map_err(|error| format!("reading the `{stage}` receipt failed: {error}"))?
        .ok_or_else(|| format!("the `{stage}` receipt {hash} is not retained in this ledger"))?;
    let text = String::from_utf8(bytes)
        .map_err(|_| format!("the `{stage}` receipt bytes are not UTF-8"))?;
    let value = JsonValue::parse(&text)
        .map_err(|error| format!("the `{stage}` receipt does not parse: {error}"))?;
    Ok((hash.to_string(), value))
}

fn required_str(value: &JsonValue, key: &str, context: &str) -> Result<String, String> {
    value
        .str_field(key)
        .map(str::to_string)
        .ok_or_else(|| format!("{context} is missing string field `{key}`"))
}

fn required_number(value: &JsonValue, key: &str, context: &str) -> Result<String, String> {
    value
        .get(key)
        .and_then(JsonValue::number_raw)
        .map(str::to_string)
        .ok_or_else(|| format!("{context} is missing numeric field `{key}`"))
}

fn required_array<'a>(
    value: &'a JsonValue,
    key: &str,
    context: &str,
) -> Result<&'a [JsonValue], String> {
    value
        .get(key)
        .and_then(JsonValue::as_array)
        .ok_or_else(|| format!("{context} is missing array field `{key}`"))
}

fn finite_or_zero(value: f64) -> f64 {
    if value.is_finite() { value } else { 0.0 }
}

fn diff_qoi(left: &JsonValue, right: &JsonValue) -> Result<Vec<QoiDiff>, String> {
    let rows_left = required_array(left, "qoi", "left QoI receipt")?;
    let rows_right = required_array(right, "qoi", "right QoI receipt")?;
    let mut rows = Vec::with_capacity(rows_left.len());
    for row_left in rows_left {
        let name = required_str(row_left, "name", "left QoI row")?;
        let row_right = rows_right
            .iter()
            .find(|row| row.str_field("name") == Some(name.as_str()))
            .ok_or_else(|| format!("right run has no QoI named `{name}`"))?;
        let nominal_left = required_number(row_left, "value", "left QoI row")?;
        let nominal_right = required_number(row_right, "value", "right QoI row")?;
        let value_left = row_left.f64_field("value").unwrap_or(f64::NAN);
        let value_right = row_right.f64_field("value").unwrap_or(f64::NAN);
        let delta = finite_or_zero(value_right - value_left);
        let rel_delta = if value_left == 0.0 {
            0.0
        } else {
            finite_or_zero(delta / value_left.abs())
        };
        rows.push(QoiDiff {
            unit_left: required_str(row_left, "unit", "left QoI row")?,
            unit_right: required_str(row_right, "unit", "right QoI row")?,
            color_left: required_str(row_left, "color", "left QoI row")?,
            color_right: required_str(row_right, "color", "right QoI row")?,
            identity_same: required_str(row_left, "identity", "left QoI row")?
                == required_str(row_right, "identity", "right QoI row")?,
            name,
            nominal_left,
            nominal_right,
            delta,
            rel_delta,
        });
    }
    if rows_right.len() != rows_left.len() {
        return Err(format!(
            "the runs retain different QoI counts ({} vs {})",
            rows_left.len(),
            rows_right.len()
        ));
    }
    Ok(rows)
}

fn diff_requirements(left: &JsonValue, right: &JsonValue) -> Result<Vec<RequirementDiff>, String> {
    let rows_left = required_array(left, "requirements", "left QoI receipt")?;
    let rows_right = required_array(right, "requirements", "right QoI receipt")?;
    let mut rows = Vec::with_capacity(rows_left.len());
    for row_left in rows_left {
        let id = required_str(row_left, "id", "left requirement row")?;
        let row_right = rows_right
            .iter()
            .find(|row| row.str_field("id") == Some(id.as_str()))
            .ok_or_else(|| format!("right run has no requirement `{id}`"))?;
        rows.push(RequirementDiff {
            outcome_left: required_str(row_left, "outcome", "left requirement row")?,
            outcome_right: required_str(row_right, "outcome", "right requirement row")?,
            nominal_margin_left: required_number(
                row_left,
                "nominal_margin_kelvin",
                "left requirement row",
            )?,
            nominal_margin_right: required_number(
                row_right,
                "nominal_margin_kelvin",
                "right requirement row",
            )?,
            effective_limit_left: required_number(
                row_left,
                "effective_limit_kelvin",
                "left requirement row",
            )?,
            effective_limit_right: required_number(
                row_right,
                "effective_limit_kelvin",
                "right requirement row",
            )?,
            id,
        });
    }
    Ok(rows)
}

fn budget_terms<'a>(receipt: &'a JsonValue, side: &str) -> Result<&'a [JsonValue], String> {
    let budgets = required_array(receipt, "budget", &format!("{side} QoI receipt"))?;
    let budget = budgets
        .first()
        .ok_or_else(|| format!("{side} QoI receipt carries no budget"))?;
    required_array(budget, "terms", &format!("{side} budget"))
}

fn diff_terms(left: &JsonValue, right: &JsonValue) -> Result<Vec<TermDiff>, String> {
    let terms_left = budget_terms(left, "left")?;
    let terms_right = budget_terms(right, "right")?;
    let mut rows = Vec::with_capacity(terms_left.len());
    for term_left in terms_left {
        let kind = required_str(term_left, "kind", "left budget term")?;
        let term_right = terms_right
            .iter()
            .find(|term| term.str_field("kind") == Some(kind.as_str()))
            .ok_or_else(|| format!("right run has no budget term `{kind}`"))?;
        let half_width = |term: &JsonValue| {
            term.get("upper_kelvin")
                .and_then(JsonValue::number_raw)
                .map(str::to_string)
        };
        rows.push(TermDiff {
            state_left: required_str(term_left, "state", "left budget term")?,
            state_right: required_str(term_right, "state", "right budget term")?,
            half_width_left: half_width(term_left),
            half_width_right: half_width(term_right),
            kind,
        });
    }
    Ok(rows)
}

fn diff_packs(left: &JsonValue, right: &JsonValue) -> Result<Vec<PackDiff>, String> {
    let packs_left = required_array(left, "packs", "left material receipt")?;
    let packs_right = required_array(right, "packs", "right material receipt")?;
    let mut kinds: Vec<String> = Vec::new();
    for pack in packs_left.iter().chain(packs_right) {
        let kind = required_str(pack, "kind", "card pack row")?;
        if !kinds.contains(&kind) {
            kinds.push(kind);
        }
    }
    let find =
        |packs: &[JsonValue], kind: &str| -> Result<(Option<String>, Option<String>), String> {
            let Some(pack) = packs
                .iter()
                .find(|pack| pack.str_field("kind") == Some(kind))
            else {
                return Ok((None, None));
            };
            Ok((
                Some(required_str(pack, "card", "card pack row")?),
                Some(required_str(pack, "identity", "card pack row")?),
            ))
        };
    let mut rows = Vec::with_capacity(kinds.len());
    for kind in kinds {
        let (card_left, identity_left) = find(packs_left, &kind)?;
        let (card_right, identity_right) = find(packs_right, &kind)?;
        rows.push(PackDiff {
            kind,
            card_left,
            card_right,
            identity_left,
            identity_right,
        });
    }
    Ok(rows)
}

fn differing_keys(left: &JsonValue, right: &JsonValue) -> Vec<String> {
    if left == right {
        return Vec::new();
    }
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
        if !left.iter().any(|(other, _)| other == key) {
            keys.push(key.clone());
        }
    }
    keys
}

/// Why two loaded runs could not be compared.
enum BuildError {
    /// The runs measure different things: no common QoI or requirement to
    /// diff against. Carries the reason.
    Incomparable(String),
    /// A retained receipt could not be read in the expected shape.
    Shape(String),
}

impl From<String> for BuildError {
    fn from(why: String) -> Self {
        BuildError::Shape(why)
    }
}

/// Sorted string fields of an array's rows.
fn sorted_fields(
    value: &JsonValue,
    array: &str,
    key: &str,
    side: &str,
) -> Result<Vec<String>, String> {
    let mut names: Vec<String> = required_array(value, array, &format!("{side} QoI receipt"))?
        .iter()
        .map(|row| required_str(row, key, &format!("{side} {array} row")))
        .collect::<Result<_, _>>()?;
    names.sort();
    Ok(names)
}

/// Two runs are comparable when they measure the same QoIs against the same
/// requirement identities (the identity names the QoI, class, region,
/// severity, and sources, not the limit, so a changed limit is a comparable
/// design change).
fn comparability(qoi_left: &JsonValue, qoi_right: &JsonValue) -> Result<(), BuildError> {
    let names_left = sorted_fields(qoi_left, "qoi", "name", "left")?;
    let names_right = sorted_fields(qoi_right, "qoi", "name", "right")?;
    if names_left != names_right {
        return Err(BuildError::Incomparable(format!(
            "the runs measure different QoIs ({} vs {})",
            names_left.join(","),
            names_right.join(",")
        )));
    }
    let ids_left = sorted_fields(qoi_left, "requirements", "id", "left")?;
    let ids_right = sorted_fields(qoi_right, "requirements", "id", "right")?;
    if ids_left != ids_right {
        return Err(BuildError::Incomparable(format!(
            "the runs evaluate different requirement identities ({} vs {})",
            ids_left.join(","),
            ids_right.join(",")
        )));
    }
    Ok(())
}

fn build(
    ledger: &Ledger,
    left: &LoadedExport,
    right: &LoadedExport,
) -> Result<Comparison, BuildError> {
    let (_, qoi_left) = stage_receipt(ledger, left, "qoi")?;
    let (_, qoi_right) = stage_receipt(ledger, right, "qoi")?;
    comparability(&qoi_left, &qoi_right)?;
    let (_, material_left) = stage_receipt(ledger, left, "material-resolve")?;
    let (_, material_right) = stage_receipt(ledger, right, "material-resolve")?;
    let mut stages = Vec::with_capacity(left.export.stages.len());
    for (stage, _, hash_left) in &left.export.stages {
        let (_, receipt_left) = stage_receipt(ledger, left, stage)?;
        let (hash_right, receipt_right) = stage_receipt(ledger, right, stage)?;
        stages.push(StageDiff {
            stage,
            left: hash_left.clone(),
            right: hash_right,
            differing_keys: differing_keys(&receipt_left, &receipt_right),
        });
    }
    if right.export.stages.len() != left.export.stages.len() {
        return Err(BuildError::Shape(format!(
            "the runs completed different stage counts ({} vs {})",
            left.export.stages.len(),
            right.export.stages.len()
        )));
    }
    Ok(Comparison {
        left_run: left.export.run.clone(),
        right_run: right.export.run.clone(),
        project_hash_left: left.export.project_hash.clone(),
        project_hash_right: right.export.project_hash.clone(),
        verification: left.export.verification,
        qoi: diff_qoi(&qoi_left, &qoi_right)?,
        requirements: diff_requirements(&qoi_left, &qoi_right)?,
        terms: diff_terms(&qoi_left, &qoi_right)?,
        pack_set_root_left: required_str(&material_left, "pack_set_root", "left material receipt")?,
        pack_set_root_right: required_str(
            &material_right,
            "pack_set_root",
            "right material receipt",
        )?,
        packs: diff_packs(&material_left, &material_right)?,
        stages,
    })
}

fn push_opt(out: &mut String, value: Option<&str>) {
    match value {
        Some(value) => push_json_string(out, value),
        None => out.push_str("null"),
    }
}

fn render_json(comparison: &Comparison, subject: &str) -> String {
    let changed = comparison.changed();
    let mut out = String::from("{\"schema\":");
    push_json_string(&mut out, RESULT_SCHEMA);
    out.push_str(",\"command\":\"compare\",\"status\":\"ok\",\"subject\":");
    push_json_string(&mut out, subject);
    out.push_str(",\"left_run\":");
    push_json_string(&mut out, &comparison.left_run);
    out.push_str(",\"right_run\":");
    push_json_string(&mut out, &comparison.right_run);
    out.push_str(",\"project_hash_left\":");
    push_json_string(&mut out, &comparison.project_hash_left);
    out.push_str(",\"project_hash_right\":");
    push_json_string(&mut out, &comparison.project_hash_right);
    let _ = write!(
        out,
        ",\"same_project\":{},\"changed\":{changed},\"summary\":",
        comparison.same_project()
    );
    push_json_string(&mut out, &comparison.summary());
    let _ = write!(
        out,
        ",\"qoi_count\":{},\"qoi_diffs\":[",
        comparison.qoi.len()
    );
    for (index, qoi) in comparison.qoi.iter().enumerate() {
        if index > 0 {
            out.push(',');
        }
        out.push_str("{\"name\":");
        push_json_string(&mut out, &qoi.name);
        out.push_str(",\"unit_left\":");
        push_json_string(&mut out, &qoi.unit_left);
        out.push_str(",\"unit_right\":");
        push_json_string(&mut out, &qoi.unit_right);
        let _ = write!(
            out,
            ",\"nominal_left\":{},\"nominal_right\":{},\"delta\":{},\"rel_delta\":{},\"color_left\":",
            qoi.nominal_left, qoi.nominal_right, qoi.delta, qoi.rel_delta
        );
        push_json_string(&mut out, &qoi.color_left);
        out.push_str(",\"color_right\":");
        push_json_string(&mut out, &qoi.color_right);
        let evolution = if qoi.color_left == qoi.color_right {
            "same".to_string()
        } else {
            format!("{}->{}", qoi.color_left, qoi.color_right)
        };
        out.push_str(",\"color_evolution\":");
        push_json_string(&mut out, &evolution);
        let _ = write!(
            out,
            ",\"identity_same\":{},\"classification\":\"{}\"}}",
            qoi.identity_same,
            if qoi.changed() { "changed" } else { "same" }
        );
    }
    out.push_str("],\"requirements\":[");
    for (index, requirement) in comparison.requirements.iter().enumerate() {
        if index > 0 {
            out.push(',');
        }
        out.push_str("{\"id\":");
        push_json_string(&mut out, &requirement.id);
        out.push_str(",\"outcome_left\":");
        push_json_string(&mut out, &requirement.outcome_left);
        out.push_str(",\"outcome_right\":");
        push_json_string(&mut out, &requirement.outcome_right);
        let _ = write!(
            out,
            ",\"verdict_changed\":{},\"nominal_margin_left\":{},\"nominal_margin_right\":{},\"effective_limit_left\":{},\"effective_limit_right\":{}}}",
            requirement.verdict_changed(),
            requirement.nominal_margin_left,
            requirement.nominal_margin_right,
            requirement.effective_limit_left,
            requirement.effective_limit_right
        );
    }
    out.push_str("],\"budget_terms\":[");
    for (index, term) in comparison.terms.iter().enumerate() {
        if index > 0 {
            out.push(',');
        }
        out.push_str("{\"kind\":");
        push_json_string(&mut out, &term.kind);
        out.push_str(",\"state_left\":");
        push_json_string(&mut out, &term.state_left);
        out.push_str(",\"state_right\":");
        push_json_string(&mut out, &term.state_right);
        let _ = write!(
            out,
            ",\"half_width_left_kelvin\":{},\"half_width_right_kelvin\":{},\"changed\":{}}}",
            term.half_width_left.as_deref().unwrap_or("null"),
            term.half_width_right.as_deref().unwrap_or("null"),
            term.changed()
        );
    }
    out.push_str("],\"materials\":{\"pack_set_root_left\":");
    push_json_string(&mut out, &comparison.pack_set_root_left);
    out.push_str(",\"pack_set_root_right\":");
    push_json_string(&mut out, &comparison.pack_set_root_right);
    let _ = write!(
        out,
        ",\"changed\":{},\"packs\":[",
        comparison.pack_set_root_left != comparison.pack_set_root_right
            || comparison.packs.iter().any(PackDiff::changed)
    );
    for (index, pack) in comparison.packs.iter().enumerate() {
        if index > 0 {
            out.push(',');
        }
        out.push_str("{\"kind\":");
        push_json_string(&mut out, &pack.kind);
        out.push_str(",\"card_left\":");
        push_opt(&mut out, pack.card_left.as_deref());
        out.push_str(",\"card_right\":");
        push_opt(&mut out, pack.card_right.as_deref());
        out.push_str(",\"identity_left\":");
        push_opt(&mut out, pack.identity_left.as_deref());
        out.push_str(",\"identity_right\":");
        push_opt(&mut out, pack.identity_right.as_deref());
        let _ = write!(out, ",\"changed\":{}}}", pack.changed());
    }
    out.push_str("]},\"stages\":[");
    for (index, stage) in comparison.stages.iter().enumerate() {
        if index > 0 {
            out.push(',');
        }
        out.push_str("{\"stage\":");
        push_json_string(&mut out, stage.stage);
        out.push_str(",\"left\":");
        push_json_string(&mut out, &stage.left);
        out.push_str(",\"right\":");
        push_json_string(&mut out, &stage.right);
        out.push_str(",\"status\":");
        push_json_string(&mut out, stage.status());
        out.push_str(",\"differing_keys\":[");
        for (key_index, key) in stage.differing_keys.iter().enumerate() {
            if key_index > 0 {
                out.push(',');
            }
            push_json_string(&mut out, key);
        }
        out.push_str("]}");
    }
    out.push_str("],\"authority\":\"projection-of-retained-receipts\",\"verification\":");
    push_json_string(&mut out, comparison.verification);
    out.push_str(",\"no_claim\":");
    push_json_string(&mut out, NO_CLAIM);
    out.push_str("}\n");
    out
}

fn render_text(comparison: &Comparison, subject: &str) -> String {
    let mut out = format!(
        "status=ok\ncommand=compare\nsubject={}\nleft_run={}\nright_run={}\nproject_hash_left={}\nproject_hash_right={}\nsame_project={}\nchanged={}\nsummary={}\nqoi_count={}\n",
        escape_text(subject),
        comparison.left_run,
        comparison.right_run,
        comparison.project_hash_left,
        comparison.project_hash_right,
        comparison.same_project(),
        comparison.changed(),
        escape_text(&comparison.summary()),
        comparison.qoi.len()
    );
    for (index, qoi) in comparison.qoi.iter().enumerate() {
        let _ = writeln!(
            out,
            "qoi.{index}={} {} {} -> {} {} delta={} color={}->{} classification={}",
            escape_text(&qoi.name),
            qoi.nominal_left,
            escape_text(&qoi.unit_left),
            qoi.nominal_right,
            escape_text(&qoi.unit_right),
            qoi.delta,
            escape_text(&qoi.color_left),
            escape_text(&qoi.color_right),
            if qoi.changed() { "changed" } else { "same" }
        );
    }
    for (index, requirement) in comparison.requirements.iter().enumerate() {
        let _ = writeln!(
            out,
            "requirement.{index}={} outcome={}->{} nominal_margin_kelvin={}->{}",
            escape_text(&requirement.id),
            escape_text(&requirement.outcome_left),
            escape_text(&requirement.outcome_right),
            requirement.nominal_margin_left,
            requirement.nominal_margin_right
        );
    }
    for (index, term) in comparison.terms.iter().enumerate() {
        let _ = writeln!(
            out,
            "budget_term.{index}={} state={}->{} half_width_kelvin={}->{}",
            escape_text(&term.kind),
            escape_text(&term.state_left),
            escape_text(&term.state_right),
            term.half_width_left.as_deref().unwrap_or("null"),
            term.half_width_right.as_deref().unwrap_or("null")
        );
    }
    let _ = writeln!(
        out,
        "materials.pack_set_root={}->{}",
        comparison.pack_set_root_left, comparison.pack_set_root_right
    );
    for (index, pack) in comparison.packs.iter().enumerate() {
        let _ = writeln!(
            out,
            "materials.pack.{index}={} card={}->{} changed={}",
            escape_text(&pack.kind),
            pack.card_left.as_deref().unwrap_or("absent"),
            pack.card_right.as_deref().unwrap_or("absent"),
            pack.changed()
        );
    }
    for stage in &comparison.stages {
        let _ = writeln!(
            out,
            "stage.{}={} {} -> {}",
            stage.stage,
            stage.status(),
            stage.left,
            stage.right
        );
    }
    let _ = writeln!(
        out,
        "authority=projection-of-retained-receipts\nverification={}\nno_claim={}",
        comparison.verification,
        escape_text(NO_CLAIM)
    );
    out
}

/// Execute the `compare` verb.
#[must_use]
pub fn compare_path(
    left_run: &str,
    right_run: &str,
    ledger_path: Option<&Path>,
    mode: OutputMode,
) -> CommandOutput {
    let subject = format!("{left_run}..{right_run}");
    let (ledger, dir) = match open_export_ledger(COMMAND, &subject, ledger_path, mode) {
        Ok(opened) => opened,
        Err(output) => return output,
    };
    let left = match load_export_from(COMMAND, left_run, &ledger, dir.clone(), mode) {
        Ok(loaded) => loaded,
        Err(output) => return output,
    };
    let right = match load_export_from(COMMAND, right_run, &ledger, dir, mode) {
        Ok(loaded) => loaded,
        Err(output) => return output,
    };
    let comparison = match build(&ledger, &left, &right) {
        Ok(comparison) => comparison,
        Err(BuildError::Shape(why)) => return shape_refusal(mode, &subject, why),
        Err(BuildError::Incomparable(why)) => {
            return refuse(
                mode,
                COMMAND,
                exit::REFUSED,
                "cli-compare-project-mismatch",
                &subject,
                format!(
                    "the runs belong to different projects with no common requirement to diff: {why}; {} solved project {} and {} solved project {}",
                    left.export.run,
                    left.export.project_hash,
                    right.export.run,
                    right.export.project_hash
                ),
                "compare two runs that measure the same QoIs against the same requirement identities; a changed card pack, limit, fidelity, or geometry import of one project is comparable",
            );
        }
    };
    let stdout = match mode {
        OutputMode::Json => render_json(&comparison, &subject),
        OutputMode::Text => render_text(&comparison, &subject),
    };
    CommandOutput {
        exit_code: exit::SUCCESS,
        stdout,
        stderr: String::new(),
    }
}

#[cfg(test)]
mod tests {
    use super::{BINDING_KEYS, BuildError, JsonValue, StageDiff, comparability, differing_keys};

    fn qoi_receipt(name: &str, requirement_id: &str) -> JsonValue {
        JsonValue::parse(&format!(
            "{{\"qoi\":[{{\"name\":\"{name}\"}}],\"requirements\":[{{\"id\":\"{requirement_id}\"}}]}}"
        ))
        .expect("fixture receipt parses")
    }

    #[test]
    fn comparability_admits_same_qoi_and_requirement_identity() {
        let left = qoi_receipt("temperature-max", "req-a");
        let right = qoi_receipt("temperature-max", "req-a");
        assert!(comparability(&left, &right).is_ok());
    }

    #[test]
    fn comparability_refuses_a_different_qoi_or_requirement_identity() {
        let left = qoi_receipt("temperature-max", "req-a");
        let other_qoi = qoi_receipt("pressure-drop", "req-a");
        let other_requirement = qoi_receipt("temperature-max", "req-b");
        for (right, expected) in [
            (
                other_qoi,
                "different QoIs (temperature-max vs pressure-drop)",
            ),
            (
                other_requirement,
                "different requirement identities (req-a vs req-b)",
            ),
        ] {
            match comparability(&left, &right) {
                Err(BuildError::Incomparable(why)) => assert!(why.contains(expected), "{why}"),
                Err(BuildError::Shape(why)) => panic!("shape error instead of refusal: {why}"),
                Ok(()) => panic!("incomparable runs admitted: {expected}"),
            }
        }
    }

    #[test]
    fn differing_keys_names_exactly_the_keys_whose_values_differ() {
        let left = JsonValue::parse("{\"schema\":\"s\",\"run\":\"a\",\"verified\":[1,2]}")
            .expect("parses");
        let right = JsonValue::parse("{\"schema\":\"s\",\"run\":\"b\",\"verified\":[1,2]}")
            .expect("parses");
        assert_eq!(differing_keys(&left, &left), Vec::<String>::new());
        assert_eq!(differing_keys(&left, &right), vec!["run".to_string()]);
        let extra = JsonValue::parse("{\"schema\":\"s\",\"run\":\"a\",\"verified\":[1,3],\"x\":1}")
            .expect("parses");
        assert_eq!(
            differing_keys(&left, &extra),
            vec!["verified".to_string(), "x".to_string()]
        );
    }

    #[test]
    fn stage_status_treats_only_binding_keys_as_unchanged_inputs() {
        let stage = |keys: &[&str]| StageDiff {
            stage: "import-verify",
            left: "l".to_string(),
            right: "r".to_string(),
            differing_keys: keys.iter().map(|key| (*key).to_string()).collect(),
        };
        assert_eq!(stage(&[]).status(), "unchanged (same receipt)");
        assert_eq!(
            stage(&BINDING_KEYS).status(),
            "unchanged (same inputs; differs only by binding keys)"
        );
        assert_eq!(stage(&["run", "verified"]).status(), "changed");
        assert!(stage(&["run", "verified"]).changed());
        assert!(!stage(&["run"]).changed());
    }
}
