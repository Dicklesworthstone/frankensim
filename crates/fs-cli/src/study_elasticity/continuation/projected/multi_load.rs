//! Independent operating conditions in the native constrained study.
//! The primary scenario is retained with weight one; additional right-edge
//! tractions are separate equilibria, never a summed force or probabilities.
use super::*;
use super::super::super::{number as node_number, pair as node_pair};
use fs_topols::{DesignBoxEdge, RobustAggregate, RobustLoadCase};
use fs_topols::robust_descent::{
    MultiLoadProjectedOptimizer, MultiLoadProjectedProgress, MultiLoadProjectedSettings,
    MultiLoadProjectedStage,
};

#[path = "multi_load/driver.rs"]
mod driver;
pub(super) use driver::drive;

#[derive(Debug, Clone)]
pub(super) struct LoadFamily {
    aggregate: RobustAggregate,
    max_solves: usize,
    max_recovery_solves: usize,
    additional: Vec<RobustLoadCase>,
}

fn aggregate_name(value: RobustAggregate) -> &'static str {
    match value {
        RobustAggregate::WeightedSum => "weighted-sum",
        RobustAggregate::WorstWeightedCase => "worst-weighted-case",
    }
}

pub(super) fn parse(fields: &[Node]) -> Result<Option<LoadFamily>> {
    let Some(pair) = fields.windows(2).find(|pair|
        matches!(&pair[0].kind, NodeKind::Keyword(key) if key == "load-family"))
    else { return Ok(None) };
    let fields = list(&pair[1], "load-family")?;
    if !fields.first().is_some_and(|node|
        matches!(&node.kind, NodeKind::Symbol(value) if value == "independent"))
    { return Err(malformed("load-family must declare independent operating conditions")); }
    let aggregate = match &field(fields, "aggregate")?.kind {
        NodeKind::Symbol(value) if value == "weighted-sum" => RobustAggregate::WeightedSum,
        NodeKind::Symbol(value) if value == "worst-weighted-case" => RobustAggregate::WorstWeightedCase,
        _ => return Err(malformed("load aggregate must be weighted-sum or worst-weighted-case")),
    };
    let max_solves = integer_node(field(fields, "max-solves")?, "max-solves")?;
    let max_recovery_solves = integer_node(field(fields, "max-recovery-solves")?, "max-recovery-solves")?;
    let additional = list(field(fields, "additional")?, "additional loads")?;
    if !(1..=15).contains(&additional.len()) || !(1..=100_000).contains(&max_solves)
        || max_solves < additional.len() + 1 || max_recovery_solves > 100_000
    { return Err(malformed("load-family requires 1..=15 additional cases, a complete baseline allowance within 100000 solves, and 0..=100000 recovery solves")); }
    let additional = additional.iter().map(|node| {
        let fields = list(node, "load case")?;
        if !fields.first().is_some_and(|node|
            matches!(&node.kind, NodeKind::Symbol(value) if value == "case"))
        { return Err(malformed("additional load must be a case declaration")); }
        let band = node_pair(field(fields, "band")?, "band")?;
        let traction = node_pair(field(fields, "traction-pa")?, "traction-pa")?;
        let weight = node_number(field(fields, "weight")?, "weight")?;
        RobustLoadCase::new(DesignBoxEdge::Right, band[0], band[1], traction, weight)
            .map_err(|error| malformed(&error.to_string()))
    }).collect::<Result<Vec<_>>>()?;
    Ok(Some(LoadFamily { aggregate, max_solves, max_recovery_solves, additional }))
}

impl LoadFamily {
    pub(super) fn html(&self, history: &History) -> String {
        format!("<p>{} independent right-edge operating conditions, including the primary scenario with weight one. Compliance columns report the {} objective; weights are not probabilities and forces are never summed. Every case, including zero-weight cases, must meet the unweighted sampled stress limit. Study case solves: {}/{}. Lifetime recovery case solves: {}/{}. Candidate CG and stress sampling are cancellable; baseline construction and two-family checkpoint recovery currently remain synchronous.</p>",
            history.cases.len(), aggregate_name(self.aggregate), history.solves, self.max_solves,
            history.recovery_solves, self.max_recovery_solves)
    }

    pub(super) fn canonical(&self, out: &mut String) {
        let _ = writeln!(out, "    :load-family (independent");
        let _ = writeln!(out, "      :aggregate {}", aggregate_name(self.aggregate));
        let _ = writeln!(out, "      :max-solves {}", self.max_solves);
        let _ = writeln!(out, "      :max-recovery-solves {}", self.max_recovery_solves);
        let _ = writeln!(out, "      :additional (");
        for case in &self.additional {
            let [a, b] = case.interval();
            let [x, y] = case.traction();
            let _ = writeln!(out, "        (case :band ({} {}) :traction-pa ({} {}) :weight {})",
                canonical_float(a), canonical_float(b), canonical_float(x), canonical_float(y),
                canonical_float(case.weight()));
        }
        let _ = writeln!(out, "      ))");
    }

    fn cases(&self, spec: &ElasticitySpec) -> Result<Vec<RobustLoadCase>> {
        let primary = fixture(spec);
        let mut cases = vec![RobustLoadCase::new(DesignBoxEdge::Right,
            0.5 - primary.band, 0.5 + primary.band, [0.0, -primary.load], 1.0)
            .map_err(|error| malformed(&error.to_string()))?];
        cases.extend_from_slice(&self.additional);
        Ok(cases)
    }

    fn controls(&self, policy: &Controls) -> MultiLoadProjectedSettings {
        MultiLoadProjectedSettings {
            max_candidates: policy.search.max_candidates, contraction: policy.search.contraction,
            min_relative_improvement: policy.search.min_relative_improvement,
            max_solves: self.max_solves,
        }
    }
}

fn case_json(case: &RobustLoadCase) -> String {
    let [a, b] = case.interval();
    let [x, y] = case.traction();
    format!("{{\"edge\":\"right\",\"band\":[{a:.17e},{b:.17e}],\"traction_pa\":[{x:.17e},{y:.17e}],\"weight\":{:.17e}}}", case.weight())
}
fn cases_json(cases: &[RobustLoadCase]) -> String {
    cases.iter().map(case_json).collect::<Vec<_>>().join(",")
}
fn states_json(states: &[SampledStressEvaluation]) -> String {
    format!("[{}]", states.iter().map(stress_json).collect::<Vec<_>>().join(","))
}
fn case_states(state: &fs_topols::RobustSampledStressEvaluation) -> Vec<SampledStressEvaluation> {
    (0..state.case_compliances.len()).map(|i| SampledStressEvaluation {
        compliance: state.case_compliances[i], volume: state.volume,
        sampled_max_von_mises: state.case_sampled_max_von_mises[i],
        max_location: state.case_max_locations[i], sample_count: state.case_sample_counts[i],
        snapshot: state.snapshot,
    }).collect()
}

// Presentation of one complete load family, not another PDE or stress sampler.
// Aggregate compliance and unweighted worst stress intentionally have different
// governing cases. A zero objective weight never removes a stress constraint.
fn summary(states: &[SampledStressEvaluation], cases: &[RobustLoadCase], aggregate: RobustAggregate)
    -> Result<SampledStressEvaluation>
{
    if states.is_empty() || states.len() != cases.len() {
        return Err(malformed("incomplete independent load family"));
    }
    let mut result = states[0].clone();
    result.compliance = 0.0;
    result.sample_count = 0;
    for (state, case) in states.iter().zip(cases) {
        if state.snapshot != result.snapshot || state.volume.to_bits() != result.volume.to_bits() {
            return Err(malformed("load cases do not describe the same exact geometry"));
        }
        let value = case.weight() * state.compliance;
        result.compliance = match aggregate {
            RobustAggregate::WeightedSum => result.compliance + value,
            RobustAggregate::WorstWeightedCase => result.compliance.max(value),
        };
        result.sample_count = result.sample_count.checked_add(state.sample_count)
            .ok_or_else(|| malformed("load-family sample count overflow"))?;
        if state.sampled_max_von_mises > result.sampled_max_von_mises {
            result.sampled_max_von_mises = state.sampled_max_von_mises;
            result.max_location = state.max_location;
        }
    }
    if !result.compliance.is_finite() { return Err(malformed("non-finite load aggregate")); }
    Ok(result)
}

fn hex(bytes: &[u8]) -> String {
    let mut text = String::with_capacity(2 * bytes.len());
    for byte in bytes { let _ = write!(text, "{byte:02x}"); }
    text
}
fn unhex(text: &str) -> Result<Vec<u8>> {
    // Native level <=5, 16 cases and 32 updates fit well inside this bound.
    // This is the owner's checkpoint, not a second optimizer wire format.
    if text.is_empty() || text.len() > 512 * 1024 || text.len() % 2 != 0
        || !text.bytes().all(|v| v.is_ascii_digit() || (b'a'..=b'f').contains(&v))
    { return Err(malformed("invalid or oversized multi-load checkpoint encoding")); }
    text.as_bytes().chunks_exact(2).map(|pair| {
        let nibble = |v| if v <= b'9' { v - b'0' } else { v - b'a' + 10 };
        Ok(16 * nibble(pair[0]) + nibble(pair[1]))
    }).collect()
}

pub(super) struct History {
    cases: Vec<RobustLoadCase>,
    baseline: Vec<SampledStressEvaluation>,
    accepted: Vec<Vec<SampledStressEvaluation>>,
    checkpoint: Vec<u8>,
    solves: usize,
    recovery_solves: usize,
}

impl History {
    fn capture(&mut self, owner: &MultiLoadProjectedOptimizer) {
        self.checkpoint = owner.checkpoint_bytes();
        self.solves = owner.solves_started();
    }

    pub(super) fn json_field(&self, family: &LoadFamily) -> String {
        format!(concat!(",\"load_family\":{{\"schema\":\"independent-loads-v1\",",
            "\"startup_and_recovery\":\"synchronous\",\"aggregate\":{},\"cases\":[{}],\"baseline_cases\":{},\"accepted_cases\":[{}],",
            "\"solves_started\":{},\"max_solves\":{},\"recovery_solves_used\":{},",
            "\"max_recovery_solves\":{},\"checkpoint_hex\":{}}}"),
            quoted(aggregate_name(family.aggregate)), cases_json(&self.cases), states_json(&self.baseline),
            self.accepted.iter().map(|row| states_json(row)).collect::<Vec<_>>().join(","),
            self.solves, family.max_solves, self.recovery_solves, family.max_recovery_solves,
            quoted(&hex(&self.checkpoint)))
    }

    pub(super) fn read(value: &JsonValue, baseline: &SampledStressEvaluation,
        accepted: &[SampledStressEvaluation], policy: &Controls) -> Result<Option<Self>>
    {
        let Some(family) = &policy.family else {
            if value.get("load_family").is_some() { return Err(malformed("undeclared load family")); }
            return Ok(None);
        };
        let value = value.get("load_family").ok_or_else(|| malformed("missing load-family history"))?;
        if value.str_field("schema") != Some("independent-loads-v1")
            || value.str_field("startup_and_recovery") != Some("synchronous")
            || value.str_field("aggregate") != Some(aggregate_name(family.aggregate))
            || integer(value, "max_solves")? != family.max_solves
            || integer(value, "max_recovery_solves")? != family.max_recovery_solves
        { return Err(malformed("retained load-family policy differs from the source")); }
        let array = |key| value.get(key).and_then(JsonValue::as_array)
            .ok_or_else(|| malformed("missing load-family array"));
        let rows = array("cases")?;
        if rows.len() != family.additional.len() + 1 { return Err(malformed("load-family case count changed")); }
        let cases = rows.iter().map(|value| {
            let pair = |key| -> Result<[f64; 2]> {
                let row = value.get(key).and_then(JsonValue::as_array)
                    .filter(|row| row.len() == 2).ok_or_else(|| malformed("invalid load-case pair"))?;
                Ok([row[0].as_f64().ok_or_else(|| malformed("nonnumeric load case"))?,
                    row[1].as_f64().ok_or_else(|| malformed("nonnumeric load case"))?])
            };
            if value.str_field("edge") != Some("right") { return Err(malformed("native load case must use the right edge")); }
            let band = pair("band")?;
            RobustLoadCase::new(DesignBoxEdge::Right, band[0], band[1], pair("traction_pa")?,
                value.f64_field("weight").ok_or_else(|| malformed("missing load weight"))?)
                .map_err(|error| malformed(&error.to_string()))
        }).collect::<Result<Vec<_>>>()?;
        if cases_json(&cases[1..]) != cases_json(&family.additional) || cases[0].weight().to_bits() != 1.0f64.to_bits() {
            return Err(malformed("additional loads or primary weight changed"));
        }
        let read_row = |row: &[JsonValue], expected: &SampledStressEvaluation| -> Result<Vec<SampledStressEvaluation>> {
            if row.len() != cases.len() { return Err(malformed("partially measured load family")); }
            let states = row.iter().map(|value| read_stress(value, policy)).collect::<Result<Vec<_>>>()?;
            if !same(&summary(&states, &cases, family.aggregate)?, expected) {
                return Err(malformed("load-family measurements disagree with the retained aggregate"));
            }
            Ok(states)
        };
        let baseline = read_row(array("baseline_cases")?, baseline)?;
        let rows = array("accepted_cases")?;
        if rows.len() != accepted.len() { return Err(malformed("load-family accepted count changed")); }
        let accepted = rows.iter().zip(accepted).map(|(row, expected)| {
            read_row(row.as_array().ok_or_else(|| malformed("invalid accepted family"))?, expected)
        }).collect::<Result<Vec<_>>>()?;
        let solves = integer(value, "solves_started")?;
        let recovery_solves = integer(value, "recovery_solves_used")?;
        if solves < cases.len() * (accepted.len() + 1) || solves > family.max_solves
            || recovery_solves > family.max_recovery_solves
        { return Err(malformed("invalid retained independent-solve accounting")); }
        let checkpoint = unhex(value.str_field("checkpoint_hex").ok_or_else(|| malformed("missing multi-load checkpoint"))?)?;
        Ok(Some(Self { cases, baseline, accepted, checkpoint, solves, recovery_solves }))
    }
}

// The library checks byte integrity, feasibility and exact numerical replay.
// This adapter additionally binds its immutable declaration to the native source.
fn bind(owner: &MultiLoadProjectedOptimizer, spec: &ElasticitySpec, policy: &Controls,
    family: &LoadFamily, fixed: &[(usize, f64)]) -> Result<()>
{
    let a = owner.settings();
    let b = settings(spec, spec.steps);
    let x = owner.controls();
    let y = family.controls(policy);
    let p = owner.projection_settings();
    let q = policy.area;
    let limit = owner.stress_limit().ok_or_else(|| malformed("missing family stress policy"))?;
    let floats_a = [a.volfrac, a.band_cells, a.move_cells, a.ell0, a.mu_al, a.sobolev_alpha,
        a.hole_radius_cells, a.youngs, a.poisson, x.contraction, x.min_relative_improvement,
        p.target, p.tolerance, p.max_shift, limit.max_von_mises, limit.absolute_tolerance];
    let floats_b = [b.volfrac, b.band_cells, b.move_cells, b.ell0, b.mu_al, b.sobolev_alpha,
        b.hole_radius_cells, b.youngs, b.poisson, y.contraction, y.min_relative_improvement,
        q.target, q.tolerance, q.max_shift, policy.stress.max_von_mises, policy.stress.absolute_tolerance];
    if a.level != b.level || a.iterations != b.iterations || a.nucleation_period != b.nucleation_period
        || x.max_candidates != y.max_candidates || x.max_solves != y.max_solves
        || p.max_evaluations != q.max_evaluations || owner.aggregate() != family.aggregate
        || owner.stress_restoration_reduction().is_some()
        || floats_a.into_iter().zip(floats_b).any(|(a, b)| a.to_bits() != b.to_bits())
        || cases_json(owner.load_cases()) != cases_json(&family.cases(spec)?)
        || owner.fixed_nodes().len() != fixed.len()
        || owner.fixed_nodes().iter().zip(fixed).any(|((i, a), (j, b))| i != j || a.to_bits() != b.to_bits())
    { return Err(malformed("multi-load checkpoint does not implement the original native study")); }
    Ok(())
}

#[cfg(test)]
#[path = "multi_load/tests.rs"]
mod tests;
