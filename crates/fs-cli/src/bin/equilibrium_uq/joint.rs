//! One joint event on a complete shared draw, never a product of marginal rates.
use super::*;
use fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::objective::design::constraints::{
    ConstraintSense, ResponseConstraint, ResponseQuantity,
};
use fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::objective::design::forward::ForwardConstraintResult;

pub(super) struct Event {
    requirements: Vec<ResponseConstraint>,
    tolerances: Vec<Option<f64>>,
    failures: Vec<usize>,
    ranges: Vec<[f64; 2]>,
    samples: usize,
}
impl Event {
    pub(super) fn new(problem: &EquilibriumDesign, tolerances: &[(String, f64)]) -> Result<Self, String> {
        let requirements = problem.constraints();
        if requirements.is_empty() { return Err("--all-constraints requires a nonempty authored constraint family".into()); }
        let mut admitted = vec![None; requirements.len()];
        for (name, tolerance) in tolerances {
            let i = requirements.iter().position(|r| &r.name == name)
                .ok_or_else(|| format!("unknown equality constraint {name}"))?;
            if requirements[i].sense != ConstraintSense::Equal || admitted[i].is_some()
                || !tolerance.is_finite() || *tolerance < 0.0 {
                return Err("tolerances apply once to named equalities, in finite nonnegative physical units".into());
            }
            admitted[i] = Some(*tolerance);
        }
        for (requirement, tolerance) in requirements.iter().zip(&admitted) {
            if requirement.sense == ConstraintSense::Equal && tolerance.is_none() {
                return Err(format!("declare --equality-tolerance {} in {} before sampling", requirement.name, unit(&requirement.quantity)));
            }
        }
        Ok(Self { requirements: requirements.to_vec(), tolerances: admitted,
            failures: vec![0; requirements.len()], ranges: vec![[f64::INFINITY, f64::NEG_INFINITY]; requirements.len()], samples: 0 })
    }

    pub(super) fn evaluate(&mut self, problem: &EquilibriumDesign, x: &[f64], work: &mut DesignControl,
        gate: &CancelGate) -> Result<f64, DesignError>
    {
        let result = problem.evaluate_forward_with_constraints(x, work, gate)?;
        self.record(&result.constraints)
    }

    fn record(&mut self, rows: &[ForwardConstraintResult]) -> Result<f64, DesignError> {
        if rows.len() != self.requirements.len() {
            return Err(DesignError::Invalid { what: "joint event needs every response row" });
        }
        // Complete admission before modifying any per-row statistics. Direct
        // physical comparisons cannot mistake an underflowed residual for zero.
        let mut failed = Vec::with_capacity(rows.len());
        for ((row, requirement), tolerance) in rows.iter().zip(&self.requirements).zip(&self.tolerances) {
            if !row.value.is_finite() || !row.residual.is_finite() {
                return Err(DesignError::Invalid { what: "nonfinite joint response" });
            }
            let accepted = match requirement.sense {
                ConstraintSense::AtMost => row.value <= requirement.bound,
                ConstraintSense::AtLeast => row.value >= requirement.bound,
                ConstraintSense::Equal => (row.value - requirement.bound).abs() <= tolerance.expect("admitted equality band"),
            };
            failed.push(!accepted);
        }
        self.samples = self.samples.checked_add(1).ok_or(DesignError::Invalid { what: "joint sample count overflow" })?;
        for (i, (row, failure)) in rows.iter().zip(&failed).enumerate() {
            self.failures[i] += usize::from(*failure);
            self.ranges[i][0] = self.ranges[i][0].min(row.value);
            self.ranges[i][1] = self.ranges[i][1].max(row.value);
        }
        // Feed ONE bounded observable to the original sampling/statistical
        // owner. Threshold zero means all constraints passed on this SAME draw.
        Ok(if failed.iter().any(|f| *f) { 1.0 } else { 0.0 })
    }
}

fn unit(quantity: &ResponseQuantity) -> &'static str {
    match quantity {
        ResponseQuantity::Displacement(_) | ResponseQuantity::ContactPenetration(_) => "m",
        ResponseQuantity::SpringForce(_) | ResponseQuantity::ContactForce(_) => "N",
    }
}
fn quantity(quantity: &ResponseQuantity) -> &'static str {
    match quantity {
        ResponseQuantity::Displacement(_) => "displacement",
        ResponseQuantity::SpringForce(_) => "signed-left-spring-force",
        ResponseQuantity::ContactForce(_) => "compressive-contact-force",
        ResponseQuantity::ContactPenetration(_) => "contact-penetration",
    }
}
fn sense(sense: ConstraintSense) -> &'static str {
    match sense { ConstraintSense::AtMost => "at-most", ConstraintSense::AtLeast => "at-least", ConstraintSense::Equal => "equal-with-explicit-band" }
}

pub(super) fn output(loaded: &EquilibriumDesignFile, options: &Options, result: &Outcome, event: &Event) -> String {
    let mut out = format!("{{\"schema\":\"frankensim-equilibrium-reliability-v1\",\"method\":\"{}\",\"status\":\"{}\",\"evidence\":\"Estimated\",\"no_claim\":{},\"model_blake3\":\"{}\",\"design_blake3\":\"{}\",\"seed\":\"{}\",\"unit\":\"1\",\"compliance_event\":\"all authored response constraints hold on the same parameter draw\",\"compliance_scope\":\"joint-all-authored-response-constraints\",\"unassessed_response_constraints\":0,\"physics_evaluation\":\"primal-only\",\"samples\":{},\"case_solves\":{},\"compliance_probability\":{:.17e},\"probability_of_any_failure\":{:.17e},\"compliance_standard_error\":{},\"completed_replicates\":{},\"standard_error_basis\":{},\"per_constraint_statistics\":\"descriptive same-draw marginals and empirical ranges; not simultaneous confidence bounds\",\"constraints\":[",
        if options.method == PropagationMethod::MonteCarlo { "mc" } else { "rqmc" }, result.status.label(),
        json_string(confidence_scope(options.policy)),
        loaded.model_info().input_hash.to_hex(), loaded.design_hash().to_hex(), options.seed,
        result.samples, result.work.case_solves, result.compliance_probability, result.mean_m,
        number(result.compliance_standard_error), result.completed_replicates.map_or_else(|| "null".into(), |n| n.to_string()),
        json_string(if options.method == PropagationMethod::MonteCarlo { "individual-monte-carlo-indicators; descriptive at optional stops" }
            else { "complete-independent-scramble-means" }));
    for (i, requirement) in event.requirements.iter().enumerate() {
        if i != 0 { out.push(','); }
        write!(&mut out, "{{\"name\":{},\"case\":{},\"quantity\":{},\"unit\":{},\"sense\":{},\"bound\":{:.17e},\"scale\":{:.17e},\"equality_tolerance\":{},\"failure_count\":{},\"empirical_compliance_probability\":{:.17e},\"observed_min\":{:.17e},\"observed_max\":{:.17e}}}",
            json_string(&requirement.name), json_string(&loaded.problem().load_cases()[requirement.case].name),
            json_string(quantity(&requirement.quantity)), json_string(unit(&requirement.quantity)),
            json_string(sense(requirement.sense)), requirement.bound, requirement.scale,
            number(event.tolerances[i]), event.failures[i], (event.samples-event.failures[i]) as f64 / event.samples as f64,
            event.ranges[i][0], event.ranges[i][1]).expect("String write");
    }
    out.push_str("],\"uncertainty_coordinates\":\"dimensionless x; physical p = reference + scale*x\",\"dependence\":\"independent input variables; responses need not be independent\",\"variables\":[");
    write_variables(&mut out, loaded, options);
    out.push(']');
    write_assessment(&mut out, options, result);
    out.push('}');
    out
}

#[cfg(test)]
#[path = "joint/tests.rs"]
mod tests;
