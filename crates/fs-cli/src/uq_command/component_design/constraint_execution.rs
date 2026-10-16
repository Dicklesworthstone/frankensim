//! Complete candidate-set evaluation using the existing same-executable child.
//! Each temperature functional gets its OWN peak and adjoint. This deliberately
//! reruns the physical producer; it does not guess a component's trajectory peak
//! from the final field or reuse the primary objective's gradient.
use super::*;
use thermal_constraints::{Assessment, Constraint, Envelope};

#[derive(Debug, Default, Clone, Copy)]
pub(super) struct Work {
    pub attempted: usize,
    pub completed: usize,
    pub steps: usize,
    pub solves: usize,
}
impl Work {
    pub fn begin(&mut self) -> Result<()> {
        self.attempted = self.attempted.checked_add(1).ok_or_else(|| budget("trajectory attempt count overflow"))?;
        Ok(())
    }
    pub fn record(&mut self, steps: usize, solves: usize) -> Result<()> {
        let completed = self.completed.checked_add(1).ok_or_else(|| budget("trajectory count overflow"))?;
        let steps = self.steps.checked_add(steps).ok_or_else(|| budget("trajectory step count overflow"))?;
        let solves = self.solves.checked_add(solves).ok_or_else(|| budget("trajectory solve count overflow"))?;
        self.completed = completed; self.steps = steps; self.solves = solves;
        Ok(())
    }
}

pub(super) fn evaluate(plan: &Plan, constraints: &[Constraint], values: &[f64],
    deadline: Instant, work: &mut Work) -> Result<Evaluation> {
    let base = plan.request(values)?;
    let mut primary = evaluate_one(plan, values, &base, deadline, work)?;
    if constraints.is_empty() { return Ok(primary); }
    let mut rows = vec![Assessment { constraint: Constraint::primary(&base, plan.limit)?,
        peak: primary.peak, time: primary.peak_time, slopes: primary.slopes.clone() }];
    for constraint in constraints {
        let request = constraint.request(&base)?;
        let result = evaluate_one(plan, values, &request, deadline, work)?;
        same_forward(&primary.document, &result.document)?;
        primary.steps = primary.steps.checked_add(result.steps).ok_or_else(|| budget("candidate-set step count overflow"))?;
        primary.solves = primary.solves.checked_add(result.solves).ok_or_else(|| budget("candidate-set solve count overflow"))?;
        rows.push(Assessment { constraint: constraint.clone(), peak: result.peak,
            time: result.peak_time, slopes: result.slopes });
        // Drop this full secondary result before starting the next objective.
        // Its observation and derivative survive, not duplicated volume fields.
    }
    let envelope = Envelope::new(rows)?;
    primary.slopes = envelope.slopes();
    primary.constraints = Some(envelope);
    Ok(primary)
}

fn evaluate_one(plan: &Plan, values: &[f64], request: &J,
    deadline: Instant, work: &mut Work) -> Result<Evaluation> {
    let encoded = serialize(request)?;
    if encoded.len() as u64 > MAX_BASE_BYTES { return Err(bad("expanded allocation request exceeds the cooling input byte cap")); }
    if Instant::now() >= deadline { return Err(budget("allocation deadline exhausted before the next thermal constraint")); }
    work.begin()?;
    let document = child::evaluate_document(&encoded, deadline).map_err(|e| match e {
        child::EvaluationError::Budget => budget("allocation deadline interrupted a thermal constraint; incomplete candidate set not used"),
        child::EvaluationError::Child(message) => model_failure(message),
    })?;
    check_selector(request, &document)?;
    let result = inspect(plan, values, document)?;
    work.record(result.steps, result.solves)?;
    Ok(result)
}

/// A complete trajectory of the WRONG functional cannot satisfy this limit.
fn check_selector(request: &J, document: &J) -> Result<()> {
    let expected = field(request, "objective")?;
    let actual = field(document, "objective")?;
    let selected: Vec<_> = ["max_vertices", "max_solid_temperature", "max_wall_region", "mean_wall_region"]
        .into_iter().filter(|key| expected.get(key).is_some()).collect();
    if selected.len() != 1 || actual.str_field("kind") != Some(selected[0]) {
        return Err(model_failure("thermal constraint received a different objective selector"));
    }
    if matches!(selected[0], "max_wall_region" | "mean_wall_region")
        && actual.str_field("region") != expected.get(selected[0]).and_then(J::as_str) {
        return Err(model_failure("thermal constraint received a different observed surface"));
    }
    Ok(())
}

/// Changing the observation must not change the physical final field, time
/// window, or integrated energies. This is a replay guard, not a proof of
/// continuous-time accuracy. Objective-dependent reverse work can differ.
fn same_forward(first: &J, next: &J) -> Result<()> {
    let field_values = field(first, "solid_temperatures_k")?.as_array()
        .ok_or_else(|| model_failure("missing physical field in multi-limit result"))?;
    if field_values.is_empty() || first.get("solid_temperatures_k") != next.get("solid_temperatures_k") {
        return Err(model_failure("changing the thermal observation changed the physical final field"));
    }
    for key in ["branches", "node_temperatures_k", "source_w", "robin_out_w"] {
        if first.get(key) != next.get(key) {
            return Err(model_failure(format!("thermal constraint replay changed physical {key}")));
        }
    }
    let phase = if first.get("repeated_cycles").is_some() { "repeated_cycles" } else { "transient" };
    let first = field(first, phase)?;
    let next = field(next, phase)?;
    for key in ["input_energy_j", "stored_energy_change_j", "air_energy_gain_j",
        "radiative_energy_loss_j", "elapsed_time_s", "time_s"] {
        if first.get(key) != next.get(key) {
            return Err(model_failure(format!("thermal constraint replay changed physical {key}")));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn changing_the_observation_cannot_change_the_solved_physical_experiment() {
        let a = J::parse(r#"{"solid_temperatures_k":[301,302],"transient":{"input_energy_j":10,"air_energy_gain_j":2,"sampled_peak_objective_k":305}}"#).unwrap();
        let mut b = a.clone();
        input::put(input::member_mut(&mut b,"transient").unwrap(),"sampled_peak_objective_k",number_value(302.0).unwrap()).unwrap();
        assert!(same_forward(&a,&b).is_ok());
        input::put(input::member_mut(&mut b,"transient").unwrap(),"air_energy_gain_j",number_value(3.0).unwrap()).unwrap();
        assert!(same_forward(&a,&b).is_err());
        let mut b = a.clone();
        input::put(&mut b,"solid_temperatures_k",J::parse("[301,303]").unwrap()).unwrap();
        assert!(same_forward(&a,&b).is_err());
    }
    #[test]
    fn a_successful_child_must_observe_the_requested_surface() {
        let request = J::parse(r#"{"objective":{"max_wall_region":"memory-wall","gradient":false}}"#).unwrap();
        let wrong = J::parse(r#"{"objective":{"kind":"max_wall_region","region":"cpu-wall"}}"#).unwrap();
        assert!(check_selector(&request,&wrong).is_err());
        let right = J::parse(r#"{"objective":{"kind":"max_wall_region","region":"memory-wall"}}"#).unwrap();
        assert!(check_selector(&request,&right).is_ok());
    }
}
