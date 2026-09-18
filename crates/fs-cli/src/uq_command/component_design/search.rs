//! Priority-ordered coordinate brackets over COMPLETE evaluated candidate sets.
//! Each temperature objective keeps its own limit; feasibility is their
//! intersection, and slopes belong to the active signed excess, not the hottest
//! absolute temperature. This remains a local declared-order allocation policy.
use super::*;
use constraint_execution::Work;

pub(super) struct Outcome {
    pub values: Vec<f64>,
    pub passing: Evaluation,
    pub baseline_peak: f64,
    pub completed: usize,
    pub reason: Option<String>,
    pub decisions: Vec<J>,
    pub history: Vec<J>,
    pub attempted: usize,
    pub steps: usize,
    pub solves: usize,
    pub newton_trials: usize,
    pub work: Work,
}

struct Runner<'a, F> {
    plan: &'a Plan,
    deadline: Instant,
    run: F,
    criteria: usize,
    attempted: usize,
    work: Work,
    history: Vec<J>,
    newton_trials: usize,
}
impl<F: FnMut(&[f64], &mut Work) -> Result<Evaluation>> Runner<'_, F> {
    fn poll(&self) -> Result<()> {
        if Instant::now() >= self.deadline { Err(budget("component allocation wall deadline exhausted")) }
        else { Ok(()) }
    }
    fn evaluate(&mut self, values: &[f64], proposal: &str) -> Result<Evaluation> {
        self.poll()?;
        if self.attempted >= self.plan.max_evaluations {
            return Err(budget("component allocation evaluation budget exhausted"));
        }
        let required = self.plan.planned_steps.checked_mul(self.criteria)
            .ok_or_else(|| budget("multi-limit candidate step count overflow"))?;
        let after = self.work.steps.checked_add(required)
            .ok_or_else(|| budget("allocation cumulative step count overflow"))?;
        if after > self.plan.max_total_steps {
            return Err(budget("component allocation cumulative trajectory-step budget cannot admit all thermal constraints"));
        }
        self.attempted += 1;
        let before = self.work;
        let result = (self.run)(values, &mut self.work)?;
        // Finished observations of an interrupted set remain charged, but a
        // successful set must contain exactly the declared number of solves.
        if result.steps != required || result.slopes.len() != values.len()
            || !result.peak.is_finite() || !result.peak_time.is_finite()
            || self.work.completed.checked_sub(before.completed) != Some(self.criteria)
            || self.work.steps.checked_sub(before.steps) != Some(required)
            || self.work.solves.checked_sub(before.solves) != Some(result.solves)
            || result.constraints.is_some() != (self.criteria > 1) {
            return Err(model_failure("allocation evaluator returned an incompatible complete thermal-constraint set"));
        }
        let margin = result.margin(self.plan)?;
        if proposal == "newton" { self.newton_trials += 1; }
        // Keep trial history O(controls), not O(observed vertices * criteria).
        // Full selectors and per-criterion derivatives belong to the selected
        // result only; repeating them here can multiply a large input by 4096.
        self.history.push(J::Object(vec![
            ("power_w".into(), number_array(values)?),
            ("sampled_peak_k".into(), number_value(result.peak)?),
            ("sampled_peak_time_s".into(), number_value(result.peak_time)?),
            ("maximum_temperature_excess_k".into(), number_value(margin)?),
            ("active_constraint".into(), J::Str(result.active_constraint().into())),
            ("passing".into(), J::Bool(margin <= 0.0)),
            ("proposal".into(), J::Str(proposal.into())),
            ("accepted_steps".into(), usize_value(result.steps)),
            ("solid_solves".into(), usize_value(result.solves)),
        ]));
        Ok(result)
    }
}

fn proposed(low: f64, high: f64, at: f64, excess: f64, slope: Option<f64>) -> Option<f64> {
    let slope = slope.filter(|s| s.is_finite() && *s > 0.0)?;
    let next = at - excess/slope;
    let guard = 0.1*(high-low);
    (next.is_finite() && next > low+guard && next < high-guard).then_some(next)
}

// Keep the existing scalar numerical test seam. Production goes through the
// metered path below, whose work survives a later constraint's interruption.
#[cfg(test)]
pub(super) fn allocate(plan: &Plan, deadline: Instant,
    mut run: impl FnMut(&[f64]) -> Result<Evaluation>) -> Result<Outcome> {
    allocate_with_work(plan, deadline, 1, |values, work| {
        work.begin()?;
        let result = run(values)?;
        work.record(result.steps, result.solves)?;
        Ok(result)
    })
}

pub(super) fn allocate_with_work(plan: &Plan, deadline: Instant, criteria: usize,
    run: impl FnMut(&[f64], &mut Work) -> Result<Evaluation>) -> Result<Outcome> {
    if !(1..=thermal_constraints::MAX_CONSTRAINTS).contains(&criteria) {
        return Err(bad("allocation requires 1..16 complete thermal constraints"));
    }
    let mut runner = Runner { plan, deadline, run, criteria, attempted: 0, work: Work::default(),
        history: Vec::new(), newton_trials: 0 };
    let mut values: Vec<_> = plan.axes.iter().map(|a| a.minimum).collect();
    let mut passing = runner.evaluate(&values, "minimum-baseline")?;
    if passing.margin(plan)? > 0.0 {
        return Err(model_failure(format!("declared minimum allocation exceeds thermal constraint {} by {} K; no global infeasibility conclusion",
            passing.active_constraint(), passing.margin(plan)?)));
    }
    let baseline_peak = passing.peak;
    let mut completed = 0;
    let mut decisions = Vec::new();
    let mut reason = None;
    'priorities: for (index, axis) in plan.axes.iter().enumerate() {
        let mut trial_values = values.clone();
        trial_values[index] = axis.maximum;
        let upper = match runner.evaluate(&trial_values, "upper-bound") {
            Ok(result) => result,
            Err(error) if error.code == "cooling-network-uq-budget" => { reason = Some(error.message); break; }
            Err(error) => return Err(error),
        };
        let mut failed = None;
        if upper.margin(plan)? <= 0.0 {
            values = trial_values; passing = upper;
        } else {
            let mut high = axis.maximum;
            let mut hot_peak = upper.peak;
            let mut hot = upper.margin(plan)?;
            let mut hot_constraint = upper.active_constraint().to_string();
            let mut hot_slope = upper.slopes[index];
            drop(upper);
            loop {
                if let Err(error) = runner.poll() { reason = Some(error.message); break 'priorities; }
                let low = values[index];
                let cool_margin = passing.margin(plan)?;
                if high-low <= plan.power_tolerance && -cool_margin <= plan.temperature_tolerance {
                    failed = Some((high,hot_peak,hot,hot_constraint)); break;
                }
                let cool = (low,cool_margin,passing.slopes[index]);
                let hot_row = (high,hot,hot_slope);
                let endpoints = if cool_margin.abs() < hot.abs() { [cool,hot_row] } else { [hot_row,cool] };
                let proposal = endpoints.into_iter().find_map(|(at,t,s)| proposed(low,high,at,t,s));
                let next = proposal.unwrap_or(0.5*low+0.5*high);
                if !(next > low && next < high) {
                    reason = Some("watt resolution exhausted before both requested tolerances were met".into());
                    break 'priorities;
                }
                let mut trial_values = values.clone(); trial_values[index] = next;
                let evaluated = match runner.evaluate(&trial_values, if proposal.is_some() { "newton" } else { "bisection" }) {
                    Ok(result) => result,
                    Err(error) if error.code == "cooling-network-uq-budget" => { reason = Some(error.message); break 'priorities; }
                    Err(error) => return Err(error),
                };
                let margin = evaluated.margin(plan)?;
                // Neither a missing constraint nor a failed requested adjoint
                // may be classified as a hot candidate or a passing subset.
                if margin <= 0.0 { values = trial_values; passing = evaluated; }
                else { high = next; hot_peak = evaluated.peak; hot = margin;
                    hot_constraint = evaluated.active_constraint().into(); hot_slope = evaluated.slopes[index]; }
            }
        }
        let (failed_power,failed_peak,failed_excess,failed_constraint,width) = match &failed {
            Some((power,peak,excess,name)) => (number_value(*power)?,number_value(*peak)?,number_value(*excess)?,J::Str(name.clone()),*power-values[index]),
            None => (J::Null,J::Null,J::Null,J::Null,0.0),
        };
        decisions.push(J::Object(vec![
            ("priority_index".into(), usize_value(index)),
            ("selected_power_w".into(), number_value(values[index])?),
            ("status".into(), J::Str(if failed.is_some() { "target-bracketed" } else { "upper-bound-passing" }.into())),
            ("failed_upper_power_w".into(), failed_power),
            ("failed_upper_peak_k".into(), failed_peak),
            ("failed_upper_temperature_excess_k".into(), failed_excess),
            ("failed_upper_constraint".into(), failed_constraint),
            ("power_bracket_width_w".into(), number_value(width)?),
            ("conditioning".into(), J::Str("every thermal constraint evaluated on the same allocation; earlier priorities fixed at selections; later priorities at minima; not a final-vector optimality certificate".into())),
        ]));
        completed += 1;
    }
    if reason.is_none() {
        if let Err(error) = runner.poll() { reason = Some(error.message); }
    }
    Ok(Outcome { values, passing, baseline_peak, completed, reason, decisions,
        history: runner.history, attempted: runner.attempted, steps: runner.work.steps,
        solves: runner.work.solves, newton_trials: runner.newton_trials, work: runner.work })
}
