//! Priority-ordered coordinate brackets over COMPLETE evaluated trajectories.
//! This is a declared allocation policy, not a general or global optimizer.
use super::*;

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
}

struct Runner<'a, F> {
    plan: &'a Plan,
    deadline: Instant,
    run: F,
    attempted: usize,
    steps: usize,
    solves: usize,
    history: Vec<J>,
    newton_trials: usize,
}
impl<F: FnMut(&[f64]) -> Result<Evaluation>> Runner<'_, F> {
    fn poll(&self) -> Result<()> {
        if Instant::now() >= self.deadline { Err(budget("component allocation wall deadline exhausted")) }
        else { Ok(()) }
    }
    fn evaluate(&mut self, values: &[f64], proposal: &str) -> Result<Evaluation> {
        self.poll()?;
        if self.attempted >= self.plan.max_evaluations {
            return Err(budget("component allocation evaluation budget exhausted"));
        }
        let after = self.steps.checked_add(self.plan.planned_steps)
            .ok_or_else(|| budget("allocation cumulative step count overflow"))?;
        if after > self.plan.max_total_steps {
            return Err(budget("component allocation cumulative trajectory-step budget exhausted"));
        }
        self.attempted += 1;
        let result = (self.run)(values)?;
        if result.steps != self.plan.planned_steps || result.slopes.len() != values.len()
            || !result.peak.is_finite() || !result.peak_time.is_finite() {
            return Err(model_failure("allocation evaluator returned an incompatible complete trajectory"));
        }
        self.steps = after;
        self.solves = self.solves.checked_add(result.solves).ok_or_else(|| budget("allocation solve count overflow"))?;
        if proposal == "newton" { self.newton_trials += 1; }
        self.history.push(J::Object(vec![
            ("power_w".into(), number_array(values)?),
            ("sampled_peak_k".into(), number_value(result.peak)?),
            ("sampled_peak_time_s".into(), number_value(result.peak_time)?),
            ("passing".into(), J::Bool(result.peak <= self.plan.limit)),
            ("proposal".into(), J::Str(proposal.into())),
            ("accepted_steps".into(), usize_value(result.steps)),
            ("solid_solves".into(), usize_value(result.solves)),
        ]));
        Ok(result)
    }
}

fn proposed(low: f64, high: f64, at: f64, temperature: f64, slope: Option<f64>, limit: f64) -> Option<f64> {
    let slope = slope.filter(|s| s.is_finite() && *s > 0.0)?;
    let next = at - (temperature-limit)/slope;
    let guard = 0.1*(high-low);
    (next.is_finite() && next > low+guard && next < high-guard).then_some(next)
}

pub(super) fn allocate(plan: &Plan, deadline: Instant,
    run: impl FnMut(&[f64]) -> Result<Evaluation>) -> Result<Outcome> {
    let mut runner = Runner { plan, deadline, run, attempted: 0, steps: 0, solves: 0,
        history: Vec::new(), newton_trials: 0 };
    let mut values: Vec<_> = plan.axes.iter().map(|a| a.minimum).collect();
    let mut passing = runner.evaluate(&values, "minimum-baseline")?;
    if passing.peak > plan.limit {
        return Err(model_failure("declared minimum allocation exceeds the sampled-peak limit; no global infeasibility conclusion"));
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
        if upper.peak <= plan.limit {
            values = trial_values; passing = upper;
        } else {
            // Retain only the failing scalar row, not its full fields. The
            // passing result already describes exactly the current vector.
            let mut high = axis.maximum;
            let mut hot = upper.peak;
            let mut hot_slope = upper.slopes[index];
            drop(upper);
            loop {
                if let Err(error) = runner.poll() { reason = Some(error.message); break 'priorities; }
                let low = values[index];
                if high-low <= plan.power_tolerance && plan.limit-passing.peak <= plan.temperature_tolerance {
                    failed = Some((high,hot)); break;
                }
                let cool = (low,passing.peak,passing.slopes[index]);
                let hot_row = (high,hot,hot_slope);
                let endpoints = if (cool.1-plan.limit).abs() < (hot-plan.limit).abs() {
                    [cool,hot_row]
                } else { [hot_row,cool] };
                let proposal = endpoints.into_iter().find_map(|(at,t,s)| proposed(low,high,at,t,s,plan.limit));
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
                // Never classify a model or adjoint failure as a hot design.
                if evaluated.peak <= plan.limit { values = trial_values; passing = evaluated; }
                else { high = next; hot = evaluated.peak; hot_slope = evaluated.slopes[index]; }
            }
        }
        let (failed_power,failed_peak,width) = match failed {
            Some((power,peak)) => (number_value(power)?,number_value(peak)?,power-values[index]),
            None => (J::Null,J::Null,0.0),
        };
        decisions.push(J::Object(vec![
            ("priority_index".into(), usize_value(index)),
            ("selected_power_w".into(), number_value(values[index])?),
            ("status".into(), J::Str(if failed.is_some() { "target-bracketed" } else { "upper-bound-passing" }.into())),
            ("failed_upper_power_w".into(), failed_power),
            ("failed_upper_peak_k".into(), failed_peak),
            ("power_bracket_width_w".into(), number_value(width)?),
            ("conditioning".into(), J::Str("earlier priorities fixed at selections; later priorities at minima; not a final-vector optimality certificate".into())),
        ]));
        completed += 1;
    }
    if reason.is_none() {
        if let Err(error) = runner.poll() { reason = Some(error.message); }
    }
    Ok(Outcome { values, passing, baseline_peak, completed, reason, decisions,
        history: runner.history, attempted: runner.attempted, steps: runner.steps,
        solves: runner.solves, newton_trials: runner.newton_trials })
}
