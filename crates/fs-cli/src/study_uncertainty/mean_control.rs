//! Frozen native calibration: support secants or one nominal adjoint. Neither
//! calibration solves nor their derivative reports are probability observations.

use fs_blake3::ContentHash;
use fs_exec::CancelGate;
use fs_ledger::Ledger;
use fs_package::Claim;
use fs_project::uncertainty::mean_control::{self as design, MeanControlMethod, MeanControlPolicy};
use fs_uq::{LinearControlEstimate, LinearControlVariate, UqStatus};
use fs_uq::product_copula::mean_control::CopulaLinearControlVariate;
use fs_uq::product_qmc::mean_control::{QmcLinearControlEstimate, QmcLinearControlVariate};

use super::{Execution, Failure, J, Loaded, Model, Result, Sample, fail, hash_field,
    integer, linked, optional, parse, quoted, read_rows_value, rows_json};

pub(super) const KIND: &str = "native-uncertainty-mean-control";
const SCHEMA: &str = "frankensim.cli.native-mean-control.v1";
const SCOPE: &str = "Fixed coefficients from predeclared whole-model support secants, not an adjoint or a local derivative certificate. Calibration solves are not probability observations. The controlled mean estimates the same model expectation; raw responses, quantiles and compliance remain unchanged. Descriptive standard errors do not justify optional stopping or bound solver, finite-grid, transform or physical-model error. Variance reduction is measured, never guaranteed.";
const ADJOINT_SCOPE: &str = "Coefficients from one final-field native nominal adjoint at the declared physical marginal means, fixed before sampling. A selected hottest-node functional is used even at a tie; no unique maximum derivative or derivative-error bound is claimed. Constant inputs need no coefficient solve. Calibration is not a probability observation. The centered control estimates the same model mean; raw temperatures, quantiles and compliance are unchanged. Descriptive standard errors do not justify optional stopping or bound solver, finite-grid, transform, continuum or physical-model error. Variance reduction is measured, never guaranteed.";
fn scope(method: MeanControlMethod) -> &'static str {
    match method { MeanControlMethod::CoordinateSecant => SCOPE, MeanControlMethod::NominalAdjoint => ADJOINT_SCOPE }
}
fn policy(model: &Model) -> MeanControlPolicy {
    model.bound.study().mean_control().expect("calibration policy")
}

enum Frozen {
    MonteCarlo(LinearControlVariate),
    Copula(CopulaLinearControlVariate),
    QuasiMonteCarlo(QmcLinearControlVariate),
    CopulaQuasiMonteCarlo(QmcLinearControlVariate),
}

pub(super) struct Calibration {
    pub(super) probes: Vec<Sample>,
    pub(super) attempted: usize,
    pub(super) failure: Option<String>,
    control: Option<Frozen>,
    method: MeanControlMethod,
}
fn invalid(message: impl Into<String>) -> Failure { fail("cli-uncertainty-mean-control", message) }
fn same_bits(a: &[f64], b: &[f64]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(x, y)| x.to_bits() == y.to_bits())
}
fn numbers(values: &[f64]) -> String {
    format!("[{}]", values.iter().map(ToString::to_string).collect::<Vec<_>>().join(","))
}
fn read_numbers(value: &J) -> Result<Vec<f64>> {
    value.as_array().filter(|v| v.len() <= 32)
        .ok_or_else(|| invalid("invalid coefficient array"))?
        .iter().map(|v| v.as_f64().filter(|x| x.is_finite())
            .ok_or_else(|| invalid("nonfinite coefficient"))).collect()
}

impl Calibration {
    /// Reconstruct from sealed calibration evidence only. Never inspect random
    /// observations to choose coefficients, and never repeat completed physics.
    pub(super) fn recover(model: &Model, ledger: &Ledger, prior: Option<&Loaded>) -> Result<Option<Self>> {
        let Some(policy) = model.bound.study().mean_control() else {
            if prior.is_some_and(|p| p.value.get("mean_control_state").is_some()) {
                return Err(invalid("uncontrolled study carries a calibration state"));
            }
            return Ok(None);
        };
        let mut state = Self { probes: Vec::new(), attempted: 0, failure: None,
            control: None, method: policy.method };
        let saved = if let Some(prior) = prior {
            let bytes = linked(ledger, &prior.value, "mean_control_state", KIND)?;
            let value = parse(&bytes)?;
            if value.str_field("schema") != Some(SCHEMA)
                || hash_field(&value, "model")? != model.identity()
                || value.get("failure") != Some(&J::Null)
            { return Err(invalid("calibration version/model differs or failed calibration was resumed")); }
            state.probes = read_rows_value(value.get("probes")
                .ok_or_else(|| invalid("missing calibration probes"))?)?;
            state.attempted = integer(&value, "attempted")?;
            if state.probes.len() > policy.probe_count(model.bound.study().parameters())
                || state.attempted != state.probes.len()
            { return Err(invalid("calibration counts differ from the original probe plan")); }
            for (ordinal, sample) in state.probes.iter().enumerate() {
                let point = policy.probe(model.bound.study().parameters(), ordinal)
                    .map_err(|e| invalid(e.detail))?;
                if !same_bits(&point, &sample.parameters) {
                    return Err(invalid("retained calibration point differs from its declared ordinal"));
                }
                if state.method == MeanControlMethod::CoordinateSecant {
                    model.verify_sample(ledger, sample)?;
                }
                // Nominal child's report-only project and all sealed lineage
                // are verified by nominal_coefficients in finish below.
            }
            Some(value)
        } else { None };
        state.finish(model, ledger)?;
        if let Some(saved) = saved {
            let value = saved.get("gradient").ok_or_else(|| invalid("missing frozen gradient field"))?;
            match state.gradient() {
                Some(gradient) if same_bits(gradient, &read_numbers(value)?) => {},
                None if value == &J::Null => {},
                _ => return Err(invalid("frozen coefficients differ from the sealed calibration solves")),
            }
            if state.failure.is_some() || (integer(&prior.expect("saved state").value, "samples_completed")? > 0 && !state.ready()) {
                return Err(invalid("random observations precede complete successful calibration"));
            }
        }
        Ok(Some(state))
    }

    fn finish(&mut self, model: &Model, ledger: &Ledger) -> Result<()> {
        if self.failure.is_some() || self.control.is_some() || self.pending(model) { return Ok(()); }
        let calculated = match self.method {
            MeanControlMethod::CoordinateSecant => {
                let values: Vec<_> = self.probes.iter().map(|p| p.value_k).collect();
                design::coefficients(model.bound.study().parameters(), &values).map_err(|e| invalid(e.detail))
            }
            MeanControlMethod::NominalAdjoint if self.probes.is_empty() =>
                Ok(vec![0.0; model.bound.study().parameters().len()]),
            MeanControlMethod::NominalAdjoint => model.nominal_coefficients(ledger, &self.probes[0]),
        };
        let gradient = match calculated {
            Ok(gradient) => gradient,
            Err(e) => { self.failure = Some(format!("calibration coefficient refused: {e}")); return Ok(()); }
        };
        let fresh = Execution::new(model)?;
        self.control = Some(match &fresh {
            Execution::MonteCarlo(run) => Frozen::MonteCarlo(run.freeze_linear_control_variate(&gradient).map_err(|e| invalid(e.to_string()))?),
            Execution::CopulaMonteCarlo(run) => Frozen::Copula(run.freeze_linear_control_variate(&gradient).map_err(|e| invalid(e.to_string()))?),
            Execution::QuasiMonteCarlo(run) => Frozen::QuasiMonteCarlo(run.freeze_linear_control_variate(&gradient).map_err(|e| invalid(e.to_string()))?),
            Execution::CopulaQuasiMonteCarlo(run) => Frozen::CopulaQuasiMonteCarlo(run.freeze_linear_control_variate(&gradient).map_err(|e| invalid(e.to_string()))?),
        });
        Ok(())
    }

    pub(super) fn pending(&self, model: &Model) -> bool {
        self.probes.len() < policy(model).probe_count(model.bound.study().parameters())
    }
    pub(super) fn ready(&self) -> bool { self.control.is_some() }
    fn gradient(&self) -> Option<&[f64]> {
        self.control.as_ref().map(|control| match control {
            Frozen::MonteCarlo(c) => c.gradient(), Frozen::Copula(c) => c.gradient(),
            Frozen::QuasiMonteCarlo(c) | Frozen::CopulaQuasiMonteCarlo(c) => c.gradient(),
        })
    }
    pub(super) fn next(&self, model: &Model) -> Result<Vec<f64>> {
        if self.failure.is_some() || self.ready() { return Err(invalid("calibration is not pending")); }
        if self.attempted >= policy(model).max_solves {
            return Err(invalid("original calibration solve allowance exhausted"));
        }
        policy(model).probe(model.bound.study().parameters(), self.probes.len()).map_err(|e| invalid(e.detail))
    }
    pub(super) fn sample(&self, model: &Model, ledger: &Ledger, gate: &CancelGate,
        point: &[f64], remaining_wall_s: f64) -> Result<Option<Sample>> {
        if !same_bits(&self.next(model)?, point) { return Err(invalid("calibration point differs")); }
        match self.method {
            MeanControlMethod::CoordinateSecant => model.sample(ledger, gate, point, remaining_wall_s),
            MeanControlMethod::NominalAdjoint => model.sample_nominal_adjoint(ledger, gate, point, remaining_wall_s),
        }
    }
    pub(super) fn accept(&mut self, model: &Model, ledger: &Ledger, sample: Sample) -> Result<()> {
        if !same_bits(&self.next(model)?, &sample.parameters) || !sample.value_k.is_finite() {
            return Err(invalid("completed calibration child differs from the proposed point"));
        }
        self.attempted += 1;
        self.probes.push(sample);
        self.finish(model, ledger)
    }
    pub(super) fn reject(&mut self, error: &Failure) {
        self.attempted += 1;
        self.failure = Some(error.to_string());
    }
    pub(super) fn snapshot(&self, model: &Model) -> String {
        // The retained canonical study already binds the exact method. Keep the
        // secant envelope unchanged; no caller-supplied coefficient is trusted.
        format!("{{\"schema\":{SCHEMA:?},\"model\":{},\"attempted\":{},\"probes\":{},\"failure\":{},\"gradient\":{}}}",
            quoted(&model.identity().to_hex()), self.attempted, rows_json(&self.probes),
            self.failure.as_deref().map_or_else(|| "null".into(), quoted),
            self.gradient().map_or_else(|| "null".into(), numbers))
    }
    pub(super) fn assessment(&self, execution: &Execution) -> Result<Option<Assessment>> {
        if self.failure.is_some() || execution.status() == UqStatus::Refused { return Ok(None); }
        let Some(control) = &self.control else { return Ok(None); };
        let estimate = match (execution, control) {
            (Execution::MonteCarlo(run), Frozen::MonteCarlo(c)) => run.assess_linear_control_variate(c)
                .map(|a| a.map(Estimate::MonteCarlo)).map_err(|e| invalid(e.to_string())),
            (Execution::CopulaMonteCarlo(run), Frozen::Copula(c)) => run.assess_linear_control_variate(c)
                .map(|a| a.map(Estimate::MonteCarlo)).map_err(|e| invalid(e.to_string())),
            (Execution::QuasiMonteCarlo(run), Frozen::QuasiMonteCarlo(c)) => run.assess_linear_control_variate(c)
                .map(|a| Some(Estimate::QuasiMonteCarlo(a))).map_err(|e| invalid(e.to_string())),
            (Execution::CopulaQuasiMonteCarlo(run), Frozen::CopulaQuasiMonteCarlo(c)) => run.assess_linear_control_variate(c)
                .map(|a| Some(Estimate::QuasiMonteCarlo(a))).map_err(|e| invalid(e.to_string())),
            _ => Err(invalid("frozen control belongs to a different sampling method")),
        }?;
        Ok(estimate.map(|estimate| Assessment { estimate, method: self.method }))
    }
    pub(super) fn report(&self, model: &Model, assessment: Option<&Assessment>) -> String {
        let p = model.bound.study().parameters();
        let order = p.iter().map(|p| quoted(&p.name)).collect::<Vec<_>>().join(",");
        let units = p.iter().map(|p| quoted(p.target.unit())).collect::<Vec<_>>().join(",");
        format!("{{\"method\":{},\"status\":{},\"probe_solves_planned\":{},\"probe_solves_allowed\":{},\"probe_solves_completed\":{},\"probe_solves_attempted\":{},\"parameter_order\":[{order}],\"parameter_units\":[{units}],\"parameter_means\":{},\"gradient\":{},\"calibration\":{},\"estimate\":{},\"scope\":{}}}",
            quoted(self.method.tag()),
            quoted(if self.failure.is_some() { "refused" } else if self.ready() { "frozen" } else { "calibrating" }),
            policy(model).probe_count(p), policy(model).max_solves,
            self.probes.len(), self.attempted, numbers(&design::means(p)),
            self.gradient().map_or_else(|| "null".into(), numbers), self.snapshot(model),
            assessment.map_or_else(|| "null".into(), Assessment::json), quoted(scope(self.method)))
    }
}

enum Estimate { MonteCarlo(LinearControlEstimate), QuasiMonteCarlo(QmcLinearControlEstimate) }
pub(super) struct Assessment { estimate: Estimate, method: MeanControlMethod }
impl Assessment {
    fn summary(&self) -> (Option<f64>, Option<f64>, Option<f64>, Option<f64>, Option<f64>, usize, &'static str) {
        match &self.estimate {
            Estimate::MonteCarlo(a) => (Some(a.raw_mean), Some(a.mean), a.raw_standard_error,
                a.standard_error, a.variance_ratio, a.n, "independent-monte-carlo-vectors"),
            Estimate::QuasiMonteCarlo(a) => (a.raw.as_ref().map(|v| v.mean), a.controlled.as_ref().map(|v| v.mean),
                a.raw.as_ref().and_then(|v| v.standard_error), a.controlled.as_ref().and_then(|v| v.standard_error),
                a.variance_ratio, a.samples_in_estimate, "complete-independent-scrambles"),
        }
    }
    pub(super) fn json(&self) -> String {
        let (raw, mean, raw_error, error, ratio, n, unit) = self.summary();
        let nets = match &self.estimate {
            Estimate::MonteCarlo(_) => String::new(),
            Estimate::QuasiMonteCarlo(a) => format!(",\"completed_replicates\":{},\"samples_accepted\":{},\"raw_replicate_means_k\":{},\"controlled_replicate_means_k\":{}",
                a.completed_replicates, a.samples_accepted, numbers(&a.raw_replicate_means), numbers(&a.controlled_replicate_means)),
        };
        format!("{{\"raw_mean_k\":{},\"controlled_mean_k\":{},\"raw_standard_error_k\":{},\"controlled_standard_error_k\":{},\"variance_ratio\":{},\"samples_in_estimate\":{n},\"statistical_unit\":{unit:?}{nets}}}",
            optional(raw), optional(mean), optional(raw_error), optional(error), optional(ratio))
    }
    pub(super) fn html(&self) -> String {
        let (raw, mean, raw_error, error, ratio, n, unit) = self.summary();
        format!("<h2>Frozen mean control</h2><p>{n} random responses in the mean estimate; statistical units: {unit}. Raw mean {} K; controlled mean {} K. Raw standard error {} K; controlled standard error {} K; variance ratio {}.</p><p>{}</p>",
            optional(raw), optional(mean), optional(raw_error), optional(error), optional(ratio), scope(self.method))
    }
    pub(super) fn claim(&self, result: ContentHash) -> Option<Claim> {
        let (_, mean, _, error, _, n, unit) = self.summary();
        let mean = mean?; let error = error?;
        let method = match self.method {
            MeanControlMethod::CoordinateSecant => "predeclared-native-coordinate-secant-mean-control",
            MeanControlMethod::NominalAdjoint => "predeclared-native-nominal-adjoint-mean-control",
        };
        Some(Claim::estimated("cooling.uncertainty.controlled-mean",
            format!("{mean} K; descriptive standard error {error} K on {n} raw native responses using {unit}. Result {}. {}", result.to_hex(), scope(self.method)),
            method, error))
    }
}
