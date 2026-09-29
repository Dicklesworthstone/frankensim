//! Nominal calibration reuses the native import/solve/recovery path. Its sole
//! project change requests an adjoint report; random child projects stay intact.

use super::*;
use fs_project::uncertainty::{Target, UniformParameter};

const OUTPUT: &str = "temperature-max-adjoint";

impl Model {
    /// Construct an admitted report-only view of the SAME physical model.
    /// Cloned assets remain under the existing input envelope (memory/4 each);
    /// no second file load, source substitution or independent physical solver.
    fn nominal_view(&self) -> Result<Self> {
        let mut spec = self.base.spec.clone();
        let outputs = spec.outputs.get_or_insert_with(Vec::new);
        let requests: Vec<_> = outputs.iter().filter(|row| row.name == OUTPUT).collect();
        if requests.len() > 1 || requests.first().is_some_and(|row| row.kind != "report") {
            return Err(invalid("nominal calibration requires one unambiguous adjoint report output"));
        }
        if requests.is_empty() {
            outputs.push(fs_project::spec::OutputRequest { name: OUTPUT.into(), kind: "report".into() });
        }
        let text = fs_project::print_sexpr(&spec).map_err(project_error)?;
        let base = fs_project::parse_sexpr(&text).map_err(project_error)?;
        let bound = self.bound.study().clone().bind(&base.spec).map_err(project_error)?;
        Self::from_parts(bound, base, self.geometry.clone(), self.cards.clone())
    }

    pub(crate) fn sample_nominal_adjoint(
        &self, ledger: &Ledger, gate: &CancelGate, point: &[f64], remaining_wall_s: f64,
    ) -> Result<Option<Sample>> {
        if gate.is_requested() || remaining_wall_s == 0.0 { return Ok(None); }
        self.nominal_view()?.sample(ledger, gate, point, remaining_wall_s)
    }

    /// Recover from sealed native evidence without rerunning calibration. The
    /// requested-report child has a distinct project/run identity from raw
    /// probability children, even when their physical coordinates coincide.
    pub(crate) fn nominal_coefficients(&self, ledger: &Ledger, sample: &Sample) -> Result<Vec<f64>> {
        let view = self.nominal_view()?;
        view.verify_sample(ledger, sample)?;
        let sealed = crate::solve::load_completed_run(ledger, &sample.run)
            .map_err(|e| invalid(format!("{}: {}", e.code, e.what)))?;
        let hash = sealed.stages.iter().find(|row| row.0 == "conduction")
            .and_then(|row| ContentHash::from_hex(&row.2))
            .ok_or_else(|| invalid("nominal calibration has no sealed conduction receipt"))?;
        let bytes = artifact(ledger, hash, "solve-stage-receipt", 4 * 1024 * 1024)?;
        let receipt = J::parse(utf8(&bytes)?).map_err(|e| invalid(e.to_string()))?;
        if receipt.str_field("run") != Some(sample.run.as_str()) {
            return Err(invalid("nominal adjoint belongs to a different physical run"));
        }
        let report = receipt.get("nominal_adjoint")
            .ok_or_else(|| invalid("nominal calibration did not produce its requested adjoint"))?;
        let region = &self.bound.base().requirements.as_ref().expect("bound requirement")[0].region;
        coefficients(report, self.bound.study().parameters(), sample.value_k, region)
    }
}

fn target_name(target: Target) -> &'static str {
    match target {
        Target::Power => "power",
        Target::ConvectionCoefficient => "convection-coefficient",
        Target::ConvectionTemperature => "convection-temperature",
        Target::HeatFlux => "heat-flux",
        Target::AirInletTemperature => "air-inlet-temperature",
        Target::FanSpeedRatio => "fan-speed-ratio",
    }
}

fn coefficients(report: &J, parameters: &[UniformParameter], nominal: f64, region: &str) -> Result<Vec<f64>> {
    if report.str_field("schema") != Some("frankensim.cli.nominal-adjoint.v1")
        || report.str_field("functional") != Some("selected-nodal-temperature")
        || report.str_field("region") != Some(region)
        || report.f64_field("value_k").map(f64::to_bits) != Some(nominal.to_bits())
        || report.f64_field("true_relative_residual").is_none_or(|v| !v.is_finite() || v < 0.0)
    { return Err(invalid("nominal adjoint version, functional, region, value or residual differs")); }
    let rows = array(report, "parameters", 256)?;
    parameters.iter().map(|p| {
        // A constant centered input is identically zero, even when its physical
        // derivative is unsupported. No derivative value is invented for it.
        if p.low == p.high { return Ok(0.0); }
        let matches: Vec<_> = rows.iter().filter(|row| row.str_field("target") == Some(target_name(p.target))
            && row.str_field("entity") == Some(p.entity.as_str())).collect();
        if matches.len() != 1 || matches[0].str_field("parameter_unit") != Some(p.target.unit()) {
            return Err(invalid(format!("nominal adjoint has missing, ambiguous or wrong-unit derivative for {}", p.name)));
        }
        matches[0].f64_field("derivative").filter(|v| v.is_finite())
            .ok_or_else(|| invalid(format!("nonfinite nominal coefficient for {}", p.name)))
    }).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    fn parameters() -> Vec<UniformParameter> {
        vec![UniformParameter { name:"power".into(),entity:"solid".into(),target:Target::Power,low:2.0,high:6.0 }]
    }
    fn report(rows: &str) -> J {
        J::parse(&format!(r#"{{"schema":"frankensim.cli.nominal-adjoint.v1","functional":"selected-nodal-temperature","region":"solid","value_k":301,"true_relative_residual":0.000000001,"parameters":[{rows}]}}"#)).unwrap()
    }
    const ROW: &str = r#"{"target":"power","entity":"solid","parameter_unit":"W","derivative":0.7}"#;
    #[test]
    fn nominal_adjoint_coefficients_bind_physical_target_and_units() {
        assert_eq!(coefficients(&report(ROW), &parameters(), 301.0, "solid").unwrap(),[0.7]);
        for rows in ["".to_string(), format!("{ROW},{ROW}"), ROW.replace("\"W\"","\"K\""), ROW.replace("power","heat-flux")] {
            assert!(coefficients(&report(&rows), &parameters(), 301.0, "solid").is_err());
        }
        assert!(coefficients(&report(ROW), &parameters(), 302.0, "solid").is_err());
        assert!(coefficients(&report(ROW), &parameters(), 301.0, "other").is_err());
    }
    #[test]
    fn nominal_adjoint_singletons_have_no_fabricated_physical_derivative() {
        let mut p=parameters(); p[0].high=p[0].low; p[0].target=Target::FanSpeedRatio;
        assert_eq!(coefficients(&report(""), &p, 301.0, "solid").unwrap(),[0.0]);
        p[0].high=3.0;
        assert!(coefficients(&report(""), &p, 301.0, "solid").is_err());
    }
}
