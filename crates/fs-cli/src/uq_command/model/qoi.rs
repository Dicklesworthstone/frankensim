//! Select the actual cooling observable before sampling. A transient's final
//! objective is not its peak; repeated-cycle output retains only the last
//! cycle under `transient`, so its all-cycle peak has a separate owner.
use super::*;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Qoi {
    Steady,
    TransientPeak { repeated_cycles: Option<usize> },
}

impl Qoi {
    pub(crate) fn parse(value: Option<&J>, base: &J) -> Result<Self> {
        let kind = match value {
            Some(value) => {
                object(value, &["kind"], "qoi")?;
                field(value, "kind")?.as_str()
            }
            None => Some("steady-objective"),
        };
        match (kind, base.get("transient")) {
            (Some("steady-objective"), None) => Ok(Self::Steady),
            (Some("steady-objective"), Some(_)) => Err(bad(
                "a transient base requires qoi.kind=transient-sampled-peak; its final temperature is not the trajectory peak",
            )),
            (Some("transient-sampled-peak"), Some(schedule)) => {
                if schedule.as_object().is_none() { return Err(bad("transient must be an object")); }
                if schedule.get("fan_speed_design").is_some() || schedule.get("power_design").is_some() {
                    return Err(bad("transient UQ requires a fixed schedule, not a nested fan/workload design search"));
                }
                let repeated_cycles = match schedule.get("repeat") {
                    None => None,
                    Some(repeat) => {
                        if repeat.get("until_periodic").is_some() {
                            return Err(bad("transient UQ requires a fixed cycle count, not a sample-dependent periodic stopping horizon"));
                        }
                        let cycles = integer(field(repeat, "cycles")?, "repeat.cycles", 4096)?;
                        if cycles == 0 { return Err(bad("repeat.cycles must be positive")); }
                        Some(cycles)
                    }
                };
                Ok(Self::TransientPeak { repeated_cycles })
            }
            (Some("transient-sampled-peak"), None) => Err(bad("transient-sampled-peak requires a transient base request")),
            _ => Err(bad("qoi.kind must be steady-objective or transient-sampled-peak")),
        }
    }

    /// fs-uq binds this name in its existing checkpoint plan identity. The
    /// exact base request separately binds schedule, cycle count and mesh.
    pub(crate) fn plan_name(self) -> &'static str {
        match self {
            Self::Steady => "cooling-objective-temperature-k",
            Self::TransientPeak { repeated_cycles: None } => "cooling-transient-sampled-peak-temperature-k",
            Self::TransientPeak { repeated_cycles: Some(_) } => "cooling-repeated-transient-sampled-peak-temperature-k",
        }
    }

    pub(crate) fn extract(self, document: &J) -> Result<f64> {
        if document.str_field("schema") != Some("frankensim.cooling-network.result.v1") {
            return Err(model_failure("cooling sample returned an unexpected result schema"));
        }
        let final_value = result_number(document.path(&["objective", "value_k"]), "objective.value_k")?;
        match self {
            Self::Steady => {
                if document.get("transient").is_some() || document.get("repeated_cycles").is_some() {
                    return Err(model_failure("steady UQ received a transient result"));
                }
                Ok(final_value)
            }
            Self::TransientPeak { repeated_cycles } => {
                let cycle_peak = result_number(document.path(&["transient", "sampled_peak_objective_k"]),
                    "transient.sampled_peak_objective_k")?;
                if final_value <= 0.0 || cycle_peak < final_value {
                    return Err(model_failure("transient sampled peak is below its final objective or has invalid absolute temperatures"));
                }
                match repeated_cycles {
                    None => {
                        if document.get("repeated_cycles").is_some() {
                            return Err(model_failure("single-cycle UQ received a repeated-cycle result"));
                        }
                        Ok(cycle_peak)
                    }
                    Some(expected) => {
                        let repeated = document.get("repeated_cycles")
                            .ok_or_else(|| model_failure("missing all-cycle peak; never substitute the last cycle"))?;
                        let completed = repeated.get("cycles_completed").and_then(J::number_raw)
                            .and_then(|value| value.parse::<usize>().ok());
                        if repeated.str_field("status") != Some("fixed-count-complete") || completed != Some(expected) {
                            return Err(model_failure("cooling sample did not complete its declared fixed cycle count"));
                        }
                        let peak = result_number(repeated.get("sampled_peak_objective_k"),
                            "repeated_cycles.sampled_peak_objective_k")?;
                        if peak < cycle_peak {
                            return Err(model_failure("all-cycle peak is below the last cycle's peak"));
                        }
                        Ok(peak)
                    }
                }
            }
        }
    }

    /// Both result forms use the same temporal/spatial interpretation. The
    /// spatial label comes from the existing admitted objective selector.
    pub(crate) fn render(self, spatial: &str) -> String {
        match self {
            Self::Steady => format!("{{\"kind\":{},\"unit\":\"K\"}}", quote(spatial)),
            Self::TransientPeak { repeated_cycles } => format!(
                "{{\"kind\":\"transient-sampled-peak\",\"spatial_kind\":{},\"unit\":\"K\",\"observation\":\"one-completed-trajectory\",\"temporal_scope\":\"initial-state-and-accepted-endpoints\",\"cycles\":{},\"continuous_time_bound\":false}}",
                quote(spatial), repeated_cycles.unwrap_or(1),
            ),
        }
    }

    pub(crate) fn no_claim(self) -> &'static str {
        match self {
            Self::Steady => "fixed-count empirical propagation through actual cooling-network child solves; no confidence sequence, optional-stopping guarantee, physical/model-form uncertainty bound, mesh-convergence certificate, experimental validation, or native .fsim/ledger package claim",
            Self::TransientPeak { .. } => "fixed-count empirical propagation of completed cooling trajectories under a fixed schedule and horizon; peaks include initial states and accepted endpoints, not an inter-step or continuous-time maximum; adaptive endpoint locations may differ between samples; no time-discretization, physical/model-form uncertainty, mesh-convergence, experimental-validation, or native .fsim/ledger package claim; empirical standard errors and quantiles are not confidence sequences",
        }
    }
}

fn result_number(value: Option<&J>, name: &str) -> Result<f64> {
    value.and_then(J::as_f64).filter(|v| v.is_finite())
        .ok_or_else(|| model_failure(format!("cooling sample has no finite {name}")))
}

#[cfg(test)]
mod tests {
    use super::*;
    fn selected() -> J { J::parse(r#"{"kind":"transient-sampled-peak"}"#).unwrap() }
    fn document() -> J {
        J::parse(r#"{"schema":"frankensim.cooling-network.result.v1","objective":{"value_k":301},"transient":{"sampled_peak_objective_k":350},"repeated_cycles":{"status":"fixed-count-complete","cycles_completed":3,"sampled_peak_objective_k":400}}"#).unwrap()
    }

    #[test]
    fn transient_selection_is_explicit_and_horizon_is_fixed() {
        let base = J::parse(r#"{"transient":{"repeat":{"cycles":3}}}"#).unwrap();
        assert!(Qoi::parse(None, &base).is_err());
        assert_eq!(Qoi::parse(Some(&selected()), &base).unwrap(),
            Qoi::TransientPeak { repeated_cycles: Some(3) });
        for text in [r#"{}"#, r#"{"transient":{"power_design":{}}}"#,
            r#"{"transient":{"repeat":{"until_periodic":{}}}}"#,
            r#"{"transient":{"repeat":{"cycles":0}}}"#] {
            assert!(Qoi::parse(Some(&selected()), &J::parse(text).unwrap()).is_err());
        }
    }

    #[test]
    fn all_cycle_peak_is_never_replaced_by_final_or_last_cycle_temperature() {
        let qoi = Qoi::TransientPeak { repeated_cycles: Some(3) };
        assert_eq!(qoi.extract(&document()).unwrap(), 400.0);
        assert!(Qoi::Steady.extract(&document()).is_err());
        assert!(Qoi::TransientPeak { repeated_cycles: None }.extract(&document()).is_err());
        assert!(Qoi::TransientPeak { repeated_cycles: Some(2) }.extract(&document()).is_err());
        let single = J::parse(r#"{"schema":"frankensim.cooling-network.result.v1","objective":{"value_k":301},"transient":{"sampled_peak_objective_k":350}}"#).unwrap();
        assert_eq!(Qoi::TransientPeak { repeated_cycles: None }.extract(&single).unwrap(), 350.0);
        assert!(qoi.extract(&single).is_err());
        let missing = J::parse(r#"{"schema":"frankensim.cooling-network.result.v1","objective":{"value_k":301}}"#).unwrap();
        assert!(Qoi::TransientPeak { repeated_cycles: None }.extract(&missing).is_err());
    }
}
