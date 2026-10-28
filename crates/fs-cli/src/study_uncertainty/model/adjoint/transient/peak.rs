//! Bind a frozen nominal control to the native selected sampled-peak branch.
//! Ties remain usable fixed controls; they are not unique maximum derivatives.

use crate::json_read::JsonValue as J;
use super::super::super::peak as observation;

#[derive(Debug, Clone, PartialEq)]
pub(super) struct SelectedPeak {
    pub(super) state: usize,
    pub(super) time_s: f64,
    value_k: f64,
    vertex: usize,
    ties: usize,
    fine_steps: usize,
}

impl SelectedPeak {
    pub(super) fn json(&self) -> String {
        format!(
            "{{\"objective\":\"temperature-peak\",\"value_k\":{},\"time_s\":{},\"state_index\":{},\"active_vertex\":{},\"active_ties\":{},\"unique_sampled_maximum\":{},\"control_prefix_steps\":{},\"fine_steps\":{},\"coefficient_semantics\":\"frozen-selected-branch\",\"unique_maximum_derivative_claim\":false}}",
            self.value_k, self.time_s, self.state, self.vertex, self.ties,
            self.ties == 1, self.state, self.fine_steps,
        )
    }
}

fn count(row: &J, key: &str) -> Result<usize, &'static str> {
    row.get(key).and_then(J::number_raw).and_then(|value| value.parse().ok())
        .ok_or("peak calibration has an invalid integer index or work count")
}

fn number(row: &J, key: &str) -> Result<f64, &'static str> {
    row.f64_field(key).filter(|value| value.is_finite())
        .ok_or("peak calibration has a missing or nonfinite physical value")
}

/// The caller has already verified the sealed child and its physical inputs.
/// Bind its branch to the same complete forward observation, fine row and
/// reverse prefix. Only metadata is copied; no field or trajectory is cloned.
pub(super) fn bind(
    receipt: &J, nominal: f64, region: &str, fine_steps: usize, enthalpy: bool,
) -> Result<SelectedPeak, &'static str> {
    let run = receipt.str_field("run").filter(|value| !value.is_empty())
        .ok_or("peak calibration has no sealed child identity")?;
    let time = receipt.get("transient").ok_or("peak calibration has no trajectory")?;
    let goal = time.get("temporal_error").ok_or("peak calibration has no final observation")?;
    let final_k = number(goal, "fine_final_k")?;
    let value_k = observation::temperature(receipt, run, region, final_k)?;
    let peak = time.get("sampled_peak").ok_or("peak calibration has no forward peak")?;
    let report = time.get("adjoint").ok_or("peak calibration has no native branch adjoint")?;
    let state = count(peak, "state")?;
    let time_s = number(peak, "time_s")?;
    let vertex = count(report, "active_vertex")?;
    let ties = count(report, "active_ties")?;
    let initial_key = if enthalpy {
        "dtemperature_dinitial_specific_enthalpy_k_kg_per_j"
    } else {
        "dtemperature_dinitial_temperatures"
    };
    let initial = report.get(initial_key).and_then(J::as_array)
        .ok_or("peak branch has no native vertex layout")?;
    let rows = time.get("fine").and_then(J::as_array)
        .filter(|rows| rows.len() == fine_steps)
        .ok_or("peak calibration does not retain its complete fine grid")?;
    if value_k.to_bits() != nominal.to_bits()
        || count(time, "fine_steps")? != fine_steps
        || count(report, "state_index")? != state
        || number(report, "time_s")?.to_bits() != time_s.to_bits()
        || number(report, "value_k")?.to_bits() != value_k.to_bits()
        || count(report, "max_reverse_steps")? != fine_steps
        || ties == 0 || vertex >= initial.len()
        || report.get("unique_sampled_maximum") != Some(&J::Bool(ties == 1))
    {
        return Err("nominal branch differs from the selected forward peak or its vertex layout");
    }
    let mut previous = 0.0;
    for (index, row) in rows.iter().enumerate() {
        let at = number(row, "time_s")?;
        let value = number(row, "final_region_max_k")?;
        if count(row, "step")? != index + 1 || at <= previous || value > value_k
            || (index + 1 < state && value == value_k)
            || (index + 1 == state && (at.to_bits() != time_s.to_bits()
                || value.to_bits() != value_k.to_bits()))
        {
            return Err("peak calibration selects another fine state or a later tied maximum");
        }
        previous = at;
    }
    if previous.to_bits() != number(time, "final_time_s")?.to_bits()
        || rows.last().and_then(|row| row.f64_field("final_region_max_k"))
            .map(f64::to_bits) != Some(final_k.to_bits())
    {
        return Err("peak calibration fine grid and final observation differ");
    }
    let revalidated = count(report, "revalidated_solid_endpoints")?;
    let reconstructed = if enthalpy { 0 } else { count(report, "reconstructed_solid_endpoints")? };
    if revalidated.checked_add(reconstructed) != Some(state) {
        return Err("peak adjoint did not reverse the selected physical prefix");
    }
    for name in ["material_controls", "prescribed_temperature", "prescribed_heat_flux",
        "prescribed_convection", "convection", "radiation"] {
        if let Some(controls) = report.get(name) {
            for key in ["control_steps", "air_control_steps", "natural_control_steps",
                "physical_endpoint_rechecks"] {
                if controls.get(key).is_some() && count(controls, key)? != state {
                    return Err("peak control family belongs to another trajectory prefix");
                }
            }
        }
    }
    Ok(SelectedPeak { state, time_s, value_k, vertex, ties, fine_steps })
}

#[cfg(test)]
mod tests {
    use super::*;

    const SOURCE: &str = r#"{"run":"nominal-child","transient":{
        "status":"completed","final_time_s":2,"fine_steps":4,
        "temporal_error":{"fine_final_k":300},
        "sampled_peak":{"region":"solid","value_k":320,"time_s":1,
            "state":2,"samples":5,"includes_initial":true,"time_scope":"accepted-fine-grid",
            "continuous_time_bound":false,"authority":"Estimated"},
        "fine":[{"step":1,"time_s":0.5,"final_region_max_k":310},
            {"step":2,"time_s":1,"final_region_max_k":320},
            {"step":3,"time_s":1.5,"final_region_max_k":315},
            {"step":4,"time_s":2,"final_region_max_k":300}],
        "adjoint":{"value_k":320,"time_s":1,"state_index":2,
            "active_vertex":1,"active_ties":1,"unique_sampled_maximum":true,
            "max_reverse_steps":4,"revalidated_solid_endpoints":2,
            "reconstructed_solid_endpoints":0,"dtemperature_dinitial_temperatures":[0.1,0.2],
            "prescribed_convection":{"control_steps":2}}}}"#;

    fn read(source: &str) -> Result<SelectedPeak, &'static str> {
        bind(&J::parse(source).unwrap(), 320.0, "solid", 4, false)
    }

    #[test]
    fn g0_peak_branch_binds_an_interior_maximum_and_selected_prefix() {
        let selected = read(SOURCE).unwrap();
        assert_eq!(selected.state, 2);
        assert_eq!(selected.time_s, 1.0);
        let encoded = J::parse(&selected.json()).unwrap();
        assert_eq!(encoded.f64_field("control_prefix_steps"), Some(2.0));
        assert_eq!(encoded.get("unique_sampled_maximum"), Some(&J::Bool(true)));
    }

    #[test]
    fn g0_tied_branch_remains_an_explicit_frozen_control_without_unique_derivative_claim() {
        let source = SOURCE.replace("\"active_ties\":1", "\"active_ties\":2")
            .replace("\"unique_sampled_maximum\":true", "\"unique_sampled_maximum\":false");
        let selected = read(&source).unwrap();
        let encoded = J::parse(&selected.json()).unwrap();
        assert_eq!(encoded.f64_field("active_ties"), Some(2.0));
        assert_eq!(encoded.str_field("coefficient_semantics"), Some("frozen-selected-branch"));
        assert_eq!(encoded.get("unique_maximum_derivative_claim"), Some(&J::Bool(false)));
        assert!(read(&source.replace("\"unique_sampled_maximum\":false",
            "\"unique_sampled_maximum\":true")).is_err());
    }

    #[test]
    fn g0_initial_tied_peak_has_an_empty_reverse_prefix_for_both_storage_laws() {
        let source = SOURCE.replace("\"value_k\":320,\"time_s\":1", "\"value_k\":320,\"time_s\":0")
            .replace("\"state\":2", "\"state\":0")
            .replace("\"state_index\":2", "\"state_index\":0")
            .replace("\"revalidated_solid_endpoints\":2", "\"revalidated_solid_endpoints\":0")
            .replace("\"control_steps\":2", "\"control_steps\":0")
            .replace("\"active_ties\":1", "\"active_ties\":3")
            .replace("\"unique_sampled_maximum\":true", "\"unique_sampled_maximum\":false");
        for enthalpy in [false, true] {
            let source = if enthalpy {
                source.replace("dtemperature_dinitial_temperatures",
                    "dtemperature_dinitial_specific_enthalpy_k_kg_per_j")
            } else { source.clone() };
            let selected = bind(&J::parse(&source).unwrap(), 320.0, "solid", 4, enthalpy).unwrap();
            assert_eq!(selected.state, 0);
            assert_eq!(selected.time_s, 0.0);
            assert_eq!(selected.ties, 3);
        }
    }

    #[test]
    fn g0_peak_calibration_rejects_mismatched_state_value_grid_and_reverse_work() {
        for (old, new) in [
            ("\"state_index\":2", "\"state_index\":3"),
            ("\"state_index\":2", "\"state_index\":2.0"),
            ("\"active_vertex\":1", "\"active_vertex\":2"),
            ("\"active_ties\":1", "\"active_ties\":0"),
            ("\"control_steps\":2", "\"control_steps\":4"),
            ("\"revalidated_solid_endpoints\":2", "\"revalidated_solid_endpoints\":4"),
            ("\"max_reverse_steps\":4", "\"max_reverse_steps\":2"),
            ("\"time_s\":1,\"state_index\"", "\"time_s\":1.5,\"state_index\""),
            ("\"value_k\":320,\"time_s\":1,\"state_index\"", "\"value_k\":319,\"time_s\":1,\"state_index\""),
            ("\"final_region_max_k\":310", "\"final_region_max_k\":320"),
            ("\"final_region_max_k\":315", "\"final_region_max_k\":321"),
            ("\"step\":3", "\"step\":4"),
            ("\"time_s\":1.5", "\"time_s\":0.75"),
        ] {
            assert!(SOURCE.contains(old), "{old}");
            assert!(read(&SOURCE.replace(old, new)).is_err(), "accepted {new}");
        }
    }
}
