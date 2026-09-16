//! Predeclared design coordinates applied AFTER the existing uncertainty draw.
//! A multiplier scales the sampled physical input, not the reported temperature.
use super::*;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum DesignControl { FanSpeed, WorkloadPower }

#[derive(Debug, Clone)]
pub(crate) struct DesignGrid {
    pub(crate) control: DesignControl,
    /// Strictly increasing multipliers. Fan preference is ascending; workload
    /// preference is descending. These are discrete alternatives, not a bracket.
    pub(crate) multipliers: Vec<f64>,
}

impl DesignGrid {
    pub(crate) fn parse(control: DesignControl, text: &str) -> Result<Self> {
        let mut multipliers = Vec::new();
        for word in text.split(',') {
            if multipliers.len() == 64 { return Err(bad("at most 64 design candidates are admitted")); }
            let value = word.trim().parse::<f64>()
                .map_err(|_| bad("design candidates must be comma-separated finite multipliers"))?;
            if !value.is_finite() || value < 0.0
                || (control == DesignControl::FanSpeed && value == 0.0)
                || multipliers.last().is_some_and(|previous| value <= *previous) {
                return Err(bad("design multipliers must be strictly increasing and nonnegative; fan multipliers must be positive"));
            }
            multipliers.push(value);
        }
        if multipliers.len() < 2 { return Err(bad("design selection requires at least two predeclared candidates")); }
        Ok(Self { control, multipliers })
    }

    pub(crate) fn label(&self) -> &'static str {
        match self.control {
            DesignControl::FanSpeed => "fan-speed-multiplier",
            DesignControl::WorkloadPower => "workload-power-multiplier",
        }
    }

    pub(crate) fn preference_order(&self) -> Vec<usize> {
        match self.control {
            DesignControl::FanSpeed => (0..self.multipliers.len()).collect(),
            DesignControl::WorkloadPower => (0..self.multipliers.len()).rev().collect(),
        }
    }

    /// Nominal shape/arithmetic admission before allocating executions or
    /// reserving output. The actual sampled model still performs all domain,
    /// correlation and physical acceptance checks on EVERY evaluation.
    pub(crate) fn validate(&self, base: &J, samples_per_candidate: usize) -> Result<()> {
        let total = self.multipliers.len().checked_mul(samples_per_candidate)
            .ok_or_else(|| budget("candidate sample budget overflow"))?;
        if total > MAX_PRODUCT_SAMPLES {
            return Err(bad(format!("candidate count times samples must not exceed {MAX_PRODUCT_SAMPLES}")));
        }
        if base.get("transient").is_some_and(|s| s.get("adjoint").is_some()) {
            return Err(bad("uncertainty-aware selection requires forward trajectories, not per-sample adjoints"));
        }
        for &multiplier in &self.multipliers {
            let mut candidate = base.clone();
            scale_candidate(&mut candidate, self.control, multiplier)?;
        }
        Ok(())
    }
}

impl Config {
    pub(crate) fn sample_candidate_request(&self, base: &J, values: &[f64],
        control: DesignControl, multiplier: f64) -> Result<String> {
        // Reuse target/dependence/source-total lowering rather than creating a
        // second sampler. Apply design last so a sampled speed/load is not lost.
        let sampled = self.sample_request(base, values)?;
        let mut candidate = J::parse(&sampled).map_err(|error| model_failure(error.to_string()))?;
        scale_candidate(&mut candidate, control, multiplier)?;
        render_json(&candidate)
    }
}

fn scale_candidate(root: &mut J, control: DesignControl, multiplier: f64) -> Result<()> {
    if !multiplier.is_finite() || multiplier < 0.0
        || (control == DesignControl::FanSpeed && multiplier == 0.0) {
        return Err(bad("invalid design multiplier"));
    }
    if control == DesignControl::FanSpeed {
        let fan = field(field(root, "hydraulics")?, "fan")?;
        let minimum = positive(field(fan, "min_speed_ratio")?, "fan.min_speed_ratio")?;
        let maximum = positive(field(fan, "max_speed_ratio")?, "fan.max_speed_ratio")?;
        let speed = |base: f64| -> Result<f64> {
            let value = base * multiplier;
            if !value.is_finite() || value <= 0.0 || value < minimum || value > maximum {
                return Err(model_failure("candidate fan speed leaves the declared speed domain; no clipping or redraw"));
            }
            Ok(value)
        };
        if root.get("transient").is_some() {
            for interval in array_mut_path(root, &["transient", "intervals"])? {
                let value = speed(positive(field(interval,"fan_speed_ratio")?,"interval fan speed")?)?;
                set_member_number(interval, "fan_speed_ratio", value)?;
            }
        } else {
            let value = speed(positive(field(fan,"speed_ratio")?,"fan speed")?)?;
            set_path_number(root, &["hydraulics","fan","speed_ratio"], value)?;
        }
        return Ok(());
    }
    if root.get("transient").is_some() {
        for interval in array_mut_path(root, &["transient", "intervals"])? {
            match (interval.get("power_scale").is_some(), interval.get("component_powers_w").is_some()) {
                (true, false) => multiply_member(interval, "power_scale", multiplier)?,
                (false, true) => {
                    let J::Object(powers) = member_mut(interval, "component_powers_w")?
                        else { return Err(bad("component_powers_w must be an object")); };
                    for (_, value) in powers { multiply_number(value, multiplier)?; }
                }
                _ => return Err(bad("each transient design interval requires exactly one workload form")),
            }
        }
    } else {
        let solid = member_mut(root, "solid")?;
        if solid.get("component_power").is_none() {
            return Err(bad("steady workload candidate selection requires solid.component_power"));
        }
        for component in array_mut_path(solid, &["component_power", "components"])? {
            multiply_member(component, "watts", multiplier)?;
        }
        recompute_component_total(root)?;
    }
    Ok(())
}

fn multiply_member(root: &mut J, name: &str, multiplier: f64) -> Result<()> {
    multiply_number(member_mut(root, name)?, multiplier)
}
fn multiply_number(slot: &mut J, multiplier: f64) -> Result<()> {
    let value = number(slot, "candidate workload")? * multiplier;
    if !value.is_finite() { return Err(model_failure("candidate workload overflowed")); }
    *slot = J::Number { value, raw: value.to_string() };
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn candidate_order_is_explicit_and_ambiguous_grids_refuse() {
        for text in ["", "1", "1,1", "2,1", "0,1", "1,NaN", "-1,1"] {
            assert!(DesignGrid::parse(DesignControl::FanSpeed,text).is_err());
        }
        let grid = DesignGrid::parse(DesignControl::WorkloadPower,"0, 1, 2").unwrap();
        assert_eq!(grid.preference_order(), vec![2,1,0]);
        assert_eq!(DesignGrid::parse(DesignControl::FanSpeed,"0.5,1,1.5").unwrap()
            .preference_order(), vec![0,1,2]);
    }

    #[test]
    fn schedule_scaling_changes_active_controls_once_and_keeps_zero_loads() {
        let base = J::parse(r#"{"hydraulics":{"fan":{"min_speed_ratio":0.5,"max_speed_ratio":2,"speed_ratio":1}},"transient":{"repeat":{"cycles":3},"intervals":[{"fan_speed_ratio":1,"power_scale":2},{"fan_speed_ratio":1.5,"component_powers_w":{"chip":0,"memory":4}}]}}"#).unwrap();
        let mut power = base.clone();
        scale_candidate(&mut power,DesignControl::WorkloadPower,0.25).unwrap();
        let rows = array_path(&power,&["transient","intervals"]).unwrap();
        assert_eq!(rows[0].f64_field("power_scale"),Some(0.5));
        assert_eq!(rows[1].path(&["component_powers_w","chip"]).unwrap().as_f64(),Some(0.0));
        assert_eq!(rows[1].path(&["component_powers_w","memory"]).unwrap().as_f64(),Some(1.0));
        let mut fan = base.clone();
        scale_candidate(&mut fan,DesignControl::FanSpeed,0.75).unwrap();
        assert_eq!(array_path(&fan,&["transient","intervals"]).unwrap()[1]
            .f64_field("fan_speed_ratio"),Some(1.125));
        assert_eq!(fan.path(&["hydraulics","fan","speed_ratio"]),base.path(&["hydraulics","fan","speed_ratio"]));
        assert_eq!(fan.path(&["transient","repeat"]),base.path(&["transient","repeat"]));
        assert!(scale_candidate(&mut fan,DesignControl::FanSpeed,4.0).is_err());
    }
}
