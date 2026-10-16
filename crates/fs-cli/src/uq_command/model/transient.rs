//! Uncertainty in quantities actually consumed by the transient producer.
//! A draw is made once per trajectory, not independently at every timestep.
//! Fixed schedules repeat the same draw in every cycle. Durations, footprints,
//! timestep rules, controllers and material assignments remain unchanged.
use super::*;

#[derive(Debug, Clone)]
pub(super) enum Target {
    InitialTemperature,
    HeatCapacity,
    ElementHeatCapacity(usize),
    IntervalPowerScale(usize),
    IntervalComponentPower { interval: usize, component: String },
    IntervalFanSpeed(usize),
}

impl Target {
    pub(super) fn parse(value: &J) -> Result<Option<Self>> {
        let kind = value.str_field("kind");
        let (target, allowed): (Self, &[&str]) = match kind {
            Some("initial-temperature") => (Self::InitialTemperature, &["kind"]),
            Some("volumetric-heat-capacity") => (Self::HeatCapacity, &["kind"]),
            Some("element-heat-capacity") => (
                Self::ElementHeatCapacity(integer_raw(field(value, "element")?, "element")?),
                &["kind", "element"],
            ),
            Some("interval-power-scale") => (
                Self::IntervalPowerScale(integer_raw(field(value, "interval")?, "interval")?),
                &["kind", "interval"],
            ),
            Some("interval-component-power") => (
                Self::IntervalComponentPower {
                    interval: integer_raw(field(value, "interval")?, "interval")?,
                    component: string(field(value, "component")?, "component")?,
                }, &["kind", "interval", "component"],
            ),
            Some("interval-fan-speed-ratio") => (
                Self::IntervalFanSpeed(integer_raw(field(value, "interval")?, "interval")?),
                &["kind", "interval"],
            ),
            _ => return Ok(None),
        };
        object(value, allowed, "transient uncertainty target")?;
        Ok(Some(target))
    }

    pub(super) fn name(&self) -> String {
        match self {
            Self::InitialTemperature => "transient.initial_temperature_k".into(),
            Self::HeatCapacity => "transient.volumetric_heat_capacity_j_m3_k".into(),
            Self::ElementHeatCapacity(index) => format!("transient.element_heat_capacities_j_m3_k[{index}]"),
            Self::IntervalPowerScale(index) => format!("transient.intervals[{index}].power_scale"),
            Self::IntervalComponentPower { interval, component } =>
                format!("transient.intervals[{interval}].component[{component}].watts"),
            Self::IntervalFanSpeed(index) => format!("transient.intervals[{index}].fan_speed_ratio"),
        }
    }

    pub(super) fn unit(&self) -> &'static str {
        match self {
            Self::InitialTemperature => "K",
            Self::HeatCapacity | Self::ElementHeatCapacity(_) => "J/(m3 K)",
            Self::IntervalPowerScale(_) | Self::IntervalFanSpeed(_) => "1",
            Self::IntervalComponentPower { .. } => "W",
        }
    }

    pub(super) fn allows_zero(&self) -> bool {
        matches!(self, Self::IntervalPowerScale(_) | Self::IntervalComponentPower { .. })
    }

    pub(super) fn render(&self) -> String {
        match self {
            Self::InitialTemperature => "{\"kind\":\"initial-temperature\"}".into(),
            Self::HeatCapacity => "{\"kind\":\"volumetric-heat-capacity\"}".into(),
            Self::ElementHeatCapacity(index) => format!("{{\"kind\":\"element-heat-capacity\",\"element\":{index}}}"),
            Self::IntervalPowerScale(index) => format!("{{\"kind\":\"interval-power-scale\",\"interval\":{index}}}"),
            Self::IntervalComponentPower { interval, component } => format!(
                "{{\"kind\":\"interval-component-power\",\"interval\":{interval},\"component\":{}}}", quote(component)),
            Self::IntervalFanSpeed(index) => format!("{{\"kind\":\"interval-fan-speed-ratio\",\"interval\":{index}}}"),
        }
    }

    pub(super) fn validate(&self, base: &J) -> Result<()> {
        let schedule = field(base, "transient")?;
        let selected = match self {
            Self::InitialTemperature => {
                exclusive(schedule, "initial_temperature_k", "initial_temperatures_k")?
            }
            Self::HeatCapacity => {
                exclusive(schedule, "volumetric_heat_capacity_j_m3_k", "element_heat_capacities_j_m3_k")?
            }
            Self::ElementHeatCapacity(index) => {
                let values = exclusive(schedule, "element_heat_capacities_j_m3_k", "volumetric_heat_capacity_j_m3_k")?
                    .as_array().ok_or_else(|| bad("element heat capacities must be an array"))?;
                values.get(*index).ok_or_else(|| bad("element heat-capacity index is out of range"))?
            }
            Self::IntervalPowerScale(index) => exclusive(interval(base, *index)?, "power_scale", "component_powers_w")?,
            Self::IntervalComponentPower { interval: index, component } => {
                component_location(base, component)?;
                let powers = exclusive(interval(base, *index)?, "component_powers_w", "power_scale")?;
                field(powers, component)?
            }
            Self::IntervalFanSpeed(index) => {
                field(field(base, "hydraulics")?, "fan")?;
                field(interval(base, *index)?, "fan_speed_ratio")?
            }
        };
        let nominal = number(selected, "nominal transient target")?;
        if nominal < 0.0 || (!self.allows_zero() && nominal == 0.0) {
            return Err(bad("nominal transient target leaves its admitted physical domain"));
        }
        Ok(())
    }

    pub(super) fn apply(&self, base: &mut J, value: f64) -> Result<()> {
        if !value.is_finite() || value < 0.0 || (!self.allows_zero() && value == 0.0) {
            return Err(model_failure("sampled transient target leaves its admitted physical domain"));
        }
        match self {
            Self::InitialTemperature => set_path_number(base, &["transient", "initial_temperature_k"], value),
            Self::HeatCapacity => set_path_number(base, &["transient", "volumetric_heat_capacity_j_m3_k"], value),
            Self::ElementHeatCapacity(index) => {
                let values = array_mut_path(base, &["transient", "element_heat_capacities_j_m3_k"])?;
                let slot = values.get_mut(*index).ok_or_else(|| bad("element heat-capacity index disappeared"))?;
                *slot = J::Number { value, raw: value.to_string() };
                Ok(())
            }
            Self::IntervalPowerScale(index) => set_member_number(interval_mut(base, *index)?, "power_scale", value),
            Self::IntervalComponentPower { interval: index, component } => {
                let powers = member_mut(interval_mut(base, *index)?, "component_powers_w")?;
                set_member_number(powers, component, value)
            }
            Self::IntervalFanSpeed(index) => set_member_number(interval_mut(base, *index)?, "fan_speed_ratio", value),
        }
    }
}

fn exclusive<'a>(value: &'a J, selected: &str, other: &str) -> Result<&'a J> {
    if value.get(other).is_some() {
        return Err(bad(format!("target {selected} cannot replace or coexist with {other}")));
    }
    field(value, selected)
}
fn interval(base: &J, index: usize) -> Result<&J> {
    array_path(base, &["transient", "intervals"])?.get(index)
        .ok_or_else(|| bad("transient interval index is out of range"))
}
fn interval_mut(base: &mut J, index: usize) -> Result<&mut J> {
    array_mut_path(base, &["transient", "intervals"])?.get_mut(index)
        .ok_or_else(|| bad("transient interval index disappeared"))
}

#[cfg(test)]
mod tests {
    use super::*;
    const BASE: &str = r#"{"hydraulics":{"fan":{"speed_ratio":1}},"solid":{"component_power":{"total_w":5,"components":[{"name":"chip","watts":5,"vertices":[0]}]}},"transient":{"initial_temperature_k":300,"element_heat_capacities_j_m3_k":[1000,2000],"intervals":[{"duration_s":3,"power_scale":1,"fan_speed_ratio":1},{"duration_s":7,"component_powers_w":{"chip":2},"fan_speed_ratio":1.5}]}}"#;

    #[test]
    fn sample_changes_only_active_declared_coordinates() {
        let mut base = J::parse(BASE).unwrap();
        let original = base.clone();
        for (target, value) in [
            (Target::InitialTemperature, 310.0),
            (Target::ElementHeatCapacity(1), 2500.0),
            (Target::IntervalPowerScale(0), 0.0),
            (Target::IntervalComponentPower { interval: 1, component: "chip".into() }, 9.0),
            (Target::IntervalFanSpeed(1), 1.8),
        ] {
            target.validate(&base).unwrap();
            target.apply(&mut base, value).unwrap();
            let roundtrip = Target::parse(&J::parse(&target.render()).unwrap()).unwrap().unwrap();
            assert_eq!(roundtrip.name(), target.name());
        }
        assert_eq!(base.path(&["solid"]), original.path(&["solid"]));
        assert_eq!(base.path(&["hydraulics"]), original.path(&["hydraulics"]));
        assert_eq!(base.path(&["transient", "initial_temperature_k"]).and_then(J::as_f64), Some(310.0));
        let capacities = array_path(&base, &["transient", "element_heat_capacities_j_m3_k"]).unwrap();
        assert_eq!(capacities[0].as_f64(), Some(1000.0));
        assert_eq!(capacities[1].as_f64(), Some(2500.0));
        assert_eq!(interval(&base, 0).unwrap().f64_field("power_scale"), Some(0.0));
        assert_eq!(interval(&base, 1).unwrap().path(&["component_powers_w", "chip"]).and_then(J::as_f64), Some(9.0));
        assert_eq!(interval(&base, 1).unwrap().f64_field("fan_speed_ratio"), Some(1.8));
        for index in 0..2 {
            assert_eq!(interval(&base, index).unwrap().get("duration_s"), interval(&original, index).unwrap().get("duration_s"));
        }
    }

    #[test]
    fn mode_mismatches_inactive_targets_and_bad_indices_refuse() {
        let base = J::parse(BASE).unwrap();
        for target in [Target::HeatCapacity, Target::ElementHeatCapacity(2),
            Target::IntervalPowerScale(1), Target::IntervalFanSpeed(2),
            Target::IntervalComponentPower { interval: 0, component: "chip".into() },
            Target::IntervalComponentPower { interval: 1, component: "unknown".into() }] {
            assert!(target.validate(&base).is_err(), "accepted {}", target.name());
        }
        assert!(Target::parse(&J::parse(r#"{"kind":"initial-temperature","interval":0}"#).unwrap()).is_err());
        assert!(Target::IntervalFanSpeed(0).validate(&J::parse("{}").unwrap()).is_err());
    }

    #[test]
    fn uniform_capacity_is_supported_without_erasing_element_assignments() {
        let mut base = J::parse(r#"{"transient":{"volumetric_heat_capacity_j_m3_k":2000}}"#).unwrap();
        Target::HeatCapacity.validate(&base).unwrap();
        Target::HeatCapacity.apply(&mut base, 3000.0).unwrap();
        assert_eq!(base.path(&["transient", "volumetric_heat_capacity_j_m3_k"]).and_then(J::as_f64), Some(3000.0));
        assert!(Target::ElementHeatCapacity(0).validate(&base).is_err());
    }

    #[test]
    fn nonphysical_samples_are_not_clipped_or_written() {
        let mut base = J::parse(BASE).unwrap();
        let before = base.clone();
        for value in [0.0, -1.0, f64::INFINITY, f64::NAN] {
            assert!(Target::InitialTemperature.apply(&mut base, value).is_err());
        }
        assert_eq!(base, before);
        Target::IntervalComponentPower { interval: 1, component: "chip".into() }.apply(&mut base, 0.0).unwrap();
    }
}
