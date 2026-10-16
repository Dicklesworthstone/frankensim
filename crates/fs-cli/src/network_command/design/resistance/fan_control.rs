//! Fan speed as the final coordinate of the joint cooling design, never a
//! nested sizing search. All trials are relative to the original speed.

use super::*;

#[derive(Debug)]
pub(super) struct FanControl {
    baseline: f64,
    minimum: f64,
    maximum: f64,
    pub(super) weight: f64,
}

impl FanControl {
    pub(super) fn parse(value: &J, request: &Request) -> Result<Self> {
        object(value, &["minimum_speed_ratio", "maximum_speed_ratio", "weight"], "fan_speed_control")?;
        let fan = request.fan.as_ref().ok_or_else(|| bad("fan_speed_control requires hydraulics.fan"))?;
        let control = Self {
            baseline: fan.speed_ratio,
            minimum: positive(get(value, "minimum_speed_ratio")?, "minimum_speed_ratio")?,
            maximum: positive(get(value, "maximum_speed_ratio")?, "maximum_speed_ratio")?,
            weight: positive(get(value, "weight")?, "fan_speed_control.weight")?,
        };
        if control.minimum >= control.maximum || control.minimum > control.baseline
            || control.maximum < control.baseline
        {
            return Err(bad("fan speed bounds must contain the original speed and have positive width"));
        }
        // The source-owned fan domain, not the optimizer box, admits speeds.
        fan.bank(control.minimum)?;
        fan.bank(control.maximum)?;
        let (lo, hi) = control.bounds();
        if !lo.is_finite() || !hi.is_finite() || lo >= hi
            || !det::exp(lo).is_finite() || det::exp(lo) <= 0.0
            || !det::exp(hi).is_finite()
        {
            return Err(bad("fan speed log-multiplier bounds are not representable"));
        }
        Ok(control)
    }

    pub(super) fn bounds(&self) -> (f64, f64) {
        (det::ln(self.minimum / self.baseline), det::ln(self.maximum / self.baseline))
    }

    pub(super) fn candidate(&self, request: &Request, coordinate: f64) -> Result<f64> {
        let (lo, hi) = self.bounds();
        if !coordinate.is_finite() || coordinate < lo || coordinate > hi {
            return Err(bad("fan speed coordinate is outside its declared box"));
        }
        let speed = finite(self.baseline * det::exp(coordinate))?.clamp(self.minimum, self.maximum);
        request.fan.as_ref().ok_or_else(|| bad("joint design lost its fan drive"))?.bank(speed)?;
        Ok(speed)
    }

    pub(super) fn render(&self, request: &Request, coordinate: f64, index: usize) -> Result<String> {
        let speed = request.fan.as_ref().ok_or_else(|| bad("joint design lost its fan drive"))?.speed_ratio;
        Ok(format!("{{\"coordinate_index\":{index},\"baseline_speed_ratio\":{},\"selected_speed_ratio\":{},\"multiplier\":{},\"log_multiplier\":{},\"minimum_speed_ratio\":{},\"maximum_speed_ratio\":{},\"weight\":{}}}",
            num(self.baseline)?, num(speed)?, num(speed / self.baseline)?, num(coordinate)?,
            num(self.minimum)?, num(self.maximum)?, num(self.weight)?))
    }
}
