//! Immutable updates of authored power values for physical inverse studies.
use super::{MovingHeatSource, ConductionError, invalid};
use fs_scenario::signal::TimeSignal;

impl MovingHeatSource {
    /// Copy the original constant, ramp endpoints, or table power values [W].
    /// Chebyshev coefficients are not pointwise nonnegative power controls and
    /// are deliberately refused by this positivity-preserving parameterization.
    pub fn absorbed_power_values(&self, max_values: usize) -> Result<Vec<f64>, ConductionError> {
        let count = match &self.absorbed_power_w {
            TimeSignal::Constant(_) => 1,
            TimeSignal::Ramp { .. } => 2,
            TimeSignal::Table { values, .. } => values.len(),
            TimeSignal::Chebfun(_) => return Err(invalid("power-value design supports constant/ramp/table signals; polynomial coefficients need a separate positivity policy")),
        };
        if count > max_values { return Err(invalid("authored power values exceed their explicit control allowance")); }
        let mut values = Vec::new();
        values.try_reserve_exact(count).map_err(|_| invalid("power-control allocation refused"))?;
        match &self.absorbed_power_w {
            TimeSignal::Constant(value) => values.push(value.value),
            TimeSignal::Ramp { from, to, .. } => values.extend([from.value, to.value]),
            TimeSignal::Table { values: original, .. } => values.extend_from_slice(original),
            TimeSignal::Chebfun(_) => unreachable!("rejected before allocation"),
        }
        Ok(values)
    }

    /// Rebuild a source with new nonnegative absorbed-power values, preserving
    /// original time knots, interpolation, position, footprint and spatial caps.
    /// The original declaration is unchanged on success or refusal. Zero/off
    /// values are valid; selected positive log-power controls are an L4 policy.
    pub fn with_absorbed_power_values(&self, values: &[f64], reference_time_s: f64,
        max_values: usize) -> Result<Self, ConductionError>
    {
        let original = self.absorbed_power_values(max_values)?;
        if values.len() != original.len() || values.iter().any(|v| !v.is_finite() || *v < 0.0) {
            return Err(invalid("replacement power values must preserve shape and be finite nonnegative watts"));
        }
        let mut power = self.absorbed_power_w.clone();
        match &mut power {
            TimeSignal::Constant(value) => value.value = values[0],
            TimeSignal::Ramp { from, to, .. } => { from.value = values[0]; to.value = values[1]; }
            TimeSignal::Table { values: out, .. } => out.copy_from_slice(values),
            TimeSignal::Chebfun(_) => unreachable!("rejected before replacement"),
        }
        Self::new(self.center_m.clone(), self.half_width_m, power, reference_time_s, self.budget)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bc::history::march::moving::MovingHeatBudget;
    use fs_qty::{Dims, QtyAny};
    use fs_scenario::signal::Interp;
    #[test]
    fn g0_power_design_preserves_the_path_knots_off_state_and_original_source() {
        let source = MovingHeatSource::new(std::array::from_fn(|_| TimeSignal::Constant(
            QtyAny::new(0.0, Dims([1,0,0,0,0,0])))), [1.;3],
            TimeSignal::Table { times: vec![0.,1.,2.], values: vec![2.,5.,0.],
                dims: crate::COMPONENT_POWER_DIMS, interp: Interp::Hold }, 0.,
            MovingHeatBudget { max_vertices:8,max_elements:6,min_active_vertices:1,power_relative_tolerance:1e-10 }).unwrap();
        let changed = source.with_absorbed_power_values(&[3.,7.,0.],0.,3).unwrap();
        assert_eq!(source.absorbed_power_values(3).unwrap(),vec![2.,5.,0.]);
        assert_eq!(changed.absorbed_power_values(3).unwrap(),vec![3.,7.,0.]);
        assert_eq!(source.center_m,changed.center_m);
        assert_eq!(source.half_width_m,changed.half_width_m);
        assert_eq!(changed.absorbed_power_w.eval(2.).unwrap().value,0.);
        for bad in [vec![1.],vec![3.,-1.,0.],vec![3.,f64::NAN,0.]] {
            assert!(source.with_absorbed_power_values(&bad,0.,3).is_err());
        }
        assert!(source.with_absorbed_power_values(&[3.,7.,0.],0.,2).is_err());
    }
}
