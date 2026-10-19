//! Immutable path and footprint edits for the original normalized heat source.
//!
//! Center coordinates are authored SI values, not sample positions at an
//! optimizer-selected time. Times, interpolation, power and mesh admission stay
//! unchanged. Polynomial paths are retained when untouched; editing their
//! coefficients requires a separate construction API and is explicitly refused.

use super::{MovingHeatSource, ConductionError, invalid, poll};
use fs_exec::Cx;
use fs_scenario::signal::TimeSignal;
use std::collections::BTreeSet;

/// One metre-valued design coordinate in a moving source declaration.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum MovingHeatGeometryControl {
    /// World axis 0..3 and a constant/ramp/table value index. Ramp order is
    /// [from, to]; a table retains its original value and knot order.
    CenterValue { axis: usize, value: usize },
    /// World-axis ellipsoid half width, strictly positive metres.
    HalfWidth { axis: usize },
}

impl MovingHeatSource {
    /// Read one original path value or half width without allocating a signal.
    /// Negative or zero center coordinates are meaningful and remain admitted.
    pub fn geometry_value(&self, control: MovingHeatGeometryControl)
        -> Result<f64, ConductionError>
    {
        match control {
            MovingHeatGeometryControl::HalfWidth { axis } => self.half_width_m.get(axis).copied()
                .ok_or_else(|| invalid("heat-source width axis must be 0, 1 or 2")),
            MovingHeatGeometryControl::CenterValue { axis, value } => {
                let signal = self.center_m.get(axis)
                    .ok_or_else(|| invalid("heat-source center axis must be 0, 1 or 2"))?;
                let result = match signal {
                    TimeSignal::Constant(q) => (value == 0).then_some(q.value),
                    TimeSignal::Ramp { from, to, .. } => [from.value, to.value].get(value).copied(),
                    TimeSignal::Table { values, .. } => values.get(value).copied(),
                    TimeSignal::Chebfun(_) => return Err(invalid(
                        "center-value editing supports constant/ramp/table paths, not polynomial coefficients")),
                };
                result.ok_or_else(|| invalid("heat-source center control names an absent authored value"))
            }
        }
    }

    /// Return a new source with selected metre-valued edits. Validate the entire
    /// batch before cloning; a repeated coordinate is an error, not last-wins.
    /// Untouched path values, every time knot, power and sampling budget survive.
    /// The cap covers the source's complete authored scalar controls, including
    /// untouched signals, and separately the edit count. It is not a byte cap.
    ///
    /// The source and any accepted trajectory are never mutated. The original
    /// constructor owns quantity, path, width and reference-time validation;
    /// later mesh support and power checks still occur at actual primal samples.
    pub fn with_geometry_values(&self, cx: &Cx<'_>, edits: &[(MovingHeatGeometryControl, f64)],
        reference_time_s: f64, max_controls: usize) -> Result<Self, ConductionError>
    {
        poll(cx, 0)?;
        if self.schedule_control_count()? > max_controls || edits.len() > max_controls {
            return Err(invalid("heat-source geometry exceeds its authored-control allowance"));
        }
        let mut seen = BTreeSet::new();
        for (index, &(control, value)) in edits.iter().enumerate() {
            poll(cx, index)?;
            self.geometry_value(control)?;
            if !seen.insert(control) || !value.is_finite()
                || matches!(control, MovingHeatGeometryControl::HalfWidth { .. }) && value <= 0.0
            { return Err(invalid("geometry edits must be unique, finite and strictly positive for half widths")); }
        }
        let mut center = self.center_m.clone();
        let mut width = self.half_width_m;
        for (index, &(control, replacement)) in edits.iter().enumerate() {
            poll(cx, index)?;
            match control {
                MovingHeatGeometryControl::HalfWidth { axis } => width[axis] = replacement,
                MovingHeatGeometryControl::CenterValue { axis, value } => match &mut center[axis] {
                    TimeSignal::Constant(q) => q.value = replacement,
                    TimeSignal::Ramp { from, to, .. } => {
                        if value == 0 { from.value = replacement; } else { to.value = replacement; }
                    }
                    TimeSignal::Table { values, .. } => values[value] = replacement,
                    TimeSignal::Chebfun(_) => unreachable!("polynomial edits refused before cloning"),
                },
            }
        }
        let source = Self::new(center, width, self.absorbed_power_w.clone(), reference_time_s, self.budget)?;
        poll(cx, edits.len())?;
        Ok(source)
    }
}
