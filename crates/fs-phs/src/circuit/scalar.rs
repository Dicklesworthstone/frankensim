//! Shared bounded scalar storage arithmetic for charge and flux laws.
//!
//! The conserved coordinate and its work-conjugate effort are unit-neutral
//! here. Physical facades supply the capacitor/inductor units and diagnostics.
//! One immutable piecewise-affine law owns its inverse, integral, exact work
//! chord and analytic chord derivative; none is fitted independently.

use std::sync::Arc;

pub(crate) const MAX_SCALAR_KNOTS: usize = 4096;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum CurveError {
    Knots,
    ZeroKnot,
    Polarities,
    PositiveSlope,
    Arithmetic,
    Energy,
    StoredEnergy,
    CoordinateSupport,
    EffortSupport,
    Effort,
    Coordinate,
}

fn finite(value: f64) -> Result<f64, CurveError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(CurveError::Arithmetic)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct ScalarCurve {
    knots: Arc<[(f64, f64)]>,
    slopes: Arc<[f64]>,
    energy_at_knots: Arc<[f64]>,
}

impl ScalarCurve {
    pub(crate) fn try_new(knots: Vec<(f64, f64)>) -> Result<Self, CurveError> {
        if !(3..=MAX_SCALAR_KNOTS).contains(&knots.len())
            || knots.iter().any(|(x, e)| !x.is_finite() || !e.is_finite())
        {
            return Err(CurveError::Knots);
        }
        let zero = knots
            .iter()
            .position(|&(x, e)| x == 0.0 && e == 0.0)
            .ok_or(CurveError::ZeroKnot)?;
        if zero == 0 || zero + 1 == knots.len() {
            return Err(CurveError::Polarities);
        }
        let mut slopes = Vec::with_capacity(knots.len() - 1);
        for pair in knots.windows(2) {
            let dx = finite(pair[1].0 - pair[0].0)?;
            let de = finite(pair[1].1 - pair[0].1)?;
            let slope = de / dx;
            if dx <= 0.0
                || de <= 0.0
                || !slope.is_finite()
                || slope <= 0.0
                || !slope.recip().is_finite()
                || slope.recip() <= 0.0
            {
                return Err(CurveError::PositiveSlope);
            }
            slopes.push(slope);
        }
        let mut energy_at_knots = vec![0.0; knots.len()];
        for i in zero..knots.len() - 1 {
            let area = f64::midpoint(knots[i].1, knots[i + 1].1) * (knots[i + 1].0 - knots[i].0);
            energy_at_knots[i + 1] = finite(energy_at_knots[i] + area)?;
            if area <= 0.0 || energy_at_knots[i + 1] <= 0.0 {
                return Err(CurveError::Energy);
            }
        }
        for i in (0..zero).rev() {
            let area = f64::midpoint(knots[i].1, knots[i + 1].1) * (knots[i + 1].0 - knots[i].0);
            energy_at_knots[i] = finite(energy_at_knots[i + 1] - area)?;
            if area >= 0.0 || energy_at_knots[i] <= 0.0 {
                return Err(CurveError::Energy);
            }
        }
        Ok(Self {
            knots: knots.into(),
            slopes: slopes.into(),
            energy_at_knots: energy_at_knots.into(),
        })
    }

    pub(crate) fn knots(&self) -> &[(f64, f64)] {
        &self.knots
    }

    pub(crate) fn range(&self) -> (f64, f64) {
        (self.knots[0].0, self.knots[self.knots.len() - 1].0)
    }

    fn segment(&self, coordinate: f64) -> Result<usize, CurveError> {
        let (lo, hi) = self.range();
        if !coordinate.is_finite() || coordinate < lo || coordinate > hi {
            return Err(CurveError::CoordinateSupport);
        }
        Ok(self
            .knots
            .partition_point(|&(x, _)| x <= coordinate)
            .saturating_sub(1)
            .min(self.knots.len() - 2))
    }

    pub(crate) fn effort(&self, coordinate: f64) -> Result<f64, CurveError> {
        let i = self.segment(coordinate)?;
        if coordinate == self.knots[i + 1].0 {
            return Ok(self.knots[i + 1].1);
        }
        let anchor = if coordinate < 0.0 { i + 1 } else { i };
        let value = finite(
            self.slopes[i].mul_add(coordinate - self.knots[anchor].0, self.knots[anchor].1),
        )?;
        if coordinate != 0.0 && value == 0.0 {
            return Err(CurveError::Effort);
        }
        Ok(value)
    }

    pub(crate) fn coordinate(&self, effort: f64) -> Result<f64, CurveError> {
        if !effort.is_finite()
            || effort < self.knots[0].1
            || effort > self.knots[self.knots.len() - 1].1
        {
            return Err(CurveError::EffortSupport);
        }
        let i = self
            .knots
            .partition_point(|&(_, e)| e <= effort)
            .saturating_sub(1)
            .min(self.knots.len() - 2);
        if effort == self.knots[i + 1].1 {
            return Ok(self.knots[i + 1].0);
        }
        let anchor = if effort < 0.0 { i + 1 } else { i };
        let coordinate =
            finite((effort - self.knots[anchor].1) / self.slopes[i] + self.knots[anchor].0)?;
        if effort != 0.0 && coordinate == 0.0 {
            return Err(CurveError::Coordinate);
        }
        self.segment(coordinate)?;
        Ok(coordinate)
    }

    pub(crate) fn inverse_slope(&self, coordinate: f64) -> Result<f64, CurveError> {
        Ok(self.slopes[self.segment(coordinate)?].recip())
    }

    pub(crate) fn energy(&self, coordinate: f64) -> Result<f64, CurveError> {
        let i = self.segment(coordinate)?;
        // Anchor nearer zero so a tiny negative-coordinate energy is not the
        // subtraction of almost all the energy at a large negative knot.
        let anchor = if coordinate < 0.0 { i + 1 } else { i };
        let dx = coordinate - self.knots[anchor].0;
        let effort = self.effort(coordinate)?;
        let energy = finite(
            self.energy_at_knots[anchor] + dx * f64::midpoint(self.knots[anchor].1, effort),
        )?;
        if energy < 0.0 || coordinate != 0.0 && energy == 0.0 {
            return Err(CurveError::StoredEnergy);
        }
        Ok(energy)
    }

    // Integral_base^end [e(x)-e(base)] dx. Walking outward from base adds
    // positive areas without subtracting large biased endpoint energies.
    fn relative_energy(&self, base: f64, end: f64) -> Result<f64, CurveError> {
        let mut i = self.segment(base)?;
        self.segment(end)?;
        let mut position = base;
        let mut effort_change = 0.0;
        let mut area = 0.0;
        while position != end {
            if end < position && position == self.knots[i].0 {
                i -= 1;
            }
            let next = if end > position {
                end.min(self.knots[i + 1].0)
            } else {
                end.max(self.knots[i].0)
            };
            let dx = next - position;
            let de = finite(self.slopes[i] * dx)?;
            area = finite(area + dx * (effort_change + 0.5 * de))?;
            effort_change = finite(effort_change + de)?;
            position = next;
            if position != end && end > position {
                i += 1;
            }
        }
        Ok(area)
    }

    /// Integral_before^after [e(x)-e(midpoint)] dx. The reference may be any
    /// admitted coordinate, including a represented midpoint with rounding.
    pub(crate) fn midpoint_remainder(
        &self,
        before: f64,
        after: f64,
        midpoint: f64,
    ) -> Result<f64, CurveError> {
        let a = self.segment(before)?;
        let b = self.segment(after)?;
        let m = self.segment(midpoint)?;
        if a == b && a == m {
            return finite(
                (after - before)
                    * (0.5 * self.slopes[m])
                    * ((before - midpoint) + (after - midpoint)),
            );
        }
        finite(self.relative_energy(midpoint, after)? - self.relative_energy(midpoint, before)?)
    }

    pub(crate) fn energy_change(&self, before: f64, after: f64) -> Result<f64, CurveError> {
        let midpoint = f64::midpoint(before, after);
        finite(
            self.effort(midpoint)? * (after - before)
                + self.midpoint_remainder(before, after, midpoint)?,
        )
    }

    pub(crate) fn discrete_effort(&self, before: f64, after: f64) -> Result<f64, CurveError> {
        let a = self.segment(before)?;
        let b = self.segment(after)?;
        let mid = f64::midpoint(before, after);
        let effort = self.effort(mid)?;
        if before == after {
            return Ok(effort);
        }
        let correction = if a == b {
            (0.5 * self.slopes[a]) * ((before - mid) + (after - mid))
        } else {
            self.midpoint_remainder(before, after, mid)? / (after - before)
        };
        finite(effort + correction)
    }

    pub(crate) fn discrete_effort_end_derivative(
        &self,
        before: f64,
        after: f64,
    ) -> Result<f64, CurveError> {
        let a = self.segment(before)?;
        let b = self.segment(after)?;
        if before == after || a == b {
            return finite(0.5 * self.slopes[b]);
        }
        let change = after - before;
        // Integral_a^b [e(b)-e(x)] dx / (b-a)^2, relative to e(b): do not
        // subtract nearly equal efforts or square a tiny coordinate change.
        finite((-self.midpoint_remainder(before, after, after)? / change) / change)
    }
}
