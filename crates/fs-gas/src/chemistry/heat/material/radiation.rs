//! Gray exchange eliminated on the same finite body's original caloric chart.

use super::*;

impl FiniteHeatBody<'_> {
    #[allow(clippy::too_many_arguments)]
    pub(in crate::chemistry::heat) fn radiative_stage_heat(
        self,
        fixed: Sum,
        volume: f64,
        capacity: f64,
        radiation: f64,
        temperature: f64,
        controls: IsochoricKineticsSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Sum, f64), IsochoricKineticsError> {
        Stage::new(self, fixed, volume, capacity, radiation, controls)?.heat(temperature, cancelled)
    }

    #[allow(clippy::too_many_arguments)]
    pub(in crate::chemistry::heat) fn radiative_stage_temperature_interval(
        self,
        fixed: Sum,
        volume: f64,
        capacity: f64,
        radiation: f64,
        bounds: [f64; 2],
        controls: IsochoricKineticsSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(f64, f64), IsochoricKineticsError> {
        let stage = Stage::new(self, fixed, volume, capacity, radiation, controls)?;
        let knots = self.curve.knots();
        let lower = stage.gas_bound(knots[0], bounds, true, cancelled)?;
        let upper = stage.gas_bound(knots[knots.len() - 1], bounds, false, cancelled)?;
        if lower > upper {
            return Err(IsochoricKineticsError::Invalid(
                "gas and radiative finite-body caloric domains do not overlap",
            ));
        }
        Ok((lower, upper))
    }
}

/// Numerical context for one already-owned chemical stage, not a time solver.
struct Stage<'m> {
    body: FiniteHeatBody<'m>,
    fixed: Sum,
    volume: f64,
    mass_density: f64,
    capacity: f64,
    radiation: f64,
    tolerance: f64,
    iterations: usize,
}

impl<'m> Stage<'m> {
    fn new(
        body: FiniteHeatBody<'m>,
        fixed: Sum,
        volume: f64,
        capacity: f64,
        radiation: f64,
        controls: IsochoricKineticsSolve,
    ) -> Result<Self, IsochoricKineticsError> {
        let mass_density = body.mass.value() / volume;
        let tolerance = controls.energy_tolerance_j / volume;
        if !positive(mass_density) || !positive(tolerance) {
            return Err(IsochoricKineticsError::Invalid(
                "unrepresentable radiative finite-body stage scale or tolerance",
            ));
        }
        Ok(Self {
            body,
            fixed,
            volume,
            mass_density,
            capacity,
            radiation,
            tolerance,
            iterations: controls.maximum_temperature_iterations,
        })
    }

    // Positive gas heat; positive feedback magnitudes with respect to Tg/Tb.
    fn exchange(
        &self,
        body_temperature: f64,
        gas_temperature: f64,
    ) -> Result<(f64, f64, f64), IsochoricKineticsError> {
        let (radiation, gas_feedback) =
            radiation_heat(self.radiation, body_temperature, gas_temperature)?;
        let heat = self.capacity * (body_temperature - gas_temperature) + radiation;
        let gas_feedback = self.capacity + gas_feedback;
        let body_feedback = self.capacity
            + ((self.radiation * body_temperature) * body_temperature) * (4.0 * body_temperature);
        if !heat.is_finite() || !positive(gas_feedback) || !positive(body_feedback) {
            return Err(IsochoricKineticsError::Invalid(
                "unrepresentable radiative finite-body heat or derivative",
            ));
        }
        Ok((heat, gas_feedback, body_feedback))
    }

    // Increasing in original h; use original authored knots to locate its one
    // source segment before any caloric inversion or nonlinear iteration.
    fn at_knot(
        &self,
        knot: EnthalpyPhaseKnot,
        gas_temperature: f64,
    ) -> Result<(f64, f64), IsochoricKineticsError> {
        let (increment, gas_feedback, _) = self.exchange(knot.temperature_k, gas_temperature)?;
        let mut residual = self.fixed;
        residual.add(
            (knot.specific_enthalpy_j_kg - self.body.state.specific_enthalpy_j_kg())
                * self.mass_density,
        );
        residual.add(increment);
        if !residual.value().is_finite() {
            return Err(IsochoricKineticsError::Invalid(
                "unrepresentable radiative finite-body source residual",
            ));
        }
        Ok((residual.value(), gas_feedback))
    }

    fn increment_at(&self, knot: EnthalpyPhaseKnot) -> f64 {
        let mut increment = Sum::default();
        increment.add(
            (self.body.state.specific_enthalpy_j_kg() - knot.specific_enthalpy_j_kg)
                * self.mass_density,
        );
        increment.add(-self.fixed.sum);
        increment.add(-self.fixed.correction);
        increment.value()
    }

    // Solve in heat relative to the fixed stage history. Then the original
    // body's h=h0-V*q/m is exactly the state used by the flux and its tangent;
    // a stiff transfer cannot magnify a separately rounded body-root residual.
    fn at_increment(
        &self,
        increment: f64,
        gas_temperature: f64,
    ) -> Result<(Sum, f64, f64, f64), IsochoricKineticsError> {
        let mut heat = self.fixed;
        heat.add(increment);
        let body = self.body.sample(heat.value() * self.volume)?;
        let (expected, gas_feedback, body_feedback) =
            self.exchange(body.temperature_k, gas_temperature)?;
        let slope = self
            .body
            .curve
            .temperature_derivative_at_specific_enthalpy(body.specific_enthalpy_j_kg)?;
        let derivative = 1.0 + (body_feedback / self.mass_density) * slope;
        let residual = increment - expected;
        let feedback = gas_feedback / derivative;
        if !positive(derivative) || !residual.is_finite() || !feedback.is_finite() || feedback < 0.0
        {
            return Err(IsochoricKineticsError::Invalid(
                "unrepresentable eliminated radiative finite-body derivative",
            ));
        }
        Ok((heat, feedback, residual, derivative))
    }

    fn heat(
        &self,
        gas_temperature: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Sum, f64), IsochoricKineticsError> {
        poll(cancelled)?;
        let knots = self.body.curve.knots();
        let (mut lower, mut upper) = (0, knots.len() - 1);
        if self.at_knot(knots[lower], gas_temperature)?.0 > 0.0
            || self.at_knot(knots[upper], gas_temperature)?.0 < 0.0
        {
            return Err(IsochoricKineticsError::Invalid(
                "radiative finite-body stage leaves its original caloric domain",
            ));
        }
        while upper - lower > 1 {
            poll(cancelled)?;
            let middle = lower + (upper - lower) / 2;
            if self.at_knot(knots[middle], gas_temperature)?.0 <= 0.0 {
                lower = middle;
            } else {
                upper = middle;
            }
        }
        if knots[lower].temperature_k == knots[upper].temperature_k {
            // Sharp latent storage has exactly constant Tb, so the complete
            // radiation/conduction increment is explicit in this one stage.
            let increment = self
                .exchange(knots[lower].temperature_k, gas_temperature)?
                .0;
            let (heat, feedback, residual, _) = self.at_increment(increment, gas_temperature)?;
            if residual.abs() <= self.tolerance {
                return Ok((heat, feedback));
            }
            return Err(self.not_converged(0, residual));
        }

        let mut lo = self.increment_at(knots[upper]);
        let mut hi = self.increment_at(knots[lower]);
        if !lo.is_finite() || !hi.is_finite() || lo >= hi {
            return Err(IsochoricKineticsError::Invalid(
                "unrepresentable radiative finite-body stage bracket",
            ));
        }
        let mut increment = if lo <= 0.0 && hi >= 0.0 {
            0.0
        } else {
            0.5 * lo + 0.5 * hi
        };
        let mut residual = 0.0;
        for iteration in 0..self.iterations {
            poll(cancelled)?;
            let (heat, feedback, value, derivative) =
                self.at_increment(increment, gas_temperature)?;
            residual = value;
            if residual.abs() <= self.tolerance {
                return Ok((heat, feedback));
            }
            if residual < 0.0 {
                lo = increment;
            } else {
                hi = increment;
            }
            let newton = increment - residual / derivative;
            let next = if newton > lo && newton < hi {
                newton
            } else {
                0.5 * lo + 0.5 * hi
            };
            if next == increment || next == lo || next == hi {
                return Err(self.not_converged(iteration + 1, residual));
            }
            increment = next;
        }
        Err(self.not_converged(self.iterations, residual))
    }

    // At each source endpoint, R(Tg) decreases strictly. Select the admitted
    // side of its root so the outer gas caloric solve never probes an
    // out-of-source body. Most stages need no additional root at either end.
    fn gas_bound(
        &self,
        knot: EnthalpyPhaseKnot,
        bounds: [f64; 2],
        lower_bound: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, IsochoricKineticsError> {
        poll(cancelled)?;
        let [mut lo, mut hi] = bounds;
        let at_lo = self.at_knot(knot, lo)?.0;
        let at_hi = self.at_knot(knot, hi)?.0;
        if lower_bound && at_lo <= 0.0 {
            return Ok(lo);
        }
        if !lower_bound && at_hi >= 0.0 {
            return Ok(hi);
        }
        if at_lo < 0.0 || at_hi > 0.0 {
            return Err(IsochoricKineticsError::Invalid(
                "gas and radiative finite-body caloric domains do not overlap",
            ));
        }
        let mut temperature = 0.5 * lo + 0.5 * hi;
        let mut residual = 0.0;
        for iteration in 0..self.iterations {
            poll(cancelled)?;
            let (value, feedback) = self.at_knot(knot, temperature)?;
            residual = value;
            let admitted = if lower_bound {
                residual <= 0.0
            } else {
                residual >= 0.0
            };
            if admitted && residual.abs() <= self.tolerance {
                return Ok(temperature);
            }
            if residual > 0.0 {
                lo = temperature;
            } else {
                hi = temperature;
            }
            let newton = temperature + residual / feedback;
            let next = if newton > lo && newton < hi && newton != temperature {
                newton
            } else {
                0.5 * lo + 0.5 * hi
            };
            if next == lo || next == hi {
                let inside = if lower_bound { hi } else { lo };
                let residual = self.at_knot(knot, inside)?.0;
                return if residual.abs() <= self.tolerance {
                    Ok(inside)
                } else {
                    Err(self.not_converged(iteration + 1, residual))
                };
            }
            temperature = next;
        }
        Err(self.not_converged(self.iterations, residual))
    }

    fn not_converged(&self, iterations: usize, residual: f64) -> IsochoricKineticsError {
        IsochoricKineticsError::TemperatureNoConvergence {
            iterations,
            residual_j: residual * self.volume,
        }
    }
}
