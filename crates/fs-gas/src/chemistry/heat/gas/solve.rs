//! Original-mixture caloric elimination inside one chemical stage.

use super::*;

impl FiniteGasHeatBody<'_> {
    fn exchange(
        self,
        body_t: f64,
        gas_t: f64,
        capacity: f64,
        radiation: f64,
    ) -> Result<(f64, f64, f64), IsochoricKineticsError> {
        let (heat, gas_feedback, body_feedback) = if radiation == 0.0 {
            (0.0, 0.0, 0.0)
        } else {
            let (heat, gas_feedback) = radiation_heat(radiation, body_t, gas_t)?;
            (
                heat,
                gas_feedback,
                ((radiation * body_t) * body_t) * (4.0 * body_t),
            )
        };
        let result = (
            capacity * (body_t - gas_t) + heat,
            capacity + gas_feedback,
            capacity + body_feedback,
        );
        if [result.0, result.1, result.2]
            .iter()
            .any(|v| !v.is_finite())
        {
            return Err(IsochoricKineticsError::Invalid(
                "unrepresentable finite gas exchange",
            ));
        }
        Ok(result)
    }

    #[allow(clippy::too_many_arguments)]
    pub(in crate::chemistry::heat) fn stage_heat(
        self,
        fixed: Sum,
        volume: f64,
        capacity: f64,
        radiation: f64,
        gas_t: f64,
        controls: IsochoricKineticsSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Sum, f64), IsochoricKineticsError> {
        if capacity == 0.0 && radiation == 0.0 {
            self.temperature_at_heat(fixed.value() * volume, controls, cancelled)?;
            return Ok((fixed, 0.0));
        }
        let evaluate = |body_t| {
            let (energy, body_capacity) = self.caloric(Temperature::new(body_t))?;
            let (exchange, _, body_feedback) = self.exchange(body_t, gas_t, capacity, radiation)?;
            let mut residual = fixed;
            residual.add(exchange);
            residual.add(energy / volume);
            Ok((
                residual.value() * volume,
                body_capacity + volume * body_feedback,
            ))
        };
        let body_t = temperature(
            self.bounds.map(Temperature::value),
            self.temperature().value(),
            controls,
            Side::Either,
            evaluate,
            cancelled,
        )?;
        let (energy, body_capacity) = self.caloric(Temperature::new(body_t))?;
        let (exchange, gas_feedback, body_feedback) =
            self.exchange(body_t, gas_t, capacity, radiation)?;
        let mut heat = fixed;
        heat.add(exchange);
        // Correct only the admitted root residual while retaining compensated
        // history. Then q and the body's original energy increment agree;
        // strong coupling cannot magnify an independent body-state roundoff.
        let mut residual = heat;
        residual.add(energy / volume);
        heat.add(-residual.sum);
        heat.add(-residual.correction);
        let feedback = gas_feedback / (1.0 + (volume / body_capacity) * body_feedback);
        if !heat.value().is_finite() || !feedback.is_finite() || feedback < 0.0 {
            return Err(IsochoricKineticsError::Invalid(
                "unrepresentable eliminated finite gas heat or derivative",
            ));
        }
        Ok((heat, feedback))
    }

    #[allow(clippy::too_many_arguments)]
    pub(in crate::chemistry::heat) fn stage_temperature_interval(
        self,
        fixed: Sum,
        volume: f64,
        capacity: f64,
        radiation: f64,
        bounds: [f64; 2],
        controls: IsochoricKineticsSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(f64, f64), IsochoricKineticsError> {
        if capacity == 0.0 && radiation == 0.0 {
            self.temperature_at_heat(fixed.value() * volume, controls, cancelled)?;
            return Ok((bounds[0], bounds[1]));
        }
        let mut result = bounds;
        for (index, body_t) in self.bounds.into_iter().enumerate() {
            poll(cancelled)?;
            let energy = self.caloric(body_t)?.0;
            let evaluate = |gas_t| {
                let (exchange, gas_feedback, _) =
                    self.exchange(body_t.value(), gas_t, capacity, radiation)?;
                let mut residual = fixed;
                residual.add(exchange);
                residual.add(energy / volume);
                Ok((-residual.value() * volume, gas_feedback * volume))
            };
            let at_lo = evaluate(bounds[0])?.0;
            let at_hi = evaluate(bounds[1])?.0;
            if (index == 0 && at_lo >= 0.0) || (index == 1 && at_hi <= 0.0) {
                continue;
            }
            let side = if index == 0 { Side::Upper } else { Side::Lower };
            result[index] = temperature(
                bounds,
                0.5 * bounds[0] + 0.5 * bounds[1],
                controls,
                side,
                evaluate,
                cancelled,
            )?;
        }
        if result[0] > result[1] {
            return Err(IsochoricKineticsError::Invalid(
                "reactor and finite gas source domains do not overlap",
            ));
        }
        Ok((result[0], result[1]))
    }
}

pub(super) enum Side {
    Either,
    Lower,
    Upper,
}

/// Bounded monotone caloric root: residuals are J and derivatives are J/K.
/// The side restriction preserves source-domain admission at clipping roots.
pub(super) fn temperature(
    bounds: [f64; 2],
    hint: f64,
    controls: IsochoricKineticsSolve,
    side: Side,
    mut evaluate: impl FnMut(f64) -> Result<(f64, f64), IsochoricKineticsError>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, IsochoricKineticsError> {
    poll(cancelled)?;
    let [mut lo, mut hi] = bounds;
    // Absolute reactor Joules alone cannot resolve a tiny partner inventory.
    // Reserve one eighth of the existing relative-temperature error allowance
    // for this inverse, so the outer full/half comparison remains meaningful.
    let admitted = |t: f64, r: f64, derivative: f64| {
        positive(derivative)
            && r.abs() <= controls.energy_tolerance_j
            && (r / derivative).abs() <= 0.125 * controls.relative_tolerance * t
            && match side {
                Side::Either => true,
                Side::Lower => r <= 0.0,
                Side::Upper => r >= 0.0,
            }
    };
    let (at_lo, derivative_lo) = evaluate(lo)?;
    let (at_hi, derivative_hi) = evaluate(hi)?;
    if admitted(lo, at_lo, derivative_lo) {
        return Ok(lo);
    }
    if admitted(hi, at_hi, derivative_hi) {
        return Ok(hi);
    }
    if !at_lo.is_finite() || !at_hi.is_finite() || at_lo > 0.0 || at_hi < 0.0 {
        return Err(IsochoricKineticsError::Invalid(
            "finite gas heat leaves the admitted source domain",
        ));
    }
    let mut t = hint.clamp(lo, hi);
    let mut residual = 0.0;
    for iteration in 0..controls.maximum_temperature_iterations {
        poll(cancelled)?;
        let (value, derivative) = evaluate(t)?;
        residual = value;
        if !residual.is_finite() || !positive(derivative) {
            return Err(IsochoricKineticsError::Invalid(
                "unrepresentable finite gas caloric residual or tangent",
            ));
        }
        if admitted(t, residual, derivative) {
            return Ok(t);
        }
        if residual < 0.0 {
            lo = t;
        } else {
            hi = t;
        }
        let newton = t - residual / derivative;
        let next = if newton > lo && newton < hi && newton != t {
            newton
        } else {
            0.5 * lo + 0.5 * hi
        };
        if next == lo || next == hi {
            for end in [lo, hi] {
                let (r, derivative) = evaluate(end)?;
                if admitted(end, r, derivative) {
                    return Ok(end);
                }
            }
            return Err(IsochoricKineticsError::TemperatureNoConvergence {
                iterations: iteration + 1,
                residual_j: residual,
            });
        }
        t = next;
    }
    Err(IsochoricKineticsError::TemperatureNoConvergence {
        iterations: controls.maximum_temperature_iterations,
        residual_j: residual,
    })
}
