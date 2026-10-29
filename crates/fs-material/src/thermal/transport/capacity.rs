//! One caloric kernel for source volumetric curves and assembled total capacity.
//! The owning adapters admit positive finite capacities and ordered temperatures.

use super::{ThermalTransportError, positive};

pub(super) trait CapacityKnot {
    fn temperature(&self) -> f64;
    fn capacity(&self) -> f64;
}

impl CapacityKnot for [f64; 2] {
    fn temperature(&self) -> f64 {
        self[0]
    }
    fn capacity(&self) -> f64 {
        self[1]
    }
}

/// Borrow an already admitted, positive, piecewise-affine capacity curve.
pub(super) struct CapacityCurve<'a, K> {
    knots: &'a [K],
}

impl<'a, K: CapacityKnot> CapacityCurve<'a, K> {
    pub(super) const fn new(knots: &'a [K]) -> Self {
        Self { knots }
    }

    fn span(&self) -> [f64; 2] {
        [
            self.knots[0].temperature(),
            self.knots[self.knots.len() - 1].temperature(),
        ]
    }

    fn check_temperature(&self, temperature: f64) -> Result<(), ThermalTransportError> {
        let [low, high] = self.span();
        if !temperature.is_finite() || temperature < low || temperature > high {
            return Err(invalid(
                "temperature is outside the prepared capacity interval",
            ));
        }
        Ok(())
    }

    fn segment(&self, temperature: f64) -> usize {
        self.knots
            .partition_point(|k| k.temperature() <= temperature)
            .clamp(1, self.knots.len() - 1)
            - 1
    }

    fn on_segment(&self, index: usize, temperature: f64) -> f64 {
        let a = &self.knots[index];
        let b = &self.knots[index + 1];
        if temperature == a.temperature() {
            return a.capacity();
        }
        if temperature == b.temperature() {
            return b.capacity();
        }
        let fraction = (temperature - a.temperature()) / (b.temperature() - a.temperature());
        (1.0 - fraction) * a.capacity() + fraction * b.capacity()
    }

    pub(super) fn capacity(&self, temperature: f64) -> Result<f64, ThermalTransportError> {
        self.check_temperature(temperature)?;
        Ok(self.on_segment(self.segment(temperature), temperature))
    }

    /// Integrate crossed segments directly to retain small changes near a datum.
    /// The same traversal owns energy and reversible sensible entropy.
    fn integral(
        &self,
        initial: f64,
        final_t: f64,
        entropy: bool,
    ) -> Result<f64, ThermalTransportError> {
        self.check_temperature(initial)?;
        self.check_temperature(final_t)?;
        if initial == final_t {
            return Ok(0.0);
        }
        let (low, high, sign) = if initial < final_t {
            (initial, final_t, 1.0)
        } else {
            (final_t, initial, -1.0)
        };
        let mut index = self.segment(low);
        let mut current = low;
        let mut sum = 0.0;
        let mut correction = 0.0;
        while current < high {
            let end = high.min(self.knots[index + 1].temperature());
            let a = self.on_segment(index, current);
            let b = self.on_segment(index, end);
            let amount = if entropy {
                let [cold, hot] = entropy_weights(current, end);
                if !positive(cold) || !positive(hot) {
                    return Err(invalid("sensible entropy weights are not representable"));
                }
                a * cold + b * hot
            } else {
                (0.5 * a + 0.5 * b) * (end - current)
            };
            if !positive(amount) {
                return Err(invalid(
                    "caloric increment is not positive and representable",
                ));
            }
            let corrected = amount - correction;
            let next = sum + corrected;
            if !positive(next) {
                return Err(invalid("caloric integral is not finite and representable"));
            }
            correction = (next - sum) - corrected;
            sum = next;
            current = end;
            index += 1;
        }
        Ok(sign * sum)
    }

    pub(super) fn energy_change(
        &self,
        initial: f64,
        final_t: f64,
    ) -> Result<f64, ThermalTransportError> {
        self.integral(initial, final_t, false)
    }

    pub(super) fn entropy_change(
        &self,
        initial: f64,
        final_t: f64,
    ) -> Result<f64, ThermalTransportError> {
        self.integral(initial, final_t, true)
    }

    /// Bounded segment traversal followed by a scaled quadratic inverse.
    /// No iteration, absolute-primitive subtraction, allocation or extrapolation.
    pub(super) fn temperature_after_energy_change(
        &self,
        initial: f64,
        change: f64,
    ) -> Result<f64, ThermalTransportError> {
        self.check_temperature(initial)?;
        if !change.is_finite() {
            return Err(invalid("storage change must be finite"));
        }
        if change == 0.0 {
            return Ok(initial);
        }
        let heating = change > 0.0;
        let span = self.span();
        let boundary = if heating { span[1] } else { span[0] };
        let target = change.abs();
        if let Ok(available) = self.energy_change(initial, boundary) {
            let available = available.abs();
            if target > available {
                return Err(invalid(
                    "storage target exceeds the source temperature support",
                ));
            }
            if target == available {
                return Ok(boundary);
            }
        }
        let mut current = initial;
        let mut cumulative = 0.0;
        let mut correction = 0.0;
        let mut i = if heating {
            self.segment(initial)
        } else {
            self.knots
                .partition_point(|k| k.temperature() < initial)
                .saturating_sub(1)
        };
        for _ in 0..self.knots.len() - 1 {
            let end = if heating {
                self.knots[i + 1].temperature()
            } else {
                self.knots[i].temperature()
            };
            let start_c = self.on_segment(i, current);
            let end_c = self.on_segment(i, end);
            let mean_c = 0.5 * start_c + 0.5 * end_c;
            let width = (end - current).abs();
            let amount = mean_c * width;
            if amount.is_nan() || amount <= 0.0 {
                return Err(invalid(
                    "storage inverse encountered an unrepresentable segment",
                ));
            }
            let corrected = amount - correction;
            let next = cumulative + corrected;
            if target <= next {
                if target == next {
                    return Ok(end);
                }
                let fraction_energy = if amount.is_finite() {
                    (target - cumulative) / amount
                } else {
                    ((target - cumulative) / mean_c) / width
                };
                let scale = start_c.max(end_c);
                let a = start_c / scale;
                let b = end_c / scale;
                let root = ((1.0 - fraction_energy) * a * a + fraction_energy * b * b).sqrt();
                let fraction = fraction_energy * (a + b) / (a + root);
                let temperature = current + fraction * (end - current);
                if !fraction.is_finite()
                    || !(0.0..=1.0).contains(&fraction)
                    || !temperature.is_finite()
                    || temperature == initial
                    || temperature < current.min(end)
                    || temperature > current.max(end)
                {
                    return Err(invalid(
                        "storage inverse is not representable at this scale",
                    ));
                }
                return Ok(temperature);
            }
            correction = (next - cumulative) - corrected;
            cumulative = next;
            current = end;
            if current == boundary {
                break;
            }
            if heating {
                i += 1;
            } else {
                i -= 1;
            }
        }
        Err(invalid(
            "storage inverse exceeded its bounded source traversal",
        ))
    }
}

// Positive endpoint weights for integral C(T)/T dT. A short series avoids
// cancellation in 1 - ln(1+r)/r when the temperatures are close together.
fn entropy_weights(low: f64, high: f64) -> [f64; 2] {
    let width = high - low;
    let ratio = width / low;
    let logarithm = if ratio.is_finite() {
        ratio.ln_1p()
    } else {
        high.ln() - low.ln()
    };
    let hot = if ratio < 0.125 {
        let mut sum = 0.0;
        let mut power = ratio;
        for n in 1..=18 {
            sum += power / f64::from(n + 1);
            power *= -ratio;
        }
        sum
    } else {
        1.0 - (low / width) * logarithm
    };
    [logarithm - hot, hot]
}

fn invalid(reason: &'static str) -> ThermalTransportError {
    ThermalTransportError::Invalid(reason)
}
