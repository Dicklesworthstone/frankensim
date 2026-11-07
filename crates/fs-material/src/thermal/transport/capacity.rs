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

/// Heat and exchange mean retain the solved relative travel, before absolute
/// temperatures are rounded. Capacity units determine the energy units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct CapacityExchange {
    pub temperature: f64,
    pub energy_change: f64,
    pub mean_temperature: f64,
    pub inlet_energy_gain: f64,
    pub inlet_energy_loss: f64,
    pub reservoir_energy_derivative: f64,
    pub mean_reservoir_derivative: f64,
}

/// Prepared search coordinates for repeated caloric queries from one fixed datum.
///
/// The coordinates locate an original capacity segment; they do not replace
/// C(T) by a piecewise-linear energy law. Each query still solves the original
/// quadratic integral on that segment. Datum-aligned forward queries add the
/// original local integral to a same-sign prefix. Other relative heat and entropy
/// integrations remain on `CapacityCurve`, without cached-primitive subtraction.
#[derive(Clone, Debug, PartialEq)]
pub(super) struct CapacityInverse {
    datum: f64,
    energy: Vec<f64>,
}

impl CapacityInverse {
    /// O(knots) preparation, with independent compensated sums on both sides
    /// of the datum. An unrepresentable index is only an optimization refusal:
    /// the owning law can retain the original direct inversion in that case.
    pub(super) fn new<K: CapacityKnot>(
        curve: &CapacityCurve<'_, K>,
        datum: f64,
    ) -> Result<Self, ThermalTransportError> {
        curve.check_temperature(datum)?;
        let split = curve.segment(datum);
        let mut energy = vec![0.0; curve.knots.len()];
        for heating in [false, true] {
            let mut current = datum;
            let mut sum = 0.0;
            let mut correction = 0.0;
            let count = if heating {
                curve.knots.len() - split - 1
            } else {
                split + 1
            };
            for step in 0..count {
                let index = if heating {
                    split + step + 1
                } else {
                    split - step
                };
                let segment = if heating { index - 1 } else { index };
                let end = curve.knots[index].temperature();
                if end != current {
                    let a = curve.on_segment(segment, current);
                    let b = curve.knots[index].capacity();
                    let amount = (0.5 * a + 0.5 * b) * (end - current).abs();
                    if !positive(amount) {
                        return Err(invalid("caloric search increment is not representable"));
                    }
                    let corrected = amount - correction;
                    let next = sum + corrected;
                    if !positive(next) {
                        return Err(invalid("caloric search coordinate is not representable"));
                    }
                    correction = (next - sum) - corrected;
                    sum = next;
                }
                energy[index] = if heating { sum } else { -sum };
                current = end;
            }
        }
        if energy.windows(2).any(|pair| pair[0] >= pair[1]) {
            return Err(invalid("caloric search knots are not numerically distinct"));
        }
        Ok(Self { datum, energy })
    }

    /// Forward energy from the same datum, with logarithmic knot lookup and
    /// constant source work away from the domain endpoints. The datum's own
    /// segment uses the direct integral, preserving small sensible increments.
    /// Elsewhere the anchor lies BETWEEN the datum and the target, so the
    /// prefix and local energy have the same sign and cannot cancel.
    ///
    /// The canonical endpoint integrals retain their original accumulation
    /// order. Owners that already retain those two values may return them
    /// directly. No difference between two large cached energies is taken.
    pub(super) fn energy_change_from_datum<K: CapacityKnot>(
        &self,
        curve: &CapacityCurve<'_, K>,
        temperature: f64,
    ) -> Result<f64, ThermalTransportError> {
        curve.check_temperature(temperature)?;
        let [cold, hot] = curve.span();
        let segment = curve.segment(temperature);
        if temperature == cold || temperature == hot || segment == curve.segment(self.datum) {
            return curve.energy_change(self.datum, temperature);
        }
        let anchor = if temperature > self.datum {
            segment
        } else {
            segment + 1
        };
        let local = curve.energy_change(curve.knots[anchor].temperature(), temperature)?;
        let energy = self.energy[anchor] + local;
        if !energy.is_finite() || energy == 0.0 || (energy > 0.0) != (temperature > self.datum) {
            // The index is only an optimization. Its accumulation order must
            // not refuse a source integral that the direct owner can represent.
            return curve.energy_change(self.datum, temperature);
        }
        Ok(energy)
    }

    /// Binary search plus one original analytic segment inverse, with no
    /// allocation or walk from the datum. Near a rounded knot or support bound,
    /// retain the direct owner's admission instead of clipping heat or T.
    pub(super) fn temperature<K: CapacityKnot>(
        &self,
        curve: &CapacityCurve<'_, K>,
        energy: f64,
    ) -> Result<f64, ThermalTransportError> {
        if energy == 0.0 {
            return Ok(self.datum);
        }
        let direct = || curve.temperature_after_energy_change(self.datum, energy);
        if !energy.is_finite()
            || energy <= self.energy[0]
            || energy >= self.energy[self.energy.len() - 1]
        {
            return direct();
        }
        let lower = self.energy.partition_point(|value| *value <= energy) - 1;
        let upper = lower + 1;
        if self.energy[lower] <= 0.0 && self.energy[upper] >= 0.0 {
            // The datum's own segment is already constant-work. Subtracting
            // a distant knot here can turn an unrepresentable tiny heat into
            // a one-ulp temperature change, even of the wrong sign.
            return direct();
        }
        if energy == self.energy[lower] {
            return Ok(curve.knots[lower].temperature());
        }
        // Work from the nearer endpoint to retain a small local increment.
        let anchor = if energy - self.energy[lower] <= self.energy[upper] - energy {
            lower
        } else {
            upper
        };
        let other = if anchor == lower { upper } else { lower };
        let start = curve.knots[anchor].temperature();
        let end = curve.knots[other].temperature();
        let change = energy - self.energy[anchor];
        let available = (0.5 * curve.knots[anchor].capacity()
            + 0.5 * curve.knots[other].capacity())
            * (end - start);
        // Prefix sums can round by a final ulp at source knots. Such an edge
        // must not make a supported query cross a different source interval.
        if change.abs() <= available.abs()
            && let Ok(temperature) = curve.temperature_after_energy_change(start, change)
            && temperature >= start.min(end)
            && temperature <= start.max(end)
            && temperature != self.datum
        {
            return Ok(temperature);
        }
        direct()
    }
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

    /// Integrate C(T) dT/(reservoir-T) for passive linear heat exchange.
    /// Only crossed source segments are visited; the final affine segment is
    /// inverted without freezing C or introducing another caloric chart.
    pub(super) fn temperature_after_linear_exchange(
        &self,
        initial: f64,
        reservoir: f64,
        exposure: f64,
    ) -> Result<f64, ThermalTransportError> {
        self.linear_heat_exchange(initial, reservoir, exposure)
            .map(|answer| answer.temperature)
    }

    pub(super) fn linear_heat_exchange(
        &self,
        initial: f64,
        reservoir: f64,
        exposure: f64,
    ) -> Result<CapacityExchange, ThermalTransportError> {
        self.check_temperature(initial)?;
        if !positive(reservoir) || !exposure.is_finite() || exposure < 0.0 {
            return Err(invalid(
                "linear heat exchange needs positive reservoir temperature and finite nonnegative exposure",
            ));
        }
        if exposure == 0.0 || initial == reservoir {
            let capacity = self.capacity(initial)?;
            let attenuation = exposure / capacity;
            let loss = -(-attenuation).exp_m1();
            let response = if attenuation == 0.0 {
                exposure
            } else if attenuation <= 0.5 {
                exposure * (loss / attenuation)
            } else {
                capacity * loss
            };
            let mean_response = if exposure == 0.0 {
                0.0
            } else if loss <= 0.5 {
                let [a0, a1, _, _] = exchange_moment_weights(loss);
                loss * a1 / (a0 + a1)
            } else {
                1.0 - response / exposure
            };
            return Ok(CapacityExchange {
                temperature: initial,
                energy_change: 0.0,
                mean_temperature: initial,
                inlet_energy_gain: (-attenuation).exp(),
                inlet_energy_loss: loss,
                reservoir_energy_derivative: response,
                mean_reservoir_derivative: mean_response,
            });
        }
        let heating = reservoir > initial;
        let direction = if heating { 1.0 } else { -1.0 };
        let mut index = if heating {
            self.segment(initial)
        } else {
            self.knots
                .partition_point(|k| k.temperature() < initial)
                .saturating_sub(1)
        };
        let mut current = initial;
        let mut cumulative = 0.0;
        let mut correction = 0.0;
        let (mut energy, mut energy_correction) = (0.0, 0.0);
        let (mut mean_travel, mut mean_correction) = (0.0, 0.0);
        let (mut inlet_energy_gain, mut inlet_energy_loss) = (1.0, 0.0);
        let mut reservoir_energy_derivative = 0.0;
        let (mut mean_response_numerator, mut response_correction) = (0.0, 0.0);
        for _ in 0..self.knots.len() - 1 {
            let end = if heating {
                self.knots[index + 1].temperature().min(reservoir)
            } else {
                self.knots[index].temperature().max(reservoir)
            };
            if end == current {
                break;
            }
            let segment = ExchangeSegment {
                distance: (reservoir - current).abs(),
                width: (end - current).abs(),
                initial_capacity: self.on_segment(index, current),
                end_capacity: self.on_segment(index, end),
            };
            let amount = segment.exposure(segment.width)?;
            let corrected = amount - correction;
            let next = cumulative + corrected;
            let finished = exposure <= next;
            let local_exposure = if finished {
                exposure - cumulative
            } else {
                amount
            };
            let travel = if !finished || exposure == next {
                segment.width
            } else {
                segment.inverse(local_exposure)?
            };
            let heat = (0.5 * segment.initial_capacity + 0.5 * segment.capacity(travel)) * travel;
            if !heat.is_finite() || !positive(local_exposure) {
                return Err(invalid(
                    "heat-exchange energy or exposure is not representable",
                ));
            }
            let gain = segment.inlet_energy_gain(travel, local_exposure)?;
            let loss = if segment.initial_capacity == segment.end_capacity {
                -(-local_exposure / segment.initial_capacity).exp_m1()
            } else {
                travel / segment.distance
            };
            let (local_mean, response, mean_response) =
                segment.exchange_moments(travel, local_exposure, heat, gain);
            // X-B stays a positive integral even when B rounds to X. This
            // preserves a weak wall response of the stream's Robin mean.
            compensated_add(
                &mut mean_response_numerator,
                &mut response_correction,
                loss * reservoir_energy_derivative + local_exposure * mean_response,
            );
            reservoir_energy_derivative = gain.mul_add(reservoir_energy_derivative, response);
            inlet_energy_loss = gain.mul_add(inlet_energy_loss, loss);
            inlet_energy_gain *= gain;
            if !reservoir_energy_derivative.is_finite()
                || !mean_response_numerator.is_finite()
                || reservoir_energy_derivative < 0.0
                || mean_response_numerator < 0.0
            {
                return Err(invalid("heat-exchange wall response is not representable"));
            }
            compensated_add(&mut energy, &mut energy_correction, heat);
            compensated_add(
                &mut mean_travel,
                &mut mean_correction,
                (local_exposure / exposure) * ((current - initial).abs() + local_mean),
            );
            if finished {
                let temperature = if travel == segment.width {
                    end
                } else if travel <= 0.5 * segment.distance {
                    current + direction * travel
                } else {
                    reservoir - direction * (segment.distance - travel)
                };
                if temperature == initial {
                    return Err(invalid(
                        "nonzero heat exchange is below temperature resolution",
                    ));
                }
                let mean_temperature = initial + direction * mean_travel;
                if !positive(energy)
                    || !positive(temperature)
                    || !positive(mean_temperature)
                    || mean_temperature < initial.min(temperature)
                    || mean_temperature > initial.max(temperature)
                    || temperature < current.min(end)
                    || temperature > current.max(end)
                {
                    return Err(invalid("heat-exchange state or mean is not representable"));
                }
                return Ok(CapacityExchange {
                    temperature,
                    energy_change: direction * energy,
                    mean_temperature,
                    inlet_energy_gain,
                    inlet_energy_loss,
                    reservoir_energy_derivative,
                    mean_reservoir_derivative: mean_response_numerator / exposure,
                });
            }
            correction = (next - cumulative) - corrected;
            cumulative = next;
            current = end;
            if heating {
                if index + 2 == self.knots.len() {
                    break;
                }
                index += 1;
            } else {
                if index == 0 {
                    break;
                }
                index -= 1;
            }
        }
        Err(invalid(
            "heat exchange exceeds the source temperature support",
        ))
    }

    #[cfg(test)]
    fn exchange_integral(
        &self,
        index: usize,
        initial: f64,
        final_t: f64,
        reservoir: f64,
    ) -> Result<f64, ThermalTransportError> {
        let segment = ExchangeSegment {
            distance: (reservoir - initial).abs(),
            width: (final_t - initial).abs(),
            initial_capacity: self.on_segment(index, initial),
            end_capacity: self.on_segment(index, final_t),
        };
        segment.exposure(segment.width)
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
            // A step resolved on its initial segment already has bounded
            // source support. Do not integrate the entire remaining table on
            // every local heat update: that made even an infinitesimal step
            // linear in the number of unrelated source knots. Only a crossing
            // needs the whole-span check, which also preserves the canonical
            // endpoint integral's rounding for reverse (cooling) traversal.
            if current == initial {
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
    entropy_weights_with_width(low, high, width)
}

fn entropy_weights_with_width(low: f64, high: f64, width: f64) -> [f64; 2] {
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

// Positive travel toward the reservoir, with the original affine capacity.
// Solving in travel preserves small heat that an absolute-T inverse quantizes.
struct ExchangeSegment {
    distance: f64,
    width: f64,
    initial_capacity: f64,
    end_capacity: f64,
}

impl ExchangeSegment {
    fn capacity(&self, travel: f64) -> f64 {
        let fraction = travel / self.width;
        (1.0 - fraction) * self.initial_capacity + fraction * self.end_capacity
    }

    fn exposure(&self, travel: f64) -> Result<f64, ThermalTransportError> {
        if travel == 0.0 {
            return Ok(0.0);
        }
        if travel == self.distance {
            return Ok(f64::INFINITY);
        }
        let [near, far] = entropy_weights_with_width(self.distance - travel, self.distance, travel);
        let amount = self.capacity(travel) * near + self.initial_capacity * far;
        if amount.is_nan() || amount < 0.0 {
            return Err(invalid("heat-exchange integral is not representable"));
        }
        Ok(amount)
    }

    fn inverse(&self, exposure: f64) -> Result<f64, ThermalTransportError> {
        let exponential = self.distance * -(-exposure / self.initial_capacity).exp_m1();
        if self.initial_capacity == self.end_capacity {
            if exponential.is_finite() && exponential >= 0.0 && exponential <= self.width {
                return Ok(exponential);
            }
            return Err(invalid(
                "constant-capacity exchange left its admitted segment",
            ));
        }
        let (mut lower, mut upper) = (0.0, self.width);
        let (mut low_value, mut high_value) = (0.0, self.exposure(self.width)?);
        let mut trial = if exponential > lower && exponential < upper {
            exponential
        } else {
            0.5 * upper
        };
        for iteration in 0..128 {
            let value = self.exposure(trial)?;
            if value == exposure || (value - exposure).abs() <= 4.0 * f64::EPSILON * exposure {
                return Ok(trial);
            }
            if value < exposure {
                lower = trial;
                low_value = value;
            } else {
                upper = trial;
                high_value = value;
            }
            if lower.next_up() >= upper {
                // Finite exposure can place the solution closer to equilibrium
                // than a representable travel. The exact supplied exposure still
                // owns the mean below; never integrate a rounded infinite log.
                return Ok(
                    if upper == self.distance || high_value - exposure < exposure - low_value {
                        upper
                    } else {
                        lower
                    },
                );
            }
            let newton =
                trial + ((exposure - value) / self.capacity(trial)) * (self.distance - trial);
            trial = if iteration % 2 == 0 && newton > lower && newton < upper {
                newton
            } else {
                lower + 0.5 * (upper - lower)
            };
        }
        Err(invalid(
            "heat-exchange inverse exceeded its bounded iteration budget",
        ))
    }

    // The same positive moments give the stream mean and its wall response.
    // B=dEout/dTw=Dout*integral C(T)/(Tw-T)^2 dT. For affine C,
    // B=C0*w0+C1*w1, w0=(1-z)*[-log(1-z)-z]/z and w1=z-w0.
    fn exchange_moments(&self, travel: f64, exposure: f64, heat: f64, gain: f64) -> (f64, f64, f64) {
        let heat_per_exposure = heat / exposure;
        let z = travel / self.distance;
        let end_capacity = self.capacity(travel);
        if z <= 0.5 {
            let [a0, a1, b0, b1] = exchange_moment_weights(z);
            let scale = self.initial_capacity.max(end_capacity);
            let (c0, c1) = (self.initial_capacity / scale, end_capacity / scale);
            let denominator = c0 * a0 + c1 * a1;
            let w0 = gain * a1;
            let response = exposure * ((c0 * w0 + c1 * (1.0 - w0)) / denominator);
            let mean_response = z * ((2.0 * c0 * b0 + c1 * (b1 - b0)) / denominator);
            let mean = if heat_per_exposure <= 0.75 * self.distance {
                self.distance - heat_per_exposure
            } else {
                travel * ((c0 * b0 + c1 * b1) / denominator)
            };
            (mean, response, mean_response)
        } else {
            let w0 = if gain == 0.0 { 0.0 } else { gain * (-gain.ln() - z) / z };
            let response = self.initial_capacity * w0 + end_capacity * (z - w0);
            let mean = if heat_per_exposure <= 0.75 * self.distance {
                self.distance - heat_per_exposure
            } else {
                let [a0, a1, b0, b1] = exchange_moment_weights(z);
                let scale = self.initial_capacity.max(end_capacity);
                let (c0, c1) = (self.initial_capacity / scale, end_capacity / scale);
                travel * ((c0 * b0 + c1 * b1) / (c0 * a0 + c1 * a1))
            };
            (mean, response, 1.0 - response / exposure)
        }
    }

    // dE_out/dE_in=(Tw-Tout)/(Tw-Tin). The original affine law gives
    // X=-C(Tw)*log(gain)-(C(Tw)-C(Tin))*travel/distance at rounded equilibrium.
    fn inlet_energy_gain(&self, travel: f64, exposure: f64) -> Result<f64, ThermalTransportError> {
        let gain = if self.initial_capacity == self.end_capacity {
            (-exposure / self.initial_capacity).exp()
        } else if travel < self.distance {
            // A resolved relative gap avoids subtracting opposing caloric
            // terms when a decreasing law has a very small reservoir Cp.
            (self.distance - travel) / self.distance
        } else if self.width == self.distance {
            let capacity_change = self.end_capacity - self.initial_capacity;
            let attenuation = if capacity_change < 0.0 {
                // Keep Cp(Tw) when subtracting it directly from a much larger
                // inlet capacity would erase it before the exposure cancels.
                let (head, tail) = fs_math::eft::two_sum(exposure, -self.initial_capacity);
                ((head + self.end_capacity) + tail) / self.end_capacity
            } else {
                exposure / self.end_capacity + capacity_change / self.end_capacity
            };
            let value = (-attenuation).exp();
            // The omitted relative gap changes this exponent by approximately
            // (Cp_in-Cp_wall)*gain/Cp_wall. It must be at rounding scale; an
            // unresolved saturated temperature must not invent a unit response.
            let gap_correction = (capacity_change.abs() * value) / self.end_capacity;
            if value == 1.0 || gap_correction > 16.0 * f64::EPSILON * attenuation.max(1.0) {
                return Err(invalid("saturated heat-exchange response is unresolved"));
            }
            value
        } else {
            return Err(invalid("heat-exchange response left its source interval"));
        };
        if !gain.is_finite() || !(0.0..=1.0).contains(&gain) {
            return Err(invalid("heat-exchange inlet response is not representable"));
        }
        Ok(gain)
    }
}

// Integrals of (1-y), y, y(1-y), y² over 1-z*y on [0,1].
// Their positive series preserve the weak-exchange limit without subtracting
// nearly equal logarithmic primitives; high-z forms keep the inlet weights finite.
fn exchange_moment_weights(z: f64) -> [f64; 4] {
    if z <= 0.5 {
        let (mut result, mut power) = ([0.0; 4], 1.0);
        let terms = if z < 0.125 { 18 } else { 64 };
        for n in 0..terms {
            let n = f64::from(n);
            result[0] += power / ((n + 1.0) * (n + 2.0));
            result[1] += power / (n + 2.0);
            result[2] += power / ((n + 2.0) * (n + 3.0));
            result[3] += power / (n + 3.0);
            power *= z;
        }
        result
    } else {
        let l = -(-z).ln_1p();
        let z2 = z * z;
        let z3 = z2 * z;
        [
            (z - (1.0 - z) * l) / z2,
            (l - z) / z2,
            (z - 0.5 * z2 - (1.0 - z) * l) / z3,
            (l - z - 0.5 * z2) / z3,
        ]
    }
}

fn compensated_add(sum: &mut f64, correction: &mut f64, value: f64) {
    let increment = value - *correction;
    let next = *sum + increment;
    *correction = (next - *sum) - increment;
    *sum = next;
}

fn invalid(reason: &'static str) -> ThermalTransportError {
    ThermalTransportError::Invalid(reason)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::Cell;

    struct CountedKnot<'a> {
        values: [f64; 2],
        queries: &'a Cell<usize>,
    }

    impl CapacityKnot for CountedKnot<'_> {
        fn temperature(&self) -> f64 {
            self.values[0]
        }

        fn capacity(&self) -> f64 {
            self.queries.set(self.queries.get() + 1);
            self.values[1]
        }
    }

    fn affine_exchange(
        initial: f64,
        final_t: f64,
        reservoir: f64,
        intercept: f64,
        slope: f64,
    ) -> f64 {
        // Independent primitive for C(T)=intercept+slope*T, used at ordinary
        // temperature separations where its subtractive form is well resolved.
        (intercept + slope * reservoir) * ((reservoir - initial) / (reservoir - final_t)).ln()
            - slope * (final_t - initial)
    }

    #[test]
    fn g1_inlet_gain_matches_constant_cp_and_equilibrium() {
        let knots = [[200.0, 1000.0], [400.0, 1000.0], [1000.0, 1000.0]];
        let curve = CapacityCurve::new(&knots);
        for (initial, wall) in [(300.0, 800.0), (800.0, 300.0), (350.0, 350.0)] {
            for exposure in [0.0_f64, 1e-6, 500.0, 50_000.0, 1e6] {
                let answer = curve.linear_heat_exchange(initial, wall, exposure).unwrap();
                let expected = (-exposure / 1000.0).exp();
                // Exp roundoff scales with exposure/Cp.
                let bound = 8.0 * f64::EPSILON * (1.0 + exposure / 1000.0) * expected;
                assert!((answer.inlet_energy_gain - expected).abs() <= bound);
            }
        }
    }

    #[test]
    fn g1_inlet_gain_crosses_knots_and_matches_independent_energy_perturbations() {
        let knots = [[200.0, 1200.0], [300.0, 1400.0], [400.0, 1800.0], [600.0, 1000.0]];
        let curve = CapacityCurve::new(&knots);
        for (initial, wall) in [(250.0, 550.0), (550.0, 250.0)] {
            let answer = curve.linear_heat_exchange(initial, wall, 4000.0).unwrap();
            let expected = (wall - answer.temperature) / (wall - initial);
            assert!((answer.inlet_energy_gain - expected).abs() < 2e-14);
            let inputs = [initial - 1e-3, initial + 1e-3];
            let outputs = inputs.map(|t| curve.linear_heat_exchange(t, wall, 4000.0).unwrap().temperature);
            let measured = curve.energy_change(outputs[0], outputs[1]).unwrap()
                / curve.energy_change(inputs[0], inputs[1]).unwrap();
            assert!((answer.inlet_energy_gain - measured).abs() < 2e-8);
        }
    }

    #[test]
    fn g1_affine_inlet_gain_retains_weak_and_saturated_response() {
        for (intercept, slope) in [(800.0, 2.0), (3000.0, -2.0)] {
            let knots = [200.0, 400.0, 800.0].map(|t| [t, intercept + slope * t]);
            let curve = CapacityCurve::new(&knots);
            for (initial, wall) in [(300.0, 600.0), (600.0, 300.0)] {
                for attenuation in [1e-12_f64, 50.0, 800.0] {
                    // Independent affine primitive at a prescribed relative gap.
                    let travel = (wall - initial) * -(-attenuation).exp_m1();
                    let exposure = (intercept + slope * wall) * attenuation - slope * travel;
                    let answer = curve.linear_heat_exchange(initial, wall, exposure).unwrap();
                    let expected = (-attenuation).exp();
                    assert!((answer.inlet_energy_gain - expected).abs() <=
                        if attenuation < 1.0 { 32.0 * f64::EPSILON } else { 2e-12 * expected });
                    if attenuation >= 50.0 {
                        assert_eq!(answer.temperature, wall);
                    }
                    if attenuation == 50.0 {
                        assert!(answer.inlet_energy_gain > 0.0);
                    }
                }
            }
        }
    }

    #[test]
    fn g1_constant_capacity_exchange_matches_newton_cooling_and_keeps_small_steps() {
        let knots = [[200.0, 1000.0], [1000.0, 1000.0]];
        let curve = CapacityCurve::new(&knots);
        for (initial, reservoir) in [(300.0, 800.0), (800.0, 300.0)] {
            for exposure in [1e-6_f64, 10.0, 500.0, 2000.0, 1e6] {
                let expected = reservoir + (initial - reservoir) * (-exposure / 1000.0).exp();
                let actual = curve
                    .temperature_after_linear_exchange(initial, reservoir, exposure)
                    .unwrap();
                assert!((actual - expected).abs() <= 2.0 * f64::EPSILON * expected);
                assert!(actual >= initial.min(reservoir) && actual <= initial.max(reservoir));
                assert_eq!(
                    actual,
                    curve
                        .temperature_after_linear_exchange(initial, reservoir, exposure)
                        .unwrap()
                );
            }
        }
        let initial: f64 = 300.0;
        let target = initial.next_up();
        let exposure = 1000.0 * ((target - initial) / (800.0 - target)).ln_1p();
        assert_eq!(
            curve
                .temperature_after_linear_exchange(initial, 800.0, exposure)
                .unwrap(),
            target
        );
        assert!(
            curve
                .temperature_after_linear_exchange(initial, 800.0, f64::from_bits(1))
                .is_err()
        );
    }

    #[test]
    fn g1_weak_exchange_retains_relative_heat_and_direct_mean() {
        for slope in [0.0, 2.0] {
            let intercept = if slope == 0.0 { 1000.0 } else { 800.0 };
            let knots = [
                [200.0, intercept + slope * 200.0],
                [800.0, intercept + slope * 800.0],
            ];
            let curve = CapacityCurve::new(&knots);
            for initial in [300.0, 500.0] {
                let reservoir = 400.0;
                let gap = reservoir - initial;
                let cp = intercept + slope * initial;
                for ntu in [1e-9, 1e-12] {
                    let exposure = cp * ntu;
                    // Independent second-order inverse of
                    // x=z+(1+slope*gap/cp)*z²/2+O(z³).
                    // O(x³) is below f64 relative heat resolution here.
                    let delta = gap * ntu * (1.0 - 0.5 * (1.0 + slope * gap / cp) * ntu);
                    let heat = cp * delta + 0.5 * slope * delta * delta;
                    let answer = curve
                        .linear_heat_exchange(initial, reservoir, exposure)
                        .unwrap();
                    assert!(
                        (answer.energy_change - heat).abs() <= 64.0 * f64::EPSILON * heat.abs(),
                        "slope={slope} initial={initial} ntu={ntu}: heat {} != {heat}",
                        answer.energy_change
                    );
                    assert!(
                        (answer.temperature - (initial + delta)).abs() <= f64::EPSILON * initial
                    );
                    let mean = initial + 0.5 * gap * ntu;
                    assert!((answer.mean_temperature - mean).abs() <= 2.0 * f64::EPSILON * initial);
                    assert!(
                        (answer.energy_change - exposure * (reservoir - answer.mean_temperature))
                            .abs()
                            <= 8.0 * f64::EPSILON * exposure * initial.max(reservoir)
                    );
                    assert!(answer.mean_temperature >= initial.min(answer.temperature));
                    assert!(answer.mean_temperature <= initial.max(answer.temperature));
                }
            }
        }
    }

    #[test]
    fn g1_exchange_heat_and_mean_cross_knots_and_retain_finite_saturation() {
        let knots = [
            [200.0, 1200.0],
            [300.0, 1400.0],
            [400.0, 1800.0],
            [600.0, 1000.0],
        ];
        let curve = CapacityCurve::new(&knots);
        for (initial, reservoir) in [(250.0, 550.0), (550.0, 250.0)] {
            for exposure in [500.0, 2000.0, 1e6, 1e12] {
                let answer = curve
                    .linear_heat_exchange(initial, reservoir, exposure)
                    .unwrap();
                let heat = curve.energy_change(initial, answer.temperature).unwrap();
                assert!((answer.energy_change - heat).abs() <= 32.0 * f64::EPSILON * heat.abs());
                assert!(
                    (answer.energy_change - exposure * (reservoir - answer.mean_temperature)).abs()
                        <= 8.0 * f64::EPSILON * exposure * initial.max(reservoir)
                );
                assert!(answer.mean_temperature >= initial.min(answer.temperature));
                assert!(answer.mean_temperature <= initial.max(answer.temperature));
            }
        }
    }

    #[test]
    fn g1_affine_capacity_exchange_matches_exact_heating_cooling_and_composition() {
        for (intercept, slope) in [(800.0, 2.0), (2200.0, -2.0)] {
            let knots = [
                [200.0, intercept + 200.0 * slope],
                [600.0, intercept + 600.0 * slope],
            ];
            let curve = CapacityCurve::new(&knots);
            for (initial, target, reservoir) in [(300.0, 500.0, 700.0), (550.0, 280.0, 180.0)] {
                let exposure = affine_exchange(initial, target, reservoir, intercept, slope);
                let actual = curve
                    .temperature_after_linear_exchange(initial, reservoir, exposure)
                    .unwrap();
                assert!(
                    (actual - target).abs() <= 4.0 * f64::EPSILON * target,
                    "{actual} != {target}"
                );
                let halfway = curve
                    .temperature_after_linear_exchange(initial, reservoir, 0.5 * exposure)
                    .unwrap();
                let composed = curve
                    .temperature_after_linear_exchange(halfway, reservoir, 0.5 * exposure)
                    .unwrap();
                assert!((composed - actual).abs() <= 4.0 * f64::EPSILON * target);
                let heat = curve.energy_change(initial, actual).unwrap();
                let exact_heat =
                    (intercept + 0.5 * slope * (target + initial)) * (target - initial);
                assert!((heat - exact_heat).abs() <= 16.0 * f64::EPSILON * exact_heat.abs());
            }
        }
    }

    #[test]
    fn g1_exchange_crosses_original_capacity_knots_in_both_directions() {
        let knots = [
            [200.0, 1200.0],
            [300.0, 1400.0],
            [400.0, 1800.0],
            [600.0, 1000.0],
        ];
        let curve = CapacityCurve::new(&knots);
        for (initial, target, reservoir) in [(250.0, 550.0, 700.0), (550.0, 250.0, 150.0)] {
            let cuts = if initial < target {
                [initial, 300.0, 400.0, target]
            } else {
                [initial, 400.0, 300.0, target]
            };
            let exposure: f64 = cuts
                .windows(2)
                .map(|pair| {
                    let middle = 0.5 * (pair[0] + pair[1]);
                    let (intercept, slope) = if middle < 300.0 {
                        (800.0, 2.0)
                    } else if middle < 400.0 {
                        (200.0, 4.0)
                    } else {
                        (3400.0, -4.0)
                    };
                    affine_exchange(pair[0], pair[1], reservoir, intercept, slope)
                })
                .sum();
            // Independent derivative of the original integral at fixed X.
            // The reservoir is outside the source support in both directions.
            let integral: f64 = cuts.windows(2).map(|pair| {
                let middle = 0.5 * (pair[0] + pair[1]);
                let (intercept, slope) = if middle < 300.0 {
                    (800.0, 2.0)
                } else if middle < 400.0 {
                    (200.0, 4.0)
                } else {
                    (3400.0, -4.0)
                };
                let (a, b) = (reservoir - pair[0], reservoir - pair[1]);
                (intercept + slope * reservoir) * (1.0 / b - 1.0 / a)
                    + slope * (b / a).ln()
            }).sum();
            let expected_response = (reservoir - target) * integral;
            let answer = curve.linear_heat_exchange(initial, reservoir, exposure).unwrap();
            let actual = answer.temperature;
            assert!(
                (actual - target).abs() <= 8.0 * f64::EPSILON * target,
                "{actual} != {target}"
            );
            assert!((answer.reservoir_energy_derivative - expected_response).abs() < 3e-9);
            assert!((answer.mean_reservoir_derivative - (1.0 - expected_response / exposure)).abs() < 3e-12);
            assert!((answer.inlet_energy_loss - (target - initial) / (reservoir - initial)).abs() < 3e-14);
            assert!((answer.inlet_energy_gain + answer.inlet_energy_loss - 1.0).abs() < 3e-14);
        }
        // Exact source knots select the correct outgoing interval when cooling.
        let exposure = affine_exchange(400.0, 350.0, 200.0, 200.0, 4.0);
        assert!(
            (curve
                .temperature_after_linear_exchange(400.0, 200.0, exposure)
                .unwrap()
                - 350.0)
                .abs()
                < 1e-12
        );
    }

    #[test]
    fn g1_reservoir_response_retains_weak_exchange_and_equilibrium_limits() {
        let knots = [[200.0, 1200.0], [800.0, 2400.0]]; // C(T)=800+2T.
        let curve = CapacityCurve::new(&knots);
        let (temperature, capacity) = (400.0, 1600.0);
        for attenuation in [0.0_f64, 1e-12, 1e-5, 1.0, 50.0] {
            let exposure = capacity * attenuation;
            let answer = curve.linear_heat_exchange(temperature, temperature, exposure).unwrap();
            let expected_loss = -(-attenuation).exp_m1();
            let expected_mean = if attenuation < 1e-3 {
                attenuation * (0.5 - attenuation / 6.0 + attenuation * attenuation / 24.0)
            } else {
                1.0 - expected_loss / attenuation
            };
            assert_eq!(answer.temperature, temperature);
            assert_eq!(answer.energy_change, 0.0);
            assert_eq!(answer.inlet_energy_loss, expected_loss);
            assert!((answer.reservoir_energy_derivative - capacity * expected_loss).abs()
                <= 8.0 * f64::EPSILON * (capacity * expected_loss));
            assert!((answer.mean_reservoir_derivative - expected_mean).abs()
                <= 8.0 * f64::EPSILON * expected_mean);
        }
        // The tiny derivative of the Robin mean is retained directly, even
        // when forming 1-B/X would lose its significant digits.
        for reservoir in [300.0, 500.0] {
            let exposure = capacity * 1e-12;
            let answer = curve.linear_heat_exchange(temperature, reservoir, exposure).unwrap();
            assert!((answer.reservoir_energy_derivative / exposure - 1.0).abs() < 2e-12);
            assert!((answer.mean_reservoir_derivative / 5e-13 - 1.0).abs() < 2e-11);
            assert!(answer.inlet_energy_loss > 0.0);
        }
    }

    #[test]
    fn g0_exchange_preserves_equilibrium_and_refuses_unsupported_paths() {
        let knots = [[300.0, 1000.0], [400.0, 1000.0]];
        let curve = CapacityCurve::new(&knots);
        for (initial, end, reservoir) in [(300.0, 400.0, 500.0), (400.0, 300.0, 200.0)] {
            let available = curve.exchange_integral(0, initial, end, reservoir).unwrap();
            assert_eq!(
                curve
                    .temperature_after_linear_exchange(initial, reservoir, available)
                    .unwrap(),
                end
            );
            assert!(
                curve
                    .temperature_after_linear_exchange(initial, reservoir, available.next_up())
                    .is_err()
            );
            assert!(
                curve
                    .temperature_after_linear_exchange(end, reservoir, 1.0)
                    .is_err()
            );
            assert_eq!(
                curve
                    .temperature_after_linear_exchange(initial, reservoir, 0.0)
                    .unwrap(),
                initial
            );
        }
        assert_eq!(
            curve
                .temperature_after_linear_exchange(350.0, 350.0, 1000.0)
                .unwrap(),
            350.0
        );
        for exposure in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(
                curve
                    .temperature_after_linear_exchange(350.0, 360.0, exposure)
                    .is_err()
            );
        }
        for reservoir in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(
                curve
                    .temperature_after_linear_exchange(350.0, reservoir, 1.0)
                    .is_err()
            );
        }
        for initial in [299.0, 401.0, f64::NAN, f64::INFINITY] {
            assert!(
                curve
                    .temperature_after_linear_exchange(initial, 350.0, 0.0)
                    .is_err()
            );
        }
    }

    #[test]
    fn g1_exchange_retains_small_heat_far_from_reservoir_and_stops_at_equilibrium() {
        let knots = [[200.0, 1200.0], [600.0, 2000.0]];
        let curve = CapacityCurve::new(&knots);
        let (initial, target, reservoir) = (300.0_f64, 300.0_f64.next_up(), 1e12_f64);
        assert_eq!(reservoir - initial, reservoir - target);
        let width = target - initial;
        // Simpson quadrature of Cp/(reservoir-T) across a one-ulp interval;
        // truncation is far below f64 resolution and uses no logarithmic weights.
        let value = |t: f64| (800.0 + 2.0 * t) / (reservoir - t);
        let exposure =
            width / 6.0 * (value(initial) + 4.0 * value(initial + 0.5 * width) + value(target));
        assert_eq!(
            curve
                .temperature_after_linear_exchange(initial, reservoir, exposure)
                .unwrap(),
            target
        );
        for reservoir in [250.0, 500.0] {
            let result = curve
                .temperature_after_linear_exchange(initial, reservoir, 1e12)
                .unwrap();
            assert!((result - reservoir).abs() <= f64::EPSILON * reservoir);
            assert!(result >= initial.min(reservoir) && result <= initial.max(reservoir));
        }
    }

    #[test]
    fn g1_exchange_does_not_walk_unvisited_capacity_intervals() {
        for count in [64, 4096] {
            let queries = Cell::new(0);
            let knots: Vec<_> = (0..count)
                .map(|i| CountedKnot {
                    values: [250.0 + f64::from(i), 1000.0 + 0.5 * f64::from(i)],
                    queries: &queries,
                })
                .collect();
            let curve = CapacityCurve::new(&knots);
            let initial = 250.25 + f64::from(count / 2);
            let target = initial + 0.125;
            let reservoir = initial + 0.5;
            let exposure = affine_exchange(initial, target, reservoir, 875.0, 0.5);
            let actual = curve
                .temperature_after_linear_exchange(initial, reservoir, exposure)
                .unwrap();
            assert!((actual - target).abs() <= 4.0 * f64::EPSILON * target);
            assert!(
                queries.get() <= 1024,
                "local exchange read {} capacities",
                queries.get()
            );
        }
    }

    #[test]
    fn g1_local_nonlinear_calorimetry_does_not_scan_unvisited_source_intervals() {
        for count in [64, 4096] {
            let queries = Cell::new(0);
            let knots: Vec<_> = (0..count)
                .map(|i| CountedKnot {
                    values: [250.0 + f64::from(i), 1000.0 + 0.5 * f64::from(i)],
                    queries: &queries,
                })
                .collect();
            let curve = CapacityCurve::new(&knots);
            let initial = 250.25 + f64::from(count / 2);
            let initial_capacity = 1000.0 + 0.5 * (initial - 250.0);
            for delta in [-0.125, 0.125] {
                // Independent exact integral of C(T)=1000+0.5*(T-250).
                let heat = initial_capacity * delta + 0.25 * delta * delta;
                queries.set(0);
                let final_t = curve
                    .temperature_after_energy_change(initial, heat)
                    .unwrap();
                assert_eq!(final_t, initial + delta);
                assert!(
                    queries.get() <= 16,
                    "local heat read {} capacities",
                    queries.get()
                );
                let recovered = curve.energy_change(initial, final_t).unwrap();
                assert!((recovered - heat).abs() <= 4.0 * f64::EPSILON * heat.abs());
            }
        }
    }

    #[test]
    fn g1_crossing_inverse_retains_exact_knots_support_and_signed_caloric_balance() {
        let knots = [
            [200.0, 1000.0],
            [300.0, 2000.0],
            [400.0, 3000.0],
            [500.0, 2000.0],
        ];
        let curve = CapacityCurve::new(&knots);
        for initial in [200.0, 250.0, 300.0, 350.0, 400.0, 500.0] {
            for final_t in [200.0, 275.0, 300.0, 350.0, 400.0, 425.0, 500.0] {
                let heat = curve.energy_change(initial, final_t).unwrap();
                let recovered = curve
                    .temperature_after_energy_change(initial, heat)
                    .unwrap();
                assert!((recovered - final_t).abs() <= 2.0 * f64::EPSILON * final_t);
                if final_t == 200.0 || final_t == 500.0 {
                    assert_eq!(recovered, final_t);
                }
            }
            for boundary in [200.0, 500.0] {
                let available = curve.energy_change(initial, boundary).unwrap();
                let beyond = if boundary == 500.0 {
                    available.next_up()
                } else {
                    available.next_down()
                };
                assert!(
                    curve
                        .temperature_after_energy_change(initial, beyond)
                        .is_err()
                );
            }
        }
    }

    #[test]
    fn g1_fixed_datum_inverse_uses_local_work_far_across_long_source_tables() {
        for count in [64, 4096, 262_144] {
            let queries = Cell::new(0);
            let knots: Vec<_> = (0..count)
                .map(|i| CountedKnot {
                    values: [250.0 + f64::from(i), 1000.0 + 0.5 * f64::from(i)],
                    queries: &queries,
                })
                .collect();
            let curve = CapacityCurve::new(&knots);
            let datum = 250.125 + f64::from(count / 2);
            let inverse = CapacityInverse::new(&curve, datum).unwrap();
            for target in [
                250.25,
                datum - 0.25,
                datum + 0.25,
                248.75 + f64::from(count),
            ] {
                let delta = target - datum;
                // Independent exact integral of C(T)=1000+0.5*(T-250).
                let heat = (1000.0 + 0.25 * (datum + target - 500.0)) * delta;
                queries.set(0);
                let actual = inverse.temperature(&curve, heat).unwrap();
                assert!(
                    (actual - target).abs() <= 8.0 * f64::EPSILON * target,
                    "{count} knots, {actual} != {target}"
                );
                assert!(
                    queries.get() <= 16,
                    "fixed-datum inversion read {} capacities on {count} knots",
                    queries.get()
                );
            }
        }
    }

    #[test]
    fn g1_datum_aligned_forward_calorimetry_uses_bounded_source_work() {
        for count in [64, 4096, 262_144] {
            let queries = Cell::new(0);
            let knots: Vec<_> = (0..count)
                .map(|i| CountedKnot {
                    values: [250.0 + f64::from(i), 1000.0 + 0.5 * f64::from(i)],
                    queries: &queries,
                })
                .collect();
            let curve = CapacityCurve::new(&knots);
            let datum = 250.125 + f64::from(count / 2);
            let index = CapacityInverse::new(&curve, datum).unwrap();
            for target in [250.25, datum - 0.25, datum, datum + 0.25, 248.75 + f64::from(count)] {
                // Independent integral of the source C(T)=1000+0.5*(T-250).
                let expected = (1000.0 + 0.25 * (datum + target - 500.0)) * (target - datum);
                queries.set(0);
                let actual = index.energy_change_from_datum(&curve, target).unwrap();
                assert!((actual - expected).abs() <= 8.0 * f64::EPSILON * expected.abs());
                assert!(queries.get() <= 8, "forward query read {} capacities", queries.get());
                let recovered = index.temperature(&curve, actual).unwrap();
                assert!((recovered - target).abs() <= 8.0 * f64::EPSILON * target);
            }
        }
    }

    #[test]
    fn g1_indexed_forward_heat_retains_datum_neighbors_and_canonical_endpoints() {
        let knots = [[200.0, 500.0], [300.0, 2000.0], [400.0, 4000.0], [500.0, 300.0]];
        let curve = CapacityCurve::new(&knots);
        for datum in [200.0_f64, 261.3, 300.0, 325.3, 500.0] {
            let index = CapacityInverse::new(&curve, datum).unwrap();
            for target in [200.0, 300.0, 400.0, 500.0, datum.next_down(), datum, datum.next_up()] {
                if !(200.0..=500.0).contains(&target) { continue; }
                let exact = curve.energy_change(datum, target).unwrap();
                let indexed = index.energy_change_from_datum(&curve, target).unwrap();
                if target == 200.0 || target == 500.0 || (target - datum).abs() < 1e-10 {
                    assert_eq!(indexed, exact);
                } else {
                    assert!((indexed - exact).abs() <= 8.0 * f64::EPSILON * exact.abs());
                }
                if target != datum { assert_ne!(indexed, 0.0); }
            }
            for target in [199.0, 501.0, f64::NAN, f64::INFINITY] {
                assert!(index.energy_change_from_datum(&curve, target).is_err());
            }
        }
    }

    #[test]
    fn g1_indexed_inverse_preserves_nonlinear_heating_cooling_and_source_support() {
        // Each interval deliberately has a different slope, including cooling
        // through a decreasing capacity. A linear interpolation of h knots
        // would not recover these temperatures.
        for scale in [1e-12, 1.0, 1e12] {
            let knots = [
                [200.0, 1000.0 * scale],
                [275.0, 7000.0 * scale],
                [310.0, 2000.0 * scale],
                [400.0, 4000.0 * scale],
                [600.0, 3000.0 * scale],
            ];
            let curve = CapacityCurve::new(&knots);
            for datum in [200.0, 261.75, 310.0, 511.25, 600.0] {
                let inverse = CapacityInverse::new(&curve, datum).unwrap();
                assert_eq!(inverse.temperature(&curve, 0.0).unwrap(), datum);
                for target in [
                    200.0, 231.25, 275.0, 301.5, 310.0, 367.0, 400.0, 521.0, 600.0,
                ] {
                    let heat = curve.energy_change(datum, target).unwrap();
                    let actual = inverse.temperature(&curve, heat).unwrap();
                    assert!(
                        (actual - target).abs() <= 8.0 * f64::EPSILON * target,
                        "scale={scale}, datum={datum}, {actual} != {target}"
                    );
                    if target == 200.0 || target == 600.0 {
                        assert_eq!(actual, target);
                    }
                }
                for boundary in [200.0, 600.0] {
                    let heat = curve.energy_change(datum, boundary).unwrap();
                    let beyond = if boundary == 200.0 {
                        heat.next_down()
                    } else {
                        heat.next_up()
                    };
                    assert!(inverse.temperature(&curve, beyond).is_err());
                }
                assert!(inverse.temperature(&curve, f64::NAN).is_err());
                assert!(inverse.temperature(&curve, f64::INFINITY).is_err());
            }
        }
    }

    #[test]
    fn g1_index_does_not_erase_small_heat_near_its_interior_datum() {
        let knots = [[200.0, 1000.0], [300.0, 2000.0], [400.0, 4000.0]];
        let curve = CapacityCurve::new(&knots);
        let datum: f64 = 325.0;
        let inverse = CapacityInverse::new(&curve, datum).unwrap();
        for target in [datum.next_down(), datum.next_up()] {
            let heat = curve.energy_change(datum, target).unwrap();
            assert_ne!(heat, 0.0);
            assert_eq!(inverse.temperature(&curve, heat).unwrap(), target);
        }
        assert!(inverse.temperature(&curve, f64::from_bits(1)).is_err());

        // A non-dyadic datum exposes cancellation against a distant cached
        // prefix: it must not manufacture warming or cooling from tiny heat.
        let knots = [[200.0, 500.0], [400.0, 4000.0]];
        let curve = CapacityCurve::new(&knots);
        let datum: f64 = 325.3;
        let inverse = CapacityInverse::new(&curve, datum).unwrap();
        for heat in [f64::from_bits(1), -f64::from_bits(1)] {
            assert!(inverse.temperature(&curve, heat).is_err());
        }
        for target in [datum.next_down(), datum.next_up()] {
            let heat = curve.energy_change(datum, target).unwrap();
            assert_eq!(inverse.temperature(&curve, heat).unwrap(), target);
        }
    }
}
