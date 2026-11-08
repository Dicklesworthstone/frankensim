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
    pub remaining_temperature_difference: f64,
    pub mean_log_exposure_derivative: f64,
}


/// Same exact response for C(T) dT/du = drive - feedback*(T-T_initial).
/// Drive, feedback and duration carry the caller's original compatible units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct CapacityDrivenExchange {
    pub temperature: f64,
    pub temperature_change: f64,
    pub energy_change: f64,
    pub mean_temperature: f64,
    pub mean_temperature_change: f64,
    pub inlet_energy_gain: f64,
    pub inlet_energy_loss: f64,
    pub source_gain: f64,
    pub mean_source_derivative: f64,
    pub remaining_drive: f64,
    pub mean_log_duration_derivative: f64,
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
        let direct = || {
            curve.temperature_after_energy_change_preserving_energy(self.datum, energy)
        };
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
            && let Ok(temperature) =
                curve.temperature_after_energy_change_preserving_energy(start, change)
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
        self.exchange::<false>(initial, reservoir, exposure)
    }

    /// Native energy owners retain the relative heat even when absolute T
    /// cannot move. Temperature-only callers keep their original refusal.
    pub(super) fn linear_heat_exchange_preserving_heat(
        &self,
        initial: f64,
        reservoir: f64,
        exposure: f64,
    ) -> Result<CapacityExchange, ThermalTransportError> {
        self.exchange::<true>(initial, reservoir, exposure)
    }

    fn exchange<const PRESERVE_HEAT: bool>(
        &self,
        initial: f64,
        reservoir: f64,
        exposure: f64,
    ) -> Result<CapacityExchange, ThermalTransportError> {
        if !positive(reservoir) {
            return Err(invalid("linear heat exchange needs positive reservoir temperature"));
        }
        // Keep the original bath as an endpoint hint, without multiplying its
        // temperature difference by a potentially very large exposure.
        let answer = self.exchange_driven::<PRESERVE_HEAT>(
            initial, reservoir - initial, 1.0, exposure, Some(reservoir),
        )?;
        Ok(CapacityExchange {
            temperature: answer.temperature,
            energy_change: answer.energy_change,
            mean_temperature: answer.mean_temperature,
            inlet_energy_gain: answer.inlet_energy_gain,
            inlet_energy_loss: answer.inlet_energy_loss,
            reservoir_energy_derivative: answer.source_gain,
            mean_reservoir_derivative: answer.mean_source_derivative,
            remaining_temperature_difference: answer.remaining_drive,
            mean_log_exposure_derivative: answer.mean_log_duration_derivative,
        })
    }

    /// Normalized-time source and linear-loss exchange. No absolute virtual
    /// reservoir or division of initial drive by feedback is required.
    pub(super) fn driven_heat_exchange(
        &self,
        initial: f64,
        initial_drive: f64,
        exposure: f64,
    ) -> Result<CapacityDrivenExchange, ThermalTransportError> {
        self.exchange_driven::<true>(initial, initial_drive, exposure, 1.0, None)
    }

    fn exchange_driven<const PRESERVE_HEAT: bool>(
        &self,
        initial: f64,
        initial_drive: f64,
        feedback: f64,
        duration: f64,
        reservoir: Option<f64>,
    ) -> Result<CapacityDrivenExchange, ThermalTransportError> {
        self.check_temperature(initial)?;
        if !initial_drive.is_finite()
            || !feedback.is_finite() || feedback < 0.0
            || !duration.is_finite() || duration < 0.0
        {
            return Err(invalid("driven exchange needs finite drive and nonnegative feedback and duration"));
        }
        if duration == 0.0 || initial_drive == 0.0 {
            let capacity = self.capacity(initial)?;
            let (gain, loss, response, mean_response) =
                equilibrium_response(capacity, feedback, duration);
            return Ok(CapacityDrivenExchange {
                temperature: initial,
                temperature_change: 0.0,
                energy_change: 0.0,
                mean_temperature: initial,
                mean_temperature_change: 0.0,
                inlet_energy_gain: gain,
                inlet_energy_loss: loss,
                source_gain: response,
                mean_source_derivative: mean_response,
                remaining_drive: initial_drive,
                mean_log_duration_derivative: 0.0,
            });
        }
        let heating = initial_drive > 0.0;
        let direction = if heating { 1.0 } else { -1.0 };
        let drive_magnitude = initial_drive.abs();
        let mut index = if heating {
            self.segment(initial)
        } else {
            self.knots.partition_point(|k| k.temperature() < initial).saturating_sub(1)
        };
        let mut current = initial;
        let (mut cumulative, mut correction) = (0.0, 0.0);
        let (mut energy, mut energy_correction) = (0.0, 0.0);
        let (mut mean_travel, mut mean_correction) = (0.0, 0.0);
        let (mut inlet_energy_gain, mut inlet_energy_loss) = (1.0, 0.0);
        let mut source_gain = 0.0;
        let (mut mean_response_integral, mut response_correction) = (0.0, 0.0);
        let (mut mean_duration_response, mut duration_response_correction) = (0.0, 0.0);
        for _ in 0..self.knots.len() - 1 {
            let boundary = if heating {
                self.knots[index + 1].temperature()
            } else {
                self.knots[index].temperature()
            };
            let source_width = (boundary - current).abs();
            if source_width == 0.0 {
                break;
            }
            let drive = reservoir.map_or_else(
                || (-feedback).mul_add((current - initial).abs(), drive_magnitude),
                |bath| (bath - current).abs(),
            );
            if !positive(drive) {
                return Err(invalid("source-relative driving power is not representable"));
            }
            // An infinite drive/feedback ratio is harmless: the finite original
            // source interval still bounds the travel, and z is computed directly.
            let equilibrium_width = if feedback == 0.0 { f64::INFINITY } else { drive / feedback };
            let width = source_width.min(equilibrium_width);
            if !positive(width) {
                return Err(invalid("driven exchange travel is not representable"));
            }
            let end = if width == source_width {
                boundary
            } else {
                reservoir.unwrap_or(current + direction * width)
            };
            let initial_capacity = self.on_segment(index, current);
            let boundary_capacity = self.on_segment(index, boundary);
            let fraction = width / source_width;
            let end_capacity = if width == source_width {
                boundary_capacity
            } else {
                (1.0 - fraction) * initial_capacity + fraction * boundary_capacity
            };
            let segment = ExchangeSegment {
                drive, feedback, width, initial_capacity, end_capacity,
            };
            let amount = segment.exposure(width)?;
            let corrected = amount - correction;
            let next = cumulative + corrected;
            let finished = duration <= next;
            let local_duration = if finished { duration - cumulative } else { amount };
            let travel = if !finished || duration == next {
                width
            } else {
                segment.inverse(local_duration)?
            };
            let heat = (0.5 * initial_capacity + 0.5 * segment.capacity(travel)) * travel;
            if !positive(heat) || !positive(local_duration) {
                return Err(invalid("driven exchange heat or duration is not representable"));
            }
            let gain = segment.inlet_energy_gain(travel, local_duration)?;
            let loss = if initial_capacity == end_capacity {
                -(-product_over(feedback, local_duration, initial_capacity)).exp_m1()
            } else {
                segment.fraction(travel)
            };
            let (local_mean, response, mean_response, local_mean_gap) =
                segment.exchange_moments(travel, local_duration, heat, gain);
            // The endpoint-minus-mean moment retains even sub-ULP absolute T.
            compensated_add(
                &mut mean_duration_response,
                &mut duration_response_correction,
                (cumulative / duration) * travel + (local_duration / duration) * local_mean_gap,
            );
            // (u - dE/dD)/feedback, accumulated without subtraction or
            // division by feedback. travel/drive is the exact zero-feedback limit.
            let previous_response = if initial_capacity == end_capacity {
                product_over(response, source_gain, initial_capacity)
            } else {
                product_over(travel, source_gain, drive)
            };
            compensated_add(
                &mut mean_response_integral,
                &mut response_correction,
                previous_response + local_duration * mean_response,
            );
            source_gain = gain.mul_add(source_gain, response);
            inlet_energy_loss = gain.mul_add(inlet_energy_loss, loss);
            inlet_energy_gain *= gain;
            if !source_gain.is_finite()
                || !mean_response_integral.is_finite()
                || !mean_duration_response.is_finite()
                || source_gain < 0.0 || mean_response_integral < 0.0 || mean_duration_response < 0.0
            {
                return Err(invalid("driven exchange response is not representable"));
            }
            compensated_add(&mut energy, &mut energy_correction, heat);
            compensated_add(
                &mut mean_travel,
                &mut mean_correction,
                (local_duration / duration) * ((current - initial).abs() + local_mean),
            );
            if finished {
                let temperature = if travel == width {
                    end
                } else if let Some(bath) = reservoir
                    && travel > 0.5 * drive
                {
                    bath - direction * (drive - travel)
                } else {
                    current + direction * travel
                };
                if temperature == initial && !PRESERVE_HEAT {
                    return Err(invalid("nonzero heat exchange is below temperature resolution"));
                }
                let temperature_change = direction * ((current - initial).abs() + travel);
                let mean_temperature_change = direction * mean_travel;
                let mean_temperature = initial + mean_temperature_change;
                if !positive(temperature) || !positive(mean_temperature)
                    || mean_temperature < initial.min(temperature)
                    || mean_temperature > initial.max(temperature)
                    || temperature < current.min(end) || temperature > current.max(end)
                {
                    return Err(invalid("driven exchange state or mean is not representable"));
                }
                // With no loss the exact source heat is the native input, not
                // a reconstruction from its rounded inverse temperature travel.
                let energy_change = if feedback == 0.0 {
                    initial_drive * duration
                } else {
                    direction * energy
                };
                if !energy_change.is_finite() {
                    return Err(invalid("driven exchange energy is not representable"));
                }
                return Ok(CapacityDrivenExchange {
                    temperature, temperature_change, energy_change, mean_temperature,
                    mean_temperature_change, inlet_energy_gain, inlet_energy_loss,
                    source_gain: if feedback == 0.0 { duration } else { source_gain },
                    mean_source_derivative: mean_response_integral / duration,
                    remaining_drive: initial_drive * inlet_energy_gain,
                    mean_log_duration_derivative: direction * mean_duration_response,
                });
            }
            correction = (next - cumulative) - corrected;
            cumulative = next;
            current = end;
            if heating {
                if index + 2 == self.knots.len() { break; }
                index += 1;
            } else {
                if index == 0 { break; }
                index -= 1;
            }
        }
        Err(invalid("heat exchange exceeds the source temperature support"))
    }

    /// Exposure required to traverse part of one original affine-Cp interval.
    /// Reaching the reservoir itself takes infinite exposure. This shares the
    /// forward integral used by the bounded exchange inverse.
    pub(super) fn exchange_integral(
        &self,
        index: usize,
        initial: f64,
        final_t: f64,
        reservoir: f64,
    ) -> Result<f64, ThermalTransportError> {
        let Some(upper) = index.checked_add(1).and_then(|i| self.knots.get(i)) else {
            return Err(invalid("heat-exchange interval is outside source support"));
        };
        let lower = &self.knots[index];
        if !positive(reservoir)
            || !initial.is_finite() || !final_t.is_finite()
            || initial < lower.temperature() || initial > upper.temperature()
            || final_t < lower.temperature() || final_t > upper.temperature()
            || (final_t > initial && final_t > reservoir)
            || (final_t < initial && final_t < reservoir)
        {
            return Err(invalid("heat-exchange exposure needs a supported path toward the reservoir"));
        }
        if initial == final_t { return Ok(0.0); }
        let segment = ExchangeSegment {
            drive: (reservoir - initial).abs(),
            feedback: 1.0,
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
        self.invert_energy_change::<false>(initial, change)
    }

    /// A retained native energy can resolve a change below one temperature ULP.
    /// Return its rounded T without erasing or rejecting that stored energy.
    /// The caller retains the energy coordinate and its original support check.
    pub(super) fn temperature_after_energy_change_preserving_energy(
        &self,
        initial: f64,
        change: f64,
    ) -> Result<f64, ThermalTransportError> {
        self.invert_energy_change::<true>(initial, change)
    }

    fn invert_energy_change<const PRESERVE_ENERGY: bool>(
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
                    || (temperature == initial && !PRESERVE_ENERGY)
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
    drive: f64,
    feedback: f64,
    width: f64,
    initial_capacity: f64,
    end_capacity: f64,
}

impl ExchangeSegment {
    fn capacity(&self, travel: f64) -> f64 {
        let fraction = travel / self.width;
        (1.0 - fraction) * self.initial_capacity + fraction * self.end_capacity
    }

    fn fraction(&self, travel: f64) -> f64 {
        if self.feedback > 0.0
            && travel == self.width && self.width == self.drive / self.feedback
        {
            1.0
        } else {
            product_over(self.feedback, travel, self.drive)
        }
    }

    fn exposure(&self, travel: f64) -> Result<f64, ThermalTransportError> {
        if travel == 0.0 { return Ok(0.0); }
        let z = self.fraction(travel);
        if z == 1.0 { return Ok(f64::INFINITY); }
        if !z.is_finite() || !(0.0..1.0).contains(&z) {
            return Err(invalid("driven exchange exceeded its positive driving interval"));
        }
        let end_capacity = self.capacity(travel);
        let amount = if self.feedback == 1.0 {
            // Preserve the original fast, cancellation-resistant bath integral.
            let [near, far] = entropy_weights_with_width(self.drive - travel, self.drive, travel);
            end_capacity * near + self.initial_capacity * far
        } else {
            let heat = (0.5 * self.initial_capacity + 0.5 * end_capacity) * travel;
            if self.feedback == 0.0 {
                heat / self.drive
            } else {
                let [a0, a1, _, _] = exchange_moment_weights(z);
                let scale = self.initial_capacity.max(end_capacity);
                let (c0, c1) = (self.initial_capacity / scale, end_capacity / scale);
                let shape = 2.0 * (c0 * a0 + c1 * a1) / (c0 + c1);
                product_over(heat, shape, self.drive)
            }
        };
        if amount.is_nan() || amount < 0.0 {
            return Err(invalid("heat-exchange integral is not representable"));
        }
        Ok(amount)
    }

    fn inverse(&self, duration: f64) -> Result<f64, ThermalTransportError> {
        if self.feedback == 0.0 {
            // Exact affine caloric inverse at constant applied power.
            let mean = 0.5 * self.initial_capacity + 0.5 * self.end_capacity;
            let amount = mean * self.width;
            let f = if amount.is_finite() {
                product_over(self.drive, duration, amount)
            } else {
                product_over(self.drive, duration, mean) / self.width
            };
            let scale = self.initial_capacity.max(self.end_capacity);
            let (a, b) = (self.initial_capacity / scale, self.end_capacity / scale);
            let root = ((1.0 - f) * a * a + f * b * b).sqrt();
            let travel = self.width * (f * (a + b) / (a + root));
            if travel.is_finite() && (0.0..=self.width).contains(&travel) {
                return Ok(travel);
            }
            return Err(invalid("constant-source exchange left its admitted segment"));
        }
        let attenuation = product_over(self.feedback, duration, self.initial_capacity);
        let exponential = if attenuation == 0.0 {
            product_over(self.drive, duration, self.initial_capacity)
        } else {
            product_over(self.drive, -(-attenuation).exp_m1(), self.feedback)
        };
        if self.initial_capacity == self.end_capacity {
            if exponential.is_finite() && (0.0..=self.width).contains(&exponential) {
                return Ok(exponential);
            }
            return Err(invalid("constant-capacity exchange left its admitted segment"));
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
            if value == duration || (value - duration).abs() <= 4.0 * f64::EPSILON * duration {
                return Ok(trial);
            }
            if value < duration {
                lower = trial;
                low_value = value;
            } else {
                upper = trial;
                high_value = value;
            }
            if lower.next_up() >= upper {
                return Ok(
                    if self.fraction(upper) == 1.0 || high_value - duration < duration - low_value {
                        upper
                    } else {
                        lower
                    },
                );
            }
            let remaining = (-self.feedback).mul_add(trial, self.drive);
            let newton = trial + ((duration - value) / self.capacity(trial)) * remaining;
            trial = if iteration % 2 == 0 && newton > lower && newton < upper {
                newton
            } else {
                lower + 0.5 * (upper - lower)
            };
        }
        Err(invalid("heat-exchange inverse exceeded its bounded iteration budget"))
    }

    // The existing positive affine-capacity moments in source-relative drive
    // coordinates. At zero feedback they are the exact cubic time moment of C(T).
    fn exchange_moments(
        &self, travel: f64, duration: f64, heat: f64, gain: f64,
    ) -> (f64, f64, f64, f64) {
        let z = self.fraction(travel);
        let end_capacity = self.capacity(travel);
        if z <= 0.5 {
            let [a0, a1, b0, b1] = exchange_moment_weights(z);
            let scale = self.initial_capacity.max(end_capacity);
            let (c0, c1) = (self.initial_capacity / scale, end_capacity / scale);
            let denominator = c0 * a0 + c1 * a1;
            let w0 = gain * a1;
            let response = duration * ((c0 * w0 + c1 * (1.0 - w0)) / denominator);
            let mean_response = product_over(
                travel, (2.0 * c0 * b0 + c1 * (b1 - b0)) / denominator, self.drive,
            );
            let mean = travel * ((c0 * b0 + c1 * b1) / denominator);
            let mean_gap = travel * ((c0 * (a0 - b0) + c1 * b0) / denominator);
            (mean, response, mean_response, mean_gap)
        } else {
            let w0 = if gain == 0.0 { 0.0 } else { gain * (-gain.ln() - z) / z };
            let scale = self.initial_capacity.max(end_capacity);
            let (c0, c1) = (self.initial_capacity / scale, end_capacity / scale);
            let response = product_over(scale, c0 * w0 + c1 * (z - w0), self.feedback);
            let heat_rate = heat / duration;
            let mean = if heat_rate <= 0.75 * self.drive {
                (self.drive - heat_rate) / self.feedback
            } else {
                let [a0, a1, b0, b1] = exchange_moment_weights(z);
                travel * ((c0 * b0 + c1 * b1) / (c0 * a0 + c1 * a1))
            };
            let mean_gap = (heat_rate - self.drive * gain) / self.feedback;
            (mean, response, (1.0 - response / duration) / self.feedback, mean_gap)
        }
    }

    // Physical fixed-intercept response dE_out/dE_in = driving_out/driving_in.
    // At rounded equilibrium retain the original finite-duration exponential.
    fn inlet_energy_gain(&self, travel: f64, duration: f64) -> Result<f64, ThermalTransportError> {
        let gain = if self.feedback == 0.0 {
            1.0
        } else if self.initial_capacity == self.end_capacity {
            (-product_over(self.feedback, duration, self.initial_capacity)).exp()
        } else if self.fraction(travel) < 1.0 {
            (-self.feedback).mul_add(travel, self.drive) / self.drive
        } else if self.width == self.drive / self.feedback {
            let capacity_change = self.end_capacity - self.initial_capacity;
            let exposure = self.feedback * duration;
            let attenuation = if capacity_change < 0.0 {
                let (head, tail) = fs_math::eft::two_sum(exposure, -self.initial_capacity);
                ((head + self.end_capacity) + tail) / self.end_capacity
            } else {
                exposure / self.end_capacity + capacity_change / self.end_capacity
            };
            let value = (-attenuation).exp();
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

// Nonnegative a*b/c without an avoidable intermediate product or ratio overflow.
fn product_over(a: f64, b: f64, c: f64) -> f64 {
    if a == 0.0 || b == 0.0 { return 0.0; }
    let (large, small) = if a >= b { (a, b) } else { (b, a) };
    let ratio = large / c;
    if ratio.is_finite() { ratio * small } else { (small / c) * large }
}

// Response to an infinitesimal drive around the actual equilibrium temperature.
// In particular, zero feedback has dE/dD=u and dT_mean/dD=u/(2*C), not zero.
fn equilibrium_response(capacity: f64, feedback: f64, duration: f64) -> (f64, f64, f64, f64) {
    if duration == 0.0 { return (1.0, 0.0, 0.0, 0.0); }
    let attenuation = product_over(feedback, duration, capacity);
    let gain = (-attenuation).exp();
    let loss = -(-attenuation).exp_m1();
    let response = if attenuation == 0.0 {
        duration
    } else if attenuation <= 0.5 {
        duration * (loss / attenuation)
    } else {
        product_over(capacity, loss, feedback)
    };
    let mean_response = if attenuation < 0.125 {
        let (mut sum, mut term) = (0.5, 0.5);
        for n in 1..=18 {
            term *= -attenuation / f64::from(n + 2);
            sum += term;
        }
        product_over(duration, sum, capacity)
    } else {
        (1.0 - response / duration) / feedback
    };
    (gain, loss, response, mean_response)
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
    fn g1_exposure_response_keeps_weak_and_saturated_newton_cooling() {
        let knots = [[200.0, 1000.0], [1000.0, 1000.0]];
        let curve = CapacityCurve::new(&knots);
        for (initial, reservoir) in [(300.0, 800.0), (800.0, 300.0)] {
            let difference = reservoir - initial;
            for ntu in [0.0_f64, 1e-12, 1e-5, 0.5, 2.0, 80.0, 1e18] {
                let answer = curve.linear_heat_exchange(initial, reservoir, 1000.0 * ntu).unwrap();
                let expected_remaining = difference * (-ntu).exp();
                // Independent Newton-cooling mean; its weak series avoids
                // cancellation in (1-exp(-ntu))/ntu-exp(-ntu).
                let expected_mean_derivative = difference * if ntu < 1e-3 {
                    ntu * (0.5 - ntu / 3.0 + ntu * ntu / 8.0 - ntu.powi(3) / 30.0)
                } else {
                    -(-ntu).exp_m1() / ntu - (-ntu).exp()
                };
                assert!((answer.remaining_temperature_difference - expected_remaining).abs()
                    <= 32.0 * f64::EPSILON * expected_remaining.abs());
                assert!((answer.mean_log_exposure_derivative - expected_mean_derivative).abs()
                    <= 64.0 * f64::EPSILON * expected_mean_derivative.abs());
                if ntu == 80.0 {
                    assert_eq!(answer.temperature, reservoir);
                    assert_ne!(answer.remaining_temperature_difference, 0.0);
                }
                if ntu == 1e18 {
                    assert_eq!(answer.temperature, answer.mean_temperature);
                    assert_ne!(answer.mean_log_exposure_derivative, 0.0);
                }
            }
        }
        for exposure in [0.0, 1e-9, 1000.0, 80_000.0] {
            let answer = curve.linear_heat_exchange(400.0, 400.0, exposure).unwrap();
            assert_eq!(answer.remaining_temperature_difference, 0.0);
            assert_eq!(answer.mean_log_exposure_derivative, 0.0);
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
                let answer = curve.linear_heat_exchange(initial, reservoir, exposure).unwrap();
                let actual = answer.temperature;
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
                let remaining = reservoir - target;
                let mean_derivative = exact_heat / exposure - remaining;
                assert!((answer.remaining_temperature_difference - remaining).abs()
                    <= 16.0 * f64::EPSILON * remaining.abs());
                assert!((answer.mean_log_exposure_derivative - mean_derivative).abs()
                    <= 32.0 * f64::EPSILON * mean_derivative.abs());
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
            let exact_heat: f64 = cuts.windows(2).map(|pair| {
                let middle = 0.5 * (pair[0] + pair[1]);
                let capacity = if middle < 300.0 {
                    800.0 + 2.0 * middle
                } else if middle < 400.0 {
                    200.0 + 4.0 * middle
                } else {
                    3400.0 - 4.0 * middle
                };
                capacity * (pair[1] - pair[0])
            }).sum();
            let remaining = reservoir - target;
            let mean_derivative = exact_heat / exposure - remaining;
            assert!((answer.remaining_temperature_difference - remaining).abs() < 2e-12);
            assert!((answer.mean_log_exposure_derivative - mean_derivative).abs() < 2e-12);
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
    fn g1_native_exchange_retains_heat_below_absolute_temperature_resolution() {
        let knots = [[300.0, 700.0], [475.0, 875.0], [700.0, 1100.0]];
        let curve = CapacityCurve::new(&knots);
        for (initial, reservoir) in [(300.0, 600.0), (600.0, 300.0)] {
            let exposure = 1e-20;
            assert!(curve.linear_heat_exchange(initial, reservoir, exposure).is_err());
            let result = curve
                .linear_heat_exchange_preserving_heat(initial, reservoir, exposure)
                .unwrap();
            let expected = exposure * (reservoir - initial);
            assert_eq!(result.temperature, initial);
            assert_eq!(result.mean_temperature, initial);
            assert!((result.energy_change - expected).abs() <= 8.0 * f64::EPSILON * expected.abs());
            assert!(result.inlet_energy_loss > 0.0);
            assert!(result.reservoir_energy_derivative > 0.0);
            assert!(result.mean_reservoir_derivative > 0.0);
            assert_eq!(
                curve
                    .temperature_after_energy_change_preserving_energy(initial, result.energy_change)
                    .unwrap(),
                initial
            );
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
        assert_eq!(inverse.temperature(&curve, f64::from_bits(1)).unwrap(), datum);
        assert!(curve.temperature_after_energy_change(datum, f64::from_bits(1)).is_err());

        // A non-dyadic datum exposes cancellation against a distant cached
        // prefix: it must not manufacture warming or cooling from tiny heat.
        let knots = [[200.0, 500.0], [400.0, 4000.0]];
        let curve = CapacityCurve::new(&knots);
        let datum: f64 = 325.3;
        let inverse = CapacityInverse::new(&curve, datum).unwrap();
        for heat in [f64::from_bits(1), -f64::from_bits(1)] {
            assert_eq!(inverse.temperature(&curve, heat).unwrap(), datum);
            assert!(curve.temperature_after_energy_change(datum, heat).is_err());
        }
        for target in [datum.next_down(), datum.next_up()] {
            let heat = curve.energy_change(datum, target).unwrap();
            assert_eq!(inverse.temperature(&curve, heat).unwrap(), target);
        }
    }
    fn close_driven(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:.3e}"
        );
    }

    #[test]
    fn g1_driven_exchange_matches_affine_source_and_allows_negative_virtual_equilibrium() {
        let knots = [[300.0, 350.0], [440.0, 420.0], [700.0, 550.0]];
        let curve = CapacityCurve::new(&knots);
        for (initial, target, virtual_t) in [(360.0, 600.0, 720.0), (650.0, 350.0, -200.0)] {
            // Independent integral for C(T)=200+T/2. The virtual intercept is
            // used only by this oracle; production accepts the relative drive.
            let x = affine_exchange(initial, target, virtual_t, 200.0, 0.5);
            let drive = x * (virtual_t - initial);
            let response = curve.driven_heat_exchange(initial, drive, x).unwrap();
            let heat = (200.0 + 0.25 * (initial + target)) * (target - initial);
            close_driven(response.temperature, target, 2e-10);
            close_driven(response.temperature_change, target - initial, 2e-10);
            close_driven(response.energy_change, heat, 2e-7);
            close_driven(response.mean_temperature_change, (drive - heat) / x, 2e-10);
            close_driven(response.remaining_drive, drive * (virtual_t - target) / (virtual_t - initial), 2e-7);
            close_driven(response.inlet_energy_gain, (virtual_t - target) / (virtual_t - initial), 2e-12);
            close_driven(response.mean_log_duration_derivative,
                response.temperature_change - response.mean_temperature_change, 2e-10);
        }
    }

    #[test]
    fn g1_pure_source_keeps_exact_heat_and_capacity_weighted_cubic_time_mean() {
        let knots = [[300.0, 350.0], [440.0, 420.0], [700.0, 550.0]];
        let curve = CapacityCurve::new(&knots);
        for (initial, target) in [(350.0, 650.0), (650.0, 350.0)] {
            let heat = (200.0 + 0.25 * (initial + target)) * (target - initial);
            let response = curve.driven_heat_exchange(initial, heat, 0.0).unwrap();
            let mean = (100.0 * (initial + target)
                + (initial * initial + initial * target + target * target) / 6.0)
                / (200.0 + 0.25 * (initial + target));
            assert_eq!(response.energy_change, heat);
            assert_eq!(response.inlet_energy_gain, 1.0);
            assert_eq!(response.inlet_energy_loss, 0.0);
            assert_eq!(response.source_gain, 1.0);
            close_driven(response.temperature, target, 2e-10);
            close_driven(response.mean_temperature, mean, 2e-10);
            close_driven(response.mean_source_derivative, (target - mean) / heat, 2e-14);
            assert!((mean - 0.5 * (initial + target)).abs() > 8.0);
        }
        let equilibrium = curve.driven_heat_exchange(500.0, 0.0, 0.0).unwrap();
        assert_eq!(equilibrium.source_gain, 1.0);
        assert_eq!(equilibrium.mean_source_derivative, 1.0 / 900.0);
    }

    #[test]
    fn g1_driven_exchange_keeps_weak_heat_moments_and_equilibrium_source_response() {
        let knots = [[300.0, 350.0], [700.0, 550.0]];
        let curve = CapacityCurve::new(&knots);
        let capacity = 450.0;
        for x in [0.0, 1e-20, 5.0, 200.0, 1e8] {
            let a: f64 = x / capacity;
            let (phi, psi) = if a < 1e-4 {
                (1.0 - a / 2.0 + a * a / 6.0, 0.5 - a / 6.0 + a * a / 24.0)
            } else {
                let phi = -(-a).exp_m1() / a;
                (phi, (1.0 - phi) / a)
            };
            for drive in [-1e-20, 0.0, 1e-20] {
                let response = curve.driven_heat_exchange(500.0, drive, x).unwrap();
                assert_eq!(response.temperature, 500.0);
                assert_eq!(response.mean_temperature, 500.0);
                close_driven(response.energy_change, drive * phi, 3e-14 * (drive * phi).abs());
                close_driven(response.temperature_change, drive * phi / capacity, 3e-14 * (drive * phi / capacity).abs());
                close_driven(response.mean_temperature_change, drive * psi / capacity, 3e-14 * (drive * psi / capacity).abs());
                close_driven(response.source_gain, phi, 3e-14 * phi);
                close_driven(response.mean_source_derivative, psi / capacity, 3e-14 * psi / capacity);
                close_driven(response.inlet_energy_gain, (-a).exp(), 3e-14);
            }
        }
    }

    #[test]
    fn g1_driven_exchange_does_not_erase_feedback_when_virtual_difference_overflows() {
        let knots = [[1.0, 1.0], [1.7e308, 1.0]];
        let curve = CapacityCurve::new(&knots);
        let (drive, x): (f64, f64) = (1e308, 0.1);
        assert!((drive / x).is_infinite());
        let response = curve.driven_heat_exchange(1.0, drive, x).unwrap();
        let phi = -(-x).exp_m1() / x;
        let psi = (1.0 - phi) / x;
        close_driven(response.energy_change / drive, phi, 2e-14);
        close_driven(response.mean_temperature_change / drive, psi, 2e-14);
        close_driven(response.source_gain, phi, 2e-14);
        close_driven(response.mean_source_derivative, psi, 2e-14);
        assert!(response.energy_change < 0.96 * drive);
    }

    #[test]
    fn g3_driven_exchange_retains_physical_inlet_and_source_derivatives_across_knots() {
        let knots = [[300.0, 350.0], [440.0, 420.0], [700.0, 550.0]];
        let curve = CapacityCurve::new(&knots);
        let (initial, drive, x) = (350.0, 160000.0, 150.0);
        let response = curve.driven_heat_exchange(initial, drive, x).unwrap();
        let [lo, hi] = [-0.1, 0.1].map(|d| curve.driven_heat_exchange(initial, drive + d, x).unwrap());
        close_driven((hi.energy_change - lo.energy_change) / 0.2, response.source_gain, 2e-7);
        close_driven((hi.mean_temperature - lo.mean_temperature) / 0.2, response.mean_source_derivative, 2e-9);
        let [lo, hi] = [-0.001, 0.001].map(|d| {
            let shifted = curve.driven_heat_exchange(initial + d, drive - x * d, x).unwrap();
            curve.energy_change(initial, initial + d).unwrap() + shifted.energy_change
        });
        let inlet_change = curve.energy_change(initial - 0.001, initial + 0.001).unwrap();
        close_driven((hi - lo) / inlet_change, response.inlet_energy_gain, 2e-7);
        let [lo, hi] = [-1e-4, 1e-4].map(|d| {
            curve.driven_heat_exchange(initial, drive * (1.0 + d), x * (1.0 + d)).unwrap()
        });
        close_driven((hi.mean_temperature - lo.mean_temperature) / 2e-4,
            response.mean_log_duration_derivative, 2e-6);
        let first = curve.driven_heat_exchange(initial, drive / 2.0, x / 2.0).unwrap();
        let second = curve.driven_heat_exchange(first.temperature, first.remaining_drive, x / 2.0).unwrap();
        close_driven(first.energy_change + second.energy_change, response.energy_change, 2e-7);
        close_driven(0.5 * (first.mean_temperature_change + first.temperature_change + second.mean_temperature_change),
            response.mean_temperature_change, 2e-10);
    }

    #[test]
    fn g4_driven_exchange_refuses_invalid_inputs_and_unsupported_source_trajectories() {
        let knots = [[300.0, 350.0], [440.0, 420.0], [700.0, 550.0]];
        let curve = CapacityCurve::new(&knots);
        for initial in [299.0, 701.0, f64::NAN, f64::INFINITY] {
            assert!(curve.driven_heat_exchange(initial, 0.0, 0.0).is_err());
        }
        for drive in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(curve.driven_heat_exchange(500.0, drive, 1.0).is_err());
        }
        for x in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(curve.driven_heat_exchange(500.0, 0.0, x).is_err());
        }
        for x in [0.0, 100.0] {
            assert!(curve.driven_heat_exchange(500.0, 1e9, x).is_err());
            assert!(curve.driven_heat_exchange(500.0, -1e9, x).is_err());
        }
    }

}
