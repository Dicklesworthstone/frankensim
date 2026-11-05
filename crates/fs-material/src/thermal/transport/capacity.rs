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
