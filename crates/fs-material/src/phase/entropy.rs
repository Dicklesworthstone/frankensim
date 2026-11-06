//! Relative thermal entropy on the original equilibrium enthalpy chart.

use super::{EquilibriumEnthalpyPhaseCurve, PhaseStateError};

mod caloric;

/// Finite thermal conjugate on one original equilibrium enthalpy chart.
///
/// The temperature is `delta_h / delta_s`, with the continuous equal-state
/// limit. It supplies a discrete-gradient heat effort, not a transport law or
/// an absolute entropy datum. Pressure, composition and work conventions stay
/// with the original chart and its coupled physical owner.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EnthalpyEntropyMean {
    /// Entropy-conjugate temperature [K], between the two endpoint temperatures.
    pub temperature_k: f64,
    /// Signed original-chart integral of `dh / T(h)` [J/(kg K)].
    pub specific_entropy_change_j_kg_k: f64,
    /// Partial derivative of the mean with respect to initial h [kg K/J].
    pub d_temperature_d_initial_enthalpy_kg_k_j: f64,
    /// Partial derivative of the mean with respect to final h [kg K/J].
    pub d_temperature_d_final_enthalpy_kg_k_j: f64,
    /// Original intervals visited; zero for equal enthalpies.
    pub intervals: usize,
}

/// Original affine Cp(T)/T integral for temperature-coordinate consumers.
/// The caller supplies capacities from the same admitted sensible source
/// interval; no latent plateau, extrapolation or capacity fitting occurs.
pub(crate) fn affine_cp_entropy(
    lower_temperature_k: f64,
    upper_temperature_k: f64,
    lower_capacity_j_kg_k: f64,
    upper_capacity_j_kg_k: f64,
) -> Result<f64, PhaseStateError> {
    caloric::affine_cp_entropy(
        lower_temperature_k,
        upper_temperature_k,
        [lower_capacity_j_kg_k, upper_capacity_j_kg_k],
    )
    .ok_or(PhaseStateError::InvalidCurve {
        what: "original affine heat-capacity entropy is not finite and representable",
    })
}

impl EquilibriumEnthalpyPhaseCurve {
    fn caloric_entropy_segment(
        &self,
        index: usize,
        lower_enthalpy: f64,
        upper_enthalpy: f64,
        lower_temperature: f64,
    ) -> Result<Option<caloric::CaloricEntropySegment>, PhaseStateError> {
        let Some(capacities) = self.sensible_heat_capacities_on_interval(index) else {
            return Ok(None);
        };
        let pair = &self.knots()[index..=index + 1];
        // Use the original quadratic inverse before its temperature is
        // rounded. This retains the local capacities for short enthalpy steps
        // and uses the same nearer-endpoint datum convention as T(h), without
        // adding binary searches or scanning unvisited knots.
        let capacity_at = |h| {
            super::caloric::temperature_fraction(pair[0], pair[1], capacities, h)
                .map(|fraction| super::caloric::capacity_at_fraction(capacities, fraction))
        };
        let local_capacities = [capacity_at(lower_enthalpy)?, capacity_at(upper_enthalpy)?];
        caloric::segment(
            upper_enthalpy - lower_enthalpy,
            lower_temperature,
            local_capacities,
        )
        .map(Some)
        .ok_or(PhaseStateError::InvalidCurve {
            what: "original affine heat-capacity entropy is not finite and representable",
        })
    }

    /// Evaluate the exact finite entropy-conjugate temperature and its Jacobian.
    ///
    /// For unequal enthalpies, `theta = (h1-h0) / integral_h0^h1 dh/T(h)`.
    /// The same logarithmic sensible integrals and exact latent plateaus as
    /// [`Self::specific_entropy_change_j_kg_k`] are used. Reversing endpoints
    /// preserves theta, reverses entropy, and swaps the two derivatives.
    /// No apparent heat capacity, resampling, extrapolation or phase-state
    /// hashing is introduced. All fixed-pressure/composition qualifications of
    /// that entropy integral still apply; a transport owner must independently
    /// account for heat direction, boundary entropy and any additional work.
    ///
    /// At equal enthalpies, theta is the original T(h), entropy is zero and
    /// each derivative is half the original dT/dh. At a chart corner this
    /// selects the existing right derivative (left at the upper endpoint),
    /// not a claim of a unique two-sided derivative. A latent interval has
    /// exactly its transition temperature and zero temperature derivatives.
    ///
    /// The endpoint derivatives integrate the positive kernels
    /// `theta^2 * (1-u) T'(h(u))/T(h(u))^2` and
    /// `theta^2 * u T'(h(u))/T(h(u))^2`, where
    /// `h(u) = (1-u) h0 + u h1`. This avoids subtracting nearly equal entropy
    /// quotients for small sensible steps, including steps whose temperatures
    /// round to the same floating-point value.
    ///
    /// Work is O(log N + K), with no allocation. Only K crossed original
    /// intervals are visited, once each; K must not exceed `max_intervals`.
    /// A caller bounding an entire solve passes its remaining interval budget
    /// and subtracts the returned count. Cancellation is checked before and
    /// after evaluation and every 64 intervals. The original endpoint density
    /// and domain admission is retained. Unrepresentable required arithmetic
    /// refuses without a candidate or a hidden temperature/derivative clamp.
    #[allow(clippy::too_many_lines)]
    pub fn thermodynamic_mean_temperature(
        &self,
        initial_enthalpy_j_kg: f64,
        final_enthalpy_j_kg: f64,
        max_intervals: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<EnthalpyEntropyMean, PhaseStateError> {
        if cancelled() {
            return Err(PhaseStateError::EntropyEvaluationCancelled);
        }
        let initial = self.values_at_specific_enthalpy(initial_enthalpy_j_kg)?;
        let final_state = self.values_at_specific_enthalpy(final_enthalpy_j_kg)?;
        let invalid = || PhaseStateError::InvalidCurve {
            what: "entropy-conjugate temperature or derivative is not finite and representable",
        };
        if initial_enthalpy_j_kg == final_enthalpy_j_kg {
            let slope = self.temperature_derivative_at_specific_enthalpy(initial_enthalpy_j_kg)?;
            let derivative = 0.5 * slope;
            if slope > 0.0 && derivative == 0.0 {
                return Err(invalid());
            }
            if cancelled() {
                return Err(PhaseStateError::EntropyEvaluationCancelled);
            }
            return Ok(EnthalpyEntropyMean {
                temperature_k: initial.temperature_k,
                specific_entropy_change_j_kg_k: 0.0,
                d_temperature_d_initial_enthalpy_kg_k_j: derivative,
                d_temperature_d_final_enthalpy_kg_k_j: derivative,
                intervals: 0,
            });
        }
        let forward = final_enthalpy_j_kg > initial_enthalpy_j_kg;
        let (lower, upper) = if forward {
            (initial, final_state)
        } else {
            (final_state, initial)
        };
        let lo = lower.specific_enthalpy_j_kg;
        let hi = upper.specific_enthalpy_j_kg;
        let width = hi - lo;
        if !width.is_finite() || width <= 0.0 {
            return Err(invalid());
        }
        let knots = self.knots();
        let first = knots.partition_point(|knot| knot.specific_enthalpy_j_kg <= lo) - 1;
        let last = knots.partition_point(|knot| knot.specific_enthalpy_j_kg < hi);
        let intervals = last - first;
        if intervals > max_intervals {
            return Err(PhaseStateError::InvalidCurve {
                what: "thermal-entropy original-interval budget exhausted",
            });
        }
        // Entropy, integral (T-T_lo)/T dh, and both derivative kernels scaled
        // by T_lo^2. Scaling avoids squaring absolute temperatures or tiny dh.
        let mut totals = [0.0; 4];
        let mut corrections = [0.0; 4];
        let mut varying_temperature = false;
        for (index, pair) in knots[first..=last].windows(2).enumerate() {
            if index.is_multiple_of(64) && cancelled() {
                return Err(PhaseStateError::EntropyEvaluationCancelled);
            }
            let a = lo.max(pair[0].specific_enthalpy_j_kg);
            let b = hi.min(pair[1].specific_enthalpy_j_kg);
            let ta = if a == lo {
                lower.temperature_k
            } else {
                pair[0].temperature_k
            };
            let tb = if b == hi {
                upper.temperature_k
            } else {
                pair[1].temperature_k
            };
            let dh = b - a;
            let (entropy, loss, zeroth_moment, lower_moment, upper_moment, slope) =
                if let Some(source) = self.caloric_entropy_segment(first + index, a, b, ta)? {
                    (
                        source.entropy,
                        source.loss,
                        source.zeroth_kernel,
                        source.lower_kernel,
                        source.upper_kernel,
                        source.reciprocal_mean_capacity,
                    )
                } else {
                    let (entropy, x) = segment_entropy(dh, ta, tb)?;
                    let (loss, first_moment) = logarithmic_moments(x);
                    let zeroth_moment = 1.0 / (1.0 + x);
                    let temperature_increment = pair[1].temperature_k - pair[0].temperature_k;
                    let slope = temperature_increment
                        / (pair[1].specific_enthalpy_j_kg - pair[0].specific_enthalpy_j_kg);
                    if !slope.is_finite() || (temperature_increment > 0.0 && slope == 0.0) {
                        return Err(PhaseStateError::UnrepresentableTemperatureDerivative);
                    }
                    (
                        entropy,
                        loss,
                        zeroth_moment,
                        zeroth_moment - first_moment,
                        first_moment,
                        slope,
                    )
                };
            varying_temperature |= slope > 0.0;
            let weight = dh / width;
            let relative_temperature = lower.temperature_k / ta;
            let scaled_slope = (slope * relative_temperature) * relative_temperature;
            if weight <= 0.0
                || lower_moment <= 0.0
                || upper_moment <= 0.0
                || (slope > 0.0 && scaled_slope == 0.0)
            {
                return Err(invalid());
            }
            let before = (a - lo) / width;
            let after = (hi - b) / width;
            let values = [
                entropy,
                dh * ((ta - lower.temperature_k) / ta + relative_temperature * loss),
                weight * scaled_slope * (after * zeroth_moment + weight * lower_moment),
                weight * scaled_slope * (before * zeroth_moment + weight * upper_moment),
            ];
            for component in 0..4 {
                if !values[component].is_finite() || values[component] < 0.0 {
                    return Err(invalid());
                }
                let increment = values[component] - corrections[component];
                let next = totals[component] + increment;
                corrections[component] = (next - totals[component]) - increment;
                totals[component] = next;
            }
        }
        let entropy = totals[0];
        // Equivalent to width/entropy, but a plateau is exactly T_lo and a
        // near-isothermal sensible interval retains its small temperature rise.
        let temperature = lower.temperature_k + totals[1] / entropy;
        let scale = temperature / lower.temperature_k;
        let lower_derivative = (totals[2] * scale) * scale;
        let upper_derivative = (totals[3] * scale) * scale;
        if !entropy.is_finite()
            || entropy <= 0.0
            || !temperature.is_finite()
            || temperature < lower.temperature_k
            || temperature > upper.temperature_k
            || !scale.is_finite()
            || !lower_derivative.is_finite()
            || !upper_derivative.is_finite()
            || lower_derivative < 0.0
            || upper_derivative < 0.0
            || (varying_temperature && (lower_derivative == 0.0 || upper_derivative == 0.0))
        {
            return Err(invalid());
        }
        if cancelled() {
            return Err(PhaseStateError::EntropyEvaluationCancelled);
        }
        let (initial_derivative, final_derivative) = if forward {
            (lower_derivative, upper_derivative)
        } else {
            (upper_derivative, lower_derivative)
        };
        Ok(EnthalpyEntropyMean {
            temperature_k: temperature,
            specific_entropy_change_j_kg_k: if forward { entropy } else { -entropy },
            d_temperature_d_initial_enthalpy_kg_k_j: initial_derivative,
            d_temperature_d_final_enthalpy_kg_k_j: final_derivative,
            intervals,
        })
    }

    /// Integrate `ds = dh / T(h)` between two admitted specific enthalpies.
    ///
    /// Each nonisothermal segment uses its exact logarithmic integral; a latent
    /// plateau contributes `delta_h / T_transition` without an apparent heat
    /// capacity or separate latent-energy credit. Sensible intervals generated
    /// from heat-capacity data integrate the original affine Cp(T)/T exactly;
    /// explicitly authored enthalpy charts retain their affine T(h) integral.
    /// Reversing the endpoints reverses the sign. No extrapolation is allowed.
    ///
    /// This is a relative equilibrium thermal entropy [J/(kg K)] at the chart's
    /// fixed pressure/composition convention. Absolute entropy, chemical mixing,
    /// pressure changes, hysteresis and additional work coordinates are not
    /// inferred from caloric data. A coupled owner must account for those
    /// separately and justify its enthalpy/internal-energy convention.
    ///
    /// Only crossed chart segments are visited: O(log N + K) work for N knots
    /// and K crossed segments, without allocation or phase-state hashing.
    /// Endpoint evaluation retains the original density/domain admission.
    ///
    /// Returns the original state/domain errors, or `InvalidCurve` when the
    /// integral cannot be represented. No entropy-production claim is made by
    /// this state function; the transport owner must audit its boundary entropy.
    pub fn specific_entropy_change_j_kg_k(
        &self,
        initial_enthalpy_j_kg: f64,
        final_enthalpy_j_kg: f64,
    ) -> Result<f64, PhaseStateError> {
        let initial = self.values_at_specific_enthalpy(initial_enthalpy_j_kg)?;
        let final_state = self.values_at_specific_enthalpy(final_enthalpy_j_kg)?;
        if initial_enthalpy_j_kg == final_enthalpy_j_kg {
            return Ok(0.0);
        }
        let (lower, upper) = if initial_enthalpy_j_kg < final_enthalpy_j_kg {
            (initial, final_state)
        } else {
            (final_state, initial)
        };
        let lo = lower.specific_enthalpy_j_kg;
        let hi = upper.specific_enthalpy_j_kg;
        let knots = self.knots();
        // At a knot, the lower endpoint uses the segment to its right and the
        // upper endpoint the segment to its left. The validated nonempty span
        // guarantees both indices exist, including at the chart boundaries.
        let first = knots.partition_point(|knot| knot.specific_enthalpy_j_kg <= lo) - 1;
        let last = knots.partition_point(|knot| knot.specific_enthalpy_j_kg < hi);
        let invalid = || PhaseStateError::InvalidCurve {
            what: "relative thermal entropy is not finite and representable",
        };
        let mut total = 0.0;
        let mut correction = 0.0;
        for (index, pair) in knots[first..=last].windows(2).enumerate() {
            let a = lo.max(pair[0].specific_enthalpy_j_kg);
            let b = hi.min(pair[1].specific_enthalpy_j_kg);
            // Only the two caller endpoints can lie inside a segment. Every
            // interior boundary is an original knot, so its admitted values
            // need neither another search nor a reconstructed state identity.
            let ta = if a == lo {
                lower.temperature_k
            } else {
                pair[0].temperature_k
            };
            let tb = if b == hi {
                upper.temperature_k
            } else {
                pair[1].temperature_k
            };
            let entropy =
                if let Some(source) = self.caloric_entropy_segment(first + index, a, b, ta)? {
                    source.entropy
                } else {
                    segment_entropy(b - a, ta, tb)?.0
                };
            // Fixed-order compensated summation preserves small segment terms.
            let increment = entropy - correction;
            let next = total + increment;
            correction = (next - total) - increment;
            total = next;
        }
        if !total.is_finite() || total <= 0.0 {
            return Err(invalid());
        }
        Ok(if final_enthalpy_j_kg < initial_enthalpy_j_kg {
            -total
        } else {
            total
        })
    }
}

// Keep the original entropy arithmetic shared by the scalar query and the
// finite conjugate. x is nonnegative on the admitted monotone caloric chart.
fn segment_entropy(dh: f64, ta: f64, tb: f64) -> Result<(f64, f64), PhaseStateError> {
    let x = (tb - ta) / ta;
    let ratio = if x == 0.0 { 1.0 } else { x.ln_1p() / x };
    let entropy = (dh / ta) * ratio;
    if !x.is_finite() || !entropy.is_finite() || entropy <= 0.0 {
        return Err(PhaseStateError::InvalidCurve {
            what: "relative thermal entropy is not finite and representable",
        });
    }
    Ok((entropy, x))
}

// Return 1-ln(1+x)/x and integral_0^1 u/(1+x*u)^2 du. Direct subtraction
// loses small sensible-temperature changes; these convergent series retain
// them. Twenty terms at x<=1/8 have truncation error below 1e-18.
fn logarithmic_moments(x: f64) -> (f64, f64) {
    if x <= 0.125 {
        let mut power = 1.0;
        let mut loss_over_x = 0.0;
        let mut first = 0.0;
        for n in 0..20 {
            let denominator = f64::from(n + 2);
            loss_over_x += power / denominator;
            first += power * f64::from(n + 1) / denominator;
            power *= -x;
        }
        (x * loss_over_x, first)
    } else {
        let logarithm = x.ln_1p();
        (1.0 - logarithm / x, ((logarithm - x / (1.0 + x)) / x) / x)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::phase::{
        EnthalpyPhaseKnot, HeatCapacityKnot, MAX_HEAT_CAPACITY_ENTHALPY_KNOTS, SolidLiquidPhase,
    };
    use fs_blake3::ContentHash;

    fn phase_curve() -> EquilibriumEnthalpyPhaseCurve {
        let knot = |h, t, f| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: t,
            liquid_mass_fraction: f,
            bulk_density_kg_m3: 1000.0,
        };
        EquilibriumEnthalpyPhaseCurve::try_new(
            ContentHash([0x74; 32]),
            vec![
                knot(0.0, 300.0, 0.0),
                knot(30_000.0, 600.0, 0.0),
                knot(55_000.0, 600.0, 1.0),
                knot(95_000.0, 800.0, 1.0),
            ],
        )
        .unwrap()
    }

    #[test]
    fn original_affine_capacity_controls_entropy_mean_and_analytic_jacobian() {
        for capacities in [[1000.0, 4000.0], [4000.0, 1000.0]] {
            let source = [300.0, 900.0].map(|temperature_k| HeatCapacityKnot {
                temperature_k,
                specific_heat_capacity_j_kg_k: if temperature_k == 300.0 {
                    capacities[0]
                } else {
                    capacities[1]
                },
                bulk_density_kg_m3: 1000.0,
            });
            // Deliberately retain one wide original interval: chord entropy
            // or a frozen Cp cannot satisfy this independent source oracle.
            let curve = EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity(
                ContentHash([0x7c; 32]),
                SolidLiquidPhase::Solid,
                0.0,
                &source,
                1e9,
            )
            .unwrap();
            assert_eq!(curve.knots().len(), 2);
            let slope = (capacities[1] - capacities[0]) / 600.0;
            let enthalpy = |t: f64| {
                let dt = t - 300.0;
                capacities[0] * dt + 0.5 * slope * dt * dt
            };
            let (ta, tb) = (400.0_f64, 750.0_f64);
            let (ha, hb) = (enthalpy(ta), enthalpy(tb));
            let entropy = (capacities[0] - 300.0 * slope) * (tb / ta).ln() + slope * (tb - ta);
            let temperature = (hb - ha) / entropy;
            let mean = curve
                .thermodynamic_mean_temperature(ha, hb, 1, || false)
                .unwrap();
            assert!((mean.specific_entropy_change_j_kg_k / entropy - 1.0).abs() < 1e-14);
            assert!((mean.temperature_k / temperature - 1.0).abs() < 1e-14);
            let lower = (temperature / ta - 1.0) / entropy;
            let upper = (1.0 - temperature / tb) / entropy;
            assert!((mean.d_temperature_d_initial_enthalpy_kg_k_j / lower - 1.0).abs() < 1e-13);
            assert!((mean.d_temperature_d_final_enthalpy_kg_k_j / upper - 1.0).abs() < 1e-13);
            assert_eq!(
                mean.specific_entropy_change_j_kg_k,
                curve.specific_entropy_change_j_kg_k(ha, hb).unwrap()
            );
            let reverse = curve
                .thermodynamic_mean_temperature(hb, ha, 1, || false)
                .unwrap();
            assert_eq!(reverse.temperature_k, mean.temperature_k);
            assert_eq!(
                reverse.d_temperature_d_initial_enthalpy_kg_k_j,
                mean.d_temperature_d_final_enthalpy_kg_k_j
            );
            assert_eq!(
                reverse.d_temperature_d_final_enthalpy_kg_k_j,
                mean.d_temperature_d_initial_enthalpy_kg_k_j
            );
        }
    }

    #[test]
    fn sourced_entropy_and_jacobians_cross_a_sharp_latent_plateau_without_smoothing() {
        let mut curve = phase_curve();
        curve
            .retain_sensible_heat_capacities(vec![Some([50.0, 150.0]), None, Some([100.0, 300.0])])
            .unwrap();
        let expected = -50.0 * 2.0_f64.ln() + 100.0 + 25_000.0 / 600.0
            - 500.0 * (800.0_f64 / 600.0).ln()
            + 200.0;
        let mean = curve
            .thermodynamic_mean_temperature(0.0, 95_000.0, 3, || false)
            .unwrap();
        assert!((mean.specific_entropy_change_j_kg_k - expected).abs() < 1e-12);
        let theta = 95_000.0 / expected;
        assert!((mean.temperature_k / theta - 1.0).abs() < 1e-14);
        assert!(
            (mean.d_temperature_d_initial_enthalpy_kg_k_j - (theta / 300.0 - 1.0) / expected).abs()
                < 1e-15
        );
        assert!(
            (mean.d_temperature_d_final_enthalpy_kg_k_j - (1.0 - theta / 800.0) / expected).abs()
                < 1e-15
        );
        let split = curve.specific_entropy_change_j_kg_k(0.0, 14_000.0).unwrap()
            + curve
                .specific_entropy_change_j_kg_k(14_000.0, 48_000.0)
                .unwrap()
            + curve
                .specific_entropy_change_j_kg_k(48_000.0, 95_000.0)
                .unwrap();
        assert!((split - expected).abs() < 1e-12);
        let plateau = curve
            .thermodynamic_mean_temperature(40_000.0, 50_000.0, 1, || false)
            .unwrap();
        assert_eq!(plateau.temperature_k, 600.0);
        assert_eq!(plateau.d_temperature_d_initial_enthalpy_kg_k_j, 0.0);
        assert_eq!(plateau.d_temperature_d_final_enthalpy_kg_k_j, 0.0);
    }

    #[test]
    fn source_capacity_retains_tiny_enthalpy_entropy_and_thermal_jacobians() {
        let curve = EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity(
            ContentHash([0x7d; 32]),
            SolidLiquidPhase::Solid,
            0.0,
            &[
                HeatCapacityKnot {
                    temperature_k: 300.0,
                    specific_heat_capacity_j_kg_k: 1000.0,
                    bulk_density_kg_m3: 1000.0,
                },
                HeatCapacityKnot {
                    temperature_k: 900.0,
                    specific_heat_capacity_j_kg_k: 4000.0,
                    bulk_density_kg_m3: 1000.0,
                },
            ],
            1e9,
        )
        .unwrap();
        let mean = curve
            .thermodynamic_mean_temperature(0.0, 1e-300, 1, || false)
            .unwrap();
        assert_eq!(mean.temperature_k, 300.0);
        assert_eq!(mean.specific_entropy_change_j_kg_k, 1e-300 / 300.0);
        assert_eq!(mean.d_temperature_d_initial_enthalpy_kg_k_j, 0.0005);
        assert_eq!(mean.d_temperature_d_final_enthalpy_kg_k_j, 0.0005);
    }

    #[test]
    fn caloric_and_latent_entropy_share_one_chart_and_are_path_additive() {
        let curve = phase_curve();
        let entropy = |a, b| curve.specific_entropy_change_j_kg_k(a, b).unwrap();
        let expected = 100.0 * 2.0_f64.ln() + 25_000.0 / 600.0 + 200.0 * (800.0_f64 / 600.0).ln();
        assert!((entropy(0.0, 95_000.0) - expected).abs() < 1e-12);
        assert_eq!(entropy(40_000.0, 50_000.0), 10_000.0 / 600.0);
        assert_eq!(entropy(95_000.0, 0.0), -entropy(0.0, 95_000.0));
        let split =
            entropy(0.0, 20_000.0) + entropy(20_000.0, 48_000.0) + entropy(48_000.0, 95_000.0);
        assert!((split - expected).abs() < 1e-12);
        assert_eq!(entropy(30_000.0, 30_000.0), 0.0);
        assert!(curve.specific_entropy_change_j_kg_k(-1.0, 0.0).is_err());
        assert!(curve.specific_entropy_change_j_kg_k(0.0, 95_001.0).is_err());
        assert!(curve.specific_entropy_change_j_kg_k(f64::NAN, 0.0).is_err());
    }

    #[test]
    fn nearly_isothermal_sensible_entropy_keeps_the_small_temperature_increment() {
        let ta = 300.0;
        let tb = 300.0 + 1e-10;
        let knot = |h, t| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: t,
            liquid_mass_fraction: 0.0,
            bulk_density_kg_m3: 1000.0,
        };
        let curve = EquilibriumEnthalpyPhaseCurve::try_single_phase(
            ContentHash([0x75; 32]),
            SolidLiquidPhase::Solid,
            vec![knot(0.0, ta), knot(1.0, tb)],
        )
        .unwrap();
        let entropy = curve.specific_entropy_change_j_kg_k(0.0, 1.0).unwrap();
        assert!(entropy >= 1.0 / tb && entropy <= 1.0 / ta);
        assert!((entropy * 0.5 * (ta + tb) - 1.0).abs() < 5e-16);
    }

    #[test]
    fn sensible_and_latent_subintervals_match_an_independent_entropy_primitive() {
        let curve = phase_curve();
        let primitive = |h: f64| {
            if h <= 30_000.0 {
                100.0 * (h / 30_000.0).ln_1p()
            } else if h <= 55_000.0 {
                100.0 * 2.0_f64.ln() + (h - 30_000.0) / 600.0
            } else {
                100.0 * 2.0_f64.ln()
                    + 25_000.0 / 600.0
                    + 200.0 * ((h - 55_000.0) / 120_000.0).ln_1p()
            }
        };
        // Both chart endpoints, both sides of the latent interval, and genuine
        // interior points exercise the lower-right/upper-left range convention.
        let endpoints = [
            0.0, 1000.0, 20_000.0, 30_000.0, 40_000.0, 55_000.0, 77_000.0, 95_000.0,
        ];
        for from in endpoints {
            for to in endpoints {
                let actual = curve.specific_entropy_change_j_kg_k(from, to).unwrap();
                let expected = primitive(to) - primitive(from);
                assert!(
                    (actual - expected).abs() < 1e-12,
                    "{from} -> {to}: {actual} != {expected}"
                );
            }
        }
    }

    #[test]
    fn short_steps_inside_a_large_sourced_capacity_chart_retain_caloric_entropy() {
        let capacity = 1024.0;
        let reference = 1.0e12;
        let source: Vec<_> = (0..MAX_HEAT_CAPACITY_ENTHALPY_KNOTS)
            .map(|index| HeatCapacityKnot {
                temperature_k: 250.0 + index as f64 * 0.125,
                specific_heat_capacity_j_kg_k: capacity,
                bulk_density_kg_m3: 1000.0,
            })
            .collect();
        let curve = EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity(
            ContentHash([0x76; 32]),
            SolidLiquidPhase::Solid,
            reference,
            &source,
            1.0,
        )
        .unwrap();
        assert_eq!(curve.knots().len(), MAX_HEAT_CAPACITY_ENTHALPY_KNOTS);
        for index in [
            0,
            MAX_HEAT_CAPACITY_ENTHALPY_KNOTS / 2,
            MAX_HEAT_CAPACITY_ENTHALPY_KNOTS - 2,
        ] {
            let from = reference + index as f64 * 128.0 + 32.0;
            let to = from + 64.0;
            let temperature = source[index].temperature_k + 32.0 / capacity;
            let expected = capacity * ((64.0 / capacity) / temperature).ln_1p();
            let actual = curve.specific_entropy_change_j_kg_k(from, to).unwrap();
            assert!((actual - expected).abs() <= 4.0 * f64::EPSILON * expected);
            assert_eq!(
                curve.specific_entropy_change_j_kg_k(to, from).unwrap(),
                -actual
            );
            assert_eq!(
                curve.specific_entropy_change_j_kg_k(from, from).unwrap(),
                0.0
            );
        }
    }

    #[test]
    fn entropy_endpoints_retain_state_density_and_domain_refusals() {
        let curve = EquilibriumEnthalpyPhaseCurve::try_single_phase(
            ContentHash([0x77; 32]),
            SolidLiquidPhase::Solid,
            vec![
                EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: 0.0,
                    temperature_k: 300.0,
                    liquid_mass_fraction: 0.0,
                    bulk_density_kg_m3: 1e-320,
                },
                EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: 1000.0,
                    temperature_k: 400.0,
                    liquid_mass_fraction: 0.0,
                    bulk_density_kg_m3: 1000.0,
                },
            ],
        )
        .unwrap();
        // Original source knots are admitted; an interior density cannot be
        // represented because the specific volume of the first knot overflows.
        for h in [-1.0, 500.0, 1001.0, f64::INFINITY] {
            let original = curve.state_at_specific_enthalpy(h).unwrap_err();
            assert_eq!(
                curve.specific_entropy_change_j_kg_k(h, 1000.0).unwrap_err(),
                original
            );
            assert_eq!(
                curve.specific_entropy_change_j_kg_k(0.0, h).unwrap_err(),
                original
            );
            assert_eq!(
                curve.specific_entropy_change_j_kg_k(h, h).unwrap_err(),
                original
            );
            assert_eq!(
                curve
                    .thermodynamic_mean_temperature(h, 1000.0, 1, || false)
                    .unwrap_err(),
                original
            );
            assert_eq!(
                curve
                    .thermodynamic_mean_temperature(0.0, h, 1, || false)
                    .unwrap_err(),
                original
            );
            assert_eq!(
                curve
                    .thermodynamic_mean_temperature(h, h, 0, || false)
                    .unwrap_err(),
                original
            );
        }
    }

    #[test]
    fn finite_thermal_conjugate_spans_sensible_and_latent_storage() {
        let curve = phase_curve();
        let mean = curve
            .thermodynamic_mean_temperature(0.0, 95_000.0, 3, || false)
            .unwrap();
        let entropy = 100.0 * 2.0_f64.ln() + 25_000.0 / 600.0 + 200.0 * (800.0_f64 / 600.0).ln();
        let temperature = 95_000.0 / entropy;
        let initial_derivative = (temperature / 300.0 - 1.0) / entropy;
        let final_derivative = (1.0 - temperature / 800.0) / entropy;
        assert!((mean.temperature_k - temperature).abs() <= 2.0e-13);
        assert!((mean.d_temperature_d_initial_enthalpy_kg_k_j - initial_derivative).abs() < 1e-16);
        assert!((mean.d_temperature_d_final_enthalpy_kg_k_j - final_derivative).abs() < 1e-16);
        assert_eq!(mean.intervals, 3);
        assert_eq!(
            mean.specific_entropy_change_j_kg_k,
            curve.specific_entropy_change_j_kg_k(0.0, 95_000.0).unwrap()
        );
        assert!(
            (mean.temperature_k * mean.specific_entropy_change_j_kg_k - 95_000.0).abs() < 3e-11
        );

        let reverse = curve
            .thermodynamic_mean_temperature(95_000.0, 0.0, 3, || false)
            .unwrap();
        assert_eq!(reverse.temperature_k, mean.temperature_k);
        assert_eq!(
            reverse.specific_entropy_change_j_kg_k,
            -mean.specific_entropy_change_j_kg_k
        );
        assert_eq!(
            reverse.d_temperature_d_initial_enthalpy_kg_k_j,
            mean.d_temperature_d_final_enthalpy_kg_k_j
        );
        assert_eq!(
            reverse.d_temperature_d_final_enthalpy_kg_k_j,
            mean.d_temperature_d_initial_enthalpy_kg_k_j
        );
    }

    #[test]
    fn latent_and_equal_state_means_retain_original_chart_jacobians() {
        let curve = phase_curve();
        for (from, to) in [(30_000.0, 55_000.0), (50_000.0, 40_000.0)] {
            let mean = curve
                .thermodynamic_mean_temperature(from, to, 1, || false)
                .unwrap();
            assert_eq!(mean.temperature_k, 600.0);
            assert_eq!(mean.specific_entropy_change_j_kg_k, (to - from) / 600.0);
            assert_eq!(mean.d_temperature_d_initial_enthalpy_kg_k_j, 0.0);
            assert_eq!(mean.d_temperature_d_final_enthalpy_kg_k_j, 0.0);
            assert_eq!(mean.intervals, 1);
        }
        // At the plateau's start/end use the chart's existing right derivative;
        // the final domain endpoint uses its existing left derivative.
        for (h, temperature, derivative) in [
            (10_000.0, 400.0, 0.005),
            (30_000.0, 600.0, 0.0),
            (55_000.0, 600.0, 0.0025),
            (95_000.0, 800.0, 0.0025),
        ] {
            let mean = curve
                .thermodynamic_mean_temperature(h, h, 0, || false)
                .unwrap();
            assert_eq!(mean.temperature_k, temperature);
            assert_eq!(mean.specific_entropy_change_j_kg_k, 0.0);
            assert_eq!(mean.d_temperature_d_initial_enthalpy_kg_k_j, derivative);
            assert_eq!(mean.d_temperature_d_final_enthalpy_kg_k_j, derivative);
            assert_eq!(mean.intervals, 0);
        }
    }

    #[test]
    fn tiny_sensible_steps_keep_their_nonzero_thermal_jacobian() {
        let knot = |h, t| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: t,
            liquid_mass_fraction: 0.0,
            bulk_density_kg_m3: 1000.0,
        };
        let curve = EquilibriumEnthalpyPhaseCurve::try_single_phase(
            ContentHash([0x78; 32]),
            SolidLiquidPhase::Solid,
            vec![knot(0.0, 300.0), knot(300.0, 600.0)],
        )
        .unwrap();
        // Both endpoint temperatures round to 300 K; the original sensible
        // chart still has dT/dh=1, rather than a latent plateau.
        let mean = curve
            .thermodynamic_mean_temperature(0.0, 1e-300, 1, || false)
            .unwrap();
        assert_eq!(mean.temperature_k, 300.0);
        assert_eq!(mean.specific_entropy_change_j_kg_k, 1e-300 / 300.0);
        assert_eq!(mean.d_temperature_d_initial_enthalpy_kg_k_j, 0.5);
        assert_eq!(mean.d_temperature_d_final_enthalpy_kg_k_j, 0.5);

        let ta = 300.0;
        let tb = ta + 1e-10;
        let curve = EquilibriumEnthalpyPhaseCurve::try_single_phase(
            ContentHash([0x79; 32]),
            SolidLiquidPhase::Solid,
            vec![knot(0.0, ta), knot(1.0, tb)],
        )
        .unwrap();
        let mean = curve
            .thermodynamic_mean_temperature(0.0, 1.0, 1, || false)
            .unwrap();
        let half_slope = 0.5 * (tb - ta);
        assert!(mean.temperature_k >= ta && mean.temperature_k <= tb);
        assert!((mean.temperature_k - 0.5 * (ta + tb)).abs() <= f64::EPSILON * ta);
        assert!((mean.d_temperature_d_initial_enthalpy_kg_k_j / half_slope - 1.0).abs() < 1e-12);
        assert!((mean.d_temperature_d_final_enthalpy_kg_k_j / half_slope - 1.0).abs() < 1e-12);
    }

    #[test]
    fn thermal_conjugates_preserve_the_energy_datum_and_reciprocal_heat_entropy() {
        let curve = phase_curve();
        let offset = 1e12;
        let shifted = EquilibriumEnthalpyPhaseCurve::try_new(
            ContentHash([0x7a; 32]),
            curve
                .knots()
                .iter()
                .map(|knot| EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: knot.specific_enthalpy_j_kg + offset,
                    ..*knot
                })
                .collect(),
        )
        .unwrap();
        assert_eq!(
            curve
                .thermodynamic_mean_temperature(1000.0, 77_000.0, 3, || false)
                .unwrap(),
            shifted
                .thermodynamic_mean_temperature(offset + 1000.0, offset + 77_000.0, 3, || false,)
                .unwrap(),
        );

        // One signed 10 kJ transfer between two 1 kg bodies has no separate
        // latent account. The exact finite entropy uses the two mean efforts.
        let hot = curve
            .thermodynamic_mean_temperature(95_000.0, 85_000.0, 1, || false)
            .unwrap();
        let cold = curve
            .thermodynamic_mean_temperature(0.0, 10_000.0, 1, || false)
            .unwrap();
        let entropy_change =
            hot.specific_entropy_change_j_kg_k + cold.specific_entropy_change_j_kg_k;
        let exchange_entropy = 10_000.0 * (cold.temperature_k.recip() - hot.temperature_k.recip());
        assert!(hot.temperature_k > cold.temperature_k);
        assert!(entropy_change > 0.0);
        assert!((entropy_change - exchange_entropy).abs() < 1e-14);
    }

    #[test]
    fn mean_queries_enforce_original_interval_work_and_cancellation() {
        let curve = EquilibriumEnthalpyPhaseCurve::try_single_phase(
            ContentHash([0x7b; 32]),
            SolidLiquidPhase::Solid,
            (0..=129)
                .map(|index| EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: 1000.0 * f64::from(index),
                    temperature_k: 300.0 + f64::from(index),
                    liquid_mass_fraction: 0.0,
                    bulk_density_kg_m3: 1000.0,
                })
                .collect(),
        )
        .unwrap();
        let full = curve
            .thermodynamic_mean_temperature(0.0, 129_000.0, 129, || false)
            .unwrap();
        assert_eq!(full.intervals, 129);
        assert_eq!(
            full.specific_entropy_change_j_kg_k,
            curve
                .specific_entropy_change_j_kg_k(0.0, 129_000.0)
                .unwrap()
        );
        let local = curve
            .thermodynamic_mean_temperature(64_250.0, 64_750.0, 1, || false)
            .unwrap();
        assert_eq!(local.intervals, 1);
        let mut polls = 0;
        let refused = curve.thermodynamic_mean_temperature(0.0, 129_000.0, 128, || {
            polls += 1;
            false
        });
        assert_eq!(
            polls, 1,
            "the budget must refuse before traversing any interval"
        );
        assert_eq!(
            refused,
            Err(PhaseStateError::InvalidCurve {
                what: "thermal-entropy original-interval budget exhausted",
            })
        );
        assert_eq!(
            curve.thermodynamic_mean_temperature(0.0, 0.0, 0, || true),
            Err(PhaseStateError::EntropyEvaluationCancelled)
        );
        polls = 0;
        let cancelled = curve.thermodynamic_mean_temperature(0.0, 129_000.0, 129, || {
            polls += 1;
            polls == 3
        });
        assert_eq!(polls, 3);
        assert_eq!(cancelled, Err(PhaseStateError::EntropyEvaluationCancelled));
        // A cancellation arriving after the one visited interval still refuses.
        polls = 0;
        let cancelled = curve.thermodynamic_mean_temperature(64_250.0, 64_750.0, 1, || {
            polls += 1;
            polls == 3
        });
        assert_eq!(polls, 3);
        assert_eq!(cancelled, Err(PhaseStateError::EntropyEvaluationCancelled));
    }

    #[test]
    fn unrepresentable_sensible_derivative_does_not_become_a_latent_plateau() {
        let knot = |h, t| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: t,
            liquid_mass_fraction: 0.0,
            bulk_density_kg_m3: 1000.0,
        };
        let curve = EquilibriumEnthalpyPhaseCurve::try_single_phase(
            ContentHash([0x7c; 32]),
            SolidLiquidPhase::Solid,
            vec![knot(0.0, 1.0), knot(1e308, 1.0_f64.next_up())],
        )
        .unwrap();
        assert_eq!(curve.specific_entropy_change_j_kg_k(0.0, 1.0).unwrap(), 1.0);
        for to in [0.0, 1.0] {
            assert_eq!(
                curve.thermodynamic_mean_temperature(0.0, to, 1, || false),
                Err(PhaseStateError::UnrepresentableTemperatureDerivative),
            );
        }
    }
}
