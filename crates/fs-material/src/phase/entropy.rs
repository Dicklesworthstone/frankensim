//! Relative thermal entropy on the original equilibrium enthalpy chart.

use super::{EquilibriumEnthalpyPhaseCurve, PhaseStateError};

impl EquilibriumEnthalpyPhaseCurve {
    /// Integrate `ds = dh / T(h)` between two admitted specific enthalpies.
    ///
    /// Each nonisothermal segment uses its exact logarithmic integral; a latent
    /// plateau contributes `delta_h / T_transition` without an apparent heat
    /// capacity or separate latent-energy credit. This integrates the admitted
    /// piecewise-linear chart, including one generated from heat-capacity data,
    /// not a higher-resolution source law that the chart only approximates.
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
        for pair in knots[first..=last].windows(2) {
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
            let dh = b - a;
            let relative_temperature = (tb - ta) / ta;
            let ratio = if relative_temperature == 0.0 {
                1.0
            } else {
                relative_temperature.ln_1p() / relative_temperature
            };
            let entropy = (dh / ta) * ratio;
            if !relative_temperature.is_finite() || !entropy.is_finite() || entropy <= 0.0 {
                return Err(invalid());
            }
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
        }
    }
}
