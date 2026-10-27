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
    /// Returns the original state/domain errors, or `InvalidCurve` when the
    /// integral cannot be represented. No entropy-production claim is made by
    /// this state function; the transport owner must audit its boundary entropy.
    pub fn specific_entropy_change_j_kg_k(
        &self,
        initial_enthalpy_j_kg: f64,
        final_enthalpy_j_kg: f64,
    ) -> Result<f64, PhaseStateError> {
        self.state_at_specific_enthalpy(initial_enthalpy_j_kg)?;
        self.state_at_specific_enthalpy(final_enthalpy_j_kg)?;
        if initial_enthalpy_j_kg == final_enthalpy_j_kg {
            return Ok(0.0);
        }
        let lo = initial_enthalpy_j_kg.min(final_enthalpy_j_kg);
        let hi = initial_enthalpy_j_kg.max(final_enthalpy_j_kg);
        let invalid = || PhaseStateError::InvalidCurve {
            what: "relative thermal entropy is not finite and representable",
        };
        let mut total = 0.0;
        let mut correction = 0.0;
        for pair in self.knots().windows(2) {
            let a = lo.max(pair[0].specific_enthalpy_j_kg);
            let b = hi.min(pair[1].specific_enthalpy_j_kg);
            if b <= a {
                continue;
            }
            let ta = self.state_at_specific_enthalpy(a)?.temperature_k();
            let tb = self.state_at_specific_enthalpy(b)?.temperature_k();
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
    use crate::phase::{EnthalpyPhaseKnot, SolidLiquidPhase};
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
            vec![knot(0.0, 300.0, 0.0), knot(30_000.0, 600.0, 0.0),
                knot(55_000.0, 600.0, 1.0), knot(95_000.0, 800.0, 1.0)],
        ).unwrap()
    }

    #[test]
    fn caloric_and_latent_entropy_share_one_chart_and_are_path_additive() {
        let curve = phase_curve();
        let entropy = |a, b| curve.specific_entropy_change_j_kg_k(a, b).unwrap();
        let expected = 100.0 * 2.0_f64.ln() + 25_000.0 / 600.0
            + 200.0 * (800.0_f64 / 600.0).ln();
        assert!((entropy(0.0, 95_000.0) - expected).abs() < 1e-12);
        assert_eq!(entropy(40_000.0, 50_000.0), 10_000.0 / 600.0);
        assert_eq!(entropy(95_000.0, 0.0), -entropy(0.0, 95_000.0));
        let split = entropy(0.0, 20_000.0) + entropy(20_000.0, 48_000.0)
            + entropy(48_000.0, 95_000.0);
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
        let knot = |h, t| EnthalpyPhaseKnot { specific_enthalpy_j_kg: h,
            temperature_k: t, liquid_mass_fraction: 0.0, bulk_density_kg_m3: 1000.0 };
        let curve = EquilibriumEnthalpyPhaseCurve::try_single_phase(
            ContentHash([0x75; 32]), SolidLiquidPhase::Solid,
            vec![knot(0.0, ta), knot(1.0, tb)],
        ).unwrap();
        let entropy = curve.specific_entropy_change_j_kg_k(0.0, 1.0).unwrap();
        assert!(entropy >= 1.0 / tb && entropy <= 1.0 / ta);
        assert!((entropy * 0.5 * (ta + tb) - 1.0).abs() < 5e-16);
    }
}
