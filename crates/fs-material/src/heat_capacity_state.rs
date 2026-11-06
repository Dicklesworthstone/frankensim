//! Temperature-coordinate operations on the existing sourced single-phase
//! heat-capacity chart. No second curve, source selection or integration model.

use crate::phase::{
    EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState, PhaseStateError,
    ResolvedHeatCapacityCurve,
};

impl ResolvedHeatCapacityCurve {
    /// Invert the admitted single-phase chart at an absolute temperature [K].
    /// The result retains the original material and chart identities. No
    /// extrapolation, phase transition, or independent density law is added.
    ///
    /// # Errors
    /// Refuses nonfinite/out-of-domain temperature or unrepresentable enthalpy.
    pub fn state_at_temperature(
        &self,
        temperature_k: f64,
    ) -> Result<EquilibriumPhaseState, PhaseStateError> {
        let h = enthalpy_at_temperature(self.curve(), temperature_k)?;
        self.curve().state_at_specific_enthalpy(h)
    }

    /// Specific enthalpy change `h(to) - h(from)` [J/kg] on the admitted
    /// original caloric law. Integrating its crossed intervals avoids
    /// subtracting a possibly large arbitrary reference enthalpy twice.
    /// Heating and cooling use the same law with opposite signs.
    ///
    /// Sourced sensible intervals retain and integrate their affine Cp(T).
    /// Explicitly authored enthalpy intervals retain their chord law. Source
    /// uncertainty remains separate from floating-point arithmetic.
    ///
    /// # Errors
    /// Both temperatures must lie inside the chart; a nonzero unrepresentable
    /// change refuses rather than being silently rounded to no stored energy.
    pub fn specific_enthalpy_change(&self, from_k: f64, to_k: f64) -> Result<f64, PhaseStateError> {
        enthalpy_change(self.curve(), from_k, to_k)
    }

    /// Relative thermal entropy [J/(kg K)] between two absolute temperatures.
    /// Integrates the original single-phase Cp/T law, including
    /// crossed knots, without converting either temperature to an absolute h.
    /// A large arbitrary enthalpy datum therefore cannot erase a small real
    /// entropy increment. Retained source Cp and authored chord laws share the
    /// same source support and original interval traversal.
    ///
    /// The same fixed source pressure/composition convention as the chart
    /// applies. Mixing entropy, pressure work and entropy production belong to
    /// the coupled owner. Reversal changes the sign; no extrapolation is allowed.
    ///
    /// # Errors
    /// Refuses out-of-domain temperatures or an unrepresentable nonzero integral.
    pub fn specific_entropy_change(&self, from_k: f64, to_k: f64) -> Result<f64, PhaseStateError> {
        entropy_change(self.curve(), from_k, to_k)
    }

    /// Positive secant heat capacity [J/(kg K)] for a conservative thermal
    /// linearization. Its product with `(to_k - from_k)` is the admitted chart's
    /// enthalpy change. At identical temperatures use the local source Cp or
    /// the authored right chord slope at
    /// interior knots and the left slope at the upper endpoint.
    ///
    /// # Errors
    /// The same domain and representability refusals as the chart operations.
    pub fn secant_heat_capacity(&self, from_k: f64, to_k: f64) -> Result<f64, PhaseStateError> {
        secant(self.curve(), from_k, to_k)
    }
}

fn refusal(reason: impl Into<String>) -> PhaseStateError {
    PhaseStateError::SourceQuery {
        property: "single-phase temperature".into(),
        reason: reason.into(),
    }
}

fn segment(
    curve: &EquilibriumEnthalpyPhaseCurve,
    temperature: f64,
) -> Result<usize, PhaseStateError> {
    let knots = curve.knots();
    let lo = knots[0].temperature_k;
    let hi = knots[knots.len() - 1].temperature_k;
    if !temperature.is_finite() || temperature < lo || temperature > hi {
        return Err(refusal(format!(
            "temperature {temperature} K is outside [{lo}, {hi}] K"
        )));
    }
    Ok(knots
        .partition_point(|k| k.temperature_k <= temperature)
        .min(knots.len() - 1)
        - 1)
}

fn slope(left: EnthalpyPhaseKnot, right: EnthalpyPhaseKnot) -> Result<f64, PhaseStateError> {
    let value = (right.specific_enthalpy_j_kg - left.specific_enthalpy_j_kg)
        / (right.temperature_k - left.temperature_k);
    if !value.is_finite() || value <= 0.0 {
        return Err(refusal(
            "single-phase heat capacity is not representable and positive",
        ));
    }
    Ok(value)
}

fn capacity(
    curve: &EquilibriumEnthalpyPhaseCurve,
    index: usize,
    temperature: f64,
) -> Result<f64, PhaseStateError> {
    let (left, right) = (curve.knots()[index], curve.knots()[index + 1]);
    let Some([a, b]) = curve.sensible_heat_capacities_on_interval(index) else {
        return slope(left, right);
    };
    let fraction = (temperature - left.temperature_k) / (right.temperature_k - left.temperature_k);
    let value = (1.0 - fraction) * a + fraction * b;
    if !value.is_finite() || value <= 0.0 {
        return Err(refusal("source heat capacity is not finite and positive"));
    }
    Ok(value)
}

pub(crate) fn enthalpy_at_temperature(
    curve: &EquilibriumEnthalpyPhaseCurve,
    temperature: f64,
) -> Result<f64, PhaseStateError> {
    let index = segment(curve, temperature)?;
    curve.specific_enthalpy_on_sensible_interval(index, temperature)
}

fn enthalpy_change(
    curve: &EquilibriumEnthalpyPhaseCurve,
    from: f64,
    to: f64,
) -> Result<f64, PhaseStateError> {
    let a = segment(curve, from)?;
    let b = segment(curve, to)?;
    if from == to {
        return Ok(0.0);
    }
    let (lo, hi) = (from.min(to), from.max(to));
    let mut change = 0.0;
    let mut correction = 0.0;
    for index in a.min(b)..=a.max(b) {
        let lower = curve.knots()[index].temperature_k.max(lo);
        let upper = curve.knots()[index + 1].temperature_k.min(hi);
        let width = upper - lower;
        if width > 0.0 {
            let mean = 0.5 * capacity(curve, index, lower)? + 0.5 * capacity(curve, index, upper)?;
            let increment = mean * width - correction;
            let next = change + increment;
            correction = (next - change) - increment;
            change = next;
        }
    }
    if !change.is_finite() || change <= 0.0 {
        return Err(refusal(
            "nonzero enthalpy change is not finite and representable",
        ));
    }
    Ok(if to > from { change } else { -change })
}

fn entropy_change(
    curve: &EquilibriumEnthalpyPhaseCurve,
    from: f64,
    to: f64,
) -> Result<f64, PhaseStateError> {
    let a = segment(curve, from)?;
    let b = segment(curve, to)?;
    if from == to {
        return Ok(0.0);
    }
    let (lo, hi) = (from.min(to), from.max(to));
    let mut total = 0.0;
    let mut correction = 0.0;
    for (offset, pair) in curve.knots()[a.min(b)..=a.max(b) + 1].windows(2).enumerate() {
        let lower = lo.max(pair[0].temperature_k);
        let upper = hi.min(pair[1].temperature_k);
        if upper <= lower {
            continue;
        }
        let index = a.min(b) + offset;
        let value = if curve.sensible_heat_capacities_on_interval(index).is_some() {
            crate::phase::affine_cp_entropy(
                lower,
                upper,
                capacity(curve, index, lower)?,
                capacity(curve, index, upper)?,
            )?
        } else {
            let delta = (upper - lower) / lower;
            let rounded = 1.0 + delta;
            let logarithm = if rounded == 1.0 {
                delta
            } else if rounded.is_finite() {
                // Correct rounding of 1+delta before deterministic ln; unlike
                // subtracting two logarithms this retains arbitrarily close T.
                fs_math::det::ln(rounded) * (delta / (rounded - 1.0))
            } else {
                fs_math::det::ln(upper) - fs_math::det::ln(lower)
            };
            slope(pair[0], pair[1])? * logarithm
        };
        if !value.is_finite() || value <= 0.0 {
            return Err(refusal(
                "nonzero thermal entropy is not finite and representable",
            ));
        }
        let increment = value - correction;
        let next = total + increment;
        correction = (next - total) - increment;
        total = next;
    }
    if !total.is_finite() || total <= 0.0 {
        return Err(refusal(
            "nonzero thermal entropy is not finite and representable",
        ));
    }
    Ok(if to < from { -total } else { total })
}

fn secant(
    curve: &EquilibriumEnthalpyPhaseCurve,
    from: f64,
    to: f64,
) -> Result<f64, PhaseStateError> {
    if from == to {
        let index = segment(curve, from)?;
        return capacity(curve, index, from);
    }
    let value = enthalpy_change(curve, from, to)? / (to - from);
    if !value.is_finite() || value <= 0.0 {
        return Err(refusal("enthalpy secant is not finite and positive"));
    }
    Ok(value)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::phase::{HeatCapacityKnot, SolidLiquidPhase};
    use fs_blake3::hash_bytes;

    fn curve(reference: f64) -> EquilibriumEnthalpyPhaseCurve {
        EquilibriumEnthalpyPhaseCurve::try_single_phase(
            hash_bytes(b"synthetic thermal chart"),
            SolidLiquidPhase::Solid,
            [(300.0, 0.0), (310.0, 20.0), (320.0, 60.0)]
                .map(|(t, h)| EnthalpyPhaseKnot {
                    temperature_k: t,
                    specific_enthalpy_j_kg: reference + h,
                    liquid_mass_fraction: 0.0,
                    bulk_density_kg_m3: 10.0,
                })
                .to_vec(),
        )
        .unwrap()
    }

    #[test]
    fn original_cp_integrals_ignore_chord_allowance_and_large_energy_datums() {
        // Cp(T)=2T, so delta_h=T1^2-T0^2 and delta_s=2(T1-T0).
        // A deliberately coarse boundary chart must still retain this law.
        for reference in [0.0, 1e12] {
            let c = EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity(
                hash_bytes(b"affine Cp source, no interpolation calibration"),
                SolidLiquidPhase::Solid,
                reference,
                &[
                    HeatCapacityKnot { temperature_k: 300.0, specific_heat_capacity_j_kg_k: 600.0, bulk_density_kg_m3: 1000.0 },
                    HeatCapacityKnot { temperature_k: 600.0, specific_heat_capacity_j_kg_k: 1200.0, bulk_density_kg_m3: 900.0 },
                ],
                1e6,
            ).unwrap();
            assert_eq!(c.knots().len(), 2);
            for (a, b) in [(300.0_f64, 500.0_f64), (500.0, 300.0), (350.25, 350.25_f64.next_up())] {
                let expected_h = (a + b) * (b - a);
                let expected_s = 2.0 * (b - a);
                assert!((enthalpy_change(&c, a, b).unwrap() - expected_h).abs() <= 16.0 * f64::EPSILON * expected_h.abs());
                assert!((entropy_change(&c, a, b).unwrap() - expected_s).abs() <= 16.0 * f64::EPSILON * expected_s.abs());
                assert!((secant(&c, a, b).unwrap() - (a + b)).abs() <= 16.0 * f64::EPSILON * (a + b));
            }
            assert_eq!(secant(&c, 350.25, 350.25).unwrap(), 700.5);
            assert!((enthalpy_at_temperature(&c, 350.25).unwrap() - reference - 32675.0625).abs() < 1e-8);
        }
    }

    #[test]
    fn g1_g3_chart_storage_is_signed_additive_and_reference_independent() {
        for reference in [0.0, -1000.0, 1.0e12] {
            let c = curve(reference);
            assert_eq!(enthalpy_change(&c, 305.0, 315.0).unwrap(), 30.0);
            assert_eq!(enthalpy_change(&c, 315.0, 305.0).unwrap(), -30.0);
            assert_eq!(secant(&c, 305.0, 315.0).unwrap(), 3.0);
            assert_eq!(secant(&c, 315.0, 305.0).unwrap(), 3.0);
            let parts = enthalpy_change(&c, 305.0, 310.0).unwrap()
                + enthalpy_change(&c, 310.0, 315.0).unwrap();
            assert_eq!(parts, 30.0);
            assert_eq!(
                enthalpy_at_temperature(&c, 315.0).unwrap(),
                reference + 40.0
            );
        }
    }

    #[test]
    fn g0_chart_knots_have_explicit_one_sided_capacity() {
        let c = curve(0.0);
        assert_eq!(secant(&c, 300.0, 300.0).unwrap(), 2.0);
        assert_eq!(secant(&c, 310.0, 310.0).unwrap(), 4.0);
        assert_eq!(secant(&c, 320.0, 320.0).unwrap(), 4.0);
        assert_eq!(enthalpy_change(&c, 310.0, 310.0).unwrap(), 0.0);
        for t in [299.0, 321.0, f64::NAN, f64::INFINITY] {
            assert!(enthalpy_at_temperature(&c, t).is_err());
            assert!(secant(&c, t, t).is_err());
            assert!(enthalpy_change(&c, 305.0, t).is_err());
        }
    }

    #[test]
    fn g3_chart_roundtrip_retains_state_identity() {
        let c = curve(-100.0);
        for t in [300.0, 305.0, 310.0, 319.0, 320.0] {
            let h = enthalpy_at_temperature(&c, t).unwrap();
            let state = c.state_at_specific_enthalpy(h).unwrap();
            assert_eq!(state.temperature_k(), t);
            assert_eq!(state.phase_curve_identity(), c.identity());
            assert_eq!(state.material_card_identity(), c.material_card_identity());
        }
    }
    #[test]
    fn g1_g3_temperature_entropy_preserves_small_changes_across_large_enthalpy_datums() {
        let from = 301.0;
        let to = from + 1e-8;
        let c = curve(1e12);
        assert_eq!(
            enthalpy_at_temperature(&c, from).unwrap(),
            enthalpy_at_temperature(&c, to).unwrap()
        );
        let entropy = entropy_change(&c, from, to).unwrap();
        let expected = 2.0 * ((to - from) / from).ln_1p();
        assert!(entropy > 0.0);
        assert!((entropy - expected).abs() < 1e-22);
        assert_eq!(entropy_change(&c, to, from).unwrap(), -entropy);
        assert_eq!(entropy_change(&curve(0.0), from, to).unwrap(), entropy);
        let whole = entropy_change(&c, 300.0, 320.0).unwrap();
        let split =
            entropy_change(&c, 300.0, 310.0).unwrap() + entropy_change(&c, 310.0, 320.0).unwrap();
        let exact = 2.0 * (310.0_f64 / 300.0).ln() + 4.0 * (320.0_f64 / 310.0).ln();
        assert!((whole - exact).abs() < 1e-14);
        assert!((whole - split).abs() < 1e-14);
        assert_eq!(entropy_change(&c, from, from).unwrap(), 0.0);
        for t in [299.0, 321.0, f64::NAN, f64::INFINITY] {
            assert!(entropy_change(&c, from, t).is_err());
        }
    }
}
