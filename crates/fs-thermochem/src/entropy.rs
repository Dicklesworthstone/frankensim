//! Cancellation-safe differences of the retained NASA-9 entropy laws.
//!
//! Reference entropy and ideal mixing cancel at frozen composition. Evaluate
//! those differences algebraically before rounding, rather than subtracting
//! two large absolute entropies. No path, reaction, or transport is inferred.

use crate::mixture::{
    FrozenIdealGasMixtureErrorV1, FrozenIdealGasMixtureModelV1, FrozenMixturePropertyV1,
};
use crate::{
    MassSpecificThermalQuantityV1, MolarEntropyV1, Nasa9PropertyV1, Nasa9StandardStateModelV1,
    ThermochemErrorV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K,
};
use fs_qty::{Pressure, Temperature};

/// ln(a/b), retaining small representable differences and avoiding ratio overflow.
/// Ordering the arguments also makes reversal exactly antisymmetric.
pub(crate) fn log_ratio(a: f64, b: f64) -> f64 {
    if a == b {
        return 0.0;
    }
    let (lo, hi, sign) = if a > b { (b, a, 1.0) } else { (a, b, -1.0) };
    let delta = (hi - lo) / lo;
    let rounded = 1.0 + delta;
    let value = if rounded == 1.0 {
        delta
    } else if rounded.is_finite() {
        // Correct the rounding of 1+delta (Kahan's log1p identity).
        fs_math::det::ln(rounded) * (delta / (rounded - 1.0))
    } else {
        fs_math::det::ln(hi) - fs_math::det::ln(lo)
    };
    sign * value
}

pub(crate) fn accumulate(sum: &mut f64, correction: &mut f64, value: f64) {
    let increment = value - *correction;
    let next = *sum + increment;
    *correction = (next - *sum) - increment;
    *sum = next;
}

impl Nasa9StandardStateModelV1 {
    /// Molar standard-state entropy change between two source-supported endpoints.
    ///
    /// Raw models retain original region selection and reference jumps.
    /// Explicitly reconciled models integrate Cp/T across each crossed region,
    /// retaining small changes through a continuous caloric seam. A common
    /// entropy reference cancels before floating-point addition.
    ///
    /// Raw endpoint support does not prove continuous source-path coverage.
    /// Fixed operation order, no allocation on success. Units: J/(mol K).
    pub fn molar_entropy_change(
        &self,
        initial: Temperature,
        final_temperature: Temperature,
    ) -> Result<MolarEntropyV1, ThermochemErrorV1> {
        let t0 = initial.value();
        let t1 = final_temperature.value();
        for t in [t0, t1] {
            if !t.is_finite() || t <= 0.0 {
                return Err(ThermochemErrorV1::InvalidEvaluationTemperature { bits: t.to_bits() });
            }
            if self.region_for_temperature(t).is_none() {
                return Err(ThermochemErrorV1::TemperatureOutOfRange { bits: t.to_bits() });
            }
        }
        if t0 == t1 {
            return Ok(MolarEntropyV1::new(0.0));
        }
        let (lo, hi, sign) = if t1 > t0 {
            (t0, t1, 1.0)
        } else {
            (t1, t0, -1.0)
        };
        let (first, lower) = self.region_for_temperature(lo).expect("validated lower endpoint");
        let (last, upper) = self.region_for_temperature(hi).expect("validated upper endpoint");
        let value = sign * if self.caloric_reference().is_some() && first != last {
            self.integrate_caloric_regions(first, last, lo, hi,
                |coefficients, start, end| entropy_difference(coefficients, coefficients, start, end))
        } else {
            entropy_difference(lower.coefficients, upper.coefficients, lo, hi)
        };
        if !value.is_finite() {
            return Err(ThermochemErrorV1::NonFiniteEvaluation {
                property: Nasa9PropertyV1::Entropy,
                bits: value.to_bits(),
            });
        }
        Ok(MolarEntropyV1::new(value))
    }
}

// Original polynomial-difference arithmetic, shared by raw endpoint differences
// and the source Cp/T integral on each reconciled temperature segment.
fn entropy_difference(lower: [f64; 9], upper: [f64; 9], lo: f64, hi: f64) -> f64 {
    let weights = |a: [f64; 9]| {
        [-0.5 * a[0], -a[1], a[2], a[3], a[4] / 2.0, a[5] / 3.0, a[6] / 4.0, a[8]]
    };
    let lower = weights(lower);
    let upper = weights(upper);
    let inverse_lo = 1.0 / lo;
    let inverse_hi = 1.0 / hi;
    let lo2 = lo * lo;
    let hi2 = hi * hi;
    let delta = hi - lo;
    let inverse_delta = -(delta / hi) / lo;
    let basis = [
        inverse_lo * inverse_lo, inverse_lo, fs_math::det::ln(lo),
        lo, lo2, lo2 * lo, lo2 * lo2, 1.0,
    ];
    let differences = [
        inverse_delta * (inverse_lo + inverse_hi),
        inverse_delta,
        log_ratio(hi, lo),
        delta,
        delta * (lo + hi),
        delta * (lo2 + lo * hi + hi2),
        (delta * (lo + hi)) * (lo2 + hi2),
        0.0,
    ];
    let mut sum = 0.0;
    let mut correction = 0.0;
    for i in 0..8 {
        if upper[i] != lower[i] {
            accumulate(&mut sum, &mut correction, (upper[i] - lower[i]) * basis[i]);
        }
        if upper[i] != 0.0 && differences[i] != 0.0 {
            accumulate(&mut sum, &mut correction, upper[i] * differences[i]);
        }
    }
    UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * sum
}

impl FrozenIdealGasMixtureModelV1 {
    /// Frozen-composition molar entropy change, including the pressure correction.
    ///
    /// Ideal mixing and each common species entropy reference cancel exactly.
    /// Composition, species and conventions belong to this immutable model at
    /// BOTH endpoints. This is not an entropy change for reacting or mixing gas.
    pub fn molar_entropy_change(
        &self,
        initial_temperature: Temperature,
        initial_pressure: Pressure,
        final_temperature: Temperature,
        final_pressure: Pressure,
    ) -> Result<MolarEntropyV1, FrozenIdealGasMixtureErrorV1> {
        self.entropy_difference(
            initial_temperature,
            initial_pressure,
            final_temperature,
            final_pressure,
        )
        .map(|(entropy, _)| MolarEntropyV1::new(entropy))
    }

    /// Same frozen-mixture entropy change in J/(kg K), using its own molar mass.
    pub fn specific_entropy_change(
        &self,
        initial_temperature: Temperature,
        initial_pressure: Pressure,
        final_temperature: Temperature,
        final_pressure: Pressure,
    ) -> Result<MassSpecificThermalQuantityV1, FrozenIdealGasMixtureErrorV1> {
        let (entropy, molar_mass) = self.entropy_difference(
            initial_temperature,
            initial_pressure,
            final_temperature,
            final_pressure,
        )?;
        let value = entropy / molar_mass;
        if !value.is_finite() {
            return Err(FrozenIdealGasMixtureErrorV1::NonFiniteEvaluation {
                property: FrozenMixturePropertyV1::MolarEntropy,
                bits: value.to_bits(),
            });
        }
        Ok(MassSpecificThermalQuantityV1::new(value))
    }

    fn entropy_difference(
        &self,
        initial_temperature: Temperature,
        initial_pressure: Pressure,
        final_temperature: Temperature,
        final_pressure: Pressure,
    ) -> Result<(f64, f64), FrozenIdealGasMixtureErrorV1> {
        for p in [initial_pressure.value(), final_pressure.value()] {
            if !p.is_finite() || p <= 0.0 {
                return Err(FrozenIdealGasMixtureErrorV1::InvalidPressure { bits: p.to_bits() });
            }
        }
        let mut entropy = 0.0;
        let mut correction = 0.0;
        let mut molar_mass = 0.0;
        for (component, (model, &fraction)) in self
            .components()
            .iter()
            .zip(self.mole_composition().fractions())
            .enumerate()
        {
            let change = model
                .molar_entropy_change(initial_temperature, final_temperature)
                .map_err(|source| FrozenIdealGasMixtureErrorV1::ComponentEvaluation {
                    component,
                    species: model.species().clone(),
                    source: Box::new(source),
                })?;
            accumulate(&mut entropy, &mut correction, fraction * change.value());
            molar_mass += fraction * model.molar_mass().value();
        }
        accumulate(
            &mut entropy,
            &mut correction,
            -UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K
                * log_ratio(final_pressure.value(), initial_pressure.value()),
        );
        if !entropy.is_finite() {
            return Err(FrozenIdealGasMixtureErrorV1::NonFiniteEvaluation {
                property: FrozenMixturePropertyV1::MolarEntropy,
                bits: entropy.to_bits(),
            });
        }
        if !molar_mass.is_finite() || molar_mass <= 0.0 {
            return Err(FrozenIdealGasMixtureErrorV1::InvalidMixtureMolarMass {
                bits: molar_mass.to_bits(),
            });
        }
        Ok((entropy, molar_mass))
    }
}
