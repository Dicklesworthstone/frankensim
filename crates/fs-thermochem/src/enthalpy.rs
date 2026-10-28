//! Reference-invariant sensible enthalpy changes of the retained NASA-9 laws.
//!
//! Nozzle kinetic energy depends on h(T0)-h(T), not on a formation-energy
//! offset. Subtract the polynomial laws before evaluation so that a large
//! common reference cannot erase the available energy. This does not remove
//! chemical energy when composition or source conventions change.

use crate::entropy::{accumulate, log_ratio};
use crate::mixture::{
    FrozenIdealGasMixtureErrorV1, FrozenIdealGasMixtureModelV1, FrozenMixturePropertyV1,
};
use crate::{
    MassSpecificEnergyQuantityV1, MolarEnthalpyV1, Nasa9PropertyV1, Nasa9StandardStateModelV1,
    ThermochemErrorV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K,
};
use fs_qty::Temperature;

impl Nasa9StandardStateModelV1 {
    /// Molar enthalpy change between two source-supported temperatures, J/mol.
    ///
    /// Raw models retain distinct coefficients and reference jumps. Explicitly
    /// reconciled models integrate the retained Cp chart over each crossed
    /// region, preserving small sensible changes through the continuous seam.
    /// A common formation offset cancels before rounding in either mode.
    ///
    /// Raw endpoint support does not prove support of the intervening path.
    /// Fixed operation order, exact reversal, no allocation on success.
    pub fn molar_enthalpy_change(
        &self,
        initial: Temperature,
        final_temperature: Temperature,
    ) -> Result<MolarEnthalpyV1, ThermochemErrorV1> {
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
            return Ok(MolarEnthalpyV1::new(0.0));
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
                |coefficients, start, end| enthalpy_difference(coefficients, coefficients, start, end))
        } else {
            enthalpy_difference(lower.coefficients, upper.coefficients, lo, hi)
        };
        if !value.is_finite() {
            return Err(ThermochemErrorV1::NonFiniteEvaluation {
                property: Nasa9PropertyV1::Enthalpy,
                bits: value.to_bits(),
            });
        }
        Ok(MolarEnthalpyV1::new(value))
    }
}

// Shared original polynomial-difference arithmetic. Equal coefficient arrays
// integrate Cp on one chart; distinct arrays retain raw source-reference jumps.
fn enthalpy_difference(lower: [f64; 9], upper: [f64; 9], lo: f64, hi: f64) -> f64 {
    let weights = |a: [f64; 9]| {
        [-a[0], a[1], a[2], a[3] / 2.0, a[4] / 3.0, a[5] / 4.0, a[6] / 5.0, a[7]]
    };
    let lower = weights(lower);
    let upper = weights(upper);
    let delta = hi - lo;
    let lo2 = lo * lo;
    let hi2 = hi * hi;
    let lo3 = lo2 * lo;
    let hi3 = hi2 * hi;
    let lo4 = lo2 * lo2;
    let hi4 = hi2 * hi2;
    // h/R = -a0/T + a1 ln(T) + a2 T + a3 T^2/2 + ... + a6 T^5/5 + a7.
    let basis = [
        1.0 / lo,
        if upper[1] != lower[1] { fs_math::det::ln(lo) } else { 0.0 },
        lo, lo2, lo3, lo4, lo4 * lo, 1.0,
    ];
    let differences = [
        -(delta / hi) / lo,
        if upper[1] != 0.0 { log_ratio(hi, lo) } else { 0.0 },
        delta,
        delta * (lo + hi),
        delta * (lo2 + lo * hi + hi2),
        (delta * (lo + hi)) * (lo2 + hi2),
        delta * (lo4 + lo3 * hi + lo2 * hi2 + lo * hi3 + hi4),
        0.0,
    ];
    let mut sum = 0.0;
    let mut correction = 0.0;
    for i in 0..8 {
        // Skip inactive powers before multiplication, avoiding 0*infinity.
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
    /// Frozen-composition molar enthalpy change, J/mol.
    ///
    /// The ideal-gas enthalpy depends only on temperature. Each common species
    /// formation reference cancels before the canonical compensated mixture
    /// reduction. This is not the chemical energy of reaction or mixing.
    pub fn molar_enthalpy_change(
        &self,
        initial: Temperature,
        final_temperature: Temperature,
    ) -> Result<MolarEnthalpyV1, FrozenIdealGasMixtureErrorV1> {
        self.enthalpy_difference(initial, final_temperature)
            .map(|(enthalpy, _)| MolarEnthalpyV1::new(enthalpy))
    }

    /// Same frozen-mixture enthalpy change in J/kg, using its own molar mass.
    pub fn specific_enthalpy_change(
        &self,
        initial: Temperature,
        final_temperature: Temperature,
    ) -> Result<MassSpecificEnergyQuantityV1, FrozenIdealGasMixtureErrorV1> {
        let (enthalpy, molar_mass) = self.enthalpy_difference(initial, final_temperature)?;
        let value = enthalpy / molar_mass;
        if !value.is_finite() {
            return Err(FrozenIdealGasMixtureErrorV1::NonFiniteEvaluation {
                property: FrozenMixturePropertyV1::MolarEnthalpy,
                bits: value.to_bits(),
            });
        }
        Ok(MassSpecificEnergyQuantityV1::new(value))
    }

    fn enthalpy_difference(
        &self,
        initial: Temperature,
        final_temperature: Temperature,
    ) -> Result<(f64, f64), FrozenIdealGasMixtureErrorV1> {
        let mut enthalpy = 0.0;
        let mut correction = 0.0;
        let mut molar_mass = 0.0;
        for (component, (model, &fraction)) in self
            .components()
            .iter()
            .zip(self.mole_composition().fractions())
            .enumerate()
        {
            let change = model
                .molar_enthalpy_change(initial, final_temperature)
                .map_err(|source| FrozenIdealGasMixtureErrorV1::ComponentEvaluation {
                    component,
                    species: model.species().clone(),
                    source: Box::new(source),
                })?;
            accumulate(&mut enthalpy, &mut correction, fraction * change.value());
            molar_mass += fraction * model.molar_mass().value();
        }
        if !enthalpy.is_finite() {
            return Err(FrozenIdealGasMixtureErrorV1::NonFiniteEvaluation {
                property: FrozenMixturePropertyV1::MolarEnthalpy,
                bits: enthalpy.to_bits(),
            });
        }
        if !molar_mass.is_finite() || molar_mass <= 0.0 {
            return Err(FrozenIdealGasMixtureErrorV1::InvalidMixtureMolarMass {
                bits: molar_mass.to_bits(),
            });
        }
        Ok((enthalpy, molar_mass))
    }
}
