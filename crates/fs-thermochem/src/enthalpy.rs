//! Reference-invariant sensible energy changes of the retained NASA-9 laws.
//!
//! Cancel formation references before evaluation. Internal-energy changes
//! integrate the same polynomial with cp/R replaced by cv/R = cp/R - 1,
//! rather than subtracting two nearly equal enthalpy and pressure-work terms.
//! This does not remove chemical energy when composition or sources change.

use crate::entropy::{accumulate, log_ratio};
use crate::mixture::{
    FrozenIdealGasMixtureErrorV1, FrozenIdealGasMixtureModelV1, FrozenMixturePropertyV1,
};
use crate::{
    MassSpecificEnergyQuantityV1, MolarEnthalpyV1, MolarInternalEnergyV1, Nasa9PropertyV1,
    Nasa9StandardStateModelV1, ThermochemErrorV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K,
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
        self.caloric_change(initial, final_temperature, false)
            .map(MolarEnthalpyV1::new)
    }

    /// Molar internal-energy change under the retained ideal-gas EOS, J/mol.
    ///
    /// Integrates Cv = Cp - R using the original NASA-9 coefficients, canceling
    /// the formation reference before rounding. The constant Cp/R coefficient
    /// is reduced by one before integration, so small positive Cv is not lost
    /// by subtracting separately rounded Delta h and R Delta T. Raw reference
    /// jumps and explicit caloric reconciliation have the same semantics as
    /// `molar_enthalpy_change`. Composition and source conventions stay fixed.
    pub fn molar_internal_energy_change(
        &self,
        initial: Temperature,
        final_temperature: Temperature,
    ) -> Result<MolarInternalEnergyV1, ThermochemErrorV1> {
        self.caloric_change(initial, final_temperature, true)
            .map(MolarInternalEnergyV1::new)
    }

    fn caloric_change(
        &self,
        initial: Temperature,
        final_temperature: Temperature,
        internal: bool,
    ) -> Result<f64, ThermochemErrorV1> {
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
            return Ok(0.0);
        }
        let (lo, hi, sign) = if t1 > t0 { (t0, t1, 1.0) } else { (t1, t0, -1.0) };
        let (first, lower) = self.region_for_temperature(lo).expect("validated lower endpoint");
        let (last, upper) = self.region_for_temperature(hi).expect("validated upper endpoint");
        let difference = |mut lower: [f64; 9], mut upper: [f64; 9], lo, hi| {
            if internal {
                lower[2] -= 1.0;
                upper[2] -= 1.0;
            }
            enthalpy_difference(lower, upper, lo, hi)
        };
        let value = sign * if self.caloric_reference().is_some() && first != last {
            self.integrate_caloric_regions(first, last, lo, hi,
                |coefficients, start, end| difference(coefficients, coefficients, start, end))
        } else {
            difference(lower.coefficients, upper.coefficients, lo, hi)
        };
        if !value.is_finite() {
            return Err(ThermochemErrorV1::NonFiniteEvaluation {
                property: if internal { Nasa9PropertyV1::InternalEnergy } else { Nasa9PropertyV1::Enthalpy },
                bits: value.to_bits(),
            });
        }
        Ok(value)
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
        self.caloric_difference(initial, final_temperature, false)
            .map(|(enthalpy, _)| MolarEnthalpyV1::new(enthalpy))
    }

    /// Same frozen-mixture enthalpy change in J/kg, using its own molar mass.
    pub fn specific_enthalpy_change(
        &self,
        initial: Temperature,
        final_temperature: Temperature,
    ) -> Result<MassSpecificEnergyQuantityV1, FrozenIdealGasMixtureErrorV1> {
        self.specific_caloric_change(initial, final_temperature, false)
    }

    /// Frozen-composition molar internal-energy change, J/mol.
    /// Each species integrates its own Cv before the compensated reduction.
    pub fn molar_internal_energy_change(
        &self,
        initial: Temperature,
        final_temperature: Temperature,
    ) -> Result<MolarInternalEnergyV1, FrozenIdealGasMixtureErrorV1> {
        self.caloric_difference(initial, final_temperature, true)
            .map(|(energy, _)| MolarInternalEnergyV1::new(energy))
    }

    /// Same frozen-mixture internal-energy change in J/kg.
    /// No absolute formation energies or mixture pressure-work subtraction enter
    /// the reduction. This API does not account for changing species amounts.
    pub fn specific_internal_energy_change(
        &self,
        initial: Temperature,
        final_temperature: Temperature,
    ) -> Result<MassSpecificEnergyQuantityV1, FrozenIdealGasMixtureErrorV1> {
        self.specific_caloric_change(initial, final_temperature, true)
    }

    fn specific_caloric_change(
        &self,
        initial: Temperature,
        final_temperature: Temperature,
        internal: bool,
    ) -> Result<MassSpecificEnergyQuantityV1, FrozenIdealGasMixtureErrorV1> {
        let (energy, molar_mass) = self.caloric_difference(initial, final_temperature, internal)?;
        let value = energy / molar_mass;
        if !value.is_finite() {
            return Err(FrozenIdealGasMixtureErrorV1::NonFiniteEvaluation {
                property: caloric_property(internal),
                bits: value.to_bits(),
            });
        }
        Ok(MassSpecificEnergyQuantityV1::new(value))
    }

    fn caloric_difference(
        &self,
        initial: Temperature,
        final_temperature: Temperature,
        internal: bool,
    ) -> Result<(f64, f64), FrozenIdealGasMixtureErrorV1> {
        let mut energy = 0.0;
        let mut correction = 0.0;
        let mut molar_mass = 0.0;
        for (component, (model, &fraction)) in self.components().iter()
            .zip(self.mole_composition().fractions()).enumerate()
        {
            let change = model.caloric_change(initial, final_temperature, internal)
                .map_err(|source| FrozenIdealGasMixtureErrorV1::ComponentEvaluation {
                    component,
                    species: model.species().clone(),
                    source: Box::new(source),
                })?;
            accumulate(&mut energy, &mut correction, fraction * change);
            molar_mass += fraction * model.molar_mass().value();
        }
        if !energy.is_finite() {
            return Err(FrozenIdealGasMixtureErrorV1::NonFiniteEvaluation {
                property: caloric_property(internal),
                bits: energy.to_bits(),
            });
        }
        if !molar_mass.is_finite() || molar_mass <= 0.0 {
            return Err(FrozenIdealGasMixtureErrorV1::InvalidMixtureMolarMass { bits: molar_mass.to_bits() });
        }
        Ok((energy, molar_mass))
    }
}

fn caloric_property(internal: bool) -> FrozenMixturePropertyV1 {
    if internal { FrozenMixturePropertyV1::MolarInternalEnergy } else { FrozenMixturePropertyV1::MolarEnthalpy }
}
