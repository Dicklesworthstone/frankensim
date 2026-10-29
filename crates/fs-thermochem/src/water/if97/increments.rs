//! Original-source caloric and volume differences at fixed pressure.
//!
//! Factoring differences of the original integer-power Gibbs monomials keeps
//! small heat and pressure work beside large absolute enthalpy references.
//! Both endpoints and the entire intervening path must belong to the same
//! released region. No phase transition, datum shift, or saturation entropy
//! convention is inferred here.

use super::{
    If97Error, If97Region, If97State, MassSpecificEnergyQuantityV1, MassSpecificThermalQuantityV1,
    MassSpecificVolume, R_SPECIFIC_J_KG_K, R1_PRESSURE_STAR_PA, R1_TEMPERATURE_STAR_K, R1_TERMS,
    R2_IDEAL, R2_PRESSURE_STAR_PA, R2_RESIDUAL, R2_TEMPERATURE_STAR_K,
};
use fs_qty::{Pressure, Temperature};

/// Finite changes of the original IF97 equation along one admitted isobar.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct If97IsobaricDifference {
    enthalpy_j_kg: f64,
    entropy_j_kg_k: f64,
    volume_m3_kg: f64,
    internal_energy_j_kg: f64,
}
impl If97IsobaricDifference {
    /// Sensible enthalpy change [J/kg], with constant source terms cancelled.
    pub const fn specific_enthalpy_change(&self) -> MassSpecificEnergyQuantityV1 {
        MassSpecificEnergyQuantityV1::new(self.enthalpy_j_kg)
    }
    /// Original-source entropy change [J/(kg K)], equal to integral Cp/T dT.
    pub const fn specific_entropy_change(&self) -> MassSpecificThermalQuantityV1 {
        MassSpecificThermalQuantityV1::new(self.entropy_j_kg_k)
    }
    /// Specific-volume change [m3/kg], retaining small thermal expansion.
    pub const fn specific_volume_change(&self) -> MassSpecificVolume {
        MassSpecificVolume::new(self.volume_m3_kg)
    }
    /// Internal-energy change [J/kg], from the same delta h - p delta v.
    pub const fn specific_internal_energy_change(&self) -> MassSpecificEnergyQuantityV1 {
        MassSpecificEnergyQuantityV1::new(self.internal_energy_j_kg)
    }
}

/// Evaluate changes before rounding them into absolute source properties.
///
/// The operation count is bounded by the retained 34 or 52 source monomials.
/// Reversing the path changes only the sign; equal endpoints return exact zero.
/// This is a mathematical difference of the original regional equation, not a
/// patch to either endpoint's rounded state. Consumers still own their actual
/// stored-state/inversion residual and phase-junction accounting.
pub fn isobaric_difference(
    pressure: Pressure,
    region: If97Region,
    from: Temperature,
    to: Temperature,
) -> Result<If97IsobaricDifference, If97Error> {
    // Retain all original state admission, including representability of the
    // endpoints; a finite difference cannot make an infinite source state valid.
    If97State::at_pressure_temperature(pressure, from, region)?;
    If97State::at_pressure_temperature(pressure, to, region)?;
    let (t0, t1, sign) = if from.value() <= to.value() {
        (from.value(), to.value(), 1.0)
    } else {
        (to.value(), from.value(), -1.0)
    };
    admit_isobaric_interval(pressure.value(), t0, t1, region)?;
    if t0 == t1 {
        return Ok(If97IsobaricDifference {
            enthalpy_j_kg: 0.0,
            entropy_j_kg_k: 0.0,
            volume_m3_kg: 0.0,
            internal_energy_j_kg: 0.0,
        });
    }
    let (p_star, t_star, shift, terms) = match region {
        If97Region::Region1 => (
            R1_PRESSURE_STAR_PA,
            R1_TEMPERATURE_STAR_K,
            1.222,
            &R1_TERMS[..],
        ),
        If97Region::Region2 => (
            R2_PRESSURE_STAR_PA,
            R2_TEMPERATURE_STAR_K,
            0.5,
            &R2_RESIDUAL[..],
        ),
    };
    let dt = t1 - t0;
    let tau0 = t_star / t0;
    let tau1 = t_star / t1;
    // Unlike subtracting the two rounded reduced temperatures, this retains
    // the physical change even at adjacent representable endpoint temperatures.
    let delta_tau = -tau0 * (dt / t1);
    let (y0, y1) = (tau0 - shift, tau1 - shift);
    let pi = pressure.value() / p_star;
    let x = if region == If97Region::Region1 {
        7.1 - pi
    } else {
        pi
    };
    let pressure_powers = super::powers::<33>(x, 0);
    let zero = if region == If97Region::Region1 { 43 } else { 1 };
    let powers = DifferencePowers::<61>::new(y0, y1, delta_tau, zero);
    let (mut dh, mut ds, mut dv) = (Sum::default(), Sum::default(), Sum::default());
    for term in terms {
        let i = term.i as usize;
        let k = (term.j + zero as i32) as usize;
        let coefficient = term.n * pressure_powers[i];
        let j = f64::from(term.j);
        let power_change = powers.delta[k];
        let derivative_change = powers.delta[k - 1];
        dh.add(coefficient * j * derivative_change);
        // tau*j*y^(j-1)-y^j = (j-1)*y^j + shift*j*y^(j-1).
        ds.add(coefficient * ((j - 1.0) * power_change + shift * j * derivative_change));
        if term.i != 0 {
            let product_change = match term.j {
                0 => dt,
                1 => -shift * dt,
                _ => dt * powers.upper[k] + t0 * power_change,
            };
            let pressure_sign = if region == If97Region::Region1 {
                -1.0
            } else {
                1.0
            };
            dv.add(
                pressure_sign
                    * term.n
                    * f64::from(term.i)
                    * pressure_powers[i - 1]
                    * product_change,
            );
        }
    }
    if region == If97Region::Region2 {
        let ideal = DifferencePowers::<11>::new(tau0, tau1, delta_tau, 7);
        for term in &R2_IDEAL {
            let j = f64::from(term.j);
            let k = (term.j + 7) as usize;
            dh.add(term.n * j * ideal.delta[k - 1]);
            ds.add(term.n * (j - 1.0) * ideal.delta[k]);
        }
    }
    let h = sign * R_SPECIFIC_J_KG_K * t_star * dh.value;
    let s = sign * R_SPECIFIC_J_KG_K * ds.value;
    let mut v = R_SPECIFIC_J_KG_K / p_star * dv.value;
    if region == If97Region::Region2 {
        v += R_SPECIFIC_J_KG_K * dt / pressure.value();
    }
    v *= sign;
    let u = h - pressure.value() * v;
    if ![h, s, v, u].into_iter().all(f64::is_finite) {
        return Err(If97Error::Unrepresentable("isobaric property increment"));
    }
    Ok(If97IsobaricDifference {
        enthalpy_j_kg: h,
        entropy_j_kg_k: s,
        volume_m3_kg: v,
        internal_energy_j_kg: u,
    })
}

// At fixed pressure each piece of either released region is connected in T.
// Check both one-sided limits where R2 changes its pressure-bound equation;
// the independently fitted R4 and B23 equations are not silently equated.
pub(super) fn admit_isobaric_interval(
    p: f64,
    lo: f64,
    hi: f64,
    region: If97Region,
) -> Result<(), If97Error> {
    if !lo.is_finite() || !hi.is_finite() || lo > hi {
        return Err(If97Error::InvalidInput(
            "invalid isobaric temperature interval",
        ));
    }
    super::admit(p, lo, region)?;
    super::admit(p, hi, region)?;
    if region == If97Region::Region2 {
        if lo <= 623.15 && hi > 623.15 {
            super::admit(p, 623.15, region)?;
            if p > super::boundary23_pressure_raw(623.15) {
                return Err(If97Error::InvalidInput(
                    "isobar crosses the R4/B23 source-domain gap",
                ));
            }
        }
        if lo < 863.15 && hi >= 863.15 {
            super::admit(p, 863.15, region)?;
        }
    }
    Ok(())
}

// delta(y^n)=y1*delta(y^(n-1))+delta(y)*y0^(n-1). Every
// contribution has the same sign. Reciprocal bases cover negative powers,
// retaining delta(1/y) analytically. Build each power once per query instead
// of repeating a polynomial factorization for every source monomial.
struct DifferencePowers<const N: usize> {
    upper: [f64; N],
    delta: [f64; N],
}
impl<const N: usize> DifferencePowers<N> {
    fn new(a: f64, b: f64, change: f64, zero: usize) -> Self {
        let mut result = Self {
            upper: [1.0; N],
            delta: [0.0; N],
        };
        let mut from_power = 1.0;
        for k in zero + 1..N {
            result.upper[k] = b * result.upper[k - 1];
            result.delta[k] = b * result.delta[k - 1] + change * from_power;
            from_power *= a;
        }
        let (a, b, change) = (1.0 / a, 1.0 / b, -change / (a * b));
        from_power = 1.0;
        for k in (0..zero).rev() {
            result.upper[k] = b * result.upper[k + 1];
            result.delta[k] = b * result.delta[k + 1] + change * from_power;
            from_power *= a;
        }
        result
    }
}

#[derive(Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) {
        crate::entropy::accumulate(&mut self.value, &mut self.correction, value);
    }
}

#[cfg(test)]
mod tests;
