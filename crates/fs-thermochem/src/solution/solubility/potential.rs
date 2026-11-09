//! Isothermal chemical work with the original pure solid as the reference.
//!
//! At saturation the aqueous and solid chemical potentials agree. The ideal
//! activity relation therefore gives `mu_aq - mu_s = RT ln(c/c_sat)` and,
//! at fixed solvent volume, `F - n_total*mu_s = RT n_aq [ln(c/c_sat)-1]`.
//! The omitted solvent/solid datum cancels in a closed isothermal system with
//! identical pair laws. It must not be discarded for open species transfer,
//! changing temperature, or a different solid phase. These are chemical-work
//! primitives, not an additional heat-of-solution source.
//!
//! Activity/chemical-potential convention:
//! <https://doi.org/10.1351/goldbook.A00115>.

use super::{IdealSolubilityV1, SolubilityErrorV1, positive};
use crate::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use fs_math::det::ln;
use fs_qty::{Amount, AmountConcentration, Temperature, Volume};

impl IdealSolubilityV1 {
    /// Aqueous-minus-solid chemical potential [J/mol] at the actual temperature.
    ///
    /// `None` is the negative-infinite ideal limit at an exactly absent aqueous
    /// species. It is never replaced by a concentration floor. Supersaturation
    /// has positive chemical potential and favors growth of admitted solid;
    /// this diagnostic supplies neither a surface nor nucleation kinetics.
    pub fn solid_relative_chemical_potential_j_mol(
        &self,
        temperature: Temperature,
        concentration: AmountConcentration,
    ) -> Result<Option<f64>, SolubilityErrorV1> {
        let saturation = self.saturation_concentration(temperature)?.value();
        let concentration = concentration.value();
        admit_amount(concentration)?;
        if concentration == 0.0 {
            return Ok(None);
        }
        Ok(Some(finite(
            (R * temperature.value()) * log_ratio(concentration, saturation),
        )?))
    }

    /// Solid-relative free energy [J] of the actual dissolved inventory.
    ///
    /// The solid contribution is zero in this reference, regardless of its
    /// amount. Changes may be summed over closed compartments sharing this
    /// original pair law and temperature, provided dissolved-plus-solid moles
    /// are separately conserved. This is not absolute Gibbs energy or a
    /// caloric state function across different temperatures.
    pub fn solid_relative_free_energy_j(
        &self,
        temperature: Temperature,
        solvent_volume: Volume,
        dissolved: Amount,
    ) -> Result<f64, SolubilityErrorV1> {
        self.solid_relative_free_energy_change_j(
            temperature,
            solvent_volume,
            Amount::new(0.0),
            dissolved,
        )
    }

    /// Fixed-volume, isothermal solid-relative free-energy change [J].
    ///
    /// Retains the `-RT ln(c_sat)` reference term: precipitation/dissolution
    /// changes dissolved moles, so a dissolved-only mixing-energy difference
    /// cannot cancel it. Both endpoint amounts may be zero; positive aqueous
    /// traces remain explicit. Nearby endpoints use the integral of the
    /// chemical potential about their compensated midpoint, avoiding
    /// subtraction of large energies and retaining the quadratic saturation
    /// minimum. Work is bounded, deterministic and allocation-free.
    pub fn solid_relative_free_energy_change_j(
        &self,
        temperature: Temperature,
        solvent_volume: Volume,
        before_dissolved: Amount,
        after_dissolved: Amount,
    ) -> Result<f64, SolubilityErrorV1> {
        let saturation = self.saturation_concentration(temperature)?.value();
        if !positive(solvent_volume.value()) {
            return Err(SolubilityErrorV1::Invalid(
                "positive solvent volume for solubility chemical work",
            ));
        }
        let capacity = saturation * solvent_volume.value();
        if !positive(capacity) {
            return Err(SolubilityErrorV1::Unrepresentable(
                "positive saturation capacity for solubility chemical work",
            ));
        }
        let (a, b) = (before_dissolved.value(), after_dissolved.value());
        admit_amount(a)?;
        admit_amount(b)?;
        if a == b {
            return Ok(0.0);
        }
        let (low, high) = (a.min(b), a.max(b));
        let delta = high - low;
        let relative = delta / low;
        let rt = finite(R * temperature.value())?;
        let change = if low == 0.0 {
            energy_product(rt, high, log_ratio(high, capacity) - 1.0)?
        } else if relative <= 1.0 {
            // The midpoint is stored as a high/low pair. When adjacent amounts
            // straddle its rounding boundary, the low part retains the leading
            // chemical work (including departure from exact saturation).
            let half_delta = 0.5 * delta;
            let middle = low + half_delta;
            let remainder = (low - middle) + half_delta;
            let chemical = if delta < f64::MIN_POSITIVE {
                // Scaling first retains a midpoint below the least representable
                // amount; half_delta itself can round to zero in this case.
                log_ratio(low, capacity) + relative_log(0.5 * relative)
            } else {
                log_ratio(middle, capacity) + relative_log(remainder / middle)
            };
            let z = relative / (2.0 + relative);
            let square = z * z;
            let mut power = square;
            let mut correction = 0.0;
            // Average log(1+x) on [-z,z] is
            // -sum_{j>=1} z^(2j)/(2j(2j+1)); |z|<=1/3.
            for j in 1..24 {
                correction += power / f64::from(2 * j * (2 * j + 1));
                power *= square;
            }
            energy_product(rt, delta, chemical - correction)?
        } else {
            // Use the larger endpoint as logarithm reference; neither high/low
            // nor a quotient of dimensional concentrations must be representable.
            energy_product(rt, delta, log_ratio(high, capacity) - 1.0)?
                - energy_product(rt, low, log_ratio(low, high))?
        };
        let signed = if b > a { change } else { -change };
        finite(signed)
    }
}

fn admit_amount(value: f64) -> Result<(), SolubilityErrorV1> {
    if value.is_finite() && value >= 0.0 {
        Ok(())
    } else {
        Err(SolubilityErrorV1::Invalid(
            "finite nonnegative aqueous amount or concentration required",
        ))
    }
}

fn finite(value: f64) -> Result<f64, SolubilityErrorV1> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(SolubilityErrorV1::Unrepresentable(
            "solid-relative solubility chemical work",
        ))
    }
}

// A tiny molar work can underflow before multiplication by RT even though the
// final joules are representable. Try other deterministic multiplication orders
// only when the ordinary order loses range; genuine underflow remains zero.
fn energy_product(rt: f64, amount: f64, factor: f64) -> Result<f64, SolubilityErrorV1> {
    let ordinary = (amount * factor) * rt;
    if ordinary.is_finite() && (ordinary != 0.0 || amount == 0.0 || factor == 0.0) {
        return Ok(ordinary);
    }
    for alternative in [(rt * amount) * factor, (rt * factor) * amount] {
        if alternative.is_finite() && alternative != 0.0 {
            return Ok(alternative);
        }
    }
    finite(ordinary)
}

// Same deterministic atanh identity as the aqueous owner's relative log.
// The bounded nearby branch also handles adjacent and subnormal inputs.
fn log_ratio(numerator: f64, denominator: f64) -> f64 {
    let relative = (numerator - denominator) / denominator;
    if (-0.5..=1.0).contains(&relative) {
        relative_log(relative)
    } else {
        ln(numerator) - ln(denominator)
    }
}

fn relative_log(relative: f64) -> f64 {
    let z = relative / (2.0 + relative);
    let square = z * z;
    let mut term = z;
    let mut result = z;
    for k in 1..24 {
        term *= square;
        result += term / f64::from(2 * k + 1);
    }
    2.0 * result
}

#[cfg(test)]
mod tests;
