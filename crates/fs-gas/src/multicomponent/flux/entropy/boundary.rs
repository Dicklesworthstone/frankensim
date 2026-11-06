//! Fixed exterior entropy efforts for the original mixture energy chart.
//!
//! eta=-S has v_i=(h_i,ref/T-s_i)/M_i-u^2/(2T), v_m=u/T, v_E=-1/T
//! and psi=u*sum_i(R*rho_i/M_i). A prescribed exterior therefore supplies
//! physical entropy -v_b.F+psi_b. This is the exterior thermodynamic effort on
//! the ACTUAL numerical mass, momentum and energy flux, not the advected entropy
//! of a separately reconstructed boundary state.
//!
//! The state account subtracts a fixed alpha_i from each molar entropy.
//! In that SAME basis h_i,ref/T-s_i+alpha_i simplifies to
//! u_i,ref/T-Delta(s_i^0)+R*ln(T/T*)+R*ln(rho_i/rho_i*).
//! Original source interval operations remove formation and entropy offsets
//! before contraction. The account restores alpha_i on its actual transfers.
//! Full momentum/kinetic efforts and psi remain essential for moving exterior
//! gas. Stationary molecular-bath heat/chemical efforts alone are insufficient.

use super::{Cell, Dd, Error, MAX_COMPONENTS, MixtureGas, Primitive, R, Temperature, log_ratio, poll, positive};
use crate::multicomponent::MAX_SPECIES;

/// Cached for one prescribed exterior trace and one immutable entropy basis.
/// No source inversion, heap allocation or logarithm is needed per RK stage.
#[derive(Clone)]
pub(in crate::multicomponent) struct ExteriorEntropyEffort {
    species_count: usize,
    positive_species: [bool; MAX_SPECIES],
    variables: [Dd; MAX_COMPONENTS],
    potential: Dd,
}

impl ExteriorEntropyEffort {
    /// Compatibility of the immutable exterior entropy effort, independent of
    /// accepted transfer history. Window steps clone this exact cached basis.
    pub(in crate::multicomponent) fn same_effort(&self, other: &Self) -> bool {
        self.species_count == other.species_count
            && self.positive_species == other.positive_species
            && self.potential.hi.to_bits() == other.potential.hi.to_bits()
            && self.potential.lo.to_bits() == other.potential.lo.to_bits()
            && self.variables.iter().zip(&other.variables).all(|(left, right)| {
                left.hi.to_bits() == right.hi.to_bits()
                    && left.lo.to_bits() == right.lo.to_bits()
            })
    }

    pub(in crate::multicomponent) fn new(
        gas: &MixtureGas,
        exterior: Cell<'_>,
        primitive: Primitive,
        reference_temperature: f64,
        reference_density: &[f64; MAX_SPECIES],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(cancelled)?;
        let n = gas.species_count();
        if exterior.partial_densities.len() != n
            || !gas.admissible_cell(exterior)
            || !positive(primitive.temperature)
            || !primitive.velocity.is_finite()
            || !positive(reference_temperature)
        {
            return Err(Error::EntropyDomain("invalid prescribed exterior entropy state"));
        }
        let mut result = Self {
            species_count: n,
            positive_species: [false; MAX_SPECIES],
            variables: [Dd::ZERO; MAX_COMPONENTS],
            potential: Dd::ZERO,
        };
        let t = Dd::from_f64(primitive.temperature);
        let velocity = Dd::from_f64(primitive.velocity);
        result.variables[n] = velocity / t;
        result.variables[n + 1] = -(Dd::ONE / t);
        let kinetic = (Dd::from_f64(0.5) * velocity) * result.variables[n];
        let log_t = Dd::from_f64(log_ratio(primitive.temperature, reference_temperature));
        checked(result.variables[n])?;
        checked(result.variables[n + 1])?;
        checked(kinetic)?;
        let mut gas_density = Dd::ZERO;
        for (k, model) in gas.source().models().iter().enumerate() {
            poll(cancelled)?;
            let density = exterior.partial_densities[k];
            if density == 0.0 {
                continue;
            }
            if !positive(reference_density[k]) {
                return Err(Error::EntropyDomain("invalid exterior entropy reference density"));
            }
            result.positive_species[k] = true;
            let mass = Dd::from_f64(model.molar_mass().value());
            let sensible = Dd::from_f64(model.molar_internal_energy_change(
                Temperature::new(gas.domain().minimum_temperature),
                Temperature::new(primitive.temperature),
            )?.value());
            let ds = Dd::from_f64(model.molar_entropy_change(
                Temperature::new(reference_temperature),
                Temperature::new(primitive.temperature),
            )?.value());
            let mu = sensible / t - ds + Dd::from_f64(R) * log_t
                + Dd::from_f64(R) * Dd::from_f64(log_ratio(density, reference_density[k]));
            result.variables[k] = mu / mass - kinetic;
            checked(result.variables[k])?;
            gas_density = gas_density + (Dd::from_f64(R) / mass) * Dd::from_f64(density);
        }
        result.potential = velocity * gas_density;
        checked(result.potential)?;
        poll(cancelled)?;
        Ok(result)
    }

    /// Check every actual stage before opposite stage transfers can cancel.
    /// A common-absent species has no finite logarithmic effort and must carry
    /// exactly zero numerical flux; no artificial activity floor is inserted.
    pub(in crate::multicomponent) fn admit_flux(
        &self,
        flux: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        poll(cancelled)?;
        if flux.len() != self.species_count + 2 {
            return Err(Error::EntropyDomain("prescribed exterior entropy flux dimension"));
        }
        for (k, value) in flux.iter().copied().enumerate() {
            if !value.is_finite() {
                return Err(Error::EntropyDomain("non-finite prescribed exterior entropy flux"));
            }
            if k < self.species_count && !self.positive_species[k] && value != 0.0 {
                return Err(Error::EntropyDomain("nonzero exterior flux at absent species activity"));
            }
        }
        Ok(())
    }

    /// Contract the native signed conservative end transfers, with a positive
    /// sign at the left end and negative sign at the right. signed_measure is
    /// the SAME oriented time*area integral. All axes include their actual
    /// stage weights; independently rounded physical fluxes must not replace it.
    pub(in crate::multicomponent) fn entropy_input(
        &self,
        transfer: &[f64],
        signed_measure: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Dd, Error> {
        self.admit_flux(transfer, cancelled)?;
        if !signed_measure.is_finite() {
            return Err(Error::EntropyDomain("non-finite prescribed exterior entropy measure"));
        }
        let mut input = Dd::from_f64(signed_measure) * self.potential;
        for (variable, transferred) in self.variables.iter().zip(transfer) {
            input = input - *variable * Dd::from_f64(*transferred);
        }
        checked(input)?;
        poll(cancelled)?;
        Ok(input)
    }
}

fn checked(value: Dd) -> Result<(), Error> {
    if value.hi.is_finite() && value.lo.is_finite() && value.to_f64().is_finite() {
        Ok(())
    } else {
        Err(Error::EntropyDomain("unrepresentable prescribed exterior entropy effort"))
    }
}

#[cfg(test)]
mod tests;
