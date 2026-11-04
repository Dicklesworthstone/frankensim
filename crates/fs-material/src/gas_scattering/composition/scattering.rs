//! Elastic molecular scattering from the SAME prepared dispersive optics.
//!
//! Independent molecular cross sections add: the coherent mixture's summed
//! polarizability must never be squared to obtain scattering. Each species uses
//! its original King factor and molecular volume at the prepared wavelength.
//! The phase is cross-section weighted, so different spectral polarizabilities
//! can change mixture depolarization even when each species' ratio is constant.
//! Cross section and King convention: Bodhaine et al. (1999),
//! <https://journals.ametsoc.org/view/journals/atot/16/11/1520-0426_1999_016_1854_orodc_2_0_co_2.xml>.
//! Matched phase convention: Lin et al. (2022), Eqs. 42–53,
//! <https://doi.org/10.3389/frsen.2022.880768>.

use super::*;

/// Local elastic scattering of the original partial-density reconstruction.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CompositionRayleighScattering {
    /// Original partial masses summed in kg/m³, including zero-polarizability gas.
    pub mass_density_kg_m3: f64,
    /// `sum_i N_i sigma_i`, in m⁻¹. This is scattering, never true absorption.
    pub scattering_per_m: f64,
    /// Scattering-weighted unpolarized 90° `I_parallel/I_perpendicular` ratio.
    /// Zero scattering uses the original material convention `delta = 0`.
    pub depolarization_ratio: f64,
}

impl MolecularRefractionWavelength<'_> {
    /// Evaluate exact molecular scattering and its matched effective phase.
    ///
    /// Inputs retain this preparation's original component order and density
    /// basis. Full phase/group/dilute admission is shared with `refractive_index`.
    /// Work is bounded by the same 128 components and uses no heap allocation.
    /// No wavelength interpolation or second polarizability fit is performed.
    ///
    /// `sigma_i = 24 pi³ F_i (v_i/lambda_m²)²`, with
    /// `F_i = (6+3 delta_i)/(6-7 delta_i)`. The mixture's isotropic fraction is
    /// `sum_i beta_i [3 delta_i/(2+delta_i)] / sum_i beta_i`. Its equivalent
    /// depolarization is `2 w/(3-w)`, not a composition average of the deltas.
    /// Species depolarization remains wavelength independent; anisotropy
    /// dispersion, correlated scattering, absorption and emission are absent.
    pub fn scattering(
        &self,
        partial_densities: &[f64],
    ) -> Result<CompositionRayleighScattering, GasScatteringError> {
        self.refractive_index(partial_densities)?;
        let lambda_m = self.wavelength_nm() * 1.0e-9;
        let lambda_squared = lambda_m * lambda_m;
        let mut contributions = [0.0; MAX_GAS_MIXTURE_COMPONENTS];
        let mut scattering_per_m = 0.0;
        let mut mass_density_kg_m3 = 0.0;
        for (i, (&density, component)) in partial_densities
            .iter()
            .zip(self.model.components())
            .enumerate()
        {
            if density == 0.0 {
                continue;
            }
            let mass = component.molar_mass_kg_mol();
            let local_mass = match self.density_basis() {
                PartialDensityBasis::MassKgPerM3 => density,
                // Divide the larger operand first: this avoids intermediate
                // overflow and avoidable underflow in N*M/N_A.
                PartialDensityBasis::MoleculesPerM3 => {
                    (density.max(mass) / AVOGADRO_PER_MOL) * density.min(mass)
                }
            };
            positive_result(local_mass)?;
            mass_density_kg_m3 += local_mass;
            let optical = component.optical();
            let sigma = if (380.0..=780.0).contains(&self.wavelength_nm())
                && component
                    .dispersion()
                    .is_none_or(|shape| shape.coefficients() == [0.0; 2])
            {
                // Keep original visible, nondispersive arithmetic. An explicit
                // zero-shape fit can also authorize a wider band; there the
                // same prepared molecular volume enters the formula below.
                optical.cross_section_m2(self.wavelength_nm())?
            } else {
                let volume = self.molecular_volumes_m3[i];
                let scaled = volume / lambda_squared;
                let delta = optical.depolarization_ratio();
                let king = (6.0 + 3.0 * delta) / (6.0 - 7.0 * delta);
                let pi = core::f64::consts::PI;
                let sigma = (24.0 * pi * pi * pi * king) * scaled * scaled;
                finite_nonnegative(sigma)?;
                if volume > 0.0 {
                    positive_result(sigma)?;
                }
                sigma
            };
            if sigma == 0.0 {
                continue;
            }
            let coefficient = match self.density_basis() {
                PartialDensityBasis::MoleculesPerM3 => sigma,
                PartialDensityBasis::MassKgPerM3 => {
                    let mut value = (sigma / mass) * AVOGADRO_PER_MOL;
                    if value == 0.0 {
                        // A representable result can survive a subnormal sigma
                        // even when the intermediate sigma/M rounds to zero.
                        value = (sigma * AVOGADRO_PER_MOL) / mass;
                    }
                    positive_result(value)?;
                    value
                }
            };
            let contribution = density * coefficient;
            positive_result(contribution)?;
            contributions[i] = contribution;
            scattering_per_m += contribution;
        }
        finite_nonnegative(mass_density_kg_m3)?;
        finite_nonnegative(scattering_per_m)?;
        let mut isotropic_weight = 0.0;
        if scattering_per_m > 0.0 {
            for (&contribution, component) in contributions.iter().zip(self.model.components()) {
                if contribution == 0.0 {
                    continue;
                }
                let weight = contribution / scattering_per_m;
                positive_result(weight)?;
                let delta = component.optical().depolarization_ratio();
                // Normalize before multiplying by delta, preserving the phase
                // even for an extremely dilute but representable scattering rate.
                let term = weight * (3.0 * delta / (2.0 + delta));
                if delta > 0.0 {
                    positive_result(term)?;
                }
                isotropic_weight += term;
            }
        }
        let depolarization_ratio = 2.0 * isotropic_weight / (3.0 - isotropic_weight);
        if !depolarization_ratio.is_finite() || !(0.0..6.0 / 7.0).contains(&depolarization_ratio) {
            return Err(GasScatteringError::Unrepresentable);
        }
        Ok(CompositionRayleighScattering {
            mass_density_kg_m3,
            scattering_per_m,
            depolarization_ratio,
        })
    }
}

fn positive_result(value: f64) -> Result<(), GasScatteringError> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(GasScatteringError::Unrepresentable)
    }
}

#[cfg(test)]
mod tests;
