//! Explicit molecular dispersion anchored to the original optical reference.
//!
//! The source-declared normal-dispersion law is
//! `r(lambda) = r_ref + B (lambda^-2-lambda_ref^-2)
//!                       + C (lambda^-4-lambda_ref^-4)`, with wavelengths in µm.
//! This is a Cauchy fit, not a gas preset, resonance law or accuracy certificate.
//! Gas Cauchy fits: <https://doi.org/10.1016/0010-2180(81)90124-3>.
//! Group index differentiates the full Lorentz–Lorenz mixture, as required by
//! Ciddor and Hill, <https://doi.org/10.1364/AO.38.001663>.

use super::*;

/// A sourced Cauchy refractivity fit at the original molecular number density.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolecularCauchyDispersion {
    reference: MolecularRayleigh,
    b_um2: f64,
    c_um4: f64,
    wavelength_range_nm: [f64; 2],
    source_identity: ContentHash,
    reference_inverse_square: f64,
}

impl MolecularCauchyDispersion {
    /// Admit finite nonnegative shape coefficients and an explicit source band
    /// within 380–780 nm containing the original reference. Monotonic endpoint
    /// checks cover the dilute refractivity range throughout the source band.
    /// A source identity records the caller's fit; it does not authenticate it.
    pub fn try_new(
        reference: MolecularRayleigh,
        b_um2: f64,
        c_um4: f64,
        wavelength_range_nm: [f64; 2],
        source_identity: ContentHash,
    ) -> Result<Self, GasScatteringError> {
        let [lo, hi] = wavelength_range_nm;
        wavelength(lo)?;
        wavelength(hi)?;
        if !b_um2.is_finite()
            || b_um2 < 0.0
            || !c_um4.is_finite()
            || c_um4 < 0.0
            || lo >= hi
            || !(lo..=hi).contains(&reference.reference_wavelength_nm)
            || source_identity == ContentHash([0; 32])
        {
            return Err(GasScatteringError::OpticalReference);
        }
        let lambda = reference.reference_wavelength_nm * 1.0e-3;
        let result = Self {
            reference,
            b_um2: if b_um2 == 0.0 { 0.0 } else { b_um2 },
            c_um4: if c_um4 == 0.0 { 0.0 } else { c_um4 },
            wavelength_range_nm,
            source_identity,
            reference_inverse_square: 1.0 / (lambda * lambda),
        };
        for endpoint in wavelength_range_nm {
            let lambda = endpoint * 1.0e-3;
            let inverse_square = 1.0 / (lambda * lambda);
            result.volume(inverse_square, inverse_square * inverse_square)?;
        }
        Ok(result)
    }

    pub const fn reference(&self) -> MolecularRayleigh {
        self.reference
    }
    /// Shape coefficients `(B_um2, C_um4)`; the original reference fixes offset.
    pub const fn coefficients(&self) -> [f64; 2] {
        [self.b_um2, self.c_um4]
    }
    pub const fn wavelength_range_nm(&self) -> [f64; 2] {
        self.wavelength_range_nm
    }
    pub const fn source_identity(&self) -> ContentHash {
        self.source_identity
    }

    fn volume(
        &self,
        inverse_square: f64,
        inverse_fourth: f64,
    ) -> Result<MolecularRefractionCoefficient, GasScatteringError> {
        let r = self.reference.refractivity
            + self.b_um2 * (inverse_square - self.reference_inverse_square)
            + self.c_um4
                * (inverse_fourth - self.reference_inverse_square * self.reference_inverse_square);
        if !r.is_finite() || !(0.0..=MAX_DILUTE_REFRACTIVITY).contains(&r) {
            return Err(GasScatteringError::RefractivityDomain);
        }
        let n2_minus_one = r * (2.0 + r);
        let denominator = 3.0 + n2_minus_one;
        let coefficient = (n2_minus_one / denominator) / self.reference.reference_number_density_m3;
        let slope = -2.0 * self.b_um2 * inverse_square - 4.0 * self.c_um4 * inverse_fourth;
        let log_wavelength_derivative = ((6.0 * (1.0 + r) / denominator / denominator) * slope)
            / self.reference.reference_number_density_m3;
        if !coefficient.is_finite()
            || !log_wavelength_derivative.is_finite()
            || (r > 0.0 && coefficient == 0.0)
            || ((self.b_um2 > 0.0 || self.c_um4 > 0.0) && log_wavelength_derivative == 0.0)
        {
            return Err(GasScatteringError::Unrepresentable);
        }
        Ok(MolecularRefractionCoefficient {
            coefficient,
            log_wavelength_derivative,
        })
    }

    fn hash_into(&self, hash: &mut DomainHasher) {
        hash.update(self.source_identity.as_bytes());
        for value in [
            self.b_um2,
            self.c_um4,
            self.wavelength_range_nm[0],
            self.wavelength_range_nm[1],
        ] {
            hash.update(&value.to_bits().to_le_bytes());
        }
    }
}

/// Additive molecular volume or mass coefficient at one admitted wavelength.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolecularRefractionCoefficient {
    /// `v_i` [m³/molecule] or `K_i` [m³/kg], in the prepared density basis.
    pub coefficient: f64,
    /// `lambda d(coefficient)/d lambda`, with the same units as the coefficient.
    pub log_wavelength_derivative: f64,
}

impl MolecularRefractionCoefficient {
    const ZERO: Self = Self {
        coefficient: 0.0,
        log_wavelength_derivative: 0.0,
    };
}

/// Reusable stack-owned spectral preparation; no per-evaluation allocation.
/// Scalar phase/group, gradients and Hessians consume these same coefficients.
pub struct MolecularRefractionWavelength<'a> {
    model: &'a MolecularRefractionModel,
    coefficients: [MolecularRefractionCoefficient; MAX_GAS_MIXTURE_COMPONENTS],
    wavelength_nm: f64,
    basis: PartialDensityBasis,
}

impl MolecularRefractionModel {
    /// Attach one explicit spectral source to every original component, in its
    /// exact order. References must match completely; species are never inferred
    /// from an optical fit. Reattachment replaces spectral sources while keeping
    /// the same original model identity, so identity is independent of history.
    pub fn with_dispersion(
        &self,
        rows: &[MolecularCauchyDispersion],
    ) -> Result<Self, GasScatteringError> {
        if rows.len() != self.components.len() {
            return Err(GasScatteringError::Components);
        }
        let mut result = self.clone();
        let mut range = [380.0_f64, 780.0_f64];
        let mut hash = DomainHasher::new("org.frankensim.material.molecular-cauchy-refraction.v1");
        hash.update(self.reference_identity.as_bytes());
        for (component, &row) in result.components.iter_mut().zip(rows) {
            if row.reference != component.optical {
                return Err(GasScatteringError::StateMismatch);
            }
            range[0] = range[0].max(row.wavelength_range_nm[0]);
            range[1] = range[1].min(row.wavelength_range_nm[1]);
            component.dispersion = Some(row);
            row.hash_into(&mut hash);
        }
        if range[0] >= range[1] {
            return Err(GasScatteringError::OpticalReference);
        }
        result.wavelength_range_nm = range;
        result.identity = hash.finalize();
        Ok(result)
    }

    /// Prepare all coefficients and analytic spectral derivatives once at a
    /// vacuum wavelength. Reciprocal wavelength powers are shared across rows.
    /// This reusable value can be retained across many spatial queries.
    pub fn at_wavelength(
        &self,
        wavelength_nm: f64,
        basis: PartialDensityBasis,
    ) -> Result<MolecularRefractionWavelength<'_>, GasScatteringError> {
        wavelength(wavelength_nm)?;
        if !(self.wavelength_range_nm[0]..=self.wavelength_range_nm[1]).contains(&wavelength_nm) {
            return Err(GasScatteringError::RefractivityDomain);
        }
        let lambda = wavelength_nm * 1.0e-3;
        let inverse_square = 1.0 / (lambda * lambda);
        let inverse_fourth = inverse_square * inverse_square;
        let mut coefficients = [MolecularRefractionCoefficient::ZERO; MAX_GAS_MIXTURE_COMPONENTS];
        for (out, component) in coefficients.iter_mut().zip(&self.components) {
            let Some(shape) = component
                .dispersion
                .filter(|s| s.b_um2 != 0.0 || s.c_um4 != 0.0)
            else {
                *out = MolecularRefractionCoefficient {
                    coefficient: component.coefficient(basis),
                    log_wavelength_derivative: 0.0,
                };
                continue;
            };
            let volume = shape.volume(inverse_square, inverse_fourth)?;
            *out = match basis {
                PartialDensityBasis::MoleculesPerM3 => volume,
                PartialDensityBasis::MassKgPerM3 => MolecularRefractionCoefficient {
                    coefficient: (volume.coefficient / component.molar_mass_kg_mol)
                        * AVOGADRO_PER_MOL,
                    log_wavelength_derivative: (volume.log_wavelength_derivative
                        / component.molar_mass_kg_mol)
                        * AVOGADRO_PER_MOL,
                },
            };
            if !out.coefficient.is_finite()
                || !out.log_wavelength_derivative.is_finite()
                || (volume.coefficient > 0.0 && out.coefficient == 0.0)
                || (volume.log_wavelength_derivative != 0.0 && out.log_wavelength_derivative == 0.0)
            {
                return Err(GasScatteringError::Unrepresentable);
            }
        }
        Ok(MolecularRefractionWavelength {
            model: self,
            coefficients,
            wavelength_nm,
            basis,
        })
    }
}

impl MolecularRefractionWavelength<'_> {
    pub const fn wavelength_nm(&self) -> f64 {
        self.wavelength_nm
    }
    pub const fn density_basis(&self) -> PartialDensityBasis {
        self.basis
    }
    /// Exact prepared coefficients, in the original model's component order.
    pub fn coefficients(&self) -> &[MolecularRefractionCoefficient] {
        &self.coefficients[..self.model.components.len()]
    }

    /// Evaluate density-dependent phase and group from the SAME molecular sum.
    /// `n_g = n - n_L sum_i density_i lambda dK_i/dlambda`; constituent group
    /// indices must not be substituted into the nonlinear mixing rule.
    pub fn refractive_index(
        &self,
        partial_densities: &[f64],
    ) -> Result<CompositionRefractiveIndex, GasScatteringError> {
        if partial_densities.len() != self.model.components.len() {
            return Err(GasScatteringError::Components);
        }
        let mut l = 0.0;
        let mut log_derivative_l = 0.0;
        for (&density, coefficient) in partial_densities.iter().zip(self.coefficients()) {
            if !density.is_finite() || density < 0.0 {
                return Err(GasScatteringError::RefractivityDomain);
            }
            let term = density * coefficient.coefficient;
            let slope = density * coefficient.log_wavelength_derivative;
            if (density > 0.0 && coefficient.coefficient > 0.0 && term == 0.0)
                || (density > 0.0 && coefficient.log_wavelength_derivative != 0.0 && slope == 0.0)
            {
                return Err(GasScatteringError::Unrepresentable);
            }
            l += term;
            log_derivative_l += slope;
        }
        let (phase_index, refractivity, denominator) = phase_from_lorentz_lorenz(l)?;
        let first = 1.5 / (phase_index * denominator * denominator);
        let second = first * (first * (1.0 + 8.0 * l) / (3.0 * phase_index));
        let group_index = phase_index - first * log_derivative_l;
        let group_refractivity = refractivity - first * log_derivative_l;
        if !group_index.is_finite() || group_index <= 0.0 || !group_refractivity.is_finite() {
            return Err(GasScatteringError::Unrepresentable);
        }
        Ok(CompositionRefractiveIndex {
            phase_index,
            refractivity,
            group_index,
            group_refractivity,
            lorentz_lorenz: l,
            derivative_lorentz_lorenz: first,
            second_derivative_lorentz_lorenz: second,
        })
    }
}

#[cfg(test)]
mod tests;
