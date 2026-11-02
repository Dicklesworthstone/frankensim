//! Coherent gas refraction with spatially varying molecular composition.
//!
//! Add the original molecular polarizabilities before the nonlinear phase law:
//! `L = sum_i N_i v_i = sum_i rho_i K_i`, `n² = (1+2L)/(1-L)`.
//! Partial densities include concentration variations without differentiating
//! normalized fractions or dividing by a vanishing total density. Scattering
//! still has its independent molecular cross-section owner. The same dilute,
//! constant-polarizability approximation as the parent applies; no gas equation
//! of state, optical dispersion, absorption or chemical species is inferred.
//! Mixture convention: Shevchenko et al. (2021), §3.3.1,
//! <https://doi.org/10.1016/j.ces.2020.116433>.

use super::thermochemical::ThermochemicalRayleighState;
use super::*;

/// Units of every partial density and its spatial derivatives in one field.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PartialDensityBasis {
    /// Partial mass density in kg/m³, not a mass fraction.
    MassKgPerM3,
    /// Partial molecule count per m³, not mol/m³ or a mole fraction.
    MoleculesPerM3,
}

/// One source-bound component in the original model's declared order.
#[derive(Clone, Debug, PartialEq)]
pub struct MolecularRefractionComponent {
    species: Option<String>,
    optical: MolecularRayleigh,
    molar_mass_kg_mol: f64,
    number_coefficient: f64,
    mass_coefficient: f64,
}

impl MolecularRefractionComponent {
    /// Canonical thermochemical species, or `None` for an ordered `GasSpec` row.
    pub fn species(&self) -> Option<&str> {
        self.species.as_deref()
    }
    pub const fn optical(&self) -> MolecularRayleigh {
        self.optical
    }
    pub const fn molar_mass_kg_mol(&self) -> f64 {
        self.molar_mass_kg_mol
    }
    /// `v_i` [m³/molecule] or `K_i` [m³/kg], prepared once from original inputs.
    pub const fn coefficient(&self, basis: PartialDensityBasis) -> f64 {
        match basis {
            PartialDensityBasis::MassKgPerM3 => self.mass_coefficient,
            PartialDensityBasis::MoleculesPerM3 => self.number_coefficient,
        }
    }
}

/// Prepared immutable component optics; evaluation allocates no memory.
#[derive(Clone, Debug, PartialEq)]
pub struct MolecularRefractionModel {
    components: Vec<MolecularRefractionComponent>,
    identity: ContentHash,
}

impl MolecularRefractionModel {
    /// Retain canonical species, component masses and optical assignments from
    /// the original accepted NASA-9 attachment. The reference state's source
    /// model identity is retained even though local concentrations may vary.
    /// New/foreign species require another model with their explicit sources.
    pub fn from_thermochemical(
        state: &ThermochemicalRayleighState,
    ) -> Result<Self, GasScatteringError> {
        let receipts = state
            .resolved()
            .thermodynamics()
            .receipt()
            .component_receipts();
        Self::prepare(
            state.rayleigh().model_identity(),
            state
                .optical_sources()
                .iter()
                .zip(receipts)
                .map(|(row, receipt)| {
                    (
                        Some(row.species.as_str()),
                        row.optical,
                        f64::from_bits(receipt.molar_mass_bits()),
                    )
                }),
        )
    }

    /// Retain the original ordered ideal-gas inputs after reproducing the full
    /// accepted state. This uses each original gas owner's molar mass, with no
    /// species-name preset or substitute mixture equation of state.
    pub fn from_mixture(
        model: &MolecularRayleighMixture,
        accepted: &GasState,
    ) -> Result<Self, GasScatteringError> {
        let state = model.bind_gas_state(accepted)?;
        let masses = model
            .specs
            .iter()
            .map(|spec| {
                GasState::try_new(spec, accepted.temperature, accepted.pressure)
                    .map(|gas| R_USSA_1976 / gas.specific_gas_constant)
                    .map_err(GasScatteringError::Gas)
            })
            .collect::<Result<Vec<_>, _>>()?;
        Self::prepare(
            state.model_identity(),
            model
                .optical
                .iter()
                .zip(masses)
                .map(|(&optical, mass)| (None, optical, mass)),
        )
    }

    fn prepare<'a>(
        source_model: ContentHash,
        rows: impl Iterator<Item = (Option<&'a str>, MolecularRayleigh, f64)>,
    ) -> Result<Self, GasScatteringError> {
        let mut components = Vec::new();
        let mut hash = DomainHasher::new("org.frankensim.material.molecular-refraction-model.v1");
        hash.update(source_model.as_bytes());
        for (species, optical, mass) in rows {
            if components.len() == MAX_GAS_MIXTURE_COMPONENTS || !mass.is_finite() || mass <= 0.0 {
                return Err(GasScatteringError::Components);
            }
            let number_coefficient = optical.lorentz_lorenz_volume_m3();
            let mass_coefficient = (number_coefficient / mass) * AVOGADRO_PER_MOL;
            finite_nonnegative(mass_coefficient)?;
            if number_coefficient > 0.0 && mass_coefficient == 0.0 {
                return Err(GasScatteringError::Unrepresentable);
            }
            optical.hash_into(&mut hash);
            hash.update(&mass.to_bits().to_le_bytes());
            components.push(MolecularRefractionComponent {
                species: species.map(str::to_owned),
                optical,
                molar_mass_kg_mol: mass,
                number_coefficient,
                mass_coefficient,
            });
        }
        if components.is_empty() {
            return Err(GasScatteringError::Components);
        }
        hash.update(&(components.len() as u64).to_le_bytes());
        Ok(Self {
            components,
            identity: hash.finalize(),
        })
    }

    /// Order is canonical species order for thermochemistry, original row order
    /// for `MolecularRayleighMixture`. Every field must retain this exact mapping.
    pub fn components(&self) -> &[MolecularRefractionComponent] {
        &self.components
    }
    pub const fn content_identity(&self) -> ContentHash {
        self.identity
    }

    /// Evaluate supplied nonnegative partial densities in the declared basis.
    /// An all-zero vector is vacuum; empty, missing or extra rows refuse.
    /// Density values describe the caller's actual frozen reconstruction, not
    /// another thermodynamic evaluation or an optical renormalization.
    pub fn refractive_index(
        &self,
        partial_densities: &[f64],
        basis: PartialDensityBasis,
        wavelength_nm: f64,
    ) -> Result<CompositionRefractiveIndex, GasScatteringError> {
        wavelength(wavelength_nm)?;
        if partial_densities.len() != self.components.len() {
            return Err(GasScatteringError::Components);
        }
        let mut l = 0.0;
        for (&density, component) in partial_densities.iter().zip(&self.components) {
            if !density.is_finite() || density < 0.0 {
                return Err(GasScatteringError::RefractivityDomain);
            }
            let coefficient = component.coefficient(basis);
            let term = density * coefficient;
            if density > 0.0 && coefficient > 0.0 && term == 0.0 {
                return Err(GasScatteringError::Unrepresentable);
            }
            l += term;
        }
        let (phase_index, refractivity, denominator) = phase_from_lorentz_lorenz(l)?;
        let first = 1.5 / (phase_index * denominator * denominator);
        let second = first * (first * (1.0 + 8.0 * l) / (3.0 * phase_index));
        Ok(CompositionRefractiveIndex {
            phase_index,
            refractivity,
            lorentz_lorenz: l,
            derivative_lorentz_lorenz: first,
            second_derivative_lorentz_lorenz: second,
        })
    }
}

/// Phase law and its dimensionless derivatives with respect to additive `L`.
/// Spatial consumers apply `grad(n)=n_L grad(L)` and
/// `Hess(n)=n_L Hess(L)+n_LL grad(L) grad(L)^T`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CompositionRefractiveIndex {
    pub phase_index: f64,
    pub refractivity: f64,
    pub lorentz_lorenz: f64,
    pub derivative_lorentz_lorenz: f64,
    pub second_derivative_lorentz_lorenz: f64,
}

#[cfg(test)]
mod tests;
