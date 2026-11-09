//! Frozen isothermal mechanics of the original source thermal potential.
//!
//! The exact current C(T), integrated free strain and native total/elastic
//! Mandel component domain remain with one source response. Reference density
//! belongs to the original inventory and is not refreshed from a heated-state
//! density query. This conversion adds neither an adiabatic tangent correction
//! nor geometric prestress. Temperature evolution and reversible thermal work
//! remain with the existing thermal potential and its coupled owner.

use super::{
    IsotropicThermalModuli, ResolvedThermoelastic, ResolvedThermoelasticLaw,
    ThermoelasticError, isotropic_thermal_moduli, tensor,
};
use crate::Tangent6;
use crate::elastic::acoustic::ElasticAcousticState;
use fs_blake3::{ContentHash, DomainHasher};

/// Original reference inventory for the fixed-geometry mechanical operator.
/// The caller retains the source/reference-volume meaning of these inputs.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IsothermalThermoelasticInput {
    /// Original inertial mass divided by reference volume [kg/m^3].
    /// Current thermal expansion does not create or remove reference mass.
    pub reference_density_kg_m3: f64,
    /// Identity of the original mechanical material/reference inventory.
    pub mechanical_material_identity: ContentHash,
}

/// One native thermal response and its exact isothermal mechanical view.
///
/// `acoustic()` acts on ELASTIC strain. Consumers must subtract the integrated
/// free strain from `source()` and admit BOTH total and elastic strains against
/// the original component box. `source().state(total, volume)` remains the
/// original owner of stress, entropy, internal energy and thermal capacity.
#[derive(Clone, Debug, PartialEq)]
pub struct IsothermalThermoelasticState {
    source: ResolvedThermoelastic,
    input: IsothermalThermoelasticInput,
    acoustic: ElasticAcousticState,
}

/// Exact local temperature derivatives in the original fixed Mandel frame.
///
/// The independent elastic engineering sources are affine; C(T) generally
/// is not. Original CTE tables may be piecewise affine, but this classical
/// second-derivative interface excludes their slope knots.
/// These derivatives retain the native rational/compliance-inverse response.
/// They do not differentiate a solved equilibrium, a geometric transformation,
/// a stabilization coefficient or a caloric transport operator.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermoelasticTemperatureDerivatives {
    /// dC/dT [Pa/K].
    pub stiffness_mandel_pa_k: Tangent6,
    /// d²C/dT² [Pa/K²].
    pub stiffness_mandel_pa_k2: Tangent6,
    /// d(epsilon_free)/dT = instantaneous expansion alpha [1/K].
    pub free_strain_mandel_per_k: [f64; 6],
    /// d²(epsilon_free)/dT² = d(alpha)/dT [1/K²].
    pub free_strain_mandel_per_k2: [f64; 6],
    /// Source-declared isotropic (dK/dT, dG/dT) [Pa/K].
    /// A generic or orthotropic tensor never infers this metadata.
    pub isotropic_moduli_pa_k: Option<(f64, f64)>,
    /// Source-declared isotropic (d²K/dT², d²G/dT²) [Pa/K²].
    pub isotropic_moduli_pa_k2: Option<(f64, f64)>,
}

impl ResolvedThermoelastic {
    fn isotropic_thermal_modes(
        &self,
    ) -> Result<Option<IsotropicThermalModuli>, ThermoelasticError> {
        match &self.response {
            ResolvedThermoelasticLaw::Isotropic(response) => isotropic_thermal_moduli(
                response.material.youngs,
                response.material.poisson,
                response.youngs_slope_pa_k,
                response.poisson_slope_per_k,
            )
            .map(Some),
            ResolvedThermoelasticLaw::Tensor(response) => Ok(response.isotropic_moduli),
        }
    }

    /// Expose the same exact coefficient derivatives used by `state`.
    /// No source query, finite difference, unit-strain probing or fit occurs.
    /// A continuous CTE slope knot refuses this classical second-derivative
    /// interface; its finite energy exchange remains available.
    pub fn temperature_derivatives(
        &self,
    ) -> Result<ThermoelasticTemperatureDerivatives, ThermoelasticError> {
        if !self.temperature_is_smooth {
            return Err(ThermoelasticError::Invalid(
                "classical thermal second derivative is undefined at a CTE slope knot"));
        }
        self.source_side_temperature_derivatives()
    }

    // Finite exchange needs C_T and alpha, both continuous at a CTE knot.
    // Its explicit piecewise source integral supplies the finite increments;
    // no classical second derivative is inferred from this one-sided snapshot.
    pub(super) fn source_side_temperature_derivatives(
        &self,
    ) -> Result<ThermoelasticTemperatureDerivatives, ThermoelasticError> {
        let moduli = self.isotropic_thermal_modes()?;
        let (first, second, alpha, alpha_first) = match &self.response {
            ResolvedThermoelasticLaw::Isotropic(response) => {
                let modes = moduli.expect("native isotropic source retains its moduli");
                (
                    tensor::bulk_shear_matrix(modes.bulk_first, modes.shear_first)?,
                    tensor::bulk_shear_matrix(modes.bulk_second, modes.shear_second)?,
                    [
                        response.expansion_per_k,
                        response.expansion_per_k,
                        response.expansion_per_k,
                        0.0,
                        0.0,
                        0.0,
                    ],
                    [
                        response.expansion_slope_per_k2,
                        response.expansion_slope_per_k2,
                        response.expansion_slope_per_k2,
                        0.0,
                        0.0,
                        0.0,
                    ],
                )
            }
            ResolvedThermoelasticLaw::Tensor(response) => (
                response.stiffness_derivative_mandel_pa_k,
                response.stiffness_second_derivative_mandel_pa_k2,
                response.expansion_mandel_per_k,
                response.expansion_derivative_mandel_per_k2,
            ),
        };
        if first
            .iter()
            .flatten()
            .chain(second.iter().flatten())
            .chain(alpha.iter())
            .chain(alpha_first.iter())
            .any(|value| !value.is_finite())
        {
            return Err(ThermoelasticError::Invalid(
                "native thermal coefficient derivatives are unrepresentable",
            ));
        }
        Ok(ThermoelasticTemperatureDerivatives {
            stiffness_mandel_pa_k: first,
            stiffness_mandel_pa_k2: second,
            free_strain_mandel_per_k: alpha,
            free_strain_mandel_per_k2: alpha_first,
            isotropic_moduli_pa_k: moduli.map(|modes| (modes.bulk_first, modes.shear_first)),
            isotropic_moduli_pa_k2: moduli.map(|modes| (modes.bulk_second, modes.shear_second)),
        })
    }

    /// Bind this exact isothermal source response to original reference inertia.
    ///
    /// Explicit native K/G metadata survives E/nu, K/G and E/G source schemas,
    /// including an independently oriented expansion tensor. The dense C(T)
    /// values are copied unchanged; no modulus round trip or isotropic fit is
    /// performed. Both isotropic and tensor sources keep their Mandel component
    /// domain, rather than inheriting another law's Frobenius norm domain.
    pub fn isothermal_elasticity(
        &self,
        input: IsothermalThermoelasticInput,
    ) -> Result<IsothermalThermoelasticState, ThermoelasticError> {
        if !input.reference_density_kg_m3.is_finite()
            || input.reference_density_kg_m3 <= 0.0
            || input.mechanical_material_identity == ContentHash([0; 32])
        {
            return Err(ThermoelasticError::Invalid(
                "isothermal elasticity requires positive finite reference density and a nonzero inventory identity",
            ));
        }
        let stiffness = self.stiffness_mandel_pa();
        let moduli = self.isotropic_thermal_modes()?;
        let mut identity =
            DomainHasher::new("org.frankensim.fs-material.isothermal-thermoelastic-state.v1");
        for hash in [
            self.elastic.identity(),
            self.expansion.identity(),
            input.mechanical_material_identity,
        ] {
            identity.update(hash.as_bytes());
        }
        for value in [
            input.reference_density_kg_m3,
            self.maximum_absolute_mandel_strain,
            self.temperature_k,
        ]
        .into_iter()
        .chain(self.thermal_identity_values())
        .chain(stiffness.into_iter().flatten())
        {
            identity.update(&value.to_bits().to_le_bytes());
        }
        identity.update(&[u8::from(moduli.is_some())]);
        if let Some(modes) = moduli {
            for value in [
                modes.bulk,
                modes.shear,
                modes.bulk_first,
                modes.shear_first,
                modes.bulk_second,
                modes.shear_second,
            ] {
                identity.update(&value.to_bits().to_le_bytes());
            }
        }
        let identity = identity.finalize();
        let acoustic = if let Some(modes) = moduli {
            ElasticAcousticState::from_isotropic_parts(
                stiffness,
                modes.bulk,
                modes.shear,
                input.reference_density_kg_m3,
                identity,
                self.maximum_absolute_mandel_strain,
            )?
        } else {
            ElasticAcousticState::try_new_mandel(
                stiffness,
                input.reference_density_kg_m3,
                identity,
                self.maximum_absolute_mandel_strain,
            )?
        };
        Ok(IsothermalThermoelasticState {
            source: self.clone(),
            input,
            acoustic,
        })
    }
}

impl IsothermalThermoelasticState {
    /// Complete original source response, integrated expansion and potential.
    #[must_use]
    pub const fn source(&self) -> &ResolvedThermoelastic {
        &self.source
    }

    /// Original reference density and material/inventory identity.
    #[must_use]
    pub const fn input(&self) -> &IsothermalThermoelasticInput {
        &self.input
    }

    /// Exact isothermal C(T), original reference density and elastic-strain box.
    #[must_use]
    pub const fn acoustic(&self) -> &ElasticAcousticState {
        &self.acoustic
    }
}
