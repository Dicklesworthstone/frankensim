//! Isentropic incremental elasticity from the original thermal potential.
//!
//! Per unit REFERENCE volume, let beta = d(sigma)/dT at fixed total strain
//! and c_epsilon = c_caloric - T*d2(Psi_elastic)/dT2. Then
//! ds = (c_epsilon/T) dT - beta:d(epsilon). Eliminating dT at ds=0 gives
//! C_S = C_T + (T/c_epsilon) beta (x) beta. The complete six-component mixed
//! derivative includes temperature-dependent shear and Poisson response;
//! an isotropic bulk-only correction is insufficient in a strained state.
//!
//! The caloric owner supplies its actual local derivative and reference mass
//! convention. Specific Cp, fixed-volume Cv, and total fixed-strain capacity
//! are NOT interchangeable inputs. This module adds the original elastic
//! contribution exactly once. No relaxation time or attenuation is inferred.
//!
//! This is a local reversible tangent in the original infinitesimal-strain,
//! reference-configuration model. It does not add geometric prestress terms,
//! finite-deformation acoustoelasticity, current-density rescaling, or evidence
//! that a source's static coefficients describe a particular ultrasonic band.

use super::{ResolvedThermoelastic, ThermoelasticError, finite};
use crate::Tangent6;
use crate::elastic::acoustic::ElasticAcousticState;
use fs_blake3::{ContentHash, DomainHasher};

/// Actual accepted mechanical and caloric inputs for one local linearization.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdiabaticThermoelasticInput {
    /// TOTAL strain, before subtracting source-integrated free thermal strain.
    pub total_strain_mandel: [f64; 6],
    /// Original inertial mass divided by reference volume [kg/m^3].
    /// Never replace this with a current density on an undeformed mass matrix.
    pub reference_density_kg_m3: f64,
    /// Local caloric dH/dT divided by reference volume [J/(m^3 K)].
    /// This EXCLUDES the elastic-field contribution. The coupled thermal owner
    /// must establish its meaning and source support; a generic Cp label alone
    /// does not establish this derivative for an arbitrary thermodynamic model.
    pub caloric_heat_capacity_j_m3_k: f64,
    /// Identity of the original mechanical material/reference inventory.
    pub mechanical_material_identity: ContentHash,
    /// Identity of the actual accepted caloric state, including its chart.
    pub caloric_state_identity: ContentHash,
}

/// Source-retaining isentropic tangent usable by the existing wave/solid owners.
#[derive(Clone, Debug, PartialEq)]
pub struct AdiabaticThermoelasticState {
    source: ResolvedThermoelastic,
    input: AdiabaticThermoelasticInput,
    fixed_strain_heat_capacity_j_m3_k: f64,
    stress_temperature_derivative_mandel_pa_k: [f64; 6],
    temperature_gradient_mandel_k: [f64; 6],
    acoustic: ElasticAcousticState,
}

impl ResolvedThermoelastic {
    /// Linearize this same thermal potential at fixed local total entropy.
    ///
    /// Temperature, expansion, E and nu retain their original query/receipts.
    /// A nonpositive TOTAL fixed-strain capacity refuses rather than falling
    /// back to isothermal elasticity or treating supplied Cp as Cv.
    /// The caller must establish a smooth, fully solid caloric neighborhood
    /// at the same physical temperature and fixed external conditions.
    pub fn adiabatic_acoustics(
        &self,
        input: AdiabaticThermoelasticInput,
    ) -> Result<AdiabaticThermoelasticState, ThermoelasticError> {
        let invalid = ThermoelasticError::Invalid;
        if !input.reference_density_kg_m3.is_finite()
            || input.reference_density_kg_m3 <= 0.0
            || !input.caloric_heat_capacity_j_m3_k.is_finite()
            || input.caloric_heat_capacity_j_m3_k <= 0.0
            || input.mechanical_material_identity == ContentHash([0; 32])
            || input.caloric_state_identity == ContentHash([0; 32])
        {
            return Err(invalid(
                "adiabatic elasticity requires original nonzero identities and positive finite reference density/caloric capacity",
            ));
        }
        // At V0=1 m^3 the returned extensive elastic heat capacity is its
        // per-reference-volume coefficient. The original owner admits strain.
        let potential = self.state(input.total_strain_mandel, 1.0)?;
        let fixed_strain_heat_capacity_j_m3_k =
            finite(input.caloric_heat_capacity_j_m3_k + potential.heat_capacity_j_k)?;
        if fixed_strain_heat_capacity_j_m3_k <= 0.0 {
            return Err(invalid(
                "caloric plus elastic heat capacity at fixed total strain must be positive",
            ));
        }
        let beta = potential.stress_temperature_derivative_mandel_pa_k;
        let mut tangent = self.material.stiffness();
        let mut temperature_gradient_mandel_k = [0.0; 6];
        if beta.iter().any(|value| *value != 0.0) {
            let scale = finite(self.temperature_k / fixed_strain_heat_capacity_j_m3_k)?;
            if scale <= 0.0 {
                return Err(invalid("isentropic thermal compliance underflows"));
            }
            for i in 0..6 {
                temperature_gradient_mandel_k[i] = finite(scale * beta[i])?;
                for j in 0..=i {
                    // Evaluate each pair once so exact major symmetry is
                    // preserved without projecting or averaging evidence.
                    let correction = finite(temperature_gradient_mandel_k[i] * beta[j])?;
                    let value = finite(tangent[i][j] + correction)?;
                    tangent[i][j] = value;
                    tangent[j][i] = value;
                }
            }
        }
        let mut identity =
            DomainHasher::new("org.frankensim.fs-material.adiabatic-thermoelastic-state.v1");
        for hash in [
            self.elastic.identity(),
            self.expansion.identity(),
            input.mechanical_material_identity,
            input.caloric_state_identity,
        ] {
            identity.update(hash.as_bytes());
        }
        for value in input
            .total_strain_mandel
            .into_iter()
            .chain([
                input.reference_density_kg_m3,
                input.caloric_heat_capacity_j_m3_k,
                self.maximum_absolute_mandel_strain,
                self.temperature_k,
                self.free_linear_strain,
                self.youngs_slope_pa_k,
                self.poisson_slope_per_k,
                self.expansion_per_k,
                self.expansion_slope_per_k2,
                fixed_strain_heat_capacity_j_m3_k,
            ])
            .chain(beta)
            .chain(tangent.into_iter().flatten())
        {
            identity.update(&value.to_bits().to_le_bytes());
        }
        let acoustic = ElasticAcousticState::try_new_mandel(
            tangent,
            input.reference_density_kg_m3,
            identity.finalize(),
            self.maximum_absolute_mandel_strain,
        )?;
        Ok(AdiabaticThermoelasticState {
            source: self.clone(),
            input,
            fixed_strain_heat_capacity_j_m3_k,
            stress_temperature_derivative_mandel_pa_k: beta,
            temperature_gradient_mandel_k,
            acoustic,
        })
    }
}

impl AdiabaticThermoelasticState {
    /// Original source-backed thermal state; none of its receipts are replaced.
    #[must_use]
    pub const fn source(&self) -> &ResolvedThermoelastic {
        &self.source
    }

    /// Exact reference inventory, total strain and caloric inputs.
    #[must_use]
    pub const fn input(&self) -> &AdiabaticThermoelasticInput {
        &self.input
    }

    /// Actual TOTAL heat capacity at fixed total strain [J/(m^3 K)].
    #[must_use]
    pub const fn fixed_strain_heat_capacity_j_m3_k(&self) -> f64 {
        self.fixed_strain_heat_capacity_j_m3_k
    }

    /// Mixed thermal derivative beta = d(sigma)/dT [Pa/K].
    #[must_use]
    pub const fn stress_temperature_derivative_mandel_pa_k(&self) -> &[f64; 6] {
        &self.stress_temperature_derivative_mandel_pa_k
    }

    /// Infinitesimal adiabatic dT/d(epsilon_Mandel) [K].
    /// Contract with an incremental strain; this is not a finite state update.
    #[must_use]
    pub const fn temperature_gradient_mandel_k(&self) -> &[f64; 6] {
        &self.temperature_gradient_mandel_k
    }

    /// Original fixed-temperature tangent, without adiabatic augmentation [Pa].
    #[must_use]
    pub fn isothermal_stiffness_mandel_pa(&self) -> Tangent6 {
        self.source.material.stiffness()
    }

    /// Isentropic tangent and reference density for native pressure/wave queries.
    ///
    /// Hydrostatic strain/work returned by this tangent are INCREMENTAL about
    /// the supplied reference state, not a new finite-pressure equation of state.
    #[must_use]
    pub const fn acoustic(&self) -> &ElasticAcousticState {
        &self.acoustic
    }
}
