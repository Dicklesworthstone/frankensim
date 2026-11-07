//! Owned interval coefficients for the original reciprocal thermal work account.
//!
//! Preparation resolves original sources and integrates thermal expansion once.
//! Constitutive trials below contain no source query, allocation, clock advance
//! or return map. The tensor material remains the owner of physical and numerical
//! dissipation; this account only supplies coefficient work and thermal entropy.

use super::{ThermalHistoryContribution, ThermoplasticIncrement};
use crate::plastic::energy::{EnergyError, EnergyIncrement, EnergyTolerances};
use crate::plastic::thermal::{FeedbackError, finite, product};
use crate::plastic::voce::VoceJ2;
use crate::plastic::{J2Plasticity, J2State};
use crate::tensor::{contract, trace};
use crate::{SmallStrainLaw, Voigt};

/// Source-admitted endpoint coefficients, independent of any specimen history.
/// Fields are writable only by the original source's preparation method.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct PreparedThermoplasticEnergy {
    pub(super) old_temperature_k: f64,
    pub(super) temperature_k: f64,
    pub(super) old_expansion_strain: f64,
    pub(super) expansion_strain: f64,
    pub(super) delta_expansion_strain: f64,
    pub(super) old_expansion_per_k: f64,
    pub(super) expansion_per_k: f64,
    pub(super) variable_expansion: bool,
    pub(super) old_law: J2Plasticity,
    pub(super) law: J2Plasticity,
    pub(super) youngs_temperature_slope_pa_k: f64,
    pub(super) hardening_temperature_slope_pa_k: f64,
    pub(super) kinematic_hardening_temperature_slope_pa_k: f64,
    pub(super) voce: Option<VoceJ2>,
    pub(super) saturation_temperature_slope_pa_k: f64,
    pub(super) old_saturation_pa: f64,
    pub(super) history_variable_storage: bool,
}

impl PreparedThermoplasticEnergy {
    /// Finish a native increment at this prepared old/new temperature pair.
    /// Additional histories supply only their coefficient work and entropy;
    /// their recoverable energy and losses already belong to `mechanical`.
    #[allow(clippy::too_many_arguments)]
    pub(in crate::plastic::thermal) fn finish<E>(
        &self,
        old_strain: &Voigt,
        old_state: &J2State,
        strain: &Voigt,
        mechanical: EnergyIncrement,
        history: ThermalHistoryContribution,
        tolerances: EnergyTolerances,
    ) -> Result<ThermoplasticIncrement, FeedbackError<E>> {
        if history.variable_storage != self.history_variable_storage {
            return Err(FeedbackError::Invalid(
                "thermal history changed its prepared storage convention",
            ));
        }
        let old_mechanical = mechanical_strain(
            old_strain,
            self.old_expansion_strain,
            self.old_law.elastic.strain_limit,
        )?;
        mechanical_strain(strain, self.expansion_strain, self.law.elastic.strain_limit)?;
        let old_elastic = core::array::from_fn(|i| old_mechanical[i] - old_state.plastic_strain[i]);
        let variable_moduli = self.youngs_temperature_slope_pa_k != 0.0
            || self.hardening_temperature_slope_pa_k != 0.0
            || self.kinematic_hardening_temperature_slope_pa_k != 0.0;
        let saturation_slope = self.saturation_temperature_slope_pa_k;
        let variable_storage =
            variable_moduli || saturation_slope != 0.0 || history.variable_storage;
        // Inspect the old elastic stress at its accepted temperature. This
        // never return-maps or otherwise advances the accepted plastic history.
        let old_stress = self.old_law.elastic.stress(&old_elastic, &());
        let delta_strain = core::array::from_fn(|i| strain[i] - old_strain[i]);
        let mechanical_work_j_m3 = finite(contract(&mechanical.stress_pa, &delta_strain))?;
        let mut old_entropy = product(self.old_expansion_per_k, trace(&old_stress))?;
        let mut new_entropy = product(self.expansion_per_k, trace(&mechanical.stress_pa))?;
        let mut coefficient_storage_change_j_m3 = if variable_moduli {
            let youngs = self.law.elastic.youngs;
            // Fixed nu makes C(T)=E(T)*C_unit. These derivatives use the native
            // frozen-endpoint elastic energies, preserving the original account.
            let elastic_derivative = self.youngs_temperature_slope_pa_k / youngs;
            let old_p_squared = product(old_state.alpha, old_state.alpha)?;
            let new_p_squared = product(mechanical.state.alpha, mechanical.state.alpha)?;
            old_entropy = finite(
                old_entropy
                    - product(elastic_derivative, mechanical.energy.initial.elastic_j_m3)?
                    - product(0.5 * self.hardening_temperature_slope_pa_k, old_p_squared)?,
            )?;
            new_entropy = finite(
                new_entropy
                    - product(
                        elastic_derivative,
                        mechanical.energy.final_state.elastic_j_m3,
                    )?
                    - product(0.5 * self.hardening_temperature_slope_pa_k, new_p_squared)?,
            )?;
            let mut coefficient = finite(
                product(
                    (youngs - self.old_law.elastic.youngs) / youngs,
                    mechanical.energy.initial.elastic_j_m3,
                )? + product(
                    0.5 * (self.law.hardening - self.old_law.hardening),
                    old_p_squared,
                )?,
            )?;
            if self.kinematic_hardening_temperature_slope_pa_k != 0.0 {
                // Native Prager storage remains defined as H_k reaches zero;
                // differentiate its strain shape without dividing by H_k.
                let old_plastic = crate::tensor::deviator(&old_state.plastic_strain);
                let new_plastic = crate::tensor::deviator(&mechanical.state.plastic_strain);
                let old_shape = finite(contract(&old_plastic, &old_plastic) / 3.0)?;
                let new_shape = finite(contract(&new_plastic, &new_plastic) / 3.0)?;
                old_entropy = finite(
                    old_entropy
                        - product(self.kinematic_hardening_temperature_slope_pa_k, old_shape)?,
                )?;
                new_entropy = finite(
                    new_entropy
                        - product(self.kinematic_hardening_temperature_slope_pa_k, new_shape)?,
                )?;
                coefficient = finite(
                    coefficient
                        + product(
                            self.law.kinematic_hardening - self.old_law.kinematic_hardening,
                            old_shape,
                        )?,
                )?;
            }
            coefficient
        } else {
            0.0
        };
        if saturation_slope != 0.0 {
            let law = self
                .voce
                .as_ref()
                .expect("selected thermal saturation has a native law");
            // Fixed-history unit-Q storage stays regular when endpoint Q is zero.
            let old_shape = finite(law.saturation_storage_per_pa(old_state.alpha))?;
            let new_shape = finite(law.saturation_storage_per_pa(mechanical.state.alpha))?;
            old_entropy = finite(old_entropy - product(saturation_slope, old_shape)?)?;
            new_entropy = finite(new_entropy - product(saturation_slope, new_shape)?)?;
            coefficient_storage_change_j_m3 = finite(
                coefficient_storage_change_j_m3
                    + product(law.saturation_pa() - self.old_saturation_pa, old_shape)?,
            )?;
        }
        if history.variable_storage {
            old_entropy = finite(old_entropy + history.old_entropy_j_m3_k)?;
            new_entropy = finite(new_entropy + history.new_entropy_j_m3_k)?;
            coefficient_storage_change_j_m3 =
                finite(coefficient_storage_change_j_m3 + history.coefficient_storage_change_j_m3)?;
        }
        let thermoelastic_storage_change_j_m3 = finite(
            product(self.temperature_k, new_entropy)?
                - product(self.old_temperature_k, old_entropy)?,
        )?;
        let reversible_heat_j_m3 = if !self.variable_expansion && !variable_storage {
            // Keep the constant-coefficient reduction and its original arithmetic.
            product(-self.old_temperature_k, finite(new_entropy - old_entropy)?)?
        } else {
            // This interval's exact expansion integral avoids subtracting two
            // large primitives relative to the stress-free reference temperature.
            finite(
                product(trace(&mechanical.stress_pa), self.delta_expansion_strain)?
                    - coefficient_storage_change_j_m3
                    - thermoelastic_storage_change_j_m3,
            )?
        };
        let stored_change_j_m3 =
            finite(mechanical.energy.stored_change_j_m3() + coefficient_storage_change_j_m3)?;
        let residual = finite(
            mechanical_work_j_m3
                - stored_change_j_m3
                - thermoelastic_storage_change_j_m3
                - mechanical.energy.physical_dissipation_j_m3
                - reversible_heat_j_m3
                - mechanical.energy.numerical_dissipation_j_m3,
        )?;
        if residual.abs() > tolerances.energy_j_m3 {
            return Err(FeedbackError::Mechanical(EnergyError::EnergyResidual {
                residual_j_m3: residual,
                tolerance_j_m3: tolerances.energy_j_m3,
            }));
        }
        Ok(ThermoplasticIncrement {
            mechanical,
            coefficient_storage_change_j_m3,
            mechanical_work_j_m3,
            reversible_heat_j_m3,
            thermoelastic_storage_change_j_m3,
        })
    }
}

/// Shared total/mechanical small-strain admission for borrowed and owned paths.
pub(super) fn mechanical_strain<E>(
    strain: &Voigt,
    expansion: f64,
    limit: f64,
) -> Result<Voigt, FeedbackError<E>> {
    let result = core::array::from_fn(|i| strain[i] - if i < 3 { expansion } else { 0.0 });
    for value in [strain, &result] {
        let norm_squared = contract(value, value);
        if value.iter().any(|v| !v.is_finite())
            || !norm_squared.is_finite()
            || norm_squared > limit * limit
        {
            return Err(FeedbackError::Invalid(
                "total or thermally offset strain exceeds the small-strain bound",
            ));
        }
    }
    Ok(result)
}
