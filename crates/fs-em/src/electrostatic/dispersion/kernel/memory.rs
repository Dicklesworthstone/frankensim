//! Original constitutive owners behind the shared charge/field evolution.
//!
//! This adapter contains no material update equation. Monomorphization keeps
//! Debye states at one vector per pole and Lorentz states at their original
//! polarization/current pair, without an enum or unused kinetic allocation.

use fs_material::field_storage::{
    FieldStorageError,
    dielectric::{
        bulk::ResolvedDebyeMaterial,
        lorentz::{LorentzState, PreparedLorentzMaterialMidpoint, ResolvedLorentzMaterial},
    },
};

type Result<T> = core::result::Result<T, FieldStorageError>;

/// Source extraction only: the shared field owner integrates Ohmic current
/// independently of the original relaxation or inertial polarization law.
#[cfg(feature = "resolved-leaky-dielectric")]
pub(crate) trait LeakyMemoryLaw<const D: usize> {
    type Memory: MemoryLaw<D>;

    fn into_parts(self) -> (Self::Memory, String, f64);
}

#[cfg(feature = "resolved-leaky-dielectric")]
impl<const D: usize> LeakyMemoryLaw<D>
    for fs_material::field_storage::dielectric::leaky::ResolvedLeakyDielectric
{
    type Memory = ResolvedDebyeMaterial;

    fn into_parts(self) -> (Self::Memory, String, f64) {
        Self::into_parts(self)
    }
}

#[cfg(feature = "resolved-leaky-dielectric")]
impl<const D: usize> LeakyMemoryLaw<D>
    for fs_material::field_storage::dielectric::lorentz::leaky::ResolvedLeakyLorentzMaterial
{
    type Memory = ResolvedLorentzMaterial;

    fn into_parts(self) -> (Self::Memory, String, f64) {
        Self::into_parts(self)
    }
}

pub(crate) struct Recovery<S, const D: usize> {
    pub states_after: Vec<S>,
    pub displacement_after_c_m2: [f64; D],
    pub stored_after_j_m3: f64,
    pub stored_change_j_m3: f64,
    pub heat_j_m3: f64,
}

pub(crate) trait MemoryLaw<const D: usize>: Clone {
    type State: Clone;
    type Preparation;

    fn polarization(state: &Self::State) -> &[f64; D];
    fn high_frequency_permittivity_f_m(&self) -> f64;
    fn displacement_c_m2(&self, field: [f64; D], states: &[Self::State]) -> Result<[f64; D]>;
    fn stored_energy_density_j_m3(&self, field: [f64; D], states: &[Self::State]) -> Result<f64>;
    fn prepare_memory(&self, dt_s: f64) -> Result<Self::Preparation>;
    fn midpoint_permittivity_f_m(&self, prepared: &Self::Preparation) -> Result<f64>;
    fn midpoint_history_c_m2(
        &self,
        prepared: &Self::Preparation,
        before_field: [f64; D],
        states: &[Self::State],
    ) -> Result<[f64; D]>;
    fn recover_midpoint(
        &self,
        prepared: &Self::Preparation,
        before_field: [f64; D],
        after_field: [f64; D],
        states: &[Self::State],
    ) -> Result<Recovery<Self::State, D>>;
}

impl<const D: usize> MemoryLaw<D> for ResolvedDebyeMaterial {
    type State = [f64; D];
    type Preparation = f64;

    fn polarization(state: &Self::State) -> &[f64; D] {
        state
    }
    fn high_frequency_permittivity_f_m(&self) -> f64 {
        ResolvedDebyeMaterial::high_frequency_permittivity_f_m(self)
    }
    fn displacement_c_m2(&self, field: [f64; D], states: &[Self::State]) -> Result<[f64; D]> {
        ResolvedDebyeMaterial::displacement_c_m2(self, field, states)
    }
    fn stored_energy_density_j_m3(&self, field: [f64; D], states: &[Self::State]) -> Result<f64> {
        ResolvedDebyeMaterial::stored_energy_density_j_m3(self, field, states)
    }
    fn prepare_memory(&self, dt_s: f64) -> Result<Self::Preparation> {
        Ok(dt_s)
    }
    fn midpoint_permittivity_f_m(&self, prepared: &Self::Preparation) -> Result<f64> {
        ResolvedDebyeMaterial::midpoint_permittivity_f_m(self, *prepared)
    }
    fn midpoint_history_c_m2(
        &self,
        prepared: &Self::Preparation,
        before_field: [f64; D],
        states: &[Self::State],
    ) -> Result<[f64; D]> {
        ResolvedDebyeMaterial::midpoint_history_c_m2(self, *prepared, before_field, states)
    }
    fn recover_midpoint(
        &self,
        prepared: &Self::Preparation,
        before_field: [f64; D],
        after_field: [f64; D],
        states: &[Self::State],
    ) -> Result<Recovery<Self::State, D>> {
        let step = self.step_midpoint(*prepared, before_field, after_field, states)?;
        Ok(Recovery {
            states_after: step.polarization_after_c_m2,
            displacement_after_c_m2: step.displacement_after_c_m2,
            stored_after_j_m3: step.stored_after_j_m3,
            stored_change_j_m3: step.stored_change_j_m3,
            heat_j_m3: step.heat_j_m3,
        })
    }
}

impl<const D: usize> MemoryLaw<D> for ResolvedLorentzMaterial {
    type State = LorentzState<D>;
    type Preparation = PreparedLorentzMaterialMidpoint;

    fn polarization(state: &Self::State) -> &[f64; D] {
        &state.polarization_c_m2
    }
    fn high_frequency_permittivity_f_m(&self) -> f64 {
        ResolvedLorentzMaterial::high_frequency_permittivity_f_m(self)
    }
    fn displacement_c_m2(&self, field: [f64; D], states: &[Self::State]) -> Result<[f64; D]> {
        ResolvedLorentzMaterial::displacement_c_m2(self, field, states)
    }
    fn stored_energy_density_j_m3(&self, field: [f64; D], states: &[Self::State]) -> Result<f64> {
        ResolvedLorentzMaterial::stored_energy_density_j_m3(self, field, states)
    }
    fn prepare_memory(&self, dt_s: f64) -> Result<Self::Preparation> {
        ResolvedLorentzMaterial::prepare_midpoint(self, dt_s)
    }
    fn midpoint_permittivity_f_m(&self, prepared: &Self::Preparation) -> Result<f64> {
        Ok(prepared.midpoint_permittivity_f_m())
    }
    fn midpoint_history_c_m2(
        &self,
        prepared: &Self::Preparation,
        before_field: [f64; D],
        states: &[Self::State],
    ) -> Result<[f64; D]> {
        prepared.midpoint_history_c_m2(before_field, states)
    }
    fn recover_midpoint(
        &self,
        prepared: &Self::Preparation,
        before_field: [f64; D],
        after_field: [f64; D],
        states: &[Self::State],
    ) -> Result<Recovery<Self::State, D>> {
        let step = prepared.step_midpoint(before_field, after_field, states)?;
        Ok(Recovery {
            states_after: step.states_after,
            displacement_after_c_m2: step.displacement_after_c_m2,
            stored_after_j_m3: step.stored_after_j_m3,
            stored_change_j_m3: step.stored_change_j_m3,
            heat_j_m3: step.heat_j_m3,
        })
    }
}
