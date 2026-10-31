//! Exact free-sided axial reductions of source-bound thermo-creep laws.
//!
//! The three-dimensional `ThermalCreep` owner remains authoritative for
//! expansion, entropy, reciprocal heat and the complete work account. This
//! module only eliminates the two free transverse stresses so existing fiber
//! sections and force-based members can consume that same law.
//!
//! For isotropic free strain `f(T)`, the axial mechanical coordinate is
//! `epsilon_x - f(T)`. The existing Norton or joint Norton/J2 uniaxial owner
//! supplies the exact zero-transverse-stress response; its compatible tensor
//! strain is then shifted by `f(T) I`. Energy admission is delegated back to
//! the original tensor thermocreep owner, so no second heat or return map exists.

use super::expansion::{ExpandingCreepSource, ThermalCreep, ThermalCreepState};
use crate::creep::plastic::source::{NortonJ2Source, ResolvedNortonJ2};
use crate::creep::plastic::uniaxial::{UniaxialNortonJ2, UniaxialNortonJ2State};
use crate::creep::source::{CreepSourceError, NortonCreepSource, ResolvedNortonCreep};
use crate::creep::uniaxial::{UniaxialNorton, UniaxialNortonState};
use crate::creep::{CreepError, CreepIncrement, CreepTolerances, NortonCreep};
use crate::fiber::Uniaxial;
use crate::state_point::ResolvedMaterialStatePoint;
use crate::{MaterialAdmissibility, SmallStrainLaw, Voigt};
use fs_evidence::ModelCard;

fn add_free(mut strain: Voigt, free: f64) -> Voigt {
    for value in &mut strain[..3] {
        *value += free;
    }
    strain
}

/// Accepted free-sided Norton history plus its actual source temperature.
#[derive(Clone, Debug, PartialEq)]
pub struct UniaxialThermalNortonState {
    pub mechanical: UniaxialNortonState,
    temperature_k: f64,
}
impl UniaxialThermalNortonState {
    pub const fn temperature_k(&self) -> f64 {
        self.temperature_k
    }
}

/// One uncommitted free-sided Norton thermocreep response.
#[derive(Clone, Debug, PartialEq)]
pub struct UniaxialThermalNortonResponse {
    pub strain: Voigt,
    pub stress_pa: f64,
    pub tangent_pa: f64,
    pub state: UniaxialThermalNortonState,
}

/// Source-bound free-sided Norton thermocreep using the original tensor owner.
#[derive(Clone, Debug)]
pub struct UniaxialThermalNorton {
    mechanical: UniaxialNorton,
    tensor: ThermalCreep<ResolvedNortonCreep>,
}
impl UniaxialThermalNorton {
    pub fn kinetics(&self) -> &NortonCreep {
        self.mechanical.backbone()
    }
    pub fn material(&self) -> Option<&ResolvedMaterialStatePoint> {
        self.mechanical.material()
    }
    pub fn expansion_material(&self) -> &ResolvedMaterialStatePoint {
        self.tensor.expansion_material()
    }
    pub fn tensor_law(&self) -> &ThermalCreep<ResolvedNortonCreep> {
        &self.tensor
    }

    pub fn tensor_strain(
        &self,
        axial_strain: f64,
        state: &UniaxialThermalNortonState,
    ) -> Result<Voigt, CreepError> {
        let free = self.tensor.free_strain_at_temperature(state.temperature_k)?;
        let mechanical = self
            .mechanical
            .tensor_strain(axial_strain - free, &state.mechanical)?;
        Ok(add_free(mechanical, free))
    }

    pub fn response(
        &self,
        axial_strain: f64,
        old: &UniaxialThermalNortonState,
    ) -> Result<UniaxialThermalNortonResponse, CreepError> {
        self.tensor_strain(axial_strain, old)?;
        let free = self.tensor.current_free_strain()?;
        let response = self
            .mechanical
            .response(axial_strain - free, &old.mechanical)?;
        Ok(UniaxialThermalNortonResponse {
            strain: add_free(response.strain, free),
            stress_pa: response.stress_pa,
            tangent_pa: response.tangent_pa,
            state: UniaxialThermalNortonState {
                mechanical: response.state,
                temperature_k: self.tensor.thermoelastic().temperature_k(),
            },
        })
    }

    pub fn advance_with_energy(
        &self,
        old_strain: f64,
        old: &UniaxialThermalNortonState,
        strain: f64,
        tolerances: CreepTolerances,
    ) -> Result<CreepIncrement<UniaxialThermalNortonState>, CreepError> {
        let old_tensor = self.tensor_strain(old_strain, old)?;
        let response = self.response(strain, old)?;
        let tensor_state =
            ThermalCreepState::from_parts(old.mechanical.tensor_history(), old.temperature_k);
        let checked =
            self.tensor
                .advance_with_energy(&old_tensor, &tensor_state, &response.strain, tolerances)?;
        let expected = response.state.mechanical.tensor_history();
        let mut residual = (checked.stress_pa[0] - response.stress_pa).abs();
        for stress in &checked.stress_pa[1..] {
            residual = residual.max(stress.abs());
        }
        let e = self.kinetics().elastic().youngs;
        for (actual, wanted) in checked
            .state
            .mechanical
            .creep_strain
            .iter()
            .zip(expected.creep_strain)
        {
            residual = residual.max(e * (actual - wanted).abs());
        }
        residual = residual.max(
            e * (checked.state.mechanical.equivalent_creep_strain
                - expected.equivalent_creep_strain)
                .abs(),
        );
        if checked.state.temperature_k() != response.state.temperature_k
            || !residual.is_finite()
            || residual > tolerances.stress_pa
        {
            return Err(CreepError::StressResidual {
                residual_pa: residual,
                tolerance_pa: tolerances.stress_pa,
            });
        }
        Ok(CreepIncrement {
            strain: response.strain,
            stress_pa: checked.stress_pa,
            tangent_pa: checked.tangent_pa,
            state: response.state,
            energy: checked.energy,
        })
    }
}
impl Uniaxial for UniaxialThermalNorton {
    type State = UniaxialThermalNortonState;
    fn initial_state(&self) -> Self::State {
        Self::State {
            mechanical: self.mechanical.initial_state(),
            temperature_k: self.tensor.thermoelastic().temperature_k(),
        }
    }
    fn stress(&self, strain: f64, state: &Self::State) -> f64 {
        self.response(strain, state)
            .map(|v| v.stress_pa)
            .unwrap_or(f64::NAN)
    }
    fn tangent(&self, strain: f64, state: &Self::State) -> f64 {
        self.response(strain, state)
            .map(|v| v.tangent_pa)
            .unwrap_or(f64::NAN)
    }
    fn update_state(&self, strain: f64, state: &Self::State) -> Self::State {
        self.response(strain, state)
            .map(|v| v.state)
            .unwrap_or(Self::State {
                mechanical: state.mechanical.clone(),
                temperature_k: f64::NAN,
            })
    }
    fn admissibility(&self) -> MaterialAdmissibility {
        self.tensor.admissibility()
    }
    fn card(&self) -> ModelCard {
        self.tensor.card()
    }
}

/// Accepted free-sided joint Norton/J2 history plus actual source temperature.
#[derive(Clone, Debug, PartialEq)]
pub struct UniaxialThermalNortonJ2State {
    pub mechanical: UniaxialNortonJ2State,
    temperature_k: f64,
}
impl UniaxialThermalNortonJ2State {
    pub const fn temperature_k(&self) -> f64 {
        self.temperature_k
    }
}

/// One uncommitted free-sided joint thermocreep response.
#[derive(Clone, Debug, PartialEq)]
pub struct UniaxialThermalNortonJ2Response {
    pub strain: Voigt,
    pub stress_pa: f64,
    pub tangent_pa: f64,
    pub state: UniaxialThermalNortonJ2State,
}

/// Source-bound free-sided simultaneous creep/yield thermomechanics.
#[derive(Clone, Debug)]
pub struct UniaxialThermalNortonJ2 {
    mechanical: UniaxialNortonJ2,
    tensor: ThermalCreep<ResolvedNortonJ2>,
}
impl UniaxialThermalNortonJ2 {
    pub fn kinetics(&self) -> &NortonCreep {
        self.mechanical.backbone().creep()
    }
    pub fn material(&self) -> Option<&ResolvedMaterialStatePoint> {
        self.mechanical.material()
    }
    pub fn expansion_material(&self) -> &ResolvedMaterialStatePoint {
        self.tensor.expansion_material()
    }
    pub fn tensor_law(&self) -> &ThermalCreep<ResolvedNortonJ2> {
        &self.tensor
    }
    pub fn hardening(&self) -> f64 {
        self.mechanical.backbone().plastic().hardening
    }
    pub fn yield_stress(&self) -> f64 {
        self.mechanical.backbone().plastic().yield_stress
    }

    pub fn tensor_strain(
        &self,
        axial_strain: f64,
        state: &UniaxialThermalNortonJ2State,
    ) -> Result<Voigt, CreepError> {
        let free = self.tensor.free_strain_at_temperature(state.temperature_k)?;
        let mechanical = self
            .mechanical
            .tensor_strain(axial_strain - free, &state.mechanical)?;
        Ok(add_free(mechanical, free))
    }

    pub fn response(
        &self,
        axial_strain: f64,
        old: &UniaxialThermalNortonJ2State,
    ) -> Result<UniaxialThermalNortonJ2Response, CreepError> {
        self.tensor_strain(axial_strain, old)?;
        let free = self.tensor.current_free_strain()?;
        let response = self
            .mechanical
            .response(axial_strain - free, &old.mechanical)?;
        Ok(UniaxialThermalNortonJ2Response {
            strain: add_free(response.strain, free),
            stress_pa: response.stress_pa,
            tangent_pa: response.tangent_pa,
            state: UniaxialThermalNortonJ2State {
                mechanical: response.state,
                temperature_k: self.tensor.thermoelastic().temperature_k(),
            },
        })
    }

    pub fn advance_with_energy(
        &self,
        old_strain: f64,
        old: &UniaxialThermalNortonJ2State,
        strain: f64,
        tolerances: CreepTolerances,
    ) -> Result<CreepIncrement<UniaxialThermalNortonJ2State>, CreepError> {
        let old_tensor = self.tensor_strain(old_strain, old)?;
        let response = self.response(strain, old)?;
        let tensor_state =
            ThermalCreepState::from_parts(old.mechanical.tensor_history(), old.temperature_k);
        let checked =
            self.tensor
                .advance_with_energy(&old_tensor, &tensor_state, &response.strain, tolerances)?;
        let expected = response.state.mechanical.tensor_history();
        let mut residual = (checked.stress_pa[0] - response.stress_pa).abs();
        for stress in &checked.stress_pa[1..] {
            residual = residual.max(stress.abs());
        }
        let e = self.kinetics().elastic().youngs;
        for (actual, wanted) in checked
            .state
            .mechanical
            .creep
            .creep_strain
            .iter()
            .zip(expected.creep.creep_strain)
            .chain(
                checked
                    .state
                    .mechanical
                    .plastic
                    .plastic_strain
                    .iter()
                    .zip(expected.plastic.plastic_strain),
            )
        {
            residual = residual.max(e * (actual - wanted).abs());
        }
        residual = residual.max(
            e * (checked.state.mechanical.creep.equivalent_creep_strain
                - expected.creep.equivalent_creep_strain)
                .abs(),
        );
        residual = residual.max(
            e.max(self.hardening())
                * (checked.state.mechanical.plastic.alpha - expected.plastic.alpha).abs(),
        );
        if checked.state.temperature_k() != response.state.temperature_k
            || !residual.is_finite()
            || residual > tolerances.stress_pa
        {
            return Err(CreepError::StressResidual {
                residual_pa: residual,
                tolerance_pa: tolerances.stress_pa,
            });
        }
        Ok(CreepIncrement {
            strain: response.strain,
            stress_pa: checked.stress_pa,
            tangent_pa: checked.tangent_pa,
            state: response.state,
            energy: checked.energy,
        })
    }
}
impl Uniaxial for UniaxialThermalNortonJ2 {
    type State = UniaxialThermalNortonJ2State;
    fn initial_state(&self) -> Self::State {
        Self::State {
            mechanical: self.mechanical.initial_state(),
            temperature_k: self.tensor.thermoelastic().temperature_k(),
        }
    }
    fn stress(&self, strain: f64, state: &Self::State) -> f64 {
        self.response(strain, state)
            .map(|v| v.stress_pa)
            .unwrap_or(f64::NAN)
    }
    fn tangent(&self, strain: f64, state: &Self::State) -> f64 {
        self.response(strain, state)
            .map(|v| v.tangent_pa)
            .unwrap_or(f64::NAN)
    }
    fn update_state(&self, strain: f64, state: &Self::State) -> Self::State {
        self.response(strain, state)
            .map(|v| v.state)
            .unwrap_or(Self::State {
                mechanical: state.mechanical.clone(),
                temperature_k: f64::NAN,
            })
    }
    fn admissibility(&self) -> MaterialAdmissibility {
        self.tensor.admissibility()
    }
    fn card(&self) -> ModelCard {
        self.tensor.card()
    }
}

impl<'a> ExpandingCreepSource<NortonCreepSource<'a>> {
    /// Resolve one free-sided fiber law at the actual endpoint temperature.
    pub fn at_uniaxial_temperature(
        &self,
        temperature_k: f64,
        duration_s: f64,
    ) -> Result<UniaxialThermalNorton, CreepSourceError> {
        let mechanical = UniaxialNorton::from_source(
            self.mechanical_source()
                .at_temperature(temperature_k, duration_s)?,
        );
        Ok(UniaxialThermalNorton {
            mechanical,
            tensor: self.at_temperature(temperature_k, duration_s)?,
        })
    }
}

impl<'a> ExpandingCreepSource<NortonJ2Source<'a>> {
    /// Resolve one free-sided simultaneous creep/yield law at endpoint temperature.
    pub fn at_uniaxial_temperature(
        &self,
        temperature_k: f64,
        duration_s: f64,
    ) -> Result<UniaxialThermalNortonJ2, CreepSourceError> {
        let mechanical = UniaxialNortonJ2::from_source(
            self.mechanical_source()
                .at_temperature(temperature_k, duration_s)?,
        );
        Ok(UniaxialThermalNortonJ2 {
            mechanical,
            tensor: self.at_temperature(temperature_k, duration_s)?,
        })
    }
}
