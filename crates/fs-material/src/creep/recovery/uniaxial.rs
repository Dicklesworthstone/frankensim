//! Exact free-sided axial restriction with independently recovering resistance.
//!
//! The original axial Norton return solves permanent flow once, with Young's
//! modulus and the private recovery-transformed primary scale. Its kinetic
//! memory is resistance, while signed creep and accumulated path remain the
//! original irreversible histories. The native tensor recovery owner audits
//! stress, heat and work; no separate axial heat law or recovery energy is added.

use super::source::{RecoveringNortonSource, ResolvedRecoveringNorton};
use super::{RecoveringCreepState, RecoveringNortonCreep};
use crate::creep::source::CreepSourceError;
use crate::creep::uniaxial::{UniaxialNorton, UniaxialNortonState};
use crate::creep::{CreepError, CreepIncrement, CreepTolerances, NortonCreep, Result};
use crate::fiber::Uniaxial;
use crate::state_point::ResolvedMaterialStatePoint;
use crate::{MaterialAdmissibility, SmallStrainLaw, Voigt};
use fs_evidence::ModelCard;
use std::sync::Arc;

/// Permanent free-sided creep history and its independently recovering memory.
#[derive(Clone, Debug, PartialEq)]
pub struct UniaxialRecoveringNortonState {
    native: UniaxialNortonState,
    resistance_strain: f64,
}

impl UniaxialRecoveringNortonState {
    pub const fn native(&self) -> &UniaxialNortonState {
        &self.native
    }

    pub const fn creep_strain(&self) -> f64 {
        self.native.creep_strain()
    }

    pub const fn equivalent_creep_strain(&self) -> f64 {
        self.native.equivalent_creep_strain()
    }

    /// Nonnegative kinetic resistance; recovery does not undo permanent strain.
    pub const fn resistance_strain(&self) -> f64 {
        self.resistance_strain
    }

    /// Original tensor history, including incompressible lateral creep strain.
    pub fn tensor_history(&self) -> RecoveringCreepState {
        RecoveringCreepState {
            native: self.native.tensor_history(),
            resistance_strain: self.resistance_strain,
        }
    }
}

/// One uncommitted axial candidate with its complete compatible total strain.
#[derive(Clone, Debug, PartialEq)]
pub struct UniaxialRecoveringNortonResponse {
    pub strain: Voigt,
    pub stress_pa: f64,
    pub tangent_pa: f64,
    pub state: UniaxialRecoveringNortonState,
}

/// Original physical interval and receipts, with a private prepared axial return.
#[derive(Clone, Debug, PartialEq)]
pub struct UniaxialRecoveringNorton {
    backbone: RecoveringNortonCreep,
    prepared: UniaxialNorton,
    material: Option<Arc<ResolvedMaterialStatePoint>>,
}

impl UniaxialRecoveringNorton {
    pub fn new(backbone: RecoveringNortonCreep) -> Self {
        let prepared = UniaxialNorton::new(backbone.prepared_norton().clone());
        Self {
            backbone,
            prepared,
            material: None,
        }
    }

    /// Retain all original primary/recovery claims and query coordinates.
    pub fn from_source(source: ResolvedRecoveringNorton) -> Self {
        let mut law = Self::new(source.law().clone());
        law.material = Some(Arc::new(source.material().clone()));
        law
    }

    pub const fn backbone(&self) -> &RecoveringNortonCreep {
        &self.backbone
    }

    /// Original Norton coefficients; the timestep-transformed scale stays private.
    pub const fn kinetics(&self) -> &NortonCreep {
        self.backbone.norton()
    }

    pub fn material(&self) -> Option<&ResolvedMaterialStatePoint> {
        self.material.as_deref()
    }

    /// Accepted compatible strain without consuming another physical interval.
    pub fn tensor_strain(
        &self,
        axial_strain: f64,
        state: &UniaxialRecoveringNortonState,
    ) -> Result<Voigt> {
        let strain = self.prepared.tensor_strain(axial_strain, &state.native)?;
        self.backbone
            .validate_history(&strain, &state.tensor_history())?;
        Ok(strain)
    }

    /// A single native axial root from frozen permanent and resistance histories.
    pub fn response(
        &self,
        strain: f64,
        old: &UniaxialRecoveringNortonState,
    ) -> Result<UniaxialRecoveringNortonResponse> {
        let old_tensor = old.tensor_history();
        self.backbone
            .validate_history(&old_tensor.native.creep_strain, &old_tensor)?;
        let (response, delta) =
            self.prepared
                .response_with_memory(strain, &old.native, old.resistance_strain)?;
        let resistance_strain = self.backbone.next_resistance(
            old.resistance_strain,
            delta,
            response.state.equivalent_creep_strain(),
        )?;
        Ok(UniaxialRecoveringNortonResponse {
            strain: response.strain,
            stress_pa: response.stress_pa,
            tangent_pa: response.tangent_pa,
            state: UniaxialRecoveringNortonState {
                native: response.state,
                resistance_strain,
            },
        })
    }

    /// Audit through the original tensor recovery owner. The returned tangent
    /// is the full tensor tangent; [`Self::response`] returns its free-sided
    /// axial restriction. Resistance recovery creates neither heat nor storage.
    pub fn advance_with_energy(
        &self,
        old_strain: f64,
        old: &UniaxialRecoveringNortonState,
        strain: f64,
        tolerances: CreepTolerances,
    ) -> Result<CreepIncrement<UniaxialRecoveringNortonState>> {
        let old_tensor = self.tensor_strain(old_strain, old)?;
        let response = self.response(strain, old)?;
        let increment = self.backbone.advance_with_energy(
            &old_tensor,
            &old.tensor_history(),
            &response.strain,
            tolerances,
        )?;
        let expected = response.state.tensor_history();
        let mut residual = (increment.stress_pa[0] - response.stress_pa).abs();
        for stress in &increment.stress_pa[1..] {
            residual = residual.max(stress.abs());
        }
        let e = self.kinetics().elastic().youngs;
        for (actual, expected) in increment
            .state
            .native
            .creep_strain
            .iter()
            .zip(expected.native.creep_strain)
        {
            residual = residual.max(e * (actual - expected).abs());
        }
        residual = residual
            .max(
                e * (increment.state.native.equivalent_creep_strain
                    - expected.native.equivalent_creep_strain)
                    .abs(),
            )
            .max(e * (increment.state.resistance_strain - expected.resistance_strain).abs());
        if !residual.is_finite() || residual > tolerances.stress_pa {
            return Err(CreepError::StressResidual {
                residual_pa: residual,
                tolerance_pa: tolerances.stress_pa,
            });
        }
        Ok(CreepIncrement {
            strain: response.strain,
            stress_pa: increment.stress_pa,
            tangent_pa: increment.tangent_pa,
            state: response.state,
            energy: increment.energy,
        })
    }
}

impl Uniaxial for UniaxialRecoveringNorton {
    type State = UniaxialRecoveringNortonState;

    fn initial_state(&self) -> Self::State {
        UniaxialRecoveringNortonState {
            native: self.prepared.initial_state(),
            resistance_strain: 0.0,
        }
    }

    fn stress(&self, strain: f64, state: &Self::State) -> f64 {
        self.response(strain, state)
            .map(|value| value.stress_pa)
            .unwrap_or(f64::NAN)
    }

    fn tangent(&self, strain: f64, state: &Self::State) -> f64 {
        self.response(strain, state)
            .map(|value| value.tangent_pa)
            .unwrap_or(f64::NAN)
    }

    fn update_state(&self, strain: f64, state: &Self::State) -> Self::State {
        self.response(strain, state)
            .map(|value| value.state)
            .unwrap_or_else(|_| UniaxialRecoveringNortonState {
                native: state.native.clone(),
                resistance_strain: f64::NAN,
            })
    }

    fn admissibility(&self) -> MaterialAdmissibility {
        self.backbone.admissibility()
    }

    fn card(&self) -> ModelCard {
        self.backbone.card()
    }
}

impl RecoveringNortonSource<'_> {
    /// Resolve one source-bounded free-sided axial interval without aging history.
    pub fn at_uniaxial_temperature(
        &self,
        temperature_k: f64,
        duration_s: f64,
    ) -> core::result::Result<UniaxialRecoveringNorton, CreepSourceError> {
        self.at_temperature(temperature_k, duration_s)
            .map(UniaxialRecoveringNorton::from_source)
    }
}
