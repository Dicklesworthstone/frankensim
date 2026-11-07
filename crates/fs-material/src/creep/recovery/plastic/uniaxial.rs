//! Exact free-sided recovery, permanent creep and isotropic plasticity.
//!
//! The original simultaneous axial owner solves
//! `q + E (dc+dp) = |E (epsilon-epsilon_c_old-epsilon_p_old)|` using the private
//! recovery-transformed primary scale and independent old resistance memory.
//! Both signed permanent strains, both accumulated paths and the solved creep
//! increment retain their original owners. The native recovering tensor law
//! audits elastic/hardening storage, separate physical heats and numerical loss.
//! No additional axial constitutive root or resistance energy is introduced.

use super::source::{RecoveringNortonJ2Source, ResolvedRecoveringNortonJ2};
use super::{RecoveringNortonJ2, RecoveringNortonJ2State};
use crate::creep::plastic::uniaxial::{UniaxialNortonJ2, UniaxialNortonJ2State};
use crate::creep::recovery::validate_resistance;
use crate::creep::source::CreepSourceError;
use crate::creep::{CreepError, CreepIncrement, CreepTolerances, NortonCreep, Result};
use crate::fiber::Uniaxial;
use crate::state_point::ResolvedMaterialStatePoint;
use crate::{MaterialAdmissibility, SmallStrainLaw, Voigt};
use fs_evidence::ModelCard;
use std::sync::Arc;

/// Both native permanent axial histories and independent primary resistance.
#[derive(Clone, Debug, PartialEq)]
pub struct UniaxialRecoveringNortonJ2State {
    native: UniaxialNortonJ2State,
    resistance_strain: f64,
}

impl UniaxialRecoveringNortonJ2State {
    pub const fn native(&self) -> &UniaxialNortonJ2State {
        &self.native
    }

    pub const fn creep_strain(&self) -> f64 {
        self.native.creep_strain()
    }

    pub const fn plastic_strain(&self) -> f64 {
        self.native.plastic_strain()
    }

    pub const fn equivalent_creep_strain(&self) -> f64 {
        self.native.equivalent_creep_strain()
    }

    pub const fn equivalent_plastic_strain(&self) -> f64 {
        self.native.equivalent_plastic_strain()
    }

    /// Kinetic resistance is bounded by permanent creep accumulation alone.
    pub const fn resistance_strain(&self) -> f64 {
        self.resistance_strain
    }

    /// Original tensor history with both incompressible lateral flow tensors.
    pub fn tensor_history(&self) -> RecoveringNortonJ2State {
        RecoveringNortonJ2State {
            native: self.native.tensor_history(),
            resistance_strain: self.resistance_strain,
        }
    }
}

/// One uncommitted axial candidate and its compatible complete tensor strain.
#[derive(Clone, Debug, PartialEq)]
pub struct UniaxialRecoveringNortonJ2Response {
    pub strain: Voigt,
    pub stress_pa: f64,
    pub tangent_pa: f64,
    pub state: UniaxialRecoveringNortonJ2State,
}

/// Original physical interval and complete optional source receipts, with one
/// privately prepared exact free-sided joint return.
#[derive(Clone, Debug, PartialEq)]
pub struct UniaxialRecoveringNortonJ2 {
    backbone: RecoveringNortonJ2,
    prepared: UniaxialNortonJ2,
    material: Option<Arc<ResolvedMaterialStatePoint>>,
}

impl UniaxialRecoveringNortonJ2 {
    pub fn new(backbone: RecoveringNortonJ2) -> Self {
        let prepared = UniaxialNortonJ2::new(backbone.prepared_norton_j2().clone());
        Self {
            backbone,
            prepared,
            material: None,
        }
    }

    /// Retain every original creep, recovery, yield and hardening claim.
    pub fn from_source(source: ResolvedRecoveringNortonJ2) -> Self {
        let mut law = Self::new(source.law().clone());
        law.material = Some(Arc::new(source.material().clone()));
        law
    }

    pub const fn backbone(&self) -> &RecoveringNortonJ2 {
        &self.backbone
    }

    /// Original Norton inputs; no timestep-scaled primary coefficient escapes.
    pub const fn kinetics(&self) -> &NortonCreep {
        self.backbone.norton_j2().creep()
    }

    pub fn material(&self) -> Option<&ResolvedMaterialStatePoint> {
        self.material.as_deref()
    }

    /// Reconstruct accepted free-sided strain without advancing either flow.
    pub fn tensor_strain(
        &self,
        axial_strain: f64,
        state: &UniaxialRecoveringNortonJ2State,
    ) -> Result<Voigt> {
        let strain = self.prepared.tensor_strain(axial_strain, &state.native)?;
        self.backbone
            .validate_history(&strain, &state.tensor_history())?;
        Ok(strain)
    }

    /// One original coupled axial root from frozen permanent/resistance memory.
    pub fn response(
        &self,
        strain: f64,
        old: &UniaxialRecoveringNortonJ2State,
    ) -> Result<UniaxialRecoveringNortonJ2Response> {
        validate_resistance(old.resistance_strain, old.equivalent_creep_strain())?;
        let (response, delta) =
            self.prepared
                .response_with_memory(strain, &old.native, old.resistance_strain)?;
        let resistance_strain = self.backbone.next_resistance(
            old.resistance_strain,
            delta,
            response.state.equivalent_creep_strain(),
        )?;
        Ok(UniaxialRecoveringNortonJ2Response {
            strain: response.strain,
            stress_pa: response.stress_pa,
            tangent_pa: response.tangent_pa,
            state: UniaxialRecoveringNortonJ2State {
                native: response.state,
                resistance_strain,
            },
        })
    }

    /// Audit the accepted candidate through the original joint tensor owner.
    /// Returned tangent is the full tensor tangent; [`Self::response`] supplies
    /// its free-sided axial restriction. Only original physical heat is retained.
    pub fn advance_with_energy(
        &self,
        old_strain: f64,
        old: &UniaxialRecoveringNortonJ2State,
        strain: f64,
        tolerances: CreepTolerances,
    ) -> Result<CreepIncrement<UniaxialRecoveringNortonJ2State>> {
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
        for (actual, wanted) in increment
            .state
            .native
            .creep
            .creep_strain
            .iter()
            .zip(expected.native.creep.creep_strain)
            .chain(
                increment
                    .state
                    .native
                    .plastic
                    .plastic_strain
                    .iter()
                    .zip(expected.native.plastic.plastic_strain),
            )
        {
            residual = residual.max(e * (actual - wanted).abs());
        }
        residual = residual
            .max(
                e * (increment.state.native.creep.equivalent_creep_strain
                    - expected.native.creep.equivalent_creep_strain)
                    .abs(),
            )
            .max(
                e.max(self.backbone.norton_j2().plastic().hardening)
                    * (increment.state.native.plastic.alpha - expected.native.plastic.alpha).abs(),
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

impl Uniaxial for UniaxialRecoveringNortonJ2 {
    type State = UniaxialRecoveringNortonJ2State;

    fn initial_state(&self) -> Self::State {
        UniaxialRecoveringNortonJ2State {
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
            .unwrap_or_else(|_| UniaxialRecoveringNortonJ2State {
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

impl RecoveringNortonJ2Source<'_> {
    /// Resolve one source-bounded axial interval without advancing specimen history.
    pub fn at_uniaxial_temperature(
        &self,
        temperature_k: f64,
        duration_s: f64,
    ) -> core::result::Result<UniaxialRecoveringNortonJ2, CreepSourceError> {
        self.at_temperature(temperature_k, duration_s)
            .map(UniaxialRecoveringNortonJ2::from_source)
    }
}
