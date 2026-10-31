//! Free-sided axial restriction of the native Norton/Arrhenius creep law.
//!
//! An axial fiber has uniaxial stress, including elastic lateral contraction
//! and traceless creep flow. Eliminating the transverse stress gives
//! `q + E dt rate(q) = |E (epsilon - epsilon_c_old)|`. This module uses the
//! same bounded logarithmic kinetic solve as the 3-D owner, with E replacing
//! 3G. The checked work account is evaluated by that original tensor owner.
//! Every trial represents one declared physical interval from frozen history;
//! repeated stress/tangent queries do not commit time or history.

use super::source::ResolvedNortonCreep;
use super::{
    CreepError, CreepIncrement, CreepState, CreepTolerances, NortonCreep, Result, require,
};
use crate::fiber::Uniaxial;
use crate::state_point::ResolvedMaterialStatePoint;
use crate::{MaterialAdmissibility, SmallStrainLaw, Voigt};
use fs_evidence::ModelCard;
use fs_math::det;
use std::sync::Arc;

/// Exact free-sided restriction, optionally retaining the original source
/// bundle. No independent axial creep coefficient or heat fraction is fitted.
#[derive(Clone, Debug, PartialEq)]
pub struct UniaxialNorton {
    backbone: NortonCreep,
    material: Option<Arc<ResolvedMaterialStatePoint>>,
    log_step_scale: f64,
}

/// Signed axial creep and accumulated path length on the uniaxial subspace.
/// General tensor histories belong to the original three-dimensional owner.
#[derive(Clone, Debug, PartialEq)]
pub struct UniaxialNortonState {
    creep_strain: f64,
    equivalent_creep_strain: f64,
}
impl UniaxialNortonState {
    pub const fn creep_strain(&self) -> f64 {
        self.creep_strain
    }
    pub const fn equivalent_creep_strain(&self) -> f64 {
        self.equivalent_creep_strain
    }
    /// Original tensor-Voigt history with incompressible lateral creep strain.
    pub fn tensor_history(&self) -> CreepState {
        let c = self.creep_strain;
        CreepState {
            creep_strain: [c, -0.5 * c, -0.5 * c, 0.0, 0.0, 0.0],
            equivalent_creep_strain: self.equivalent_creep_strain,
        }
    }
}

/// Candidate response with the complete compatible total strain. Its axial
/// tangent is condensed at zero transverse stress, not the 3-D C_xxxx entry.
#[derive(Clone, Debug, PartialEq)]
pub struct UniaxialNortonResponse {
    pub strain: Voigt,
    pub stress_pa: f64,
    pub tangent_pa: f64,
    pub state: UniaxialNortonState,
}

impl UniaxialNorton {
    /// Restrict an already admitted native physical-step law.
    pub fn new(backbone: NortonCreep) -> Self {
        let log_step_scale = det::ln(backbone.elastic.youngs)
            + backbone.log_duration
            + backbone.log_rate_at_reference_stress
            - backbone.log_reference_stress;
        Self {
            backbone,
            material: None,
            log_step_scale,
        }
    }

    /// Retain all original material identities, query conventions and receipts.
    /// Temperature/duration changes should resolve a new native source step and
    /// preserve the accepted axial state, rather than recreating virgin history.
    pub fn from_source(source: ResolvedNortonCreep) -> Self {
        let mut law = Self::new(source.law().clone());
        law.material = Some(Arc::new(source.material().clone()));
        law
    }

    pub const fn backbone(&self) -> &NortonCreep {
        &self.backbone
    }
    pub fn material(&self) -> Option<&ResolvedMaterialStatePoint> {
        self.material.as_deref()
    }

    /// Compatible accepted tensor strain without advancing the creep law.
    pub fn tensor_strain(&self, axial_strain: f64, state: &UniaxialNortonState) -> Result<Voigt> {
        let transverse = -self.backbone.elastic.poisson * (axial_strain - state.creep_strain)
            - 0.5 * state.creep_strain;
        let strain = [axial_strain, transverse, transverse, 0.0, 0.0, 0.0];
        self.backbone
            .validate_history(&strain, &state.tensor_history())?;
        Ok(strain)
    }

    /// One implicit axial candidate from accepted history. Full tensor strain,
    /// elastic strain and creep strain retain the original admission ceiling.
    pub fn response(
        &self,
        strain: f64,
        old: &UniaxialNortonState,
    ) -> Result<UniaxialNortonResponse> {
        let old_tensor = old.tensor_history();
        self.backbone
            .validate_history(&old_tensor.creep_strain, &old_tensor)?;
        let trial = self.backbone.elastic.youngs * (strain - old.creep_strain);
        require(trial.is_finite(), "axial trial stress is not representable")?;
        let relaxed = self.backbone.relaxation(trial.abs(), self.log_step_scale)?;
        let state = UniaxialNortonState {
            creep_strain: old.creep_strain + trial.signum() * relaxed.delta_alpha,
            equivalent_creep_strain: old.equivalent_creep_strain + relaxed.delta_alpha,
        };
        let tensor = self.tensor_strain(strain, &state)?;
        self.backbone.validate_history(&tensor, &old_tensor)?;
        let stress_pa = trial.signum() * relaxed.q;
        let tangent_pa = self.backbone.elastic.youngs * relaxed.radial;
        require(
            stress_pa.is_finite() && tangent_pa.is_finite(),
            "axial stress or tangent is not representable",
        )?;
        Ok(UniaxialNortonResponse {
            strain: tensor,
            stress_pa,
            tangent_pa,
            state,
        })
    }

    /// Audit the candidate through the original full-tensor work/heat owner.
    /// Returned tangent remains the full 3-D tangent; `response` supplies the
    /// condensed axial tangent used by fiber equilibrium. All energy densities
    /// are per original material volume, including physical and numerical loss.
    pub fn advance_with_energy(
        &self,
        old_strain: f64,
        old_state: &UniaxialNortonState,
        strain: f64,
        tolerances: CreepTolerances,
    ) -> Result<CreepIncrement> {
        let old_tensor = self.tensor_strain(old_strain, old_state)?;
        let response = self.response(strain, old_state)?;
        let increment = self.backbone.advance_with_energy(
            &old_tensor,
            &old_state.tensor_history(),
            &response.strain,
            tolerances,
        )?;
        let mut residual = (increment.stress_pa[0] - response.stress_pa).abs();
        for stress in &increment.stress_pa[1..] {
            residual = residual.max(stress.abs());
        }
        let e = self.backbone.elastic.youngs;
        let native_history = response.state.tensor_history();
        for (actual, expected) in increment
            .state
            .creep_strain
            .iter()
            .zip(native_history.creep_strain)
        {
            residual = residual.max(e * (actual - expected).abs());
        }
        residual = residual.max(
            e * (increment.state.equivalent_creep_strain - response.state.equivalent_creep_strain)
                .abs(),
        );
        if !residual.is_finite() || residual > tolerances.stress_pa {
            return Err(CreepError::StressResidual {
                residual_pa: residual,
                tolerance_pa: tolerances.stress_pa,
            });
        }
        Ok(increment)
    }
}

impl Uniaxial for UniaxialNorton {
    type State = UniaxialNortonState;
    fn initial_state(&self) -> Self::State {
        UniaxialNortonState {
            creep_strain: 0.0,
            equivalent_creep_strain: 0.0,
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
            .unwrap_or(UniaxialNortonState {
                creep_strain: f64::NAN,
                equivalent_creep_strain: f64::NAN,
            })
    }
    fn admissibility(&self) -> MaterialAdmissibility {
        self.backbone.admissibility()
    }
    fn card(&self) -> ModelCard {
        self.backbone.card()
    }
}
