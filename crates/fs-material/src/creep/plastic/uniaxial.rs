//! Exact free-sided fibers with simultaneous permanent creep and plasticity.
//!
//! Eliminating transverse stress gives
//! `q + E (delta_c + delta_p) = |E (epsilon - epsilon_c_old - epsilon_p_old)|`.
//! The original combined scalar update supplies both increments and its exact
//! derivative, using E in place of 3G. Each signed axial history maps to its
//! own traceless tensor `[e, -e/2, -e/2, 0, 0, 0]`; accumulated paths remain
//! distinct through unloading and reversal. The native full-tensor owner audits
//! stored elastic/hardening energy, both physical heats and numerical loss.
//!
//! All trial calls borrow accepted history without advancing time. This is the
//! same bounded small-strain isotropic model, not a finite-strain beam law.

use super::source::ResolvedNortonJ2;
use super::{CoupledStiffness, NortonJ2, NortonJ2State};
use crate::creep::{CreepError, CreepIncrement, CreepState, CreepTolerances, Result, require};
use crate::fiber::Uniaxial;
use crate::plastic::J2State;
use crate::state_point::ResolvedMaterialStatePoint;
use crate::{MaterialAdmissibility, SmallStrainLaw, Voigt};
use fs_evidence::ModelCard;
use std::sync::Arc;

#[derive(Clone, Debug, PartialEq)]
enum Backbone {
    Authored(NortonJ2),
    Sourced(Arc<ResolvedNortonJ2>),
}

/// Native free-sided reduction retaining the complete optional source bundle.
/// Cloning a source-bound fiber shares its immutable law and original receipts.
#[derive(Clone, Debug, PartialEq)]
pub struct UniaxialNortonJ2 {
    backbone: Backbone,
    stiffness: CoupledStiffness,
}

/// Separate signed permanent strains and accumulated equivalent paths.
#[derive(Clone, Debug, PartialEq)]
pub struct UniaxialNortonJ2State {
    creep_strain: f64,
    plastic_strain: f64,
    equivalent_creep_strain: f64,
    equivalent_plastic_strain: f64,
}

impl UniaxialNortonJ2State {
    pub const fn creep_strain(&self) -> f64 {
        self.creep_strain
    }
    pub const fn plastic_strain(&self) -> f64 {
        self.plastic_strain
    }
    pub const fn equivalent_creep_strain(&self) -> f64 {
        self.equivalent_creep_strain
    }
    pub const fn equivalent_plastic_strain(&self) -> f64 {
        self.equivalent_plastic_strain
    }

    /// The two original native histories, with their incompressible lateral flow.
    pub fn tensor_history(&self) -> NortonJ2State {
        let c = self.creep_strain;
        let p = self.plastic_strain;
        NortonJ2State {
            creep: CreepState {
                creep_strain: [c, -0.5 * c, -0.5 * c, 0.0, 0.0, 0.0],
                equivalent_creep_strain: self.equivalent_creep_strain,
            },
            plastic: J2State {
                plastic_strain: [p, -0.5 * p, -0.5 * p, 0.0, 0.0, 0.0],
                alpha: self.equivalent_plastic_strain,
            },
        }
    }
}

/// Candidate complete tensor strain and condensed uniaxial response.
#[derive(Clone, Debug, PartialEq)]
pub struct UniaxialNortonJ2Response {
    pub strain: Voigt,
    pub stress_pa: f64,
    pub tangent_pa: f64,
    pub state: UniaxialNortonJ2State,
}

impl UniaxialNortonJ2 {
    /// Restrict an admitted native law; no new constitutive coefficients are fitted.
    pub fn new(backbone: NortonJ2) -> Self {
        let stiffness = Self::axial_stiffness(&backbone);
        Self {
            backbone: Backbone::Authored(backbone),
            stiffness,
        }
    }

    /// Keep the original union of creep, elasticity and hardening source claims.
    /// A new temperature/duration resolves a new law while keeping accepted state.
    pub fn from_source(source: ResolvedNortonJ2) -> Self {
        let stiffness = Self::axial_stiffness(source.law());
        Self {
            backbone: Backbone::Sourced(Arc::new(source)),
            stiffness,
        }
    }

    fn axial_stiffness(backbone: &NortonJ2) -> CoupledStiffness {
        CoupledStiffness::new(
            backbone.creep(),
            backbone.plastic().hardening,
            backbone.creep().elastic().youngs,
        )
    }

    pub fn backbone(&self) -> &NortonJ2 {
        match &self.backbone {
            Backbone::Authored(law) => law,
            Backbone::Sourced(source) => source.law(),
        }
    }

    pub fn source(&self) -> Option<&ResolvedNortonJ2> {
        match &self.backbone {
            Backbone::Authored(_) => None,
            Backbone::Sourced(source) => Some(source.as_ref()),
        }
    }

    pub fn material(&self) -> Option<&ResolvedMaterialStatePoint> {
        self.source().map(ResolvedNortonJ2::material)
    }

    /// Reconstruct accepted free-sided strain without performing another update.
    pub fn tensor_strain(&self, axial_strain: f64, state: &UniaxialNortonJ2State) -> Result<Voigt> {
        let elastic = axial_strain - state.creep_strain - state.plastic_strain;
        let transverse = -self.backbone().creep().elastic().poisson * elastic
            - 0.5 * state.creep_strain
            - 0.5 * state.plastic_strain;
        let strain = [axial_strain, transverse, transverse, 0.0, 0.0, 0.0];
        self.backbone()
            .validate_history(&strain, &state.tensor_history())?;
        Ok(strain)
    }

    /// One simultaneous candidate from frozen history, including Poisson
    /// contraction and both lateral inelastic strains. Tangent is E*dq/dq_trial.
    pub fn response(
        &self,
        strain: f64,
        old: &UniaxialNortonJ2State,
    ) -> Result<UniaxialNortonJ2Response> {
        let backbone = self.backbone();
        let old_tensor = old.tensor_history();
        backbone.validate_state(&old_tensor)?;
        let e = backbone.creep().elastic().youngs;
        let trial = e * (strain - old.creep_strain - old.plastic_strain);
        require(trial.is_finite(), "axial trial stress is not representable")?;
        let relaxed =
            backbone.relaxation(trial.abs(), old.equivalent_plastic_strain, &self.stiffness)?;
        let state = UniaxialNortonJ2State {
            creep_strain: old.creep_strain + trial.signum() * relaxed.delta_c,
            plastic_strain: old.plastic_strain + trial.signum() * relaxed.delta_p,
            equivalent_creep_strain: old.equivalent_creep_strain + relaxed.delta_c,
            equivalent_plastic_strain: old.equivalent_plastic_strain + relaxed.delta_p,
        };
        let tensor = self.tensor_strain(strain, &state)?;
        backbone.validate_history(&tensor, &old_tensor)?;
        let stress_pa = trial.signum() * relaxed.q;
        let tangent_pa = e * relaxed.radial;
        require(
            stress_pa.is_finite() && tangent_pa.is_finite(),
            "axial stress or tangent is not representable",
        )?;
        Ok(UniaxialNortonJ2Response {
            strain: tensor,
            stress_pa,
            tangent_pa,
            state,
        })
    }

    /// Audit one candidate through the original full-tensor energy owner.
    /// Returned tangent is the 3-D tangent; `response` supplies its free-sided
    /// axial condensation. No accepted strain, history or source is mutated.
    pub fn advance_with_energy(
        &self,
        old_strain: f64,
        old_state: &UniaxialNortonJ2State,
        strain: f64,
        tolerances: CreepTolerances,
    ) -> Result<CreepIncrement<NortonJ2State>> {
        let old_tensor = self.tensor_strain(old_strain, old_state)?;
        let response = self.response(strain, old_state)?;
        let backbone = self.backbone();
        let increment = backbone.advance_with_energy(
            &old_tensor,
            &old_state.tensor_history(),
            &response.strain,
            tolerances,
        )?;
        let mut residual = (increment.stress_pa[0] - response.stress_pa).abs();
        for stress in &increment.stress_pa[1..] {
            residual = residual.max(stress.abs());
        }
        let e = backbone.creep().elastic().youngs;
        let expected = response.state.tensor_history();
        for (actual, wanted) in increment
            .state
            .creep
            .creep_strain
            .iter()
            .zip(expected.creep.creep_strain)
            .chain(
                increment
                    .state
                    .plastic
                    .plastic_strain
                    .iter()
                    .zip(expected.plastic.plastic_strain),
            )
        {
            residual = residual.max(e * (actual - wanted).abs());
        }
        residual = residual.max(
            e * (increment.state.creep.equivalent_creep_strain
                - expected.creep.equivalent_creep_strain)
                .abs(),
        );
        residual = residual.max(
            e.max(backbone.plastic().hardening)
                * (increment.state.plastic.alpha - expected.plastic.alpha).abs(),
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

impl Uniaxial for UniaxialNortonJ2 {
    type State = UniaxialNortonJ2State;
    fn initial_state(&self) -> Self::State {
        UniaxialNortonJ2State {
            creep_strain: 0.0,
            plastic_strain: 0.0,
            equivalent_creep_strain: 0.0,
            equivalent_plastic_strain: 0.0,
        }
    }
    fn stress(&self, strain: f64, old: &Self::State) -> f64 {
        self.response(strain, old)
            .map(|v| v.stress_pa)
            .unwrap_or(f64::NAN)
    }
    fn tangent(&self, strain: f64, old: &Self::State) -> f64 {
        self.response(strain, old)
            .map(|v| v.tangent_pa)
            .unwrap_or(f64::NAN)
    }
    fn update_state(&self, strain: f64, old: &Self::State) -> Self::State {
        self.response(strain, old)
            .map(|v| v.state)
            .unwrap_or(UniaxialNortonJ2State {
                creep_strain: f64::NAN,
                plastic_strain: f64::NAN,
                equivalent_creep_strain: f64::NAN,
                equivalent_plastic_strain: f64::NAN,
            })
    }
    fn admissibility(&self) -> MaterialAdmissibility {
        self.backbone().admissibility()
    }
    fn card(&self) -> ModelCard {
        self.backbone().card()
    }
}
