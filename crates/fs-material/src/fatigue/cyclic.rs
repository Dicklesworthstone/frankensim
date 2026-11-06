//! Fatigue observation of original, source-bound cyclic plasticity histories.
//!
//! J2/Prager, Voce and one-through-eight-branch Chaboche laws share the same
//! physical stress recovery: the original elastic map of total mechanical
//! strain minus accepted plastic strain. Backstresses translate the yield
//! surface; they must never be subtracted from the physical fatigue stress.
//!
//! This sealed interface exposes accepted state, not a second return map. The
//! native mechanical owner remains responsible for advancing each complete
//! history, converging equilibrium and admitting its fixed source conditions.
//! It supplies no thermal, viscous or creep strain decomposition and does not
//! transfer a calibration to a different temperature or hardening model.

use core::fmt::Debug;

use super::critical_plane::TensorFatigueSample;
use crate::plastic::J2Plasticity;
use crate::plastic::chaboche::source::{ResolvedChabocheBranches, ResolvedChabocheJ2};
use crate::plastic::source::ResolvedJ2Plasticity;
use crate::plastic::voce::source::ResolvedVoceJ2;
use crate::state_point::ResolvedMaterialStatePoint;
use crate::{MaterialError, SmallStrainLaw, Voigt};

mod sealed {
    pub trait Sealed {}
}

/// Original fixed-condition metal laws admitted by native fatigue consumers.
///
/// Constitutive parameters, source receipts and all internal variables remain
/// owned by the native law. Consumers may inspect its common elastic/initial-
/// yield backbone for calibration binding, but must advance the complete law;
/// the backbone alone cannot replace Voce saturation or Chaboche recovery.
/// Static dispatch retains each actual state type, including every independent
/// backstress tensor, without erased histories or per-sample source cloning.
pub trait SourceBoundFatigueLaw: SmallStrainLaw + Clone + Debug + sealed::Sealed {
    /// The full immutable source bundle, including all nonlinear branch claims.
    fn fatigue_source(&self) -> &ResolvedMaterialStatePoint;

    /// Original elasticity and initial yield reference, never evolved strength.
    fn fatigue_backbone(&self) -> &J2Plasticity;

    /// Accepted signed tensor plastic strain; backstress is a different state.
    fn accepted_plastic_strain<'a>(&self, state: &'a Self::State) -> &'a Voigt;

    /// Validate native finite/strain/history bounds without a constitutive step.
    /// Chaboche admission includes every branch's deviatoric and reachable
    /// magnitude bounds. Equilibrium and instantaneous yield convergence remain
    /// the native solver's responsibility; no arbitrary yield tolerance is
    /// introduced by this observation interface.
    fn validate_fatigue_history(
        &self,
        total_mechanical_strain: &Voigt,
        state: &Self::State,
    ) -> Result<(), MaterialError>;

    /// Recover physical stress from one already accepted material point.
    ///
    /// Useful for reductions such as shells that retain full accepted strain
    /// and history but store only integrated stress resultants. A solver that
    /// already retains its accepted tensor stress should use that tensor and
    /// call `validate_fatigue_history` directly. Neither operation updates or
    /// collapses history, evaluates a yield-surface effective stress, or calls
    /// the algorithmic `SmallStrainLaw::stress` on an accepted plastic state.
    fn accepted_fatigue_sample(
        &self,
        total_mechanical_strain: &Voigt,
        state: &Self::State,
    ) -> Result<TensorFatigueSample, MaterialError> {
        self.validate_fatigue_history(total_mechanical_strain, state)?;
        let plastic = self.accepted_plastic_strain(state);
        let elastic = core::array::from_fn(|i| total_mechanical_strain[i] - plastic[i]);
        let stress_pa = self.fatigue_backbone().elastic.stress(&elastic, &());
        if stress_pa.iter().any(|value| !value.is_finite()) {
            return Err(MaterialError::State {
                what: "accepted fatigue stress is unrepresentable".into(),
            });
        }
        Ok(TensorFatigueSample {
            stress_pa,
            total_mechanical_strain: *total_mechanical_strain,
        })
    }
}

impl sealed::Sealed for ResolvedJ2Plasticity {}
impl SourceBoundFatigueLaw for ResolvedJ2Plasticity {
    fn fatigue_source(&self) -> &ResolvedMaterialStatePoint {
        self.material()
    }
    fn fatigue_backbone(&self) -> &J2Plasticity {
        self.law()
    }
    fn accepted_plastic_strain<'a>(&self, state: &'a Self::State) -> &'a Voigt {
        &state.plastic_strain
    }
    fn validate_fatigue_history(
        &self,
        strain: &Voigt,
        state: &Self::State,
    ) -> Result<(), MaterialError> {
        self.law()
            .validate_overstress_history(strain, state)
            .map_err(history_error)
    }
}

impl sealed::Sealed for ResolvedVoceJ2 {}
impl SourceBoundFatigueLaw for ResolvedVoceJ2 {
    fn fatigue_source(&self) -> &ResolvedMaterialStatePoint {
        self.material()
    }
    fn fatigue_backbone(&self) -> &J2Plasticity {
        self.law().backbone()
    }
    fn accepted_plastic_strain<'a>(&self, state: &'a Self::State) -> &'a Voigt {
        &state.plastic_strain
    }
    fn validate_fatigue_history(
        &self,
        strain: &Voigt,
        state: &Self::State,
    ) -> Result<(), MaterialError> {
        self.law()
            .backbone()
            .validate_overstress_history(strain, state)
            .map_err(history_error)
    }
}

impl sealed::Sealed for ResolvedChabocheJ2 {}
impl SourceBoundFatigueLaw for ResolvedChabocheJ2 {
    fn fatigue_source(&self) -> &ResolvedMaterialStatePoint {
        self.material()
    }
    fn fatigue_backbone(&self) -> &J2Plasticity {
        self.law().backbone().backbone()
    }
    fn accepted_plastic_strain<'a>(&self, state: &'a Self::State) -> &'a Voigt {
        &state.plastic.plastic_strain
    }
    fn validate_fatigue_history(
        &self,
        strain: &Voigt,
        state: &Self::State,
    ) -> Result<(), MaterialError> {
        self.law()
            .validate_history(strain, state)
            .map_err(history_error)
    }
}

impl<const N: usize> sealed::Sealed for ResolvedChabocheBranches<N> {}
impl<const N: usize> SourceBoundFatigueLaw for ResolvedChabocheBranches<N> {
    fn fatigue_source(&self) -> &ResolvedMaterialStatePoint {
        self.material()
    }
    fn fatigue_backbone(&self) -> &J2Plasticity {
        self.law().backbone().backbone()
    }
    fn accepted_plastic_strain<'a>(&self, state: &'a Self::State) -> &'a Voigt {
        &state.plastic.plastic_strain
    }
    fn validate_fatigue_history(
        &self,
        strain: &Voigt,
        state: &Self::State,
    ) -> Result<(), MaterialError> {
        self.law()
            .validate_history(strain, state)
            .map_err(history_error)
    }
}

// Native reductions can borrow the same source owner without duplicating its
// material bundle. A reference cannot opt an unsupported law into this trait.
impl<L: SourceBoundFatigueLaw> sealed::Sealed for &L {}
impl<L: SourceBoundFatigueLaw> SourceBoundFatigueLaw for &L {
    fn fatigue_source(&self) -> &ResolvedMaterialStatePoint {
        L::fatigue_source(*self)
    }
    fn fatigue_backbone(&self) -> &J2Plasticity {
        L::fatigue_backbone(*self)
    }
    fn accepted_plastic_strain<'a>(&self, state: &'a Self::State) -> &'a Voigt {
        L::accepted_plastic_strain(*self, state)
    }
    fn validate_fatigue_history(
        &self,
        strain: &Voigt,
        state: &Self::State,
    ) -> Result<(), MaterialError> {
        L::validate_fatigue_history(*self, strain, state)
    }
}

fn history_error(error: impl core::fmt::Display) -> MaterialError {
    MaterialError::State {
        what: format!("native cyclic fatigue history: {error}"),
    }
}
