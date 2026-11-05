//! Work-conjugate stress averages on the original straight strain increment.
//!
//! The native cohesive owner integrates its complete vector traction, including
//! irreversible first initiation, nonproportional mixed-mode loading, rupture,
//! unloading and normal reclosure. The positive in-plane elastic remainder is
//! averaged exactly. For fixed mode tensors `A_i` and width `h`,
//! `mean_stress = C_parallel * (old_strain + strain)/2 + sum A_i * mean_t_i`.
//! Its tensor double contraction with the strain change equals the original
//! volumetric path work, independently of the endpoint stress.
//!
//! The optional tangent differentiates this average with respect to the NEW
//! strain at fixed accepted history: `C_parallel/2 + h sum A_i Kmean_ij A_j`.
//! Tensor-shear columns carry their factor of two. First-onset BK sensitivities
//! stay with the original cohesive derivative; it explicitly declines a
//! classical derivative at unresolved nonsmooth/degenerate paths.
//!
//! These allocation-free point operations provide the force/work primitive for
//! a dynamic or discrete-gradient owner. They do not advance time, supply mass,
//! add numerical damping, resolve unseen strain reversals, or establish a
//! finite-element time-integration scheme. Coefficients, plane and width stay
//! fixed over each supplied strain path.

use super::{
    BrittleCrackBand, BrittleError, BrittleState, BrittleStep, add_outer, contract, finite,
    finite_tangent, finite_voigt,
};
use crate::cohesive::CohesiveStep;
use crate::{MaterialAdmissibility, SmallStrainLaw, Tangent6, Voigt};
use fs_evidence::ModelCard;

/// Original endpoint and a distinct work-conjugate stress for its whole path.
#[derive(Clone, Debug, PartialEq)]
pub struct BrittlePathStep<'a> {
    /// Same endpoint stress, tangents, history and energies as `trial`.
    pub endpoint: BrittleStep<'a>,
    /// Analytic stress average on the affine TOTAL tensor-strain path, Pa.
    pub mean_stress_pa: Voigt,
    /// `mean_stress : (new_strain - old_strain)`, J/m³.
    pub path_work_j_m3: f64,
    /// Path-mean work minus the original material work, J/m³. This signed
    /// arithmetic defect is neither fracture dissipation nor generated heat.
    pub path_work_defect_j_m3: f64,
}

/// Borrowed work-conjugate force adapter for the native small-strain assembly.
///
/// `strain` is the NEW endpoint and `State` is the original accepted
/// [`BrittleState`]. The force response is the exact path-mean stress, whereas
/// `update_state` retains the original endpoint history. This adapter neither
/// averages history nor supplies time integration, mass or a Newton fallback.
///
/// A classical derivative can be absent even for a valid zero-length path.
/// The infallible trait reports this with a nonfinite tangent; the checked
/// [`Self::trial`] preserves the native `Option`. A nonlinear consumer must
/// explicitly admit its own numerical or semismooth fallback. No endpoint
/// tangent is substituted for the derivative of the mean force.
#[derive(Clone, Copy, Debug)]
pub struct BrittleMeanSmallStrain<'a> {
    owner: &'a BrittleCrackBand,
}

impl BrittleCrackBand {
    /// Borrow the original law for assembly of its exact path-mean stress.
    /// See [`BrittleMeanSmallStrain`] for the explicit nonsmooth-tangent boundary.
    #[must_use]
    pub fn mean_small_strain_law(&self) -> BrittleMeanSmallStrain<'_> {
        BrittleMeanSmallStrain { owner: self }
    }

    /// Original endpoint trial plus the exact vector stress average. This
    /// force-only entry avoids evaluating the optional path-mean derivative.
    ///
    /// One native cohesive trial supplies both the endpoint and the average;
    /// no damage evolution or path quadrature is repeated. Owner identity,
    /// both Frobenius strain limits and the original material-energy checks
    /// apply. The mean's work must additionally agree with endpoint energy
    /// accounting within `settings().energy_tolerance_j_m3`.
    pub fn trial_with_mean_stress<'a>(
        &'a self,
        strain: &Voigt,
        state: &BrittleState<'a>,
    ) -> Result<BrittlePathStep<'a>, BrittleError> {
        let jump = self.jump_for_trial(strain, state)?;
        let (step, mean) = state
            .point
            .trial_with_mean_traction(jump, self.area_energy_tolerance)?;
        self.finish_path(strain, state, step, mean)
    }

    /// Same original trial and average, plus
    /// `d(mean_stress)/d(new_tensor_strain)` in Pa at fixed accepted history.
    /// This generally nonsymmetric derivative differs from the returned
    /// endpoint tangent, even for an intact elastic step (whose mean tangent
    /// is HALF the elastic stiffness).
    ///
    /// `None` preserves the native cohesive owner's explicit refusal of a
    /// classical derivative at a kink or unresolved crossing sensitivity. A
    /// consumer must use its admitted fallback; no endpoint or fabricated zero
    /// tangent is substituted. The physical trial/mean can still be valid.
    pub fn trial_with_mean_stress_tangent<'a>(
        &'a self,
        strain: &Voigt,
        state: &BrittleState<'a>,
    ) -> Result<(BrittlePathStep<'a>, Option<Tangent6>), BrittleError> {
        let jump = self.jump_for_trial(strain, state)?;
        let (step, mean, local_tangent) = state
            .point
            .trial_with_mean_traction_tangent(jump, self.area_energy_tolerance)?;
        let path = self.finish_path(strain, state, step, mean)?;
        let tangent = local_tangent
            .map(|local| self.path_tangent_from_modes(local))
            .transpose()?;
        Ok((path, tangent))
    }

    fn finish_path<'a>(
        &'a self,
        strain: &Voigt,
        state: &BrittleState<'a>,
        step: CohesiveStep<'a>,
        mean_traction_pa: [f64; 3],
    ) -> Result<BrittlePathStep<'a>, BrittleError> {
        let midpoint = core::array::from_fn(|i| 0.5 * state.strain[i] + 0.5 * strain[i]);
        let mut mean_stress_pa = self.parallel_stress(&midpoint);
        for (mode, traction) in self.modes.iter().zip(mean_traction_pa) {
            for i in 0..6 {
                mean_stress_pa[i] += mode[i] * traction;
            }
        }
        finite_voigt(&mean_stress_pa)?;
        let endpoint = self.finish_trial(strain, state, step)?;
        let change = core::array::from_fn(|i| strain[i] - state.strain[i]);
        let path_work_j_m3 = finite(contract(&mean_stress_pa, &change))?;
        let path_work_defect_j_m3 = finite(path_work_j_m3 - endpoint.energy.work)?;
        if path_work_defect_j_m3.abs() > self.settings.energy_tolerance_j_m3 {
            return Err(BrittleError::EnergyDefect(path_work_defect_j_m3));
        }
        Ok(BrittlePathStep {
            endpoint,
            mean_stress_pa,
            path_work_j_m3,
            path_work_defect_j_m3,
        })
    }

    fn path_tangent_from_modes(&self, local: [[f64; 3]; 3]) -> Result<Tangent6, BrittleError> {
        let mut tangent = self
            .parallel_tangent
            .map(|row| row.map(|value| 0.5 * value));
        for (i, row) in local.iter().enumerate() {
            for (j, value) in row.iter().enumerate() {
                add_outer(
                    &mut tangent,
                    &self.modes[i],
                    &self.modes[j],
                    finite(self.settings.band_width_m * value)?,
                );
            }
        }
        finite_tangent(&tangent)?;
        Ok(tangent)
    }
}

impl<'a> BrittleMeanSmallStrain<'a> {
    /// One original cohesive trial supplies endpoint history, path-mean stress,
    /// work audit and its optional classical derivative at fixed old history.
    /// Physical refusals and valid paths without a derivative remain distinct.
    pub fn trial(
        &self,
        strain: &Voigt,
        state: &BrittleState<'a>,
    ) -> Result<(BrittlePathStep<'a>, Option<Tangent6>), BrittleError> {
        self.owner.trial_with_mean_stress_tangent(strain, state)
    }
}

impl<'a> SmallStrainLaw for BrittleMeanSmallStrain<'a> {
    type State = BrittleState<'a>;

    fn initial_state(&self) -> Self::State {
        self.owner.initial_state()
    }

    fn stress(&self, strain: &Voigt, state: &Self::State) -> Voigt {
        self.owner
            .trial_with_mean_stress(strain, state)
            .map(|path| path.mean_stress_pa)
            .unwrap_or([f64::NAN; 6])
    }

    fn tangent(&self, strain: &Voigt, state: &Self::State) -> Tangent6 {
        self.trial(strain, state)
            .ok()
            .and_then(|(_, tangent)| tangent)
            .unwrap_or([[f64::NAN; 6]; 6])
    }

    fn update_state(&self, strain: &Voigt, state: &Self::State) -> Self::State {
        self.owner.small_strain_law().update_state(strain, state)
    }

    fn admissibility(&self) -> MaterialAdmissibility {
        self.owner.small_strain_law().admissibility()
    }

    fn card(&self) -> ModelCard {
        let mut card = self.owner.small_strain_law().card();
        card.name = "material.brittle-fixed-plane-path-mean".into();
        card.assumptions.push(
            "exact stress average on one straight strain increment at fixed accepted history; original endpoint state and work retained".into(),
        );
        card.known_failures.push(
            "native path-mean derivative can be absent at held contact or damage switches, including zero-length Newton starts; consumer owns an explicitly admitted fallback".into(),
        );
        card.known_failures.push(
            "mean force stress is not the endpoint stress; no time integration, mass or intermediate strain reversals are supplied".into(),
        );
        card.assumptions.sort_unstable();
        card.known_failures.sort_unstable();
        card
    }
}
