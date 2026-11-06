//! Plane-stress evolution with the original finite-viscoelastic history owner.
//!
//! A transverse solve borrows one accepted history and one physical interval.
//! Every local trial therefore sees the same viscous metrics and `dt`; a Newton
//! iteration is never a material time step. Only after transverse equilibrium
//! succeeds does the original material owner propose the updated history and
//! its physical work/storage/loss accounts. The accepted input stays unchanged.
//!
//! Exact incompressibility uses the original remainder history and native split
//! admission. It carries the pressure multiplier in the physical Piola stress,
//! with no bulk penalty or invented pressure heat. As in the original return
//! map, loading work is the frozen-metric loading increment; a spatial integrator
//! owns its separate difference from endpoint nodal work.

/// Prescribed natural-strain events with retained histories and transverse balance.
pub mod reference_change;
pub use reference_change::{
    ViscoIncompressiblePlaneStressReferenceChange, ViscoPlaneStressReferenceChange,
};

use super::{
    ViscoHyperelastic, ViscoRemainderState, ViscoRemainderStepResponse, ViscoState,
    ViscoStepResponse,
};
use crate::MaterialError;
use crate::hyper::incompressible_plane_stress::{
    IncompressiblePlaneStressState, incompressible_plane_stress_state_with_seed,
};
use crate::hyper::plane_stress::{
    PlaneStressSettings, PlaneStressState, plane_stress_state_with_seed,
};

/// Proposed finite-compliance plane-stress history and its converged reduction.
#[derive(Clone, Debug)]
pub struct ViscoPlaneStressStepResponse {
    /// Complete original material update, including physical stored energy,
    /// loading work and nonnegative branch relaxation loss.
    pub material: ViscoStepResponse,
    /// Reduced endpoint potential and in-plane Piola stress. The potential is
    /// the incremental value function, not the storage in `material.state`.
    /// The full physical first Piola stress is `material.first_piola_pa`.
    pub plane_stress: PlaneStressState,
}

/// Proposed exactly incompressible plane-stress remainder history.
#[derive(Clone, Debug)]
pub struct ViscoIncompressiblePlaneStressStepResponse {
    /// Original bulk-free material update. Its Piola stress is the remainder
    /// gradient and deliberately excludes the pressure multiplier.
    pub material: ViscoRemainderStepResponse,
    /// Physical endpoint stress, including `full_piola` and compression-positive
    /// `pressure_pa` in `P = P_remainder - p cof(F)`. Its energy field is an
    /// incremental potential; physical storage/work/loss remain in `material`.
    pub plane_stress: IncompressiblePlaneStressState,
}

impl ViscoHyperelastic {
    /// Propose one physical interval at the requested in-plane deformation.
    ///
    /// The original finite-compliance potential relaxes all three transverse
    /// tractions, then the original return map advances the retained history at
    /// precisely that full deformation and the same `dt_s`. `dt_s = 0` retains
    /// the existing elastic-jump convention without relaxation loss.
    ///
    /// An explicit `[F02,F12,F22]` can reuse a previously converged column. With
    /// `None`, the accepted transverse shears and natural-volume thickness
    /// `det(G)/det(F_in_plane)` seed the solve. This thickness is only a guess;
    /// finite-compliance relaxation remains free to change the actual volume.
    /// The shared reducer owns seed fallback, stability, domains and budgets.
    /// No partial history is published if reduction or final material accounting
    /// fails. The returned proposal is committed only by its caller.
    pub fn step_plane_stress_response(
        &self,
        accepted: &ViscoState,
        in_plane: &[f64; 4],
        dt_s: f64,
        initial_transverse: Option<[f64; 3]>,
        settings: PlaneStressSettings,
    ) -> Result<ViscoPlaneStressStepResponse, MaterialError> {
        let frozen = self.freeze(accepted, dt_s)?;
        let seed = initial_transverse.unwrap_or_else(|| {
            let determinant = in_plane[0] * in_plane[3] - in_plane[1] * in_plane[2];
            [
                accepted.deformation_gradient[2],
                accepted.deformation_gradient[5],
                self.reference_deformation().volume_ratio() / determinant,
            ]
        });
        let plane_stress = plane_stress_state_with_seed(&frozen, in_plane, seed, settings)?;
        let material = self.step_response(accepted, &plane_stress.deformation_gradient, dt_s)?;
        Ok(ViscoPlaneStressStepResponse {
            material,
            plane_stress,
        })
    }

    /// Propose one exact-volume interval using the native bulk-free history.
    ///
    /// `det(F)=det(G)` is imposed during every trial, independently of the finite
    /// bulk coefficient. The two transverse shears are relaxed and pressure is
    /// recovered from the normal traction equation. The existing remainder
    /// owner's split restrictions apply; this method does not turn a coupled
    /// Neo-Hookean full-energy history into an admitted remainder history.
    ///
    /// The default seed retains `[F02,F12]` from the accepted deformation. An
    /// explicit pair may instead reuse the spatial solver's final local state.
    /// Thickness is always computed from the new in-plane area and reference
    /// volume, never copied from the old state. All trials and the final original
    /// return map use the same accepted history and `dt_s`; errors leave it
    /// unchanged. Pressure contributes no fabricated branch dissipation.
    pub fn step_incompressible_plane_stress_response(
        &self,
        accepted: &ViscoRemainderState,
        in_plane: &[f64; 4],
        dt_s: f64,
        initial_shear: Option<[f64; 2]>,
        settings: PlaneStressSettings,
    ) -> Result<ViscoIncompressiblePlaneStressStepResponse, MaterialError> {
        let frozen = self.freeze_remainder(accepted, dt_s)?;
        let seed = initial_shear.unwrap_or([
            accepted.deformation_gradient[2],
            accepted.deformation_gradient[5],
        ]);
        let plane_stress =
            incompressible_plane_stress_state_with_seed(&frozen, in_plane, seed, settings)?;
        let material =
            self.step_remainder_response(accepted, &plane_stress.deformation_gradient, dt_s)?;
        Ok(ViscoIncompressiblePlaneStressStepResponse {
            material,
            plane_stress,
        })
    }
}
