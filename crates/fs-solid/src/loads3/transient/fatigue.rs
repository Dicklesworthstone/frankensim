//! Fatigue from the actual accepted inertial trajectory, including vibration
//! after all external forcing has stopped. Native Newmark mechanics, stress
//! recovery and the shared rainflow/Miner owner remain unchanged.
//!
//! The wrapper owns the native run: its initial state and every successful
//! mechanical step are observed exactly once. There is no mutable run escape,
//! externally supplied stress, or reconfiguration API that could skip samples
//! or replace the bound mesh/material/temperature during an observation.
//! Cloning checkpoints both mechanics and unresolved cycles together. Matrices
//! remain shared; trial storage is bounded by state size and history limits.
//!
//! This is the same fixed-temperature, elastic, fixed-plane empirical S-N
//! baseline as `crate::fatigue`. Sampling must resolve stress reversals through
//! time refinement. It is not a multiaxial, plastic, thermal-fatigue, fracture,
//! stiffness-degradation or validated component-life model. Damage is not heat.

use super::response::TransientResponse3;
use super::{TetTransientRun3, TransientStep3};
use crate::fatigue::{
    TetFatigueError, TetFatigueGauge, TetFatigueHistoryConfig, TetFatigueObserver,
};
use crate::loads3::TetLoadCase3;
use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_material::fatigue::{FatigueSummary, SnCurve};

/// One accepted native inertial interval and its fatigue observations.
#[derive(Debug)]
pub struct TetTransientFatigueStep3 {
    /// Original integrator's solve and discrete work/energy report.
    pub mechanics: TransientStep3,
    /// Original endpoint stress, inertia, reactions and physical elastic energy.
    pub response: TransientResponse3,
    /// Fixed-plane signed normal stresses [Pa], in the original gauge order.
    pub signed_stress_pa: Vec<f64>,
    /// Open-history damage; unresolved final half cycles are included by `finish`.
    pub fatigue: Vec<FatigueSummary>,
}

/// An inertial trajectory and its complete, source-bound fatigue history.
/// No accepted step can bypass the shared observer. A failed mechanical solve,
/// late gauge refusal, stress-domain error or cancellation commits neither side.
#[derive(Clone)]
pub struct TetTransientFatigueRun3<'g, 'm, 'p, 's> {
    run: TetTransientRun3<'g, 'm, 'p, 's>,
    observation: TetFatigueObserver,
}

impl<'g, 'm, 'p, 's> TetTransientRun3<'g, 'm, 'p, 's> {
    /// Attach fatigue BEFORE the first mechanical advance and immediately
    /// sample the real initial stress. A thermal state is mandatory even for a
    /// zero-expansion specimen: its actual absolute temperature and material
    /// binding must match each calibration. No temperature is invented.
    ///
    /// The consumed run retains its original matrices, supports, forces,
    /// displacement, velocity, acceleration and native budgets. Inertial steps
    /// are never replaced by static equilibrium. Use `Clone` on the returned
    /// wrapper to continue a checkpoint with its unresolved cycles intact.
    pub fn with_fatigue(
        self,
        mesh_frame: ContentHash,
        gauges: Vec<TetFatigueGauge>,
        config: TetFatigueHistoryConfig,
        cx: &Cx<'_>,
    ) -> Result<TetTransientFatigueRun3<'g, 'm, 'p, 's>, TetFatigueError> {
        cx.checkpoint().map_err(|_| TetFatigueError::Cancelled)?;
        if self.state.step != 0 {
            return Err(TetFatigueError::Invalid(
                "fatigue must begin at the original initial state; past stresses cannot be reconstructed",
            ));
        }
        let thermal = self
            .thermal_strains
            .as_ref()
            .ok_or(TetFatigueError::Invalid(
                "fatigue requires the native run's explicit frozen temperature/material state",
            ))?;
        let mut observation = TetFatigueObserver::new(
            self.problem,
            thermal.field(),
            mesh_frame,
            gauges,
            config,
            cx,
        )?;
        let response =
            self.recover_response(mesh_frame, config.maximum_absolute_mandel_strain, cx)?;
        observation.accept(self.state.time_s, &response.stress, cx)?;
        Ok(TetTransientFatigueRun3 {
            run: self,
            observation,
        })
    }
}

impl<'g, 'm, 'p, 's> TetTransientFatigueRun3<'g, 'm, 'p, 's> {
    /// Read-only original trajectory. An independently cloned run cannot be
    /// attached back to this observation, overwrite its state or skip a sample.
    pub fn mechanics(&self) -> &TetTransientRun3<'g, 'm, 'p, 's> {
        &self.run
    }

    /// Original immutable calibration, including its caller-declared scope.
    pub fn curve(&self, gauge: usize) -> Option<&SnCurve> {
        self.observation.curve(gauge)
    }

    pub fn summaries(&self) -> Vec<FatigueSummary> {
        self.observation.summaries()
    }

    /// Actual accepted stress/reactions with this observer's frozen frame and
    /// strain bound. This read does not add another fatigue sample.
    pub fn recover_response(&self, cx: &Cx<'_>) -> Result<TransientResponse3, TetFatigueError> {
        Ok(self.run.recover_response(
            self.observation.mesh_frame(),
            self.observation.maximum_strain(),
            cx,
        )?)
    }

    /// Advance stationary-support mechanics and every gauge as one transaction.
    pub fn advance(
        &mut self,
        load_next: &TetLoadCase3<'_, '_>,
        h_s: f64,
        cx: &Cx<'_>,
    ) -> Result<TetTransientFatigueStep3, TetFatigueError> {
        self.observation
            .check_next(self.run.state.time_s + h_s, cx)?;
        let mut candidate = self.run.clone();
        let mechanics = candidate.advance(load_next, h_s, cx)?;
        self.accept_candidate(candidate, mechanics, cx)
    }

    /// Advance the original explicitly moving support endpoints. Support
    /// velocities/accelerations, their inertia and work stay in the native
    /// owner; the observed stress includes their actual mechanical effects.
    pub fn advance_with_supports(
        &mut self,
        load_next: &TetLoadCase3<'_, '_>,
        prescribed_dofs_m: &[(usize, f64)],
        h_s: f64,
        cx: &Cx<'_>,
    ) -> Result<TetTransientFatigueStep3, TetFatigueError> {
        self.observation
            .check_next(self.run.state.time_s + h_s, cx)?;
        let mut candidate = self.run.clone();
        let mechanics = candidate.advance_with_supports(load_next, prescribed_dofs_m, h_s, cx)?;
        self.accept_candidate(candidate, mechanics, cx)
    }

    fn accept_candidate(
        &mut self,
        candidate: TetTransientRun3<'g, 'm, 'p, 's>,
        mechanics: TransientStep3,
        cx: &Cx<'_>,
    ) -> Result<TetTransientFatigueStep3, TetFatigueError> {
        let response = candidate.recover_response(
            self.observation.mesh_frame(),
            self.observation.maximum_strain(),
            cx,
        )?;
        let (signed_stress_pa, fatigue) =
            self.observation
                .accept(candidate.state.time_s, &response.stress, cx)?;
        // Observer acceptance is the last fallible operation. Native mechanics
        // and all gauges therefore publish together, including cancellation.
        self.run = candidate;
        Ok(TetTransientFatigueStep3 {
            mechanics,
            response,
            signed_stress_pa,
            fatigue,
        })
    }

    /// Close the WHOLE observation once, not a processing chunk. Refusal leaves
    /// all histories open and mechanics untouched; success prevents more steps.
    pub fn finish(&mut self, cx: &Cx<'_>) -> Result<Vec<FatigueSummary>, TetFatigueError> {
        self.observation.finish(cx)
    }
}
