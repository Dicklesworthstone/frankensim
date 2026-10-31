//! Fatigue from the actual accepted inertial trajectory, including vibration
//! after all external forcing has stopped. Native Newmark mechanics and stress
//! recovery feed the same S-N and Findley owners as the quasistatic path.
//!
//! The wrapper owns the native run: its initial state and every successful
//! mechanical step are observed exactly once. There is no mutable run escape,
//! externally supplied stress or reconfiguration that can skip samples or
//! replace the bound mesh/material/temperature. Cloning checkpoints mechanics
//! and all unresolved histories together; matrices remain shared.
//!
//! Findley retains the in-plane shear-vector trajectory, searches declared
//! fixed planes and reports calibrated fatigue-limit utilization for a closed
//! window. Its result is not a life or damage estimate. Sampling must resolve
//! reversals and rotating stress paths. Both observers retain their fixed-state
//! elastic applicability; no plastic response, damage heat or stiffness loss
//! is inferred from this observation.

use super::response::TransientResponse3;
use super::{TetTransientRun3, TransientStep3};
use crate::fatigue::{
    TetFatigueError, TetFatigueGauge, TetFatigueHistoryConfig, TetFatigueObserver, TetFindleyGauge,
};
use crate::loads3::TetLoadCase3;
use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_material::fatigue::critical_plane::findley::{FindleyCalibration, FindleyCycleResult};
use fs_material::fatigue::{FatigueSummary, SnCurve};
use fs_material::state_point::ResolvedMaterialStatePoint;

/// One accepted native inertial interval and its S-N observations. Findley
/// window results remain available on the owning trajectory.
#[derive(Debug)]
pub struct TetTransientFatigueStep3 {
    /// Original integrator's solve and discrete work/energy report.
    pub mechanics: TransientStep3,
    /// Original endpoint stress, inertia, reactions and physical elastic energy.
    pub response: TransientResponse3,
    /// Fixed-plane signed normal stresses [Pa], in the original S-N gauge order.
    pub signed_stress_pa: Vec<f64>,
    /// Open-history damage; unresolved final half cycles are included by `finish`.
    pub fatigue: Vec<FatigueSummary>,
}

/// An inertial trajectory and its complete, source-bound fatigue history.
/// No accepted step can bypass the shared observer. A failed mechanical solve,
/// late gauge refusal, closure failure or cancellation commits neither side.
#[derive(Clone)]
pub struct TetTransientFatigueRun3<'g, 'm, 'p, 's> {
    run: TetTransientRun3<'g, 'm, 'p, 's>,
    observation: TetFatigueObserver,
}

impl<'g, 'm, 'p, 's> TetTransientRun3<'g, 'm, 'p, 's> {
    /// Attach S-N fatigue BEFORE the first mechanical advance and immediately
    /// sample the real initial stress. A thermal state is mandatory even for a
    /// zero-expansion specimen: temperature and material binding must match.
    pub fn with_fatigue(
        self,
        mesh_frame: ContentHash,
        gauges: Vec<TetFatigueGauge>,
        config: TetFatigueHistoryConfig,
        cx: &Cx<'_>,
    ) -> Result<TetTransientFatigueRun3<'g, 'm, 'p, 's>, TetFatigueError> {
        self.with_fatigue_and_findley(mesh_frame, gauges, Vec::new(), config, cx)
    }

    /// Attach S-N gauges, Findley gauges or both before the first advance.
    /// Initial displacement, velocity, acceleration and frozen thermal state
    /// remain the original native run's; initial stress is observed immediately.
    /// No rest or zero-stress sample is fabricated. Clone the returned wrapper
    /// to retain native mechanics and every unresolved cycle at a checkpoint.
    pub fn with_fatigue_and_findley(
        self,
        mesh_frame: ContentHash,
        gauges: Vec<TetFatigueGauge>,
        findley: Vec<TetFindleyGauge>,
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
            findley,
            config,
            cx,
        )?;
        let response =
            self.recover_response(mesh_frame, config.maximum_absolute_mandel_strain, cx)?;
        observation.accept(self.state.time_s, &response.stress, false, cx)?;
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

    /// Original immutable S-N calibration, including its declared scope.
    pub fn curve(&self, gauge: usize) -> Option<&SnCurve> {
        self.observation.curve(gauge)
    }

    pub fn summaries(&self) -> Vec<FatigueSummary> {
        self.observation.summaries()
    }

    pub fn observed_samples(&self) -> u64 {
        self.observation.observed_samples()
    }

    pub fn findley_results(&self) -> Vec<Option<FindleyCycleResult>> {
        self.observation.findley_results()
    }

    pub fn findley_calibration(&self, gauge: usize) -> Option<&FindleyCalibration> {
        self.observation.findley_calibration(gauge)
    }

    pub fn findley_material_source(&self, gauge: usize) -> Option<&ResolvedMaterialStatePoint> {
        self.observation.findley_material_source(gauge)
    }

    /// Close only the current accepted Findley window, preserving S-N history.
    pub fn close_findley_cycle(
        &mut self,
        cx: &Cx<'_>,
    ) -> Result<Vec<FindleyCycleResult>, TetFatigueError> {
        self.observation.close_findley_cycle(cx)
    }

    /// Discard the current Findley window and begin at the last actual accepted
    /// stress. Native state, time, S-N residue and global samples are unchanged.
    pub fn restart_findley_cycle(&mut self, cx: &Cx<'_>) -> Result<(), TetFatigueError> {
        self.observation.restart_findley_cycle(cx)
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
        self.advance_inner(load_next, None, h_s, false, cx)
    }

    /// Accept the endpoint and close Findley in the SAME native transaction.
    /// Late closure or circle-work refusal publishes neither mechanics nor any
    /// S-N/Findley history, allowing refinement or a corrected endpoint retry.
    pub fn advance_and_close_findley_cycle(
        &mut self,
        load_next: &TetLoadCase3<'_, '_>,
        h_s: f64,
        cx: &Cx<'_>,
    ) -> Result<TetTransientFatigueStep3, TetFatigueError> {
        self.advance_inner(load_next, None, h_s, true, cx)
    }

    /// Advance the original explicitly moving support endpoints. Support
    /// velocities/accelerations, their inertia and work stay in the native
    /// owner; observed stress includes their actual mechanical effects.
    pub fn advance_with_supports(
        &mut self,
        load_next: &TetLoadCase3<'_, '_>,
        prescribed_dofs_m: &[(usize, f64)],
        h_s: f64,
        cx: &Cx<'_>,
    ) -> Result<TetTransientFatigueStep3, TetFatigueError> {
        self.advance_inner(load_next, Some(prescribed_dofs_m), h_s, false, cx)
    }

    /// Moving-support endpoint and Findley closure are accepted together.
    pub fn advance_with_supports_and_close_findley_cycle(
        &mut self,
        load_next: &TetLoadCase3<'_, '_>,
        prescribed_dofs_m: &[(usize, f64)],
        h_s: f64,
        cx: &Cx<'_>,
    ) -> Result<TetTransientFatigueStep3, TetFatigueError> {
        self.advance_inner(load_next, Some(prescribed_dofs_m), h_s, true, cx)
    }

    fn advance_inner(
        &mut self,
        load_next: &TetLoadCase3<'_, '_>,
        prescribed_dofs_m: Option<&[(usize, f64)]>,
        h_s: f64,
        close_findley: bool,
        cx: &Cx<'_>,
    ) -> Result<TetTransientFatigueStep3, TetFatigueError> {
        self.observation
            .check_next(self.run.state.time_s + h_s, cx)?;
        if close_findley {
            self.observation.require_findley(cx)?;
        }
        let mut candidate = self.run.clone();
        let mechanics = match prescribed_dofs_m {
            Some(supports) => candidate.advance_with_supports(load_next, supports, h_s, cx)?,
            None => candidate.advance(load_next, h_s, cx)?,
        };
        let response = candidate.recover_response(
            self.observation.mesh_frame(),
            self.observation.maximum_strain(),
            cx,
        )?;
        let (signed_stress_pa, fatigue) =
            self.observation
                .accept(candidate.state.time_s, &response.stress, close_findley, cx)?;
        // Observer acceptance (including requested closure) is the last fallible
        // operation. Native mechanics and every gauge publish together.
        self.run = candidate;
        Ok(TetTransientFatigueStep3 {
            mechanics,
            response,
            signed_stress_pa,
            fatigue,
        })
    }

    /// Close the WHOLE observation, including S-N endpoint half cycles and every
    /// Findley window. Refusal leaves all histories and native motion untouched;
    /// success prevents more steps. Findley results remain available by accessor.
    pub fn finish(&mut self, cx: &Cx<'_>) -> Result<Vec<FatigueSummary>, TetFatigueError> {
        self.observation.finish(cx)
    }
}
