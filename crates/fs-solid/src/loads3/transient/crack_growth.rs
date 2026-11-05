//! Mode-I Paris observation on the original accepted inertial trajectory.
//!
//! The same bounded chronological cycle observer as the quasistatic owner
//! receives initial and accepted native stresses, including free vibration.
//! No zero-load static replacement, skipped native step, mutable-run escape or
//! external stress trace can alter its history. All fixed-temperature,
//! nominal-stress, bounded-geometry and simple-cycle restrictions remain in force.

use super::response::TransientResponse3;
use super::{TetTransientRun3, TransientStep3};
use crate::fatigue::crack_growth::{
    TetCrackGrowthConfig, TetCrackGrowthError, TetCrackGrowthGauge, TetCrackGrowthObserver,
};
use crate::loads3::TetLoadCase3;
use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_material::fatigue::crack_growth::{
    CrackGeometry, CrackGrowthForecast, CrackGrowthState, CrackGrowthStep, ParisLaw,
};
use fs_material::state_point::ResolvedMaterialStatePoint;

#[derive(Debug)]
pub struct TetTransientCrackGrowthStep3 {
    pub mechanics: TransientStep3,
    pub response: TransientResponse3,
    pub normal_stress_pa: Vec<f64>,
    pub cracks: Vec<CrackGrowthState>,
    pub completed_cycles: Option<Vec<CrackGrowthStep>>,
}

/// Cloning retains mechanical motion, open cycle and crack length together.
/// Failure commits neither candidate mechanics nor any observer history.
#[derive(Clone)]
pub struct TetTransientCrackGrowthRun3<'g, 'm, 'p, 's> {
    run: TetTransientRun3<'g, 'm, 'p, 's>,
    observation: TetCrackGrowthObserver,
}

impl<'g, 'm, 'p, 's> TetTransientRun3<'g, 'm, 'p, 's> {
    /// Attach before the first native advance and observe the ACTUAL initial
    /// stress immediately. Begin at a physical cycle extremum. An explicit
    /// frozen temperature/material state is required, including at zero strain.
    pub fn with_crack_growth(
        self,
        mesh_frame: ContentHash,
        gauges: Vec<TetCrackGrowthGauge>,
        config: TetCrackGrowthConfig,
        cx: &Cx<'_>,
    ) -> Result<TetTransientCrackGrowthRun3<'g, 'm, 'p, 's>, TetCrackGrowthError> {
        cx.checkpoint()
            .map_err(|_| TetCrackGrowthError::Cancelled)?;
        if self.state.step != 0 {
            return Err(TetCrackGrowthError::Invalid(
                "crack growth must begin at the original initial state; prior cycles cannot be reconstructed",
            ));
        }
        let thermal = self
            .thermal_strains
            .as_ref()
            .ok_or(TetCrackGrowthError::Invalid(
                "crack growth requires the native run's explicit frozen temperature/material state",
            ))?;
        let mut observation = TetCrackGrowthObserver::new(
            self.problem,
            thermal.field(),
            mesh_frame,
            gauges,
            config,
            cx,
        )?;
        let response = self.recover_response(mesh_frame, observation.maximum_strain(), cx)?;
        observation.accept(self.state.time_s, &response.stress, false, cx)?;
        Ok(TetTransientCrackGrowthRun3 {
            run: self,
            observation,
        })
    }
}

impl<'g, 'm, 'p, 's> TetTransientCrackGrowthRun3<'g, 'm, 'p, 's> {
    pub fn mechanics(&self) -> &TetTransientRun3<'g, 'm, 'p, 's> {
        &self.run
    }
    pub fn states(&self) -> Vec<CrackGrowthState> {
        self.observation.states()
    }
    pub fn observed_samples(&self) -> u64 {
        self.observation.observed_samples()
    }
    pub fn law(&self, gauge: usize) -> Option<&ParisLaw> {
        self.observation.law(gauge)
    }
    pub fn geometry(&self, gauge: usize) -> Option<&CrackGeometry> {
        self.observation.geometry(gauge)
    }
    pub fn material_source(&self, gauge: usize) -> Option<&ResolvedMaterialStatePoint> {
        self.observation.material_source(gauge)
    }

    /// Conditional exposure forecast from a fully sampled, still unconsumed
    /// cycle. Identical future extrema are a caller assumption, not a prediction
    /// of the dynamics. No native step, cycle window or crack state changes.
    pub fn cycles_to_length(
        &self,
        gauge: usize,
        target_length_m: f64,
        cx: &Cx<'_>,
    ) -> Result<CrackGrowthForecast, TetCrackGrowthError> {
        self.observation
            .cycles_to_length(gauge, target_length_m, cx)
    }

    /// Recover current native response without adding an observer sample.
    pub fn recover_response(&self, cx: &Cx<'_>) -> Result<TransientResponse3, TetCrackGrowthError> {
        Ok(self.run.recover_response(
            self.observation.mesh_frame(),
            self.observation.maximum_strain(),
            cx,
        )?)
    }

    pub fn advance(
        &mut self,
        load: &TetLoadCase3<'_, '_>,
        h_s: f64,
        cx: &Cx<'_>,
    ) -> Result<TetTransientCrackGrowthStep3, TetCrackGrowthError> {
        self.advance_inner(load, None, h_s, false, cx)
    }

    /// The endpoint and its completed physical cycle publish together.
    pub fn advance_and_close_cycle(
        &mut self,
        load: &TetLoadCase3<'_, '_>,
        h_s: f64,
        cx: &Cx<'_>,
    ) -> Result<TetTransientCrackGrowthStep3, TetCrackGrowthError> {
        self.advance_inner(load, None, h_s, true, cx)
    }

    /// Moving supports retain their original inertia and work in native motion.
    pub fn advance_with_supports(
        &mut self,
        load: &TetLoadCase3<'_, '_>,
        prescribed_dofs_m: &[(usize, f64)],
        h_s: f64,
        cx: &Cx<'_>,
    ) -> Result<TetTransientCrackGrowthStep3, TetCrackGrowthError> {
        self.advance_inner(load, Some(prescribed_dofs_m), h_s, false, cx)
    }

    pub fn advance_with_supports_and_close_cycle(
        &mut self,
        load: &TetLoadCase3<'_, '_>,
        prescribed_dofs_m: &[(usize, f64)],
        h_s: f64,
        cx: &Cx<'_>,
    ) -> Result<TetTransientCrackGrowthStep3, TetCrackGrowthError> {
        self.advance_inner(load, Some(prescribed_dofs_m), h_s, true, cx)
    }

    fn advance_inner(
        &mut self,
        load: &TetLoadCase3<'_, '_>,
        prescribed_dofs_m: Option<&[(usize, f64)]>,
        h_s: f64,
        close_cycle: bool,
        cx: &Cx<'_>,
    ) -> Result<TetTransientCrackGrowthStep3, TetCrackGrowthError> {
        self.observation
            .check_next(self.run.state.time_s + h_s, cx)?;
        let mut candidate = self.run.clone();
        let mechanics = match prescribed_dofs_m {
            Some(supports) => candidate.advance_with_supports(load, supports, h_s, cx)?,
            None => candidate.advance(load, h_s, cx)?,
        };
        let response = candidate.recover_response(
            self.observation.mesh_frame(),
            self.observation.maximum_strain(),
            cx,
        )?;
        let (normal_stress_pa, completed_cycles) =
            self.observation
                .accept(candidate.state.time_s, &response.stress, close_cycle, cx)?;
        self.run = candidate;
        Ok(TetTransientCrackGrowthStep3 {
            mechanics,
            response,
            normal_stress_pa,
            cracks: self.states(),
            completed_cycles,
        })
    }

    pub fn close_cycle(
        &mut self,
        cx: &Cx<'_>,
    ) -> Result<Vec<CrackGrowthStep>, TetCrackGrowthError> {
        self.observation.close_cycle(cx)
    }

    pub fn finish(&mut self, cx: &Cx<'_>) -> Result<Vec<CrackGrowthState>, TetCrackGrowthError> {
        self.observation.finish(cx)
    }
}
