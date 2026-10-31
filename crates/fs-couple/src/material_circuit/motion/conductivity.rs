//! Nonlinear bulk conduction in the original reciprocal field/body descriptor.
//!
//! Changing conductor material changes discharge, field force and actual motion.
//! The shared conductivity owner resolves E/J at that SAME descriptor trajectory;
//! there is no stationary-circuit prediction followed by a mechanical correction.

use super::{FieldMotionCandidate, FieldMotionInput, FieldMotionStep, MaterialFieldMotion};
use crate::material_circuit::conductivity::{ConductivityDriver, ConductivitySample,
    ConductivitySolve, ConductorResponse, solve_native};
use crate::material_circuit::{MaterialCircuitError, MaterialCircuitSolve, MaterialResistor,
    ResistorHeat};

/// Original field/body work and the material response at its solved E/J.
pub struct FieldConductivityStep {
    /// Authoritative native electrical, mechanical, support and Joule account.
    pub field: FieldMotionStep,
    /// Complete source evaluations at actual solved branch fields.
    pub evaluated: Vec<ConductivitySample>,
    /// All original descriptor solves, including discarded constitutive trials.
    pub iterations: usize,
    /// Original bundle refreshes, excluding the initial seed-admission replay.
    pub source_evaluations: usize,
}

/// Joint charge/flux, gap and momentum, uncommitted until all sources agree.
pub struct FieldConductivityCandidate {
    motion: MaterialFieldMotion,
    step: FieldConductivityStep,
}
impl FieldConductivityCandidate {
    /// Proposed original coupled state; field and motion cannot be split.
    pub const fn motion(&self) -> &MaterialFieldMotion { &self.motion }
    /// Native energy transfer and separate constitutive residuals.
    pub const fn step(&self) -> &FieldConductivityStep { &self.step }
    /// Publish this endpoint alongside any finite thermal receivers.
    pub fn into_parts(self) -> (MaterialFieldMotion, FieldConductivityStep) {
        (self.motion, self.step)
    }
}

impl ConductivityDriver for FieldMotionCandidate {
    fn resistor_heats(&self) -> &[ResistorHeat] { &self.step().resistors }
}

impl MaterialFieldMotion {
    /// Close scalar conductor response, circuit storage and field/body motion.
    /// Every material iteration starts at the SAME accepted charge/flux, gap,
    /// momentum and support history. The original descriptor supplies branch
    /// voltage/current; the original pinned material bundle supplies resistance.
    /// Both resistance and constitutive power residuals must meet their gates.
    /// Native Joule work is authoritative and counted once; iteration defects
    /// never become heat or an independently applied force.
    ///
    /// Electrical axes are explicit existing SI source coordinates. All other
    /// material coordinates stay fixed here. Conductor geometry is stationary;
    /// gap motion retains the original uniform, constant-area field reduction,
    /// fixed field coefficient, guided inertia and optional conservative support.
    /// Series contacts require their own jointly refreshed interface law and
    /// refuse this bulk-only extension. No motional conductor transport, Hall
    /// law, breakdown, hysteresis, field entropy or convergence guarantee is added.
    pub fn conductivity_candidate(
        &self, sources: &[ConductorResponse<'_>], input: FieldMotionInput<'_>,
        controls: MaterialCircuitSolve, solve: ConductivitySolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<FieldConductivityCandidate, MaterialCircuitError> {
        self.conductivity_with_materials_candidate(self.circuit().resistors().to_vec(),
            sources, input, controls, solve, cancelled)
    }

    /// The same E/J-motion closure at a proposed thermal/material context.
    /// Seeds must retain each original complete source bundle, geometry and
    /// constitutive basis. Their temperature/pressure context is retained while
    /// actual E/J is resolved again. An outer thermal owner may use native heat
    /// to close its own enthalpy residual, then publish all histories together.
    #[allow(clippy::too_many_arguments)]
    pub fn conductivity_with_materials_candidate(
        &self, seeds: Vec<MaterialResistor>, sources: &[ConductorResponse<'_>],
        input: FieldMotionInput<'_>, controls: MaterialCircuitSolve,
        solve: ConductivitySolve, cancelled: impl FnMut() -> bool,
    ) -> Result<FieldConductivityCandidate, MaterialCircuitError> {
        let trial = solve_native(self.circuit(), seeds, sources, Vec::new(), None, controls, solve,
            |applied, _, check| self.step_candidate(applied, input, controls, check), cancelled)?;
        let (motion, field) = trial.candidate.into_parts();
        Ok(FieldConductivityCandidate { motion, step: FieldConductivityStep { field,
            evaluated: trial.evaluated, iterations: trial.iterations,
            source_evaluations: trial.source_evaluations } })
    }
}
