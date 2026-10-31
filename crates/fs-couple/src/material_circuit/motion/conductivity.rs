//! Nonlinear bulk/contact conduction in the original reciprocal field/body descriptor.
//!
//! Changing conductor material changes discharge, field force and actual motion.
//! The shared conductivity owner resolves E/J at that SAME descriptor trajectory;
//! there is no stationary-circuit prediction followed by a mechanical correction.

use super::{FieldMotionCandidate, FieldMotionInput, FieldMotionStep, MaterialFieldMotion};
use crate::material_circuit::conductivity::{ConductivityDriver, ConductivitySample,
    ConductivitySolve, ConductorResponse, solve_native};
use crate::material_circuit::conductivity::contact::{ContactConductivitySample, ContactResponse};
use crate::material_circuit::contacts::{ContactHeat, SeriesElectricalContact};
use crate::material_circuit::{MaterialCircuitError, MaterialCircuitSolve, MaterialResistor,
    ResistorHeat};

/// Original field/body work and the material response at its solved E/J.
pub struct FieldConductivityStep {
    /// Authoritative native electrical, mechanical, support and Joule account.
    pub field: FieldMotionStep,
    /// Complete source evaluations at actual solved branch fields.
    pub evaluated: Vec<ConductivitySample>,
    /// Complete interface sources at their OWN solved drops and current densities.
    /// Empty when no electrical operating-response bindings are supplied.
    pub evaluated_contacts: Vec<ContactConductivitySample>,
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
    fn contact_heats(&self) -> &[ContactHeat] { &self.step().contacts }
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
        self.conductivity_with_contact_response_candidate(seeds, Vec::new(), sources,
            None, input, controls, solve, cancelled)
    }

    /// Close bulk E/J response with electrically independent series contacts in
    /// the SAME reciprocal field/body step. Original contact cards, selected
    /// claims, areas and branch assignments are retained; their supplied thermal
    /// context remains fixed within this inner iteration. Operating-dependent
    /// contacts require [`Self::conductivity_with_contact_response_candidate`].
    #[allow(clippy::too_many_arguments)]
    pub fn conductivity_with_contacts_candidate(
        &self, seeds: Vec<MaterialResistor>, contacts: Vec<SeriesElectricalContact>,
        sources: &[ConductorResponse<'_>], input: FieldMotionInput<'_>,
        controls: MaterialCircuitSolve, solve: ConductivitySolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<FieldConductivityCandidate, MaterialCircuitError> {
        self.conductivity_with_contact_response_candidate(seeds, contacts, sources,
            None, input, controls, solve, cancelled)
    }

    /// Jointly close bulk E/J, individual contact V/I/J and reciprocal field
    /// motion using the original native descriptor. Every constitutive trial
    /// restarts from the same accepted charge/flux, gap, momentum and support.
    /// Contact drops are excluded from bulk field queries. Every contact uses
    /// its own area and explicitly oriented source coordinates; no branch-wide
    /// substitute, stationary circuit predictor or mechanical correction is used.
    ///
    /// `Some` requires one binding per original interface, including fixed laws.
    /// `None` admits only electrically independent contact sources. Resistance
    /// and constitutive-power residuals gate both bulk and contacts. Their native
    /// losses remain separate and are the only Joule heat; constitutive defects
    /// do not become heat. Actual material responses remain inspectable.
    ///
    /// This heats stationary conductors/interfaces around the moving lossless
    /// gap. Conductor/contact geometry, field coefficients and body inertia stay
    /// fixed. Optional support retains its frozen thermal state. Motional
    /// conductor transport, interface storage, arcing, opening contacts, thermal
    /// field entropy and live thermoelastic supports need their own owners.
    #[allow(clippy::too_many_arguments)]
    pub fn conductivity_with_contact_response_candidate(
        &self, seeds: Vec<MaterialResistor>, contacts: Vec<SeriesElectricalContact>,
        sources: &[ConductorResponse<'_>], contact_sources: Option<&[ContactResponse<'_>]>,
        input: FieldMotionInput<'_>, controls: MaterialCircuitSolve,
        solve: ConductivitySolve, cancelled: impl FnMut() -> bool,
    ) -> Result<FieldConductivityCandidate, MaterialCircuitError> {
        let trial = solve_native(self.circuit(), seeds, sources, contacts, contact_sources,
            controls, solve, |applied, interfaces, check| {
                self.step_with_contacts_candidate(applied, interfaces, input, controls, check)
            }, cancelled)?;
        let (motion, field) = trial.candidate.into_parts();
        Ok(FieldConductivityCandidate { motion, step: FieldConductivityStep { field,
            evaluated: trial.evaluated, evaluated_contacts: trial.evaluated_contacts,
            iterations: trial.iterations,
            source_evaluations: trial.source_evaluations } })
    }
}
