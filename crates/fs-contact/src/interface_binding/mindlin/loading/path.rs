//! Source-admitted simultaneous Hertz indentation and compliant scalar shear.
//!
//! Normal load, footprint and reversible normal energy come from the original
//! material/geometry owner. The shared analytic tribology owner resolves shear
//! history and body recoil; this adapter supplies no independent coefficient or
//! normal-force override. Temperature and all admitted source conditions remain
//! those of the immutable isothermal contact.

/// Pressure-weighted irreversible sliding on this same original admitted path.
#[path = "path/sliding.rs"]
pub mod sliding;

use super::{SourceHertzMindlinContact, SourceHertzMindlinState, invalid, match_original};
use crate::interface_binding::mindlin::SourceMindlinError;
use crate::normal_patch::PointNormalPatchReceipt;
use fs_tribo::mindlin::loading::path::{MindlinPathError, MindlinPathStep};
use fs_tribo::{TangentialSolveConfig, TangentialStep, TangentialStepInput};

/// One affine indentation path and released scalar tangential body increment.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SourceMindlinPathInput {
    /// Absolute endpoint indentation in the original circular geometry [m].
    pub indentation_m: f64,
    /// Released motion along the attachment's fixed first tangent [m].
    pub free_slip_increment_m: f64,
    /// Physical path duration, admitted against the original rate envelope [s].
    pub dt_s: f64,
    /// Mechanical displacement response to resisting average shear effort [m/N].
    /// Zero prescribes the shear path; positive values resolve body recoil.
    pub compliance_m_per_n: f64,
    /// Original bounded constitutive solve controls.
    pub solve: TangentialSolveConfig,
}

/// Accepted-source candidate containing separate normal and shear work ports.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceMindlinPathStep {
    /// Actual shear increment, average effort, retained memory, storage and heat.
    pub tangential: TangentialStep<SourceHertzMindlinState>,
    /// Original endpoint load, pressure, footprint, authority and uncertainty.
    pub normal_response: PointNormalPatchReceipt,
    /// Signed work from the original conservative Hertz potential [J].
    pub normal_work_j: f64,
    /// Change in that same normal potential [J], separate from shear storage.
    pub normal_stored_change_j: f64,
}

/// Original source admission or constitutive path refusal, including cancellation.
#[derive(Clone, Debug, PartialEq)]
pub enum SourceMindlinPathError {
    Source(SourceMindlinError),
    Path(MindlinPathError),
}

impl From<SourceMindlinError> for SourceMindlinPathError {
    fn from(error: SourceMindlinError) -> Self {
        Self::Source(error)
    }
}
impl From<MindlinPathError> for SourceMindlinPathError {
    fn from(error: MindlinPathError) -> Self {
        Self::Path(error)
    }
}
impl core::fmt::Display for SourceMindlinPathError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source-backed Hertz--Mindlin path: {self:?}")
    }
}
impl std::error::Error for SourceMindlinPathError {}

impl SourceHertzMindlinContact {
    /// Resolve simultaneous normal loading and shear with optional body recoil.
    ///
    /// Indentation and the solved shear displacement are affine over the
    /// interval. The original normal evaluator derives the endpoint force;
    /// normal force itself follows the nonlinear Hertz law. The mechanical
    /// owner supplies the scalar compliance, giving
    /// `actual_slip + compliance * average_effort = free_slip`.
    ///
    /// This retains the original source identities, frozen temperature and
    /// complete shear memory through loading, unloading and separation. It
    /// returns a candidate only: no source or accepted state is mutated on
    /// refusal, cancellation, or successful evaluation. Normal dynamics,
    /// changing material conditions and rotating tangential frames require
    /// their own mechanical/history-transfer owners.
    pub fn advance_coupled_path(
        &self,
        state: &SourceHertzMindlinState,
        input: &SourceMindlinPathInput,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SourceMindlinPathStep, SourceMindlinPathError> {
        let (response, native_input) = self.admit_coupled_path(state, input, cancelled)?;
        let path = self
            .law
            .advance_coupled_path(&state.constitutive, &native_input, cancelled)?;
        self.bind_coupled_path(state, input, response, path, cancelled)
    }

    // Shared source admission and native inputs; neither entry repeats a solve.
    fn admit_coupled_path(
        &self,
        state: &SourceHertzMindlinState,
        input: &SourceMindlinPathInput,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(PointNormalPatchReceipt, TangentialStepInput), SourceMindlinPathError> {
        if cancelled() {
            return Err(MindlinPathError::Cancelled.into());
        }
        if !self.owns_state(state) {
            return Err(
                invalid("simultaneous Mindlin path requires original source identities").into(),
            );
        }
        let normal = self.source.normal_source();
        normal
            .validate_kinematics(state.indentation_m, 0.0, input.dt_s)
            .map_err(SourceMindlinError::from)?;
        let rate = (input.indentation_m - state.indentation_m) / input.dt_s;
        normal
            .validate_kinematics(state.indentation_m, rate, input.dt_s)
            .map_err(SourceMindlinError::from)?;
        let response = normal
            .evaluate(
                "source-mindlin/simultaneous-normal-shear-path",
                input.indentation_m,
                rate,
                input.dt_s,
            )
            .map_err(SourceMindlinError::from)?;
        let native_input = TangentialStepInput {
            normal_n: response.normal_force_n,
            slip_increment_m: [input.free_slip_increment_m, 0.0],
            dt_s: input.dt_s,
            contact_temperature_k: Some(self.temperature_kelvin()),
            nominal_pressure_pa: None,
            compliance_m_per_n: input.compliance_m_per_n,
            solve: input.solve,
        };
        Ok((response, native_input))
    }

    fn bind_coupled_path(
        &self,
        state: &SourceHertzMindlinState,
        input: &SourceMindlinPathInput,
        response: PointNormalPatchReceipt,
        path: MindlinPathStep,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SourceMindlinPathStep, SourceMindlinPathError> {
        match_original(&path.tangential.state, &response)?;
        if cancelled() {
            return Err(MindlinPathError::Cancelled.into());
        }
        let change = response.reversible_energy_j - state.normal_energy_j;
        let candidate = path.tangential;
        Ok(SourceMindlinPathStep {
            tangential: TangentialStep {
                traction_n: candidate.traction_n,
                slip_increment_m: candidate.slip_increment_m,
                state: SourceHertzMindlinState {
                    normal_model_identity: state.normal_model_identity,
                    friction_state_identity: state.friction_state_identity,
                    constitutive: candidate.state,
                    indentation_m: input.indentation_m,
                    normal_energy_j: response.reversible_energy_j,
                },
                regime: candidate.regime,
                stored_change_j: candidate.stored_change_j,
                dissipated_j: candidate.dissipated_j,
                closure_residual_j: candidate.closure_residual_j,
                compatibility_residual_m: candidate.compatibility_residual_m,
                iterations: candidate.iterations,
            },
            normal_response: response,
            normal_work_j: change,
            normal_stored_change_j: change,
        })
    }
}
