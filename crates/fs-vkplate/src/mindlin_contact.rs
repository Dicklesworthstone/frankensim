//! A sourced, fixed-preload tangential contact on the plate body's full response.
//!
//! The footprint is one work-conjugate scalar shear coordinate. Its state lives
//! on the midpoint displacement grid: an interval advances it by
//! `du = T^T (w_next - w_previous) / 2`, with duration `k`. The returned average
//! traction therefore does exactly the work used by the centered body step.
//! All normal, dry-friction and retained shear reactions solve together.
//!
//! The original source supplies the held circular Hertz patch, elastic moduli,
//! coefficient and temperature. The normal loading is an externally maintained
//! boundary condition; this attachment does not solve its normal deformation or
//! introduce an additional normal reaction. The footprint must represent shear
//! along that patch's fixed first tangent. Two-axis/nonproportional slip,
//! changing preload, contact migration and thermal evolution are outside this
//! model. In particular it cannot replace a distributed hi-hat rim by assumption.

use fs_contact::interface_binding::mindlin::{
    SourceMindlinContact, SourceMindlinError, SourceMindlinState,
};
use fs_tribo::{TangentialLaw, TangentialSolveConfig, TangentialStepInput, TriboError};

use crate::friction_island::Controls;

/// Last accepted interval's retained contact work, in original source units.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct MindlinContactReport {
    /// Average effort entering the contact [N]; bodies receive its negative.
    pub traction_n: f64,
    /// Endpoint effort stored in the original constitutive checkpoint [N].
    pub endpoint_force_n: f64,
    pub slip_increment_m: f64,
    pub stored_change_j: f64,
    /// Microslip/gross-slip heat from the original law, not a numerical deficit.
    pub dissipated_j: f64,
    /// Signed work on the bodies [J]; positive when elastic memory returns work.
    pub mechanical_work_j: f64,
    /// `stored_change + dissipated + mechanical_work` [J].
    pub work_residual_j: f64,
}

/// An immutable source/footprint with its original retained shear checkpoint.
/// Clone this together with the instrument to preserve contact memory.
#[derive(Clone, Debug)]
pub struct MindlinContact {
    source: SourceMindlinContact,
    footprint: Vec<(usize, usize, f64)>,
    state: SourceMindlinState,
    solve: TangentialSolveConfig,
    report: MindlinContactReport,
    dissipated_j: f64,
}

#[derive(Clone)]
pub(crate) struct Candidate {
    state: SourceMindlinState,
    pub report: MindlinContactReport,
    pub residual: f64,
    dissipated_j: f64,
}

impl MindlinContact {
    /// Bind an explicit source checkpoint and numerical/force budgets.
    /// For initially unloaded shear, pass `source.initial_state()` explicitly.
    /// The caller declares this midpoint configuration as the checkpoint's
    /// physical reference; no rebasing of loaded material history is inferred.
    pub fn new(
        source: SourceMindlinContact,
        footprint: Vec<(usize, usize, f64)>,
        state: SourceMindlinState,
        solve: TangentialSolveConfig,
    ) -> Result<Self, SourceMindlinError> {
        if footprint.iter().any(|t| !t.2.is_finite()) || !footprint.iter().any(|t| t.2 != 0.0) {
            return Err(SourceMindlinError::Invalid(
                "Mindlin shear footprint must be finite and nonempty",
            ));
        }
        let contact = Self {
            source,
            footprint,
            state,
            solve,
            report: MindlinContactReport::default(),
            dissipated_j: 0.0,
        };
        // Reuse the source's identity, temperature, load and numerical admission.
        // A checkpoint's endpoint effort need not equal the next accepted mean
        // reaction, so assembly does not impose a fictitious zero-slip step.
        contact.source.trial_prescribed_for_solve(
            &contact.state,
            &contact.input(
                0.0,
                contact.source.normal_source().envelope().maximum_step_s,
            ),
        )?;
        Ok(contact)
    }

    pub fn source(&self) -> &SourceMindlinContact {
        &self.source
    }
    pub fn state(&self) -> &SourceMindlinState {
        &self.state
    }
    pub fn footprint(&self) -> &[(usize, usize, f64)] {
        &self.footprint
    }
    pub fn report(&self) -> MindlinContactReport {
        self.report
    }
    pub fn stored_energy_j(&self) -> f64 {
        self.source.stored_energy_j(&self.state)
    }
    /// Heat accumulated by this attachment after its supplied checkpoint [J].
    pub fn dissipated_j(&self) -> f64 {
        self.dissipated_j
    }

    fn input(&self, increment_m: f64, k: f64) -> TangentialStepInput {
        TangentialStepInput {
            normal_n: self.source.normal_response().normal_force_n,
            slip_increment_m: [increment_m, 0.0],
            dt_s: k,
            contact_temperature_k: Some(self.source.temperature_kelvin()),
            nominal_pressure_pa: None,
            compliance_m_per_n: 0.0,
            solve: self.solve,
        }
    }

    /// Force and derivative with respect to the full centered nodal path.
    /// Every trial borrows the original history; it allocates no new memory.
    /// Only the converged body reaction is subject to the accepted-force cap:
    /// the free-motion response can be larger while the coupled root is valid.
    pub(crate) fn trial(
        &self,
        centered_path: f64,
        k: f64,
    ) -> Result<(f64, f64), SourceMindlinError> {
        let trial = self
            .source
            .trial_prescribed_for_solve(&self.state, &self.input(0.5 * centered_path, k))?;
        Ok((trial.average_force_n, 0.5 * trial.average_tangent_n_per_m))
    }

    /// Construct history once from the actual rounded nodal increment. The
    /// caller publishes this only after every contact and body check succeeds.
    pub(crate) fn candidate(
        &self,
        increment: f64,
        force: f64,
        k: f64,
        controls: Controls,
    ) -> Result<Candidate, SourceMindlinError> {
        if !force.is_finite() {
            return Err(SourceMindlinError::Invalid(
                "Mindlin applied reaction must be finite",
            ));
        }
        if force.abs() > self.solve.maximum_force_n {
            return Err(SourceMindlinError::Tribo(
                TriboError::OutsideApplicability {
                    field: "maximum_force_n",
                    value: force.abs(),
                    minimum: 0.0,
                    maximum: self.solve.maximum_force_n,
                },
            ));
        }
        let step = self
            .source
            .advance(&self.state, &self.input(increment, k))
            .map_err(SourceMindlinError::Tribo)?;
        let before = self.stored_energy_j();
        let after = self.source.stored_energy_j(&step.state);
        let force_scale = controls.force_tolerance
            + controls.relative_tolerance * force.abs().max(step.traction_n[0].abs());
        let slip_scale = controls.displacement_tolerance
            + controls.relative_tolerance * increment.abs().max(step.slip_increment_m[0].abs());
        let mechanical_work_j = -force * increment;
        let work_residual_j = step.stored_change_j + step.dissipated_j + mechanical_work_j;
        let work_scale = force_scale * increment.abs()
            + controls.displacement_tolerance * force.abs().max(step.traction_n[0].abs())
            + controls.relative_tolerance
                * before
                    .max(after)
                    .max(step.dissipated_j)
                    .max(mechanical_work_j.abs())
            + controls.force_tolerance * controls.displacement_tolerance;
        let dissipated_j = self.dissipated_j + step.dissipated_j;
        if [
            before,
            after,
            force_scale,
            slip_scale,
            work_scale,
            mechanical_work_j,
            work_residual_j,
            dissipated_j,
            step.stored_change_j,
            step.dissipated_j,
            step.closure_residual_j,
            step.compatibility_residual_m[0],
            step.slip_increment_m[0],
        ]
        .iter()
        .any(|x| !x.is_finite())
            || before < 0.0
            || after < 0.0
            || step.dissipated_j < 0.0
            || force_scale <= 0.0
            || slip_scale <= 0.0
            || work_scale <= 0.0
            || step.traction_n[1] != 0.0
            || step.slip_increment_m[1] != 0.0
        {
            return Err(SourceMindlinError::Invalid(
                "Mindlin candidate has invalid storage, work or scalar kinematics",
            ));
        }
        let residual = ((force - step.traction_n[0]).abs() / force_scale)
            .max((step.slip_increment_m[0] - increment).abs() / slip_scale)
            .max(step.compatibility_residual_m[0].abs() / slip_scale)
            .max(work_residual_j.abs() / work_scale)
            .max(step.closure_residual_j.abs() / work_scale)
            .max((after - before - step.stored_change_j).abs() / work_scale);
        if !residual.is_finite() {
            return Err(SourceMindlinError::Invalid(
                "non-finite Mindlin acceptance residual",
            ));
        }
        let report = MindlinContactReport {
            traction_n: force,
            endpoint_force_n: step.state.endpoint_force_n(),
            slip_increment_m: increment,
            stored_change_j: step.stored_change_j,
            dissipated_j: step.dissipated_j,
            mechanical_work_j,
            work_residual_j,
        };
        Ok(Candidate {
            state: step.state,
            report,
            residual,
            dissipated_j,
        })
    }

    pub(crate) fn publish(&mut self, candidate: Candidate) {
        self.state = candidate.state;
        self.report = candidate.report;
        self.dissipated_j = candidate.dissipated_j;
    }
}
