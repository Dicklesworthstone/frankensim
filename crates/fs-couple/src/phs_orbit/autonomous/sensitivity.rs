//! Physical parameter changes of complete self-excited pHS cycles.

use super::super::{
    PhsFeedbackParameterAction, PhsFeedbackStateJacobian, PhsOrbitError, PhsStorageHessianAction,
    sideband::PhsStorageParameterAction, variational_feedback::feedback_tangent_problem,
};
use super::AutonomousPortFeedback;
use fs_orbit::sensitivity::{
    AutonomousMeanObjective, AutonomousOrbitSensitivity, ForcedMeanObjectiveControls,
    ForcedSensitivityControls, ForcedSensitivityError, PreparedAutonomousSensitivity,
};
use fs_orbit::{HbOrbit, HbScales};
use fs_phs::{PhsError, PhsStructureParameterDirection, PortHamiltonian};

/// Complete direct parameter direction in the unchanged physical state chart.
/// Every source callback is deterministic/read-only; memory remains explicit.
pub struct PhsAutonomousParameterDirection<'a> {
    pub structure: PhsStructureParameterDirection<'a>,
    /// `partial_p grad H` at fixed original state, including all memory terms.
    /// None holds the complete storage law fixed.
    pub storage: Option<&'a PhsStorageParameterAction<'a>>,
    /// Required TOTAL fixed-state feedback input derivative. It receives both
    /// original conjugate flows and their exact direct parameter derivative.
    /// Explicitly write zeros when the original law has no direct dependence.
    pub feedback: &'a PhsFeedbackParameterAction<'a>,
}

#[derive(Clone, Debug, PartialEq)]
pub enum PhsAutonomousSensitivityError {
    Source(PhsOrbitError),
    Physical(PhsError),
    Sensitivity(ForcedSensitivityError),
}
impl core::fmt::Display for PhsAutonomousSensitivityError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "pHS autonomous harmonic sensitivity refused: {self:?}")
    }
}
impl std::error::Error for PhsAutonomousSensitivityError {}
impl From<PhsOrbitError> for PhsAutonomousSensitivityError {
    fn from(value: PhsOrbitError) -> Self {
        Self::Source(value)
    }
}
impl From<PhsError> for PhsAutonomousSensitivityError {
    fn from(value: PhsError) -> Self {
        Self::Physical(value)
    }
}
impl From<ForcedSensitivityError> for PhsAutonomousSensitivityError {
    fn from(value: ForcedSensitivityError) -> Self {
        Self::Sensitivity(value)
    }
}
type Result<T> = core::result::Result<T, PhsAutonomousSensitivityError>;

mod adjoint;
pub use adjoint::PreparedPhsAutonomousAdjoint;

/// One source-bound exact tangent/factor for many physical parameter directions.
///
/// Borrows the original storage, J/R/G and external feedback unchanged. The
/// generic fs-orbit owner solves both frequency and every Fourier coefficient;
/// the shared complete feedback Jacobian supplies the original state chain.
pub struct PreparedPhsAutonomousSensitivity<'a> {
    system: &'a PortHamiltonian,
    feedback: &'a dyn AutonomousPortFeedback,
    prepared: PreparedAutonomousSensitivity,
}

impl<'a> PreparedPhsAutonomousSensitivity<'a> {
    /// Prepare from the original storage-owned complete continuous tangent.
    /// `phase_reference` uses the original complete packed autonomous HB chart.
    /// All original states, including constitutive/controller memories, remain.
    pub fn new(
        system: &'a PortHamiltonian,
        feedback: &'a dyn AutonomousPortFeedback,
        feedback_jacobian: &PhsFeedbackStateJacobian<'_>,
        orbit: &HbOrbit,
        phase_reference: &[f64],
        controls: ForcedSensitivityControls,
        scales: &HbScales,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        Self::prepare(
            system,
            feedback,
            feedback_jacobian,
            orbit,
            phase_reference,
            controls,
            scales,
            None,
            cancelled,
        )
    }

    /// The same source with an explicitly authored COMPLETE storage Hessian.
    /// The feedback state Jacobian already includes its flow/Hessian chain.
    pub fn new_analytic(
        system: &'a PortHamiltonian,
        feedback: &'a dyn AutonomousPortFeedback,
        feedback_jacobian: &PhsFeedbackStateJacobian<'_>,
        orbit: &HbOrbit,
        phase_reference: &[f64],
        controls: ForcedSensitivityControls,
        scales: &HbScales,
        hessian: &PhsStorageHessianAction<'_>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        Self::prepare(
            system,
            feedback,
            feedback_jacobian,
            orbit,
            phase_reference,
            controls,
            scales,
            Some(hessian),
            cancelled,
        )
    }

    fn prepare(
        system: &'a PortHamiltonian,
        feedback: &'a dyn AutonomousPortFeedback,
        feedback_jacobian: &PhsFeedbackStateJacobian<'_>,
        orbit: &HbOrbit,
        phase_reference: &[f64],
        mut controls: ForcedSensitivityControls,
        scales: &HbScales,
        hessian: Option<&PhsStorageHessianAction<'_>>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        if cancelled() {
            return Err(ForcedSensitivityError::Cancelled.into());
        }
        let (n, m) = (system.state_dim() as u128, system.port_dim() as u128);
        // Complete physical-feedback tangent plus all source scratch. This is
        // reserved in addition to the numerical owner's bounded AFT/LU payload.
        let bytes = 8_u128.saturating_mul(
            2_u128
                .saturating_mul(n.saturating_mul(n))
                .saturating_add(n.saturating_mul(m))
                .saturating_add(6 * n)
                .saturating_add(6 * m),
        );
        if bytes > controls.maximum_payload_bytes as u128 {
            return Err(ForcedSensitivityError::Budget {
                quantity: "physical feedback adapter payload",
                requested: bytes,
                cap: controls.maximum_payload_bytes as u128,
            }
            .into());
        }
        controls.maximum_payload_bytes -= bytes as usize;
        let problem = feedback_tangent_problem(system, feedback, feedback_jacobian, hessian)?;
        let result = PreparedAutonomousSensitivity::new(
            &problem,
            orbit,
            phase_reference,
            controls,
            scales,
            |point, out| {
                // This original adapter's port IS s I. No general frequency
                // load or descriptor is silently replaced by an ODE port.
                let d = system.state_dim();
                out.fill(0.0);
                for harmonic in 1..=point.harmonics {
                    for state in 0..d {
                        let re = d + (harmonic - 1) * 2 * d + state;
                        let im = re + d;
                        out[re] = -(harmonic as f64) * point.coefficients[im];
                        out[im] = harmonic as f64 * point.coefficients[re];
                    }
                }
                out.iter().all(|v| v.is_finite())
            },
            &mut cancelled,
        );
        if problem.source_refused() {
            return Err(PhsOrbitError::SourceRefused.into());
        }
        Ok(Self {
            system,
            feedback,
            prepared: result?,
        })
    }

    /// Original numerical diagnostics and complete residual-direction access.
    pub fn prepared(&self) -> &PreparedAutonomousSensitivity {
        &self.prepared
    }

    /// Mean observation of the complete original self-excited cycle. The last
    /// callback output is the complete frequency partial at fixed phase and
    /// state, including `f_omega - (t/omega)*f_t` for physical-time observations.
    /// The numerical owner checks both value and full coefficient/frequency
    /// covector on its two original AFT grids. Direct physical parameter
    /// dependence remains separate from the solved waveform/frequency chain.
    pub fn mean_objective(
        &mut self,
        controls: ForcedMeanObjectiveControls,
        observation: impl FnMut(f64, &[f64], &mut f64, &mut [f64], &mut f64) -> bool,
        cancelled: impl FnMut() -> bool,
    ) -> Result<AutonomousMeanObjective> {
        self.prepared
            .mean_objective(controls, observation, cancelled)
            .map_err(Into::into)
    }

    /// Mean original Hamiltonian and exact full-state gradient, including all
    /// constitutive memory. Its explicit frequency partial is zero; frequency
    /// changes still affect the solved waveform through the bordered adjoint.
    /// A material-dependent H retains its separate explicit parameter partial.
    pub fn mean_stored_energy(
        &mut self,
        controls: ForcedMeanObjectiveControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<AutonomousMeanObjective> {
        let system = self.system;
        self.mean_objective(
            controls,
            |_time, state, value, gradient, frequency| {
                *value = system.hamiltonian(state);
                system.storage().gradient(state, gradient);
                *frequency = 0.0;
                value.is_finite() && gradient.iter().all(|v| v.is_finite())
            },
            cancelled,
        )
    }

    /// Differentiate stiffness/storage, loss, interconnection, reciprocal port
    /// maps and feedback parameters in one physical direction.
    ///
    /// At fixed state, `e=grad H`, `h_p=partial_p grad H`,
    /// `y_p=G^T h_p+G_p^T e`, and the authored feedback callback gives the TOTAL
    /// `u_p` including its direct parameter and `y_p` chain exactly once. The
    /// original source-bound owner evaluates
    /// `(J-R)h_p+(J_p-R_p)e+G_p u+G u_p` with compensated contractions. The
    /// cached total state Jacobian owns the separate unknown `x_p` chain.
    ///
    /// The same generic phase-bordered factor supplies omega_p and every
    /// original state/memory Fourier derivative. The source is never perturbed,
    /// refitted, retuned or reduced. Output is per physical parameter unit;
    /// parameter_scale affects conditioning and residual acceptance only.
    ///
    /// This is a local finite-order sensitivity at fixed state/phase charts;
    /// it does not prove stability, infinite harmonic resolution or physical
    /// constitutive validity. The shared explicit work budget excludes source
    /// callbacks; their count is bounded by the two declared AFT grids. Failure
    /// or cancellation publishes no partial result and leaves the factor reusable.
    pub fn solve_parameter(
        &mut self,
        parameter_scale: f64,
        direction: PhsAutonomousParameterDirection<'_>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AutonomousOrbitSensitivity> {
        if cancelled() {
            return Err(ForcedSensitivityError::Cancelled.into());
        }
        with_parameter_direction(self.system, self.feedback, direction, |action| {
            self.prepared
                .solve_island_direction(parameter_scale, action, &mut cancelled)
        })
    }
}

// Both complete waveform derivatives and scalar objective adjoints traverse
// this same original J/R/G/storage/feedback chain, with no duplicated law.
fn with_parameter_direction<T>(
    system: &PortHamiltonian,
    feedback: &dyn AutonomousPortFeedback,
    direction: PhsAutonomousParameterDirection<'_>,
    run: impl FnOnce(
        &mut dyn FnMut(f64, &[f64], &mut [f64]) -> bool,
    ) -> core::result::Result<T, ForcedSensitivityError>,
) -> Result<T> {
    let physical = system.prepare_parameter_direction(direction.structure)?;
    let (n, m) = (system.state_dim(), system.port_dim());
    let mut effort = vec![f64::NAN; n];
    let mut partial_effort = vec![0.0; n];
    let mut flows = vec![0.0; m];
    let mut partial_flows = vec![0.0; m];
    let mut input = vec![f64::NAN; m];
    let mut partial_input = vec![f64::NAN; m];
    let mut source_refused = false;
    let mut physical_error = None;
    let result = run(&mut |_time, state, out| {
        effort.fill(f64::NAN);
        system.storage().gradient(state, &mut effort);
        if effort.iter().any(|v| !v.is_finite())
            || system.output_from_effort_into(&effort, &mut flows).is_err()
        {
            source_refused = true;
            return false;
        }
        input.fill(f64::NAN);
        if !feedback.input(state, &flows, &mut input) || input.iter().any(|v| !v.is_finite()) {
            source_refused = true;
            return false;
        }
        if let Some(storage) = direction.storage {
            partial_effort.fill(f64::NAN);
            if !storage(state, &mut partial_effort) || partial_effort.iter().any(|v| !v.is_finite())
            {
                return false;
            }
        }
        if let Err(error) =
            physical.output_from_effort_into(&effort, &partial_effort, &mut partial_flows)
        {
            physical_error = Some(error);
            return false;
        }
        partial_input.fill(f64::NAN);
        if !(direction.feedback)(state, &flows, &partial_flows, &mut partial_input)
            || partial_input.iter().any(|v| !v.is_finite())
        {
            return false;
        }
        match physical.vector_field_from_effort_into(
            &effort,
            &partial_effort,
            &input,
            &partial_input,
            out,
        ) {
            Ok(()) => true,
            Err(error) => {
                physical_error = Some(error);
                false
            }
        }
    });
    if source_refused {
        return Err(PhsOrbitError::SourceRefused.into());
    }
    if let Some(error) = physical_error {
        return Err(error.into());
    }
    result.map_err(Into::into)
}

#[cfg(test)]
mod tests;
