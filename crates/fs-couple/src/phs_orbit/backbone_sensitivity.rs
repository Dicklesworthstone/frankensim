//! Original material, interconnection and load derivatives of a backbone.

use super::{
    PhsOrbitError, PhsStorageHessianAction, autonomous_problem, bad,
    sideband::PhsStorageParameterAction, tangent::StorageTangent,
};
use fs_orbit::sensitivity::{
    AutonomousOrbitSensitivity, BackboneSensitivityAnchor, BackboneSensitivityDirection,
    ForcedSensitivityControls, ForcedSensitivityError, PreparedBackboneSensitivity,
};
use fs_orbit::{HbOrbit, HbScales};
use fs_phs::{PhsError, PhsStructureParameterDirection, PortHamiltonian};

/// Complete direct physical parameter changes at fixed original state.
///
/// This is a direction, not a replacement source: the unchanged Hamiltonian,
/// J/R/G and bias still supply every base evaluation and state tangent. All
/// explicit material and coupled coordinates remain in the original chart.
#[derive(Clone, Copy, Default)]
pub struct PhsBackboneParameterDirection<'a> {
    pub structure: PhsStructureParameterDirection<'a>,
    /// Complete `partial_p grad H(x)` at FIXED x. None fixes the storage law.
    /// Do not include `Hxx*x_p`; the cached original tangent already owns it.
    pub storage: Option<&'a PhsStorageParameterAction<'a>>,
    /// Derivative of every constant physical port input. None fixes the bias.
    pub bias: Option<&'a [f64]>,
}

#[derive(Clone, Debug, PartialEq)]
pub enum PhsBackboneSensitivityError {
    Source(PhsOrbitError),
    Physical(PhsError),
    Sensitivity(ForcedSensitivityError),
}
impl core::fmt::Display for PhsBackboneSensitivityError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "pHS backbone sensitivity refused: {self:?}")
    }
}
impl std::error::Error for PhsBackboneSensitivityError {}
impl From<PhsOrbitError> for PhsBackboneSensitivityError {
    fn from(value: PhsOrbitError) -> Self {
        Self::Source(value)
    }
}
impl From<PhsError> for PhsBackboneSensitivityError {
    fn from(value: PhsError) -> Self {
        Self::Physical(value)
    }
}
impl From<ForcedSensitivityError> for PhsBackboneSensitivityError {
    fn from(value: ForcedSensitivityError) -> Self {
        Self::Sensitivity(value)
    }
}
type Result<T> = core::result::Result<T, PhsBackboneSensitivityError>;

/// Frequency and every original harmonic's parameter derivative on a backbone.
///
/// The original pHS is borrowed unchanged. Its complete storage tangent and
/// `s I` port feed fs-orbit's exact AFT linearization. That owner retains every
/// physical balance row, plus full-waveform phase and selected first-harmonic
/// amplitude closures, in one reusable rectangular Householder QR factor.
///
/// A passive source need not have a conservative periodic family. No resistance
/// is removed and no balancing drive is invented: the base and differentiated
/// original equations must actually be compatible. A new loss direction can
/// therefore refuse even when its J/R/G derivatives are structurally valid.
/// Compatibility is measured numerical evidence, not a proof of conservation,
/// infinite harmonic resolution, branch identity or Floquet stability.
pub struct PreparedPhsBackboneSensitivity<'a> {
    system: &'a PortHamiltonian,
    bias: &'a [f64],
    prepared: PreparedBackboneSensitivity,
}

impl<'a> PreparedPhsBackboneSensitivity<'a> {
    /// Prepare with the original storage-owned COMPLETE continuous Hessian.
    ///
    /// `anchor.coordinate` is an original physical state; no first-coordinate
    /// assumption or state permutation is introduced. The complete packed
    /// autonomous `phase_reference` fixes one phase chart for all directions.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        system: &'a PortHamiltonian,
        bias: &'a [f64],
        orbit: &HbOrbit,
        anchor: BackboneSensitivityAnchor,
        phase_reference: &[f64],
        controls: ForcedSensitivityControls,
        scales: &HbScales,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        Self::prepare(
            system,
            bias,
            orbit,
            anchor,
            phase_reference,
            controls,
            scales,
            None,
            cancelled,
        )
    }

    /// The same original source with an authored complete storage Hessian.
    /// Unavailable, incomplete or nonfinite tangent entries refuse; no numerical
    /// probing or quadratic substitute replaces a declared constitutive law.
    #[allow(clippy::too_many_arguments)]
    pub fn new_analytic(
        system: &'a PortHamiltonian,
        bias: &'a [f64],
        orbit: &HbOrbit,
        anchor: BackboneSensitivityAnchor,
        phase_reference: &[f64],
        controls: ForcedSensitivityControls,
        scales: &HbScales,
        hessian: &PhsStorageHessianAction<'_>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        Self::prepare(
            system,
            bias,
            orbit,
            anchor,
            phase_reference,
            controls,
            scales,
            Some(hessian),
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn prepare(
        system: &'a PortHamiltonian,
        bias: &'a [f64],
        orbit: &HbOrbit,
        anchor: BackboneSensitivityAnchor,
        phase_reference: &[f64],
        mut controls: ForcedSensitivityControls,
        scales: &HbScales,
        hessian: Option<&PhsStorageHessianAction<'_>>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        if cancelled() {
            return Err(ForcedSensitivityError::Cancelled.into());
        }
        let (n, m) = (system.state_dim(), system.port_dim());
        if n == 0 || bias.len() != m || bias.iter().any(|v| !v.is_finite()) {
            return Err(bad("constant bias must cover every original physical port").into());
        }
        // Reserve the original storage tangent/cache and all adapter scratch
        // before constructing it, separately from the numerical owner's budget.
        let (nd, md) = (n as u128, m as u128);
        let bytes = 8_u128.saturating_mul(
            2_u128
                .saturating_mul(nd.saturating_mul(nd))
                .saturating_add(8 * nd)
                .saturating_add(4 * md),
        );
        if bytes > controls.maximum_payload_bytes as u128 {
            return Err(ForcedSensitivityError::Budget {
                quantity: "physical backbone adapter payload",
                requested: bytes,
                cap: controls.maximum_payload_bytes as u128,
            }
            .into());
        }
        controls.maximum_payload_bytes -= bytes as usize;
        let mut problem = autonomous_problem(system, bias)?;
        if let Some(action) = hessian {
            problem.tangent = StorageTangent::hessian(n, action);
        }
        let result = PreparedBackboneSensitivity::new(
            &problem,
            orbit,
            anchor,
            phase_reference,
            controls,
            scales,
            |point, out| {
                // This complete original ODE has precisely sI and a constant
                // bias. At fixed Fourier phase, only i*n*omega*X changes.
                out.fill(0.0);
                for harmonic in 1..=point.harmonics {
                    for state in 0..n {
                        let re = n + (harmonic - 1) * 2 * n + state;
                        let im = re + n;
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
            bias,
            prepared: result?,
        })
    }

    /// Shared exact-factor diagnostics and complete residual-direction access.
    pub fn prepared(&self) -> &PreparedBackboneSensitivity {
        &self.prepared
    }

    /// Differentiate original material, interconnection, resistance, physical
    /// port maps and constant loads, plus an independently authored amplitude.
    ///
    /// For `e=grad H` and `e_p=partial_p grad H`, the original physical owner
    /// applies `(J-R)e_p+(J_p-R_p)e+G_p u+G u_p` in one compensated contraction.
    /// The cached complete state tangent supplies the separate `x_p` chain;
    /// fs-orbit solves every coefficient and `omega_p` together. J_p must be
    /// skew and R_p symmetric, as admitted by the original parameter owner.
    ///
    /// Set `direction.first_harmonic_amplitude_derivative=0` for a physical
    /// parameter at fixed amplitude. Set it to one with the default zero source
    /// direction to obtain the nonlinear backbone slope and harmonic changes.
    /// `parameter_scale` affects only conditioning and residual acceptance;
    /// the result remains per physical parameter unit. Both original AFT grids
    /// must satisfy the differentiated balances, phase and amplitude constraint.
    ///
    /// No model perturbation, refitting or new nonlinear orbit solve occurs.
    /// Source callbacks are read-only and bounded by the owner's two grids.
    /// Failure/cancellation publishes nothing and leaves the factor reusable.
    pub fn solve_parameter(
        &mut self,
        direction: BackboneSensitivityDirection,
        parameter: PhsBackboneParameterDirection<'_>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AutonomousOrbitSensitivity> {
        if cancelled() {
            return Err(ForcedSensitivityError::Cancelled.into());
        }
        let (n, m) = (self.system.state_dim(), self.system.port_dim());
        if parameter
            .bias
            .is_some_and(|values| values.len() != m || values.iter().any(|v| !v.is_finite()))
        {
            return Err(bad("bias derivative must cover every original physical port").into());
        }
        let physical = self
            .system
            .prepare_parameter_direction(parameter.structure)?;
        let mut effort = vec![f64::NAN; n];
        let mut partial_effort = vec![0.0; n];
        let zero_bias = vec![0.0; m];
        let partial_bias = parameter.bias.unwrap_or(&zero_bias);
        let mut source_refused = false;
        let mut physical_error = None;
        let result = self.prepared.solve_island_direction(
            direction,
            |_time, state, out| {
                effort.fill(f64::NAN);
                self.system.storage().gradient(state, &mut effort);
                if effort.iter().any(|v| !v.is_finite()) {
                    source_refused = true;
                    return false;
                }
                if let Some(storage) = parameter.storage {
                    partial_effort.fill(f64::NAN);
                    if !storage(state, &mut partial_effort)
                        || partial_effort.iter().any(|v| !v.is_finite())
                    {
                        return false;
                    }
                }
                match physical.vector_field_from_effort_into(
                    &effort,
                    &partial_effort,
                    self.bias,
                    partial_bias,
                    out,
                ) {
                    Ok(()) => true,
                    Err(error) => {
                        physical_error = Some(error);
                        false
                    }
                }
            },
            &mut cancelled,
        );
        if source_refused {
            return Err(PhsOrbitError::SourceRefused.into());
        }
        if let Some(error) = physical_error {
            return Err(error.into());
        }
        result.map_err(Into::into)
    }
}

#[cfg(test)]
mod tests;
