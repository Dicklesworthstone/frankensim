//! Frequency sensitivity of the changing original forced periodic response.

use super::super::{
    ObservationBudget, PhsPortSidebandCheck, minimum_observation_work, observe_with_budget,
};
use super::{
    C64, FactoredPhsPortSidebands, OrbitProblem, PhsPortSidebandError, PreparedPhsPortSidebands,
    Result, SidebandResponse, poll,
};

/// Local response change for a unit relative change of the forcing frequency.
///
/// These coefficients are `omega * dX/domega` and `omega * dY/domega`, at fixed
/// forcing phase and fixed authored DC/complex input coefficients. They retain
/// the original state/output units; divide by `omega_rad_s` for derivatives per
/// rad/s. They are derivatives of the changing nonlinear BASE response, unlike
/// a fixed-base small-probe carrier derivative.
///
/// Both signed Fourier indices participate, with zero carrier. A real waveform
/// is `Re(sum_k Z_k exp(i*k*theta))`: the positive and negative entries are the
/// conjugate Fourier coefficients, without an additional factor of two. DC can
/// change under nonlinear mixing even though the clock derivative's DC forcing
/// vanishes. Phase-zero states and every explicit material-memory coordinate
/// remain in their original order. The derivative is at fixed `theta=omega*t`,
/// not at fixed physical time; the latter also needs `t * dx/dt`.
#[derive(Clone, Debug)]
pub struct PhsForcedFrequencyTangent {
    pub omega_rad_s: f64,
    pub state_harmonics: usize,
    pub output_harmonics: usize,
    pub state_dimension: usize,
    pub port_count: usize,
    /// Sideband-major, then original state: `(k+state_harmonics)*n+j`.
    /// Residuals check the differentiated original equations and omitted bands.
    pub state_log_frequency_derivative: SidebandResponse,
    /// `(k+output_harmonics)*ports+p`, in original physical output units.
    pub port_log_frequency_derivative: Vec<C64>,
    /// Original nonlinear base balance, independently rechecked at preparation.
    pub base_residual: f64,
    pub tangent_sampling_defect: f64,
    /// Independent physical-output checks in the original output units.
    pub observation_trace: Vec<PhsPortSidebandCheck>,
    pub hessian_actions: usize,
}

impl PreparedPhsPortSidebands<'_, '_> {
    /// Factor at zero carrier and differentiate the complete forced base orbit.
    /// See [`FactoredPhsPortSidebands::forced_frequency_tangent`] for its fixed
    /// source/phase contract, full-state equation and resolution requirements.
    ///
    /// # Errors
    /// Autonomous sources refuse before factorization. All original source,
    /// singularity, residual, omitted-band, observation and cancellation errors
    /// remain active; no partial derivative is returned.
    pub fn forced_frequency_tangent(
        &mut self,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PhsForcedFrequencyTangent> {
        poll(&mut cancelled)?;
        if self.source.problem.autonomous() {
            return Err(PhsPortSidebandError::Invalid(
                "a forced frequency tangent requires an externally fixed base clock",
            ));
        }
        self.factor(0.0, &mut cancelled)?
            .forced_frequency_tangent(cancelled)
    }
}

impl FactoredPhsPortSidebands<'_, '_> {
    /// Differentiate a forced nonlinear branch using this SAME zero-carrier LU.
    ///
    /// In the original complete pHS ODE, write `omega*x_theta=f(x)+G*u(theta)`.
    /// Holding the physical model and authored input phasors fixed gives
    /// `(omega*d/dtheta-Df) z = -omega*x_theta`, with `z=d x/d ln(omega)`.
    /// The right-hand side is `-i*k*omega*X_k` at every signed base harmonic.
    /// fs-orbit already owns this exact periodic tangent, solve, residual
    /// correction and omitted-equation checks. No second Jacobian, perturbed
    /// frequency, fitted resonance or mass-normalized replacement is introduced.
    ///
    /// The output derivative is the original `G^T Hxx(x_base) z`, observed on
    /// independently refined grids. Native nonlinear storage, explicit memory,
    /// and driven instantaneous feedback retain their complete derivatives.
    /// Fixed phasors mean relative phase, amplitude and DC input stay fixed:
    /// physically delayed/dispersive forcing must be realized in the original
    /// state, rather than silently treated as a frequency-independent input.
    ///
    /// All active base harmonics must fit in the retained response order.
    /// Preparation rechecks base balance and tangent sampling; this operation
    /// keeps the full retained/omitted residual gates and physical output limits.
    /// Log-frequency normalization leaves these tolerances in their existing
    /// state/equation/output units. Derivative resolution still needs harmonic
    /// refinement; a small base residual alone does not resolve the tangent.
    ///
    /// This is a local finite-harmonic sensitivity, not a branch selector,
    /// Floquet classification, forward-error bound or continuum certificate.
    /// At a fold a frequency graph may cease to exist: the singular native
    /// operator refuses without regularization or switching branches. Near a
    /// fold a large derivative and small backward error do not bound sensitivity.
    /// Work/payload use the source's original per-operation caps; the clock RHS
    /// is charged before the existing physical observation. Sources remain fixed.
    ///
    /// # Errors
    /// Nonzero carrier, autonomous source, omitted active base coefficients,
    /// unrepresentable arithmetic, budget/resolution/source refusal or
    /// cancellation. No incomplete derivative is published.
    pub fn forced_frequency_tangent(
        &mut self,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PhsForcedFrequencyTangent> {
        poll(&mut cancelled)?;
        if self.carrier_angular_frequency() != 0.0 || self.source.problem.autonomous() {
            return Err(PhsPortSidebandError::Invalid(
                "forced frequency tangents require a zero-carrier factor and an external base clock",
            ));
        }
        let original = self.source.problem.original();
        let (n, m) = (original.system.state_dim(), original.system.port_dim());
        let controls = &self.source.controls;
        let h = controls.sideband.harmonics;
        let base = self.source.orbit;
        if base
            .coeffs
            .iter()
            .any(|row| row.iter().skip(h + 1).any(|z| z.re != 0.0 || z.im != 0.0))
        {
            return Err(PhsPortSidebandError::Invalid(
                "the frequency tangent must retain every active original base harmonic",
            ));
        }
        let mut budget = ObservationBudget {
            actions: 0,
            work: 16 * n as u128 * (2 * h + 1) as u128,
            retained_bytes: 0,
        };
        if budget.work + minimum_observation_work(n, m, controls.initial_samples)
            > controls.max_observation_work as u128
        {
            return Err(PhsPortSidebandError::Budget(
                "frequency tangent and physical observation exceed the arithmetic cap",
            ));
        }
        let state_log_frequency_derivative = {
            // Source/response admission bounds these extents. The existing
            // adapter payload includes both this RHS and the returned response.
            let mut rhs = vec![C64::ZERO; n * (2 * h + 1)];
            for (j, row) in base.coeffs.iter().enumerate() {
                poll(&mut cancelled)?;
                for (k, &x) in row.iter().enumerate().skip(1).take(h) {
                    let rate = k as f64 * base.omega;
                    let value = C64::new(rate * x.im, -rate * x.re);
                    if !value.re.is_finite() || !value.im.is_finite() {
                        return Err(PhsPortSidebandError::Invalid(
                            "original clock derivative is not representable",
                        ));
                    }
                    rhs[(h + k) * n + j] = value;
                    rhs[(h - k) * n + j] = value.conj();
                }
            }
            self.factor.solve(&rhs, &mut cancelled)?
        };
        let (port_log_frequency_derivative, observation_trace, hessian_actions) =
            observe_with_budget(
                original,
                base,
                &state_log_frequency_derivative.coefficients,
                controls,
                &controls.envelope_tolerances,
                &mut budget,
                &mut cancelled,
            )?;
        poll(&mut cancelled)?;
        Ok(PhsForcedFrequencyTangent {
            omega_rad_s: base.omega,
            state_harmonics: h,
            output_harmonics: controls.output_harmonics,
            state_dimension: n,
            port_count: m,
            state_log_frequency_derivative,
            port_log_frequency_derivative,
            base_residual: self.base_residual,
            tangent_sampling_defect: self.tangent_sampling_defect,
            observation_trace,
            hessian_actions,
        })
    }
}
