//! Complete parameter derivatives of the original periodically pumped system.

use super::{
    Adapter, MAX_PARAMETRIC_HARMONICS, ParametricError, ParametricHamiltonian,
    PeriodicQuadraticModulation, Result,
};
use crate::phs_orbit::{PhsForcedParameterDirection, sensitivity::ode_clock_column};
use fs_orbit::sensitivity::{
    ForcedOrbitSensitivity, ForcedSensitivityControls, ForcedSensitivityDirection,
    ForcedSensitivityError, PreparedForcedSensitivity,
};
use fs_orbit::{HbOrbit, HbScales};

mod direction;
mod power;
pub use power::{ParametricPowerObjective, PreparedParametricPowerAdjoint};

/// One complete physical parameter direction at unchanged state and phase.
/// The same clock drives the original input and Hamiltonian modulation.
pub struct ParametricParameterDirection<'a> {
    /// Original J/R/G, H0 and input Fourier-coefficient derivatives. The storage
    /// callback writes partial_p grad H0; the modulation contribution is added
    /// separately, exactly once, from `modulation` below.
    pub physical: PhsForcedParameterDirection<'a>,
    /// Derivative of the complete original symmetric Fourier energy matrices,
    /// at FIXED PHASE. None holds Q fixed. Dimension, clock and retained
    /// bandwidth must match the prepared source; signed directions are allowed.
    pub modulation: Option<&'a PeriodicQuadraticModulation>,
}

/// Original time-dependent Hamiltonian tangent/factor and cached full-state
/// trajectories. The source, pump, drive and any authored H0 Hessian stay bound
/// unchanged for every waveform and objective parameter query.
pub struct PreparedParametricSensitivity<'s, 'a> {
    source: &'s ParametricHamiltonian<'a>,
    prepared: PreparedForcedSensitivity,
}

impl From<ForcedSensitivityError> for ParametricError {
    fn from(value: ForcedSensitivityError) -> Self {
        Self::Sensitivity(value)
    }
}
impl From<fs_phs::PhsError> for ParametricError {
    fn from(value: fs_phs::PhsError) -> Self {
        Self::Physical(value)
    }
}

impl<'a> ParametricHamiltonian<'a> {
    /// Prepare exact implicit derivatives of the ORIGINAL pumped HB response.
    /// The full tangent is `(J-R)*(H0_xx+Q(theta))`, including every nonlinear
    /// mode and explicit memory state. The existing fs-orbit owner reevaluates
    /// the base and both independent AFT Jacobians before factoring once.
    ///
    /// The complete original H0 tangent is required; use `new_analytic` when
    /// its author supplies an explicit action. No quadratic replacement,
    /// finite-difference fallback or new Fourier/Newton/factor owner is used.
    ///
    /// Query scratch for the largest admitted modulation direction is reserved
    /// before numerical preparation. Already-owned source storage and callback
    /// work retain the generic sensitivity contract's exclusions; source
    /// dimensions/bandwidth remain subject to this module's hard limits.
    /// Nonlinear internal resistance is currently refused: this parameter
    /// direction contract does not yet supply the independent original D_p.
    /// Omitting that chain would produce an incomplete physical derivative.
    pub fn prepare_sensitivity(
        &self,
        orbit: &HbOrbit,
        mut controls: ForcedSensitivityControls,
        scales: &HbScales,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedParametricSensitivity<'_, 'a>> {
        if cancelled() {
            return Err(ForcedSensitivityError::Cancelled.into());
        }
        if self.dissipation.is_some() {
            return Err(ParametricError::Invalid(
                "nonlinear internal resistance requires a complete original dissipation parameter direction",
            ));
        }
        let harmonics = orbit
            .coeffs
            .first()
            .map_or(0, |row| row.len().saturating_sub(1));
        self.admit_order(harmonics)?;
        if orbit.omega != self.angular_frequency() || !self.tangent.available() {
            return Err(ParametricError::Invalid(
                "parametric sensitivity requires the response clock and complete H0 tangent",
            ));
        }
        let (n, m) = (
            self.system.state_dim() as u128,
            self.system.port_dim() as u128,
        );
        let bytes = 8 * (2 * n * n + 3 * n + m)
            + 16 * ((harmonics / self.periods).min(MAX_PARAMETRIC_HARMONICS) as u128 + 1);
        if bytes > controls.maximum_payload_bytes as u128 {
            return Err(ForcedSensitivityError::Budget {
                quantity: "parametric physical direction payload",
                requested: bytes,
                cap: controls.maximum_payload_bytes as u128,
            }
            .into());
        }
        controls.maximum_payload_bytes -= bytes as usize;
        self.source_refused.set(false);
        let prepared =
            PreparedForcedSensitivity::new(&Adapter(self), orbit, controls, scales, &mut cancelled);
        if matches!(prepared, Err(ForcedSensitivityError::Cancelled)) {
            return Err(ForcedSensitivityError::Cancelled.into());
        }
        if self.source_refused.get() {
            return Err(ParametricError::SourceRefused);
        }
        let prepared = prepared?;
        if cancelled() {
            return Err(ForcedSensitivityError::Cancelled.into());
        }
        Ok(PreparedParametricSensitivity {
            source: self,
            prepared,
        })
    }
}

impl PreparedParametricSensitivity<'_, '_> {
    /// Original factor, complete-state scales and independent base diagnostics.
    pub fn prepared(&self) -> &PreparedForcedSensitivity {
        &self.prepared
    }

    /// Differentiate all original Fourier states under simultaneous material,
    /// matrix, pump, drive and prescribed-clock changes. At fixed phase,
    /// `ep=partial_p grad H0 + Qp*x` and the original pHS owner supplies
    /// `(J-R)ep+(Jp-Rp)e+Gp*u+G*up`. The original sI clock column supplies
    /// `i*n*Omega_p*X_n`; Q and the drive gain no extra fixed-time phase term.
    /// Omega_p is the RESPONSE-chart clock derivative. On an over_periods(m)
    /// chart the original physical source clock therefore changes by m*Omega_p;
    /// authored pump/input coefficient directions retain their SOURCE clock.
    ///
    /// These are fixed-order implicit derivatives, not a Floquet stability or
    /// omitted-harmonic certificate. Source refusal/cancellation publishes no
    /// partial waveform and leaves the original factor reusable.
    pub fn solve_parameter(
        &mut self,
        chart: ForcedSensitivityDirection,
        direction: ParametricParameterDirection<'_>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ForcedOrbitSensitivity> {
        if cancelled() {
            return Err(ForcedSensitivityError::Cancelled.into());
        }
        let n = self.source.system.state_dim();
        direction::with_samples(
            self.source,
            self.prepared.harmonics(),
            direction,
            |samples| {
                self.prepared.solve_split_direction(
                    chart,
                    |time, state, out| samples.borrow_mut().field(time, state, out),
                    |point, out| {
                        ode_clock_column(n, chart.angular_frequency_derivative, point, out)
                    },
                    &mut cancelled,
                )
            },
        )
    }
}
