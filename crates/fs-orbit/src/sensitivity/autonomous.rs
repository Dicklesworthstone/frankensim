//! Frequency and full-waveform derivatives of an isolated self-excited cycle.

use super::{
    ForcedOrbitSensitivity, ForcedSensitivityControls, ForcedSensitivityDirection,
    PreparedForcedSensitivity, Result,
};
use crate::{HbOrbit, HbScales, OrbitProblem};
use fs_math::c64::C64;

/// The unchanged original autonomous balance in its fixed physical phase chart.
/// No phase row, frequency slot or normalization appears in the coefficients.
#[derive(Clone, Copy, Debug)]
pub struct AutonomousSensitivityPoint<'a> {
    pub angular_frequency: f64,
    pub harmonics: usize,
    /// Differentiate this original AFT quadrature, including its centering.
    pub aft_samples: usize,
    /// All DC states, then real and imaginary state blocks per positive harmonic.
    pub coefficients: &'a [f64],
}

/// Local physical parameter derivative with the authored phase reference fixed.
#[derive(Clone, Debug)]
pub struct AutonomousOrbitSensitivity {
    pub omega: f64,
    /// Solved frequency derivative in radians/time per physical parameter unit.
    pub angular_frequency_derivative: f64,
    /// Complete `dX_n/dp` at fixed phase, in original state/parameter units.
    /// These are Fourier half-amplitudes, including every DC and memory state.
    /// For absolute-time derivatives add the phase drift `t*omega_p*dx/dphase`.
    pub coeffs: Vec<Vec<C64>>,
    /// Componentwise backward error of the complete bordered linear solve.
    pub relative_residual: f64,
    /// Maximum two-grid physical differentiated-balance norm, combined with
    /// the normalized fixed-phase defect per declared parameter-scale increment.
    pub derivative_residual: f64,
}

impl From<ForcedOrbitSensitivity> for AutonomousOrbitSensitivity {
    fn from(value: ForcedOrbitSensitivity) -> Self {
        Self {
            omega: value.omega,
            angular_frequency_derivative: value.angular_frequency_derivative,
            coeffs: value.coeffs,
            relative_residual: value.relative_residual,
            derivative_residual: value.derivative_residual,
        }
    }
}

/// One original AFT factor with a frequency unknown and fixed phase border.
///
/// Reuses the forced sensitivity owner's original coupled Jacobian, scaling,
/// two-grid checks, resource admission, LU and residual refinement. The phase
/// row is the same full-waveform Parseval hyperplane used by autonomous HB.
/// It solves `[R_X R_omega; phase_X 0] [X_p,omega_p] = [-R_p,0]` without
/// deleting a balance row, fixing amplitude or regularizing a neutral mode.
pub struct PreparedAutonomousSensitivity {
    prepared: PreparedForcedSensitivity,
}

impl PreparedAutonomousSensitivity {
    /// Prepare at an unchanged complete source and accepted physical waveform.
    ///
    /// `phase_reference` is a complete packed autonomous seed, including omega,
    /// in the same original state chart and retained order. It remains fixed for
    /// every parameter direction. The base must satisfy both the original
    /// balances and this hyperplane; a DC-only phase reference refuses.
    ///
    /// `frequency_column` writes the COMPLETE `partial R/partial omega` with
    /// coefficients and physical parameters held fixed, in canonical real
    /// balance units per radian/time. It includes the actual frequency-port and
    /// AFT-time-chart dependence; no general port is presumed to be `s I`.
    /// Called once on each original AFT grid, it must fill every row. Refusal,
    /// nonfinite entries and unavailable exact island tangents never fall back
    /// to probes. The phase border has zero frequency derivative.
    ///
    /// This is a local fixed-order derivative. Independent order refinement is
    /// required for omitted harmonics. Singular/nonisolated cycles and excessive
    /// conditioning refuse; existence, branch uniqueness, primitive period and
    /// Floquet stability are not established. Callback work is outside the shared
    /// explicit numerical work budget, as it is for forced sensitivities.
    pub fn new<P: OrbitProblem>(
        problem: &P,
        orbit: &HbOrbit,
        phase_reference: &[f64],
        controls: ForcedSensitivityControls,
        scales: &HbScales,
        mut frequency_column: impl FnMut(AutonomousSensitivityPoint<'_>, &mut [f64]) -> bool,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        Ok(Self {
            prepared: PreparedForcedSensitivity::prepare(
                problem,
                orbit,
                controls,
                scales,
                Some(phase_reference),
                Some(&mut frequency_column),
                cancelled,
            )?,
        })
    }

    pub fn dimension(&self) -> usize {
        self.prepared.dimension()
    }
    pub fn harmonics(&self) -> usize {
        self.prepared.harmonics()
    }
    pub fn angular_frequency(&self) -> f64 {
        self.prepared.angular_frequency()
    }
    pub fn base_residual(&self) -> f64 {
        self.prepared.base_residual
    }
    pub fn tangent_sampling_defect(&self) -> f64 {
        self.prepared.tangent_sampling_defect
    }
    pub fn condition_estimate(&self) -> f64 {
        self.prepared.condition_estimate
    }

    /// Differentiate any complete source-parameter residual, with all original
    /// coefficients AND omega held fixed in this callback. The cached bordered
    /// solve supplies the unknown frequency response exactly once. Include all
    /// physical parameter-dependent port, material, source and memory terms;
    /// do not add `R_X*X_p` or `R_omega*omega_p` to this direct partial.
    ///
    /// `parameter_scale` is a positive physical increment used for conditioning
    /// and absolute residual acceptance only; output is per physical unit.
    /// The unchanged phase reference has zero direct parameter derivative.
    pub fn solve(
        &self,
        parameter_scale: f64,
        mut residual_direction: impl FnMut(AutonomousSensitivityPoint<'_>, &mut [f64]) -> bool,
        cancelled: impl FnMut() -> bool,
    ) -> Result<AutonomousOrbitSensitivity> {
        self.prepared
            .solve(
                direction(parameter_scale),
                |point, out| {
                    residual_direction(
                        AutonomousSensitivityPoint {
                            angular_frequency: point.angular_frequency,
                            harmonics: point.harmonics,
                            aft_samples: point.aft_samples,
                            coefficients: point.coefficients,
                        },
                        out,
                    )
                },
                cancelled,
            )
            .map(Into::into)
    }

    /// Differentiate a time-local island with its original linear frequency
    /// port fixed. The callback is the COMPLETE direct parameter derivative at
    /// unchanged state, not its state chain. Cached original waveforms and the
    /// same centered AFT owner project the source on both grids; no original
    /// island, Hessian, phase or frequency-column callback is repeated.
    pub fn solve_island_direction(
        &mut self,
        parameter_scale: f64,
        island_direction: impl FnMut(f64, &[f64], &mut [f64]) -> bool,
        cancelled: impl FnMut() -> bool,
    ) -> Result<AutonomousOrbitSensitivity> {
        self.prepared
            .solve_island_direction(direction(parameter_scale), island_direction, cancelled)
            .map(Into::into)
    }
}

fn direction(parameter_scale: f64) -> ForcedSensitivityDirection {
    ForcedSensitivityDirection {
        parameter: 0.0,
        parameter_scale,
        // Internally the autonomous frequency is solved as an additional state.
        angular_frequency_derivative: 0.0,
    }
}
