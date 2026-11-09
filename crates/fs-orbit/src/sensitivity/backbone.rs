//! Fixed-amplitude derivatives of a conservative periodic family.

use super::{
    AutonomousOrbitSensitivity, AutonomousSensitivityPoint, ForcedSensitivityControls,
    ForcedSensitivityDirection, ForcedSensitivityError, Lu, PreparedForcedSensitivity, Result,
    invalid, lu, positive, scaled,
};
use crate::{
    BackboneContinuationError, HbOrbit, HbScales, OrbitError, OrbitProblem, Packing,
    backbone::continuation::linear, scaling::AdmittedScales,
};
use fs_math::{dd::Dd, det};

/// A physical chart on the conservative family, with no state reordering.
#[derive(Clone, Copy, Debug)]
pub struct BackboneSensitivityAnchor {
    /// Original state whose first-harmonic peak fixes the amplitude.
    pub coordinate: usize,
    /// Positive `2*abs(X_1[coordinate])` in the original coordinate's units.
    /// This is neither the turning-point amplitude nor the full waveform peak.
    pub first_harmonic_amplitude: f64,
}

/// An original parameter direction with an independently authored amplitude law.
#[derive(Clone, Copy, Debug)]
pub struct BackboneSensitivityDirection {
    /// Positive parameter increment used for numerical scaling only.
    pub parameter_scale: f64,
    /// Physical `dA/dp`. Use zero for material/geometry changes at fixed A,
    /// or one with a zero source direction to differentiate with respect to A.
    pub first_harmonic_amplitude_derivative: f64,
}
impl BackboneSensitivityDirection {
    fn original(self) -> Result<ForcedSensitivityDirection> {
        if !positive(self.parameter_scale) || !self.first_harmonic_amplitude_derivative.is_finite()
        {
            return Err(invalid(
                "positive parameter scale and finite amplitude derivative required",
            ));
        }
        Ok(ForcedSensitivityDirection {
            parameter: 0.0,
            parameter_scale: self.parameter_scale,
            angular_frequency_derivative: 0.0,
        })
    }
}

/// Frequency and complete waveform derivatives along a conservative family.
///
/// All original HB rows remain present. A fixed full-waveform phase row and
/// a physical amplitude row close the unknown coefficients and frequency:
/// `[R_X R_omega; phase_X 0; A_X 0] [X_p,omega_p] = [-R_p,0,A_p]`.
/// The original backbone-continuation Householder QR factors this rectangular
/// system once. No normal equations, removed balance row, artificial damping,
/// singular-cycle LU or finite-difference tangent is introduced.
///
/// Conservation leaves a redundant balance relation; consistency is therefore
/// checked, not presumed. In particular a new dissipative direction generally
/// admits no nonzero fixed-amplitude conservative orbit and must refuse.
pub struct PreparedBackboneSensitivity {
    prepared: PreparedForcedSensitivity,
}
impl PreparedBackboneSensitivity {
    /// Prepare the exact original AFT tangent and frequency column on two grids.
    /// The source must be autonomous and unchanged during subsequent queries.
    /// `phase_reference` is a complete packed autonomous seed in the same state
    /// chart/order, including frequency. The base must satisfy that fixed phase
    /// and the authored amplitude as well as every original balance row.
    ///
    /// `frequency_column` declares the complete partial R/partial omega at fixed
    /// physical parameters and coefficients. It includes the actual port and AFT
    /// time chart; no `s I` port is inferred. It must fill every original row.
    /// Controls retain the shared payload/work, independent residual and tangent
    /// sampling bounds. The condition estimate concerns the column-equilibrated
    /// QR R factor and is not a certified rank or forward-error bound.
    /// Relative solve error uses maximum row defect divided by maximum row
    /// action; the absolute two-grid physical residual gate remains independent.
    ///
    /// These local fixed-order derivatives neither prove conservation of an
    /// arbitrary source nor establish existence, harmonic-tail convergence,
    /// primitive period or Floquet stability. An amplitude fold, additional
    /// neutral family or incompatible parameter direction can refuse.
    #[allow(clippy::too_many_arguments)]
    pub fn new<P: OrbitProblem>(
        problem: &P,
        orbit: &HbOrbit,
        anchor: BackboneSensitivityAnchor,
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
                Some(anchor),
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

    /// Differentiate a complete original residual partial at fixed coefficients
    /// and omega. Include original material, port, geometry and source changes;
    /// omit state/frequency chain terms, which the factor solves for once.
    /// The phase reference stays fixed. Output includes EVERY original DC and
    /// harmonic coefficient in state/parameter units, plus d omega/dp.
    pub fn solve(
        &self,
        direction: BackboneSensitivityDirection,
        mut residual_direction: impl FnMut(AutonomousSensitivityPoint<'_>, &mut [f64]) -> bool,
        cancelled: impl FnMut() -> bool,
    ) -> Result<AutonomousOrbitSensitivity> {
        self.prepared
            .solve_with_amplitude(
                direction.original()?,
                direction.first_harmonic_amplitude_derivative,
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

    /// Differentiate the original local island with its frequency port fixed.
    /// Cached full physical waveforms and the same centered AFT projection form
    /// the two independent parameter columns. No nonlinear orbit solve or source
    /// state/frequency-tangent callback repeats. A zero callback with dA/dp=1
    /// gives the backbone slope and all harmonic changes with amplitude.
    pub fn solve_island_direction(
        &mut self,
        direction: BackboneSensitivityDirection,
        island_direction: impl FnMut(f64, &[f64], &mut [f64]) -> bool,
        cancelled: impl FnMut() -> bool,
    ) -> Result<AutonomousOrbitSensitivity> {
        self.prepared
            .solve_island_with_amplitude(
                direction.original()?,
                direction.first_harmonic_amplitude_derivative,
                island_direction,
                cancelled,
            )
            .map(Into::into)
    }
}

pub(super) struct AmplitudeClosure {
    indices: [usize; 2],
    row: [f64; 2],
    amplitude: f64,
    pub(super) base_defect: f64,
}
impl AmplitudeClosure {
    pub(super) fn new(
        anchor: BackboneSensitivityAnchor,
        pack: &Packing,
        coefficients: &[f64],
        scales: &AdmittedScales,
    ) -> Result<Self> {
        if anchor.coordinate >= pack.d || !positive(0.5 * anchor.first_harmonic_amplitude) {
            return Err(invalid(
                "original coordinate and positive representable amplitude required",
            ));
        }
        let indices = [
            pack.re_index(1, anchor.coordinate),
            pack.im_index(1, anchor.coordinate),
        ];
        let [x, y] = indices.map(|i| coefficients[i]);
        let radius = det::hypot(x, y);
        if !positive(radius) {
            return Err(invalid("nonzero original fundamental required"));
        }
        let half = 0.5 * anchor.first_harmonic_amplitude;
        let row = [
            scaled(x / radius, scales.columns[indices[0]], half)?,
            scaled(y / radius, scales.columns[indices[1]], half)?,
        ];
        Ok(Self {
            indices,
            row,
            amplitude: anchor.first_harmonic_amplitude,
            base_defect: radius / half - 1.0,
        })
    }
    pub(super) fn write_row(&self, out: &mut [f64]) {
        for k in 0..2 {
            out[self.indices[k]] = self.row[k];
        }
    }
    pub(super) fn rhs(&self, derivative: f64, parameter_scale: f64) -> Result<f64> {
        scaled(derivative, parameter_scale, self.amplitude)
    }
    pub(super) fn derivative_defect(
        &self,
        physical: &[f64],
        columns: &[f64],
        derivative: f64,
        parameter_scale: f64,
    ) -> Result<f64> {
        let mut residual = Dd::from_f64(-self.rhs(derivative, parameter_scale)?);
        for k in 0..2 {
            let col = self.indices[k];
            residual = residual
                + Dd::from_f64(self.row[k])
                    * Dd::from_f64(scaled(physical[col], parameter_scale, columns[col])?);
        }
        Ok(residual.to_f64())
    }
}

pub(super) enum SensitivityFactor {
    Square(Lu),
    Rectangular(linear::Admitted),
}
impl SensitivityFactor {
    pub(super) fn new(
        matrix: &[f64],
        n: usize,
        maximum_condition: f64,
        rectangular: bool,
    ) -> Result<(Self, f64)> {
        if rectangular {
            let factor = linear::admit(&mut matrix.to_vec(), n + 1, n, maximum_condition)
                .map_err(qr_error)?;
            let condition = factor.condition;
            Ok((Self::Rectangular(factor), condition))
        } else {
            let factor = lu(matrix, n).map_err(|_| ForcedSensitivityError::Singular)?;
            let condition = factor.condition_1(matrix);
            Ok((Self::Square(factor), condition))
        }
    }
    pub(super) fn solve(&self, rhs: &[f64]) -> Result<Vec<f64>> {
        match self {
            Self::Square(factor) => {
                let mut solution = rhs.to_vec();
                factor.solve(&mut solution);
                Ok(solution)
            }
            Self::Rectangular(factor) => factor.solve(rhs).map_err(qr_error),
        }
    }
}
fn qr_error(error: BackboneContinuationError) -> ForcedSensitivityError {
    match error {
        BackboneContinuationError::IllConditioned { estimate, maximum } => {
            ForcedSensitivityError::IllConditioned {
                estimate,
                limit: maximum,
            }
        }
        BackboneContinuationError::Numerical(OrbitError::SingularJacobian) => {
            ForcedSensitivityError::Singular
        }
        BackboneContinuationError::Numerical(error) => ForcedSensitivityError::Orbit(error),
        _ => invalid("original rectangular sensitivity QR refused"),
    }
}
