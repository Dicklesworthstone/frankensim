//! Implicit objective gradients through the original forced harmonic balance.

use super::{
    ForcedParameterPoint, ForcedSensitivityDirection, ForcedSensitivityError,
    PreparedForcedSensitivity, Result, SensitivityFactor, admit_direction, invalid, poll, positive,
    project_direction, scaled,
};
use fs_math::{dd::Dd, det};

/// One physical objective derivative, without a full waveform-direction solve.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ForcedObjectiveGradient {
    /// `q_p - lambda^T R_p`, in original objective/parameter units. The explicit
    /// term includes any authored clock dependence of the objective itself.
    pub derivative: f64,
    /// L2 change of the complete residual direction between the two AFT grids,
    /// normalized by equation units and the declared parameter increment.
    pub direction_sampling_defect: f64,
    /// Difference of the two contracted objective derivatives, normalized by
    /// objective units and the declared parameter increment.
    pub gradient_sampling_defect: f64,
}

/// A single objective's adjoint of a prepared forced HB response.
///
/// Solve `R_X^T lambda = q_X` once, then contract each complete physical
/// parameter direction in O(number of retained real coefficients). The existing
/// partial-pivot LU supplies its transpose solve; no transpose matrix, new
/// factor, orbit solve or parameter probe is needed. The borrow keeps the
/// weights paired with their unchanged physical base and cached AFT grids.
///
/// These are fixed-order implicit derivatives. The original source and the
/// objective derivative must remain unchanged. Harmonic refinement, objective
/// differentiability, source/clock equivalence and physical validity remain
/// caller obligations, not consequences of a small dual residual.
pub struct PreparedForcedAdjoint<'a> {
    prepared: &'a mut PreparedForcedSensitivity,
    weights: Vec<f64>,
    objective_scale: f64,
    /// Componentwise backward error of the conditioned transpose solve.
    pub relative_residual: f64,
    /// Maximum original dual-equation defect on the two independent grids,
    /// with state/objective units applied after physical contraction.
    pub dual_residual: f64,
}

impl PreparedForcedSensitivity {
    /// Prepare the derivative of one REAL objective with respect to EVERY
    /// retained real Fourier coordinate, at fixed physical forcing frequency.
    ///
    /// `objective_gradient` has the canonical packing: all DC states, followed
    /// by all real states then all imaginary states for each positive harmonic.
    /// It is the ordinary real-coordinate gradient, not complex amplitudes or
    /// a Parseval-weighted covector. For mean square of one real state,
    /// `q = X_0^2 + 2*sum(|X_n|^2)`, its entries are `2*X_0`, `4*Re X_n`, and
    /// `4*Im X_n`. Include every state, including material memory, even when its
    /// objective derivative is zero. Positive `objective_scale` conditions only
    /// the computation and its absolute acceptance; returned weights retain
    /// original objective/equation units.
    ///
    /// The existing sensitivity payload/work admission covers this O(N)
    /// additional storage and bounded transpose-solve corrections. The original
    /// maximum_relative_residual and maximum_derivative_residual govern the
    /// backward error and independent normalized dual equations, respectively.
    /// Refusal and cancellation publish no adjoint. Autonomous/conservative
    /// phase closures are outside this forced entry.
    pub fn prepare_adjoint(
        &mut self,
        objective_gradient: &[f64],
        objective_scale: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedForcedAdjoint<'_>> {
        poll(&mut cancelled)?;
        let n = self.pack.coeff_len();
        if self.pack.omega_slot.is_some()
            || self.phase.is_some()
            || self.amplitude.is_some()
            || objective_gradient.len() != n
            || objective_gradient.iter().any(|v| !v.is_finite())
            || !positive(objective_scale)
        {
            return Err(invalid(
                "forced adjoint needs a complete finite real objective gradient and positive objective scale",
            ));
        }
        let SensitivityFactor::Square(factor) = &self.factor else {
            return Err(invalid("forced adjoint requires the original square LU"));
        };
        let rhs = objective_gradient
            .iter()
            .zip(&self.scales.columns)
            .map(|(&value, &unit)| scaled(value, unit, objective_scale))
            .collect::<Result<Vec<_>>>()?;
        let mut solution = rhs.clone();
        factor.solve_transpose(&mut solution);
        let mut correction = vec![0.0; n];
        let mut relative_residual = f64::INFINITY;
        for iteration in 0..=3 {
            poll(&mut cancelled)?;
            if solution.iter().any(|v| !v.is_finite()) {
                return Err(invalid("finite forced objective adjoint required"));
            }
            relative_residual = transpose_residual(
                &self.normalized,
                &solution,
                &rhs,
                &mut correction,
                &mut cancelled,
            )?;
            if relative_residual <= self.controls.maximum_relative_residual || iteration == 3 {
                break;
            }
            factor.solve_transpose(&mut correction);
            for (value, &change) in solution.iter_mut().zip(&correction) {
                *value += change;
            }
        }
        let weights = solution
            .iter()
            .zip(&self.scales.rows)
            .map(|(&value, &unit)| scaled(value, objective_scale, unit))
            .collect::<Result<Vec<_>>>()?;
        let mut dual_residual: f64 = 0.0;
        for matrix in [&self.coarse_jacobian, &self.fine_jacobian] {
            let mut norm = 0.0;
            for col in 0..n {
                poll(&mut cancelled)?;
                let mut defect = Dd::from_f64(-objective_gradient[col]);
                for row in 0..n {
                    defect =
                        defect + Dd::from_f64(matrix[row * n + col]) * Dd::from_f64(weights[row]);
                }
                norm = det::hypot(
                    norm,
                    scaled(defect.to_f64(), self.scales.columns[col], objective_scale)?,
                );
            }
            dual_residual = dual_residual.max(norm);
        }
        if !relative_residual.is_finite()
            || !dual_residual.is_finite()
            || relative_residual > self.controls.maximum_relative_residual
            || dual_residual > self.controls.maximum_derivative_residual
        {
            return Err(ForcedSensitivityError::DerivativeUnresolved {
                relative_residual,
                derivative_residual: dual_residual,
            });
        }
        poll(&mut cancelled)?;
        Ok(PreparedForcedAdjoint {
            prepared: self,
            weights,
            objective_scale,
            relative_residual,
            dual_residual,
        })
    }
}

impl PreparedForcedAdjoint<'_> {
    /// Original objective/equation weights, in canonical real residual order.
    /// No extra complex-conjugation or harmonic multiplicity applies when
    /// contracting these weights with a packed real residual direction.
    pub fn residual_weights(&self) -> &[f64] {
        &self.weights
    }

    /// Contract a COMPLETE physical residual direction with the cached adjoint.
    /// The callback uses the same convention as PreparedForcedSensitivity::solve:
    /// differentiate original forcing, material, frequency port and clock at
    /// FIXED Fourier coordinates. `explicit_derivative` is the objective's own
    /// partial at fixed coordinates, including any explicit frequency term.
    /// Each of the two AFT grids is evaluated once; neither invokes a linear
    /// solve. Direction and contracted-gradient sampling changes must both meet
    /// the prepared maximum_derivative_residual in their normalized units.
    pub fn gradient(
        &self,
        direction: ForcedSensitivityDirection,
        explicit_derivative: f64,
        mut residual_direction: impl FnMut(ForcedParameterPoint<'_>, &mut [f64]) -> bool,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ForcedObjectiveGradient> {
        admit_objective_direction(direction, explicit_derivative, &mut cancelled)?;
        let mut columns = Vec::with_capacity(2);
        for samples in [self.prepared.coarse.samples, self.prepared.fine.samples] {
            poll(&mut cancelled)?;
            let mut column = vec![f64::NAN; self.weights.len()];
            let point = ForcedParameterPoint {
                parameter: direction.parameter,
                angular_frequency: self.prepared.omega,
                harmonics: self.prepared.pack.n_harm,
                aft_samples: samples,
                coefficients: &self.prepared.coefficients,
            };
            if !residual_direction(point, &mut column) || column.iter().any(|v| !v.is_finite()) {
                return Err(ForcedSensitivityError::DirectionRefused { samples });
            }
            columns.push(column);
        }
        self.contract(
            direction,
            explicit_derivative,
            &columns[0],
            &columns[1],
            &mut cancelled,
        )
    }

    /// Project a complete local-island parameter partial using the SAME cached
    /// waveforms and centered AFT owner as forward sensitivities. Clock and
    /// linear port must be fixed; other directions use gradient(). There are
    /// no baseline source/tangent calls and no per-parameter linear solves.
    pub fn gradient_island_direction(
        &mut self,
        direction: ForcedSensitivityDirection,
        explicit_derivative: f64,
        mut island_direction: impl FnMut(f64, &[f64], &mut [f64]) -> bool,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ForcedObjectiveGradient> {
        admit_objective_direction(direction, explicit_derivative, &mut cancelled)?;
        if direction.angular_frequency_derivative != 0.0 {
            return Err(invalid(
                "island-only objective direction requires a fixed physical clock",
            ));
        }
        let p = &mut self.prepared;
        let coarse = project_direction(
            &p.pack,
            p.omega,
            &mut p.coarse,
            &mut island_direction,
            &mut cancelled,
        )?;
        let fine = project_direction(
            &p.pack,
            p.omega,
            &mut p.fine,
            &mut island_direction,
            &mut cancelled,
        )?;
        self.contract(
            direction,
            explicit_derivative,
            &coarse,
            &fine,
            &mut cancelled,
        )
    }

    fn contract(
        &self,
        direction: ForcedSensitivityDirection,
        explicit: f64,
        coarse: &[f64],
        fine: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ForcedObjectiveGradient> {
        let mut values = [Dd::from_f64(explicit); 2];
        let mut direction_sampling_defect = 0.0;
        for row in 0..self.weights.len() {
            poll(cancelled)?;
            let weight = Dd::from_f64(self.weights[row]);
            values[0] = values[0] - weight * Dd::from_f64(coarse[row]);
            values[1] = values[1] - weight * Dd::from_f64(fine[row]);
            let difference = Dd::from_f64(coarse[row]) - Dd::from_f64(fine[row]);
            direction_sampling_defect = det::hypot(
                direction_sampling_defect,
                scaled(
                    difference.to_f64(),
                    direction.parameter_scale,
                    self.prepared.scales.rows[row],
                )?,
            );
        }
        let gradient_sampling_defect = scaled(
            (values[0] - values[1]).to_f64(),
            direction.parameter_scale,
            self.objective_scale,
        )?
        .abs();
        if !direction_sampling_defect.is_finite()
            || direction_sampling_defect > self.prepared.controls.maximum_derivative_residual
            || gradient_sampling_defect > self.prepared.controls.maximum_derivative_residual
        {
            return Err(ForcedSensitivityError::DerivativeUnresolved {
                relative_residual: self.relative_residual,
                derivative_residual: direction_sampling_defect.max(gradient_sampling_defect),
            });
        }
        let derivative = values[1].to_f64();
        if !derivative.is_finite() {
            return Err(invalid("physical objective gradient must be representable"));
        }
        poll(cancelled)?;
        Ok(ForcedObjectiveGradient {
            derivative,
            direction_sampling_defect,
            gradient_sampling_defect,
        })
    }
}

fn admit_objective_direction(
    direction: ForcedSensitivityDirection,
    explicit: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<()> {
    poll(cancelled)?;
    admit_direction(direction)?;
    if !explicit.is_finite() {
        return Err(invalid("explicit objective derivative must be finite"));
    }
    Ok(())
}

// Contract the transpose in place: no second dense matrix or factor is stored.
fn transpose_residual(
    matrix: &[f64],
    solution: &[f64],
    rhs: &[f64],
    correction: &mut [f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64> {
    let n = solution.len();
    let mut maximum: f64 = 0.0;
    for col in 0..n {
        poll(cancelled)?;
        let mut residual = Dd::from_f64(-rhs[col]);
        let mut magnitude = Dd::from_f64(rhs[col].abs());
        for row in 0..n {
            let term = Dd::from_f64(matrix[row * n + col]) * Dd::from_f64(solution[row]);
            residual = residual + term;
            magnitude = magnitude + term.abs();
        }
        correction[col] = -residual.to_f64();
        let relative = if magnitude == Dd::ZERO {
            0.0
        } else {
            (residual / magnitude).to_f64().abs()
        };
        if !relative.is_finite() || !correction[col].is_finite() {
            return Err(invalid("forced adjoint transpose residual must be finite"));
        }
        maximum = maximum.max(relative);
    }
    Ok(maximum)
}
