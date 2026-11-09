//! Prepare the full passive Kelvin/Maxwell operator without forming its Gram.

use super::*;
use fs_la::factor::{qr, svd_jacobi};
use fs_math::det;

impl PreparedTensorKelvinStressControl {
    /// Bind this history's immutable material and factor its complete coupled
    /// constitutive operator once. Preparation does not advance the history.
    ///
    /// # Errors
    /// Missing or zero Kelvin viscosity, exceeded dimension/storage/work
    /// budgets, unresolved original factors, unrepresentable rates, cancellation.
    #[allow(clippy::too_many_lines)] // One bounded preparation of the original full operator.
    pub fn new(
        history: &TensorPronyHistory,
        budget: TensorKelvinStressControlBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(&mut cancelled)?;
        let model = &history.model;
        let n = 6 + 6 * model.branches().len();
        let storage = n
            .checked_mul(n)
            .and_then(|x| x.checked_mul(16))
            .and_then(|x| x.checked_add(128 * n + 1024));
        let work = n
            .checked_pow(3)
            .and_then(|x| x.checked_mul(1024))
            .and_then(|x| x.checked_add(100 * n * n + 10_000));
        if n > budget.max_state_coordinates
            || storage.is_none_or(|s| s > budget.max_dense_scalars)
            || work.is_none_or(|w| {
                w > budget.max_preparation_work || w > fs_la::eigen::MAX_EIGEN_UNPOLLED_SCALAR_WORK
            })
        {
            return Err(invalid(
                "Kelvin stress-control preparation exceeds its capacity or work budget",
            ));
        }
        let kelvin_factor = model
            .kv_factor
            .as_ref()
            .ok_or_else(|| invalid("stress control requires a nonzero Kelvin viscosity tensor"))?;
        // The original PSD owner leaves zero columns for missing dissipative
        // directions. Rotation preserves those columns; do not promote a
        // numerically tiny QR pivot of a rank-deficient tensor to viscosity.
        let columns: Vec<_> = (0..6)
            .filter(|&j| kelvin_factor.iter().any(|row| row[j] != 0.0))
            .collect();
        if columns.is_empty() {
            return Err(invalid(
                "zero Kelvin viscosity uses the Maxwell stress controller",
            ));
        }
        if columns.len() < 6 {
            return Self::prepare_semidefinite(history, &columns, &mut cancelled);
        }
        let viscosity_scale = kelvin_factor
            .iter()
            .flatten()
            .fold(0.0_f64, |s, v| s.max(v.abs()));
        if !viscosity_scale.is_finite() || viscosity_scale <= 0.0 {
            return Err(invalid("unrepresentable Kelvin factor scale"));
        }
        let mut rows = [0.0; 36];
        for i in 0..6 {
            for j in 0..6 {
                rows[6 * i + j] = kelvin_factor[j][i] / viscosity_scale;
            }
        }
        let viscosity_qr = qr(&rows, 6, 6);
        poll(&mut cancelled)?;
        if (0..6).any(|i| !viscosity_qr.r(i, i).is_finite() || viscosity_qr.r(i, i) == 0.0) {
            return Err(invalid("Kelvin viscosity factor is numerically singular"));
        }
        let mut prepared = Self {
            model: Arc::clone(model),
            viscosity_qr: Some(viscosity_qr),
            viscosity_scale,
            constrained: None,
            branch_rate_roots: Vec::new(),
            rates: Vec::new(),
            left_modes: Vec::new(),
            strain_modes: Vec::new(),
            memory_modes: Vec::new(),
        };
        let mut factor = zeros(n * n)?;
        for column in 0..6 {
            let transformed =
                prepared.solve_viscosity_root(model.equilibrium_gram[column], true)?;
            for row in 0..6 {
                factor[row * n + column] = transformed[row];
            }
        }
        for (j, branch) in model.branches().iter().enumerate() {
            poll(&mut cancelled)?;
            let root = checked(
                1.0 / det::sqrt(branch.tau_s),
                "unrepresentable Maxwell rate root",
            )?;
            if root <= 0.0 {
                return Err(invalid("Maxwell rate root underflowed"));
            }
            prepared.branch_rate_roots.push(root);
            let b = branch.resolved_gram_rows(&model.factors[j]);
            for (axis, row) in b.iter().enumerate() {
                let column = 6 + 6 * j + axis;
                let transformed = prepared.solve_viscosity_root(*row, true)?;
                for row in 0..6 {
                    factor[row * n + column] = transformed[row];
                }
                factor[column * n + column] = root;
            }
        }
        prepared.prepare_modes(factor, n, None, &mut cancelled)?;
        Ok(prepared)
    }

    // The same native loss factor and exact scalar ramps serve both viscosity
    // ranks. `energy_basis` embeds constrained dynamic coordinates into the
    // complete original equilibrium/Maxwell state; definite viscosity is the
    // identity case and keeps its established arithmetic.
    #[allow(clippy::too_many_lines)] // One native factor preparation and original-source checks.
    pub(super) fn prepare_modes(
        &mut self,
        mut factor: Vec<f64>,
        n: usize,
        energy_basis: Option<&[f64]>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<()> {
        let scale = factor.iter().fold(0.0_f64, |s, v| s.max(v.abs()));
        if !scale.is_finite() || scale <= 0.0 {
            return Err(invalid("unrepresentable Kelvin relaxation scale"));
        }
        for value in &mut factor {
            *value /= scale;
        }
        poll(cancelled)?;
        let decomposition = svd_jacobi(&factor, n, n);
        poll(cancelled)?;
        let tolerance = 8192.0 * f64::EPSILON * n as f64;
        let memory = 6 * self.model.branches().len();
        self.memory_modes = zeros(memory * n)?;
        for mode in 0..n {
            poll(cancelled)?;
            let singular = decomposition.sigma[mode];
            if !singular.is_finite() || singular <= 0.0 {
                return Err(invalid("Kelvin relaxation singular value is unresolved"));
            }
            for row in 0..n {
                let mut actual = Dd::ZERO;
                for col in 0..n {
                    actual = actual
                        + Dd::from_f64(factor[row * n + col])
                            * Dd::from_f64(decomposition.v[col * n + mode]);
                }
                let residual = (actual
                    - Dd::from_f64(singular) * Dd::from_f64(decomposition.u[row * n + mode]))
                .to_f64();
                if !residual.is_finite() || (residual / singular).abs() > tolerance {
                    return Err(invalid("Kelvin relaxation factor residual is unresolved"));
                }
            }
            for other in 0..n {
                for vectors in [&decomposition.u, &decomposition.v] {
                    let mut product = Dd::ZERO;
                    for row in 0..n {
                        product = product
                            + Dd::from_f64(vectors[row * n + mode])
                                * Dd::from_f64(vectors[row * n + other]);
                    }
                    let residual = product.to_f64() - if mode == other { 1.0 } else { 0.0 };
                    if !residual.is_finite() || residual.abs() > tolerance {
                        return Err(invalid("Kelvin relaxation vectors are unresolved"));
                    }
                }
            }
            let rate_root = checked(scale * singular, "unrepresentable Kelvin modal rate root")?;
            let rate = checked(
                rate_root * rate_root,
                "unrepresentable Kelvin modal relaxation rate",
            )?;
            if rate <= 0.0 {
                return Err(invalid("Kelvin modal relaxation rate underflowed"));
            }
            self.rates.push(rate);
            // F^T U_i = sqrt(lambda_i) V_i makes this equal to
            // sqrt(lambda_i) G_eq^-1 V_i[0..6], without inverting G_eq.
            let strain_mode = if let Some(constraint) = &self.constrained {
                let mut direction = [Dd::ZERO; 6];
                for (row, loss) in constraint.strain_loss.iter().enumerate() {
                    for axis in 0..6 {
                        direction[axis] = direction[axis]
                            + Dd::from_f64(loss[axis])
                                * Dd::from_f64(decomposition.u[row * n + mode]);
                    }
                }
                let result = direction.map(Dd::to_f64);
                if !finite(&result) {
                    return Err(invalid(
                        "unrepresentable constrained Kelvin modal direction",
                    ));
                }
                result
            } else {
                self.solve_viscosity_root(
                    core::array::from_fn(|i| decomposition.u[i * n + mode]),
                    false,
                )?
            };
            self.strain_modes.push(strain_mode);
            for row in 0..memory {
                let direction = if let Some(basis) = energy_basis {
                    let mut value = Dd::ZERO;
                    for column in 0..n {
                        value = value
                            + Dd::from_f64(basis[(row + 6) * n + column])
                                * Dd::from_f64(decomposition.v[column * n + mode]);
                    }
                    value.to_f64()
                } else {
                    decomposition.v[(row + 6) * n + mode]
                };
                self.memory_modes[row * n + mode] = checked(
                    rate_root * direction,
                    "unrepresentable Maxwell modal rate direction",
                )?;
            }
        }
        self.left_modes = decomposition.u;
        poll(cancelled)?;
        Ok(())
    }
}
