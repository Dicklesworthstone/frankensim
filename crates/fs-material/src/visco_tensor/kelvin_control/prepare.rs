//! Prepare the full passive Kelvin/Maxwell operator without forming its Gram.

use super::*;
use fs_la::factor::{qr, svd_jacobi};
use fs_math::det;

impl PreparedTensorKelvinStressControl {
    /// Bind this history's immutable material and factor its complete coupled
    /// constitutive operator once. Preparation does not advance the history.
    ///
    /// # Errors
    /// Missing or singular Kelvin viscosity, exceeded dimension/storage/work
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
            .ok_or_else(|| invalid("stress control requires an SPD Kelvin viscosity tensor"))?;
        // The original PSD owner leaves zero columns for missing dissipative
        // directions. Rotation preserves those columns; do not promote a
        // numerically tiny QR pivot of a rank-deficient tensor to viscosity.
        if (0..6).any(|j| kelvin_factor.iter().all(|row| row[j] == 0.0)) {
            return Err(invalid(
                "singular Kelvin viscosity requires a constrained stress controller",
            ));
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
            viscosity_qr,
            viscosity_scale,
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
        let scale = factor.iter().fold(0.0_f64, |s, v| s.max(v.abs()));
        if !scale.is_finite() || scale <= 0.0 {
            return Err(invalid("unrepresentable Kelvin relaxation scale"));
        }
        for value in &mut factor {
            *value /= scale;
        }
        poll(&mut cancelled)?;
        let decomposition = svd_jacobi(&factor, n, n);
        poll(&mut cancelled)?;
        let tolerance = 8192.0 * f64::EPSILON * n as f64;
        prepared.memory_modes = zeros((n - 6) * n)?;
        for mode in 0..n {
            poll(&mut cancelled)?;
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
            prepared.rates.push(rate);
            // F^T U_i = sqrt(lambda_i) V_i makes this equal to
            // sqrt(lambda_i) G_eq^-1 V_i[0..6], without inverting G_eq.
            prepared.strain_modes.push(prepared.solve_viscosity_root(
                core::array::from_fn(|i| decomposition.u[i * n + mode]),
                false,
            )?);
            for row in 0..n - 6 {
                prepared.memory_modes[row * n + mode] = checked(
                    rate_root * decomposition.v[(row + 6) * n + mode],
                    "unrepresentable Maxwell modal rate direction",
                )?;
            }
        }
        prepared.left_modes = decomposition.u;
        poll(&mut cancelled)?;
        Ok(prepared)
    }
}
