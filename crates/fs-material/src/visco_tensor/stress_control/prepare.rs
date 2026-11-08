//! Positive compliance preparation; all numerical factors use fs-la owners.

use super::*;
use fs_la::factor::{cholesky, qr, svd_jacobi};
use fs_math::det;

impl PreparedTensorStressControl {
    /// Bind a reusable controller to this history's original immutable model.
    /// The history itself is not advanced or changed by preparation.
    ///
    /// # Errors
    /// Nonzero Kelvin viscosity, capacity/work exhaustion, unresolved native
    /// factors, unrepresentable source scales, or cancellation.
    #[allow(clippy::too_many_lines)] // One admitted factorization of the original compliance law.
    pub fn new(
        history: &TensorPronyHistory,
        budget: TensorStressControlBudget,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        Self::prepare(history, [true; 6], budget, cancelled)
    }

    #[allow(clippy::too_many_lines)] // Same positive factor for full and partially constrained creep.
    pub(super) fn prepare(
        history: &TensorPronyHistory,
        stress_components: [bool; 6],
        budget: TensorStressControlBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(&mut cancelled)?;
        let model = &history.model;
        if model
            .kv_pa_s()
            .is_some_and(|c| c.iter().flatten().any(|&v| v != 0.0))
        {
            return Err(invalid(
                "tensor stress control does not admit nonzero Kelvin viscosity",
            ));
        }
        let n = model.branches().len() * 6;
        let storage = n
            .checked_mul(n)
            .and_then(|x| x.checked_mul(16))
            .and_then(|x| x.checked_add(128 * n + 1024));
        let work = n
            .checked_pow(3)
            .and_then(|x| x.checked_mul(1024))
            .and_then(|x| x.checked_add(100 * n * n + 10_000));
        if n > budget.max_memory_coordinates
            || storage.is_none_or(|s| s > budget.max_dense_scalars)
            || work.is_none_or(|w| {
                w > budget.max_preparation_work || w > fs_la::eigen::MAX_EIGEN_UNPOLLED_SCALAR_WORK
            })
        {
            return Err(invalid(
                "tensor stress-control preparation exceeds its capacity or work budget",
            ));
        }
        let mut c0 = [[Dd::ZERO; 6]; 6];
        for (i, row) in c0.iter_mut().enumerate() {
            for (j, value) in row.iter_mut().enumerate() {
                for factor in &model.equilibrium_gram {
                    *value = *value + Dd::from_f64(factor[i]) * Dd::from_f64(factor[j]);
                }
            }
        }
        for (branch, fallback) in model.branches().iter().zip(&model.factors) {
            let b = branch.resolved_gram_rows(fallback);
            for (i, row) in c0.iter_mut().enumerate() {
                for (j, value) in row.iter_mut().enumerate() {
                    for factor in &b {
                        *value = *value + Dd::from_f64(factor[i]) * Dd::from_f64(factor[j]);
                    }
                }
            }
        }
        let c0 = c0.map(|r| r.map(Dd::to_f64));
        if c0.iter().flatten().any(|v| !v.is_finite()) {
            return Err(invalid("unrepresentable instantaneous tensor stiffness"));
        }
        let axes: Vec<_> = (0..6).filter(|&i| stress_components[i]).collect();
        let m = axes.len();
        let scale = if m == 0 {
            1.0
        } else {
            axes.iter()
                .flat_map(|&i| axes.iter().map(move |&j| c0[i][j]))
                .fold(0.0_f64, |a, b| a.max(b.abs()))
        };
        if !scale.is_finite() || scale <= 0.0 {
            return Err(invalid("unrepresentable instantaneous tensor stiffness"));
        }
        let mut flat = [0.0; 36];
        for (i, &axis_i) in axes.iter().enumerate() {
            for (j, &axis_j) in axes.iter().enumerate() {
                flat[m * i + j] = c0[axis_i][axis_j] / scale;
            }
        }
        let instantaneous = if m == 0 {
            None
        } else {
            Some(cholesky(&flat[..m * m], m).map_err(|_| {
                invalid("controlled instantaneous tensor stiffness is not numerically SPD")
            })?)
        };
        let mut prepared = Self {
            model: Arc::clone(model),
            instantaneous,
            instantaneous_stiffness: c0,
            stress_components,
            stiffness_scale: scale,
            branch_rate_roots: Vec::new(),
            rates: Vec::new(),
            modes: Vec::new(),
            forcing: Vec::new(),
        };
        if n == 0 {
            poll(&mut cancelled)?;
            return Ok(prepared);
        }
        // Rotated models retain a nontriangular original Gram factor. QR it
        // directly; treating its diagonal as a Cholesky factor would change
        // the physical compliance after a material-frame rotation.
        let eq_scale = axes
            .iter()
            .flat_map(|&i| model.equilibrium_gram.iter().map(move |row| row[i]))
            .fold(0.0_f64, |a, b| a.max(b.abs()));
        let mut eq_rows = [0.0; 36];
        for i in 0..6 {
            for (j, &axis) in axes.iter().enumerate() {
                eq_rows[m * i + j] = model.equilibrium_gram[i][axis] / eq_scale;
            }
        }
        let eq_factor = (m > 0).then(|| qr(&eq_rows[..6 * m], 6, m));
        let tau_min = model
            .branches()
            .iter()
            .map(|b| b.tau_s)
            .fold(f64::INFINITY, f64::min);
        let rate_scale = checked(1.0 / tau_min, "unrepresentable fastest material rate")?;
        if rate_scale <= 0.0 {
            return Err(invalid("unresolved fastest material rate"));
        }
        let mut normalized_roots = zeros(n)?;
        let mut rows = vec![[0.0; 6]; n];
        let mut compliance_rows = vec![[0.0; 6]; n];
        let mut qr_scale = 1.0_f64;
        for (branch_index, branch) in model.branches().iter().enumerate() {
            poll(&mut cancelled)?;
            let root = checked(
                1.0 / det::sqrt(branch.tau_s),
                "unrepresentable branch rate root",
            )?;
            let normalized = det::sqrt(tau_min / branch.tau_s);
            if root <= 0.0 || normalized <= 0.0 {
                return Err(invalid("material rate range exceeds the creep transform"));
            }
            prepared.branch_rate_roots.push(root);
            let b = branch.resolved_gram_rows(&model.factors[branch_index]);
            for (coordinate, original_row) in b.iter().enumerate() {
                let row = 6 * branch_index + coordinate;
                normalized_roots[row] = normalized;
                rows[row] = *original_row;
                // P = B_F R_eq,FF^-1; retain original six-row equilibrium
                // factors restricted to the stress-controlled columns F.
                let mut p = [0.0; 6];
                for (i, &axis) in axes.iter().enumerate() {
                    let eq_factor = eq_factor
                        .as_ref()
                        .ok_or_else(|| invalid("missing controlled equilibrium factor"))?;
                    let mut value = Dd::from_f64(original_row[axis] / eq_scale);
                    for (j, &coefficient) in p.iter().take(i).enumerate() {
                        value = value - Dd::from_f64(eq_factor.r(j, i)) * Dd::from_f64(coefficient);
                    }
                    p[i] = checked(
                        value.to_f64() / eq_factor.r(i, i),
                        "unrepresentable material compliance factor",
                    )?;
                    qr_scale = qr_scale.max(p[i].abs());
                }
                compliance_rows[row] = p;
            }
        }
        // [I; P^T] = Q Rj gives I + P P^T = Rj^T Rj. Retain I in
        // this rectangular factor rather than adding it to a huge Gram matrix.
        let mut augmented = zeros((n + m) * n)?;
        for row in 0..n {
            augmented[row * n + row] = 1.0 / qr_scale;
            for axis in 0..m {
                augmented[(n + axis) * n + row] = compliance_rows[row][axis] / qr_scale;
            }
        }
        poll(&mut cancelled)?;
        let factor = qr(&augmented, n + m, n);
        poll(&mut cancelled)?;
        let mut relaxation_factor = zeros(n * n)?;
        // W = Rj^-T sqrt(R/rate_scale), and A/rate_scale = W^T W.
        for col in 0..n {
            for row in 0..n {
                let diagonal = factor.r(row, row);
                if !diagonal.is_finite() || diagonal == 0.0 {
                    return Err(invalid("creep compliance QR is numerically singular"));
                }
                let mut value = Dd::from_f64(if row == col {
                    normalized_roots[col] / qr_scale
                } else {
                    0.0
                });
                for k in 0..row {
                    value = value
                        - Dd::from_f64(factor.r(k, row))
                            * Dd::from_f64(relaxation_factor[k * n + col]);
                }
                relaxation_factor[row * n + col] = checked(
                    value.to_f64() / diagonal,
                    "unrepresentable positive creep factor",
                )?;
            }
        }
        poll(&mut cancelled)?;
        let decomposition = svd_jacobi(&relaxation_factor, n, n);
        poll(&mut cancelled)?;
        let tolerance = 8192.0 * f64::EPSILON * n as f64;
        for mode in 0..n {
            let singular = decomposition.sigma[mode];
            if !singular.is_finite() || singular <= 0.0 {
                return Err(invalid("creep relaxation singular value is unresolved"));
            }
            for row in 0..n {
                let mut actual = Dd::ZERO;
                for col in 0..n {
                    actual = actual
                        + Dd::from_f64(relaxation_factor[row * n + col])
                            * Dd::from_f64(decomposition.v[col * n + mode]);
                }
                let residual = (actual
                    - Dd::from_f64(singular) * Dd::from_f64(decomposition.u[row * n + mode]))
                .to_f64();
                if !residual.is_finite() || (residual / singular).abs() > tolerance {
                    return Err(invalid("creep factor SVD residual is unresolved"));
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
                        return Err(invalid("creep relaxation vectors are unresolved"));
                    }
                }
            }
            let rate = checked(
                singular * (singular * rate_scale),
                "unrepresentable creep relaxation rate",
            )?;
            if rate <= 0.0 {
                return Err(invalid("creep relaxation rate underflowed"));
            }
            prepared.rates.push(rate);
            let mut forcing = [Dd::ZERO; 6];
            for (row, original_row) in rows.iter().enumerate() {
                for axis in 0..6 {
                    forcing[axis] = forcing[axis]
                        + Dd::from_f64(decomposition.v[row * n + mode])
                            * Dd::from_f64(prepared.branch_rate_roots[row / 6])
                            * Dd::from_f64(original_row[axis]);
                }
            }
            let h = forcing.map(Dd::to_f64);
            let controlled = prepared.compliance(h)?;
            let coupled = apply(&c0, &controlled);
            // G_F = C0_FF^-1 H_F^T; G_K = H_K^T - C0_KF G_F.
            // The latter includes reaction work under imposed strain, without
            // discarding any original branch coordinates or coupling terms.
            prepared.forcing.push(core::array::from_fn(|axis| {
                if stress_components[axis] {
                    controlled[axis]
                } else {
                    h[axis] - coupled[axis]
                }
            }));
        }
        prepared.modes = decomposition.v;
        poll(&mut cancelled)?;
        Ok(prepared)
    }
}
