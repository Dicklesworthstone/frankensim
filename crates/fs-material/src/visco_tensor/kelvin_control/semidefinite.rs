//! Algebraic elastic directions of the original semidefinite Kelvin tensor.
//!
//! With z=[G_eq e;q], C=[G_eq;B], and N spanning ker(D), A=C N imposes
//! A^T z=N^T stress. Its compliance is J=N(A^T A)^-1 N^T. If Q spans
//! ker(A^T), the square loss factor FQ has rank(D)+6*branches coordinates:
//! F=[H^-T P^T C^T; sqrt(R)], D=P H^T H P^T. The original loss g obeys
//! g'=-FQ(FQ)^T g-K^T stress', while e'=J stress'-K g. Preparation uses
//! the original nonzero Gram columns and native orthogonal factors, never an
//! eigenvalue cutoff, a diagonal regularizer, or a difference of PSD matrices.

use super::*;
use fs_la::factor::qr;
use fs_math::det;

#[derive(Clone, Debug)]
pub(super) struct ConstrainedKelvin {
    /// Columns of P H^-1: turn original stress residuals into physical loss.
    pub(super) range: Vec<Voigt>,
    null: Vec<Voigt>,
    /// J=W W^T; retaining its root avoids squaring a dimensional compliance.
    compliance_root: Vec<Voigt>,
    /// Columns of K, in original Kelvin then Maxwell loss-coordinate order.
    pub(super) strain_loss: Vec<Voigt>,
    tolerance: f64,
}

impl ConstrainedKelvin {
    pub(super) fn compliance(&self, rhs: &Voigt) -> Result<Voigt> {
        let mut answer = [Dd::ZERO; 6];
        for root in &self.compliance_root {
            let value = dd_dot(root, rhs);
            for axis in 0..6 {
                answer[axis] = answer[axis] + Dd::from_f64(root[axis]) * value;
            }
        }
        let answer = answer.map(Dd::to_f64);
        if !finite(&answer) {
            return Err(invalid("unrepresentable Kelvin nullspace compliance"));
        }
        Ok(answer)
    }

    pub(super) fn admit_residual(&self, residual: &[Dd; 6], magnitude: &Voigt) -> Result<()> {
        for direction in &self.null {
            let mut value = Dd::ZERO;
            let mut scale = 0.0;
            for axis in 0..6 {
                value = value + Dd::from_f64(direction[axis]) * residual[axis];
                scale += direction[axis].abs() * magnitude[axis];
            }
            if !scale.is_finite() || value.abs().to_f64() > self.tolerance * scale {
                return Err(invalid(
                    "original Kelvin algebraic stress constraint is unresolved",
                ));
            }
        }
        Ok(())
    }
}

fn triangular(qr: &Qr, size: usize, rhs: Voigt, transpose: bool, scale: f64) -> Result<Voigt> {
    let mut answer = [0.0; 6];
    for index in 0..size {
        let i = if transpose { index } else { size - 1 - index };
        let mut value = Dd::from_f64(rhs[i] / scale);
        if transpose {
            for (j, &v) in answer.iter().take(i).enumerate() {
                value = value - Dd::from_f64(qr.r(j, i)) * Dd::from_f64(v);
            }
        } else {
            for (j, &v) in answer.iter().enumerate().take(size).skip(i + 1) {
                value = value - Dd::from_f64(qr.r(i, j)) * Dd::from_f64(v);
            }
        }
        answer[i] = checked(
            value.to_f64() / qr.r(i, i),
            "unrepresentable constrained Kelvin factor solve",
        )?;
    }
    Ok(answer)
}

fn admit_pivots(qr: &Qr, size: usize) -> Result<()> {
    if (0..size).any(|i| !qr.r(i, i).is_finite() || qr.r(i, i) == 0.0) {
        return Err(invalid(
            "original Kelvin rank or elastic constraint is unresolved",
        ));
    }
    Ok(())
}

impl PreparedTensorKelvinStressControl {
    // The public entry has already admitted the full original dimension and
    // conservative native-factor storage/work ceilings before this branch.
    #[allow(clippy::too_many_lines)]
    pub(super) fn prepare_semidefinite(
        history: &TensorPronyHistory,
        columns: &[usize],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        let model = &history.model;
        let original = 6 + 6 * model.branches().len();
        let rank = columns.len();
        let algebraic = 6 - rank;
        let dynamic = original - algebraic;
        let tolerance = 8192.0 * f64::EPSILON * original as f64;
        let source = model.kv_factor.as_ref().expect("admitted Kelvin factor");
        let scale = source.iter().flatten().fold(0.0_f64, |s, v| s.max(v.abs()));
        if !scale.is_finite() || scale <= 0.0 {
            return Err(invalid("unrepresentable semidefinite Kelvin factor scale"));
        }
        let mut thin = zeros(6 * rank)?;
        for axis in 0..6 {
            for (j, &column) in columns.iter().enumerate() {
                thin[axis * rank + j] = source[axis][column] / scale;
            }
        }
        poll(cancelled)?;
        let range_qr = qr(&thin, 6, rank);
        poll(cancelled)?;
        admit_pivots(&range_qr, rank)?;
        let mut basis = Vec::with_capacity(6);
        for column in 0..6 {
            let mut vector = [0.0; 6];
            vector[column] = 1.0;
            range_qr.apply_q(&mut vector);
            if !finite(&vector) {
                return Err(invalid("unrepresentable Kelvin range/null basis"));
            }
            basis.push(vector);
        }
        for i in 0..6 {
            for j in 0..6 {
                let expected = if i == j { 1.0 } else { 0.0 };
                if (dd_dot(&basis[i], &basis[j]).to_f64() - expected).abs() > tolerance {
                    return Err(invalid("Kelvin range/null orthogonality is unresolved"));
                }
            }
        }
        let mut constraint = ConstrainedKelvin {
            range: Vec::with_capacity(rank),
            null: basis[rank..].to_vec(),
            compliance_root: Vec::with_capacity(algebraic),
            strain_loss: Vec::with_capacity(dynamic),
            tolerance,
        };
        for column in 0..rank {
            let mut rhs = [0.0; 6];
            rhs[column] = 1.0;
            let coefficients = triangular(&range_qr, rank, rhs, true, scale)?;
            let direction = core::array::from_fn(|axis| {
                (0..rank)
                    .fold(Dd::ZERO, |value, j| {
                        value + Dd::from_f64(basis[j][axis]) * Dd::from_f64(coefficients[j])
                    })
                    .to_f64()
            });
            if !finite(&direction) {
                return Err(invalid("unrepresentable Kelvin range inverse root"));
            }
            constraint.range.push(direction);
        }
        // Establish the exact original rank, including weak positive channels.
        // Numerical failure refuses this law; it never changes that source rank.
        for (j, &column) in columns.iter().enumerate() {
            let source_column = core::array::from_fn(|axis| source[axis][column]);
            for (i, direction) in constraint.range.iter().enumerate() {
                let expected = if i == j { 1.0 } else { 0.0 };
                if (dd_dot(direction, &source_column).to_f64() - expected).abs() > tolerance {
                    return Err(invalid(
                        "original Kelvin inverse-root reconstruction is unresolved",
                    ));
                }
            }
            for direction in &constraint.null {
                let residual = dd_dot(direction, &source_column).abs().to_f64();
                let magnitude = (0..6)
                    .map(|i| (direction[i] * source_column[i]).abs())
                    .sum::<f64>();
                if !magnitude.is_finite() || residual > tolerance * magnitude {
                    return Err(invalid(
                        "original Kelvin nullspace reconstruction is unresolved",
                    ));
                }
            }
        }
        let mut state_rows = Vec::with_capacity(original);
        state_rows.extend(model.equilibrium_gram);
        let mut roots = Vec::with_capacity(model.branches().len());
        for (j, branch) in model.branches().iter().enumerate() {
            poll(cancelled)?;
            let root = checked(
                1.0 / det::sqrt(branch.tau_s),
                "unrepresentable Maxwell rate root",
            )?;
            if root <= 0.0 {
                return Err(invalid("Maxwell rate root underflowed"));
            }
            roots.push(root);
            state_rows.extend(branch.resolved_gram_rows(&model.factors[j]));
        }
        let mut a = zeros(original * algebraic)?;
        for (row, c) in state_rows.iter().enumerate() {
            for (column, direction) in constraint.null.iter().enumerate() {
                a[row * algebraic + column] = checked(
                    dd_dot(c, direction).to_f64(),
                    "unrepresentable Kelvin elastic constraint",
                )?;
            }
        }
        let elastic_scale = a.iter().fold(0.0_f64, |s, v| s.max(v.abs()));
        if !elastic_scale.is_finite() || elastic_scale <= 0.0 {
            return Err(invalid("unrepresentable Kelvin elastic-constraint scale"));
        }
        for value in &mut a {
            *value /= elastic_scale;
        }
        poll(cancelled)?;
        let elastic_qr = qr(&a, original, algebraic);
        poll(cancelled)?;
        admit_pivots(&elastic_qr, algebraic)?;
        for column in 0..algebraic {
            let mut rhs = [0.0; 6];
            rhs[column] = 1.0;
            let coefficients = triangular(&elastic_qr, algebraic, rhs, false, elastic_scale)?;
            let root = core::array::from_fn(|axis| {
                (0..algebraic)
                    .fold(Dd::ZERO, |value, j| {
                        value
                            + Dd::from_f64(constraint.null[j][axis]) * Dd::from_f64(coefficients[j])
                    })
                    .to_f64()
            });
            if !finite(&root) {
                return Err(invalid("unrepresentable Kelvin elastic compliance root"));
            }
            constraint.compliance_root.push(root);
        }
        for direction in &constraint.range {
            let mut instantaneous = [Dd::ZERO; 6];
            for c in &state_rows {
                let value = dd_dot(c, direction);
                for axis in 0..6 {
                    instantaneous[axis] = instantaneous[axis] + Dd::from_f64(c[axis]) * value;
                }
            }
            let correction = constraint.compliance(&instantaneous.map(Dd::to_f64))?;
            let loss = core::array::from_fn(|i| {
                (Dd::from_f64(direction[i]) - Dd::from_f64(correction[i])).to_f64()
            });
            constraint.strain_loss.push(loss);
        }
        for (row, c) in state_rows.iter().enumerate().skip(6) {
            let response = constraint.compliance(c)?;
            let loss = response.map(|v| -v * roots[(row - 6) / 6]);
            if !finite(&loss) {
                return Err(invalid(
                    "unrepresentable constrained Maxwell loss direction",
                ));
            }
            constraint.strain_loss.push(loss);
        }
        let mut energy_basis = zeros(original * dynamic)?;
        let mut factor = zeros(dynamic * dynamic)?;
        for column in 0..dynamic {
            poll(cancelled)?;
            let mut vector = zeros(original)?;
            vector[algebraic + column] = 1.0;
            elastic_qr.apply_q(&mut vector);
            let mut force = [Dd::ZERO; 6];
            for (row, c) in state_rows.iter().enumerate() {
                energy_basis[row * dynamic + column] = vector[row];
                for axis in 0..6 {
                    force[axis] = force[axis] + Dd::from_f64(c[axis]) * Dd::from_f64(vector[row]);
                }
            }
            let force = force.map(Dd::to_f64);
            for (row, direction) in constraint.range.iter().enumerate() {
                factor[row * dynamic + column] = checked(
                    dd_dot(direction, &force).to_f64(),
                    "unrepresentable constrained Kelvin loss factor",
                )?;
            }
            for row in 6..original {
                factor[(rank + row - 6) * dynamic + column] = checked(
                    roots[(row - 6) / 6] * vector[row],
                    "unrepresentable constrained Maxwell loss factor",
                )?;
            }
        }
        let mut prepared = Self {
            model: Arc::clone(model),
            viscosity_qr: None,
            viscosity_scale: scale,
            constrained: Some(constraint),
            branch_rate_roots: roots,
            rates: Vec::new(),
            left_modes: Vec::new(),
            strain_modes: Vec::new(),
            memory_modes: Vec::new(),
        };
        prepared.prepare_modes(factor, dynamic, Some(&energy_basis), cancelled)?;
        Ok(prepared)
    }
}
