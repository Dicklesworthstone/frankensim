//! Transfer factor coordinates at fixed strain without erasing Maxwell memory.

use super::*;

fn recoverable_factor(
    stiffness: &Tangent6,
    rows: &Tangent6,
    relative_tolerance: f64,
) -> Result<fs_la::factor::Cholesky> {
    let matrix: Vec<_> = stiffness.iter().flatten().copied().collect();
    let factor = fs_la::factor::cholesky(&matrix, 6).map_err(|_| {
        invalid("changed Maxwell factors must preserve all physical spring strain directions")
    })?;
    // A rounded Gram matrix can give a singular B a tiny positive Cholesky
    // pivot. Check every factor-space unit vector, even at zero actual memory.
    for column in 0..6 {
        let mut inverse_column = rows[column];
        factor.solve(&mut inverse_column);
        let recovered_column = apply(rows, &inverse_column);
        for (row, value) in recovered_column.iter().enumerate() {
            let expected = if row == column { 1.0 } else { 0.0 };
            if !value.is_finite() || (value - expected).abs() > relative_tolerance {
                return Err(invalid(
                    "Maxwell factor has unresolved or singular physical strain directions",
                ));
            }
        }
    }
    Ok(factor)
}

impl ProjectedTensorProny {
    /// Re-express each spring strain in another law's factor coordinates.
    /// Both projections must have identical strain columns, reference volume
    /// and branch order. The source-binding owner must additionally establish
    /// that corresponding branches describe the SAME physical internal variable.
    ///
    /// For z=sqrt(V) B e, recover e by C e=B^T z/sqrt(V), C=B^T B,
    /// then form z_new=sqrt(V) B_new e. Changed singular factors refuse at BOTH
    /// endpoints: old hidden strain cannot be reconstructed, and a new nullspace
    /// would erase strain required by later continuation. Exactly unchanged
    /// factors need no inversion, including rank-deficient branches.
    /// No relaxation, strain increment, heat or clock advancement occurs here.
    pub fn transfer_factor_memory(
        &self,
        next: &Self,
        memory: &[f64],
        relative_tolerance: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<f64>> {
        poll(&mut cancelled)?;
        if self.columns != next.columns
            || self.volume_m3.to_bits() != next.volume_m3.to_bits()
            || self.model.branches.len() != next.model.branches.len()
            || memory.len() != 6 * self.model.branches.len()
            || memory.iter().any(|x| !x.is_finite())
            || !relative_tolerance.is_finite()
            || !(0.0..=1e-8).contains(&relative_tolerance)
            || relative_tolerance == 0.0
        {
            return Err(invalid(
                "memory transfer requires identical kinematics, complete ordered branches and a bounded relative tolerance",
            ));
        }
        let root = det::sqrt(self.volume_m3);
        let mut output = Vec::with_capacity(memory.len());
        for (index, (old, new)) in self
            .model
            .branches
            .iter()
            .zip(&next.model.branches)
            .enumerate()
        {
            poll(&mut cancelled)?;
            let old_rows = old.resolved_gram_rows(&self.model.factors[index]);
            let new_rows = new.resolved_gram_rows(&next.model.factors[index]);
            let z = &memory[6 * index..6 * index + 6];
            if old_rows == new_rows {
                output.extend_from_slice(z);
                continue;
            }
            let factor = recoverable_factor(&old.c_pa, &old_rows, relative_tolerance)?;
            recoverable_factor(&new.c_pa, &new_rows, relative_tolerance)?;
            let mut strain: [f64; 6] =
                core::array::from_fn(|i| (0..6).map(|j| old_rows[j][i] * z[j] / root).sum());
            factor.solve(&mut strain);
            if !finite(&strain) {
                return Err(invalid("Maxwell spring-strain recovery overflow"));
            }
            let recovered = apply(&old_rows, &strain).map(|x| root * x);
            let scale = z.iter().fold(0.0_f64, |s, x| s.max(x.abs()));
            let residual = recovered
                .iter()
                .zip(z)
                .fold(0.0_f64, |s, (a, b)| s.max((a - b).abs()));
            if !finite(&recovered) || !residual.is_finite() || residual > relative_tolerance * scale
            {
                return Err(invalid(
                    "Maxwell spring-strain reconstruction exceeds its relative allowance",
                ));
            }
            let transformed = apply(&new_rows, &strain).map(|x| root * x);
            if !finite(&transformed) {
                return Err(invalid("transferred Maxwell factor coordinates overflow"));
            }
            output.extend_from_slice(&transformed);
        }
        poll(&mut cancelled)?;
        Ok(output)
    }
}
