//! Transfer factor coordinates at fixed strain without erasing Maxwell memory.

use super::*;

/// Original point's mechanical thermodynamic contributions. Entropy and U
/// are signed contributions relative to its unstressed caloric material.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct TensorThermalState {
    pub free_energy_j: f64,
    pub entropy_j_k: f64,
    pub internal_energy_j: f64,
}

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

fn recover_strain(
    stiffness: &Tangent6,
    rows: &Tangent6,
    z: &[f64],
    root: f64,
    tolerance: f64,
) -> Result<Voigt> {
    let factor = recoverable_factor(stiffness, rows, tolerance)?;
    let mut strain = core::array::from_fn(|i| (0..6).map(|j| rows[j][i] * z[j] / root).sum());
    factor.solve(&mut strain);
    if !finite(&strain) {
        return Err(invalid("Maxwell spring-strain recovery overflow"));
    }
    let recovered = apply(rows, &strain).map(|x| root * x);
    let scale = z.iter().fold(0.0_f64, |s, x| s.max(x.abs()));
    let residual = recovered
        .iter()
        .zip(z)
        .fold(0.0_f64, |s, (a, b)| s.max((a - b).abs()));
    if !finite(&recovered) || !residual.is_finite() || residual > tolerance * scale {
        return Err(invalid(
            "Maxwell spring-strain reconstruction exceeds its relative allowance",
        ));
    }
    Ok(strain)
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
    /// If free strain changes, subtract B_new*(free_new-free_old) from each
    /// new factor coordinate. Thus the original VISCOUS strain survives the
    /// thermal strain change; keeping spring strain fixed would erase thermal
    /// stress. No relaxation, dissipative heat or clock advancement occurs here.
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
        let free_change =
            core::array::from_fn(|i| next.free_strain_mandel[i] - self.free_strain_mandel[i]);
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
            let shift = apply(&new_rows, &free_change).map(|x| root * x);
            if old_rows == new_rows {
                if free_change.iter().all(|x| *x == 0.0) {
                    output.extend_from_slice(z);
                } else {
                    let transformed: Vec<_> = z.iter().zip(shift).map(|(x, d)| x - d).collect();
                    if transformed.iter().any(|v| !v.is_finite()) {
                        return Err(invalid("thermal memory transfer overflow"));
                    }
                    output.extend(transformed);
                }
                continue;
            }
            let strain = recover_strain(&old.c_pa, &old_rows, z, root, relative_tolerance)?;
            recoverable_factor(&new.c_pa, &new_rows, relative_tolerance)?;
            let transformed =
                core::array::from_fn::<_, 6, _>(|i| root * dot(&new_rows[i], &strain) - shift[i]);
            if !finite(&transformed) {
                return Err(invalid("transferred Maxwell factor coordinates overflow"));
            }
            output.extend_from_slice(&transformed);
        }
        poll(&mut cancelled)?;
        Ok(output)
    }

    /// Evaluate psi, S=-psi_T and U=psi+T*S using this ORIGINAL projected
    /// strain map and complete volume-scaled factor memory. Derivatives must
    /// come from the same source law and reference frame as this projection;
    /// the source/coupling owner establishes that identity.
    ///
    /// Every varying branch recovers physical elastic strain through the same
    /// checked factor inverse used by native memory transfer. A varying singular
    /// branch refuses even at zero actual memory. An unchanged singular branch
    /// contributes its original free energy and zero thermal entropy.
    pub fn thermal_state(
        &self,
        derivative: &source::thermal::TensorTemperatureDerivative,
        displacement: &[f64],
        memory: &[f64],
        temperature_k: f64,
        relative_memory_tolerance: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<TensorThermalState> {
        self.thermal_state_with_expansion(
            derivative,
            displacement,
            memory,
            temperature_k,
            [0.0; 6],
            relative_memory_tolerance,
            cancelled,
        )
    }

    /// The same complete potential, with a source-admitted free-strain
    /// derivative in this point's Mandel frame. At fixed total/viscous strains,
    /// S includes V*alpha:sigma for equilibrium AND every Maxwell spring.
    /// Constant singular branches need no inverse: B*alpha dotted with their
    /// original factor memory retains all energetic directions exactly.
    #[allow(clippy::too_many_arguments)]
    pub fn thermal_state_with_expansion(
        &self,
        derivative: &source::thermal::TensorTemperatureDerivative,
        displacement: &[f64],
        memory: &[f64],
        temperature_k: f64,
        alpha_mandel: Voigt,
        relative_memory_tolerance: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<TensorThermalState> {
        poll(&mut cancelled)?;
        if memory.len() != 6 * self.model.branches.len()
            || derivative.branches_pa_k.len() != self.model.branches.len()
            || memory.iter().any(|v| !v.is_finite())
            || !(derivative.interval_k[0]..=derivative.interval_k[1]).contains(&temperature_k)
            || !relative_memory_tolerance.is_finite()
            || relative_memory_tolerance <= 0.0
            || relative_memory_tolerance > 1e-8
            || !finite(&alpha_mandel)
        {
            return Err(invalid(
                "thermal tensor state needs complete memory, bounded temperature and reconstruction tolerance",
            ));
        }
        let total = mandel(&self.strain(displacement)?);
        let strain = core::array::from_fn(|i| total[i] - self.free_strain_mandel[i]);
        let n = self.coordinates();
        let mut free_energy_j = 0.0;
        for i in 0..n {
            for j in 0..n {
                free_energy_j +=
                    0.5 * displacement[i] * self.equilibrium[i * n + j] * displacement[j];
            }
        }
        if self.free_strain_mandel.iter().any(|x| *x != 0.0) {
            // Evaluate physical elastic strain directly; subtracting nearly
            // equal affine terms would lose free-expansion energy accuracy.
            free_energy_j = 0.5 * self.volume_m3 * dot(&strain, &apply(&self.model.c_eq, &strain));
        }
        let mut entropy_j_k =
            -0.5 * self.volume_m3 * dot(&strain, &apply(&derivative.equilibrium_pa_k, &strain));
        entropy_j_k += self.volume_m3 * dot(&alpha_mandel, &apply(&self.model.c_eq, &strain));
        let root = det::sqrt(self.volume_m3);
        for (index, (branch, slope)) in self
            .model
            .branches
            .iter()
            .zip(&derivative.branches_pa_k)
            .enumerate()
        {
            poll(&mut cancelled)?;
            let z = &memory[6 * index..6 * index + 6];
            free_energy_j += 0.5 * z.iter().map(|v| v * v).sum::<f64>();
            if alpha_mandel.iter().any(|v| *v != 0.0) {
                let rows = branch.resolved_gram_rows(&self.model.factors[index]);
                let rate = apply(&rows, &alpha_mandel);
                entropy_j_k += root * rate.iter().zip(z).map(|(a, b)| a * b).sum::<f64>();
            }
            if slope.iter().flatten().any(|v| *v != 0.0) {
                let rows = branch.resolved_gram_rows(&self.model.factors[index]);
                let q = recover_strain(&branch.c_pa, &rows, z, root, relative_memory_tolerance)?;
                if alpha_mandel.iter().any(|v| *v != 0.0) {
                    // The combined capacity bound assumes ||q|| <= 2L.
                    // Use the SAME native Mandel-norm admission as total and
                    // equilibrium elastic strain, scaled by one half.
                    self.model.validate_mandel_strain(&q.map(|v| 0.5 * v))
                        .map_err(|_| invalid("affine expansion Maxwell spring strain exceeds its 2L thermal envelope"))?;
                }
                entropy_j_k -= 0.5 * self.volume_m3 * dot(&q, &apply(slope, &q));
            }
        }
        let internal_energy_j = temperature_k.mul_add(entropy_j_k, free_energy_j);
        if !free_energy_j.is_finite() || !entropy_j_k.is_finite() || !internal_energy_j.is_finite()
        {
            return Err(invalid(
                "tensor mechanical thermodynamic contribution overflow",
            ));
        }
        poll(&mut cancelled)?;
        Ok(TensorThermalState {
            free_energy_j,
            entropy_j_k,
            internal_energy_j,
        })
    }
}
