//! Transfer factor coordinates at fixed strain without erasing Maxwell memory.

use super::*;

/// Original point's mechanical thermodynamic contributions. Entropy and U
/// are signed contributions relative to its unstressed caloric material.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct TensorThermalState {
    pub free_energy_j: f64,
    pub entropy_j_k: f64,
    pub internal_energy_j: f64,
    /// Mechanical dU/dT at fixed total and viscous strains [J/K]. The
    /// caloric material contributes separately; only their total must be positive.
    pub heat_capacity_j_k: f64,
}

// A right inverse on the ORIGINAL nonzero Gram rows. Normalize each row
// independently before solving: a weak positive spring direction is never
// classified as zero from its magnitude relative to another direction.
// Only exactly zero original rows are absent. Redundant or numerically
// unresolved nonzero rows refuse instead of inventing a spectral cutoff.
struct FactorRecovery {
    rows: Tangent6,
    scales: Voigt,
    dual: Tangent6,
    projector: Tangent6,
    rank: usize,
}
impl FactorRecovery {
    fn new(original: &Tangent6, tolerance: f64) -> Result<Self> {
        let scales = original.map(|row| row.iter().fold(0.0_f64, |s, x| s.max(x.abs())));
        let active: Vec<_> = (0..6).filter(|&i| scales[i] != 0.0).collect();
        let rows = core::array::from_fn(|i| {
            if scales[i] == 0.0 { [0.0; 6] } else { original[i].map(|x| x / scales[i]) }
        });
        let rank = active.len();
        let mut dual = [[0.0; 6]; 6];
        if rank != 0 {
            let gram: Vec<_> = active.iter().flat_map(|&i|
                active.iter().map(move |&j| dot(&rows[i], &rows[j]))).collect();
            let factor = fs_la::factor::cholesky(&gram, rank).map_err(|_| {
                invalid("Maxwell nonzero factor rows must have independently recoverable directions")
            })?;
            for (column, &index) in active.iter().enumerate() {
                let mut weights = vec![0.0; rank];
                weights[column] = 1.0;
                factor.solve(&mut weights);
                dual[index] = core::array::from_fn(|i| active.iter().zip(&weights)
                    .map(|(&row, weight)| rows[row][i] * weight).sum());
                for &row in &active {
                    let value = dot(&rows[row], &dual[index]);
                    let expected = if row == index { 1.0 } else { 0.0 };
                    if !value.is_finite() || (value - expected).abs() > tolerance {
                        return Err(invalid("Maxwell factor row inverse exceeds its relative allowance"));
                    }
                }
            }
        }
        let projector: Tangent6 = core::array::from_fn(|i| core::array::from_fn(|j|
            active.iter().map(|&row| dual[row][i] * rows[row][j]).sum()));
        for i in 0..6 {
            for j in 0..6 {
                let value = projector[i][j];
                if !value.is_finite() || value.abs() > 1.0 + tolerance
                    || (value - projector[j][i]).abs() > tolerance {
                    return Err(invalid("Maxwell physical factor range is numerically unresolved"));
                }
            }
        }
        Ok(Self { rows, scales, dual, projector, rank })
    }

    fn same_range(&self, other: &Self, tolerance: f64) -> Result<()> {
        if self.rank != other.rank || self.projector.iter().flatten()
            .zip(other.projector.iter().flatten()).any(|(a, b)| (a-b).abs() > tolerance) {
            return Err(invalid("changed Maxwell factors must preserve the same physical strain range"));
        }
        Ok(())
    }

    // Entropy must not depend on any unobserved strain in the original
    // factor's nullspace. Normalize every derivative row independently too.
    fn admit_derivative(&self, derivative: &Tangent6, tolerance: f64) -> Result<()> {
        for row in derivative {
            if !finite(row) { return Err(invalid("nonfinite Maxwell thermal derivative")); }
            let scale = row.iter().fold(0.0_f64, |s, x| s.max(x.abs()));
            if scale == 0.0 { continue; }
            let direction = row.map(|x| x / scale);
            let represented = apply(&self.projector, &direction);
            if represented.iter().zip(direction).any(|(a, b)| (a-b).abs() > tolerance) {
                return Err(invalid("Maxwell thermal derivative requires an unobserved strain direction"));
            }
        }
        Ok(())
    }

    fn strain(&self, z: &[f64], root: f64, tolerance: f64) -> Result<Voigt> {
        let mut coordinates = [0.0; 6];
        for i in 0..6 {
            if self.scales[i] == 0.0 {
                if z[i] != 0.0 {
                    return Err(invalid("Maxwell memory lies outside its original factor range"));
                }
            } else {
                let denominator = root * self.scales[i];
                coordinates[i] = if denominator.is_finite() && denominator != 0.0 {
                    z[i] / denominator
                } else { (z[i] / root) / self.scales[i] };
                if !coordinates[i].is_finite() || (z[i] != 0.0 && coordinates[i] == 0.0) {
                    return Err(invalid("Maxwell observable strain is not representable"));
                }
            }
        }
        let strain = core::array::from_fn(|i|
            (0..6).map(|row| self.dual[row][i] * coordinates[row]).sum());
        let recovered = apply(&self.rows, &strain);
        let scale = coordinates.iter().fold(0.0_f64, |s, x| s.max(x.abs()));
        let residual = recovered.iter().zip(coordinates)
            .fold(0.0_f64, |s, (a, b)| s.max((a-b).abs()));
        if !finite(&strain) || !finite(&recovered) || !residual.is_finite()
            || residual > tolerance * scale {
            return Err(invalid("Maxwell spring-strain reconstruction exceeds its relative allowance"));
        }
        Ok(strain)
    }
}

impl ProjectedTensorProny {
    /// Re-express each spring strain in another law's factor coordinates.
    /// Both projections must have identical strain columns, reference volume
    /// and branch order. The source-binding owner must additionally establish
    /// that corresponding branches describe the SAME physical internal variable.
    ///
    /// For z=sqrt(V) B e, recover the observable strain through the original
    /// independently normalized nonzero factor rows, then form sqrt(V) B_new e.
    /// Rank-deficient factors may change when their physical row spaces agree.
    /// A changed rank/nullspace or unresolved dependent rows refuse: hidden
    /// strain cannot be invented or erased for later continuation. Every nonzero
    /// row survives, however weak; there is no eigenvalue cutoff. Exactly
    /// unchanged factors need no inversion, including rank-deficient branches.
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
            let old_recovery = FactorRecovery::new(&old_rows, relative_tolerance)?;
            let new_recovery = FactorRecovery::new(&new_rows, relative_tolerance)?;
            old_recovery.same_range(&new_recovery, relative_tolerance)?;
            let strain = old_recovery.strain(z, root, relative_tolerance)?;
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
    /// Every varying branch recovers observable elastic strain through the same
    /// checked factor inverse used by native memory transfer. Rank-deficient
    /// branches are admitted only when their thermal derivative also acts within
    /// that physical range; derivatives requiring hidden strain refuse even at
    /// zero memory. An unchanged branch needs no strain reconstruction.
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

    /// The same complete potential, with a source-admitted constant free-strain
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
        cancelled: impl FnMut() -> bool,
    ) -> Result<TensorThermalState> {
        let expansion = TensorThermalExpansion::try_new(
            temperature_k, alpha_mandel, [0.0; 6], derivative.interval_k,
        )?;
        self.thermal_state_with_expansion_derivatives(
            derivative, displacement, memory, temperature_k, &expansion,
            relative_memory_tolerance, cancelled,
        )
    }

    /// Complete original potential for an admitted affine alpha(T).
    /// At fixed total/viscous strain, e_T=-alpha and e_TT=-alpha'.
    /// Constant singular Maxwell branches evaluate alpha' work directly from
    /// retained factors, and admit its observable 2L envelope without inversion.
    #[allow(clippy::too_many_arguments)]
    pub fn thermal_state_with_expansion_derivatives(
        &self,
        derivative: &source::thermal::TensorTemperatureDerivative,
        displacement: &[f64],
        memory: &[f64],
        temperature_k: f64,
        expansion: &TensorThermalExpansion,
        relative_memory_tolerance: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<TensorThermalState> {
        poll(&mut cancelled)?;
        expansion.admit(derivative, temperature_k)?;
        let alpha_mandel = expansion.alpha_mandel_per_k();
        let beta_mandel = expansion.alpha_slope_mandel_per_k2();
        if memory.len() != 6 * self.model.branches.len()
            || derivative.branches_pa_k.len() != self.model.branches.len()
            || memory.iter().any(|v| !v.is_finite())
            || !(derivative.interval_k[0]..=derivative.interval_k[1]).contains(&temperature_k)
            || derivative
                .resolved_temperature_k
                .is_some_and(|t| t.to_bits() != temperature_k.to_bits())
            || !relative_memory_tolerance.is_finite()
            || relative_memory_tolerance <= 0.0
            || relative_memory_tolerance > 1e-8
            || !finite(&alpha_mandel)
        {
            return Err(invalid(
                "thermal tensor state needs complete memory, matching derivative temperature and reconstruction tolerance",
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
        let slope_strain = apply(&derivative.equilibrium_pa_k, &strain);
        let mut entropy_j_k = -0.5 * self.volume_m3 * dot(&strain, &slope_strain);
        let has_expansion = alpha_mandel.iter().any(|v| *v != 0.0);
        let changing_expansion = beta_mandel.iter().any(|v| *v != 0.0);
        // At fixed total/viscous strain, e_T=-alpha and e_TT=-beta:
        // psi_TT/V = .5 e:C_TT:e - 2 alpha:C_T:e + alpha:C:alpha - beta:C:e.
        // The affine Maxwell branches below have no C_TT contribution.
        let mut curvature = if derivative
            .equilibrium_pa_k2
            .iter()
            .flatten()
            .any(|v| *v != 0.0)
        {
            0.5 * dot(&strain, &apply(&derivative.equilibrium_pa_k2, &strain))
        } else {
            0.0
        };
        if has_expansion {
            entropy_j_k += self.volume_m3 * dot(&alpha_mandel, &apply(&self.model.c_eq, &strain));
            curvature += dot(&alpha_mandel, &apply(&self.model.c_eq, &alpha_mandel))
                - 2.0 * dot(&alpha_mandel, &slope_strain);
        }
        if changing_expansion {
            curvature -= dot(&beta_mandel, &apply(&self.model.c_eq, &strain));
        }
        let root = det::sqrt(self.volume_m3);
        for (index, (branch, slope)) in self
            .model
            .branches
            .iter()
            .zip(derivative.branches_pa_k.iter())
            .enumerate()
        {
            poll(&mut cancelled)?;
            let z = &memory[6 * index..6 * index + 6];
            free_energy_j += 0.5 * z.iter().map(|v| v * v).sum::<f64>();
            if has_expansion {
                let rows = branch.resolved_gram_rows(&self.model.factors[index]);
                let rate = apply(&rows, &alpha_mandel);
                entropy_j_k += root * rate.iter().zip(z).map(|(a, b)| a * b).sum::<f64>();
                curvature += dot(&alpha_mandel, &apply(&branch.c_pa, &alpha_mandel));
            }
            if changing_expansion {
                let rows = branch.resolved_gram_rows(&self.model.factors[index]);
                let rate = apply(&rows, &beta_mandel);
                let observable = rate.iter().zip(z).map(|(a, b)| a * b).sum::<f64>() / root;
                let direction = apply(&branch.c_pa, &beta_mandel);
                let allowed = THERMAL_EXPANSION_ENVELOPE_GUARD
                    * (2.0 * self.model.maximum_strain_norm()) * det::sqrt(dot(&direction, &direction));
                if !observable.is_finite() || !allowed.is_finite() || observable.abs() > allowed {
                    return Err(invalid("thermal expansion Maxwell curvature exceeds its observable 2L envelope"));
                }
                // beta:C:q = (B beta).z/sqrt(V), including singular branches.
                curvature -= observable;
            }
            if slope.iter().flatten().any(|v| *v != 0.0) {
                let rows = branch.resolved_gram_rows(&self.model.factors[index]);
                let recovery = FactorRecovery::new(&rows, relative_memory_tolerance)?;
                recovery.admit_derivative(slope, relative_memory_tolerance)?;
                let q = recovery.strain(z, root, relative_memory_tolerance)?;
                if has_expansion || changing_expansion {
                    // The combined capacity bound assumes ||q|| <= 2L.
                    // Use the SAME native Mandel-norm admission as total and
                    // equilibrium elastic strain, scaled by one half.
                    self.model.validate_mandel_strain(&q.map(|v| 0.5 * v))
                        .map_err(|_| invalid("affine expansion Maxwell spring strain exceeds its 2L thermal envelope"))?;
                    curvature -= 2.0 * dot(&alpha_mandel, &apply(slope, &q));
                }
                entropy_j_k -= 0.5 * self.volume_m3 * dot(&q, &apply(slope, &q));
            }
        }
        let internal_energy_j = temperature_k.mul_add(entropy_j_k, free_energy_j);
        let heat_capacity_j_k = if curvature == 0.0 {
            0.0
        } else {
            -temperature_k * (self.volume_m3 * curvature)
        };
        if !free_energy_j.is_finite()
            || !entropy_j_k.is_finite()
            || !internal_energy_j.is_finite()
            || !curvature.is_finite()
            || !heat_capacity_j_k.is_finite()
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
            heat_capacity_j_k,
        })
    }
}
