//! Bulk stress from the retained Q1 field and the assembly's reference tensor.
//!
//! Voigt entries are `[xx, yy, zz, xy, yz, xz]`, with physical shear stresses
//! (no engineering-strain factor in the output). These are numerical samples,
//! not recovered nodal stresses or bounds on a continuum maximum. Ghost and
//! Nitsche stabilization affect the solved field but are not material stress.
use super::material::ReferenceMaterial3;
use super::*;

/// One original bulk quadrature point; point order is cell order followed by
/// that cell's retained rule order. No geometry is regenerated or reclassified.
#[derive(Debug, Clone, PartialEq)]
pub struct BulkStressPoint3 {
    /// Index in the operator's `scales()` and `volumes()` arrays.
    pub cell: usize,
    pub position: [f64; 3],
    /// Positive reference-volume quadrature weight, not normalized.
    pub weight: f64,
    /// `C_reference : eps(u)` for the original admitted material.
    pub reference_stress: [f64; 6],
    /// Current physical stress: `scale[cell] * reference_stress`.
    pub stress: [f64; 6],
}

/// One affine displacement family, sampled on the SAME retained bulk points.
/// For coefficients with Euclidean norm at most one, the reference stress is
/// `baseline[p].reference_stress + sum_j a[j]*basis_reference_stresses[j][p]`.
/// Source-bound material strain domains are admitted over this entire ball,
/// not by treating coordinate directions as independent physical states.
#[derive(Debug, Clone, PartialEq)]
pub struct AffineBulkStress3 {
    /// Offset-field stresses and retained quadrature metadata, including scales.
    pub baseline: Vec<BulkStressPoint3>,
    /// Field-major: original basis order, then exactly the baseline point order.
    /// These are unscaled constitutive directions, with physical shear entries.
    pub basis_reference_stresses: Vec<Vec<[f64; 6]>>,
}

/// Partial derivatives of a scalar functional of the physical stress samples.
/// These do not include the effect of equilibrium changing with material.
#[derive(Debug, Clone, PartialEq)]
pub struct BulkStressPullback3 {
    /// `dJ/du` in the operator's independent displacement coordinates.
    pub displacement: Vec<f64>,
    /// `dJ/dscale` at FIXED u, one value per retained active cell.
    pub scales: Vec<f64>,
}

fn poll(checkpoint: &mut impl FnMut() -> ControlFlow<()>) -> Result<(), ElasticityError3> {
    if checkpoint().is_break() {
        Err(ElasticityError3::Cancelled)
    } else {
        Ok(())
    }
}

fn admit_displacement(
    u: &[f64],
    n: usize,
    checkpoint: &mut impl FnMut() -> ControlFlow<()>,
) -> Result<(), ElasticityError3> {
    poll(checkpoint)?;
    if u.len() != n {
        return Err(ElasticityError3::Invalid(
            "invalid bulk stress displacement shape",
        ));
    }
    for chunk in u.chunks(256) {
        poll(checkpoint)?;
        if chunk.iter().any(|v| !v.is_finite()) {
            return Err(ElasticityError3::Invalid(
                "nonfinite bulk stress displacement",
            ));
        }
    }
    Ok(())
}

fn stress_basis(material: &ReferenceMaterial3, gradients: &[[f64; 3]; 8]) -> [[f64; 6]; 24] {
    std::array::from_fn(|i| material.basis_stress(gradients[i / 3], i % 3))
}

fn apply_stress_basis(basis: &[[f64; 6]; 24], local: &[f64; 24]) -> [f64; 6] {
    let mut value = [0.0; 6];
    for i in 0..24 {
        for k in 0..6 {
            value[k] += basis[i][k] * local[i];
        }
    }
    value
}

fn stress(material: &ReferenceMaterial3, gradients: &[[f64; 3]; 8], local: &[f64; 24]) -> [f64; 6] {
    apply_stress_basis(&stress_basis(material, gradients), local)
}

fn mandel_strain(gradients: &[[f64; 3]; 8], local: &[f64; 24]) -> [f64; 6] {
    let mut gradient = [[0.0; 3]; 3];
    for a in 0..8 {
        for c in 0..3 {
            for d in 0..3 {
                gradient[c][d] += local[3 * a + c] * gradients[a][d];
            }
        }
    }
    let shear = core::f64::consts::FRAC_1_SQRT_2;
    [
        gradient[0][0],
        gradient[1][1],
        gradient[2][2],
        shear * (gradient[0][1] + gradient[1][0]),
        shear * (gradient[1][2] + gradient[2][1]),
        shear * (gradient[2][0] + gradient[0][2]),
    ]
}

// Admit the actual displacement strain, never the unit basis functions or
// arbitrary adjoint vectors. The reference VJP does not evaluate a primal.
fn admit_material_strain(
    material: &ReferenceMaterial3,
    gradients: &[[f64; 3]; 8],
    local: &[f64; 24],
) -> Result<(), ElasticityError3> {
    let Some(source) = material.source_state() else {
        return Ok(());
    };
    let strain = mandel_strain(gradients, local);
    if strain
        .iter()
        .any(|v| !v.is_finite() || v.abs() > source.maximum_absolute_mandel_strain())
    {
        return Err(ElasticityError3::Invalid(
            "bulk strain exceeds retained material domain",
        ));
    }
    Ok(())
}

fn admit_affine_material_strain(
    material: &ReferenceMaterial3,
    gradients: &[[f64; 3]; 8],
    offset: &[f64; 24],
    basis: impl IntoIterator<Item = [f64; 24]>,
    checkpoint: &mut impl FnMut() -> ControlFlow<()>,
) -> Result<(), ElasticityError3> {
    let Some(source) = material.source_state() else {
        return Ok(());
    };
    let baseline = mandel_strain(gradients, offset);
    let mut norm_scale = [0.0_f64; 6];
    let mut norm_sum = [0.0_f64; 6];
    for local in basis {
        poll(checkpoint)?;
        let strain = mandel_strain(gradients, &local);
        for k in 0..6 {
            let value = strain[k].abs();
            if !value.is_finite() {
                return Err(ElasticityError3::Invalid("nonfinite affine bulk strain"));
            }
            // Scaled sum of squares avoids squaring a large or tiny direction.
            if value > norm_scale[k] {
                let ratio = norm_scale[k] / value;
                norm_sum[k] = 1.0 + norm_sum[k] * ratio * ratio;
                norm_scale[k] = value;
            } else if value != 0.0 {
                let ratio = value / norm_scale[k];
                norm_sum[k] += ratio * ratio;
            }
        }
    }
    for k in 0..6 {
        // This is the exact component support formula in real arithmetic;
        // floating-point evaluation is numerical admission, not an enclosure.
        let radius = norm_scale[k] * norm_sum[k].sqrt();
        let maximum = baseline[k].abs() + radius;
        if !maximum.is_finite() || maximum > source.maximum_absolute_mandel_strain() {
            return Err(ElasticityError3::Invalid(
                "affine bulk strain ball exceeds retained material domain",
            ));
        }
    }
    Ok(())
}

impl CutElasticity3 {
    /// Exact `z^T (dK/dscale_c) u`, including each cell's original Nitsche
    /// block and its material-weighted share of each ghost term. No current scale factor is
    /// applied. Direct contraction avoids cancellation in polarization when
    /// primal and adjoint magnitudes differ substantially.
    pub fn scale_bilinear_forms(
        &self,
        z: &[f64],
        u: &[f64],
        mut checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<Vec<f64>, ElasticityError3> {
        self.admit_stress_state(Some(u), &mut checkpoint)?;
        self.admit_stress_state(Some(z), &mut checkpoint)?;
        let mut result = vec![0.0; self.cells.len()];
        for (id, cell) in self.cells.iter().enumerate() {
            poll(&mut checkpoint)?;
            let local = self.stress_local(cell, u);
            let dual = self.stress_local(cell, z);
            for i in 0..24 {
                let applied: f64 = cell.stiffness[i]
                    .iter()
                    .zip(&local)
                    .map(|(k, u)| k * u)
                    .sum();
                result[id] += dual[i] * applied;
            }
        }
        for face in &self.ghosts {
            poll(&mut checkpoint)?;
            let mut value = 0.0;
            for c in 0..3 {
                let jump = |x: &[f64]| -> f64 {
                    face.nodes
                        .iter()
                        .zip(&face.jump)
                        .map(|(&n, &j)| if self.fixed[n] { 0.0 } else { j * x[3 * n + c] })
                        .sum()
                };
                value += jump(z) * jump(u);
            }
            for (&cell, weight) in face.cells.iter().zip(face.scale_weights()) {
                result[cell] += weight * value;
            }
        }
        if !result.iter().all(|v| v.is_finite()) {
            return Err(ElasticityError3::Invalid(
                "bilinear scale contraction overflow",
            ));
        }
        poll(&mut checkpoint)?;
        Ok(result)
    }

    fn stress_point_count(
        &self,
        checkpoint: &mut impl FnMut() -> ControlFlow<()>,
    ) -> Result<usize, ElasticityError3> {
        let mut count = 0usize;
        for cell in &self.cells {
            poll(checkpoint)?;
            count = count
                .checked_add(cell.rules.bulk().len())
                .ok_or(ElasticityError3::Invalid(
                    "bulk stress point count overflow",
                ))?;
        }
        Ok(count)
    }

    // Adaptive reconstruction must admit this output family BEFORE allocating
    // its expanded displacement vectors. Reuse the same point/field accounting.
    pub(super) fn affine_stress_point_count(
        &self,
        basis_count: usize,
        max_point_fields: usize,
        checkpoint: &mut impl FnMut() -> ControlFlow<()>,
    ) -> Result<usize, ElasticityError3> {
        let fields = basis_count.checked_add(1).ok_or(ElasticityError3::Invalid(
            "affine bulk stress field count overflow",
        ))?;
        let count = self.stress_point_count(checkpoint)?;
        if count
            .checked_mul(fields)
            .is_none_or(|n| n > max_point_fields)
        {
            return Err(ElasticityError3::Invalid(
                "affine bulk stress point-field allowance exhausted",
            ));
        }
        Ok(count)
    }

    fn admit_stress_state(
        &self,
        u: Option<&[f64]>,
        checkpoint: &mut impl FnMut() -> ControlFlow<()>,
    ) -> Result<(), ElasticityError3> {
        poll(checkpoint)?;
        if let Some(u) = u {
            admit_displacement(u, self.n(), checkpoint)?;
        }
        for chunk in self.scales.chunks(256) {
            poll(checkpoint)?;
            if chunk
                .iter()
                .any(|v| !v.is_finite() || *v <= 0.0 || *v > 1.0)
            {
                return Err(ElasticityError3::Invalid(
                    "invalid bulk stress material scale",
                ));
            }
        }
        Ok(())
    }

    fn stress_local(&self, cell: &Cell3, u: &[f64]) -> [f64; 24] {
        std::array::from_fn(|i| {
            let node = cell.nodes[i / 3];
            if self.fixed[node] {
                0.0
            } else {
                u[3 * node + i % 3]
            }
        })
    }

    /// Check the retained source-material domain for one physical displacement.
    /// No stress samples or displacement-family buffers are allocated.
    pub fn admit_material_domain(
        &self,
        u: &[f64],
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<(), ElasticityError3> {
        self.admit_affine_material_domain(u, &[], checkpoint)
    }

    /// Admit every field `offset + sum_j a[j] * basis[j]`, `||a||_2 <= 1`,
    /// against the retained source strain ceilings at original bulk points.
    /// Shapes, finite coefficients and scales are checked even for legacy laws;
    /// no ceiling is invented for a law without a retained source domain.
    ///
    /// This streams cell-local fields with constant auxiliary storage and no
    /// heap allocation. It shares the stable component support calculation with
    /// `affine_bulk_stress`, but neither evaluates nor stores constitutive stress.
    /// Cancellation is polled while checking inputs, at cells and points, within
    /// each basis family and before publication. Operator state is unchanged.
    /// The pointwise floating-point admission is not a continuum certificate.
    pub fn admit_affine_material_domain(
        &self,
        offset: &[f64],
        basis: &[&[f64]],
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<(), ElasticityError3> {
        self.admit_affine_material_domain_mapped(
            offset,
            basis,
            self.n(),
            |cell, field| self.stress_local(cell, field),
            checkpoint,
        )
    }

    // The adaptive owner supplies the exact constrained local reconstruction;
    // this native traversal still owns material validity and point ordering.
    pub(super) fn admit_affine_material_domain_mapped(
        &self,
        offset: &[f64],
        basis: &[&[f64]],
        n: usize,
        mut local: impl FnMut(&Cell3, &[f64]) -> [f64; 24],
        mut checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<(), ElasticityError3> {
        self.admit_stress_state(None, &mut checkpoint)?;
        admit_displacement(offset, n, &mut checkpoint)?;
        for &field in basis {
            admit_displacement(field, n, &mut checkpoint)?;
        }
        for cell in &self.cells {
            poll(&mut checkpoint)?;
            if cell.material.source_state().is_none() {
                continue;
            }
            let local_offset = local(cell, offset);
            for &(position, _) in cell.rules.bulk() {
                poll(&mut checkpoint)?;
                let (_, gradients) = q1(cell.bounds, position);
                admit_affine_material_strain(
                    &cell.material,
                    &gradients,
                    &local_offset,
                    basis.iter().map(|field| local(cell, *field)),
                    &mut checkpoint,
                )?;
            }
        }
        poll(&mut checkpoint)?;
        Ok(())
    }

    /// Evaluate the actual bulk quadrature stress with each cell's retained
    /// tensor and current cell scales. Clamped coefficients are masked exactly
    /// as in the operator. `max_points` admits the output allocation before any
    /// point is evaluated. Every point and publication boundary polls the caller.
    /// Source-bound samples must remain inside the admitted Mandel strain
    /// ceiling. No state changes on success, refusal or cancellation.
    pub fn bulk_stress(
        &self,
        u: &[f64],
        max_points: usize,
        mut checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<Vec<BulkStressPoint3>, ElasticityError3> {
        self.admit_stress_state(Some(u), &mut checkpoint)?;
        let count = self.stress_point_count(&mut checkpoint)?;
        if count > max_points {
            return Err(ElasticityError3::Invalid(
                "bulk stress point allowance exhausted",
            ));
        }
        let mut points = Vec::with_capacity(count);
        for (id, cell) in self.cells.iter().enumerate() {
            let local = self.stress_local(cell, u);
            for &(position, weight) in cell.rules.bulk() {
                poll(&mut checkpoint)?;
                let (_, gradients) = q1(cell.bounds, position);
                admit_material_strain(&cell.material, &gradients, &local)?;
                let reference_stress = stress(&cell.material, &gradients, &local);
                let stress = reference_stress.map(|s| self.scales[id] * s);
                if !reference_stress
                    .iter()
                    .chain(&stress)
                    .all(|s| s.is_finite())
                {
                    return Err(ElasticityError3::Invalid("bulk stress overflow"));
                }
                points.push(BulkStressPoint3 {
                    cell: id,
                    position,
                    weight,
                    reference_stress,
                    stress,
                });
            }
        }
        poll(&mut checkpoint)?;
        Ok(points)
    }

    /// Sample one affine field family in a single retained-quadrature traversal.
    /// `max_point_fields` bounds `point_count * (basis.len()+1)` before any
    /// output or local-family allocation. Baseline metadata and every direction
    /// use identical points, masks, tensors and stress arithmetic. Constitutive
    /// basis stresses are prepared once per point and reused for all fields.
    ///
    /// Each source Mandel-strain component must satisfy
    /// `abs(epsilon_offset) + l2(epsilon_basis) <= source_ceiling`. Independent
    /// coordinate-state checks would miss correlated adverse combinations.
    /// Legacy operators retain their existing lack of a source-domain ceiling.
    /// No point family escapes on refusal/cancellation; operator state is unchanged.
    pub fn affine_bulk_stress(
        &self,
        offset: &[f64],
        basis: &[&[f64]],
        max_point_fields: usize,
        mut checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<AffineBulkStress3, ElasticityError3> {
        self.admit_stress_state(Some(offset), &mut checkpoint)?;
        let count =
            self.affine_stress_point_count(basis.len(), max_point_fields, &mut checkpoint)?;
        for &field in basis {
            self.admit_stress_state(Some(field), &mut checkpoint)?;
        }
        let mut result = AffineBulkStress3 {
            baseline: Vec::with_capacity(count),
            basis_reference_stresses: (0..basis.len())
                .map(|_| Vec::with_capacity(count))
                .collect(),
        };
        let mut local_basis = Vec::with_capacity(basis.len());
        for (id, cell) in self.cells.iter().enumerate() {
            poll(&mut checkpoint)?;
            let local_offset = self.stress_local(cell, offset);
            local_basis.clear();
            for &field in basis {
                poll(&mut checkpoint)?;
                local_basis.push(self.stress_local(cell, field));
            }
            for &(position, weight) in cell.rules.bulk() {
                poll(&mut checkpoint)?;
                let (_, gradients) = q1(cell.bounds, position);
                admit_affine_material_strain(
                    &cell.material,
                    &gradients,
                    &local_offset,
                    local_basis.iter().copied(),
                    &mut checkpoint,
                )?;
                let constitutive = stress_basis(&cell.material, &gradients);
                let reference_stress = apply_stress_basis(&constitutive, &local_offset);
                let stress = reference_stress.map(|s| self.scales[id] * s);
                if !reference_stress
                    .iter()
                    .chain(&stress)
                    .all(|s| s.is_finite())
                {
                    return Err(ElasticityError3::Invalid("affine bulk stress overflow"));
                }
                for (column, local) in result.basis_reference_stresses.iter_mut().zip(&local_basis)
                {
                    poll(&mut checkpoint)?;
                    let direction = apply_stress_basis(&constitutive, local);
                    if direction.iter().any(|s| !s.is_finite()) {
                        return Err(ElasticityError3::Invalid("affine bulk stress overflow"));
                    }
                    column.push(direction);
                }
                result.baseline.push(BulkStressPoint3 {
                    cell: id,
                    position,
                    weight,
                    reference_stress,
                    stress,
                });
            }
        }
        poll(&mut checkpoint)?;
        Ok(result)
    }

    /// Exact VJP of the PHYSICAL stresses returned by `bulk_stress`.
    /// `derivatives[i]` is dJ/dsigma_i in the exact retained point order;
    /// quadrature weights are NOT inserted here. Include them in the supplied
    /// cotangents when differentiating a volume integral. The direct scale
    /// derivative holds u fixed. An equilibrium adjoint must separately include
    /// the operator's complete scale derivative, including stabilization.
    pub fn bulk_stress_pullback(
        &self,
        u: &[f64],
        derivatives: &[[f64; 6]],
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<BulkStressPullback3, ElasticityError3> {
        self.stress_pullback(Some(u), derivatives, checkpoint)
    }

    /// Exact VJP of `reference_stress`, without the current stiffness scale.
    /// This is the linear map `(C_reference B)^T`, useful for explicitly
    /// declared stress-relaxation models. It has no direct scale derivative.
    /// Point ordering and unweighted cotangent conventions match the physical
    /// pullback. The reference material is the admitted card, not unit Young's
    /// modulus. Neither entry point differentiates geometry or quadrature.
    pub fn reference_bulk_stress_pullback(
        &self,
        derivatives: &[[f64; 6]],
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<Vec<f64>, ElasticityError3> {
        Ok(self
            .stress_pullback(None, derivatives, checkpoint)?
            .displacement)
    }

    fn stress_pullback(
        &self,
        u: Option<&[f64]>,
        derivatives: &[[f64; 6]],
        mut checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<BulkStressPullback3, ElasticityError3> {
        self.admit_stress_state(u, &mut checkpoint)?;
        if derivatives.len() != self.stress_point_count(&mut checkpoint)? {
            return Err(ElasticityError3::Invalid(
                "one stress cotangent per bulk point required",
            ));
        }
        for d in derivatives {
            poll(&mut checkpoint)?;
            if !d.iter().all(|v| v.is_finite()) {
                return Err(ElasticityError3::Invalid("nonfinite bulk stress cotangent"));
            }
        }
        let mut displacement = vec![0.0; self.n()];
        let mut scales = if u.is_some() {
            vec![0.0; self.cells.len()]
        } else {
            Vec::new()
        };
        let mut index = 0;
        for (id, cell) in self.cells.iter().enumerate() {
            let local = u.map(|u| self.stress_local(cell, u));
            let factor = if u.is_some() { self.scales[id] } else { 1.0 };
            for &(position, _) in cell.rules.bulk() {
                poll(&mut checkpoint)?;
                let (_, gradients) = q1(cell.bounds, position);
                let d = &derivatives[index];
                index += 1;
                if let Some(local) = &local {
                    admit_material_strain(&cell.material, &gradients, local)?;
                    let value = stress(&cell.material, &gradients, local);
                    scales[id] += d.iter().zip(value).map(|(a, b)| a * b).sum::<f64>();
                }
                for i in 0..24 {
                    let node = cell.nodes[i / 3];
                    if self.fixed[node] {
                        continue;
                    }
                    let b = cell.material.basis_stress(gradients[i / 3], i % 3);
                    displacement[3 * node + i % 3] +=
                        factor * d.iter().zip(b).map(|(a, b)| a * b).sum::<f64>();
                }
            }
        }
        if !displacement.iter().chain(&scales).all(|v| v.is_finite()) {
            return Err(ElasticityError3::Invalid("bulk stress pullback overflow"));
        }
        poll(&mut checkpoint)?;
        Ok(BulkStressPullback3 {
            displacement,
            scales,
        })
    }
}
