//! Bulk stress from the retained Q1 field and the assembly's reference tensor.
//!
//! Voigt entries are `[xx, yy, zz, xy, yz, xz]`, with physical shear stresses
//! (no engineering-strain factor in the output). These are numerical samples,
//! not recovered nodal stresses or bounds on a continuum maximum. B-bar uses
//! the same projected dilatation as bulk assembly; its stress is constitutive
//! stress of that assumed strain, not the unprojected displacement strain.
//! Ghost and Nitsche stabilization affect the solved field but are not material stress.
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
    /// `C_reference : eps_formulation(u)` for the original admitted material,
    /// including the cell's B-bar dilatation projection when selected.
    pub reference_stress: [f64; 6],
    /// Current physical stress: `scale[cell] * reference_stress`.
    pub stress: [f64; 6],
}

/// Strain, stress and elastic storage at one retained physical bulk point.
/// Both strains use tensor shear in `[xx, yy, zz, xy, yz, zx]` order.
#[derive(Debug, Clone, PartialEq)]
pub struct BulkStrainStressPoint3 {
    /// Stress, material scaling and the original integration point/weight.
    pub response: BulkStressPoint3,
    /// Actual displacement strain `sym(grad(u))`, before B-bar projection.
    pub strain: [f64; 6],
    /// Strain used by the constitutive law. Equals `strain` for displacement
    /// elements; B-bar replaces only its trace by the physical cell average.
    pub constitutive_strain: [f64; 6],
    /// `sigma : constitutive_strain / 2`, per reference volume. Excludes
    /// ghost/Nitsche stabilization and boundary work; includes current scaling.
    pub energy_density: f64,
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

fn stress_basis(cell: &Cell3, gradients: &[[f64; 3]; 8]) -> [[f64; 6]; 24] {
    let mut basis = std::array::from_fn(|i| {
        if i < 3 {
            [0.0; 6]
        } else {
            cell.basis_stress(gradients[i / 3], i / 3, i % 3)
        }
    });
    // Every formulation differentiates u_a-u_0. Retain its exact algebraic
    // transpose too, so rigid translation cannot create constitutive stress
    // or a spurious material-domain exit through cancellation of large terms.
    for component in 0..3 {
        for k in 0..6 {
            basis[component][k] = -(1..8)
                .map(|corner| basis[3 * corner + component][k])
                .sum::<f64>();
        }
    }
    basis
}

fn apply_stress_basis(basis: &[[f64; 6]; 24], local: &[f64; 24]) -> [f64; 6] {
    let mut value = [0.0; 6];
    for i in 3..24 {
        let displacement = local[i] - local[i % 3];
        for k in 0..6 {
            value[k] += basis[i][k] * displacement;
        }
    }
    value
}

fn stress(cell: &Cell3, gradients: &[[f64; 3]; 8], local: &[f64; 24]) -> [f64; 6] {
    apply_stress_basis(&stress_basis(cell, gradients), local)
}

pub(super) fn mandel_strain(gradients: &[[f64; 3]; 8], local: &[f64; 24]) -> [f64; 6] {
    let mut gradient = [[0.0; 3]; 3];
    for a in 1..8 {
        for c in 0..3 {
            let displacement = local[3 * a + c] - local[c];
            for d in 0..3 {
                gradient[c][d] += displacement * gradients[a][d];
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
    cell: &Cell3,
    gradients: &[[f64; 3]; 8],
    local: &[f64; 24],
) -> Result<(), ElasticityError3> {
    let Some(source) = cell.material.source_state() else {
        return Ok(());
    };
    let strain = mandel_strain(gradients, local);
    let projected = cell.project_mandel_strain(strain, local);
    admit_strains(source.maximum_absolute_mandel_strain(), &strain, &projected)
}

fn admit_strains(
    ceiling: f64,
    strain: &[f64; 6],
    projected: &[f64; 6],
) -> Result<(), ElasticityError3> {
    if strain
        .iter()
        .any(|v| !v.is_finite() || v.abs() > ceiling)
    {
        return Err(ElasticityError3::Invalid(
            "bulk strain exceeds retained material domain",
        ));
    }
    if projected.iter().any(|v| !v.is_finite() || v.abs() > ceiling)
    {
        return Err(ElasticityError3::Invalid(
            "projected bulk strain exceeds retained material domain",
        ));
    }
    Ok(())
}

fn admit_affine_material_strain(
    cell: &Cell3,
    gradients: &[[f64; 3]; 8],
    offset: &[f64; 24],
    basis: impl IntoIterator<Item = [f64; 24]>,
    checkpoint: &mut impl FnMut() -> ControlFlow<()>,
) -> Result<(), ElasticityError3> {
    let Some(source) = cell.material.source_state() else {
        return Ok(());
    };
    // Admit both the actual Q1 kinematics and the assumed constitutive strain.
    // Projection must not conceal a displacement-domain exit, and the projected
    // normal components can themselves exceed the original component ceiling.
    let projected = cell.kinematics.is_projected();
    let family_count = if projected { 2 } else { 1 };
    let baseline = [
        mandel_strain(gradients, offset),
        if projected {
            cell.projected_mandel_strain(gradients, offset)
        } else {
            [0.0; 6]
        },
    ];
    let mut norm_scale = [[0.0_f64; 6]; 2];
    let mut norm_sum = [[0.0_f64; 6]; 2];
    for local in basis {
        poll(checkpoint)?;
        let strains = [
            mandel_strain(gradients, &local),
            if projected {
                cell.projected_mandel_strain(gradients, &local)
            } else {
                [0.0; 6]
            },
        ];
        for family in 0..family_count {
            for k in 0..6 {
                let value = strains[family][k].abs();
                if !value.is_finite() {
                    return Err(ElasticityError3::Invalid("nonfinite affine bulk strain"));
                }
                // Scaled sum of squares avoids squaring a large or tiny direction.
                if value > norm_scale[family][k] {
                    let ratio = norm_scale[family][k] / value;
                    norm_sum[family][k] = 1.0 + norm_sum[family][k] * ratio * ratio;
                    norm_scale[family][k] = value;
                } else if value != 0.0 {
                    let ratio = value / norm_scale[family][k];
                    norm_sum[family][k] += ratio * ratio;
                }
            }
        }
    }
    for family in 0..family_count {
        for k in 0..6 {
            // This is the exact component support formula in real arithmetic;
            // floating-point evaluation is numerical admission, not an enclosure.
            let radius = norm_scale[family][k] * norm_sum[family][k].sqrt();
            let maximum = baseline[family][k].abs() + radius;
            if !maximum.is_finite() || maximum > source.maximum_absolute_mandel_strain() {
                return Err(ElasticityError3::Invalid(
                    "affine bulk strain ball exceeds retained material domain",
                ));
            }
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
    /// against the retained source strain ceilings throughout each active Q1 cell.
    /// Shapes, finite coefficients and scales are checked even for legacy laws;
    /// no ceiling is invented for a law without a retained source domain.
    /// B-bar checks both the actual displacement strain and its projected strain.
    ///
    /// This streams cell-local fields with constant auxiliary storage and no
    /// heap allocation. It shares the stable component support calculation with
    /// `affine_bulk_stress`, but neither evaluates nor stores constitutive stress.
    /// Every Q1 strain component is multiaffine in the three cell coordinates.
    /// Its coefficient-ball support `abs(offset) + norm(basis)` is separately
    /// convex, so checking the eight vertices covers the whole cell in real
    /// arithmetic. B-bar's constant dilatation projection preserves this property.
    /// The check includes the fictitious extension of cut cells and may therefore
    /// conservatively refuse a field whose physical subdomain stays admissible.
    /// Cancellation is polled while checking inputs, at cells and vertices, within
    /// each basis family and before publication. Operator state is unchanged.
    /// This floating-point discrete-field bound is not an interval enclosure or
    /// a continuum-error certificate.
    pub fn admit_affine_material_domain(
        &self,
        offset: &[f64],
        basis: &[&[f64]],
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<(), ElasticityError3> {
        self.admit_affine_material_domain_mapped(
            Some(offset),
            basis,
            self.n(),
            |cell, field| self.stress_local(cell, field),
            checkpoint,
        )
    }

    /// Admit every field `sum_j a[j] * basis[j]`, `||a||_2 <= 1`, against
    /// the retained source strain ceilings. Uses the same streaming whole-ball
    /// owner as affine admission with an implicit zero local offset; no global
    /// zero displacement vector, stress output or other heap storage is allocated.
    pub fn admit_load_space_material_domain(
        &self,
        basis: &[&[f64]],
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<(), ElasticityError3> {
        self.admit_affine_material_domain_mapped(
            None,
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
        offset: Option<&[f64]>,
        basis: &[&[f64]],
        n: usize,
        mut local: impl FnMut(&Cell3, &[f64]) -> [f64; 24],
        mut checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<(), ElasticityError3> {
        self.admit_stress_state(None, &mut checkpoint)?;
        if let Some(offset) = offset {
            admit_displacement(offset, n, &mut checkpoint)?;
        }
        for &field in basis {
            admit_displacement(field, n, &mut checkpoint)?;
        }
        for cell in &self.cells {
            poll(&mut checkpoint)?;
            if cell.material.source_state().is_none() {
                continue;
            }
            let local_offset = offset.map_or([0.0; 24], |field| local(cell, field));
            let (lo, hi) = (cell.bounds.lo(), cell.bounds.hi());
            // Integration nodes alone miss material-domain exits near edges and
            // corners. Multi-affinity reduces the complete Q1-cell bound to these
            // eight vertices, independently of cut-rule size or quadrature depth.
            for corner in 0..8 {
                poll(&mut checkpoint)?;
                let position = std::array::from_fn(|axis| {
                    if (corner >> axis) & 1 == 0 {
                        lo[axis]
                    } else {
                        hi[axis]
                    }
                });
                let (_, gradients) = q1(cell.bounds, position);
                admit_affine_material_strain(
                    cell,
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
                admit_material_strain(cell, &gradients, &local)?;
                let reference_stress = stress(cell, &gradients, &local);
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

    /// Recover actual strain, constitutive strain, stress and elastic storage
    /// together, on the original physical bulk quadrature. Tensor shear is
    /// retained explicitly; B-bar's assumed strain is never mislabeled as the
    /// displacement strain. Both obey the retained source-material domain.
    ///
    /// Only the returned point array is allocated. `max_points` admits that
    /// array before evaluation. Input, cell, point and publication boundaries
    /// poll cancellation. These are numerical samples, not bounds between
    /// quadrature points or a finite-strain constitutive model.
    pub fn bulk_strain_stress(
        &self,
        u: &[f64],
        max_points: usize,
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<Vec<BulkStrainStressPoint3>, ElasticityError3> {
        self.bulk_strain_stress_mapped(
            u,
            self.n(),
            max_points,
            |cell, field| self.stress_local(cell, field),
            checkpoint,
        )
    }

    // Adaptive observations reconstruct only the current cell, rather than
    // expanding a second global displacement vector just to recover strains.
    pub(super) fn bulk_strain_stress_mapped(
        &self,
        u: &[f64],
        n: usize,
        max_points: usize,
        mut local: impl FnMut(&Cell3, &[f64]) -> [f64; 24],
        mut checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<Vec<BulkStrainStressPoint3>, ElasticityError3> {
        self.admit_stress_state(None, &mut checkpoint)?;
        admit_displacement(u, n, &mut checkpoint)?;
        let count = self.stress_point_count(&mut checkpoint)?;
        if count > max_points {
            return Err(ElasticityError3::Invalid(
                "bulk strain/stress point allowance exhausted",
            ));
        }
        let mut points = Vec::with_capacity(count);
        for (id, cell) in self.cells.iter().enumerate() {
            poll(&mut checkpoint)?;
            let local = local(cell, u);
            for &(position, weight) in cell.rules.bulk() {
                poll(&mut checkpoint)?;
                let (_, gradients) = q1(cell.bounds, position);
                let strain = mandel_strain(&gradients, &local);
                let constitutive = cell.project_mandel_strain(strain, &local);
                if let Some(source) = cell.material.source_state() {
                    admit_strains(
                        source.maximum_absolute_mandel_strain(),
                        &strain,
                        &constitutive,
                    )?;
                }
                let strain = fs_material::tensor::from_mandel(&strain);
                let constitutive_strain = fs_material::tensor::from_mandel(&constitutive);
                let reference_stress = stress(cell, &gradients, &local);
                let stress = reference_stress.map(|s| self.scales[id] * s);
                let energy_density =
                    0.5 * fs_material::tensor::contract(&constitutive_strain, &stress);
                if !strain
                    .iter()
                    .chain(&constitutive_strain)
                    .chain(&reference_stress)
                    .chain(&stress)
                    .all(|v| v.is_finite())
                    || !energy_density.is_finite()
                    || energy_density < 0.0
                {
                    return Err(ElasticityError3::Invalid(
                        "unrepresentable bulk strain/stress response",
                    ));
                }
                points.push(BulkStrainStressPoint3 {
                    response: BulkStressPoint3 {
                        cell: id,
                        position,
                        weight,
                        reference_stress,
                        stress,
                    },
                    strain,
                    constitutive_strain,
                    energy_density,
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
                    cell,
                    &gradients,
                    &local_offset,
                    local_basis.iter().copied(),
                    &mut checkpoint,
                )?;
                let constitutive = stress_basis(cell, &gradients);
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
    /// This is the linear map `(C_reference B_formulation)^T`, useful for explicitly
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
                let constitutive = stress_basis(cell, &gradients);
                if let Some(local) = &local {
                    admit_material_strain(cell, &gradients, local)?;
                    let value = apply_stress_basis(&constitutive, local);
                    scales[id] += d.iter().zip(value).map(|(a, b)| a * b).sum::<f64>();
                }
                for i in 0..24 {
                    let node = cell.nodes[i / 3];
                    if self.fixed[node] {
                        continue;
                    }
                    let b = constitutive[i];
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
