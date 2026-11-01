//! Prepared consistent native-density inertia on the retained 3-D cut geometry.
//!
//! Each cell contributes `occupancy[cell] * integral rho0 N^T N dV`.
//! Scalar node-pair coefficients are shared by all three displacement components.
//! Adaptive preparation uses the existing homogeneous constraint transpose.
//! Stiffness/ghost penalties never manufacture inertia. Positive semidefinite
//! quadrature is preserved, but positive definiteness or a continuum mass bound
//! is not inferred; exact void and strong clamps can leave zero rows.

use std::collections::{BTreeMap, btree_map::Entry};
use std::ops::ControlFlow;
use std::sync::Arc;

use fs_material::elastic::acoustic::ElasticAcousticState;

use super::material::ReferenceMaterial3;
use super::{CutElasticity3, ElasticityError3, q1};

fn poll(checkpoint: &mut impl FnMut() -> ControlFlow<()>) -> Result<(), ElasticityError3> {
    if checkpoint().is_break() {
        Err(ElasticityError3::Cancelled)
    } else {
        Ok(())
    }
}

fn invalid(message: &'static str) -> ElasticityError3 {
    ElasticityError3::Invalid(message)
}

/// Cached physical consistent mass, linear in independent cell occupancy.
///
/// The immutable binding uses the same original geometry/material owner as
/// prepared body loads. Rebuilding an otherwise identical operator requires
/// preparing its own mass. Moving the operator or changing stiffness scales
/// preserves this binding. Geometry coordinates must be in metres and native
/// density in kg/m³; the matrix coefficients then have units of kilograms.
pub struct MaterialMass3 {
    binding: Arc<()>,
    adaptive: bool,
    n: usize,
    fixed: Vec<bool>,
    materials: Vec<ReferenceMaterial3>,
    offsets: Vec<usize>,
    pairs: Vec<[usize; 2]>,
    values: Vec<f64>,
}

impl core::fmt::Debug for MaterialMass3 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("MaterialMass3")
            .field("displacement_dofs", &self.n)
            .field("cells", &self.cells())
            .field("scalar_terms", &self.values.len())
            .finish_non_exhaustive()
    }
}

impl MaterialMass3 {
    /// Displacement dimension in the original operator's ordering.
    #[must_use]
    pub const fn n(&self) -> usize {
        self.n
    }

    /// Number of independent physical cell occupancies.
    #[must_use]
    pub fn cells(&self) -> usize {
        self.materials.len()
    }

    /// Strong homogeneous node clamps; each mask entry covers three components.
    /// Their mass rows are zero, not artificial identity rows.
    #[must_use]
    pub fn fixed(&self) -> &[bool] {
        &self.fixed
    }

    /// Retained upper-triangular scalar node-pair entries, before replication
    /// across the three vector components.
    #[must_use]
    pub fn scalar_terms(&self) -> usize {
        self.values.len()
    }

    /// Original source state; `None` only when the cell index is out of range.
    #[must_use]
    pub fn cell_material_state(&self, cell: usize) -> Option<&ElasticAcousticState> {
        self.materials.get(cell)?.source_state()
    }

    fn admit_occupancy(
        &self,
        occupancy: &[f64],
        checkpoint: &mut impl FnMut() -> ControlFlow<()>,
    ) -> Result<(), ElasticityError3> {
        poll(checkpoint)?;
        if occupancy.len() != self.cells() {
            return Err(invalid("one material-mass occupancy per active cell required"));
        }
        for (cell, &fraction) in occupancy.iter().enumerate() {
            if cell % 64 == 0 {
                poll(checkpoint)?;
            }
            if !fraction.is_finite() || !(0.0..=1.0).contains(&fraction) {
                return Err(invalid("material-mass occupancy must lie in [0,1]"));
            }
        }
        Ok(())
    }

    fn admit_field(
        &self,
        field: &[f64],
        checkpoint: &mut impl FnMut() -> ControlFlow<()>,
    ) -> Result<(), ElasticityError3> {
        if field.len() != self.n {
            return Err(invalid("material-mass displacement shape mismatch"));
        }
        for (dof, value) in field.iter().enumerate() {
            if dof % 192 == 0 {
                poll(checkpoint)?;
            }
            if !value.is_finite() {
                return Err(invalid("nonfinite material-mass displacement"));
            }
        }
        Ok(())
    }

    /// Apply consistent mass using physical fractions in `[0,1]`.
    /// Exact void has zero inertia; the fractions are independent of stiffness.
    pub fn apply(
        &self,
        occupancy: &[f64],
        x: &[f64],
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<Vec<f64>, ElasticityError3> {
        let mut output = vec![0.0; self.n];
        self.apply_into(occupancy, x, &mut output, checkpoint)?;
        Ok(output)
    }

    /// Allocation-free mass action for iterative modal or dynamic solves.
    ///
    /// Shapes and finite inputs are admitted before writing `output`.
    /// An arithmetic failure or cancellation during accumulation may leave
    /// partial output, which the caller must discard.
    pub fn apply_into(
        &self,
        occupancy: &[f64],
        x: &[f64],
        output: &mut [f64],
        mut checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<(), ElasticityError3> {
        self.admit_occupancy(occupancy, &mut checkpoint)?;
        self.admit_field(x, &mut checkpoint)?;
        if output.len() != self.n {
            return Err(invalid("material-mass output shape mismatch"));
        }
        output.fill(0.0);
        for (cell, &fraction) in occupancy.iter().enumerate() {
            poll(&mut checkpoint)?;
            if fraction == 0.0 {
                continue;
            }
            for entry in self.offsets[cell]..self.offsets[cell + 1] {
                if entry % 64 == 0 {
                    poll(&mut checkpoint)?;
                }
                let [a, b] = self.pairs[entry];
                let value = fraction * self.values[entry];
                for component in 0..3 {
                    output[3 * a + component] += value * x[3 * b + component];
                    if a != b {
                        output[3 * b + component] += value * x[3 * a + component];
                    }
                }
            }
        }
        for (dof, value) in output.iter().enumerate() {
            if dof % 192 == 0 {
                poll(&mut checkpoint)?;
            }
            if !value.is_finite() {
                return Err(invalid("material-mass action overflow"));
            }
        }
        poll(&mut checkpoint)
    }

    /// Exact diagonal of the consistent matrix, without row-sum lumping.
    pub fn diagonal(
        &self,
        occupancy: &[f64],
        mut checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<Vec<f64>, ElasticityError3> {
        self.admit_occupancy(occupancy, &mut checkpoint)?;
        let mut diagonal = vec![0.0; self.n];
        for (cell, &fraction) in occupancy.iter().enumerate() {
            poll(&mut checkpoint)?;
            for entry in self.offsets[cell]..self.offsets[cell + 1] {
                if entry % 64 == 0 {
                    poll(&mut checkpoint)?;
                }
                let [a, b] = self.pairs[entry];
                if a == b {
                    let value = fraction * self.values[entry];
                    for component in 0..3 {
                        diagonal[3 * a + component] += value;
                    }
                }
            }
        }
        for (dof, value) in diagonal.iter().enumerate() {
            if dof % 192 == 0 {
                poll(&mut checkpoint)?;
            }
            if !value.is_finite() || *value < 0.0 {
                return Err(invalid("material-mass diagonal overflow"));
            }
        }
        poll(&mut checkpoint)?;
        Ok(diagonal)
    }

    /// Exact derivative of `left^T M(occupancy) right` with respect to physical
    /// cell occupancy. No occupancy, SIMP, filter or projection factor is applied.
    /// Mode vectors are algebraic directions, so no physical strain amplitude
    /// is inferred from their arbitrary normalization.
    pub fn occupancy_pullback(
        &self,
        left: &[f64],
        right: &[f64],
        mut checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<Vec<f64>, ElasticityError3> {
        poll(&mut checkpoint)?;
        self.admit_field(left, &mut checkpoint)?;
        self.admit_field(right, &mut checkpoint)?;
        let mut work = vec![0.0; self.cells()];
        for (cell, total) in work.iter_mut().enumerate() {
            poll(&mut checkpoint)?;
            for entry in self.offsets[cell]..self.offsets[cell + 1] {
                if entry % 64 == 0 {
                    poll(&mut checkpoint)?;
                }
                let [a, b] = self.pairs[entry];
                let value = self.values[entry];
                for component in 0..3 {
                    *total += value * left[3 * a + component] * right[3 * b + component];
                    if a != b {
                        *total += value * left[3 * b + component] * right[3 * a + component];
                    }
                }
            }
            if !total.is_finite() {
                return Err(invalid("material-mass occupancy pullback overflow"));
            }
        }
        poll(&mut checkpoint)?;
        Ok(work)
    }

    pub(super) fn validate_for(
        &self,
        raw: &CutElasticity3,
        adaptive: bool,
        n: usize,
    ) -> Result<(), ElasticityError3> {
        if !Arc::ptr_eq(&self.binding, &raw.material_body_binding)
            || self.adaptive != adaptive
            || self.n != n
            || self.cells() != raw.cells.len()
        {
            return Err(invalid(
                "material mass belongs to another operator or displacement space",
            ));
        }
        Ok(())
    }

    // Integrate once in the native Q1 space, then reduce each cell's symmetric
    // block with the actual homogeneous reconstruction T. No full global matrix
    // or componentwise duplicate of the scalar node-pair storage is required.
    pub(super) fn prepare(
        raw: &CutElasticity3,
        rows: Option<&[Vec<(usize, f64)>]>,
        fixed: &[bool],
        max_terms: usize,
        mut checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<Self, ElasticityError3> {
        poll(&mut checkpoint)?;
        let n = fixed
            .len()
            .checked_mul(3)
            .ok_or(invalid("material-mass displacement size overflow"))?;
        if rows.is_some_and(|r| r.len() != raw.nodes.len())
            || (rows.is_none() && fixed.len() != raw.nodes.len())
        {
            return Err(invalid("material-mass constraint shape mismatch"));
        }
        for cell in &raw.cells {
            poll(&mut checkpoint)?;
            if cell.material.source_state().is_none() {
                return Err(invalid("material mass requires native material density"));
            }
        }
        let allocation = || invalid("material-mass allocation failed");
        let offsets_len = raw
            .cells
            .len()
            .checked_add(1)
            .ok_or(invalid("material-mass cell count overflow"))?;
        let mut result = Self {
            binding: Arc::clone(&raw.material_body_binding),
            adaptive: rows.is_some(),
            n,
            fixed: Vec::new(),
            materials: Vec::new(),
            offsets: Vec::new(),
            pairs: Vec::new(),
            values: Vec::new(),
        };
        result.fixed.try_reserve(fixed.len()).map_err(|_| allocation())?;
        result.fixed.extend_from_slice(fixed);
        result.materials.try_reserve(raw.cells.len()).map_err(|_| allocation())?;
        result.offsets.try_reserve(offsets_len).map_err(|_| allocation())?;
        result.offsets.push(0);
        for cell in &raw.cells {
            poll(&mut checkpoint)?;
            let density = cell.material.source_state().expect("native density admitted").density_kg_m3();
            let mut local = [[0.0; 8]; 8];
            for &(position, weight) in cell.rules.bulk() {
                poll(&mut checkpoint)?;
                let (shape, _) = q1(cell.bounds, position);
                let weight = density * weight;
                if !weight.is_finite() || weight < 0.0 {
                    return Err(invalid("material-mass quadrature weight overflow"));
                }
                for a in 0..8 {
                    for b in a..8 {
                        local[a][b] += weight * shape[a] * shape[b];
                    }
                }
            }
            for a in 0..8 {
                for b in a..8 {
                    if !local[a][b].is_finite() {
                        return Err(invalid("material-mass integration overflow"));
                    }
                    local[b][a] = local[a][b];
                }
            }
            let mut block = BTreeMap::<[usize; 2], f64>::new();
            for a in 0..8 {
                let identity_a = [(cell.nodes[a], 1.0)];
                let row_a = rows.map_or(identity_a.as_slice(), |r| r[cell.nodes[a]].as_slice());
                for b in 0..8 {
                    poll(&mut checkpoint)?;
                    let identity_b = [(cell.nodes[b], 1.0)];
                    let row_b = rows.map_or(identity_b.as_slice(), |r| r[cell.nodes[b]].as_slice());
                    for &(i, wi) in row_a {
                        if i >= fixed.len() || !wi.is_finite() {
                            return Err(invalid("invalid material-mass constraint row"));
                        }
                        if fixed[i] || wi == 0.0 {
                            continue;
                        }
                        for &(j, wj) in row_b {
                            poll(&mut checkpoint)?;
                            if j >= fixed.len() || !wj.is_finite() {
                                return Err(invalid("invalid material-mass constraint row"));
                            }
                            if fixed[j] || wj == 0.0 || i > j {
                                continue;
                            }
                            let value = wi * local[a][b] * wj;
                            if !value.is_finite() {
                                return Err(invalid("constrained material-mass overflow"));
                            }
                            if value == 0.0 {
                                continue;
                            }
                            let potential = result.values.len().checked_add(block.len())
                                .and_then(|n| n.checked_add(1))
                                .ok_or(invalid("material-mass term count overflow"))?;
                            let target = match block.entry([i, j]) {
                                Entry::Occupied(entry) => entry.into_mut(),
                                Entry::Vacant(entry) => {
                                    if potential > max_terms {
                                        return Err(invalid("material-mass term allowance exhausted"));
                                    }
                                    entry.insert(0.0)
                                }
                            };
                            *target += value;
                            if !target.is_finite() {
                                return Err(invalid("constrained material-mass overflow"));
                            }
                        }
                    }
                }
            }
            result.pairs.try_reserve(block.len()).map_err(|_| allocation())?;
            result.values.try_reserve(block.len()).map_err(|_| allocation())?;
            for (pair, value) in block {
                poll(&mut checkpoint)?;
                if value != 0.0 {
                    result.pairs.push(pair);
                    result.values.push(value);
                }
            }
            result.materials.push(cell.material.clone());
            result.offsets.push(result.values.len());
        }
        poll(&mut checkpoint)?;
        Ok(result)
    }
}

impl CutElasticity3 {
    /// Cache native-density consistent inertia on the retained cut quadrature.
    /// `max_terms` bounds upper-triangular scalar node-pair slots, including
    /// the current cell's assembly map, before vector-component replication.
    pub fn prepare_material_mass(
        &self,
        max_terms: usize,
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<MaterialMass3, ElasticityError3> {
        MaterialMass3::prepare(self, None, &self.fixed, max_terms, checkpoint)
    }

    /// Require the exact original geometry/material/displacement-space owner.
    pub fn validate_material_mass(&self, mass: &MaterialMass3) -> Result<(), ElasticityError3> {
        mass.validate_for(self, false, 3 * self.nodes.len())
    }
}
