//! Prepared body forces from native reference density and a fixed acceleration.
//!
//! Each cell contributes a sparse force column
//! `B_c = integral rho0_c N^T acceleration dV`. Physical occupancy gives `B r`;
//! its exact transpose gives `B^T z` for design sensitivities. Native cut rules,
//! material states and displacement constraints are consumed once in preparation.
//! Neither stiffness interpolation nor an assumed unit density enters this map.
//! Coordinates are in metres, acceleration in m/s² and assembled forces in N.

use std::collections::{BTreeMap, btree_map::Entry};
use std::ops::ControlFlow;
use std::sync::Arc;

use fs_material::elastic::acoustic::ElasticAcousticState;

use super::material::ReferenceMaterial3;
use super::{CutElasticity3, ElasticityError3, q1};

#[derive(Clone, Copy, PartialEq, Eq)]
pub(super) enum BodyLoadSpace3 {
    Cartesian,
    Adaptive,
}

/// Geometry-prepared acceleration load in one operator's displacement space.
///
/// Columns follow active-cell order and already eliminate strong homogeneous
/// clamps. Adaptive columns include the original hanging-node transpose. All
/// actions are sparse, allocate only their output, and never call acceleration
/// or regenerate quadrature. Source states remain attached to their columns.
///
/// An immutable owner token binds geometry, sources, clamps and displacement
/// ordering. Use the operator's `validate_material_body_load` before attaching
/// this value to a study. Independently rebuilt operators require preparation
/// on their own geometry, even when their dimensions match. This in-memory
/// ownership identity is not a content digest or scientific provenance claim.
pub struct MaterialBodyLoad3 {
    binding: Arc<()>,
    space: BodyLoadSpace3,
    n: usize,
    materials: Vec<ReferenceMaterial3>,
    offsets: Vec<usize>,
    nodes: Vec<usize>,
    forces: Vec<[f64; 3]>,
}

impl core::fmt::Debug for MaterialBodyLoad3 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("MaterialBodyLoad3")
            .field("displacement_dofs", &self.n)
            .field("cells", &self.cells())
            .field("scalar_terms", &(3 * self.nodes.len()))
            .finish_non_exhaustive()
    }
}

fn poll(checkpoint: &mut impl FnMut() -> ControlFlow<()>) -> Result<(), ElasticityError3> {
    if checkpoint().is_break() {
        Err(ElasticityError3::Cancelled)
    } else {
        Ok(())
    }
}

fn allocation() -> ElasticityError3 {
    ElasticityError3::Invalid("material body-load allocation failed")
}

impl MaterialBodyLoad3 {
    /// Number of displacement entries in the original independent ordering.
    #[must_use]
    pub const fn n(&self) -> usize {
        self.n
    }

    /// Number of physical occupancy entries, in original active-cell order.
    #[must_use]
    pub fn cells(&self) -> usize {
        self.materials.len()
    }

    /// Original immutable source state; `None` only for an invalid cell index.
    #[must_use]
    pub fn cell_material_state(&self, cell: usize) -> Option<&ElasticAcousticState> {
        self.materials.get(cell)?.source_state()
    }

    /// Apply `B occupancy` in the retained displacement space.
    ///
    /// Occupancy is the physical remaining fraction in `[0,1]`, after any
    /// filter/projection and protected-region overrides. It is independent of
    /// SIMP stiffness scales: exact void has no material body force.
    pub fn load(
        &self,
        occupancy: &[f64],
        mut checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<Vec<f64>, ElasticityError3> {
        poll(&mut checkpoint)?;
        if occupancy.len() != self.cells() {
            return Err(ElasticityError3::Invalid(
                "one body-load occupancy per active cell required",
            ));
        }
        for (cell, &fraction) in occupancy.iter().enumerate() {
            if cell % 64 == 0 {
                poll(&mut checkpoint)?;
            }
            if !fraction.is_finite() || !(0.0..=1.0).contains(&fraction) {
                return Err(ElasticityError3::Invalid(
                    "material body-load occupancy must lie in [0,1]",
                ));
            }
        }
        let mut rhs = vec![0.0; self.n];
        for (cell, &fraction) in occupancy.iter().enumerate() {
            poll(&mut checkpoint)?;
            for entry in self.offsets[cell]..self.offsets[cell + 1] {
                if entry % 64 == 0 {
                    poll(&mut checkpoint)?;
                }
                let node = self.nodes[entry];
                for component in 0..3 {
                    rhs[3 * node + component] += fraction * self.forces[entry][component];
                }
            }
        }
        for (dof, value) in rhs.iter().enumerate() {
            if dof % 192 == 0 {
                poll(&mut checkpoint)?;
            }
            if !value.is_finite() {
                return Err(ElasticityError3::Invalid(
                    "material body-load action overflow",
                ));
            }
        }
        poll(&mut checkpoint)?;
        Ok(rhs)
    }

    /// Exact `B^T dual`: derivative of `dual^T load(occupancy)` with respect to
    /// physical cell occupancy. Filter/projection derivatives belong outside;
    /// no occupancy or stiffness multiplier is applied a second time.
    pub fn occupancy_pullback(
        &self,
        dual: &[f64],
        mut checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<Vec<f64>, ElasticityError3> {
        poll(&mut checkpoint)?;
        if dual.len() != self.n {
            return Err(ElasticityError3::Invalid(
                "material body-load dual shape mismatch",
            ));
        }
        for (dof, value) in dual.iter().enumerate() {
            if dof % 192 == 0 {
                poll(&mut checkpoint)?;
            }
            if !value.is_finite() {
                return Err(ElasticityError3::Invalid(
                    "nonfinite material body-load dual",
                ));
            }
        }
        let mut work = vec![0.0; self.cells()];
        for (cell, value) in work.iter_mut().enumerate() {
            poll(&mut checkpoint)?;
            for entry in self.offsets[cell]..self.offsets[cell + 1] {
                if entry % 64 == 0 {
                    poll(&mut checkpoint)?;
                }
                let node = self.nodes[entry];
                for component in 0..3 {
                    *value += self.forces[entry][component] * dual[3 * node + component];
                }
            }
            if !value.is_finite() {
                return Err(ElasticityError3::Invalid(
                    "material body-load pullback overflow",
                ));
            }
        }
        poll(&mut checkpoint)?;
        Ok(work)
    }

    pub(super) fn validate_for(
        &self,
        raw: &CutElasticity3,
        space: BodyLoadSpace3,
        n: usize,
    ) -> Result<(), ElasticityError3> {
        if !Arc::ptr_eq(&self.binding, &raw.material_body_binding)
            || self.space != space
            || self.n != n
            || self.cells() != raw.cells.len()
        {
            return Err(ElasticityError3::Invalid(
                "material body load belongs to another operator or displacement space",
            ));
        }
        Ok(())
    }

    // One geometry/material integrator and one sparse transpose construction.
    // None means identity raw-node rows; adaptive callers pass their exact T.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn prepare(
        raw: &CutElasticity3,
        rows: Option<&[Vec<(usize, f64)>]>,
        fixed: &[bool],
        space: BodyLoadSpace3,
        acceleration: &dyn Fn([f64; 3]) -> [f64; 3],
        max_terms: usize,
        mut checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<Self, ElasticityError3> {
        poll(&mut checkpoint)?;
        let n = fixed.len().checked_mul(3).ok_or(ElasticityError3::Invalid(
            "material body-load displacement size overflow",
        ))?;
        if rows.is_some_and(|r| r.len() != raw.nodes.len())
            || (rows.is_none() && fixed.len() != raw.nodes.len())
        {
            return Err(ElasticityError3::Invalid(
                "material body-load constraint shape mismatch",
            ));
        }
        // Refuse missing density before invoking any caller acceleration law.
        for cell in &raw.cells {
            poll(&mut checkpoint)?;
            if cell.material.source_state().is_none() {
                return Err(ElasticityError3::Invalid(
                    "material body load requires native material density",
                ));
            }
        }
        let offset_count = raw
            .cells
            .len()
            .checked_add(1)
            .ok_or(ElasticityError3::Invalid(
                "material body-load cell count overflow",
            ))?;
        let mut result = Self {
            binding: Arc::clone(&raw.material_body_binding),
            space,
            n,
            materials: Vec::new(),
            offsets: Vec::new(),
            nodes: Vec::new(),
            forces: Vec::new(),
        };
        result
            .materials
            .try_reserve(raw.cells.len())
            .map_err(|_| allocation())?;
        result
            .offsets
            .try_reserve(offset_count)
            .map_err(|_| allocation())?;
        result.offsets.push(0);
        for cell in &raw.cells {
            poll(&mut checkpoint)?;
            let density = cell
                .material
                .source_state()
                .expect("native density admitted")
                .density_kg_m3();
            let mut local = [[0.0; 3]; 8];
            for &(position, weight) in cell.rules.bulk() {
                poll(&mut checkpoint)?;
                let acceleration = acceleration(position);
                poll(&mut checkpoint)?;
                if acceleration.iter().any(|a| !a.is_finite()) {
                    return Err(ElasticityError3::Invalid("nonfinite material acceleration"));
                }
                let force = acceleration.map(|a| density * a);
                if force.iter().any(|f| !f.is_finite()) {
                    return Err(ElasticityError3::Invalid(
                        "material body-force density overflow",
                    ));
                }
                let (shape, _) = q1(cell.bounds, position);
                for node in 0..8 {
                    for component in 0..3 {
                        local[node][component] += weight * shape[node] * force[component];
                    }
                }
            }
            if local.iter().flatten().any(|v| !v.is_finite()) {
                return Err(ElasticityError3::Invalid(
                    "material body-load integration overflow",
                ));
            }
            let mut column = BTreeMap::<usize, [f64; 3]>::new();
            for (node, force) in cell.nodes.iter().zip(local) {
                poll(&mut checkpoint)?;
                if force == [0.0; 3] {
                    continue;
                }
                let identity = [(*node, 1.0)];
                let row = rows.map_or(identity.as_slice(), |r| r[*node].as_slice());
                for &(master, factor) in row {
                    poll(&mut checkpoint)?;
                    if master >= fixed.len() || !factor.is_finite() {
                        return Err(ElasticityError3::Invalid(
                            "invalid material body-load constraint row",
                        ));
                    }
                    if fixed[master] || factor == 0.0 {
                        continue;
                    }
                    let potential_nodes = result
                        .nodes
                        .len()
                        .checked_add(column.len())
                        .and_then(|n| n.checked_add(1))
                        .ok_or(ElasticityError3::Invalid(
                            "material body-load term count overflow",
                        ))?;
                    let target = match column.entry(master) {
                        Entry::Occupied(entry) => entry.into_mut(),
                        Entry::Vacant(entry) => {
                            if potential_nodes.checked_mul(3).is_none_or(|n| n > max_terms) {
                                return Err(ElasticityError3::Invalid(
                                    "material body-load term allowance exhausted",
                                ));
                            }
                            entry.insert([0.0; 3])
                        }
                    };
                    for component in 0..3 {
                        target[component] += factor * force[component];
                        if !target[component].is_finite() {
                            return Err(ElasticityError3::Invalid(
                                "constrained material body-load overflow",
                            ));
                        }
                    }
                }
            }
            // The allowance includes the current column's assembly slots;
            // check it before any map insertion or retained-vector growth.
            result
                .nodes
                .try_reserve(column.len())
                .map_err(|_| allocation())?;
            result
                .forces
                .try_reserve(column.len())
                .map_err(|_| allocation())?;
            for (node, force) in column {
                poll(&mut checkpoint)?;
                if force != [0.0; 3] {
                    result.nodes.push(node);
                    result.forces.push(force);
                }
            }
            result.materials.push(cell.material.clone());
            result.offsets.push(result.nodes.len());
        }
        poll(&mut checkpoint)?;
        Ok(result)
    }
}

impl CutElasticity3 {
    /// Prepare native-density self-weight or another fixed reference acceleration.
    /// The pure acceleration callback is evaluated at retained reference points
    /// only during preparation; later actions retain the same environmental law.
    /// `max_terms` bounds retained scalar slots (three per sparse node record),
    /// including the current cell's assembly map. Missing density, exhausted
    /// allowance, nonfinite arithmetic or cancellation publishes no stencil.
    pub fn prepare_material_body_load(
        &self,
        acceleration: &dyn Fn([f64; 3]) -> [f64; 3],
        max_terms: usize,
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<MaterialBodyLoad3, ElasticityError3> {
        MaterialBodyLoad3::prepare(
            self,
            None,
            &self.fixed,
            BodyLoadSpace3::Cartesian,
            acceleration,
            max_terms,
            checkpoint,
        )
    }

    /// Require the original immutable geometry/material/displacement owner.
    /// Changing stiffness scales or moving the operator preserves this binding.
    pub fn validate_material_body_load(
        &self,
        load: &MaterialBodyLoad3,
    ) -> Result<(), ElasticityError3> {
        load.validate_for(self, BodyLoadSpace3::Cartesian, 3 * self.nodes.len())
    }
}
