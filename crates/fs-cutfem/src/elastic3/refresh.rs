//! Transactional material reassembly on the original reference geometry.
//!
//! Only constitutive values change. Retained quadrature, support selection,
//! ghost stencils and displacement coordinates are neither rebuilt nor copied.

use super::*;
use std::mem::size_of;
use std::sync::Arc;

/// Explicit native work/storage bounds before material callbacks or allocation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MaterialRefreshLimits3 {
    /// Maximum original active cells, each queried exactly once.
    pub max_cells: usize,
    /// Sum of retained bulk points, selected embedded-support points, and ghost
    /// points visited by the refresh. No new geometric quadrature is produced.
    pub max_points: usize,
    /// Conservative fixed native staging payload: replacement cell blocks and
    /// kinematics, one fixed-size native tensor cache per cell, and ghost
    /// coefficients. This is NOT a total heap cap. It excludes original dynamic
    /// thermal/source payloads and receipts, caller-owned material resolution/
    /// query allocations, allocator bookkeeping and fixed kernel stack.
    /// Existing geometry is borrowed.
    pub max_scratch_bytes: usize,
}

impl Default for MaterialRefreshLimits3 {
    fn default() -> Self {
        Self {
            max_cells: 32_768,
            max_points: 1_000_000,
            max_scratch_bytes: 256 * 1024 * 1024,
        }
    }
}

struct CellRefresh3 {
    material: ReferenceMaterial3,
    kinematics: CellKinematics3,
    stiffness: [[f64; 24]; 24],
}

fn invalid(reason: &'static str) -> ElasticityError3 {
    ElasticityError3::Invalid(reason)
}

impl CutElasticity3 {
    fn admit_material_refresh(
        &self,
        limits: MaterialRefreshLimits3,
        checkpoint: &mut impl FnMut() -> Result<(), ElasticityError3>,
    ) -> Result<(), ElasticityError3> {
        checkpoint()?;
        if self.cells.len() > limits.max_cells {
            return Err(invalid("material refresh exceeds the cell allowance"));
        }
        if self.reference_cell_masses_kg.is_none() {
            return Err(invalid("material refresh requires original native reference density"));
        }
        let mut points = self.ghosts.len();
        for (index, cell) in self.cells.iter().enumerate() {
            checkpoint()?;
            if cell.material.source_state().is_none() {
                return Err(invalid("material refresh requires original native reference density"));
            }
            let surface_points = self.embedded.as_ref()
                .map_or(0, |support| support.selected_points(index).len());
            points = points.checked_add(cell.rules.bulk().len())
                .and_then(|count| count.checked_add(surface_points))
                .ok_or(invalid("material refresh point count overflow"))?;
            if points > limits.max_points {
                return Err(invalid("material refresh exceeds the retained-point allowance"));
            }
        }
        // Charge a separate fixed native cache even for a uniform provider,
        // which shares one Arc payload. Arc counters and the new owner token
        // are included. Dynamic source payloads and allocator overhead are not.
        let cell_bytes = size_of::<CellRefresh3>()
            .checked_add(size_of::<material::PreparedMaterial3>())
            .and_then(|bytes| bytes.checked_add(2 * size_of::<usize>()))
            .ok_or(invalid("material refresh workspace overflow"))?;
        let bytes = self.cells.len().checked_mul(cell_bytes)
            .and_then(|bytes| self.ghosts.len().checked_mul(size_of::<[f64; 2]>())
                .and_then(|ghosts| bytes.checked_add(ghosts)))
            .and_then(|bytes| bytes.checked_add(2 * size_of::<usize>()))
            .ok_or(invalid("material refresh workspace overflow"))?;
        if bytes > limits.max_scratch_bytes {
            return Err(invalid("material refresh exceeds the staging-byte allowance"));
        }
        checkpoint()
    }

    /// Replace the complete native constitutive field on this SAME reference
    /// geometry. Reassembles bulk, B-bar, selected Nitsche and both ghost sides
    /// from the original retained rules and stencils; no SDF, patch, clamp or
    /// geometric quadrature callback is invoked. Each active-cell material is
    /// resolved once and retained, with no material queries in later applies.
    ///
    /// This is a fixed-reference thermal/material update: every cell's new
    /// native reference-density bits must equal its original density. A current
    /// thermal density cannot change the mass of an undeformed closed specimen.
    /// Legacy states without native inventory refuse. Stiffness-design scales,
    /// physical inventory, geometry, supports and formulation remain unchanged.
    /// B-bar retains its original physical-volume projection and requires the
    /// new source to remain explicitly isotropic; no isotropic fit is inferred.
    ///
    /// All dependent values are staged before publication. Any returned error
    /// or cancellation preserves the old operator and ALL of its preparations.
    /// A successful update rotates ownership: mass, body-load, eigenstrain and
    /// prestress preparations must be recreated even when their dimensions or
    /// reference density match. Source/material domains are never enlarged.
    ///
    /// This commits an OPERATOR update, not a complete heat/mechanical step.
    /// It does not advance temperature, replace a source query, differentiate
    /// C(T), or claim atomicity for a later equilibrium/thermal solve.
    pub fn refresh_material_field(
        &mut self,
        field: ElasticMaterialField3<'_>,
        limits: MaterialRefreshLimits3,
        mut checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<(), ElasticityError3> {
        let mut poll = || {
            if checkpoint().is_break() { Err(ElasticityError3::Cancelled) } else { Ok(()) }
        };
        self.admit_material_refresh(limits, &mut poll)?;
        let provider = MaterialProvider3::from_field(field)?;
        poll()?;
        let mut cells = Vec::with_capacity(self.cells.len());
        for (index, cell) in self.cells.iter().enumerate() {
            poll()?;
            let material = provider.at(cell.bounds)?;
            poll()?;
            let original = cell.material.source_state()
                .expect("original native density was admitted");
            let refreshed = material.source_state()
                .ok_or(invalid("material refresh requires native reference density"))?;
            if original.density_kg_m3().to_bits() != refreshed.density_kg_m3().to_bits() {
                return Err(invalid("material refresh cannot change reference density"));
            }
            let kinematics = cell.kinematics.with_material(&material)?;
            let mut stiffness = assemble_cell_stiffness_controlled(
                cell.bounds, &cell.rules, &material, kinematics, &mut poll)?;
            if let Some(support) = &self.embedded {
                support.reassemble_cell_stiffness(index, cell, &material, kinematics,
                    &mut stiffness, &mut poll)?;
            }
            cells.push(CellRefresh3 { material, kinematics, stiffness });
        }
        let mut ghost_weights = Vec::with_capacity(self.ghosts.len());
        for ghost in &self.ghosts {
            poll()?;
            let weights = ghost.cells.map(|index| {
                let cell = &cells[index];
                ghost.geometry.weight(cell.kinematics.ghost_modulus(&cell.material))
            });
            if weights.iter().any(|weight| !weight.is_finite() || *weight <= 0.0) {
                return Err(invalid("material refresh ghost coefficient overflow"));
            }
            ghost_weights.push(weights);
        }
        let binding = Arc::new(());
        poll()?;
        // The final publication has no callbacks, allocations or fallible work.
        // Geometry, scales and reference inventory retain their exact storage.
        for (cell, replacement) in self.cells.iter_mut().zip(cells) {
            cell.material = replacement.material;
            cell.kinematics = replacement.kinematics;
            cell.stiffness = replacement.stiffness;
        }
        for (ghost, weights) in self.ghosts.iter_mut().zip(ghost_weights) {
            ghost.weights = weights;
        }
        self.material_body_binding = binding;
        Ok(())
    }
}
