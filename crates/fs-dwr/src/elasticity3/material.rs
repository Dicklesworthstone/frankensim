//! Native material-body work on the two actual displacement spaces.
use super::*;
use fs_cutfem::elastic3::MaterialBodyLoad3;
use fs_cutfem::elastic3::adaptive::enrichment::CellResidual3;

/// Prepared native body forces and physical occupancies on both grids.
///
/// Each stencil belongs to its exact operator. Fine physical occupancy must be
/// copied from the coarse parent, independently of the inherited stiffness
/// scale; it must not be refiltered or inferred by inverting SIMP. Source-state
/// inheritance is already checked by `AdaptiveTransfer3`.
///
/// Prepare both stencils from the same pure reference acceleration law. This
/// caller obligation cannot be inferred from two sampled force maps. No nodal
/// prolongation, new quadrature or assumed unit material density is used here.
#[derive(Clone, Copy)]
pub struct MaterialBodyLoads3<'a> {
    pub coarse: &'a MaterialBodyLoad3,
    pub fine: &'a MaterialBodyLoad3,
    pub coarse_occupancy: &'a [f64],
    pub fine_occupancy: &'a [f64],
}

impl MaterialBodyLoads3<'_> {
    pub(super) fn admit(self, transfer: &AdaptiveTransfer3<'_>,
        checkpoint: &mut impl FnMut() -> ControlFlow<()>) -> Result<(), GoalError3> {
        poll(checkpoint)?;
        transfer.coarse().validate_material_body_load(self.coarse)?;
        transfer.fine().validate_material_body_load(self.fine)?;
        if self.coarse_occupancy.len() != self.coarse.cells()
            || self.fine_occupancy.len() != self.fine.cells() {
            return Err(GoalError3::Invalid("one physical body-load occupancy per active cell required"));
        }
        for (i, &rho) in self.coarse_occupancy.iter().enumerate() {
            if i % 64 == 0 { poll(checkpoint)?; }
            if !rho.is_finite() || !(0.0..=1.0).contains(&rho) {
                return Err(GoalError3::Invalid("body-load occupancy must lie in [0,1]"));
            }
        }
        for (i, (&rho, &parent)) in self.fine_occupancy.iter().zip(transfer.parents()).enumerate() {
            if i % 64 == 0 { poll(checkpoint)?; }
            if rho.to_bits() != self.coarse_occupancy[parent].to_bits() {
                return Err(GoalError3::Invalid("enriched physical body-load occupancy is not parent-inherited"));
            }
        }
        poll(checkpoint)
    }

    pub(super) fn add_forces(self, coarse: &mut [f64], fine: &mut [f64],
        checkpoint: &mut impl FnMut() -> ControlFlow<()>) -> Result<(), GoalError3> {
        for (load, occupancy, rhs) in [
            (self.coarse, self.coarse_occupancy, coarse),
            (self.fine, self.fine_occupancy, fine),
        ] {
            let force = load.load(occupancy, &mut *checkpoint)?;
            if force.len() != rhs.len() { return Err(GoalError3::Invalid("material body-force shape mismatch")); }
            for (i, (r, f)) in rhs.iter_mut().zip(force).enumerate() {
                if i % 192 == 0 { poll(checkpoint)?; }
                *r += f;
                if !r.is_finite() { return Err(GoalError3::Invalid("material equilibrium RHS overflow")); }
            }
        }
        poll(checkpoint)
    }

    pub(super) fn add_fine_work(self, rows: &mut [CellResidual3], test: &[f64],
        checkpoint: &mut impl FnMut() -> ControlFlow<()>) -> Result<(), GoalError3> {
        if rows.len() != self.fine_occupancy.len() {
            return Err(GoalError3::Invalid("material body residual shape mismatch"));
        }
        let work = self.fine.occupancy_pullback(test, &mut *checkpoint)?;
        for (i, ((row, value), &rho)) in rows.iter_mut().zip(work).zip(self.fine_occupancy).enumerate() {
            if i % 64 == 0 { poll(checkpoint)?; }
            row.load += rho*value;
            if !row.load.is_finite() || !row.residual().is_finite() {
                return Err(GoalError3::Invalid("material body residual overflow"));
            }
        }
        poll(checkpoint)
    }
}
