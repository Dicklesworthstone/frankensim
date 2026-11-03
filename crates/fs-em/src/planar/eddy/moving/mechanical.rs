//! Planar conduction-current forces with exact P1 torque and virtual work.
//!
//! Reuses the spatial induction load integrator with three in-plane shape
//! functions. Each uniform extrusion is represented at its actual mid-depth,
//! so arbitrary receiving-frame wrenches retain the full force lever arm.
//! These volume forces exclude magnetization and interface traction; do not
//! also apply an overlapping enclosing Maxwell-stress load to the same body.
//! The original fixed-occupancy induction scope and prescribed motion remain.

use super::{MovingEddyCandidate, MovingEddyCell};
use crate::planar::eddy::EddySection;
use crate::planar::eddy::nonlinear::moving::NonlinearMovingEddyCandidate;
use crate::spatial::vector_potential::induction::mechanical::load::{CellLoad, integrate_mapped};
use fs_rep_mesh::TriComplex2;

pub use crate::spatial::vector_potential::induction::mechanical::{
    SpatialLorentzLoadError as PlanarLorentzLoadError,
    SpatialLorentzLoadLimits as PlanarLorentzLoadLimits, SpatialLorentzLoads as PlanarLorentzLoads,
};

impl MovingEddyCandidate<'_> {
    /// Extract the SAME midpoint current/field stage as the induction work.
    ///
    /// Selection names distinct original conducting triangles, in any order.
    /// Nodal forces use original vertex order and planar x/y axes, embedded
    /// at z = thickness/2 in the declared 0..thickness extrusion. Pairing them
    /// with [vx, vy, 0] recovers the selected native mechanical power exactly
    /// up to floating-point reduction. Use `wrench` for receiving-frame and
    /// moment-origin conversion; multiply by this step's duration for impulse.
    ///
    /// This transfer does not advance the field or mechanical state. A coupled
    /// receiver must use the SAME stage and close its actual mechanical work.
    pub fn lorentz_loads(
        &self,
        selected_triangles: &[usize],
        limits: PlanarLorentzLoadLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PlanarLorentzLoads, PlanarLorentzLoadError> {
        let problem = self.state().problem();
        recover(
            problem.mesh(),
            problem.thickness_m(),
            &self.step().cells,
            selected_triangles,
            limits,
            cancelled,
        )
    }
}

impl NonlinearMovingEddyCandidate<'_> {
    /// Extract the accepted ENDPOINT current/field stage, including admitted
    /// saturation, fixed recoil and explicit solid-current drives.
    ///
    /// Geometry and selection follow [`MovingEddyCandidate::lorentz_loads`].
    /// This endpoint load pairs with the nonlinear update's own duration and
    /// prescribed velocity. Backward-Euler numerical dissipation remains
    /// separate from both Lorentz work and physical Joule heat.
    pub fn lorentz_loads(
        &self,
        selected_triangles: &[usize],
        limits: PlanarLorentzLoadLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PlanarLorentzLoads, PlanarLorentzLoadError> {
        let problem = self.state().problem();
        recover(
            problem.mesh(),
            problem.thickness_m(),
            &self.step().cells,
            selected_triangles,
            limits,
            cancelled,
        )
    }
}

fn recover(
    mesh: &TriComplex2,
    thickness_m: f64,
    cells: &[MovingEddyCell],
    selected_triangles: &[usize],
    limits: PlanarLorentzLoadLimits,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<PlanarLorentzLoads, PlanarLorentzLoadError> {
    use PlanarLorentzLoadError::{Cancelled, Invalid};
    if cancelled() {
        return Err(Cancelled);
    }
    if selected_triangles.is_empty()
        || selected_triangles.len() > limits.max_cells
        || mesh.vertices().len() > limits.max_vertices
    {
        return Err(Invalid("empty or over-budget planar Lorentz selection"));
    }
    // Native conductor binding sorts records by original triangle. Preserve
    // sparse body selection without a full conductor scan per mechanical body.
    let mut selected = Vec::with_capacity(selected_triangles.len());
    for &triangle in selected_triangles {
        if cancelled() {
            return Err(Cancelled);
        }
        let record = cells
            .binary_search_by_key(&triangle, |record| record.joule.source.triangle)
            .ok()
            .and_then(|index| cells.get(index))
            .ok_or(Invalid(
                "Lorentz selection requires original conducting triangles",
            ))?;
        let vertices = *mesh
            .faces()
            .get(triangle)
            .ok_or(Invalid("Lorentz triangle is outside its original mesh"))?;
        selected.push((
            triangle,
            vertices,
            record.joule.volume_m3,
            record.lorentz_force_density_n_m3.map(|f| [f[0], f[1], 0.0]),
        ));
    }
    selected.sort_unstable_by_key(|entry| entry.0);
    if selected.windows(2).any(|pair| pair[0].0 == pair[1].0) {
        return Err(Invalid("Lorentz selection contains a repeated triangle"));
    }
    integrate_mapped(
        mesh.vertices().len(),
        |node| {
            mesh.vertices()
                .get(node)
                .map(|point| [point[0], point[1], 0.5 * thickness_m])
        },
        selected
            .iter()
            .map(|(_, vertices, volume, density)| CellLoad {
                vertices: *vertices,
                volume_m3: *volume,
                force_density_n_m3: density,
            }),
        limits,
        cancelled,
    )
}
