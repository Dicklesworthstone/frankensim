//! Static magnetic force and torque from the computed exterior field.
//!
//! Integrate T.n with T = B tensor H - (B.H) I/2 over a CLOSED contour in
//! one homogeneous, current-free linear exterior medium. The integration
//! region encloses the target and an exterior collar; its artificial contour
//! is not the body's mechanical surface. The result is a net wrench, never a
//! guessed local distribution of nodal/contact loads.
//!
//! Reference: J. L. Kirtley, MIT 6.685, Class Notes 1, section 4.2,
//! https://web.mit.edu/6.685/www/chapter1.pdf. The stress includes magnetic
//! contrast inside the contour without assigning a separate force law to each
//! material. Magnetostriction, exterior-medium pressure response, electric
//! stress, moving geometry and end fields are outside this planar model.

use fs_ga::{Vec3, Wrench};
use fs_rep_mesh::TriComplex2;

use crate::{
    planar::{MagneticError, PlanarField, poll},
    transfer::PlanarFrame,
};

/// Bounds on the original mesh traversal and temporary contour storage.
#[derive(Clone, Copy, Debug)]
pub struct MaxwellForceLimits {
    /// Maximum number of original mesh triangles, not only selected cells.
    pub max_cells: usize,
    /// Maximum number of original mesh edges.
    pub max_edges: usize,
}

impl PlanarField<'_> {
    /// Integrate the net magnetic load on a selected region, in SI units.
    ///
    /// `enclosed_cells` names original triangles covering the target plus a
    /// surrounding collar. `TriComplex2::trace_for_faces` owns their boundary:
    /// selected-selected edges cancel and disconnected components/holes retain
    /// their orientation. Every contour edge must have TWO incident cells of
    /// the SAME resolved induced-linear exterior material state, with zero
    /// impressed current and no permanent remanence ON the contour. Enclosed
    /// recoil magnets are supported after the field owner's operating checks.
    /// This refuses contours on material interfaces, current-carrying cells or
    /// the exterior mesh boundary instead of choosing an ambiguous field side.
    /// All contour cells must share that one exterior snapshot. The physical
    /// interpretation assumes a nondeforming exterior without magnetostriction;
    /// a material name or matching permeability alone does not establish it.
    ///
    /// The original exterior-side P1 field is constant along each edge, so
    /// force and its first moment are integrated exactly for THAT discrete
    /// field. Mesh/contour refinement is still required for physical accuracy.
    /// No smoothing across source/material interfaces is performed. Extrusion
    /// thickness is the original metric; invariant end-cap contributions cancel.
    ///
    /// `frame` places the magnetic coordinates in the desired mechanical frame.
    /// Torque is about `moment_origin_m` expressed in that SAME frame. The
    /// returned canonical `fs_ga::Wrench` can feed native articulated dynamics
    /// directly when this frame/origin is the receiving body's frame/origin.
    /// This is an instantaneous static load: displacement requires a new field
    /// solve and a coupled owner for source/field/mechanical work. Do not add an
    /// overlapping Lorentz load to this same enclosed body.
    ///
    /// Fixed mesh-edge order makes the result independent of selection order.
    /// Cancellation, invalid geometry/source correspondence and nonfinite
    /// arithmetic publish no partial wrench and leave the field unchanged.
    pub fn maxwell_wrench(
        &self,
        enclosed_cells: &[usize],
        frame: &PlanarFrame,
        moment_origin_m: [f64; 3],
        limits: MaxwellForceLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Wrench, MagneticError> {
        poll(cancelled)?;
        let mesh = self.problem().mesh();
        let currents = self
            .impressed_current_density_a_m2()
            .ok_or(MagneticError::Invalid(
                "Maxwell force requires an original static current snapshot",
            ))?;
        if enclosed_cells.is_empty()
            || enclosed_cells.len() > mesh.faces().len()
            || mesh.faces().len() > limits.max_cells
            || mesh.edges().len() > limits.max_edges
            || !moment_origin_m.iter().all(|v| v.is_finite())
        {
            return Err(MagneticError::Invalid(
                "invalid magnetic force region, budget or moment origin",
            ));
        }
        let trace = mesh
            .trace_for_faces(enclosed_cells.iter().copied())
            .map_err(|_| {
                MagneticError::Invalid("magnetic force region has repeated or invalid cells")
            })?;
        poll(cancelled)?;
        if trace.edges.is_empty() {
            return Err(MagneticError::Invalid(
                "magnetic force requires a nonempty closed contour",
            ));
        }
        let mut selected = vec![false; mesh.faces().len()];
        for &cell in enclosed_cells {
            selected[cell] = true;
        }
        // Only selected/unselected edges are queried below. For these edges
        // the admitted manifold mesh has exactly one unselected neighbour.
        let mut outside = vec![None; mesh.edges().len()];
        for (cell, face) in mesh.faces().iter().enumerate() {
            poll(cancelled)?;
            if !selected[cell] {
                for (a, b) in [(face[0], face[1]), (face[1], face[2]), (face[2], face[0])] {
                    let edge = mesh
                        .edge_index(a, b)
                        .ok_or(MagneticError::Invalid("magnetic contour edge is absent"))?;
                    outside[edge] = Some(cell);
                }
            }
        }
        let materials = self.problem().materials();
        let recoil = self.problem().recoil_materials();
        let exterior = &materials[trace.edges[0].source_face];
        let mut sum = [0.0; 6];
        let mut correction = [0.0; 6];
        for edge in &trace.edges {
            poll(cancelled)?;
            let inner = edge.source_face;
            let outer = outside[edge.global_edge].ok_or(MagneticError::Invalid(
                "magnetic contour must lie strictly inside the mesh",
            ))?;
            if recoil[inner].is_some() || recoil[outer].is_some() {
                return Err(MagneticError::Invalid(
                    "Maxwell force contour requires an induced-linear exterior without remanence",
                ));
            }
            if materials[inner] != *exterior
                || materials[outer] != *exterior
                || currents[inner] != 0.0
                || currents[outer] != 0.0
            {
                return Err(MagneticError::Invalid(
                    "magnetic contour needs one homogeneous current-free material on both sides",
                ));
            }
            let [a, b] = edge.oriented_vertices.map(|v| mesh.vertices()[v as usize]);
            let orientation = triangle_orientation(mesh, inner)?;
            if triangle_orientation(mesh, outer)? != orientation {
                return Err(MagneticError::Invalid(
                    "magnetic contour neighbour is folded onto its enclosed side",
                ));
            }
            // Clockwise and counterclockwise stored complexes produce the
            // same geometric outward normal. Keep normal*area unnormalized.
            let depth = self.thickness_m() * orientation;
            let normal_area = [depth * (b[1] - a[1]), depth * (a[0] - b[0])];
            if normal_area == [0.0; 2] || !normal_area.iter().all(|v| v.is_finite()) {
                return Err(MagneticError::Invalid(
                    "magnetic contour area is not representable",
                ));
            }
            let flux = self.magnetic_flux_density_t()[outer];
            let h = self.magnetic_field_a_m()[outer];
            let pressure = 0.5 * flux[0].mul_add(h[0], flux[1] * h[1]);
            let h_normal_area = h[0].mul_add(normal_area[0], h[1] * normal_area[1]);
            let local_force = [
                flux[0].mul_add(h_normal_area, -pressure * normal_area[0]),
                flux[1].mul_add(h_normal_area, -pressure * normal_area[1]),
                0.0,
            ];
            let force = frame.vector_to_global(local_force);
            let center = frame.point_to_global([
                0.5 * a[0] + 0.5 * b[0],
                0.5 * a[1] + 0.5 * b[1],
                0.5 * self.thickness_m(),
            ]);
            let arm = [0, 1, 2].map(|i| center[i] - moment_origin_m[i]);
            let torque = [
                arm[1].mul_add(force[2], -arm[2] * force[1]),
                arm[2].mul_add(force[0], -arm[0] * force[2]),
                arm[0].mul_add(force[1], -arm[1] * force[0]),
            ];
            if !force
                .iter()
                .chain(&center)
                .chain(&arm)
                .chain(&torque)
                .all(|v| v.is_finite())
            {
                return Err(MagneticError::Invalid(
                    "magnetic force or moment is not representable",
                ));
            }
            // Opposing contour tractions can be much larger than the net load.
            // Compensate their deterministic sum without modifying the field.
            for (i, value) in torque.into_iter().chain(force).enumerate() {
                let corrected = value - correction[i];
                let next = sum[i] + corrected;
                correction[i] = (next - sum[i]) - corrected;
                sum[i] = next;
                if !sum[i].is_finite() || !correction[i].is_finite() {
                    return Err(MagneticError::Invalid(
                        "integrated magnetic load is not representable",
                    ));
                }
            }
        }
        poll(cancelled)?;
        Ok(Wrench::new(
            Vec3::new(sum[0], sum[1], sum[2]),
            Vec3::new(sum[3], sum[4], sum[5]),
        ))
    }
}

fn triangle_orientation(mesh: &TriComplex2, cell: usize) -> Result<f64, MagneticError> {
    let [p, q, r] = mesh.faces()[cell].map(|v| mesh.vertices()[v as usize]);
    let determinant = (q[0] - p[0]).mul_add(r[1] - p[1], -(q[1] - p[1]) * (r[0] - p[0]));
    if !determinant.is_finite() || determinant == 0.0 {
        return Err(MagneticError::Invalid(
            "magnetic contour orientation is not representable",
        ));
    }
    Ok(determinant.signum())
}
