//! Three-dimensional Maxwell-stress transfer into the native mechanical wrench.
//!
//! A union of original tetrahedra declares a closed control surface surrounding
//! the target and an exterior collar. Only the original exterior-side field is
//! integrated. The stress is `T = B tensor H - (B.H) I/2` in one homogeneous,
//! induced-linear, current-free exterior; no stress is inferred inside a magnet
//! or nonlinear material. Reference: Kirtley, MIT 6.685, Class Notes 1, §4.2,
//! https://web.mit.edu/6.685/www/chapter1.pdf.
//!
//! These are instantaneous forces for the solved fixed geometry. The result
//! does not advance motion or supply a source/field/mechanical energy exchange.
//! Material deformation, magnetostriction and exterior-medium pressure response
//! are excluded. An overlapping Lorentz load must not be added to the same
//! enclosed body. Exact quadrature of the discrete field does not bound its
//! physical error: field and surface refinement remain necessary. This is not
//! a virtual-work comparison, an equilibrated stress, or a Verified force bound.
//! The `eggshell` child adds a prepared distributed virtual-displacement
//! functional with independent coenergy and force/torque accuracy tests.

pub mod eggshell;

use super::{
    Face, SpatialMagneticError, SpatialMagneticField, SpatialMagneticMaterial, cross, dot, finite,
    invalid, poll, sub,
};
use fs_ga::{Se3, Vec3, Wrench};
use fs_rep_mesh::TetComplex;

/// Bounds on the actual prepared mesh, selection and surface traversal.
#[derive(Clone, Copy, Debug)]
pub struct SpatialMaxwellForceLimits {
    pub max_cells: usize,
    pub max_faces: usize,
}

impl SpatialMagneticField<'_> {
    /// Net Maxwell force [N] and torque [N m] on the selected target/collar.
    ///
    /// `enclosed_cells` contains distinct original tetrahedra. Their complete
    /// boundary is derived from the PREPARED face adjacency: selected-selected
    /// faces cancel and every remaining normal points from the selected volume
    /// into its exterior. Holes and disconnected selected regions retain this
    /// orientation. Every surface face must lie strictly INSIDE the mesh, with
    /// two incident cells of the same original induced-linear material snapshot.
    /// All surface cells must share that snapshot; a matching permeability or
    /// material name alone is insufficient. Recoil magnets may be enclosed but
    /// cannot touch the integration surface.
    ///
    /// The caller admits the exterior stress law above (typically air/vacuum),
    /// without material striction or pressure forces. Constant P1 exterior B/H
    /// give constant face traction, so triangle area and centroid integrate
    /// both force and torque exactly for that discrete field, without smoothing
    /// across interfaces or treating its flux jumps as equilibrated.
    ///
    /// `placement` is the canonical SE(3) pose from magnetic mesh coordinates
    /// into the receiving mechanical frame, with translation in metres.
    /// `moment_origin_m` is expressed in that SAME receiving frame. The returned
    /// canonical `fs_ga::Wrench` uses `[torque, force]`, matching native `fs-mbd`
    /// body-wrench inputs when frame and origin are the receiving body's.
    /// Selection order cannot change the fixed canonical face reduction order.
    /// Cancellation or refusal returns no partial wrench and mutates no field.
    pub fn maxwell_wrench(
        &self,
        enclosed_cells: &[usize],
        placement: &Se3,
        moment_origin_m: [f64; 3],
        limits: SpatialMaxwellForceLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Wrench, SpatialMagneticError> {
        let problem = self.problem();
        let mut exterior = None;
        integrate(
            SurfaceGeometry {
                mesh: &problem.mesh,
                positions_m: &problem.positions_m,
                faces: &problem.faces,
            },
            enclosed_cells,
            placement,
            moment_origin_m,
            limits,
            |inner, outer| {
                let (SpatialMagneticMaterial::Linear(a), SpatialMagneticMaterial::Linear(b)) =
                    (&problem.materials[inner], &problem.materials[outer])
                else {
                    return Err(invalid(
                        "3D Maxwell surface requires induced-linear material on both sides",
                    ));
                };
                let source = *exterior.get_or_insert(a);
                if a != source || b != source {
                    return Err(invalid(
                        "3D Maxwell surface requires one homogeneous original exterior source",
                    ));
                }
                Ok((
                    self.flux_density_t()[outer],
                    self.field_strength_a_m()[outer],
                ))
            },
            cancelled,
        )
    }
}

impl super::nonlinear::SpatialNonlinearField<'_> {
    /// Maxwell load on an enclosed nonlinear/recoil target through its linear
    /// exterior collar. Saturated-core fields are retained exactly as solved;
    /// the nonlinear material itself is never assigned a Maxwell stress law.
    ///
    /// Selection, pose, origin, deterministic reduction and refinement limits
    /// have the same contract as `SpatialMagneticField::maxwell_wrench`.
    /// Both neighbors of every surface face must be original `Linear` cells,
    /// all bound to the same complete material snapshot.
    pub fn maxwell_wrench(
        &self,
        enclosed_cells: &[usize],
        placement: &Se3,
        moment_origin_m: [f64; 3],
        limits: SpatialMaxwellForceLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Wrench, SpatialMagneticError> {
        use super::nonlinear::SpatialNonlinearMaterial;
        let problem = self.problem();
        let mut exterior = None;
        integrate(
            SurfaceGeometry {
                mesh: problem.mesh(),
                positions_m: problem.positions_m(),
                faces: &problem.geometry.faces,
            },
            enclosed_cells,
            placement,
            moment_origin_m,
            limits,
            |inner, outer| {
                let (SpatialNonlinearMaterial::Linear(a), SpatialNonlinearMaterial::Linear(b)) =
                    (&problem.materials()[inner], &problem.materials()[outer])
                else {
                    return Err(invalid(
                        "3D Maxwell surface requires induced-linear material on both sides",
                    ));
                };
                let source = *exterior.get_or_insert(a);
                if a != source || b != source {
                    return Err(invalid(
                        "3D Maxwell surface requires one homogeneous original exterior source",
                    ));
                }
                Ok((
                    self.flux_density_t()[outer],
                    self.field_strength_a_m()[outer],
                ))
            },
            cancelled,
        )
    }
}

// A sealed borrowed view of already admitted geometry. Other native 3D field
// owners reuse this integration; they must admit their own original exterior
// source/current snapshots before supplying the OUTER cell's physical B/H.
// No public arbitrary-field constructor can bypass the field/source owner.
pub(super) struct SurfaceGeometry<'a> {
    pub(super) mesh: &'a TetComplex,
    pub(super) positions_m: &'a [[f64; 3]],
    pub(super) faces: &'a [Face],
}

#[allow(clippy::too_many_arguments)]
pub(super) fn integrate(
    geometry: SurfaceGeometry<'_>,
    enclosed_cells: &[usize],
    placement: &Se3,
    moment_origin_m: [f64; 3],
    limits: SpatialMaxwellForceLimits,
    mut fields: impl FnMut(usize, usize) -> Result<([f64; 3], [f64; 3]), SpatialMagneticError>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Wrench, SpatialMagneticError> {
    poll(cancelled)?;
    if enclosed_cells.is_empty()
        || enclosed_cells.len() > geometry.mesh.tets.len()
        || geometry.mesh.tets.len() > limits.max_cells
        || geometry.faces.len() > limits.max_faces
        || moment_origin_m.iter().any(|v| !v.is_finite())
    {
        return Err(invalid(
            "invalid 3D Maxwell region, mesh budget or torque origin",
        ));
    }
    let mut selected = vec![false; geometry.mesh.tets.len()];
    let mut first = usize::MAX;
    for &cell in enclosed_cells {
        poll(cancelled)?;
        if cell >= selected.len() || selected[cell] {
            return Err(invalid(
                "3D Maxwell region has invalid or repeated tetrahedra",
            ));
        }
        selected[cell] = true;
        first = first.min(cell);
    }
    // Moment reduction near the region avoids translating large opposing
    // face tractions separately. The final moment shift acts on the NET force.
    let reference = geometry.positions_m[geometry.mesh.tets[first][0] as usize];
    let mut sum = [0.0; 6];
    let mut correction = [0.0; 6];
    let mut surface_faces = 0;
    for (face_id, face) in geometry.faces.iter().enumerate() {
        poll(cancelled)?;
        let a = face.cells[0];
        let Some(&b) = face.cells.get(1) else {
            if selected[a] {
                return Err(invalid(
                    "3D Maxwell surface must lie strictly inside the original mesh",
                ));
            }
            continue;
        };
        let (inner, outer, sign) = match (selected[a], selected[b]) {
            (true, false) => (a, b, 1.0),
            (false, true) => (b, a, -1.0),
            _ => continue,
        };
        surface_faces += 1;
        let (field, dual) = fields(inner, outer)?;
        if field.iter().chain(&dual).any(|v| !v.is_finite()) {
            return Err(invalid("3D Maxwell exterior field must be finite"));
        }
        let area = face.outward_normal.map(|n| sign * face.area_m2 * n);
        let pressure = finite(0.5 * dot(field, dual))?;
        let normal_dual = finite(dot(dual, area))?;
        let force = std::array::from_fn(|j| field[j].mul_add(normal_dual, -pressure * area[j]));
        let points = geometry.mesh.faces[face_id].map(|v| geometry.positions_m[v as usize]);
        let offset = points.map(|p| sub(p, reference));
        let centroid =
            std::array::from_fn(|j| offset[0][j] / 3.0 + offset[1][j] / 3.0 + offset[2][j] / 3.0);
        let torque = cross(centroid, force);
        for (i, value) in torque.into_iter().chain(force).enumerate() {
            finite(value)?;
            // Opposing Maxwell tractions can greatly exceed their resultant.
            let corrected = finite(value - correction[i])?;
            let next = finite(sum[i] + corrected)?;
            correction[i] = finite((next - sum[i]) - corrected)?;
            sum[i] = next;
        }
    }
    if surface_faces == 0 {
        return Err(invalid(
            "3D Maxwell force requires a nonempty closed surface",
        ));
    }
    let frame_error = |_| invalid("3D Maxwell wrench frame transformation is unrepresentable");
    let force = placement
        .transform_vector(Vec3::new(sum[3], sum[4], sum[5]))
        .map_err(frame_error)?;
    let torque = placement
        .transform_vector(Vec3::new(sum[0], sum[1], sum[2]))
        .map_err(frame_error)?;
    let point = placement
        .transform_point(Vec3::new(reference[0], reference[1], reference[2]))
        .map_err(frame_error)?;
    let arm = Vec3::new(
        point.x - moment_origin_m[0],
        point.y - moment_origin_m[1],
        point.z - moment_origin_m[2],
    );
    let shift = arm.cross(force);
    let result = Wrench::new(
        Vec3::new(torque.x + shift.x, torque.y + shift.y, torque.z + shift.z),
        force,
    );
    for value in result.to_array() {
        finite(value)?;
    }
    poll(cancelled)?;
    Ok(result)
}
