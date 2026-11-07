//! Shared geometric Maxwell-stress integral; constitutive and source admission
//! stays with each field owner. Nonlinear materials supply their own isotropic
//! pressure potential; geometry and first moments have one shared owner.

use fs_ga::{Vec3, Wrench};
use fs_rep_mesh::TriComplex2;

use super::MaxwellForceLimits;
use crate::transfer::PlanarFrame;

pub(crate) enum ContourError {
    Invalid(&'static str),
    Cancelled,
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), ContourError> {
    if cancelled() {
        Err(ContourError::Cancelled)
    } else {
        Ok(())
    }
}

/// Integrate `a tensor b - (a.b) I/2` on the exterior side of each contour
/// edge. The caller must admit both incident cells before supplying that pair.
#[allow(clippy::too_many_arguments)]
pub(crate) fn integrate(
    mesh: &TriComplex2,
    thickness_m: f64,
    enclosed_cells: &[usize],
    frame: &PlanarFrame,
    moment_origin_m: [f64; 3],
    limits: MaxwellForceLimits,
    mut fields: impl FnMut(usize, usize, [u32; 2]) -> Result<([f64; 2], [f64; 2]), ContourError>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Wrench, ContourError> {
    integrate_with_pressure(
        mesh,
        thickness_m,
        enclosed_cells,
        frame,
        moment_origin_m,
        limits,
        |inner, outer, vertices| {
            let (field, dual) = fields(inner, outer, vertices)?;
            let pressure = 0.5 * field[0].mul_add(dual[0], field[1] * dual[1]);
            Ok((field, dual, pressure))
        },
        cancelled,
    )
}

/// Integrate `a tensor b - pressure I` using the original closed geometric
/// contour. The constitutive owner supplies the pressure [J/m3 = Pa] and
/// admits the material/mechanical convention on both incident cells.
#[allow(clippy::too_many_arguments)]
pub(crate) fn integrate_with_pressure(
    mesh: &TriComplex2,
    thickness_m: f64,
    enclosed_cells: &[usize],
    frame: &PlanarFrame,
    moment_origin_m: [f64; 3],
    limits: MaxwellForceLimits,
    mut fields: impl FnMut(usize, usize, [u32; 2]) -> Result<([f64; 2], [f64; 2], f64), ContourError>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Wrench, ContourError> {
    integrate_traction(
        mesh,
        thickness_m,
        enclosed_cells,
        frame,
        moment_origin_m,
        limits,
        |inner, outer, vertices, normal_area| {
            let (field, dual, pressure) = fields(inner, outer, vertices)?;
            if !pressure.is_finite() {
                return Err(ContourError::Invalid(
                    "Maxwell pressure is not representable",
                ));
            }
            let normal_dual = dual[0].mul_add(normal_area[0], dual[1] * normal_area[1]);
            Ok([
                field[0].mul_add(normal_dual, -pressure * normal_area[0]),
                field[1].mul_add(normal_dual, -pressure * normal_area[1]),
            ])
        },
        cancelled,
    )
}

/// Apply the caller's admitted stress directly to the oriented contour area.
/// This permits a material owner to retain stable principal stresses without
/// subtracting a large isotropic pressure from a large dyadic contribution.
#[allow(clippy::too_many_arguments)]
pub(crate) fn integrate_traction(
    mesh: &TriComplex2,
    thickness_m: f64,
    enclosed_cells: &[usize],
    frame: &PlanarFrame,
    moment_origin_m: [f64; 3],
    limits: MaxwellForceLimits,
    mut traction: impl FnMut(usize, usize, [u32; 2], [f64; 2]) -> Result<[f64; 2], ContourError>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Wrench, ContourError> {
    poll(cancelled)?;
    if enclosed_cells.is_empty()
        || enclosed_cells.len() > mesh.faces().len()
        || mesh.faces().len() > limits.max_cells
        || mesh.edges().len() > limits.max_edges
        || !moment_origin_m.iter().all(|v| v.is_finite())
        || !thickness_m.is_finite()
        || thickness_m <= 0.0
    {
        return Err(ContourError::Invalid(
            "invalid Maxwell force region, budget, thickness or moment origin",
        ));
    }
    let trace = mesh
        .trace_for_faces(enclosed_cells.iter().copied())
        .map_err(|_| ContourError::Invalid("Maxwell force region has repeated or invalid cells"))?;
    poll(cancelled)?;
    if trace.edges.is_empty() {
        return Err(ContourError::Invalid(
            "Maxwell force requires a nonempty closed contour",
        ));
    }
    let mut selected = vec![false; mesh.faces().len()];
    for &cell in enclosed_cells {
        selected[cell] = true;
    }
    let mut outside = vec![None; mesh.edges().len()];
    for (cell, face) in mesh.faces().iter().enumerate() {
        poll(cancelled)?;
        if !selected[cell] {
            for (a, b) in [(face[0], face[1]), (face[1], face[2]), (face[2], face[0])] {
                let edge = mesh
                    .edge_index(a, b)
                    .ok_or(ContourError::Invalid("Maxwell contour edge is absent"))?;
                outside[edge] = Some(cell);
            }
        }
    }
    // Reduce about a nearby point in mesh coordinates. Forming every global
    // center first can round away physical lever arms under a large common
    // translation, and shifting opposing tractions separately amplifies error.
    let vertex = mesh.vertices()[trace.edges[0].oriented_vertices[0] as usize];
    let reference = [vertex[0], vertex[1], 0.5 * thickness_m];
    let mut sum = [0.0; 6];
    let mut correction = [0.0; 6];
    for edge in &trace.edges {
        poll(cancelled)?;
        let inner = edge.source_face;
        let outer = outside[edge.global_edge].ok_or(ContourError::Invalid(
            "Maxwell contour must lie strictly inside the mesh",
        ))?;
        let [a, b] = edge.oriented_vertices.map(|v| mesh.vertices()[v as usize]);
        let orientation = triangle_orientation(mesh, inner)?;
        if triangle_orientation(mesh, outer)? != orientation {
            return Err(ContourError::Invalid(
                "Maxwell contour neighbour is folded onto its enclosed side",
            ));
        }
        // The signed triangle orientation makes this geometric normal outward
        // for either coherent mesh winding; leave normal times area unnormalized.
        let depth = thickness_m * orientation;
        let normal_area = [depth * (b[1] - a[1]), depth * (a[0] - b[0])];
        if normal_area == [0.0; 2] || !normal_area.iter().all(|v| v.is_finite()) {
            return Err(ContourError::Invalid(
                "Maxwell contour area is not representable",
            ));
        }
        let local_force = traction(inner, outer, edge.oriented_vertices, normal_area)?;
        let force = [local_force[0], local_force[1], 0.0];
        let arm = Vec3::new(
            0.5 * (a[0] - reference[0]) + 0.5 * (b[0] - reference[0]),
            0.5 * (a[1] - reference[1]) + 0.5 * (b[1] - reference[1]),
            0.0,
        );
        let moment = arm.cross(Vec3::new(force[0], force[1], force[2]));
        let torque = [moment.x, moment.y, moment.z];
        if !force
            .iter()
            .chain(&[arm.x, arm.y, arm.z])
            .chain(&torque)
            .all(|v| v.is_finite())
        {
            return Err(ContourError::Invalid(
                "Maxwell force or moment is not representable",
            ));
        }
        // The individual opposing tractions can greatly exceed their net load.
        for (i, value) in torque.into_iter().chain(force).enumerate() {
            let corrected = value - correction[i];
            let next = sum[i] + corrected;
            correction[i] = (next - sum[i]) - corrected;
            sum[i] = next;
            if !sum[i].is_finite() || !correction[i].is_finite() {
                return Err(ContourError::Invalid(
                    "integrated Maxwell load is not representable",
                ));
            }
        }
    }
    poll(cancelled)?;
    let force = frame.vector_to_global([sum[3], sum[4], sum[5]]);
    let torque = frame.vector_to_global([sum[0], sum[1], sum[2]]);
    let offset = frame.vector_to_global(reference);
    let translation = frame.point_to_global([0.0; 3]);
    // Subtract the receiving origin before adding the small mesh offset.
    let arm = [0, 1, 2].map(|i| (translation[i] - moment_origin_m[i]) + offset[i]);
    let force = Vec3::new(force[0], force[1], force[2]);
    let torque = Vec3::new(torque[0], torque[1], torque[2]);
    let result = Wrench::new(
        torque + Vec3::new(arm[0], arm[1], arm[2]).cross(force),
        force,
    );
    if !result.to_array().iter().all(|v| v.is_finite()) {
        return Err(ContourError::Invalid(
            "integrated Maxwell load is not representable",
        ));
    }
    Ok(result)
}

pub(crate) fn triangle_orientation(mesh: &TriComplex2, cell: usize) -> Result<f64, ContourError> {
    let [p, q, r] = mesh.faces()[cell].map(|v| mesh.vertices()[v as usize]);
    let determinant = (q[0] - p[0]).mul_add(r[1] - p[1], -(q[1] - p[1]) * (r[0] - p[0]));
    if !determinant.is_finite() || determinant == 0.0 {
        return Err(ContourError::Invalid(
            "Maxwell contour orientation is not representable",
        ));
    }
    Ok(determinant.signum())
}
