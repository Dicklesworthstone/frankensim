//! Complete-revolution magnetic loads from the original meridional field.
//!
//! For poloidal B and azimuthal J, integration over phi cancels every radial
//! force vector and intrinsic torque. The axial surface load is
//! `Fz = 2*pi*integral r*(Bz*Hr*n_r + (Bz*Hz-Br*Hr)*n_z/2) ds`.
//! A meridional edge has oriented normal times length `(dz,-dr)`. Affine B
//! makes the radially weighted stress cubic along that edge: two-point Gauss
//! integration is exact for this discrete field, including sloping contours.
//! The free-current volume load is `Fz = -integral J_phi*B_r dV`.
//!
//! These follow the Lorentz and Maxwell stress descriptions in J. L. Kirtley,
//! MIT 6.685, Class Notes 1, sections 4.1--4.2:
//! <https://web.mit.edu/6.685/www/chapter1.pdf>.
//! Their agreement is a numerical consistency check, not physical validation
//! or an error enclosure. A computed contour field need not be equilibrated;
//! physical loads still require mesh/contour refinement. No magnetic pressure
//! is relabeled as a radial net force or as a local hoop-stress distribution.

use super::{AxisymmetricMagneticField, finite, nonlinear::AxisymmetricNonlinearField, poll};
use crate::{
    force::{
        MaxwellForceLimits,
        contour::{ContourError, triangle_orientation},
    },
    planar::{MagneticError, nonlinear::IsotropicMagneticMaterial},
    transfer::PlanarFrame,
};
use fs_ga::{Vec3, Wrench};
use fs_material::{
    field_storage::linear_magnetic_permeability, state_point::ResolvedMaterialStatePoint,
};
use fs_rep_mesh::{Metric2, TriComplex2};

impl AxisymmetricMagneticField<'_> {
    /// Net magnetic wrench on a complete body of revolution enclosed by cells.
    ///
    /// The original selected-cell trace owns the closed swept surface. Every
    /// positive-area edge must lie strictly inside the mesh, with the SAME
    /// source-resolved linear material and zero impressed current on both sides.
    /// All such edges share one exterior source snapshot. Axis edges have zero
    /// swept area and require no fictitious negative-radius neighbour. The
    /// exterior is assumed mechanically inert, without magnetostriction or a
    /// field-dependent pressure law. This is not inferred from a material name.
    ///
    /// `frame` places the meridional [r,z] plane: its first direction is radial
    /// at phi=0, its SECOND direction is the axis, and its origin lies on that
    /// axis at z=0. `moment_origin_m` is expressed in the same receiving frame.
    /// A full 2*pi sweep is required. Net force is axial and intrinsic torque
    /// vanishes; the returned torque is the moment of that force about the
    /// requested origin. Radial tractions can still produce internal hoop load.
    ///
    /// The native `fs_ga::Wrench` can drive a body using this same frame/origin.
    /// This is an instantaneous static estimate; displacement needs a new field
    /// solve and a coupled work account. Do not add the impressed-current load
    /// for the same enclosed body. No partial load is returned on refusal/cancel.
    pub fn maxwell_wrench(
        &self,
        enclosed_cells: &[usize],
        frame: &PlanarFrame,
        moment_origin_m: [f64; 3],
        limits: MaxwellForceLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Wrench, MagneticError> {
        let mut exterior = None;
        integrate_maxwell(
            self.problem().mesh(),
            self.cell_vertex_flux_density_t(),
            self.impressed_current_density_a_m2(),
            enclosed_cells,
            frame,
            moment_origin_m,
            limits,
            |inner, outer| {
                admit_exterior(
                    &mut exterior,
                    &self.problem().materials()[inner],
                    &self.problem().materials()[outer],
                )
            },
            cancelled,
        )
    }

    /// Integrate ONLY the original impressed-current Lorentz load J_phi cross B.
    ///
    /// `current_cells` names original triangles once each. This can be the total
    /// load on a nonmagnetic coil; it omits magnetization, material-interface,
    /// striction and bound-current forces. A magnetic core or permanent magnet
    /// requires the surrounding Maxwell contour for its complete net load.
    /// The whole accepted field, including self-field, enters the integral.
    ///
    /// The complete-revolution, frame, moment-origin and static-work conventions
    /// are those of `maxwell_wrench`. Exact radial triangle moments integrate
    /// this original cellwise-current/affine-B load without point sampling.
    pub fn impressed_current_wrench(
        &self,
        current_cells: &[usize],
        frame: &PlanarFrame,
        moment_origin_m: [f64; 3],
        limits: MaxwellForceLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Wrench, MagneticError> {
        integrate_current(
            self.problem().mesh(),
            self.cell_vertex_flux_density_t(),
            self.impressed_current_density_a_m2(),
            current_cells,
            frame,
            moment_origin_m,
            limits,
            cancelled,
        )
    }
}

impl AxisymmetricNonlinearField<'_> {
    /// Net load on an enclosed nonlinear core, recoil magnet or current winding.
    ///
    /// This uses the SAME cylindrical stress integral, source-current checks,
    /// full 2*pi sweep and [r,z] frame/origin convention as the linear field's
    /// `maxwell_wrench`. Reversible saturation and admitted fixed recoil can be
    /// anywhere INSIDE the contour; both cells at every positive-area contour
    /// edge must instead share one original current-free induced-linear exterior
    /// source. No nonlinear material stress or magnetostriction law is invented.
    /// Original whole-cell source and recoil admission remains with the field.
    pub fn maxwell_wrench(
        &self,
        enclosed_cells: &[usize],
        frame: &PlanarFrame,
        moment_origin_m: [f64; 3],
        limits: MaxwellForceLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Wrench, MagneticError> {
        let mut exterior = None;
        integrate_maxwell(
            self.problem().mesh(),
            self.cell_vertex_flux_density_t(),
            self.impressed_current_density_a_m2(),
            enclosed_cells,
            frame,
            moment_origin_m,
            limits,
            |inner, outer| {
                let (IsotropicMagneticMaterial::Linear(a), IsotropicMagneticMaterial::Linear(b)) = (
                    &self.problem().materials()[inner],
                    &self.problem().materials()[outer],
                ) else {
                    return Err(MagneticError::Invalid(
                        "axisymmetric stress contour requires an induced-linear exterior",
                    ));
                };
                admit_exterior(&mut exterior, a, b)
            },
            cancelled,
        )
    }

    /// Original free-current J_phi cross B load, with the linear field's exact
    /// moment, full-revolution and frame conventions. This is ONLY the impressed
    /// current contribution; it cannot replace the complete force on a nonlinear
    /// core or permanent magnet. Do not add it to an overlapping Maxwell load.
    pub fn impressed_current_wrench(
        &self,
        current_cells: &[usize],
        frame: &PlanarFrame,
        moment_origin_m: [f64; 3],
        limits: MaxwellForceLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Wrench, MagneticError> {
        integrate_current(
            self.problem().mesh(),
            self.cell_vertex_flux_density_t(),
            self.impressed_current_density_a_m2(),
            current_cells,
            frame,
            moment_origin_m,
            limits,
            cancelled,
        )
    }
}

fn admit_exterior<'a>(
    reference: &mut Option<(&'a ResolvedMaterialStatePoint, f64)>,
    inner: &'a ResolvedMaterialStatePoint,
    outer: &'a ResolvedMaterialStatePoint,
) -> Result<f64, MagneticError> {
    if inner != outer || reference.is_some_and(|(source, _)| source != inner) {
        return Err(MagneticError::Invalid(
            "axisymmetric stress contour requires one original exterior material snapshot",
        ));
    }
    if let Some((_, mu)) = reference {
        return Ok(*mu);
    }
    let mu = linear_magnetic_permeability(inner).map_err(MagneticError::Material)?;
    *reference = Some((inner, mu));
    Ok(mu)
}

fn admission(
    mesh: &TriComplex2,
    selected_cells: &[usize],
    moment_origin_m: [f64; 3],
    limits: MaxwellForceLimits,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<bool>, MagneticError> {
    poll(cancelled)?;
    if mesh.metric()
        != (Metric2::Axisymmetric {
            angular_span: core::f64::consts::TAU,
        })
        || selected_cells.is_empty()
        || selected_cells.len() > mesh.faces().len()
        || mesh.faces().len() > limits.max_cells
        || mesh.edges().len() > limits.max_edges
        || moment_origin_m.iter().any(|v| !v.is_finite())
    {
        return Err(MagneticError::Invalid(
            "finite origin, bounded nonempty region and complete 2*pi sweep required",
        ));
    }
    let mut selected = vec![false; mesh.faces().len()];
    for &cell in selected_cells {
        poll(cancelled)?;
        if cell >= selected.len() || selected[cell] {
            return Err(MagneticError::Invalid(
                "axisymmetric load region has repeated or invalid cells",
            ));
        }
        selected[cell] = true;
    }
    Ok(selected)
}

#[allow(clippy::too_many_arguments)]
fn integrate_maxwell(
    mesh: &TriComplex2,
    flux: &[[[f64; 2]; 3]],
    currents: &[f64],
    enclosed_cells: &[usize],
    frame: &PlanarFrame,
    moment_origin_m: [f64; 3],
    limits: MaxwellForceLimits,
    mut exterior: impl FnMut(usize, usize) -> Result<f64, MagneticError>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Wrench, MagneticError> {
    let selected = admission(mesh, enclosed_cells, moment_origin_m, limits, cancelled)?;
    let trace = mesh
        .trace_for_faces(enclosed_cells.iter().copied())
        .map_err(|_| MagneticError::Invalid("axisymmetric stress needs a closed mesh trace"))?;
    poll(cancelled)?;
    let mut outside = vec![None; mesh.edges().len()];
    for (cell, face) in mesh.faces().iter().enumerate() {
        poll(cancelled)?;
        if !selected[cell] {
            for (a, b) in [(face[0], face[1]), (face[1], face[2]), (face[2], face[0])] {
                let edge = mesh
                    .edge_index(a, b)
                    .ok_or(MagneticError::Invalid("missing axisymmetric contour edge"))?;
                outside[edge] = Some(cell);
            }
        }
    }
    let mut force = 0.0;
    let mut correction = 0.0;
    let mut active = false;
    // Nodes on [0,1]; the two weights are 1/2. The radial surface measure
    // multiplies quadratic stress, so midpoint-only integration is insufficient.
    let offset = 0.5 / 3.0_f64.sqrt();
    for edge in &trace.edges {
        poll(cancelled)?;
        let [a, b] = edge.oriented_vertices.map(|v| mesh.vertices()[v as usize]);
        if a[0] == 0.0 && b[0] == 0.0 {
            continue;
        }
        active = true;
        let inner = edge.source_face;
        let outer = outside[edge.global_edge].ok_or(MagneticError::Invalid(
            "positive-area axisymmetric stress contour must lie inside the mesh",
        ))?;
        let orientation = triangle_orientation(mesh, inner).map_err(contour_error)?;
        if triangle_orientation(mesh, outer).map_err(contour_error)? != orientation {
            return Err(MagneticError::Invalid(
                "axisymmetric stress contour neighbour is folded onto the enclosed side",
            ));
        }
        if currents[inner] != 0.0 || currents[outer] != 0.0 {
            return Err(MagneticError::Invalid(
                "axisymmetric stress contour requires zero current on both sides",
            ));
        }
        let mu = exterior(inner, outer)?;
        let face = mesh.faces()[outer];
        let ends = edge
            .oriented_vertices
            .map(|v| face.iter().position(|&w| v == w).map(|i| flux[outer][i]));
        let [Some(ba), Some(bb)] = ends else {
            return Err(MagneticError::Invalid(
                "axisymmetric stress field lacks its original edge vertices",
            ));
        };
        let normal = [
            finite(orientation * (b[1] - a[1]))?,
            finite(orientation * (a[0] - b[0]))?,
        ];
        for t in [0.5 - offset, 0.5 + offset] {
            let radius = finite((1.0 - t).mul_add(a[0], t * b[0]))?;
            let br = finite((1.0 - t).mul_add(ba[0], t * bb[0]))?;
            let bz = finite((1.0 - t).mul_add(ba[1], t * bb[1]))?;
            let hr = finite(br / mu)?;
            let hz = finite(bz / mu)?;
            let shear = finite(bz * hr)?;
            let axial = finite(0.5 * bz.mul_add(hz, -br * hr))?;
            let traction = finite(shear.mul_add(normal[0], axial * normal[1]))?;
            let weighted = finite((core::f64::consts::PI * radius) * traction)?;
            accumulate(&mut force, &mut correction, weighted)?;
        }
    }
    if !active {
        return Err(MagneticError::Invalid(
            "axisymmetric stress needs a positive-area closed surface",
        ));
    }
    poll(cancelled)?;
    axial_wrench(force, frame, moment_origin_m)
}

#[allow(clippy::too_many_arguments)]
fn integrate_current(
    mesh: &TriComplex2,
    flux: &[[[f64; 2]; 3]],
    currents: &[f64],
    current_cells: &[usize],
    frame: &PlanarFrame,
    moment_origin_m: [f64; 3],
    limits: MaxwellForceLimits,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Wrench, MagneticError> {
    let selected = admission(mesh, current_cells, moment_origin_m, limits, cancelled)?;
    let mut force = 0.0;
    let mut correction = 0.0;
    // Original triangle order makes sums independent of the selection order.
    for (cell, face) in mesh.faces().iter().enumerate() {
        poll(cancelled)?;
        if !selected[cell] || currents[cell] == 0.0 {
            continue;
        }
        let points = face.map(|v| mesh.vertices()[v as usize]);
        let mean_radius = finite(points.iter().map(|p| p[0] / 3.0).sum())?;
        let volume = mesh.face_measure(cell).ok_or(MagneticError::Invalid(
            "missing original swept current volume",
        ))?;
        if mean_radius <= 0.0 || volume <= 0.0 {
            return Err(MagneticError::Invalid(
                "positive swept current volume required",
            ));
        }
        // Exact integral r*B_r dA, with affine r and B_r: the normalized
        // vertex weights are (1 + r_i/sum(r_j))/4, never the centroid weights.
        let mut mean_br = 0.0;
        for i in 0..3 {
            let weight = 0.25 * (1.0 + (points[i][0] / mean_radius) / 3.0);
            mean_br = finite(weight.mul_add(flux[cell][i][0], mean_br))?;
        }
        let value = finite(-currents[cell] * finite(volume * mean_br)?)?;
        accumulate(&mut force, &mut correction, value)?;
    }
    poll(cancelled)?;
    axial_wrench(force, frame, moment_origin_m)
}

fn contour_error(error: ContourError) -> MagneticError {
    match error {
        ContourError::Invalid(reason) => MagneticError::Invalid(reason),
        ContourError::Cancelled => MagneticError::Cancelled,
    }
}

fn accumulate(sum: &mut f64, correction: &mut f64, value: f64) -> Result<(), MagneticError> {
    let adjusted = finite(value - *correction)?;
    let next = finite(*sum + adjusted)?;
    *correction = finite((next - *sum) - adjusted)?;
    *sum = next;
    Ok(())
}

fn axial_wrench(
    axial_force_n: f64,
    frame: &PlanarFrame,
    moment_origin_m: [f64; 3],
) -> Result<Wrench, MagneticError> {
    let force = frame.vector_to_global([0.0, axial_force_n, 0.0]);
    let axis_origin = frame.point_to_global([0.0; 3]);
    let arm = [0, 1, 2].map(|i| axis_origin[i] - moment_origin_m[i]);
    let torque = [
        arm[1].mul_add(force[2], -arm[2] * force[1]),
        arm[2].mul_add(force[0], -arm[0] * force[2]),
        arm[0].mul_add(force[1], -arm[1] * force[0]),
    ];
    if force
        .iter()
        .chain(&arm)
        .chain(&torque)
        .any(|v| !v.is_finite())
    {
        return Err(MagneticError::Invalid(
            "axisymmetric force or moment is not representable",
        ));
    }
    Ok(Wrench::new(
        Vec3::new(torque[0], torque[1], torque[2]),
        Vec3::new(force[0], force[1], force[2]),
    ))
}
