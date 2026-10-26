//! Electromagnetic force and torque from computed exterior fields.
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
//! material. The electric owner uses the same contour, normals, first moments
//! and compensated reduction with its own D/E source admission. Striction,
//! exterior-medium pressure response, moving geometry and end fields remain
//! outside this planar model.

use fs_ga::Wrench;

pub(crate) mod contour;
mod induction;
use contour::ContourError;

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
    /// A current derivative on a magnetized problem is not an actual operating
    /// field and refuses; unbiased current responses remain usable.
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
        if self.is_incremental_response() && self.problem().has_permanent_magnets() {
            return Err(MagneticError::Invalid(
                "Maxwell force requires the actual operating field including permanent-magnet bias",
            ));
        }
        let currents = self
            .impressed_current_density_a_m2()
            .ok_or(MagneticError::Invalid(
                "Maxwell force requires an original static current snapshot",
            ))?;
        let materials = self.problem().materials();
        let recoil = self.problem().recoil_materials();
        let mut exterior = None;
        contour::integrate(
            self.problem().mesh(),
            self.thickness_m(),
            enclosed_cells,
            frame,
            moment_origin_m,
            limits,
            |inner, outer, _| {
                if recoil[inner].is_some() || recoil[outer].is_some() {
                    return Err(ContourError::Invalid(
                        "Maxwell force contour requires an induced-linear exterior without remanence",
                    ));
                }
                let source = *exterior.get_or_insert(&materials[inner]);
                if materials[inner] != *source
                    || materials[outer] != *source
                    || currents[inner] != 0.0
                    || currents[outer] != 0.0
                {
                    return Err(ContourError::Invalid(
                        "magnetic contour needs one homogeneous current-free material on both sides",
                    ));
                }
                Ok((self.magnetic_flux_density_t()[outer], self.magnetic_field_a_m()[outer]))
            },
            cancelled,
        ).map_err(|error| match error {
            ContourError::Invalid(reason) => MagneticError::Invalid(reason),
            ContourError::Cancelled => MagneticError::Cancelled,
        })
    }
}
