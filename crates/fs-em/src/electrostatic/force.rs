//! Electric stress shares the magnetic owner's geometric contour integral.

use fs_ga::Wrench;
use fs_rep_mesh::Metric2;

use super::{ElectricField, ElectrostaticError, poll};
use crate::{
    force::{
        MaxwellForceLimits,
        contour::{self, ContourError},
    },
    transfer::PlanarFrame,
};

impl ElectricField<'_> {
    /// Integrate the electric Maxwell stress `D tensor E - (D.E) I/2` over a
    /// closed contour enclosing the selected cells and their exterior collar.
    /// The result is a native mechanical `Wrench`, with force in N and torque
    /// in N m about `moment_origin_m` in the supplied frame.
    ///
    /// Both cells touching each contour edge must share ONE complete resolved
    /// linear isotropic dielectric snapshot and have zero original free-charge
    /// source. Contours touching prescribed electrode vertices, dielectric
    /// interfaces, folded cells or the mesh boundary refuse. Enclosed charge
    /// and contrasting dielectrics are allowed: their computed field supplies
    /// the force without a separate material-specific force coefficient.
    ///
    /// The exterior is assumed nondeforming with no electrostriction or fluid
    /// pressure response. A permittivity value alone does not establish these
    /// mechanical assumptions. The artificial contour is not the body's real
    /// surface: this is a NET load, not an inferred solid nodal traction map.
    ///
    /// Original exterior-side P1 fields, actual extrusion depth, orientation,
    /// first moments and compensated sums use the magnetic force owner's same
    /// implementation. Physical accuracy still needs mesh/contour refinement.
    /// End fields are absent from the invariant planar model.
    ///
    /// This is an instantaneous load. Moving the object requires a new field
    /// and explicit source/field/mechanical work accounting; do not reuse frozen
    /// capacitance as an energy-conserving actuator or add an overlapping rho E
    /// force to the same enclosed body. Failure publishes no partial load.
    pub fn maxwell_wrench(
        &self,
        enclosed_cells: &[usize],
        frame: &PlanarFrame,
        moment_origin_m: [f64; 3],
        limits: MaxwellForceLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Wrench, ElectrostaticError> {
        poll(cancelled)?;
        let Metric2::Planar { thickness } = self.problem().mesh().metric() else {
            return Err(ElectrostaticError::Invalid(
                "electric force requires the original planar metric",
            ));
        };
        let materials = self.problem().materials();
        let mut exterior = None;
        contour::integrate(
            self.problem().mesh(),
            thickness,
            enclosed_cells,
            frame,
            moment_origin_m,
            limits,
            |inner, outer, vertices| {
                let source = *exterior.get_or_insert(&materials[inner]);
                if materials[inner] != *source || materials[outer] != *source {
                    return Err(ContourError::Invalid(
                        "electric contour needs one homogeneous dielectric material on both sides",
                    ));
                }
                if self.free_charge_density_c_m3()[inner] != 0.0
                    || self.free_charge_density_c_m3()[outer] != 0.0
                    || vertices
                        .iter()
                        .any(|&v| self.nodal_boundary_charge_c()[v as usize].is_some())
                {
                    return Err(ContourError::Invalid(
                        "electric contour must avoid free-charge sources and electrode vertices",
                    ));
                }
                Ok((
                    self.displacement_c_m2()[outer],
                    self.electric_field_v_m()[outer],
                ))
            },
            cancelled,
        )
        .map_err(|error| match error {
            ContourError::Invalid(reason) => ElectrostaticError::Invalid(reason),
            ContourError::Cancelled => ElectrostaticError::Cancelled,
        })
    }
}
