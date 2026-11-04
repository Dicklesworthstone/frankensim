//! Native electric force through the shared three-dimensional stress integral.

use super::{SpatialDielectricLaw, SpatialElectricField, SpatialElectrostaticError};
use crate::spatial::{
    SpatialMagneticError,
    force::{SpatialMaxwellForceLimits, SurfaceGeometry, integrate},
};
use fs_ga::{Se3, Wrench};

impl SpatialElectricField<'_> {
    /// Integrate `T = D tensor E - (D.E) I/2` over the closed boundary of the
    /// selected original tetrahedra, returning force [N] and torque [N m].
    /// The canonical `placement` maps mesh coordinates to the mechanical frame;
    /// `moment_origin_m` is expressed in that same receiving frame.
    ///
    /// Every surface face lies strictly inside the mesh. Both neighboring
    /// tetrahedra must have the SAME complete original isotropic dielectric
    /// snapshot, zero prescribed volume charge and zero impressed polarization.
    /// That one source must extend
    /// around the whole surface. Conservatively, neither neighboring tetrahedron
    /// may contain an electrode, imposed-displacement or charged-sheet vertex. Charge and
    /// different materials may be enclosed beyond this clear exterior collar.
    ///
    /// The exterior medium is nondeforming, without electrostriction or fluid
    /// pressure forces. Permittivity alone does not establish those mechanical
    /// assumptions. The existing geometric owner integrates the exterior-side
    /// P1 stress, physical areas and centroids in canonical face order, with
    /// compensated force/moment sums and the native SE(3) wrench transformation.
    /// Field and surface refinement remain necessary; this is not a continuum
    /// force bound, equilibrated stress or solid nodal traction distribution.
    ///
    /// This is an instantaneous load on the solved fixed geometry. Advancing
    /// motion needs a new field and explicit electrical/storage/mechanical work
    /// accounting. Do not also add overlapping volume-charge forces to the same
    /// enclosed body. Cancellation or refusal publishes no partial wrench.
    pub fn maxwell_wrench(
        &self,
        enclosed_cells: &[usize],
        placement: &Se3,
        moment_origin_m: [f64; 3],
        limits: SpatialMaxwellForceLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Wrench, SpatialElectrostaticError> {
        let problem = self.problem();
        let geometry = &problem.geometry;
        let mut exterior = None;
        integrate(
            SurfaceGeometry {
                mesh: &geometry.mesh,
                positions_m: &geometry.positions_m,
                faces: &geometry.faces,
            },
            enclosed_cells,
            placement,
            moment_origin_m,
            limits,
            |inner, outer| {
                if !matches!(problem.laws[inner], SpatialDielectricLaw::Isotropic)
                    || !matches!(problem.laws[outer], SpatialDielectricLaw::Isotropic)
                {
                    return Err(SpatialMagneticError::Invalid(
                        "3D electric stress requires an explicitly isotropic exterior dielectric",
                    ));
                }
                let source = *exterior.get_or_insert(&problem.materials[inner]);
                if problem.materials[inner] != *source || problem.materials[outer] != *source {
                    return Err(SpatialMagneticError::Invalid(
                        "3D electric stress needs one homogeneous original exterior source",
                    ));
                }
                if self.free_charge_density_c_m3[inner] != 0.0
                    || self.free_charge_density_c_m3[outer] != 0.0
                    || self.impressed_polarization_c_m2.as_ref().is_some_and(|p| {
                        p[inner] != [0.0; 3] || p[outer] != [0.0; 3]
                    })
                    || geometry.mesh.tets[inner]
                        .iter()
                        .chain(&geometry.mesh.tets[outer])
                        .any(|&v| {
                            self.nodal_boundary_charge_c[v as usize].is_some()
                                || self.charged_sheet_vertices.get(v as usize) == Some(&true)
                        })
                {
                    return Err(SpatialMagneticError::Invalid(
                        "3D electric stress collar must avoid charge, impressed polarization and source vertices",
                    ));
                }
                Ok((
                    self.displacement_c_m2[outer],
                    self.electric_field_v_m[outer],
                ))
            },
            cancelled,
        )
        .map_err(SpatialElectrostaticError::from)
    }
}
