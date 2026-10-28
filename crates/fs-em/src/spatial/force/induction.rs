//! Mechanical loads from the actual accepted A–phi induction history.
//!
//! Both reversible-linear and saturating transient owners reuse the original
//! closed-surface and prepared P1 eggshell functionals. B and H are recovered
//! from retained A and its actual constitutive source; no static re-solve may
//! erase eddy-current memory after an impressed source is switched off.
//!
//! Every integration cell must be an original homogeneous induced-linear
//! exterior without a conductor or impressed current. Conducting, saturating
//! and recoil targets may be enclosed. The conservative no-conductor collar
//! admission does not mistake a momentarily zero current for insulating air.
//! These are instantaneous magnetic wrenches, not an update of mechanical
//! state or a debit of mechanical work. Do not add overlapping Lorentz loads.
//! Displacement-current, moving-geometry and exterior-medium striction effects
//! remain outside the underlying magnetoquasistatic model. Mesh/collar
//! refinement is still required; exact discrete integration is not an error bound.

use super::super::{
    SpatialGeometry, SpatialMagneticError, SpatialMagneticMaterial, invalid,
    nonlinear::SpatialNonlinearMaterial,
    poll,
    vector_potential::induction::{
        SpatialInduction, SpatialInductionField,
        nonlinear::{SpatialNonlinearInduction, SpatialNonlinearInductionError},
    },
};
use super::{SpatialMaxwellForceLimits, SurfaceGeometry, eggshell, integrate};
use eggshell::PreparedSpatialEggshell;
use fs_ga::{Se3, Wrench};
use fs_material::{conductor::ResolvedConductorState, state_point::ResolvedMaterialStatePoint};

// Private: arbitrary detached B/H samples cannot implement the public load API.
// The same integration/admission code serves both original transient owners.
trait Owner {
    type Error: From<SpatialMagneticError>;
    fn geometry(&self) -> &SpatialGeometry;
    fn linear_source(&self, cell: usize) -> Option<&ResolvedMaterialStatePoint>;
    fn conductors(&self) -> &[Option<ResolvedConductorState>];
    fn current(&self) -> &[[f64; 3]];
    fn recover(
        &self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialInductionField, Self::Error>;
}
impl Owner for SpatialInduction<'_> {
    type Error = SpatialMagneticError;
    fn geometry(&self) -> &SpatialGeometry {
        &self.problem().geometry
    }
    fn linear_source(&self, cell: usize) -> Option<&ResolvedMaterialStatePoint> {
        match &self.problem().materials()[cell] {
            SpatialMagneticMaterial::Linear(source) => Some(source),
            SpatialMagneticMaterial::Recoil(_) => None,
        }
    }
    fn conductors(&self) -> &[Option<ResolvedConductorState>] {
        self.conductors()
    }
    fn current(&self) -> &[[f64; 3]] {
        self.impressed_current_density_a_m2()
    }
    fn recover(
        &self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialInductionField, Self::Error> {
        self.field(cancelled)
    }
}
impl Owner for SpatialNonlinearInduction<'_> {
    type Error = SpatialNonlinearInductionError;
    fn geometry(&self) -> &SpatialGeometry {
        &self.problem().geometry
    }
    fn linear_source(&self, cell: usize) -> Option<&ResolvedMaterialStatePoint> {
        match &self.problem().materials()[cell] {
            SpatialNonlinearMaterial::Linear(source) => Some(source),
            _ => None,
        }
    }
    fn conductors(&self) -> &[Option<ResolvedConductorState>] {
        self.conductors()
    }
    fn current(&self) -> &[[f64; 3]] {
        self.impressed_current_density_a_m2()
    }
    fn recover(
        &self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialInductionField, Self::Error> {
        self.field(cancelled)
    }
}

fn admit<'a>(
    owner: &'a impl Owner,
    cell: usize,
    exterior: &mut Option<&'a ResolvedMaterialStatePoint>,
) -> Result<(), SpatialMagneticError> {
    if owner.conductors()[cell].is_some() || owner.current()[cell].iter().any(|&x| x != 0.0) {
        return Err(invalid(
            "induction Maxwell collar must be nonconducting and free of impressed current",
        ));
    }
    let source = owner.linear_source(cell).ok_or_else(|| {
        invalid("induction Maxwell collar requires induced-linear exterior material")
    })?;
    if source != *exterior.get_or_insert(source) {
        return Err(invalid(
            "induction Maxwell collar requires one complete original exterior source",
        ));
    }
    Ok(())
}

fn surface<O: Owner>(
    owner: &O,
    enclosed_cells: &[usize],
    placement: &Se3,
    moment_origin_m: [f64; 3],
    limits: SpatialMaxwellForceLimits,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Wrench, O::Error> {
    poll(cancelled)?;
    let geometry = owner.geometry();
    // Reject a full-field recovery outside the user's declared mesh budget.
    if geometry.mesh.tets.len() > limits.max_cells || geometry.faces.len() > limits.max_faces {
        return Err(invalid("induction Maxwell field recovery exceeds mesh budget").into());
    }
    let field = owner.recover(cancelled)?;
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
            admit(owner, inner, &mut exterior)?;
            admit(owner, outer, &mut exterior)?;
            Ok((
                field.flux_density_t()[outer],
                field.field_strength_a_m()[outer],
            ))
        },
        cancelled,
    )
    .map_err(Into::into)
}

fn distributed<O: Owner>(
    owner: &O,
    shell: &PreparedSpatialEggshell<'_>,
    placement: &Se3,
    moment_origin_m: [f64; 3],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Wrench, O::Error> {
    let field = owner.recover(cancelled)?;
    let geometry = owner.geometry();
    let mut exterior = None;
    eggshell::apply(
        shell,
        &geometry.mesh,
        &geometry.positions_m,
        placement,
        moment_origin_m,
        |cell| {
            admit(owner, cell, &mut exterior)?;
            Ok((
                field.flux_density_t()[cell],
                field.field_strength_a_m()[cell],
            ))
        },
        cancelled,
    )
    .map_err(Into::into)
}

macro_rules! load_api {
    ($owner:ty, $error:ty) => {
        impl $owner {
            /// Instantaneous Maxwell wrench from retained transient B/H. The
            /// complete closed surface must have source-identical insulating
            /// linear material on both sides. Pose and moment origin follow
            /// `SpatialMagneticField::maxwell_wrench`; no static solve occurs.
            pub fn maxwell_wrench(
                &self,
                enclosed_cells: &[usize],
                placement: &Se3,
                moment_origin_m: [f64; 3],
                limits: SpatialMaxwellForceLimits,
                cancelled: &mut impl FnMut() -> bool,
            ) -> Result<Wrench, $error> {
                surface(
                    self,
                    enclosed_cells,
                    placement,
                    moment_origin_m,
                    limits,
                    cancelled,
                )
            }

            /// Reuse the original geometry's prepared distributed Maxwell
            /// functional on this accepted induction state. Every transition
            /// or fractional plateau must be source-identical insulating
            /// linear exterior. The weight-one target may carry eddy current.
            /// The returned canonical wrench does not advance mechanics.
            pub fn eggshell_wrench(
                &self,
                shell: &PreparedSpatialEggshell<'_>,
                placement: &Se3,
                moment_origin_m: [f64; 3],
                cancelled: &mut impl FnMut() -> bool,
            ) -> Result<Wrench, $error> {
                distributed(self, shell, placement, moment_origin_m, cancelled)
            }
        }
    };
}
load_api!(SpatialInduction<'_>, SpatialMagneticError);
load_api!(
    SpatialNonlinearInduction<'_>,
    SpatialNonlinearInductionError
);
