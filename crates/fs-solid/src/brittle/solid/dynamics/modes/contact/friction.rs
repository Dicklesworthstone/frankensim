//! Conditional sticking modes inside the native static friction ellipse.
//!
//! Strictly interior physical preload and negligible tangential relative speed
//! select two exact tangent constraints. The original normal constraints, bulk
//! frozen stiffness and reference inertia remain unchanged. Friction does not
//! become a tangent spring or modal damping. A zero/zero ellipse leaves tangent
//! motion free; a single zero coefficient retains the native sustained-law
//! refusal instead of acquiring fictitious traction through regularization.
//!
//! The supplied preload is a held branch assumption, not an equilibrium proof.
//! Shapes and reaction increments use the same physical coordinates as the
//! original contact owner. No finite oscillation amplitude, future stick/slip
//! transition, heat, finite sliding or material-source authority is inferred.

/// Original material-interface sources on the held static-friction branch.
pub mod source;

use super::{
    BrittleModeConstraint, FrozenBrittleContactControls, FrozenBrittleContactError,
    FrozenBrittleContactModes, admit_counts, assemble_frozen_pencil, dot, normal_constraints,
    project_contact_pencil,
};
use crate::brittle::solid::checkpoint;
use crate::brittle::solid::dynamics::BrittleTetMotion;
use crate::brittle::solid::dynamics::fracture::contact::friction::{
    BrittleFaceFriction, BrittleFaceFrictionFrame, EllipticFriction,
};
use crate::brittle::solid::fracture::ReleasedBrittleFace;
use crate::inertia::P1ReferenceMass;
use crate::linear3::TetAssemblyBudget;
use crate::loads3::TetLoadGeometry3;
use fs_contact::sustained::elliptic::EllipticForceCoordinates;
use fs_exec::Cx;
use fs_material::brittle::{BrittleCrackBand, BrittleState};

/// Native normal admission plus an explicit fractional distance from sliding.
#[derive(Clone, Copy, Debug)]
pub struct FrozenBrittleStickingControls {
    /// The row cap counts ALL scalar constraints: normal and both tangents.
    pub contact: FrozenBrittleContactControls,
    /// In [0,1). Transformed tangent load / native disk capacity must be
    /// STRICTLY less than 1-margin; zero still excludes the sliding boundary.
    pub minimum_friction_margin: f64,
}

impl Default for FrozenBrittleStickingControls {
    fn default() -> Self {
        Self {
            contact: FrozenBrittleContactControls::default(),
            minimum_friction_margin: 1e-6,
        }
    }
}

impl BrittleTetMotion<'_> {
    /// Observe loaded sticking contact through exact native trace constraints.
    ///
    /// `friction` and `preload_n` follow retained face order. Each preload is
    /// [vertex][normal,tangent1,tangent2], in the native frame formed by that
    /// face's tangent direction. Open and tied pairs require all three forces
    /// to be zero and do not evaluate inactive friction inputs. Loaded pairs
    /// with two positive coefficients must rest tangentially and lie strictly
    /// inside the original ellipse; its boundary and sliding states refuse.
    /// Two zero coefficients require zero tangent preload and retain free slip.
    ///
    /// The caller holds the contact state and separately establishes equilibrium
    /// when equilibrium frequencies are required. History and time never advance.
    pub fn frozen_sticking_contact_modes(
        &self,
        friction: &[BrittleFaceFriction],
        preload_n: &[[[f64; 3]; 3]],
        fixed_dofs: &[usize],
        controls: FrozenBrittleStickingControls,
        cx: &Cx<'_>,
    ) -> Result<FrozenBrittleContactModes, FrozenBrittleContactError> {
        assemble_sticking_contact_modes(
            self.solid.geometry,
            self.solid.budget,
            self.solid.materials.iter().map(|source| source.law()),
            &self.solid.solution.element_states,
            &self.solid.solution.displacement_m,
            &self.velocity_mps,
            &self.solid.released_faces,
            &self.mass,
            friction,
            preload_n,
            fixed_dofs,
            controls,
            cx,
        )
    }
}

/// Thermal and mechanical histories use identical physical contact reduction.
#[allow(clippy::too_many_arguments)]
pub(crate) fn assemble_sticking_contact_modes<'m>(
    geometry: &TetLoadGeometry3<'_>,
    budget: TetAssemblyBudget,
    materials: impl ExactSizeIterator<Item = &'m BrittleCrackBand>,
    states: &[BrittleState<'_>],
    displacement_m: &[[f64; 3]],
    velocity_mps: &[[f64; 3]],
    released_faces: &[ReleasedBrittleFace],
    mass: &P1ReferenceMass,
    friction: &[BrittleFaceFriction],
    preload_n: &[[[f64; 3]; 3]],
    fixed_dofs: &[usize],
    controls: FrozenBrittleStickingControls,
    cx: &Cx<'_>,
) -> Result<FrozenBrittleContactModes, FrozenBrittleContactError> {
    assemble_sticking_contact_modes_with(
        geometry,
        budget,
        materials,
        states,
        displacement_m,
        velocity_mps,
        released_faces,
        mass,
        friction.len(),
        preload_n,
        fixed_dofs,
        controls,
        cx,
        |face, normal, _velocity| {
            Ok((
                friction[face].vertices[normal.row.vertex],
                friction[face].tangent_direction,
            ))
        },
    )
}

/// Only a loaded original vertex requests coefficients. Source callers keep
/// their candidates private until the entire pencil and spectrum are accepted.
#[allow(clippy::too_many_arguments)]
pub(crate) fn assemble_sticking_contact_modes_with<'m, E>(
    geometry: &TetLoadGeometry3<'_>,
    budget: TetAssemblyBudget,
    materials: impl ExactSizeIterator<Item = &'m BrittleCrackBand>,
    states: &[BrittleState<'_>],
    displacement_m: &[[f64; 3]],
    velocity_mps: &[[f64; 3]],
    released_faces: &[ReleasedBrittleFace],
    mass: &P1ReferenceMass,
    friction_faces: usize,
    preload_n: &[[[f64; 3]; 3]],
    fixed_dofs: &[usize],
    controls: FrozenBrittleStickingControls,
    cx: &Cx<'_>,
    mut friction_at: impl FnMut(
        usize,
        BrittleModeConstraint,
        [f64; 3],
    ) -> Result<(EllipticFriction, [f64; 3]), E>,
) -> Result<FrozenBrittleContactModes, E>
where
    E: From<FrozenBrittleContactError>,
{
    use FrozenBrittleContactError as Error;
    checkpoint(cx).map_err(Error::from)?;
    let full = geometry
        .nodes_m()
        .len()
        .checked_mul(3)
        .ok_or(Error::Invalid("sticking modal DOF overflow"))?;
    let free = full
        .checked_sub(fixed_dofs.len())
        .ok_or(Error::Invalid("too many fixed sticking modal DOFs"))?;
    admit_counts(full, free, 0, controls.contact)?;
    if !controls.minimum_friction_margin.is_finite()
        || !(0.0..1.0).contains(&controls.minimum_friction_margin)
        || displacement_m.len() != geometry.nodes_m().len()
        || velocity_mps.len() != geometry.nodes_m().len()
        || preload_n.len() != released_faces.len()
        || friction_faces != released_faces.len()
    {
        return Err(Error::Invalid("sticking modal source dimensions or margin").into());
    }
    let mut constraints = Vec::new();
    for (face_index, face) in released_faces.iter().enumerate() {
        checkpoint(cx).map_err(Error::from)?;
        let preload = preload_n[face_index];
        // A stack-sized native normal admission avoids copying all inactive
        // faces into a new preload buffer. The cumulative scalar cap below
        // still includes every admitted normal and tangent row.
        let normals = normal_constraints(
            displacement_m,
            velocity_mps,
            core::slice::from_ref(face),
            &[preload.map(|force| force[0])],
            controls.contact,
            cx,
        )?;
        for (vertex, force) in preload.iter().enumerate() {
            if !normals.iter().any(|row| row.row.vertex == vertex)
                && force.iter().any(|value| *value != 0.0)
            {
                return Err(Error::Contact {
                    face: face.id,
                    vertex,
                    reason: "open or tied contact cannot carry any preload",
                }
                .into());
            }
        }
        let mut frame = None;
        for normal in normals {
            checkpoint(cx).map_err(Error::from)?;
            let vertex = normal.row.vertex;
            let fail = |reason| Error::Contact {
                face: face.id,
                vertex,
                reason,
            };
            let force = preload[vertex];
            let [a, b] = normal.row.nodes;
            let relative =
                core::array::from_fn(|axis| velocity_mps[b][axis] - velocity_mps[a][axis]);
            // Even a frictionless active point needs one scalar row. Check
            // that minimum before a source callback allocates its snapshot.
            admit_counts(full, free, constraints.len() + 1, controls.contact)?;
            let (friction, tangent_direction) = friction_at(face_index, normal, relative)?;
            let coordinates = EllipticForceCoordinates::prepare(friction)
                .map_err(|_| fail("native sustained friction ellipse refused"))?;
            let q = coordinates
                .force_to_solver(force)
                .map_err(|_| fail("native preload force coordinates refused"))?;
            let sticking = coordinates.coefficient() > 0.0;
            let additional = if sticking { 3 } else { 1 };
            let total = constraints
                .len()
                .checked_add(additional)
                .ok_or(Error::Invalid("sticking constraint count overflow"))?;
            admit_counts(full, free, total, controls.contact)?;
            if !sticking {
                if force[1] != 0.0 || force[2] != 0.0 {
                    return Err(fail("zero friction cannot carry tangent preload").into());
                }
                constraints.push(normal);
                continue;
            }
            let capacity = coordinates.coefficient() * force[0];
            let magnitude = q[1].hypot(q[2]);
            let utilization = magnitude / capacity;
            if !capacity.is_finite()
                || capacity <= 0.0
                || !magnitude.is_finite()
                || !utilization.is_finite()
                || utilization >= 1.0 - controls.minimum_friction_margin
            {
                return Err(fail(
                    "preload is not representably inside the static friction ellipse",
                )
                .into());
            }
            let local_frame = match frame {
                Some(frame) => frame,
                None => {
                    let built =
                        BrittleFaceFrictionFrame::from_unit_normal(face.normal, tangent_direction)
                            .map_err(|_| fail("native friction frame refused"))?;
                    frame = Some(built);
                    built
                }
            };
            let speed1 = dot(local_frame.tangent1, relative);
            let speed2 = dot(local_frame.tangent2, relative);
            if !speed1.is_finite()
                || !speed2.is_finite()
                || speed1.hypot(speed2) > controls.contact.velocity_tolerance_mps
            {
                return Err(fail("sticking contact must be tangentially resting").into());
            }
            constraints.push(normal);
            constraints.push(BrittleModeConstraint {
                row: normal.row,
                direction: local_frame.tangent1,
            });
            constraints.push(BrittleModeConstraint {
                row: normal.row,
                direction: local_frame.tangent2,
            });
        }
    }
    let pencil = assemble_frozen_pencil(
        geometry,
        budget,
        materials,
        states,
        displacement_m,
        mass,
        fixed_dofs,
        controls.contact.assembly,
        cx,
    )
    .map_err(Error::from)?;
    project_contact_pencil(pencil, constraints, controls.contact, cx).map_err(E::from)
}
