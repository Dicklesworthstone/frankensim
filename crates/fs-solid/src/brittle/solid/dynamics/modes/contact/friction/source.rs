//! Original material-interface sources for conditional sticking vibration.
//!
//! The existing source sampler refreshes the explicitly selected static claims
//! at each loaded pair's actual tangent speed. No slip, distance, temperature,
//! damage, physical time or body state advances. Observed source candidates are
//! returned only with a completely admitted original contact spectrum.

use super::{FrozenBrittleStickingControls, assemble_sticking_contact_modes_with};
use crate::brittle::solid::dynamics::BrittleTetMotion;
use crate::brittle::solid::dynamics::fracture::contact::friction::source::{
    BrittleFaceFrictionSample, BrittleFaceFrictionSource, BrittleFaceSourcedFrictionError,
    sample_source, validate_source_bindings,
};
use crate::brittle::solid::dynamics::fracture::contact::friction::{
    BrittleFaceFrictionFrame, EllipticFriction,
};
use crate::brittle::solid::dynamics::modes::contact::{
    FrozenBrittleContactError, FrozenBrittleContactModes, dot,
};
use crate::brittle::solid::fracture::ReleasedBrittleFace;
use crate::inertia::P1ReferenceMass;
use crate::linear3::TetAssemblyBudget;
use crate::loads3::TetLoadGeometry3;
use fs_contact::interface_binding::thermal_friction::elliptic::impulse::SourceImpulseCoefficientRole;
use fs_exec::Cx;
use fs_matdb::MaterialStateId;
use fs_material::brittle::{BrittleCrackBand, BrittleState};

/// One complete conditional spectrum and exactly the loaded source snapshots.
#[derive(Debug)]
pub struct FrozenBrittleSourcedStickingModes {
    pub contact: FrozenBrittleContactModes,
    /// Original face/vertex identities survive constraint alias reduction.
    /// These are coefficient observations, not a physical history advance.
    pub samples: Vec<BrittleFaceFrictionSample>,
}

#[derive(Debug)]
pub enum FrozenBrittleSourcedStickingError {
    Contact(FrozenBrittleContactError),
    Source(BrittleFaceSourcedFrictionError),
}

impl From<FrozenBrittleContactError> for FrozenBrittleSourcedStickingError {
    fn from(error: FrozenBrittleContactError) -> Self {
        Self::Contact(error)
    }
}
impl From<BrittleFaceSourcedFrictionError> for FrozenBrittleSourcedStickingError {
    fn from(error: BrittleFaceSourcedFrictionError) -> Self {
        Self::Source(error)
    }
}
impl core::fmt::Display for FrozenBrittleSourcedStickingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source-bound brittle sticking modes refused: {self:?}")
    }
}
impl std::error::Error for FrozenBrittleSourcedStickingError {}

impl BrittleTetMotion<'_> {
    /// Observe a held sticking spectrum using the original interface claims.
    ///
    /// Sources follow retained face and vertex order, name the actual ordered
    /// material states, and explicitly select the Static coefficient role.
    /// Only loaded pairs sample their source, at their own measured tangent
    /// speed in the declared frame. Accepted path, slip, texture and finite-step
    /// temperature histories remain unchanged; external source conditions are
    /// held. No normal pressure is inferred from the supplied nodal preload.
    ///
    /// Native static-ellipse admission, original mass and fixed-damage stiffness
    /// determine the same exact constrained modes as the explicit-law route.
    /// Neither the preload nor its source establishes mechanical equilibrium
    /// or guarantees finite-amplitude sticking. A late source, modal, budget or
    /// cancellation refusal returns no spectrum or source snapshots.
    pub fn frozen_sourced_sticking_contact_modes(
        &self,
        sources: &[BrittleFaceFrictionSource<'_>],
        preload_n: &[[[f64; 3]; 3]],
        fixed_dofs: &[usize],
        controls: FrozenBrittleStickingControls,
        cx: &Cx<'_>,
    ) -> Result<FrozenBrittleSourcedStickingModes, FrozenBrittleSourcedStickingError> {
        assemble_sourced_sticking_contact_modes(
            self.solid.geometry,
            self.solid.budget,
            self.solid.materials.iter().map(|source| source.law()),
            |element| {
                self.solid
                    .materials
                    .get(element)
                    .map(|source| source.material().material())
            },
            &self.solid.solution.element_states,
            &self.solid.solution.displacement_m,
            &self.velocity_mps,
            &self.solid.released_faces,
            &self.mass,
            sources,
            preload_n,
            fixed_dofs,
            controls,
            cx,
        )
    }
}

/// Mechanical and thermal motions share binding, sampling and modal admission.
#[allow(clippy::too_many_arguments)]
pub(crate) fn assemble_sourced_sticking_contact_modes<'m>(
    geometry: &TetLoadGeometry3<'_>,
    budget: TetAssemblyBudget,
    materials: impl ExactSizeIterator<Item = &'m BrittleCrackBand>,
    material: impl Fn(usize) -> Option<&'m MaterialStateId>,
    states: &[BrittleState<'_>],
    displacement_m: &[[f64; 3]],
    velocity_mps: &[[f64; 3]],
    released_faces: &[ReleasedBrittleFace],
    mass: &P1ReferenceMass,
    sources: &[BrittleFaceFrictionSource<'_>],
    preload_n: &[[[f64; 3]; 3]],
    fixed_dofs: &[usize],
    controls: FrozenBrittleStickingControls,
    cx: &Cx<'_>,
) -> Result<FrozenBrittleSourcedStickingModes, FrozenBrittleSourcedStickingError> {
    validate_source_bindings(released_faces, material, sources, cx)?;
    for source in sources {
        cx.checkpoint()
            .map_err(|_| BrittleFaceSourcedFrictionError::Cancelled)?;
        if source
            .vertices
            .iter()
            .any(|vertex| vertex.coefficient_role != SourceImpulseCoefficientRole::Static)
        {
            return Err(BrittleFaceSourcedFrictionError::Invalid(
                "sticking modes require the explicit original Static coefficient role",
            )
            .into());
        }
    }
    let mut samples = Vec::new();
    let mut frame = None;
    let contact = assemble_sticking_contact_modes_with(
        geometry,
        budget,
        materials,
        states,
        displacement_m,
        velocity_mps,
        released_faces,
        mass,
        sources.len(),
        preload_n,
        fixed_dofs,
        controls,
        cx,
        |face_index, normal, relative_velocity| {
            let source = sources[face_index];
            let local = match frame {
                Some((index, value)) if index == face_index => value,
                _ => {
                    let value = BrittleFaceFrictionFrame::from_unit_normal(
                        normal.direction,
                        source.tangent_direction,
                    )
                    .map_err(BrittleFaceSourcedFrictionError::from)?;
                    frame = Some((face_index, value));
                    value
                }
            };
            let velocity = [
                dot(local.normal, relative_velocity),
                dot(local.tangent1, relative_velocity),
                dot(local.tangent2, relative_velocity),
            ];
            let sample =
                sample_source(source.vertices[normal.row.vertex], normal.row, velocity, cx)?;
            let friction = EllipticFriction {
                coefficients_mu: sample.sample.coefficients_mu,
                principal_axis: sample.sample.state.constitutive_state().principal_axis(),
            };
            samples.push(sample);
            Ok::<_, FrozenBrittleSourcedStickingError>((friction, source.tangent_direction))
        },
    )?;
    Ok(FrozenBrittleSourcedStickingModes { contact, samples })
}
