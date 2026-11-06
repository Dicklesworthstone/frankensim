//! Explicit material-interface friction on native released-face impacts.
//!
//! Each active paired vertex samples its original source at actual incoming
//! slip in the declared face frame. Source candidates and the mechanical
//! successor are returned together only after the complete event is accepted.

use super::*;
use fs_contact::interface_binding::thermal_friction::elliptic::impulse::{
    SourceFrictionConditions, SourceImpulseCoefficientRole, SourceImpulseError, SourceImpulseSample,
};
use fs_contact::interface_binding::thermal_friction::elliptic::{
    SourceThermalEllipticCoulomb, SourceThermalEllipticState,
};
use fs_matdb::MaterialStateId;

/// An explicit source and accepted history for one original paired vertex.
#[derive(Clone, Copy, Debug)]
pub struct BrittleVertexFrictionSource<'s> {
    pub law: &'s SourceThermalEllipticCoulomb,
    pub state: &'s SourceThermalEllipticState,
    /// Complete externally supplied conditions. Declared speed coordinates use
    /// actual incoming slip; accumulated distance stays at its accepted value.
    pub conditions: &'s SourceFrictionConditions,
    /// One fixed static or kinetic ellipse for this instantaneous event.
    pub coefficient_role: SourceImpulseCoefficientRole,
}

/// Source ownership follows retained face order and original vertex order.
#[derive(Clone, Copy, Debug)]
pub struct BrittleFaceFrictionSource<'s> {
    pub face_id: u64,
    /// The source's retained texture coordinates must use the frame constructed
    /// from this direction and the retained face normal. No texture is inferred.
    pub tangent_direction: [f64; 3],
    pub vertices: [BrittleVertexFrictionSource<'s>; 3],
}

/// Accepted coefficient/source candidate, identified independently of row order.
#[derive(Clone, Debug)]
pub struct BrittleFaceFrictionSample {
    pub row: BrittleFaceImpactRow,
    pub sample: SourceImpulseSample,
}

/// One accepted body event and exactly its active original source candidates.
#[derive(Debug)]
pub struct BrittleFaceSourcedFrictionStep<'a> {
    pub impact: BrittleFaceFrictionStep<'a>,
    /// Inactive and tied vertices are not sampled and retain their old history.
    pub samples: Vec<BrittleFaceFrictionSample>,
}

#[derive(Debug)]
pub enum BrittleFaceSourcedFrictionError {
    Impact(BrittleFaceFrictionError),
    Invalid(&'static str),
    Binding {
        face: u64,
        reason: &'static str,
    },
    Source {
        row: BrittleFaceImpactRow,
        source: SourceImpulseError,
    },
    Cancelled,
}
impl From<BrittleFaceFrictionError> for BrittleFaceSourcedFrictionError {
    fn from(error: BrittleFaceFrictionError) -> Self {
        Self::Impact(error)
    }
}
impl core::fmt::Display for BrittleFaceSourcedFrictionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "source-bound released brittle face impact refused: {self:?}"
        )
    }
}
impl std::error::Error for BrittleFaceSourcedFrictionError {}

impl<'a> BrittleTetMotion<'a> {
    /// Resolve native directional friction from explicitly selected interfaces.
    ///
    /// Each source must name the actual ordered material states on its face.
    /// Native source ownership retains selected claims, validity, authority,
    /// accepted slip, path, transported texture and finite-step thermal history.
    /// A fixed coefficient role is explicit; this event does not invent a
    /// static-to-kinetic switch, pressure from impulse, finite travel or heat.
    ///
    /// Every active row uses its own actual preimpact relative speed. Inactive
    /// and tied rows never evaluate a source. No input history is mutated: a
    /// source, geometry, cancellation, mobility or late work refusal returns no
    /// body or source successor. The caller retains all accepted inputs.
    pub fn impact_released_faces_sourced_frictional(
        &self,
        restitution: &[f64],
        sources: &[BrittleFaceFrictionSource<'_>],
        held_velocity_dofs: &[usize],
        origin_m: [f64; 3],
        controls: BrittleFaceFrictionControls,
        cx: &Cx<'_>,
    ) -> Result<BrittleFaceSourcedFrictionStep<'a>, BrittleFaceSourcedFrictionError> {
        let (response, samples) = released_face_sourced_frictional_impact(
            ReleasedFaceImpactBody {
                reference_nodes_m: self.solid.geometry.nodes_m(),
                displacement_m: &self.solid.solution.displacement_m,
                released_faces: &self.solid.released_faces,
                mass: &self.mass,
                velocity_mps: &self.velocity_mps,
                kinetic_energy_j: self.kinetic_energy_j,
            },
            |element| {
                self.solid
                    .materials
                    .get(element)
                    .map(|source| source.material().material())
            },
            restitution,
            sources,
            held_velocity_dofs,
            origin_m,
            controls,
            cx,
        )?;
        let impact = BrittleFaceFrictionStep::from_response(self, response, cx)?;
        Ok(BrittleFaceSourcedFrictionStep { impact, samples })
    }
}

/// Both thermal and isothermal owners retain their actual material identities.
pub(crate) fn released_face_sourced_frictional_impact<'m>(
    body: ReleasedFaceImpactBody<'_>,
    material: impl Fn(usize) -> Option<&'m MaterialStateId>,
    restitution: &[f64],
    sources: &[BrittleFaceFrictionSource<'_>],
    held_velocity_dofs: &[usize],
    origin_m: [f64; 3],
    controls: BrittleFaceFrictionControls,
    cx: &Cx<'_>,
) -> Result<
    (ReleasedFaceFrictionResponse, Vec<BrittleFaceFrictionSample>),
    BrittleFaceSourcedFrictionError,
> {
    let poll = || {
        cx.checkpoint()
            .map_err(|_| BrittleFaceSourcedFrictionError::Cancelled)
    };
    poll()?;
    if sources.len() != body.released_faces.len() {
        return Err(BrittleFaceSourcedFrictionError::Invalid(
            "one explicit source binding per retained face required",
        ));
    }
    for (face, binding) in body.released_faces.iter().zip(sources) {
        poll()?;
        if face.id != binding.face_id {
            return Err(BrittleFaceSourcedFrictionError::Binding {
                face: face.id,
                reason: "source face IDs must match retained input order",
            });
        }
        for vertex in &binding.vertices {
            for principal in vertex.law.principal_sources() {
                let actual = principal.interface().surface_materials();
                if (0..2).any(|side| material(face.elements[side]) != Some(&actual[side])) {
                    return Err(BrittleFaceSourcedFrictionError::Binding {
                        face: face.id,
                        reason: "ordered interface material states must match incident elements",
                    });
                }
            }
        }
    }
    let directions: Vec<_> = sources
        .iter()
        .map(|source| source.tangent_direction)
        .collect();
    let mut samples = Vec::new();
    let response = released_face_frictional_impact_with(
        body,
        restitution,
        &directions,
        held_velocity_dofs,
        origin_m,
        controls,
        cx,
        |face_index, row, _frame, velocity| {
            let binding = sources[face_index].vertices[row.vertex];
            let sample = binding
                .law
                .sample_preimpact_coefficients(
                    binding.state,
                    [velocity[1], velocity[2]],
                    binding.conditions,
                    binding.coefficient_role,
                    &mut || cx.checkpoint().is_err(),
                )
                .map_err(|source| BrittleFaceSourcedFrictionError::Source { row: *row, source })?;
            let friction = EllipticFriction {
                coefficients_mu: sample.coefficients_mu,
                principal_axis: sample.state.constitutive_state().principal_axis(),
            };
            samples.push(BrittleFaceFrictionSample { row: *row, sample });
            Ok::<_, BrittleFaceSourcedFrictionError>(friction)
        },
    )?;
    Ok((response, samples))
}
