//! Native primal-ray traversal through a continuous P1 thermal mesh.
//!
//! Faces join only through the original vertex IDs. One reference glass and
//! thermo-optic law apply throughout; shared nodal temperatures therefore give
//! continuous n, although grad(n) can jump. Such a face changes curvature but
//! introduces no Fresnel interface or instantaneous direction jump. Native
//! adaptive integration evaluates the original element on each side of the
//! face, including the numerical arrival/departure collar. No arc length,
//! phase/group integral, or absorption is omitted at a transition.
//!
//! Preparation stores adjacency once. Propagation keeps only two prepared
//! elements, allocates nothing, and reuses the native boundary and ODE owners.
//! The two-cell collar excludes every other face of both elements. Elsewhere
//! every stage stays in the current element's native geometric free ball.
//!
//! This does not locate an arbitrary launch, certify global mesh embedding,
//! cross duplicated contact nodes, or supply other scene-object geometry.
//! Native local ODE error estimates remain local numerical estimates. Bundle
//! traversal is deliberately absent: a gradient jump also needs the correct
//! crossing-time variation (saltation), not merely piecewise interior Hessians.

use super::{
    ConductionElementBoundary, P1ConductionElement, P1ConductionTemperatureSnapshot,
    TemperatureGradientField, TemperatureGradientSample, TemperatureJet,
    TransientTemperatureFieldError,
};
use crate::charts::Ray;
use crate::charts::proximity::{ProximitySide, triangle_proximity};
use crate::gradient_index::event::{
    BoundaryTermination, GradientDeparture, GradientEventOptions,
    propagate_to_boundary_with_control,
};
use crate::gradient_index::thermal::TemperatureField;
use crate::gradient_index::{GradientIndexError, GradientIndexPropagation};
use crate::material_state::thermo_optic::{MaterialThermoOpticField, ThermoOpticFieldError};
use core::cell::Cell;
use fs_exec::Cx;
use fs_geom::{Point3, Vec3};
use fs_math::dd::Dd;

mod transition;
use transition::{MeshBoundary, Transition};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct ElementFace {
    element: usize,
    local_face: usize,
}

/// One accepted original temperature snapshot and its face-ID adjacency.
/// Preparation is O(tetrahedra * log(faces)) and allocates O(tetrahedra + faces).
/// The original mesh must represent the intended conforming body; this owner
/// checks shared-face incidence and orientation, not global self-intersection.
pub struct P1ConductionMesh<'a> {
    snapshot: P1ConductionTemperatureSnapshot<'a>,
    neighbours: Vec<[Option<ElementFace>; 4]>,
    face_ids: Vec<[usize; 4]>,
}

impl<'a> P1ConductionMesh<'a> {
    pub fn try_new(
        cx: &Cx<'_>,
        snapshot: P1ConductionTemperatureSnapshot<'a>,
    ) -> Result<Self, P1MeshError> {
        let complex = snapshot.mesh().complex();
        let mut owners: Vec<Option<ElementFace>> = vec![None; complex.faces.len()];
        let mut neighbours = vec![[None; 4]; complex.tets.len()];
        let mut face_ids = vec![[0; 4]; complex.tets.len()];
        for (element, tet) in complex.tets.iter().enumerate() {
            if element % 128 == 0 {
                cx.checkpoint().map_err(|_| GradientIndexError::Cancelled)?;
            }
            for local_face in 0..4 {
                let mut vertices = [0; 3];
                let mut i = 0;
                for (j, &vertex) in tet.iter().enumerate() {
                    if j != local_face {
                        vertices[i] = vertex;
                        i += 1;
                    }
                }
                vertices.sort_unstable();
                let face = complex.faces.binary_search(&vertices).map_err(|_| {
                    P1MeshError::InvalidTopology {
                        element,
                        local_face,
                    }
                })?;
                face_ids[element][local_face] = face;
                let current = ElementFace {
                    element,
                    local_face,
                };
                if let Some(other) = owners[face] {
                    if neighbours[other.element][other.local_face].is_some() {
                        return Err(P1MeshError::InvalidTopology {
                            element,
                            local_face,
                        });
                    }
                    let previous = snapshot.element(other.element)?;
                    let next = snapshot.element(element)?;
                    let side = triangle_proximity(
                        next.positions[local_face],
                        previous.faces_m()[other.local_face],
                        other.local_face as u32,
                    )
                    .map_err(|_| P1MeshError::InvalidTopology {
                        element,
                        local_face,
                    })?
                    .1
                    .side;
                    // A topological match is not permission to enter an
                    // overlapping tetrahedron on the same geometric side.
                    if side != ProximitySide::Positive {
                        return Err(P1MeshError::InvalidTopology {
                            element,
                            local_face,
                        });
                    }
                    neighbours[other.element][other.local_face] = Some(current);
                    neighbours[element][local_face] = Some(other);
                } else {
                    owners[face] = Some(current);
                }
            }
        }
        Ok(Self {
            snapshot,
            neighbours,
            face_ids,
        })
    }

    pub const fn snapshot(&self) -> P1ConductionTemperatureSnapshot<'a> {
        self.snapshot
    }

    /// Propagate the original thermo-optic binding from its selected element.
    /// The binding must borrow this exact mesh and nodal slice. Its original
    /// material laws, spectral domains and snapshot identity remain in force.
    ///
    /// An exterior or duplicated contact-node face returns an ordinary native
    /// `Surface` termination. Its feature ID is the original mesh's canonical
    /// face index; `options.boundary_id` identifies the enclosing mesh. A
    /// distance-limited endpoint can lie inside a transition collar and reports
    /// the original element actually sampled there. No endpoint is moved onto
    /// a face and no part of the requested arc length is skipped.
    #[allow(clippy::too_many_arguments)]
    pub fn propagate(
        &self,
        cx: &Cx<'_>,
        binding: &MaterialThermoOpticField<'_, P1ConductionElement<'a>>,
        ray: Ray,
        transverse_axis: Vec3,
        wavelength_nm: f64,
        max_distance_m: f64,
        options: P1MeshOptions,
    ) -> Result<P1MeshPropagation, P1MeshError> {
        options.event.validate()?;
        if options.max_internal_faces > MAX_INTERNAL_FACES {
            return Err(GradientIndexError::InvalidOptions.into());
        }
        let start = *binding.temperature_source();
        let source = start.snapshot();
        if !core::ptr::eq(source.mesh(), self.snapshot.mesh())
            || !core::ptr::eq(
                source.nodal_temperature_k(),
                self.snapshot.nodal_temperature_k(),
            )
        {
            return Err(P1MeshError::ForeignSnapshot);
        }
        let cursor = MeshTemperature {
            current: Cell::new(start),
            previous: Cell::new(None),
        };
        let original = MaterialThermoOpticField::try_new(
            binding.reference_glass(),
            binding.law(),
            &cursor,
            binding.snapshot_identity(),
        )?;
        let prepared = original.prepare(wavelength_nm)?;
        let mut boundary = MeshBoundary::new(
            start.boundary_oracle(options.boundary_id, options.event.surface_tolerance_m)?,
        );
        let mut departure: Option<GradientDeparture> = None;
        let mut current_ray = ray;
        let mut current_axis = transverse_axis;
        let mut total: Option<GradientIndexPropagation> = None;
        // Accumulate only between element segments. The original hot ODE
        // arithmetic is unchanged; small collars survive long paths.
        let mut sums = [Dd::ZERO; 4]; // arc, phase, group, absorption
        let mut boundary_evaluations = 0;
        let mut internal_faces = 0;
        let mut contact_resolution_m = 0.0_f64;
        loop {
            cx.checkpoint().map_err(|_| GradientIndexError::Cancelled)?;
            let mut event_options = options.event;
            if let Some(previous) = total {
                event_options.integration.max_steps = event_options
                    .integration
                    .max_steps
                    .checked_sub(previous.accepted_steps + previous.rejected_steps)
                    .filter(|&n| n > 0)
                    .ok_or(GradientIndexError::StepBudgetExceeded)?;
                event_options.max_boundary_queries = event_options
                    .max_boundary_queries
                    .checked_sub(boundary_evaluations)
                    .filter(|&n| n > 0)
                    .ok_or(GradientIndexError::BoundaryQueryBudgetExceeded)?;
            }
            let remaining = (Dd::from_f64(max_distance_m) - sums[0]).to_f64();
            let segment = propagate_to_boundary_with_control(
                &prepared,
                &mut boundary,
                current_ray,
                current_axis,
                wavelength_nm,
                remaining,
                departure,
                event_options,
                || cx.checkpoint().is_ok(),
            )?;
            let p = segment.propagation;
            let increments = [
                p.distance_m,
                p.optical_path_length_m,
                p.group_optical_path_length_m.unwrap_or(0.0),
                p.optical_depth,
            ];
            for (sum, value) in sums.iter_mut().zip(increments) {
                *sum = *sum + Dd::from_f64(value);
                if !sum.hi.is_finite() || !sum.lo.is_finite() {
                    return Err(GradientIndexError::NumericalRange.into());
                }
            }
            boundary_evaluations += segment.boundary_evaluations;
            contact_resolution_m = contact_resolution_m.max(segment.contact_resolution_m);
            let mut propagation = p;
            if let Some(previous) = total {
                propagation.initial_sample = previous.initial_sample;
                propagation.accepted_steps += previous.accepted_steps;
                propagation.rejected_steps += previous.rejected_steps;
                propagation.field_evaluations += previous.field_evaluations;
                propagation.max_normalized_error_estimate = propagation
                    .max_normalized_error_estimate
                    .max(previous.max_normalized_error_estimate);
                if previous.group_optical_path_length_m.is_some()
                    != p.group_optical_path_length_m.is_some()
                {
                    return Err(GradientIndexError::GroupIndexUnavailable.into());
                }
            }
            propagation.distance_m = sums[0].to_f64();
            propagation.optical_path_length_m = sums[1].to_f64();
            propagation.group_optical_path_length_m =
                p.group_optical_path_length_m.map(|_| sums[2].to_f64());
            propagation.optical_depth = sums[3].to_f64();
            let reached_distance = segment.termination == BoundaryTermination::DistanceLimit
                || propagation.distance_m >= max_distance_m;
            let termination = if reached_distance {
                // The final native segment integrated precisely its remaining
                // requested arc. Preserve the public request's exact f64 endpoint.
                propagation.distance_m = max_distance_m;
                Some(BoundaryTermination::DistanceLimit)
            } else if let BoundaryTermination::Surface(mut witness) = segment.termination {
                let active = cursor.current.get();
                let local_face = witness.feature_id;
                let Some(next) = self.neighbours[active.element][local_face] else {
                    witness.feature_id = self.face_ids[active.element][local_face];
                    return Ok(P1MeshPropagation {
                        propagation,
                        termination: BoundaryTermination::Surface(witness),
                        element_index: active.element,
                        internal_faces,
                        boundary_evaluations,
                        contact_resolution_m,
                    });
                };
                if internal_faces >= options.max_internal_faces {
                    return Err(P1MeshError::FaceBudgetExceeded);
                }
                let neighbour = self.snapshot.element(next.element)?;
                let transition =
                    Transition::new(active, neighbour, next.local_face, p.ray.origin, witness);
                departure = Some(transition.departure());
                boundary = MeshBoundary::with_transition(
                    neighbour
                        .boundary_oracle(options.boundary_id, options.event.surface_tolerance_m)?,
                    transition,
                );
                cursor.previous.set(Some(active));
                cursor.current.set(neighbour);
                internal_faces += 1;
                None
            } else {
                unreachable!("all native boundary terminations were handled")
            };
            if let Some(termination) = termination {
                return Ok(P1MeshPropagation {
                    element_index: cursor.sample(p.ray.origin)?.0,
                    propagation,
                    termination,
                    internal_faces,
                    boundary_evaluations,
                    contact_resolution_m,
                });
            }
            total = Some(propagation);
            current_ray = p.ray;
            current_axis = p.transverse_axis;
        }
    }
}

/// Hard work ceiling, independent of the native attempted-step/query ceilings.
pub const MAX_INTERNAL_FACES: u32 = 1_000_000;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct P1MeshOptions {
    /// Native tolerances and WHOLE-path step/query budgets, including collars.
    pub event: GradientEventOptions,
    /// Maximum resolved internal-face continuations. A distance limit inside
    /// the collar can stop before mathematically crossing the last such face.
    /// Zero permits propagation only inside the selected original tetrahedron.
    pub max_internal_faces: u32,
    /// Identity returned for a boundary of the enclosing original mesh.
    pub boundary_id: usize,
}

impl Default for P1MeshOptions {
    fn default() -> Self {
        Self {
            event: GradientEventOptions::default(),
            max_internal_faces: 4096,
            boundary_id: 0,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct P1MeshPropagation {
    pub propagation: GradientIndexPropagation,
    pub termination: BoundaryTermination,
    /// Original element containing the returned integrated endpoint.
    pub element_index: usize,
    /// Resolved internal-face continuations, including a final partial collar.
    pub internal_faces: u32,
    pub boundary_evaluations: u32,
    /// Maximum native numerical contact-resolution scale; not a global bound.
    pub contact_resolution_m: f64,
}

// Routing is private per propagation. Its snapshot and pointwise temperature
// never change. Queries do not change ownership, and the native event guard
// confines the previous-cell fallback to the admitted two-cell departure.
struct MeshTemperature<'a> {
    current: Cell<P1ConductionElement<'a>>,
    previous: Cell<Option<P1ConductionElement<'a>>>,
}

impl MeshTemperature<'_> {
    fn sample(&self, p: Point3) -> Result<(usize, TemperatureJet), GradientIndexError> {
        let current = self.current.get();
        match current.sample_temperature_jet(p) {
            Ok(jet) => Ok((current.element, jet)),
            Err(TransientTemperatureFieldError::OutsideElement) => {
                let previous = self
                    .previous
                    .get()
                    .ok_or(GradientIndexError::OutsideFieldDomain)?;
                previous
                    .sample_temperature_jet(p)
                    .map(|jet| (previous.element, jet))
                    .map_err(TransientTemperatureFieldError::into_gradient_error)
            }
            Err(error) => Err(error.into_gradient_error()),
        }
    }
}
impl TemperatureField for MeshTemperature<'_> {
    fn temperature_k(&self, p: Point3) -> Result<f64, GradientIndexError> {
        Ok(self.sample(p)?.1.sample.temperature_k)
    }
}
impl TemperatureGradientField for MeshTemperature<'_> {
    fn temperature_gradient(
        &self,
        p: Point3,
    ) -> Result<TemperatureGradientSample, GradientIndexError> {
        Ok(self.sample(p)?.1.sample)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum P1MeshError {
    Temperature(TransientTemperatureFieldError),
    Material(ThermoOpticFieldError),
    Propagation(GradientIndexError),
    ForeignSnapshot,
    InvalidTopology { element: usize, local_face: usize },
    FaceBudgetExceeded,
}

impl From<TransientTemperatureFieldError> for P1MeshError {
    fn from(error: TransientTemperatureFieldError) -> Self {
        Self::Temperature(error)
    }
}
impl From<ThermoOpticFieldError> for P1MeshError {
    fn from(error: ThermoOpticFieldError) -> Self {
        Self::Material(error)
    }
}
impl From<GradientIndexError> for P1MeshError {
    fn from(error: GradientIndexError) -> Self {
        Self::Propagation(error)
    }
}
impl core::fmt::Display for P1MeshError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "transient P1 mesh refraction: {self:?}")
    }
}
impl core::error::Error for P1MeshError {}

#[cfg(test)]
mod tests;
