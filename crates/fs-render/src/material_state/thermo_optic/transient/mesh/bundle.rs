//! Ray bundles through original continuous-temperature element faces.
//!
//! The original launch-angle columns and native optical integrals are retained
//! across the whole mesh. Each gradient jump uses a temporary coordinate normal
//! to its original face, with the two derivatives of physical arc length carried
//! through the same adaptive stages. Returning to fixed total arc length gives
//! the crossing-time (saltation) correction. No neighboring-ray trace, central
//! pretrace, temperature smoothing or zero-Hessian substitution is performed.
//!
//! As in the primal owner, vertex IDs determine adjacency, one material law
//! covers the frozen temperature field, and contacts/exterior faces stop. The
//! finite departure neighborhood and transverse admission apply throughout the
//! collar. A roundoff-unresolved endpoint on a genuine gradient-jump plane
//! refuses a derivative; an endpoint resolved on either side remains usable.
//! Receiving-plane caustic/irradiance admission stays with the native bundle.

use super::*;
use crate::gradient_index::bundle::event::{
    BundleEventOptions, propagate_seeded_bundle_across_continuous_plane_with_control,
    propagate_seeded_bundle_to_boundary_with_control,
};
use crate::gradient_index::bundle::{
    BundleOptions, RayBundlePropagation, SecondOrderIndexField, propagate_bundle_with_control,
};

/// Whole-mesh geometry/work controls plus native angular-column tolerances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct P1MeshBundleOptions {
    pub mesh: P1MeshOptions,
    /// Metres/radian, also controlling the face-coordinate arc-length columns.
    pub position_derivative_tolerance_m: f64,
    pub direction_derivative_tolerance: f64,
}
impl Default for P1MeshBundleOptions {
    fn default() -> Self {
        let bundle = BundleOptions::default();
        Self {
            mesh: P1MeshOptions::default(),
            position_derivative_tolerance_m: bundle.position_derivative_tolerance_m,
            direction_derivative_tolerance: bundle.direction_derivative_tolerance,
        }
    }
}

/// Original launch-angle family, with derivatives at equal TOTAL physical arc.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct P1MeshBundlePropagation {
    pub bundle: RayBundlePropagation,
    pub termination: BoundaryTermination,
    pub element_index: usize,
    /// Internal coordinate continuations; a partial collar need not have
    /// reached the mathematical face. These introduce no Fresnel event weight.
    pub internal_faces: u32,
    pub boundary_evaluations: u32,
    pub contact_resolution_m: f64,
}

impl<'a> P1ConductionMesh<'a> {
    /// Propagate two original orthonormal launch-angle variations through the
    /// same frozen mesh as [`Self::propagate`]. Every Hessian comes from the
    /// original one-sided P1 reconstruction and native material chain rule.
    /// A changed normal gradient has its crossing-time correction, even when
    /// both one-sided temperature Hessians are zero.
    #[allow(clippy::too_many_arguments)]
    pub fn propagate_bundle(
        &self,
        cx: &Cx<'_>,
        binding: &MaterialThermoOpticField<'_, P1ConductionElement<'a>>,
        ray: Ray,
        transverse_axis: Vec3,
        wavelength_nm: f64,
        max_distance_m: f64,
        options: P1MeshBundleOptions,
    ) -> Result<P1MeshBundlePropagation, P1MeshError> {
        options.mesh.event.validate()?;
        if !max_distance_m.is_finite()
            || max_distance_m < 0.0
            || options.mesh.max_internal_faces > MAX_INTERNAL_FACES
        {
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
        self.propagate_bundle_prepared(
            cx, &cursor, &prepared, ray, transverse_axis, wavelength_nm, max_distance_m, options,
        )
    }

    // The original variational/face-coordinate transport is shared by glass
    // and water; each prepared law supplies its complete nonlinear index jet.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn propagate_bundle_prepared<F: SecondOrderIndexField>(
        &self,
        cx: &Cx<'_>,
        cursor: &MeshTemperature<'a>,
        prepared: &F,
        ray: Ray,
        transverse_axis: Vec3,
        wavelength_nm: f64,
        max_distance_m: f64,
        options: P1MeshBundleOptions,
    ) -> Result<P1MeshBundlePropagation, P1MeshError> {
        let start = cursor.current.get();
        // The native zero-length launch establishes the immutable angular basis
        // and validates the source/ray. Every later call retains those columns.
        let mut bundle = propagate_bundle_with_control(
            prepared,
            ray,
            transverse_axis,
            wavelength_nm,
            0.0,
            BundleOptions {
                integration: options.mesh.event.integration,
                position_derivative_tolerance_m: options.position_derivative_tolerance_m,
                direction_derivative_tolerance: options.direction_derivative_tolerance,
            },
            || cx.checkpoint().is_ok(),
        )?;
        let mut boundary = MeshBoundary::new(start.boundary_oracle(
            options.mesh.boundary_id,
            options.mesh.event.surface_tolerance_m,
        )?);
        let mut departure = None;
        let mut distance = Dd::ZERO;
        let mut internal_faces = 0;
        let mut boundary_evaluations = 0_u32;
        let mut contact_resolution_m = 0.0_f64;
        loop {
            cx.checkpoint().map_err(|_| GradientIndexError::Cancelled)?;
            let past = *bundle.propagation();
            let mut event = options.mesh.event;
            event.integration.max_steps = event
                .integration
                .max_steps
                .checked_sub(past.accepted_steps + past.rejected_steps)
                .filter(|&n| n > 0)
                .ok_or(GradientIndexError::StepBudgetExceeded)?;
            event.max_boundary_queries = event
                .max_boundary_queries
                .checked_sub(boundary_evaluations)
                .filter(|&n| n > 0)
                .ok_or(GradientIndexError::BoundaryQueryBudgetExceeded)?;
            let native = BundleEventOptions {
                boundary: event,
                position_derivative_tolerance_m: options.position_derivative_tolerance_m,
                direction_derivative_tolerance: options.direction_derivative_tolerance,
            };
            let remaining = (Dd::from_f64(max_distance_m) - distance).to_f64();
            let seed = bundle.continuation();
            let (result, in_coordinate_collar) = if let Some(departure) = departure {
                propagate_seeded_bundle_across_continuous_plane_with_control(
                    prepared,
                    &mut boundary,
                    &seed,
                    remaining,
                    departure,
                    native,
                    || cx.checkpoint().is_ok(),
                )?
            } else {
                (
                    propagate_seeded_bundle_to_boundary_with_control(
                        prepared,
                        &mut boundary,
                        &seed,
                        remaining,
                        None,
                        native,
                        || cx.checkpoint().is_ok(),
                    )?,
                    false,
                )
            };
            distance = distance + Dd::from_f64(result.segment.distance_m);
            if !distance.hi.is_finite() || !distance.lo.is_finite() {
                return Err(GradientIndexError::NumericalRange.into());
            }
            bundle = result.bundle;
            let p = *bundle.propagation();
            boundary_evaluations += result.boundary_evaluations;
            contact_resolution_m = contact_resolution_m.max(result.contact_resolution_m);
            if result.termination == BoundaryTermination::DistanceLimit
                || distance.to_f64() >= max_distance_m
            {
                if in_coordinate_collar {
                    cx.checkpoint().map_err(|_| GradientIndexError::Cancelled)?;
                    if boundary_evaluations >= options.mesh.event.max_boundary_queries {
                        return Err(GradientIndexError::BoundaryQueryBudgetExceeded.into());
                    }
                    boundary_evaluations += 1;
                    boundary.admit_bundle_endpoint(p.ray.origin)?;
                }
                return Ok(P1MeshBundlePropagation {
                    bundle,
                    termination: BoundaryTermination::DistanceLimit,
                    element_index: cursor.sample(p.ray.origin)?.0,
                    internal_faces,
                    boundary_evaluations,
                    contact_resolution_m,
                });
            }
            let BoundaryTermination::Surface(mut witness) = result.termination else {
                unreachable!("all native boundary terminations were handled")
            };
            let active = cursor.current.get();
            let local_face = witness.feature_id;
            let Some(next) = self.neighbours[active.element][local_face] else {
                witness.feature_id = self.face_ids[active.element][local_face];
                return Ok(P1MeshBundlePropagation {
                    bundle,
                    termination: BoundaryTermination::Surface(witness),
                    element_index: active.element,
                    internal_faces,
                    boundary_evaluations,
                    contact_resolution_m,
                });
            };
            if internal_faces >= options.mesh.max_internal_faces {
                return Err(P1MeshError::FaceBudgetExceeded);
            }
            let neighbour = self.snapshot.element(next.element)?;
            let transition =
                Transition::new(active, neighbour, next.local_face, p.ray.origin, witness);
            departure = Some(transition.departure());
            boundary = MeshBoundary::with_transition(
                neighbour.boundary_oracle(
                    options.mesh.boundary_id,
                    options.mesh.event.surface_tolerance_m,
                )?,
                transition,
            );
            cursor.previous.set(Some(active));
            cursor.current.set(neighbour);
            internal_faces += 1;
        }
    }
}

#[cfg(test)]
mod tests;
