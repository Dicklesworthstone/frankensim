//! Ordinary-scene bridge to the existing two-interface Snell connector.
//! This owner supplies original geometry/visibility and forward s/p power;
//! source calibration and spectral quadrature stay with photometry.

use super::*;
use crate::polarization::Stokes;

pub(crate) struct SlabConnection<'a> {
    scene: &'a Scene,
    slab: ParallelSlab,
    geometry: geometry::PreparedGeometry,
    source: Point3,
    receiver: IlluminanceProbe,
    options: TransportOptions,
}

pub(crate) struct SlabTransfer {
    pub(crate) launch_direction: Vec3,
    pub(crate) area_per_solid_angle: f64,
    pub(crate) power_transmission: f64,
    pub(crate) unoccluded: bool,
}

impl<'a> SlabConnection<'a> {
    pub(crate) fn try_new(
        cx: &Cx<'_>,
        scene: &'a Scene,
        primitive: usize,
        source: Point3,
        receiver: IlluminanceProbe,
        options: TransportOptions,
    ) -> Result<Self, PhotometryError> {
        cx.checkpoint()?;
        if !options.ray_epsilon_m.is_finite()
            || options.ray_epsilon_m <= 0.0
            || !options.max_distance_m.is_finite()
            || options.max_distance_m <= options.ray_epsilon_m
            || options.max_outputs == 0
        {
            return Err(PhotometryError::InvalidOptions);
        }
        for (index, object) in scene.primitives.iter().enumerate() {
            cx.checkpoint()?;
            if matches!(object.shape, Shape::AnimatedInstance(_)) {
                return Err(PhotometryError::UnsupportedPrimitive {
                    primitive: index,
                    reason: "slab photometry requires a static geometry snapshot",
                });
            }
        }
        let boundary = scene
            .primitives
            .get(primitive)
            .ok_or(PhotometryError::InvalidPrimitive { primitive })?;
        let delta = receiver.position_m().delta_from(source);
        let distance = norm(delta)?;
        if distance <= options.ray_epsilon_m {
            return Err(PhotometryError::SourceTooClose { source: 0 });
        }
        if distance > options.max_distance_m {
            return Err(PhotometryError::DistanceLimitExceeded);
        }
        let media = initial_medium::InitialMedia::try_new(scene, cx)?;
        for endpoint in [source, receiver.position_m()] {
            if media.stack(scene, cx, endpoint, None)?.len() != 0 {
                return Err(slab_nee_refusal(
                    primitive,
                    "source and receiver must both be in vacuum",
                )
                .into());
            }
        }
        // Select the authored glass, not an opaque blocker on the unbent chord.
        // Only the subsequently solved physical segments determine visibility.
        let chord = Ray {
            origin: source,
            dir: delta.scale(1.0 / distance),
        };
        let first =
            intersect_primitive(primitive, boundary, cx, &chord, None, Some(distance), None)?
                .ok_or_else(|| {
                    slab_nee_refusal(primitive, "source/receiver chord misses the selected slab")
                })?;
        let slab = discover_parallel_slab_with_probe(scene, &chord, first, |probe| {
            intersect_primitive(
                primitive,
                boundary,
                cx,
                probe,
                None,
                Some(options.max_distance_m),
                None,
            )
        })?;
        // Endpoint stacks skip bodies whose AABBs cannot contain the endpoint.
        // Query the selected body's actual interior so two open parallel sheets
        // cannot acquire a homogeneous-glass volume from their face geometry.
        let interior = slab
            .entry_reference
            .offset(slab.axis.scale(0.5 * slab.thickness_m));
        if !media.boundary_contains(cx, primitive, interior, None)? {
            return Err(slab_nee_refusal(
                primitive,
                "selected slab has no admitted closed interior",
            )
            .into());
        }
        if slab.entry_reference.delta_from(source).dot(slab.axis) <= 0.0
            || receiver
                .position_m()
                .delta_from(slab.exit_reference)
                .dot(slab.axis)
                <= 0.0
        {
            return Err(
                slab_nee_refusal(primitive, "source and receiver do not bracket the slab").into(),
            );
        }
        Ok(Self {
            scene,
            slab,
            geometry: geometry::PreparedGeometry::try_new(&scene.primitives, cx)?,
            source,
            receiver,
            options,
        })
    }

    pub(crate) fn sample(
        &self,
        cx: &Cx<'_>,
        wavelength: f64,
    ) -> Result<SlabTransfer, PhotometryError> {
        cx.checkpoint()?;
        let slab = self.slab;
        let refusal = |reason| slab_nee_refusal(slab.boundary_primitive, reason);
        let index = slab
            .glass
            .ior()
            .eval(wavelength)
            .map_err(TracerError::from)?;
        let path = slab_connection_to_point(
            self.source,
            self.receiver.position_m(),
            self.receiver.normal(),
            slab.axis,
            slab.thickness_m,
            1.0 / index,
        )
        .map_err(refusal)?;
        if self.receiver.normal().dot(path.incident_direction) >= 0.0 {
            return Err(refusal("selected slab image reaches the back of the receiver").into());
        }
        let cosine = path.incident_direction.dot(slab.axis);
        let inside_cosine = path.internal_direction.dot(slab.axis);
        let normal_distance = self
            .receiver
            .position_m()
            .delta_from(self.source)
            .dot(slab.axis);
        let distance =
            (normal_distance - slab.thickness_m) / cosine + slab.thickness_m / inside_cosine;
        if !distance.is_finite() || distance > self.options.max_distance_m {
            return Err(PhotometryError::DistanceLimitExceeded);
        }
        let mut result = SlabTransfer {
            launch_direction: path.incident_direction,
            area_per_solid_angle: path
                .light_area_jacobian
                .ok_or(PhotometryError::NonFiniteResult)?,
            power_transmission: 0.0,
            unoccluded: false,
        };
        let trace = |ray: &Ray, limit| {
            self.geometry
                .intersect_scene(&self.scene.primitives, cx, ray, limit)
        };
        let entry_ray = Ray {
            origin: self.source,
            dir: path.incident_direction,
        };
        let Some(entry) = trace(&entry_ray, distance)? else {
            return Err(refusal("selected refracted entry misses the finite slab").into());
        };
        if entry.primitive_index != slab.boundary_primitive {
            if primitive_is_dielectric(self.scene, entry.primitive_index) {
                return Err(
                    refusal("another dielectric interrupts the selected slab image").into(),
                );
            }
            return Ok(result);
        }
        let entry_frame = surface_frame(&entry.hit, &entry_ray)?;
        let tolerance = SLAB_CONNECTION_REL_TOLERANCE * slab.thickness_m.max(RAY_EPS);
        if !entry_frame.entering
            || entry_frame.geometric.dot(slab.axis) > -1.0 + SLAB_PARALLEL_COSINE_TOLERANCE
            || entry
                .hit
                .point
                .delta_from(slab.entry_reference)
                .dot(slab.axis)
                .abs()
                > tolerance
        {
            return Err(
                refusal("selected image does not enter the original parallel broad face").into(),
            );
        }
        let surface = dielectric_surface(self.scene, slab.boundary_primitive)?;
        let entry_optics = surface
            .interface(wavelength, false)
            .map_err(TracerError::from)?
            .fresnel_polarized(cosine.clamp(0.0, 1.0), 1.0, index)
            .map_err(TracerError::from)?;
        let internal_direction = crate::dielectric::refract(
            path.incident_direction.scale(-1.0),
            entry_frame.geometric,
            1.0,
            index,
            entry_optics.scalar.transmitted_cosine,
        )
        .map_err(TracerError::from)?;
        if !parallel_directions(internal_direction, path.internal_direction) {
            return Err(
                refusal("native entry refraction disagrees with the slab connection").into(),
            );
        }
        let inside_ray = Ray {
            origin: entry.spawn_boundary_ray(entry_frame.geometric, internal_direction),
            dir: internal_direction,
        };
        let Some(exit) = trace(&inside_ray, distance)? else {
            return Err(refusal("selected image escapes before the slab exit").into());
        };
        if exit.primitive_index != slab.boundary_primitive {
            if primitive_is_dielectric(self.scene, exit.primitive_index) {
                return Err(refusal("nested dielectric interrupts the slab image").into());
            }
            return Ok(result);
        }
        let exit_frame = surface_frame(&exit.hit, &inside_ray)?;
        if exit_frame.entering
            || exit_frame.geometric.dot(slab.axis) < 1.0 - SLAB_PARALLEL_COSINE_TOLERANCE
            || exit
                .hit
                .point
                .delta_from(slab.exit_reference)
                .dot(slab.axis)
                .abs()
                > tolerance
        {
            return Err(
                refusal("selected image does not exit the original parallel broad face").into(),
            );
        }
        let exit_optics = surface
            .interface(wavelength, true)
            .map_err(TracerError::from)?
            .fresnel_polarized(inside_cosine.clamp(0.0, 1.0), index, 1.0)
            .map_err(TracerError::from)?;
        if exit_optics.scalar.total_internal_reflection {
            return Err(refusal("selected slab image has no propagating exit").into());
        }
        let outgoing = crate::dielectric::refract(
            internal_direction.scale(-1.0),
            exit_frame.geometric.scale(-1.0),
            index,
            1.0,
            exit_optics.scalar.transmitted_cosine,
        )
        .map_err(TracerError::from)?;
        if !parallel_directions(outgoing, path.incident_direction) {
            return Err(refusal("native exit does not restore the exterior direction").into());
        }
        let to_receiver = self.receiver.position_m().delta_from(exit.hit.point);
        let remaining = to_receiver.dot(outgoing);
        let residual = vec_sub(to_receiver, outgoing.scale(remaining)).norm();
        if remaining <= self.options.ray_epsilon_m
            || !residual.is_finite()
            || residual > SLAB_CONNECTION_REL_TOLERANCE * distance.max(RAY_EPS)
        {
            return Err(
                refusal("native refracted path does not reach the selected receiver").into(),
            );
        }
        let exit_ray = Ray {
            origin: exit.spawn_boundary_ray(exit_frame.geometric, outgoing),
            dir: outgoing,
        };
        let limit = self
            .receiver
            .position_m()
            .delta_from(exit_ray.origin)
            .dot(outgoing)
            - self.options.ray_epsilon_m;
        if limit <= 0.0 {
            return Err(refusal("receiver lies within the exit boundary exclusion").into());
        }
        if let Some(blocker) = trace(&exit_ray, limit)? {
            if primitive_is_dielectric(self.scene, blocker.primitive_index) {
                return Err(
                    refusal("another dielectric interrupts the outgoing slab image").into(),
                );
            }
            return Ok(result);
        }
        // Parallel interfaces share the s/p plane. Preserve the first face's
        // induced Q through the second, instead of multiplying mean T values.
        let power = exit_optics
            .transmission
            .apply(entry_optics.transmission.apply(Stokes::unpolarized(1.0)))
            .i;
        let beer = medium_transmittance(
            Some(slab.glass),
            wavelength,
            exit.hit.point.delta_from(entry.hit.point).norm(),
        )?;
        result.power_transmission = power * beer;
        if !result.power_transmission.is_finite() || result.power_transmission < 0.0 {
            return Err(PhotometryError::NonFiniteResult);
        }
        result.unoccluded = true;
        Ok(result)
    }
}
