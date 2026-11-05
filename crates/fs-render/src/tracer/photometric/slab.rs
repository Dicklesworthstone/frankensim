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
    log_power_transmission: f64,
    pub(crate) unoccluded: bool,
}

impl SlabTransfer {
    pub(crate) fn irradiance(&self, intensity: f64) -> Result<f64, PhotometryError> {
        if !intensity.is_finite() || intensity < 0.0 {
            return Err(PhotometryError::NonFiniteResult);
        }
        if intensity == 0.0 || self.log_power_transmission == f64::NEG_INFINITY {
            return Ok(0.0);
        }
        let spreading = intensity / self.area_per_solid_angle;
        let ordinary = spreading * self.power_transmission;
        // Retain the established arithmetic in its ordinary range. Resort to
        // log scaling only when an intermediate would lose representable light.
        let value = if spreading.is_finite()
            && spreading >= f64::MIN_POSITIVE
            && self.power_transmission >= f64::MIN_POSITIVE
            && ordinary.is_finite()
            && ordinary > 0.0
        {
            ordinary
        } else {
            det::exp(
                det::ln(intensity) - det::ln(self.area_per_solid_angle)
                    + self.log_power_transmission,
            )
        };
        if value.is_finite() && value >= 0.0 {
            Ok(value)
        } else {
            Err(PhotometryError::NonFiniteResult)
        }
    }
}

/// One unfolded inverse solve; visibility is always evaluated after folding.
pub(crate) struct SlabImage {
    wavelength: f64,
    round_trips: u32,
    index: f64,
    path: SlabConnectionGeometry,
    receiver_cosine: f64,
    internal_depth_m: f64,
    distance_m: f64,
}

fn log_sum_exp(a: f64, b: f64) -> f64 {
    let largest = a.max(b);
    if largest.is_infinite() {
        largest
    } else {
        largest + det::ln(det::exp(a - largest) + det::exp(b - largest))
    }
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
        let image = self.image(cx, wavelength, 0)?;
        // Keep the original selected-image admission. The history owner instead
        // treats back-facing/grazing orders as zero and continues later orders.
        if image.receiver_cosine <= 1.0e-12 {
            return Err(slab_nee_refusal(
                self.slab.boundary_primitive,
                "selected slab image reaches the back or grazing edge of the receiver",
            )
            .into());
        }
        self.sample_image(cx, &image)
    }

    pub(crate) fn image(
        &self,
        cx: &Cx<'_>,
        wavelength: f64,
        round_trips: u32,
    ) -> Result<SlabImage, PhotometryError> {
        cx.checkpoint()?;
        let slab = self.slab;
        let refusal = |reason| slab_nee_refusal(slab.boundary_primitive, reason);
        let index = slab
            .glass
            .ior()
            .eval(wavelength)
            .map_err(TracerError::from)?;
        // Unfold each complete internal return into two extra slab widths.
        // Translate the receiver too, so the exterior normal distance remains
        // the physical one. Translation leaves the receiving Jacobian intact.
        let extra_depth = 2.0 * f64::from(round_trips) * slab.thickness_m;
        let internal_depth_m = slab.thickness_m + extra_depth;
        let path = slab_connection_to_point(
            self.source,
            self.receiver
                .position_m()
                .offset(slab.axis.scale(extra_depth)),
            slab.axis,
            slab.axis,
            internal_depth_m,
            1.0 / index,
        )
        .map_err(refusal)?;
        let receiver_cosine = -self.receiver.normal().dot(path.incident_direction);
        let cosine = path.incident_direction.dot(slab.axis);
        let inside_cosine = path.internal_direction.dot(slab.axis);
        let normal_distance = self
            .receiver
            .position_m()
            .delta_from(self.source)
            .dot(slab.axis);
        let distance_m =
            (normal_distance - slab.thickness_m) / cosine + internal_depth_m / inside_cosine;
        if !distance_m.is_finite() {
            return Err(PhotometryError::NonFiniteResult);
        }
        Ok(SlabImage {
            wavelength,
            round_trips,
            index,
            path,
            receiver_cosine,
            internal_depth_m,
            distance_m,
        })
    }

    /// Analytic omitted-image bound, not a floating-point error enclosure.
    pub(crate) fn omitted_bound(
        &self,
        cx: &Cx<'_>,
        image: &SlabImage,
        peak_intensity: f64,
    ) -> Result<f64, PhotometryError> {
        cx.checkpoint()?;
        if peak_intensity == 0.0 {
            return Ok(0.0);
        }
        let exterior = self
            .receiver
            .position_m()
            .delta_from(self.source)
            .dot(self.slab.axis)
            - self.slab.thickness_m;
        let cosine = image.path.incident_direction.dot(self.slab.axis);
        let width = exterior + image.internal_depth_m / image.index;
        let increment = 2.0 * self.slab.thickness_m / image.index;
        let extinction = self
            .slab
            .glass
            .absorption_ref()
            .extinction_per_m(image.wavelength)
            .map_err(TracerError::from)?;
        // Keep the source scale attached until the last exponential. Tiny
        // reflection/Beer factors can still yield a finite bound for a bright
        // source; multiplying already-underflowed factors would lose it.
        let log_scale = det::ln(peak_intensity) - extinction * image.internal_depth_m;
        // Future launch angles decrease. Their two transverse ray-differential
        // widths give J_m >= cos(theta_N)^2 [H + (2m+1)h/n]^2, including a
        // tilted receiver. Passivity and a decreasing-series integral then
        // bound arbitrary coatings without sampling an angular maximum.
        let mut log_bound = if increment > 0.0 {
            let first = -2.0 * det::ln(width);
            let second = -det::ln(width) - det::ln(increment);
            log_scale - 2.0 * det::ln(cosine) + log_sum_exp(first, second)
        } else {
            f64::INFINITY
        };
        let surface = dielectric_surface(self.scene, self.slab.boundary_primitive)?;
        if surface.coating().is_none() {
            let optics = surface
                .interface(image.wavelength, false)
                .map_err(TracerError::from)?
                .fresnel_polarized(cosine.clamp(0.0, 1.0), 1.0, image.index)
                .map_err(TracerError::from)?;
            let pure_s = Stokes {
                i: 1.0,
                q: 1.0,
                u: 0.0,
                v: 0.0,
            };
            let reflection = optics.reflection.apply(pure_s).i.clamp(0.0, 1.0);
            let transmission = optics.transmission.apply(pure_s).i.clamp(0.0, 1.0);
            if reflection == 0.0 && image.round_trips != 0 {
                return Ok(0.0);
            }
            // Bare Rs grows with angle and bounds Rp. Every future history
            // therefore has at most Rs(theta_N)^(2m) reflected power. Native
            // Ts retains 1-Rs^2 even when Rs rounds to one near grazing.
            let denominator = transmission * (2.0 - transmission);
            let log_power = if image.round_trips == 0 {
                0.0
            } else {
                2.0 * f64::from(image.round_trips) * det::ln(reflection)
            };
            let minimum_width = exterior.max(cosine * width);
            let bare = log_scale - 2.0 * det::ln(minimum_width) + log_power - det::ln(denominator);
            log_bound = log_bound.min(bare);
        }
        if log_bound.is_nan() {
            return Err(PhotometryError::NonFiniteResult);
        }
        // Never report a physically positive tail as exactly zero just because
        // its ordinary floating-point evaluation underflowed.
        Ok(det::exp(log_bound).max(f64::from_bits(1)))
    }

    pub(crate) fn sample_image(
        &self,
        cx: &Cx<'_>,
        image: &SlabImage,
    ) -> Result<SlabTransfer, PhotometryError> {
        cx.checkpoint()?;
        let slab = self.slab;
        let refusal = |reason| slab_nee_refusal(slab.boundary_primitive, reason);
        let SlabImage {
            wavelength,
            index,
            path,
            distance_m: distance,
            ..
        } = *image;
        let cosine = path.incident_direction.dot(slab.axis);
        let inside_cosine = path.internal_direction.dot(slab.axis);
        let mut result = SlabTransfer {
            launch_direction: path.incident_direction,
            area_per_solid_angle: 1.0,
            power_transmission: 0.0,
            log_power_transmission: f64::NEG_INFINITY,
            unoccluded: false,
        };
        if image.receiver_cosine <= 0.0 {
            return Ok(result);
        }
        result.area_per_solid_angle = path
            .light_area_jacobian
            .ok_or(PhotometryError::NonFiniteResult)?
            * cosine
            / image.receiver_cosine;
        if !result.area_per_solid_angle.is_finite() || result.area_per_solid_angle <= 0.0 {
            return Err(PhotometryError::NonFiniteResult);
        }
        if distance > self.options.max_distance_m {
            return Err(PhotometryError::DistanceLimitExceeded);
        }
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
        let mut inside_ray = Ray {
            origin: entry.spawn_boundary_ray(entry_frame.geometric, internal_direction),
            dir: internal_direction,
        };
        let exit_optics = surface
            .interface(wavelength, true)
            .map_err(TracerError::from)?
            .fresnel_polarized(inside_cosine.clamp(0.0, 1.0), index, 1.0)
            .map_err(TracerError::from)?;
        if exit_optics.scalar.total_internal_reflection {
            return Err(refusal("selected slab image has no propagating exit").into());
        }
        let mut state = entry_optics.transmission.apply(Stokes::unpolarized(1.0));
        let mut previous_point = entry.hit.point;
        let mut internal_distance = 0.0;
        let final_interface = 2 * u64::from(image.round_trips);
        let mut interaction = 0;
        let (exit, exit_frame) = loop {
            cx.checkpoint()?;
            let Some(hit) = trace(&inside_ray, distance)? else {
                return Err(refusal("selected image escapes before the slab exit").into());
            };
            if hit.primitive_index != slab.boundary_primitive {
                if primitive_is_dielectric(self.scene, hit.primitive_index) {
                    return Err(refusal("nested dielectric interrupts the slab image").into());
                }
                return Ok(result);
            }
            let frame = surface_frame(&hit.hit, &inside_ray)?;
            let far_face = interaction % 2 == 0;
            let (reference, sign) = if far_face {
                (slab.exit_reference, 1.0)
            } else {
                (slab.entry_reference, -1.0)
            };
            if frame.entering
                || sign * frame.geometric.dot(slab.axis) < 1.0 - SLAB_PARALLEL_COSINE_TOLERANCE
                || hit.hit.point.delta_from(reference).dot(slab.axis).abs() > tolerance
            {
                return Err(
                    refusal("selected history misses its original parallel broad face").into(),
                );
            }
            internal_distance += hit.hit.point.delta_from(previous_point).norm();
            if !internal_distance.is_finite() {
                return Err(PhotometryError::NonFiniteResult);
            }
            if interaction == final_interface {
                break (hit, frame);
            }
            // The original surface is seen from inside at both return faces.
            // Its reversed coating, s/p powers and native specular direction
            // remain bound to each physical interaction, never the unfolded ray.
            state = exit_optics.reflection.apply(state);
            let reflected = crate::dielectric::reflect(inside_ray.dir.scale(-1.0), frame.geometric);
            inside_ray = Ray {
                origin: hit.spawn_boundary_ray(frame.geometric, reflected),
                dir: reflected,
            };
            previous_point = hit.hit.point;
            interaction += 1;
        };
        let outgoing = crate::dielectric::refract(
            inside_ray.dir.scale(-1.0),
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
        let physical_distance =
            entry.hit.point.delta_from(self.source).norm() + internal_distance + to_receiver.norm();
        if !physical_distance.is_finite() || physical_distance > self.options.max_distance_m {
            return Err(PhotometryError::DistanceLimitExceeded);
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
        // Parallel interfaces share s/p. Unpolarized entry can produce Q but
        // no U/V; all repeated reflection and exit transmission retain it.
        let power = exit_optics.transmission.apply(state).i;
        let beer = medium_transmittance(Some(slab.glass), wavelength, internal_distance)?;
        result.power_transmission = power * beer;
        if !result.power_transmission.is_finite() || result.power_transmission < 0.0 {
            return Err(PhotometryError::NonFiniteResult);
        }
        if result.power_transmission >= f64::MIN_POSITIVE {
            result.log_power_transmission = det::ln(result.power_transmission);
            result.unoccluded = true;
            return Ok(result);
        }
        // Preserve positive throughput that only becomes representable after
        // multiplying a bright calibrated source. The same common-plane s/p
        // history permits a closed log product of the native power matrices.
        let channel = |q| {
            let polarization = Stokes {
                i: 1.0,
                q,
                u: 0.0,
                v: 0.0,
            };
            let entry = entry_optics.transmission.apply(polarization).i;
            let exit = exit_optics.transmission.apply(polarization).i;
            let reflection = exit_optics.reflection.apply(polarization).i;
            if entry == 0.0 || exit == 0.0 || (image.round_trips > 0 && reflection == 0.0) {
                return f64::NEG_INFINITY;
            }
            let returns = if image.round_trips == 0 {
                0.0
            } else {
                2.0 * f64::from(image.round_trips) * det::ln(reflection)
            };
            det::ln(entry) + det::ln(exit) + returns
        };
        let extinction = slab
            .glass
            .absorption_ref()
            .extinction_per_m(wavelength)
            .map_err(TracerError::from)?;
        result.log_power_transmission = log_sum_exp(channel(1.0), channel(-1.0))
            - core::f64::consts::LN_2
            - extinction * internal_distance;
        if result.log_power_transmission.is_nan() {
            return Err(PhotometryError::NonFiniteResult);
        }
        result.unoccluded = true;
        Ok(result)
    }
}
