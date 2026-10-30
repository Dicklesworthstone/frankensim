//! Calibrated mesh-area proposals and their complementary source-hit weights.
//!
//! Geometry, BSDFs, phase laws, endpoint visibility and physical transmittance
//! follow the shared tracer. The proposal retains actual emitter points: a
//! hidden face on the same primitive is not another density for a visible hit.

use super::*;
use crate::photometry::area::{EmissionSidedness, PhotometricAreaEmitter};
use crate::photometry::source::SourceError;
use core::ops::Range;

const AREA_NEE_DOMAIN: u32 = 0x6172_6561;
const PROBE_AREA_NEE_DOMAIN: u32 = 0x6172_7072;

struct Emitter<'a> {
    source: &'a PhotometricAreaEmitter,
    faces: Range<usize>,
    luminance_cd_m2: f64,
    projector: Option<ProjectorPlacement>,
    probability: f64,
    cumulative_probability: f64,
}

struct ProjectorPlacement {
    world_to_optical: crate::instances::RigidTransform,
    clear_flux_per_area_lm_m2: f64,
}

struct Face {
    vertices: [Point3; 3],
    normal: Vec3,
    area_m2: f64,
    probability: f64,
    cumulative_probability: f64,
}

/// Admitted immutable emitters, retaining native world-space triangle witnesses.
pub(in crate::tracer) struct AreaLighting<'a> {
    emitters: Vec<Emitter<'a>>,
    faces: Vec<Face>,
    // Sorted primitive indices keep source-hit lookup independent of mesh size.
    primitive_emitters: Vec<(usize, usize)>,
    active: bool,
}

#[derive(Clone, Copy)]
struct AreaSample {
    emitter: usize,
    point: Point3,
    pdf_solid_angle: f64,
}

/// The actual receiver response and continuation density at one path vertex.
#[derive(Clone, Copy)]
pub(in crate::tracer) enum Receiver<'a> {
    Surface {
        material: &'a Material,
        frame: SurfaceFrame,
        outgoing: Vec3,
        incident_medium: Option<DielectricGlass>,
        boundary: Option<BoundaryMedia>,
        brushed: Option<&'a crate::surface_detail::ShadingFrame>,
        intersection: &'a SceneIntersection,
        incoming_ray: &'a Ray,
    },
    Phase {
        phase: PhaseFunction,
        incoming_travel: Vec3,
    },
    Probe {
        normal: Vec3,
    },
}

impl<'a> AreaLighting<'a> {
    pub(in crate::tracer) const fn empty() -> Self {
        Self {
            emitters: Vec::new(),
            faces: Vec::new(),
            primitive_emitters: Vec::new(),
            active: false,
        }
    }

    pub(in crate::tracer) fn try_new(
        scene: &Scene,
        sources: &'a [PhotometricAreaEmitter],
        cx: &Cx<'_>,
    ) -> Result<Self, PhotometryError> {
        cx.checkpoint()?;
        let mut result = Self::empty();
        result
            .emitters
            .try_reserve_exact(sources.len())
            .map_err(|_| PhotometryError::AllocationRefused)?;
        result
            .primitive_emitters
            .try_reserve_exact(sources.len())
            .map_err(|_| PhotometryError::AllocationRefused)?;
        let maximum_flux = sources
            .iter()
            .map(PhotometricAreaEmitter::flux_lm)
            .fold(0.0, f64::max);
        let mut flux_sum = 0.0;
        for (emitter_index, source) in sources.iter().enumerate() {
            cx.checkpoint()?;
            let primitive_index = source.primitive_index();
            let primitive =
                scene
                    .primitives
                    .get(primitive_index)
                    .ok_or(PhotometryError::InvalidPrimitive {
                        primitive: primitive_index,
                    })?;
            if primitive.emission.is_some() {
                return Err(PhotometryError::UnsupportedPrimitive {
                    primitive: primitive_index,
                    reason: "calibrated area emission cannot share legacy relative emission",
                });
            }
            if matches!(primitive.material, Material::Dielectric { .. }) {
                return Err(PhotometryError::UnsupportedPrimitive {
                    primitive: primitive_index,
                    reason: "calibrated surface emission requires an opaque emitting boundary",
                });
            }
            let (mesh, placement) = match &primitive.shape {
                Shape::Mesh(mesh) => (mesh, None),
                Shape::Instance(instance) => match instance.geometry() {
                    SharedGeometry::Mesh(mesh) => (mesh.as_ref(), Some(instance.transform())),
                    SharedGeometry::Chart(_) => {
                        return Err(PhotometryError::UnsupportedPrimitive {
                            primitive: primitive_index,
                            reason: "calibrated area emission requires native mesh faces",
                        });
                    }
                },
                Shape::Chart(_) | Shape::AnimatedInstance(_) => {
                    return Err(PhotometryError::UnsupportedPrimitive {
                        primitive: primitive_index,
                        reason: "calibrated area emission requires static native mesh faces",
                    });
                }
            };
            // Optics live in the primitive's original body frame. Check the
            // pupil in that frame, before a possibly large world translation
            // can hide a nonplanar aperture in coordinate rounding.
            let primitive_to_optical = source
                .projector()
                .map(|projector| projector.frame().inverse())
                .transpose()
                .map_err(|_| PhotometryError::Source(SourceError::NumericalRange))?;
            let world_to_optical = source
                .projector()
                .map(|projector| {
                    placement
                        .unwrap_or_else(crate::instances::RigidTransform::identity)
                        .compose(projector.frame())?
                        .inverse()
                })
                .transpose()
                .map_err(|_| PhotometryError::Source(SourceError::NumericalRange))?;
            let mut pupil_abs_xy_m = [0.0_f64; 2];
            result
                .faces
                .try_reserve(mesh.triangles.len())
                .map_err(|_| PhotometryError::AllocationRefused)?;
            let first_face = result.faces.len();
            let mut area_sum = 0.0;
            let mut triangle_keys = Vec::new();
            triangle_keys
                .try_reserve_exact(mesh.triangles.len())
                .map_err(|_| PhotometryError::AllocationRefused)?;
            for triangle in &mesh.triangles {
                cx.checkpoint()?;
                let mut vertices = [Point3::new(0.0, 0.0, 0.0); 3];
                let mut optical_vertices = vertices;
                for ((point, optical_point), &index) in
                    vertices.iter_mut().zip(&mut optical_vertices).zip(triangle)
                {
                    let vertex = mesh.vertices.get(index as usize).ok_or(
                        PhotometryError::InvalidPrimitive {
                            primitive: primitive_index,
                        },
                    )?;
                    let local = Point3::new(vertex[0], vertex[1], vertex[2]);
                    if let (Some(projector), Some(to_optical)) =
                        (source.projector(), primitive_to_optical)
                    {
                        *optical_point = to_optical.transform_point(local);
                        let roundoff_scale_m = vertex
                            .iter()
                            .chain(projector.frame().translation_m().iter())
                            .map(|value| value.abs())
                            .fold(0.0_f64, f64::max);
                        projector.validate_aperture_point(*optical_point, roundoff_scale_m)?;
                        pupil_abs_xy_m[0] = pupil_abs_xy_m[0].max(optical_point.x.abs());
                        pupil_abs_xy_m[1] = pupil_abs_xy_m[1].max(optical_point.y.abs());
                    }
                    *point = placement.map_or(local, |transform| transform.transform_point(local));
                    if [point.x, point.y, point.z]
                        .iter()
                        .any(|value| !value.is_finite())
                    {
                        return Err(PhotometryError::InvalidPrimitive {
                            primitive: primitive_index,
                        });
                    }
                }
                if let Some(projector) = source.projector() {
                    projector.validate_aperture_normal(cross(
                        optical_vertices[1].delta_from(optical_vertices[0]),
                        optical_vertices[2].delta_from(optical_vertices[0]),
                    ))?;
                }
                triangle_keys.push(triangle_key(vertices));
                let perpendicular = cross(
                    vertices[1].delta_from(vertices[0]),
                    vertices[2].delta_from(vertices[0]),
                );
                let twice_area =
                    norm(perpendicular).map_err(|_| PhotometryError::InvalidPrimitive {
                        primitive: primitive_index,
                    })?;
                let area_m2 = 0.5 * twice_area;
                let next = area_sum + area_m2;
                if !(area_m2 > 0.0 && next.is_finite() && next > area_sum) {
                    return Err(PhotometryError::Source(SourceError::NumericalRange));
                }
                area_sum = next;
                result.faces.push(Face {
                    vertices,
                    normal: perpendicular.scale(1.0 / twice_area),
                    area_m2,
                    probability: 0.0,
                    cumulative_probability: area_sum,
                });
            }
            if area_sum == 0.0 {
                return Err(PhotometryError::InvalidPrimitive {
                    primitive: primitive_index,
                });
            }
            // Coincident duplicate faces would count flux twice while the
            // native first-hit geometry exposes only one of those faces.
            triangle_keys.sort_unstable();
            if triangle_keys.windows(2).any(|pair| pair[0] == pair[1]) {
                return Err(PhotometryError::InvalidPrimitive {
                    primitive: primitive_index,
                });
            }
            let mut previous = 0.0;
            for face in &mut result.faces[first_face..] {
                cx.checkpoint()?;
                face.cumulative_probability /= area_sum;
                face.probability = face.cumulative_probability - previous;
                if !(face.probability > 0.0 && face.probability.is_finite()) {
                    return Err(PhotometryError::Source(SourceError::NumericalRange));
                }
                previous = face.cumulative_probability;
            }
            let sides = match source.sidedness() {
                EmissionSidedness::Front => 1.0,
                EmissionSidedness::Both => 2.0,
            };
            let luminance = source.flux_lm() / area_sum / PI / sides;
            let projector = world_to_optical.map(|world_to_optical| ProjectorPlacement {
                world_to_optical,
                clear_flux_per_area_lm_m2: source.clear_flux_lm() / area_sum,
            });
            let maximum_luminance = if let Some(optics) = source.projector() {
                let scale = source.clear_flux_lm() / area_sum;
                let maximum = scale * optics.maximum_luminance_factor(pupil_abs_xy_m)?;
                if !scale.is_finite()
                    || !maximum.is_finite()
                    || (source.flux_lm() > 0.0 && maximum <= 0.0)
                {
                    return Err(PhotometryError::Source(SourceError::NumericalRange));
                }
                maximum
            } else {
                luminance
            };
            if !luminance.is_finite()
                || (source.flux_lm() > 0.0 && luminance <= 0.0)
                || source
                    .spectrum()
                    .bins_w_per_lumen_nm()
                    .iter()
                    .any(|power| !(maximum_luminance * *power).is_finite())
            {
                return Err(PhotometryError::Source(SourceError::NumericalRange));
            }
            let weight = if maximum_flux > 0.0 {
                source.flux_lm() / maximum_flux
            } else {
                0.0
            };
            let next = flux_sum + weight;
            if !next.is_finite() || (source.flux_lm() > 0.0 && next <= flux_sum) {
                return Err(PhotometryError::Source(SourceError::NumericalRange));
            }
            flux_sum = next;
            result.emitters.push(Emitter {
                source,
                faces: first_face..result.faces.len(),
                luminance_cd_m2: luminance,
                projector,
                probability: 0.0,
                cumulative_probability: flux_sum,
            });
            result
                .primitive_emitters
                .push((primitive_index, emitter_index));
        }
        result
            .primitive_emitters
            .sort_unstable_by_key(|(primitive, _)| *primitive);
        if result
            .primitive_emitters
            .windows(2)
            .any(|pair| pair[0].0 == pair[1].0)
        {
            return Err(PhotometryError::InvalidOptions);
        }
        if flux_sum > 0.0 {
            let mut previous = 0.0;
            for emitter in &mut result.emitters {
                cx.checkpoint()?;
                emitter.cumulative_probability /= flux_sum;
                emitter.probability = emitter.cumulative_probability - previous;
                if emitter.source.flux_lm() > 0.0 && emitter.probability <= 0.0 {
                    return Err(PhotometryError::Source(SourceError::NumericalRange));
                }
                previous = emitter.cumulative_probability;
            }
            result.active = true;
        }
        cx.checkpoint()?;
        Ok(result)
    }

    pub(in crate::tracer) fn has_emission(&self) -> bool {
        self.active
    }

    /// Spectral emission at the actual hit, with its complementary NEE weight.
    #[allow(clippy::too_many_arguments)]
    pub(in crate::tracer) fn hit_radiance(
        &self,
        cx: &Cx<'_>,
        strategy: DirectStrategy,
        primitive_index: usize,
        origin: Point3,
        point: Point3,
        geometric_normal: Vec3,
        outgoing: Vec3,
        previous: Option<PreviousBsdf>,
        wavelengths: &[f64; PACKET],
    ) -> Result<[f64; PACKET], TracerError> {
        let Some(emitter_index) = self.emitter_for_primitive(primitive_index) else {
            return Ok([0.0; PACKET]);
        };
        let emitter = &self.emitters[emitter_index];
        if emitter.source.flux_lm() == 0.0
            || !emits_towards(emitter.source.sidedness(), geometric_normal, outgoing)
        {
            return Ok([0.0; PACKET]);
        }
        let luminance = emitter.luminance(point, outgoing)?;
        if luminance == 0.0 {
            return Ok([0.0; PACKET]);
        }
        let weight = if previous.is_none()
            || previous
                .is_some_and(|previous| previous.pdf.is_delta() || previous.smooth_slab.is_some())
            || strategy == DirectStrategy::BsdfOnly
        {
            1.0
        } else {
            let pdf = self.point_pdf(cx, emitter_index, origin, point)?;
            if pdf > 0.0 {
                emissive_hit_weight(strategy, previous, Some(pdf))
            } else {
                // No source proposal covers this numerical endpoint.
                1.0
            }
        };
        let mut radiance = [0.0; PACKET];
        for (lane, &wavelength) in wavelengths.iter().enumerate() {
            radiance[lane] =
                luminance * emitter.source.spectrum().value_w_per_lumen_nm(wavelength) * weight;
        }
        Ok(radiance)
    }

    /// One area proposal, paired with the original continuation technique.
    /// Probe results have spectral irradiance units; other receivers return
    /// radiance after their original BSDF or phase response.
    #[allow(clippy::too_many_arguments)]
    pub(in crate::tracer) fn direct_contribution(
        &self,
        scene: &Scene,
        fixtures: &FixtureLighting<'_>,
        cx: &Cx<'_>,
        settings: &Settings,
        depth: u32,
        medium: Option<PathMedium<'_, '_>>,
        point: Point3,
        receiver: Receiver<'_>,
        wavelengths: &[f64; PACKET],
        rng: &PathRng,
    ) -> Result<[f64; PACKET], TracerError> {
        let mut result = [0.0; PACKET];
        if !self.active
            || settings.strategy == DirectStrategy::BsdfOnly
            || matches!(
                receiver,
                Receiver::Surface {
                    material: Material::Dielectric { surface, .. },
                    ..
                } if surface.is_delta()
            )
        {
            // Smooth interfaces have no continuous source proposal. Their
            // native reflected/refracted continuations retain full hit weight,
            // without inspecting an unrelated sampled source connection.
            return Ok(result);
        }
        cx.checkpoint()?;
        let domain = if matches!(receiver, Receiver::Probe { .. }) {
            PROBE_AREA_NEE_DOMAIN
        } else {
            AREA_NEE_DOMAIN
        };
        let variates = philox4x32_10([rng.pixel, rng.sample, depth, domain], rng.key).map(u32_unit);
        let Some(sample) = self.sample(scene, cx, point, receiver, variates[0], variates[1])?
        else {
            return Ok(result);
        };
        let (direction, distance) = fixtures.admit_source_segment(
            point,
            sample
                .point
                .delta_from(point)
                .scale(1.0 / norm(sample.point.delta_from(point))?),
            Some(sample.point),
        )?;
        if !receiver.admits(direction) {
            return Ok(result);
        }
        let emitter = &self.emitters[sample.emitter];
        let luminance = emitter.luminance(sample.point, direction.scale(-1.0))?;
        if luminance == 0.0 {
            return Ok(result);
        }
        let opaque_pdf = receiver.opaque_pdf(direction)?;
        if !opaque_pdf.is_finite() || opaque_pdf < 0.0 {
            return Err(TracerError::InvalidInput);
        }
        if !matches!(receiver, Receiver::Phase { .. }) && distance <= fixtures.options.ray_epsilon_m
        {
            // This is a geometric resolution limit, not a point-light
            // singularity or a reason to clamp the area-source energy.
            return Err(TracerError::InvalidInput);
        }
        let shadow_origin = match receiver {
            Receiver::Phase { .. } => point,
            Receiver::Surface {
                frame,
                boundary: Some(_),
                ..
            } => dielectric_spawn_origin(point, frame.geometric, direction),
            Receiver::Surface {
                frame: SurfaceFrame {
                    oriented: normal, ..
                },
                ..
            }
            | Receiver::Probe { normal } => {
                let shifted = point.offset(direction.scale(fixtures.options.ray_epsilon_m));
                let displacement = shifted.delta_from(point);
                if [shifted.x, shifted.y, shifted.z]
                    .iter()
                    .any(|value| !value.is_finite())
                    || !displacement.dot(direction).is_finite()
                    || displacement.dot(direction) <= 0.0
                    || !displacement.dot(normal).is_finite()
                    || displacement.dot(normal) <= 0.0
                {
                    return Err(TracerError::InvalidInput);
                }
                shifted
            }
        };
        if !visible_endpoint(
            scene,
            &fixtures.geometry,
            cx,
            point,
            shadow_origin,
            sample.point,
            emitter.source.primitive_index(),
        )? {
            return Ok(result);
        }
        // The material stack owns the homogeneous glass interior. Ambient
        // haze is displaced there, including at a rough interface whose NEE
        // direction lies on the transmitted side.
        let shadow_medium = receiver.shadow_medium(direction);
        let ambient_medium = medium.filter(|_| shadow_medium.is_none());
        let transmittance = if let Some(glass) = shadow_medium {
            let mut values = [1.0; PACKET];
            for (lane, &wavelength) in wavelengths.iter().enumerate() {
                values[lane] = medium_transmittance(Some(glass), wavelength, distance)?;
            }
            values
        } else if let Some(medium) = ambient_medium {
            medium.transmittance(cx, point, direction, distance, wavelengths)?
        } else {
            [1.0; PACKET]
        };
        let survival =
            ambient_medium.map_or(1.0, |medium| medium.escape_probability(&transmittance));
        // A probe launches its incoming ray even when max_depth=1. Ordinary
        // terminal vertices do not evaluate another continuation source hit.
        let terminal =
            !matches!(receiver, Receiver::Probe { .. }) && depth + 1 == settings.max_depth;
        for (lane, &wavelength) in wavelengths.iter().enumerate() {
            // At a dispersive boundary the native continuation fans out by
            // wavelength. Its actual lane PDF, including event probability,
            // is the competitor rather than a hero-wavelength approximation.
            let (response, pdf) = receiver.response_and_pdf(direction, wavelength, opaque_pdf)?;
            let competing_pdf = pdf * survival;
            if !competing_pdf.is_finite() || competing_pdf < 0.0 {
                return Err(TracerError::InvalidInput);
            }
            let weight = if terminal {
                1.0
            } else {
                match settings.strategy {
                    DirectStrategy::Mis => {
                        balance_heuristic(1, sample.pdf_solid_angle, 1, competing_pdf)
                    }
                    DirectStrategy::PowerMis => {
                        power_heuristic(1, sample.pdf_solid_angle, 1, competing_pdf)
                    }
                    DirectStrategy::NeeOnly => 1.0,
                    DirectStrategy::BsdfOnly => return Ok(result),
                }
            };
            result[lane] = response
                * luminance
                * emitter.source.spectrum().value_w_per_lumen_nm(wavelength)
                * transmittance[lane]
                / sample.pdf_solid_angle
                * weight;
            if !result[lane].is_finite() || result[lane] < 0.0 {
                return Err(TracerError::InvalidInput);
            }
        }
        cx.checkpoint()?;
        Ok(result)
    }

    fn emitter_for_primitive(&self, primitive: usize) -> Option<usize> {
        self.primitive_emitters
            .binary_search_by_key(&primitive, |(index, _)| *index)
            .ok()
            .map(|index| self.primitive_emitters[index].1)
    }

    fn sample(
        &self,
        scene: &Scene,
        cx: &Cx<'_>,
        origin: Point3,
        receiver: Receiver<'_>,
        u1: f64,
        u2: f64,
    ) -> Result<Option<AreaSample>, TracerError> {
        let emitter_index = self
            .emitters
            .partition_point(|emitter| emitter.cumulative_probability <= u1);
        let emitter = self
            .emitters
            .get(emitter_index)
            .ok_or(TracerError::InvalidInput)?;
        let previous = emitter.cumulative_probability - emitter.probability;
        let local_u = ((u1 - previous) / emitter.probability).clamp(0.0, 1.0_f64.next_down());
        let faces = &self.faces[emitter.faces.clone()];
        let face_index = faces.partition_point(|face| face.cumulative_probability <= local_u);
        let face = faces.get(face_index).ok_or(TracerError::InvalidInput)?;
        if let Receiver::Surface {
            intersection,
            incoming_ray,
            ..
        } = receiver
            && intersection.primitive_index == emitter.source.primitive_index()
        {
            // A ray-derived hit can lie an ulp off its actual triangle plane.
            // Do not turn that error into a grazing self-lighting connection:
            // coplanar surface patches have exactly zero view factor. Retain
            // the native face witness instead of a cutoff on shadow angles or
            // distances; probes and other primitives never take this shortcut.
            let triangle_index = match intersection.instance_hit {
                Some(instance_hit) => match instance_hit.surface_feature {
                    InstanceSurfaceFeature::MeshTriangle { triangle_index, .. } => triangle_index,
                    InstanceSurfaceFeature::ChartUnavailable => {
                        return Err(TracerError::InvalidInput);
                    }
                },
                None => {
                    let Shape::Mesh(mesh) = &scene.primitives[intersection.primitive_index].shape
                    else {
                        return Err(TracerError::InvalidInput);
                    };
                    mesh.intersect_surface_with_cx_bounded(cx, incoming_ray, f64::MAX)?
                        .ok_or(TracerError::InvalidInput)?
                        .triangle_index
                }
            };
            let receiver_face = faces
                .get(triangle_index as usize)
                .ok_or(TracerError::InvalidInput)?;
            if face_index == triangle_index as usize || face.shares_plane(receiver_face) {
                return Ok(None);
            }
        }
        let previous = face.cumulative_probability - face.probability;
        let local_u = ((local_u - previous) / face.probability).clamp(0.0, 1.0_f64.next_down());
        let radius = local_u.sqrt();
        let point = face.vertices[0]
            .offset(
                face.vertices[1]
                    .delta_from(face.vertices[0])
                    .scale(radius * (1.0 - u2)),
            )
            .offset(
                face.vertices[2]
                    .delta_from(face.vertices[0])
                    .scale(radius * u2),
            );
        let outgoing = origin.delta_from(point);
        if !emits_towards(emitter.source.sidedness(), face.normal, outgoing) {
            return Ok(None);
        }
        let pdf_solid_angle = self.point_pdf(cx, emitter_index, origin, point)?;
        if !(pdf_solid_angle > 0.0 && pdf_solid_angle.is_finite()) {
            return Err(TracerError::InvalidInput);
        }
        Ok(Some(AreaSample {
            emitter: emitter_index,
            point,
            pdf_solid_angle,
        }))
    }

    /// Density for this actual surface endpoint. Only faces containing the
    /// point contribute; farther intersections on the same ray are excluded.
    fn point_pdf(
        &self,
        cx: &Cx<'_>,
        emitter_index: usize,
        origin: Point3,
        point: Point3,
    ) -> Result<f64, TracerError> {
        let emitter = &self.emitters[emitter_index];
        let displacement = origin.delta_from(point);
        let distance = norm(displacement)?;
        let outgoing = displacement.scale(1.0 / distance);
        let mut density = 0.0;
        for face in &self.faces[emitter.faces.clone()] {
            cx.checkpoint()?;
            if !face.contains(point)
                || !emits_towards(emitter.source.sidedness(), face.normal, outgoing)
            {
                continue;
            }
            let cosine = face.normal.dot(outgoing).abs();
            if cosine == 0.0 {
                continue;
            }
            density += emitter.probability
                * face.probability
                * (distance / face.area_m2)
                * (distance / cosine);
        }
        if !density.is_finite() || density < 0.0 {
            return Err(TracerError::InvalidInput);
        }
        Ok(density)
    }
}

impl Emitter<'_> {
    /// Source hits and NEE must evaluate the same pupil-position/direction
    /// law. The area proposal remains a density over geometric endpoints;
    /// a blocked gobo ray contributes zero without changing that density.
    fn luminance(&self, point: Point3, outgoing: Vec3) -> Result<f64, TracerError> {
        let Some(placement) = &self.projector else {
            return Ok(self.luminance_cd_m2);
        };
        let optics = self.source.projector().ok_or(TracerError::InvalidInput)?;
        let factor = optics
            .luminance_per_clear_flux_per_area(
                placement.world_to_optical.transform_point(point),
                placement.world_to_optical.transform_vector(outgoing),
            )
            .map_err(|_| TracerError::InvalidInput)?;
        let luminance = placement.clear_flux_per_area_lm_m2 * factor;
        if !luminance.is_finite() || luminance < 0.0 {
            return Err(TracerError::InvalidInput);
        }
        Ok(luminance)
    }
}

impl Receiver<'_> {
    fn admits(self, direction: Vec3) -> bool {
        match self {
            Self::Surface {
                material: Material::Dielectric { surface, .. },
                frame,
                ..
            } => !surface.is_delta() && frame.oriented.dot(direction) != 0.0,
            Self::Surface { frame, .. } => frame.oriented.dot(direction) > 0.0,
            Self::Probe { normal } => normal.dot(direction) > 0.0,
            Self::Phase { .. } => true,
        }
    }

    fn shadow_medium(self, direction: Vec3) -> Option<DielectricGlass> {
        match self {
            Self::Surface {
                frame,
                incident_medium,
                boundary,
                ..
            } => boundary.map_or(incident_medium, |boundary| {
                if frame.oriented.dot(direction) > 0.0 {
                    boundary.incident
                } else {
                    boundary.transmitted
                }
            }),
            Self::Phase { .. } | Self::Probe { .. } => None,
        }
    }

    /// Opaque, phase and probe densities are shared by the packet. Dielectric
    /// response evaluation supplies its wavelength-dependent density instead.
    fn opaque_pdf(self, direction: Vec3) -> Result<f64, TracerError> {
        match self {
            Self::Surface {
                material,
                frame,
                outgoing,
                brushed,
                ..
            } => Ok(if frame.oriented.dot(direction) > 0.0 {
                bsdf_pdf(material, frame.oriented, outgoing, direction, brushed)
            } else {
                0.0
            }),
            Self::Phase {
                phase,
                incoming_travel,
            } => Ok(phase.eval_cos(incoming_travel.dot(direction))?),
            Self::Probe { normal } => Ok(normal.dot(direction).max(0.0) / PI),
        }
    }

    fn response_and_pdf(
        self,
        direction: Vec3,
        wavelength: f64,
        opaque_pdf: f64,
    ) -> Result<(f64, f64), TracerError> {
        match self {
            Self::Surface {
                material,
                frame,
                outgoing,
                incident_medium,
                boundary,
                brushed,
                ..
            } => {
                if let Material::Dielectric { surface, .. } = material {
                    let Some(alpha) = surface.roughness_alpha() else {
                        return Ok((0.0, 0.0));
                    };
                    let boundary = boundary.ok_or(TracerError::InvalidInput)?;
                    let evaluation = evaluate_rough_dielectric(
                        frame.oriented,
                        outgoing,
                        direction,
                        medium_ior(boundary.incident, wavelength)?,
                        medium_ior(boundary.transmitted, wavelength)?,
                        alpha,
                    )?;
                    Ok((
                        frame.oriented.dot(direction).abs() * evaluation.value,
                        evaluation.pdf,
                    ))
                } else {
                    let cosine = frame.oriented.dot(direction).max(0.0);
                    Ok((
                        cosine
                            * opaque_bsdf_eval(
                                material,
                                frame.oriented,
                                outgoing,
                                direction,
                                wavelength,
                                incident_medium,
                                brushed,
                            )?,
                        opaque_pdf,
                    ))
                }
            }
            Self::Phase {
                phase,
                incoming_travel,
            } => Ok((phase.eval_cos(incoming_travel.dot(direction))?, opaque_pdf)),
            Self::Probe { normal } => Ok((normal.dot(direction).max(0.0), opaque_pdf)),
        }
    }
}

impl Face {
    /// Compare native supporting planes, not a rounded path hit. The bound
    /// covers cancellation in the scalar triple product, with no metric or
    /// angular epsilon that could erase a separate nearly grazing surface.
    fn shares_plane(&self, other: &Self) -> bool {
        let on_plane = |vertices: [Point3; 3], point: Point3| {
            let vectors = [
                vertices[1].delta_from(vertices[0]),
                vertices[2].delta_from(vertices[0]),
                point.delta_from(vertices[0]),
            ];
            let scale = vectors
                .iter()
                .flat_map(|vector| [vector.x.abs(), vector.y.abs(), vector.z.abs()])
                .fold(0.0, f64::max);
            if scale == 0.0 {
                return true;
            }
            let [a, b, p] = vectors.map(|vector| vector.scale(1.0 / scale));
            let products = [
                [a.y * b.z, a.z * b.y, p.x],
                [a.z * b.x, a.x * b.z, p.y],
                [a.x * b.y, a.y * b.x, p.z],
            ];
            let determinant = products
                .iter()
                .map(|[positive, negative, coordinate]| (positive - negative) * coordinate)
                .sum::<f64>();
            let permanent = products
                .iter()
                .map(|[positive, negative, coordinate]| {
                    (positive.abs() + negative.abs()) * coordinate.abs()
                })
                .sum::<f64>();
            determinant.is_finite()
                && permanent.is_finite()
                && determinant.abs() <= 16.0 * f64::EPSILON * permanent
        };
        other
            .vertices
            .iter()
            .all(|&point| on_plane(self.vertices, point))
            && self
                .vertices
                .iter()
                .all(|&point| on_plane(other.vertices, point))
    }

    // Same scaled planar/barycentric admission used by native thermal area
    // proposals, with winding retained for one-sided calibrated emission.
    fn contains(&self, point: Point3) -> bool {
        let relative = point.delta_from(self.vertices[0]);
        let scale = self.vertices.iter().fold(
            point.x.abs().max(point.y.abs()).max(point.z.abs()),
            |scale, vertex| {
                scale
                    .max(vertex.x.abs())
                    .max(vertex.y.abs())
                    .max(vertex.z.abs())
            },
        );
        let plane_error = relative.dot(self.normal).abs();
        if !plane_error.is_finite() || plane_error > 128.0 * f64::EPSILON * scale {
            return false;
        }
        let project = |vector: Vec3| {
            if self.normal.x.abs() >= self.normal.y.abs()
                && self.normal.x.abs() >= self.normal.z.abs()
            {
                [vector.y, vector.z]
            } else if self.normal.y.abs() >= self.normal.z.abs() {
                [vector.x, vector.z]
            } else {
                [vector.x, vector.y]
            }
        };
        let a = project(self.vertices[1].delta_from(self.vertices[0]));
        let b = project(self.vertices[2].delta_from(self.vertices[0]));
        let p = project(relative);
        let scale = a[0].abs().max(a[1].abs()).max(b[0].abs()).max(b[1].abs());
        let a = a.map(|value| value / scale);
        let b = b.map(|value| value / scale);
        let p = p.map(|value| value / scale);
        let determinant = a[0] * b[1] - a[1] * b[0];
        let u = (p[0] * b[1] - p[1] * b[0]) / determinant;
        let v = (a[0] * p[1] - a[1] * p[0]) / determinant;
        let tolerance = 128.0 * f64::EPSILON;
        u.is_finite()
            && v.is_finite()
            && u >= -tolerance
            && v >= -tolerance
            && u + v <= 1.0 + tolerance
    }
}

fn emits_towards(sidedness: EmissionSidedness, normal: Vec3, outgoing: Vec3) -> bool {
    let cosine = normal.dot(outgoing);
    match sidedness {
        EmissionSidedness::Front => cosine > 0.0,
        EmissionSidedness::Both => cosine != 0.0,
    }
}

fn triangle_key(vertices: [Point3; 3]) -> [[u64; 3]; 3] {
    let mut key = vertices.map(|point| {
        [point.x, point.y, point.z].map(|value| if value == 0.0 { 0 } else { value.to_bits() })
    });
    key.sort_unstable();
    key
}

fn visible_endpoint(
    scene: &Scene,
    geometry: &geometry::PreparedGeometry,
    cx: &Cx<'_>,
    source: Point3,
    shadow_origin: Point3,
    target: Point3,
    primitive_index: usize,
) -> Result<bool, TracerError> {
    let displacement = target.delta_from(shadow_origin);
    let distance = norm(displacement)?;
    let ray = Ray {
        origin: shadow_origin,
        dir: displacement.scale(1.0 / distance),
    };
    let Some(hit) = geometry.intersect_scene(&scene.primitives, cx, &ray, f64::MAX)? else {
        return Ok(false);
    };
    let scale = [source.x, source.y, source.z, target.x, target.y, target.z]
        .into_iter()
        .map(f64::abs)
        .fold(norm(target.delta_from(source))?, f64::max);
    // A metric epsilon here would leak nearby hidden faces of the same body.
    let tolerance = 128.0 * f64::EPSILON * scale;
    Ok(hit.primitive_index == primitive_index
        && hit.hit.point.delta_from(target).norm() <= tolerance)
}
