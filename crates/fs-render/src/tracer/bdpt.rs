//! Bidirectional path-space transport for the spectral tracer.
//!
//! This module follows Veach's path-space construction: camera and finite-
//! emitter subpaths retain reciprocal area-measure densities, every admissible
//! `(s,t)` vertex-connection strategy is evaluated, and the strategies are
//! combined with balance-heuristic multiple importance sampling. Infinite
//! environment illumination remains a disjoint camera-subpath/NEE estimator;
//! no finite launch surface is invented for it.

use super::transport::{TransportMode, refractive_transport_factor};
use super::*;

/// Bit-affecting semantics of the opt-in bidirectional path integrator.
/// Version 2 corrects the global sample multiplicity and estimator
/// normalization of light-subpath strategies that splat into the raster.
pub const BIDIRECTIONAL_TRACER_SEMANTICS_VERSION: u32 = 2;

const CAMERA_WALK_DOMAIN: u32 = 0x4244_4341;
const LIGHT_WALK_DOMAIN: u32 = 0x4244_4c49;
const CONNECTION_DOMAIN: u32 = 0x4244_434e;

/// Counts of evaluated and nonzero path-space strategies.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct BidirectionalStrategyStats {
    /// Candidate `(s,t)` strategies within the configured depth bound.
    pub evaluated: u64,
    /// Candidates carrying a finite, positive contribution before MIS.
    pub nonzero: u64,
    /// Light-subpath contributions projected into a raster pixel (`t=1`).
    pub camera_splats: u64,
}

/// A complete fixed-SPP bidirectional render and direct strategy evidence.
#[derive(Debug, Clone, PartialEq)]
pub struct BidirectionalRenderOutput {
    /// Raw spectral CIE XYZ film. Values are summed over `spp_done`.
    pub film: Film,
    /// Strategy counters for diagnosing unsupported or zero-contribution paths.
    pub strategies: BidirectionalStrategyStats,
}

/// Execution evidence for sample-block-parallel bidirectional rendering.
///
/// Ordinary camera-path rendering owns disjoint image tiles. BDPT cannot use
/// that ownership rule because `t=1` light subpaths splat into arbitrary
/// pixels. Instead, each logical task owns a fixed sample-index block over the
/// complete raster and the completed blocks are reduced in ascending sample
/// order. The block shape is semantics, not a worker-count decision.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BidirectionalExecutionReport {
    /// Raster geometry associated with the film. This is not the executor's
    /// sample-block plan.
    pub image_layout: RenderTileLayout,
    /// Number of fixed logical sample blocks dispatched.
    pub sample_blocks: u64,
    /// Maximum samples in a logical block.
    pub samples_per_block: u32,
    /// Caller-requested worker ceiling.
    pub requested_workers: usize,
    /// Workers actually admitted by the executor.
    pub workers: usize,
    /// Bytes retained by all private sample-block films.
    pub staging_film_bytes: u64,
    /// Conservative concurrent per-source splat scratch admission.
    pub splat_scratch_envelope_bytes: u64,
    /// Shared Sobol direction-state payload charged for the run.
    pub sampler_state_bytes: u64,
    /// Setup and admission wall time.
    pub setup_ns: u64,
    /// Throughput-lane wall time, including drain on failure.
    pub traversal_ns: u64,
    /// Sum of sample-block compute time across workers.
    pub block_compute_ns: u64,
    /// Time spent publishing completed blocks in canonical order.
    pub publication_ns: u64,
    /// Conservative worker-capacity time not accounted for by block compute.
    pub idle_worker_ns: u64,
    /// Executor scheduling, completion, and cancellation diagnostics.
    pub executor: RunReport,
    /// Exact operation-memory admission trace after transient charges release.
    pub memory: LeaseReceipt,
}

/// A complete sample-block-parallel BDPT render and its execution evidence.
#[derive(Debug, Clone, PartialEq)]
pub struct BidirectionalRenderExecutionOutput {
    /// Fully published raw spectral CIE XYZ film.
    pub film: Film,
    /// Direct path-strategy counters.
    pub strategies: BidirectionalStrategyStats,
    /// Deterministic execution and memory evidence.
    pub report: BidirectionalExecutionReport,
}

const BIDIRECTIONAL_SAMPLE_BLOCK: u32 = 8;
const BIDIRECTIONAL_SAMPLE_KERNEL: &str = "fs-render/bidirectional-sample-block-v1";

#[derive(Clone, Copy)]
enum VertexKind {
    Camera,
    Light {
        light_index: usize,
        primitive_index: usize,
        emission: (LiftedSpectrum, f64),
    },
    Surface {
        primitive_index: usize,
    },
}

#[derive(Clone, Copy)]
struct Vertex {
    kind: VertexKind,
    point: Point3,
    geometric_normal: Option<Vec3>,
    /// Direction from this vertex to its predecessor.
    wo: Option<Vec3>,
    beta: f64,
    pdf_fwd: f64,
    pdf_rev: f64,
    delta: bool,
    positive_medium: Option<MediumEntry>,
    negative_medium: Option<MediumEntry>,
}

impl Vertex {
    fn camera(point: Point3) -> Self {
        Self {
            kind: VertexKind::Camera,
            point,
            geometric_normal: None,
            wo: None,
            beta: 1.0,
            pdf_fwd: 1.0,
            pdf_rev: 0.0,
            delta: false,
            positive_medium: None,
            negative_medium: None,
        }
    }

    fn light(sample: crate::lighting::RectEmissionSample, wavelength_nm: f64) -> Self {
        let (spectrum, scale) = sample.emission;
        Self {
            kind: VertexKind::Light {
                light_index: sample.light_index,
                primitive_index: sample.primitive_index,
                emission: sample.emission,
            },
            point: sample.point,
            geometric_normal: Some(sample.normal),
            wo: None,
            beta: spectrum.eval(wavelength_nm) * scale,
            pdf_fwd: sample.pdf_position_area,
            pdf_rev: 0.0,
            delta: false,
            positive_medium: None,
            negative_medium: None,
        }
    }

    fn surface(
        scene: &Scene,
        intersection: SceneIntersection,
        ray: Ray,
        frame: SurfaceFrame,
        stack: &MediumStack,
        beta: f64,
        pdf_fwd: f64,
    ) -> Self {
        let primitive_index = intersection.primitive_index;
        let primitive = &scene.primitives[primitive_index];
        let (positive_medium, negative_medium) = match primitive.material {
            Material::Dielectric { glass, .. } => {
                let positive = if frame.entering {
                    stack.last().copied()
                } else {
                    stack
                        .len()
                        .checked_sub(2)
                        .and_then(|index| stack.get(index))
                        .copied()
                };
                (
                    positive,
                    Some(MediumEntry {
                        boundary_primitive: primitive_index,
                        glass,
                    }),
                )
            }
            Material::Lambertian { .. } | Material::Ggx { .. } | Material::Conductor { .. } => {
                let medium = stack.last().copied();
                (medium, medium)
            }
        };
        Self {
            kind: VertexKind::Surface { primitive_index },
            point: intersection.hit.point,
            geometric_normal: Some(frame.geometric),
            wo: Some(ray.dir.scale(-1.0)),
            beta,
            pdf_fwd,
            pdf_rev: 0.0,
            delta: false,
            positive_medium,
            negative_medium,
        }
    }

    fn is_surface(self) -> bool {
        matches!(self.kind, VertexKind::Surface { .. })
    }

    fn is_connectible(self, scene: &Scene) -> bool {
        match self.kind {
            VertexKind::Camera | VertexKind::Light { .. } => true,
            VertexKind::Surface { primitive_index } => {
                match scene.primitives[primitive_index].material {
                    Material::Dielectric { surface, .. } => !surface.is_delta(),
                    Material::Lambertian { .. }
                    | Material::Ggx { .. }
                    | Material::Conductor { .. } => true,
                }
            }
        }
    }

    fn medium_toward(self, direction: Vec3) -> Option<MediumEntry> {
        let Some(normal) = self.geometric_normal else {
            return None;
        };
        if normal.dot(direction) >= 0.0 {
            self.positive_medium
        } else {
            self.negative_medium
        }
    }

    fn emitted_radiance(self, scene: &Scene, wavelength_nm: f64) -> f64 {
        let emission = match self.kind {
            VertexKind::Light { emission, .. } => Some(emission),
            VertexKind::Surface { primitive_index } => scene.primitives[primitive_index].emission,
            VertexKind::Camera => None,
        };
        emission.map_or(0.0, |(spectrum, scale)| {
            spectrum.eval(wavelength_nm) * scale
        })
    }
}

struct CameraEscape {
    beta: f64,
    origin: Point3,
    direction: Vec3,
    previous_pdf_solid_angle: f64,
    previous_delta: bool,
}

struct CameraSubpath {
    vertices: Vec<Vertex>,
    escape: Option<CameraEscape>,
}

/// One cross-pixel camera splat staged under logical, scheduling-independent
/// identity. A source pixel publishes these records only after all of its
/// samples finish; source pixels themselves publish in ascending order.
#[derive(Clone, Copy)]
struct SplatRecord {
    target_pixel: u32,
    sample: u32,
    strategy_s: usize,
    xyz: [f64; 3],
}

fn publish_source_splats(film: &mut Film, splats: &mut Vec<SplatRecord>) {
    // This is the renderer analogue of fs-sparse COO canonical assembly:
    // arrival order is irrelevant because logical sample identity owns the
    // reduction order. The current serial source-pixel loop bounds scratch;
    // a parallel tile implementation can retain the same per-source batches
    // and merge source ranges in ascending order.
    splats.sort_by_key(|record| (record.target_pixel, record.sample, record.strategy_s));
    for splat in splats.drain(..) {
        add_scaled_xyz(&mut film.xyz[splat.target_pixel as usize], splat.xyz);
    }
}

#[derive(Clone, Copy)]
struct SubpathRng<'a> {
    pixel: u32,
    sample: u32,
    dimension: u32,
    domain: u32,
    key: [u32; 2],
    sobol: Option<&'a Sobol>,
}

impl<'a> SubpathRng<'a> {
    fn new(
        settings: &Settings,
        pixel: u32,
        sample: u32,
        domain: u32,
        sobol: Option<&'a Sobol>,
    ) -> Self {
        Self {
            pixel,
            sample,
            dimension: 0,
            domain,
            key: [
                (settings.seed & 0xffff_ffff) as u32,
                (settings.seed >> 32) as u32,
            ],
            sobol,
        }
    }

    fn next4(&mut self) -> [f64; 4] {
        if let Some(sobol) = self.sobol {
            // Each logical path block owns a separate nested Owen tree. This
            // avoids data-dependent dimension consumption while retaining the
            // sample index as the Sobol point index. Camera walk, light walk,
            // and every `(s,t)` connector use disjoint domains.
            let seed_words = philox4x32_10(
                [self.pixel, self.domain, self.dimension, 0x4244_514d],
                self.key,
            );
            let scramble_seed = u64::from(seed_words[0]) | (u64::from(seed_words[1]) << 32);
            let mut values = [0.0; PATH_QMC_DIMENSIONS];
            sobol.point_with_owen_scramble(self.sample, &mut values, scramble_seed);
            self.dimension = self.dimension.wrapping_add(1);
            return [values[0], values[1], values[2], values[3]];
        }
        let words = philox4x32_10(
            [self.pixel, self.sample, self.dimension, self.domain],
            self.key,
        );
        self.dimension = self.dimension.wrapping_add(1);
        words.map(u32_unit)
    }
}

#[derive(Clone, Copy)]
struct ScatterSample {
    direction: Vec3,
    weight: f64,
    pdf_fwd_solid_angle: f64,
    pdf_rev_solid_angle: f64,
    delta: bool,
    event: Option<DielectricEvent>,
}

fn same_medium(left: Option<MediumEntry>, right: Option<MediumEntry>) -> bool {
    match (left, right) {
        (None, None) => true,
        (Some(left), Some(right)) => {
            left.boundary_primitive == right.boundary_primitive && left.glass == right.glass
        }
        (None, Some(_)) | (Some(_), None) => false,
    }
}

fn vector_between(from: Point3, to: Point3) -> Option<(Vec3, f64)> {
    let displacement = to.delta_from(from);
    let distance_squared = displacement.dot(displacement);
    if !distance_squared.is_finite() || distance_squared <= 0.0 {
        return None;
    }
    let distance = distance_squared.sqrt();
    Some((displacement.scale(1.0 / distance), distance_squared))
}

fn convert_density(source: &Vertex, pdf_solid_angle: f64, next: &Vertex) -> f64 {
    if !pdf_solid_angle.is_finite() || pdf_solid_angle <= 0.0 {
        return 0.0;
    }
    let Some((direction, distance_squared)) = vector_between(source.point, next.point) else {
        return 0.0;
    };
    let cosine = next
        .geometric_normal
        .map_or(1.0, |normal| normal.dot(direction.scale(-1.0)).abs());
    pdf_solid_angle * cosine / distance_squared
}

fn oriented_normal(vertex: &Vertex, outgoing: Vec3) -> Option<Vec3> {
    let normal = vertex.geometric_normal?;
    Some(if normal.dot(outgoing) >= 0.0 {
        normal
    } else {
        normal.scale(-1.0)
    })
}

fn surface_scattering(
    scene: &Scene,
    vertex: &Vertex,
    predecessor_direction: Vec3,
    next_direction: Vec3,
    wavelength_nm: f64,
    mode: TransportMode,
) -> Result<f64, TracerError> {
    let VertexKind::Surface { primitive_index } = vertex.kind else {
        return Ok(0.0);
    };
    let material = &scene.primitives[primitive_index].material;
    let normal =
        oriented_normal(vertex, predecessor_direction).ok_or(TracerError::MissingNormal)?;
    match material {
        Material::Lambertian { .. } | Material::Ggx { .. } | Material::Conductor { .. } => {
            let medium = vertex
                .medium_toward(predecessor_direction)
                .map(|entry| entry.glass);
            opaque_bsdf_eval(
                material,
                normal,
                predecessor_direction,
                next_direction,
                wavelength_nm,
                medium,
            )
        }
        Material::Dielectric { surface, .. } => {
            let Some(alpha) = surface.roughness_alpha() else {
                return Ok(0.0);
            };
            let incident = vertex
                .medium_toward(predecessor_direction)
                .map(|entry| entry.glass);
            let transmitted = vertex
                .medium_toward(predecessor_direction.scale(-1.0))
                .map(|entry| entry.glass);
            let eta_i = medium_ior(incident, wavelength_nm)?;
            let eta_t = medium_ior(transmitted, wavelength_nm)?;
            let evaluation = evaluate_rough_dielectric(
                normal,
                predecessor_direction,
                next_direction,
                eta_i,
                eta_t,
                alpha,
            )?;
            let mut value = evaluation.value;
            if mode == TransportMode::Importance
                && evaluation.event == DielectricEvent::Transmission
            {
                let radiance_factor =
                    refractive_transport_factor(TransportMode::Radiance, eta_i, eta_t)
                        .ok_or(TracerError::InvalidInput)?;
                value /= radiance_factor;
            }
            Ok(value)
        }
    }
}

fn surface_pdf_solid_angle(
    scene: &Scene,
    vertex: &Vertex,
    predecessor_direction: Vec3,
    next_direction: Vec3,
    wavelength_nm: f64,
) -> Result<f64, TracerError> {
    let VertexKind::Surface { primitive_index } = vertex.kind else {
        return Ok(0.0);
    };
    let material = &scene.primitives[primitive_index].material;
    let normal =
        oriented_normal(vertex, predecessor_direction).ok_or(TracerError::MissingNormal)?;
    match material {
        Material::Lambertian { .. } | Material::Ggx { .. } | Material::Conductor { .. } => Ok(
            bsdf_pdf(material, normal, predecessor_direction, next_direction),
        ),
        Material::Dielectric { surface, .. } => {
            let Some(alpha) = surface.roughness_alpha() else {
                return Ok(0.0);
            };
            let incident = vertex
                .medium_toward(predecessor_direction)
                .map(|entry| entry.glass);
            let transmitted = vertex
                .medium_toward(predecessor_direction.scale(-1.0))
                .map(|entry| entry.glass);
            let eta_i = medium_ior(incident, wavelength_nm)?;
            let eta_t = medium_ior(transmitted, wavelength_nm)?;
            Ok(evaluate_rough_dielectric(
                normal,
                predecessor_direction,
                next_direction,
                eta_i,
                eta_t,
                alpha,
            )?
            .pdf)
        }
    }
}

fn sample_surface(
    scene: &Scene,
    vertex: &Vertex,
    wavelength_nm: f64,
    mode: TransportMode,
    random: [f64; 4],
) -> Result<Option<ScatterSample>, TracerError> {
    let VertexKind::Surface { primitive_index } = vertex.kind else {
        return Ok(None);
    };
    let wo = vertex.wo.ok_or(TracerError::InvalidInput)?;
    let material = &scene.primitives[primitive_index].material;
    let normal = oriented_normal(vertex, wo).ok_or(TracerError::MissingNormal)?;
    match material {
        Material::Lambertian { .. } | Material::Ggx { .. } | Material::Conductor { .. } => {
            let Some((wi, pdf_fwd)) = bsdf_sample(material, normal, wo, random[0], random[1])
            else {
                return Ok(None);
            };
            let cosine = normal.dot(wi).max(0.0);
            let incident_medium = vertex.medium_toward(wo).map(|entry| entry.glass);
            let value = opaque_bsdf_eval(material, normal, wo, wi, wavelength_nm, incident_medium)?;
            let pdf_rev = bsdf_pdf(material, normal, wi, wo);
            let weight = value * cosine / pdf_fwd;
            Ok(
                (weight.is_finite() && weight > 0.0).then_some(ScatterSample {
                    direction: wi,
                    weight,
                    pdf_fwd_solid_angle: pdf_fwd,
                    pdf_rev_solid_angle: pdf_rev,
                    delta: false,
                    event: None,
                }),
            )
        }
        Material::Dielectric { surface, .. } => {
            let incident = vertex.medium_toward(wo).map(|entry| entry.glass);
            let transmitted = vertex
                .medium_toward(wo.scale(-1.0))
                .map(|entry| entry.glass);
            let eta_i = medium_ior(incident, wavelength_nm)?;
            let eta_t = medium_ior(transmitted, wavelength_nm)?;
            if let Some(alpha) = surface.roughness_alpha() {
                let Some(sample) = sample_rough_dielectric(
                    normal, wo, eta_i, eta_t, alpha, random[0], random[1], random[2],
                )?
                else {
                    return Ok(None);
                };
                let mut weight = sample.radiance_weight;
                if mode == TransportMode::Importance
                    && sample.event == DielectricEvent::Transmission
                {
                    let radiance_factor =
                        refractive_transport_factor(TransportMode::Radiance, eta_i, eta_t)
                            .ok_or(TracerError::InvalidInput)?;
                    weight /= radiance_factor;
                }
                let pdf_rev = if sample.delta {
                    0.0
                } else {
                    surface_pdf_solid_angle(scene, vertex, sample.direction, wo, wavelength_nm)?
                };
                Ok(Some(ScatterSample {
                    direction: sample.direction,
                    weight,
                    pdf_fwd_solid_angle: if sample.delta { 0.0 } else { sample.pdf },
                    pdf_rev_solid_angle: pdf_rev,
                    delta: sample.delta,
                    event: Some(sample.event),
                }))
            } else {
                let sample = sample_smooth_dielectric(normal, wo, eta_i, eta_t, random[2])?;
                let weight = match mode {
                    TransportMode::Radiance => sample.radiance_weight,
                    TransportMode::Importance => 1.0,
                };
                Ok(Some(ScatterSample {
                    direction: sample.direction,
                    weight,
                    pdf_fwd_solid_angle: 0.0,
                    pdf_rev_solid_angle: 0.0,
                    delta: true,
                    event: Some(sample.event),
                }))
            }
        }
    }
}

fn update_stack_after_sample(
    scene: &Scene,
    vertex: &Vertex,
    stack: &mut MediumStack,
    sample: ScatterSample,
) -> Result<(), TracerError> {
    if sample.event != Some(DielectricEvent::Transmission) {
        return Ok(());
    }
    let VertexKind::Surface { primitive_index } = vertex.kind else {
        return Ok(());
    };
    let Material::Dielectric { glass, .. } = scene.primitives[primitive_index].material else {
        return Ok(());
    };
    let geometric = vertex.geometric_normal.ok_or(TracerError::MissingNormal)?;
    let wo = vertex.wo.ok_or(TracerError::InvalidInput)?;
    let entering = geometric.dot(wo) > 0.0;
    let boundary = boundary_media(primitive_index, glass, entering, stack)?;
    apply_medium_transition(stack, boundary.transition)
}

#[allow(clippy::too_many_arguments)]
fn random_walk(
    scene: &Scene,
    cx: &Cx<'_>,
    ray_time: Option<&PathTime>,
    wavelength_nm: f64,
    mode: TransportMode,
    mut ray: Ray,
    mut beta: f64,
    mut pdf_fwd_solid_angle: f64,
    max_vertices: usize,
    vertices: &mut Vec<Vertex>,
    rng: &mut SubpathRng<'_>,
) -> Result<Option<CameraEscape>, TracerError> {
    let mut stack = MediumStack::new();
    let mut segment_origin = ray.origin;
    while vertices.len() < max_vertices {
        cx.checkpoint()?;
        let Some(intersection) = intersect(scene, cx, &ray, ray_time)? else {
            if stack.last().is_some() {
                return Err(TracerError::InvalidInput);
            }
            return Ok((mode == TransportMode::Radiance).then_some(CameraEscape {
                beta,
                origin: segment_origin,
                direction: ray.dir,
                previous_pdf_solid_angle: pdf_fwd_solid_angle,
                // A primary camera ray has no competing light-sampling
                // strategy at a surface. Treat the sensor endpoint like a
                // delta predecessor for the environment MIS decision.
                previous_delta: vertices.len() == 1
                    || vertices.last().is_some_and(|vertex| vertex.delta),
            }));
        };

        if let Some(active) = stack.last() {
            let distance_m = intersection.hit.point.delta_from(segment_origin).norm();
            beta *= active
                .glass
                .absorption()
                .transmittance(wavelength_nm, distance_m)?;
        }
        if !beta.is_finite() || beta <= 0.0 {
            break;
        }
        let frame = surface_frame(&intersection.hit, &ray)?;
        let placeholder = Vertex {
            kind: VertexKind::Surface {
                primitive_index: intersection.primitive_index,
            },
            point: intersection.hit.point,
            geometric_normal: Some(frame.geometric),
            wo: Some(ray.dir.scale(-1.0)),
            beta,
            pdf_fwd: 0.0,
            pdf_rev: 0.0,
            delta: false,
            positive_medium: None,
            negative_medium: None,
        };
        let pdf_area = convert_density(
            vertices.last().ok_or(TracerError::InvalidInput)?,
            pdf_fwd_solid_angle,
            &placeholder,
        );
        let vertex = Vertex::surface(scene, intersection, ray, frame, &stack, beta, pdf_area);
        vertices.push(vertex);

        if scene.primitives[intersection.primitive_index]
            .emission
            .is_some()
        {
            break;
        }
        if vertices.len() >= max_vertices {
            break;
        }

        let Some(sample) = sample_surface(
            scene,
            vertices.last().unwrap(),
            wavelength_nm,
            mode,
            rng.next4(),
        )?
        else {
            break;
        };
        let current_index = vertices.len() - 1;
        vertices[current_index].delta = sample.delta;
        if current_index > 0 {
            let reverse_area = convert_density(
                &vertices[current_index],
                sample.pdf_rev_solid_angle,
                &vertices[current_index - 1],
            );
            vertices[current_index - 1].pdf_rev = reverse_area;
        }
        update_stack_after_sample(scene, &vertices[current_index], &mut stack, sample)?;
        beta *= sample.weight;
        pdf_fwd_solid_angle = sample.pdf_fwd_solid_angle;
        segment_origin = vertices[current_index].point;
        ray = Ray {
            origin: dielectric_spawn_origin(
                vertices[current_index].point,
                vertices[current_index]
                    .geometric_normal
                    .ok_or(TracerError::MissingNormal)?,
                sample.direction,
            ),
            dir: sample.direction,
        };
    }
    Ok(None)
}

#[allow(clippy::too_many_arguments)]
fn generate_camera_subpath(
    scene: &Scene,
    cx: &Cx<'_>,
    settings: &Settings,
    physical_camera: &PhysicalCamera,
    ray: Ray,
    ray_time: Option<&PathTime>,
    wavelength_nm: f64,
    pixel: u32,
    sample: u32,
    path_sobol: Option<&Sobol>,
) -> Result<CameraSubpath, TracerError> {
    let raster = physical_camera
        .pinhole_raster_sample(ray.origin.offset(ray.dir), settings.width, settings.height)?
        .ok_or(TracerError::InvalidInput)?;
    if raster.pixel != pixel {
        return Err(TracerError::InvalidInput);
    }
    let mut vertices = Vec::with_capacity(settings.max_depth as usize + 2);
    vertices.push(Vertex::camera(ray.origin));
    let mut rng = SubpathRng::new(settings, pixel, sample, CAMERA_WALK_DOMAIN, path_sobol);
    let escape = random_walk(
        scene,
        cx,
        ray_time,
        wavelength_nm,
        TransportMode::Radiance,
        ray,
        1.0,
        raster.pdf_solid_angle,
        settings.max_depth as usize + 2,
        &mut vertices,
        &mut rng,
    )?;
    Ok(CameraSubpath { vertices, escape })
}

fn generate_light_subpath(
    scene: &Scene,
    lighting: &AdmittedLighting<'_>,
    cx: &Cx<'_>,
    settings: &Settings,
    ray_time: Option<&PathTime>,
    wavelength_nm: f64,
    pixel: u32,
    sample: u32,
    path_sobol: Option<&Sobol>,
) -> Result<Vec<Vertex>, TracerError> {
    let mut rng = SubpathRng::new(settings, pixel, sample, LIGHT_WALK_DOMAIN, path_sobol);
    let random = rng.next4();
    let Some(emission) =
        lighting.sample_rectangle_emission(random[0], random[1], random[2], random[3])
    else {
        return Ok(Vec::new());
    };
    let mut vertices = Vec::with_capacity(settings.max_depth as usize + 1);
    vertices.push(Vertex::light(emission, wavelength_nm));
    let (spectrum, scale) = emission.emission;
    let emitted = spectrum.eval(wavelength_nm) * scale;
    let cosine = emission.normal.dot(emission.direction).abs();
    let beta = emitted * cosine / (emission.pdf_position_area * emission.pdf_direction_solid_angle);
    let ray = Ray {
        origin: dielectric_spawn_origin(emission.point, emission.normal, emission.direction),
        dir: emission.direction,
    };
    let _ = random_walk(
        scene,
        cx,
        ray_time,
        wavelength_nm,
        TransportMode::Importance,
        ray,
        beta,
        emission.pdf_direction_solid_angle,
        settings.max_depth as usize + 1,
        &mut vertices,
        &mut rng,
    )?;
    Ok(vertices)
}

fn path_time_for_sample(
    scene: &Scene,
    cx: &Cx<'_>,
    settings: &Settings,
    shutter: ShutterInterval,
    pixel: u32,
    sample: u32,
) -> Result<PathTime, TracerError> {
    let normalized = shutter.sample_for_stream(settings.seed, u64::from(pixel), u64::from(sample));
    let absolute_time_s = shutter.time_at(normalized);
    let mut cached_animated = std::array::from_fn(|_| None);
    let mut cached_count = 0;
    for (primitive_index, primitive) in scene.primitives.iter().enumerate() {
        if let Shape::AnimatedInstance(instance) = &primitive.shape {
            let Some(slot) = cached_animated.get_mut(cached_count) else {
                break;
            };
            *slot = Some(CachedAnimatedInstance {
                primitive_index,
                instance: instance.instance_at(cx, absolute_time_s)?,
            });
            cached_count += 1;
        }
    }
    Ok(PathTime {
        interval: shutter,
        normalized,
        cached_animated,
    })
}

fn vertex_pdf(
    scene: &Scene,
    lighting: &AdmittedLighting<'_>,
    camera: &PhysicalCamera,
    current: &Vertex,
    previous: Option<&Vertex>,
    next: &Vertex,
    wavelength_nm: f64,
    width: u32,
    height: u32,
) -> Result<f64, TracerError> {
    let Some((direction, _)) = vector_between(current.point, next.point) else {
        return Ok(0.0);
    };
    let pdf_solid_angle = match current.kind {
        VertexKind::Camera => camera
            .pinhole_raster_sample(next.point, width, height)?
            .map_or(0.0, |sample| sample.pdf_solid_angle),
        VertexKind::Light { light_index, .. } => lighting
            .rectangle_emission_pdfs(light_index, direction)
            .map_or(0.0, |(_, direction_pdf)| direction_pdf),
        VertexKind::Surface { .. } => {
            let Some(previous) = previous else {
                return Ok(0.0);
            };
            let Some((predecessor_direction, _)) = vector_between(current.point, previous.point)
            else {
                return Ok(0.0);
            };
            surface_pdf_solid_angle(
                scene,
                current,
                predecessor_direction,
                direction,
                wavelength_nm,
            )?
        }
    };
    Ok(convert_density(current, pdf_solid_angle, next))
}

fn light_identity(lighting: &AdmittedLighting<'_>, vertex: &Vertex) -> Option<usize> {
    match vertex.kind {
        VertexKind::Light { light_index, .. } => Some(light_index),
        VertexKind::Surface { primitive_index } => {
            lighting.rect_index_for_primitive(primitive_index)
        }
        VertexKind::Camera => None,
    }
}

fn pdf_light_origin(lighting: &AdmittedLighting<'_>, light: &Vertex, next: &Vertex) -> f64 {
    let Some(light_index) = light_identity(lighting, light) else {
        return 0.0;
    };
    let Some((direction, _)) = vector_between(light.point, next.point) else {
        return 0.0;
    };
    lighting
        .rectangle_emission_pdfs(light_index, direction)
        .map_or(0.0, |(position_pdf, _)| position_pdf)
}

fn pdf_light_direction(lighting: &AdmittedLighting<'_>, light: &Vertex, next: &Vertex) -> f64 {
    let Some(light_index) = light_identity(lighting, light) else {
        return 0.0;
    };
    let Some((direction, _)) = vector_between(light.point, next.point) else {
        return 0.0;
    };
    let Some((_, direction_pdf)) = lighting.rectangle_emission_pdfs(light_index, direction) else {
        return 0.0;
    };
    convert_density(light, direction_pdf, next)
}

fn endpoint_primitive(vertex: &Vertex) -> Option<usize> {
    match vertex.kind {
        VertexKind::Light {
            primitive_index, ..
        }
        | VertexKind::Surface { primitive_index } => Some(primitive_index),
        VertexKind::Camera => None,
    }
}

fn segment_transmittance(
    scene: &Scene,
    cx: &Cx<'_>,
    ray_time: Option<&PathTime>,
    from: &Vertex,
    to: &Vertex,
    wavelength_nm: f64,
) -> Result<f64, TracerError> {
    let Some((direction, distance_squared)) = vector_between(from.point, to.point) else {
        return Ok(0.0);
    };
    let reverse = direction.scale(-1.0);
    let from_medium = from.medium_toward(direction);
    let to_medium = to.medium_toward(reverse);
    if !same_medium(from_medium, to_medium) {
        return Ok(0.0);
    }
    let origin = from.geometric_normal.map_or(from.point, |normal| {
        dielectric_spawn_origin(from.point, normal, direction)
    });
    let ray = Ray {
        origin,
        dir: direction,
    };
    let target_primitive = endpoint_primitive(to);
    let visible = match intersect(scene, cx, &ray, ray_time)? {
        None => target_primitive.is_none(),
        Some(hit) => {
            let distance = distance_squared.sqrt();
            target_primitive == Some(hit.primitive_index)
                && hit.hit.point.delta_from(to.point).norm() <= 8.0 * RAY_EPS + 2.0e-8 * distance
        }
    };
    if !visible {
        return Ok(0.0);
    }
    let distance = distance_squared.sqrt();
    medium_transmittance(
        from_medium.map(|entry| entry.glass),
        wavelength_nm,
        distance,
    )
}

fn environment_visible(
    scene: &Scene,
    cx: &Cx<'_>,
    ray_time: Option<&PathTime>,
    from: &Vertex,
    direction: Vec3,
) -> Result<bool, TracerError> {
    if from.medium_toward(direction).is_some() {
        return Ok(false);
    }
    let origin = from.geometric_normal.map_or(from.point, |normal| {
        dielectric_spawn_origin(from.point, normal, direction)
    });
    Ok(intersect(
        scene,
        cx,
        &Ray {
            origin,
            dir: direction,
        },
        ray_time,
    )?
    .is_none())
}

fn geometry_term(
    scene: &Scene,
    cx: &Cx<'_>,
    ray_time: Option<&PathTime>,
    left: &Vertex,
    right: &Vertex,
    wavelength_nm: f64,
) -> Result<f64, TracerError> {
    let Some((direction, distance_squared)) = vector_between(left.point, right.point) else {
        return Ok(0.0);
    };
    let left_cosine = left
        .geometric_normal
        .map_or(1.0, |normal| normal.dot(direction).abs());
    let right_cosine = right
        .geometric_normal
        .map_or(1.0, |normal| normal.dot(direction.scale(-1.0)).abs());
    if left_cosine <= 0.0 || right_cosine <= 0.0 {
        return Ok(0.0);
    }
    Ok(left_cosine * right_cosine / distance_squared
        * segment_transmittance(scene, cx, ray_time, left, right, wavelength_nm)?)
}

#[allow(clippy::too_many_arguments)]
fn finite_mis_weight(
    scene: &Scene,
    lighting: &AdmittedLighting<'_>,
    camera: &PhysicalCamera,
    light_vertices: &[Vertex],
    camera_vertices: &[Vertex],
    sampled: Option<Vertex>,
    s: usize,
    t: usize,
    wavelength_nm: f64,
    width: u32,
    height: u32,
) -> Result<f64, TracerError> {
    if s + t == 2 {
        return Ok(1.0);
    }
    let mut light = light_vertices.to_vec();
    let mut camera_path = camera_vertices.to_vec();
    if s == 1 {
        light[0] = sampled.ok_or(TracerError::InvalidInput)?;
    } else if t == 1 {
        camera_path[0] = sampled.ok_or(TracerError::InvalidInput)?;
    }

    let qs_index = s.checked_sub(1);
    let pt_index = t.checked_sub(1);
    let qs_minus_index = s.checked_sub(2);
    let pt_minus_index = t.checked_sub(2);
    if let Some(index) = pt_index {
        camera_path[index].delta = false;
    }
    if let Some(index) = qs_index {
        light[index].delta = false;
    }

    if let Some(pt) = pt_index {
        camera_path[pt].pdf_rev = if let Some(qs) = qs_index {
            vertex_pdf(
                scene,
                lighting,
                camera,
                &light[qs],
                qs_minus_index.map(|index| &light[index]),
                &camera_path[pt],
                wavelength_nm,
                width,
                height,
            )?
        } else {
            let Some(pt_minus) = pt_minus_index else {
                return Ok(1.0);
            };
            pdf_light_origin(lighting, &camera_path[pt], &camera_path[pt_minus])
        };
    }
    if let Some(pt_minus) = pt_minus_index {
        camera_path[pt_minus].pdf_rev = if let Some(qs) = qs_index {
            vertex_pdf(
                scene,
                lighting,
                camera,
                &camera_path[pt_index.unwrap()],
                Some(&light[qs]),
                &camera_path[pt_minus],
                wavelength_nm,
                width,
                height,
            )?
        } else {
            pdf_light_direction(
                lighting,
                &camera_path[pt_index.unwrap()],
                &camera_path[pt_minus],
            )
        };
    }
    if let Some(qs) = qs_index {
        let pt = pt_index.ok_or(TracerError::InvalidInput)?;
        light[qs].pdf_rev = vertex_pdf(
            scene,
            lighting,
            camera,
            &camera_path[pt],
            pt_minus_index.map(|index| &camera_path[index]),
            &light[qs],
            wavelength_nm,
            width,
            height,
        )?;
    }
    if let Some(qs_minus) = qs_minus_index {
        light[qs_minus].pdf_rev = vertex_pdf(
            scene,
            lighting,
            camera,
            &light[qs_index.unwrap()],
            pt_index.map(|index| &camera_path[index]),
            &light[qs_minus],
            wavelength_nm,
            width,
            height,
        )?;
    }

    let remap_zero = |pdf: f64| if pdf == 0.0 { 1.0 } else { pdf };
    // A t=1 strategy is sampled once for every source pixel and splatted into
    // the pixel selected by the camera endpoint.  All other strategies are
    // sampled once for the target pixel.  Multi-sample MIS therefore uses the
    // global t=1 sample count in addition to the path-density ratio.  Omitting
    // this factor makes t=1 both misweighted and resolution dependent.
    let mut sum_ratio = 0.0;
    let mut ratio = 1.0;
    for index in (1..t).rev() {
        ratio *= remap_zero(camera_path[index].pdf_rev) / remap_zero(camera_path[index].pdf_fwd);
        if !camera_path[index].delta && !camera_path[index - 1].delta {
            sum_ratio += ratio_with_sample_multiplicity(ratio, index, t, width, height)?;
        }
    }
    ratio = 1.0;
    for index in (0..s).rev() {
        ratio *= remap_zero(light[index].pdf_rev) / remap_zero(light[index].pdf_fwd);
        let preceding_delta = index > 0 && light[index - 1].delta;
        if !light[index].delta && !preceding_delta {
            let alternative_t = t + (s - index);
            sum_ratio += ratio_with_sample_multiplicity(ratio, alternative_t, t, width, height)?;
        }
    }
    let weight = 1.0 / (1.0 + sum_ratio);
    if weight.is_finite() && (0.0..=1.0).contains(&weight) {
        Ok(weight)
    } else {
        Err(TracerError::InvalidInput)
    }
}

fn ratio_with_sample_multiplicity(
    density_ratio: f64,
    alternative_t: usize,
    current_t: usize,
    width: u32,
    height: u32,
) -> Result<f64, TracerError> {
    let alternative_samples = strategy_sample_multiplicity(alternative_t, width, height)?;
    let current_samples = strategy_sample_multiplicity(current_t, width, height)?;
    Ok(density_ratio * alternative_samples / current_samples)
}

fn strategy_sample_multiplicity(t: usize, width: u32, height: u32) -> Result<f64, TracerError> {
    if t != 1 {
        return Ok(1.0);
    }
    let pixels = width.checked_mul(height).ok_or(TracerError::InvalidInput)?;
    if pixels == 0 {
        return Err(TracerError::InvalidInput);
    }
    Ok(f64::from(pixels))
}

fn scalar_to_xyz(value: f64, wavelength_nm: f64) -> [f64; 3] {
    let weight = value * (LAMBDA_MAX - LAMBDA_MIN) / y_integral();
    [
        weight * cie_x(wavelength_nm),
        weight * cie_y(wavelength_nm),
        weight * cie_z(wavelength_nm),
    ]
}

fn add_scaled_xyz(target: &mut [f64; 3], value: [f64; 3]) {
    target[0] += value[0];
    target[1] += value[1];
    target[2] += value[2];
}

#[allow(clippy::too_many_arguments)]
fn connect_s1(
    scene: &Scene,
    lighting: &AdmittedLighting<'_>,
    camera: &PhysicalCamera,
    cx: &Cx<'_>,
    ray_time: Option<&PathTime>,
    settings: &Settings,
    light_vertices: &[Vertex],
    camera_vertices: &[Vertex],
    t: usize,
    wavelength_nm: f64,
    random: [f64; 4],
) -> Result<Option<f64>, TracerError> {
    let pt = &camera_vertices[t - 1];
    if !pt.is_connectible(scene) || !pt.is_surface() {
        return Ok(None);
    }
    let Some(light_sample) = lighting.sample(pt.point, random[0], random[1]) else {
        return Ok(None);
    };
    let predecessor = &camera_vertices[t - 2];
    let predecessor_direction = vector_between(pt.point, predecessor.point)
        .map(|(direction, _)| direction)
        .ok_or(TracerError::InvalidInput)?;
    match light_sample {
        LightSample::Environment(sample) => {
            let f = surface_scattering(
                scene,
                pt,
                predecessor_direction,
                sample.direction,
                wavelength_nm,
                TransportMode::Radiance,
            )?;
            let cosine = pt
                .geometric_normal
                .map_or(1.0, |normal| normal.dot(sample.direction).abs());
            if f <= 0.0
                || cosine <= 0.0
                || !environment_visible(scene, cx, ray_time, pt, sample.direction)?
            {
                return Ok(None);
            }
            let bsdf_pdf = surface_pdf_solid_angle(
                scene,
                pt,
                predecessor_direction,
                sample.direction,
                wavelength_nm,
            )?;
            let weight = balance_heuristic(1, sample.pdf_solid_angle, 1, bsdf_pdf);
            let (spectrum, scale) = sample.emission;
            let value = pt.beta * f * cosine * spectrum.eval(wavelength_nm) * scale
                / sample.pdf_solid_angle
                * weight;
            Ok((value.is_finite() && value > 0.0).then_some(value))
        }
        LightSample::Rectangle(sample) => {
            let direction = sample.point.delta_from(pt.point);
            let distance_squared = direction.dot(direction);
            if !distance_squared.is_finite() || distance_squared <= 0.0 {
                return Ok(None);
            }
            let direction = direction.scale(1.0 / distance_squared.sqrt());
            let f = surface_scattering(
                scene,
                pt,
                predecessor_direction,
                direction,
                wavelength_nm,
                TransportMode::Radiance,
            )?;
            let cosine = pt
                .geometric_normal
                .map_or(1.0, |normal| normal.dot(direction).abs());
            if f <= 0.0 || cosine <= 0.0 {
                return Ok(None);
            }
            let emitted = sample.emission.0.eval(wavelength_nm) * sample.emission.1;
            let position_pdf = lighting
                .rectangle_emission_pdfs(sample.light_index, direction.scale(-1.0))
                .map_or(0.0, |(pdf, _)| pdf);
            let sampled = Vertex {
                kind: VertexKind::Light {
                    light_index: sample.light_index,
                    primitive_index: sample.primitive_index,
                    emission: sample.emission,
                },
                point: sample.point,
                geometric_normal: Some(sample.normal),
                wo: None,
                beta: emitted / sample.pdf_solid_angle,
                pdf_fwd: position_pdf,
                pdf_rev: 0.0,
                delta: false,
                positive_medium: None,
                negative_medium: None,
            };
            let transmittance =
                segment_transmittance(scene, cx, ray_time, pt, &sampled, wavelength_nm)?;
            if transmittance <= 0.0 {
                return Ok(None);
            }
            let mis = finite_mis_weight(
                scene,
                lighting,
                camera,
                light_vertices,
                camera_vertices,
                Some(sampled),
                1,
                t,
                wavelength_nm,
                settings.width,
                settings.height,
            )?;
            let value =
                pt.beta * f * cosine * emitted / sample.pdf_solid_angle * transmittance * mis;
            Ok((value.is_finite() && value > 0.0).then_some(value))
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn connect_t1(
    scene: &Scene,
    lighting: &AdmittedLighting<'_>,
    camera: &PhysicalCamera,
    cx: &Cx<'_>,
    ray_time: Option<&PathTime>,
    settings: &Settings,
    light_vertices: &[Vertex],
    camera_vertices: &[Vertex],
    s: usize,
    wavelength_nm: f64,
) -> Result<Option<(u32, f64)>, TracerError> {
    let qs = &light_vertices[s - 1];
    if !qs.is_connectible(scene) || !qs.is_surface() {
        return Ok(None);
    }
    let Some(raster) = camera.pinhole_raster_sample(qs.point, settings.width, settings.height)?
    else {
        return Ok(None);
    };
    let camera_vertex = Vertex {
        pdf_fwd: 0.0,
        ..Vertex::camera(camera.eye())
    };
    let predecessor = &light_vertices[s - 2];
    let predecessor_direction = vector_between(qs.point, predecessor.point)
        .map(|(direction, _)| direction)
        .ok_or(TracerError::InvalidInput)?;
    let direction_to_camera = raster.direction_from_camera.scale(-1.0);
    let f = surface_scattering(
        scene,
        qs,
        predecessor_direction,
        direction_to_camera,
        wavelength_nm,
        TransportMode::Importance,
    )?;
    if f <= 0.0 {
        return Ok(None);
    }
    let camera_importance_area = convert_density(&camera_vertex, raster.pdf_solid_angle, qs);
    if camera_importance_area <= 0.0 {
        return Ok(None);
    }
    let transmittance =
        segment_transmittance(scene, cx, ray_time, qs, &camera_vertex, wavelength_nm)?;
    if transmittance <= 0.0 {
        return Ok(None);
    }
    let mis = finite_mis_weight(
        scene,
        lighting,
        camera,
        light_vertices,
        camera_vertices,
        Some(camera_vertex),
        s,
        1,
        wavelength_nm,
        settings.width,
        settings.height,
    )?;
    // The complete raster launches width*height independent light subpaths
    // for every sample index.  A target pixel receives arbitrary splats from
    // that global pool, so each splat carries the ordinary 1/n estimator
    // normalization for the t=1 strategy.
    let light_path_samples = strategy_sample_multiplicity(1, settings.width, settings.height)?;
    let value = qs.beta * f * camera_importance_area * transmittance * mis / light_path_samples;
    Ok((value.is_finite() && value > 0.0).then_some((raster.pixel, value)))
}

#[allow(clippy::too_many_arguments)]
fn connect_general(
    scene: &Scene,
    lighting: &AdmittedLighting<'_>,
    camera: &PhysicalCamera,
    cx: &Cx<'_>,
    ray_time: Option<&PathTime>,
    settings: &Settings,
    light_vertices: &[Vertex],
    camera_vertices: &[Vertex],
    s: usize,
    t: usize,
    wavelength_nm: f64,
) -> Result<Option<f64>, TracerError> {
    let qs = &light_vertices[s - 1];
    let pt = &camera_vertices[t - 1];
    if !qs.is_connectible(scene)
        || !pt.is_connectible(scene)
        || !qs.is_surface()
        || !pt.is_surface()
    {
        return Ok(None);
    }
    let qs_predecessor = &light_vertices[s - 2];
    let pt_predecessor = &camera_vertices[t - 2];
    let qs_wo = vector_between(qs.point, qs_predecessor.point)
        .map(|(direction, _)| direction)
        .ok_or(TracerError::InvalidInput)?;
    let pt_wo = vector_between(pt.point, pt_predecessor.point)
        .map(|(direction, _)| direction)
        .ok_or(TracerError::InvalidInput)?;
    let Some((qs_to_pt, _)) = vector_between(qs.point, pt.point) else {
        return Ok(None);
    };
    let qs_f = surface_scattering(
        scene,
        qs,
        qs_wo,
        qs_to_pt,
        wavelength_nm,
        TransportMode::Importance,
    )?;
    let pt_f = surface_scattering(
        scene,
        pt,
        pt_wo,
        qs_to_pt.scale(-1.0),
        wavelength_nm,
        TransportMode::Radiance,
    )?;
    if qs_f <= 0.0 || pt_f <= 0.0 {
        return Ok(None);
    }
    let geometry = geometry_term(scene, cx, ray_time, qs, pt, wavelength_nm)?;
    if geometry <= 0.0 {
        return Ok(None);
    }
    let mis = finite_mis_weight(
        scene,
        lighting,
        camera,
        light_vertices,
        camera_vertices,
        None,
        s,
        t,
        wavelength_nm,
        settings.width,
        settings.height,
    )?;
    let value = qs.beta * qs_f * pt_f * pt.beta * geometry * mis;
    Ok((value.is_finite() && value > 0.0).then_some(value))
}

#[allow(clippy::too_many_arguments)]
fn evaluate_strategies(
    scene: &Scene,
    lighting: &AdmittedLighting<'_>,
    camera: &PhysicalCamera,
    cx: &Cx<'_>,
    ray_time: Option<&PathTime>,
    settings: &Settings,
    source_pixel: u32,
    sample: u32,
    wavelength_nm: f64,
    path_sobol: Option<&Sobol>,
    light_vertices: &[Vertex],
    camera_subpath: &CameraSubpath,
    film: &mut Film,
    splats: &mut Vec<SplatRecord>,
    stats: &mut BidirectionalStrategyStats,
) -> Result<(), TracerError> {
    let camera_vertices = &camera_subpath.vertices;
    if let Some(escape) = &camera_subpath.escape
        && let Some(environment) = lighting.environment_evaluation(escape.origin, escape.direction)
    {
        let weight = if escape.previous_delta {
            1.0
        } else {
            balance_heuristic(
                1,
                escape.previous_pdf_solid_angle,
                1,
                environment.pdf_solid_angle,
            )
        };
        let emitted = environment.emission.0.eval(wavelength_nm) * environment.emission.1;
        let value = escape.beta * emitted * weight;
        if value.is_finite() && value > 0.0 {
            add_scaled_xyz(
                &mut film.xyz[source_pixel as usize],
                scalar_to_xyz(value, wavelength_nm),
            );
            stats.nonzero = stats.nonzero.saturating_add(1);
        }
    }

    for t in 1..=camera_vertices.len() {
        for s in 0..=light_vertices.len() {
            let Some(depth) = s.checked_add(t).and_then(|sum| sum.checked_sub(2)) else {
                continue;
            };
            if (s == 1 && t == 1) || depth > settings.max_depth as usize {
                continue;
            }
            stats.evaluated = stats.evaluated.saturating_add(1);
            let mut connection_rng = SubpathRng::new(
                settings,
                source_pixel,
                sample,
                CONNECTION_DOMAIN,
                path_sobol,
            );
            connection_rng.dimension = ((s as u32) << 16) ^ t as u32;
            let contribution = if s == 0 {
                if t < 2 {
                    None
                } else {
                    let pt = camera_vertices[t - 1];
                    let emitted = pt.emitted_radiance(scene, wavelength_nm);
                    if emitted <= 0.0 {
                        None
                    } else {
                        let mis = finite_mis_weight(
                            scene,
                            lighting,
                            camera,
                            light_vertices,
                            camera_vertices,
                            None,
                            0,
                            t,
                            wavelength_nm,
                            settings.width,
                            settings.height,
                        )?;
                        Some((source_pixel, pt.beta * emitted * mis, false))
                    }
                }
            } else if t == 1 {
                connect_t1(
                    scene,
                    lighting,
                    camera,
                    cx,
                    ray_time,
                    settings,
                    light_vertices,
                    camera_vertices,
                    s,
                    wavelength_nm,
                )?
                .map(|(pixel, value)| (pixel, value, true))
            } else if s == 1 {
                connect_s1(
                    scene,
                    lighting,
                    camera,
                    cx,
                    ray_time,
                    settings,
                    light_vertices,
                    camera_vertices,
                    t,
                    wavelength_nm,
                    connection_rng.next4(),
                )?
                .map(|value| (source_pixel, value, false))
            } else {
                connect_general(
                    scene,
                    lighting,
                    camera,
                    cx,
                    ray_time,
                    settings,
                    light_vertices,
                    camera_vertices,
                    s,
                    t,
                    wavelength_nm,
                )?
                .map(|value| (source_pixel, value, false))
            };
            if let Some((pixel, value, splat)) = contribution
                && value.is_finite()
                && value > 0.0
            {
                let xyz = scalar_to_xyz(value, wavelength_nm);
                if splat {
                    splats.push(SplatRecord {
                        target_pixel: pixel,
                        sample,
                        strategy_s: s,
                        xyz,
                    });
                } else {
                    add_scaled_xyz(&mut film.xyz[pixel as usize], xyz);
                }
                stats.nonzero = stats.nonzero.saturating_add(1);
                if splat {
                    stats.camera_splats = stats.camera_splats.saturating_add(1);
                }
            }
        }
    }
    Ok(())
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn render_admitted_sample_range(
    scene: &Scene,
    lighting: &AdmittedLighting<'_>,
    camera: &AnimatedCamera,
    exposure: CameraExposure<'_>,
    cx: &Cx<'_>,
    settings: &Settings,
    shutter: ShutterInterval,
    from: u32,
    to: u32,
    sobol: Option<&Sobol>,
    path_sobol: Option<&Sobol>,
) -> Result<BidirectionalRenderOutput, TracerError> {
    if from > to || to > settings.spp {
        return Err(TracerError::InvalidInput);
    }
    let mut film = Film::try_new(settings.width, settings.height)?;
    let mut stats = BidirectionalStrategyStats::default();
    let key = [
        (settings.seed & 0xffff_ffff) as u32,
        (settings.seed >> 32) as u32,
    ];
    let render_cx = cx.with_stream_seed(settings.seed);
    let pixel_count = settings
        .width
        .checked_mul(settings.height)
        .ok_or(TracerError::InvalidInput)?;
    // Sample-major traversal is load-bearing. A complete light subpath is
    // launched for every `(sample, source_pixel)`, and all of its arbitrary
    // raster splats publish before the next source identity. Fixed sample
    // blocks can therefore execute in any worker order and still be merged in
    // one canonical arithmetic order.
    for sample in from..to {
        cx.checkpoint()?;
        for pixel in 0..pixel_count {
            cx.checkpoint()?;
            let mut source_splats = Vec::new();
            let (jx, jy, wavelength_sample) = pixel_dims(settings, sobol, key, pixel, sample)?;
            let wavelength_nm = LAMBDA_MIN + wavelength_sample * (LAMBDA_MAX - LAMBDA_MIN);
            let ray_time =
                path_time_for_sample(scene, &render_cx, settings, shutter, pixel, sample)?;
            let absolute_time_s = shutter.time_at(ray_time.normalized);
            let physical = camera.evaluate_exposure(&render_cx, exposure, absolute_time_s)?;
            if !physical.aperture().is_pinhole() {
                return Err(TracerError::InvalidInput);
            }
            let px = pixel % settings.width;
            let py = pixel / settings.width;
            let width = f64::from(settings.width);
            let height = f64::from(settings.height);
            let half_tan = physical.projection().vertical_half_tan();
            let aspect = width / height;
            let x_tan = (2.0 * (f64::from(px) + jx) / width - 1.0) * aspect * half_tan;
            let y_tan = (1.0 - 2.0 * (f64::from(py) + jy) / height) * half_tan;
            let ray = physical.generate_ray_from_tangent_offsets(
                &render_cx,
                x_tan,
                y_tan,
                camera_lens_sample(key, pixel, sample)?,
            )?;
            let camera_subpath = generate_camera_subpath(
                scene,
                &render_cx,
                settings,
                &physical,
                ray,
                Some(&ray_time),
                wavelength_nm,
                pixel,
                sample,
                path_sobol,
            )?;
            let light_subpath = generate_light_subpath(
                scene,
                &lighting,
                &render_cx,
                settings,
                Some(&ray_time),
                wavelength_nm,
                pixel,
                sample,
                path_sobol,
            )?;
            evaluate_strategies(
                scene,
                &lighting,
                &physical,
                &render_cx,
                Some(&ray_time),
                settings,
                pixel,
                sample,
                wavelength_nm,
                path_sobol,
                &light_subpath,
                &camera_subpath,
                &mut film,
                &mut source_splats,
                &mut stats,
            )?;
            publish_source_splats(&mut film, &mut source_splats);
        }
    }
    film.spp_done = to - from;
    Ok(BidirectionalRenderOutput {
        film,
        strategies: stats,
    })
}

fn merge_block(
    film: &mut Film,
    strategies: &mut BidirectionalStrategyStats,
    block: BidirectionalRenderOutput,
) -> Result<(), TracerError> {
    if film.width != block.film.width
        || film.height != block.film.height
        || film.xyz.len() != block.film.xyz.len()
    {
        return Err(TracerError::InvalidInput);
    }
    for (target, source) in film.xyz.iter_mut().zip(block.film.xyz) {
        add_scaled_xyz(target, source);
    }
    film.spp_done = film.spp_done.saturating_add(block.film.spp_done);
    strategies.evaluated = strategies
        .evaluated
        .saturating_add(block.strategies.evaluated);
    strategies.nonzero = strategies.nonzero.saturating_add(block.strategies.nonzero);
    strategies.camera_splats = strategies
        .camera_splats
        .saturating_add(block.strategies.camera_splats);
    Ok(())
}

/// Render a fixed-SPP cinematic frame with bidirectional path tracing.
///
/// The serial oracle and parallel executor share the same fixed eight-sample
/// reduction blocks. This keeps floating-point image bits invariant to worker
/// count while supporting `t=1` light-tracing strategies that splat into
/// arbitrary pixels. The estimator supports all current opaque and
/// smooth/rough dielectric materials, finite rectangle emitters, animated
/// geometry at sampled shutter time, and a pinhole cinematic camera. A finite
/// aperture refuses rather than pretending an optical-centre connection
/// represents the lens integral.
#[allow(clippy::too_many_arguments)]
pub fn render_cinematic_bidirectional(
    scene: &Scene,
    camera: &AnimatedCamera,
    cut_side: CutSide,
    cx: &Cx<'_>,
    settings: &Settings,
    shutter: ShutterInterval,
) -> Result<BidirectionalRenderOutput, TracerError> {
    cx.checkpoint()?;
    let exposure = camera.admit_shutter(cx, shutter, cut_side)?;
    let camera_path = CameraPath::Cinematic { camera, exposure };
    let (lighting, time_mode) = preflight_render(
        scene,
        cx,
        settings,
        None,
        0,
        settings.spp,
        Some(shutter),
        camera_path,
    )?;
    let sobol = pixel_sobol(settings.sampler, settings.seed);
    let path_sobol = path_sobol(settings.sampler);
    let mut film = Film::try_new(settings.width, settings.height)?;
    let mut strategies = BidirectionalStrategyStats::default();
    let mut from = 0;
    while from < settings.spp {
        let to = from
            .saturating_add(BIDIRECTIONAL_SAMPLE_BLOCK)
            .min(settings.spp);
        let block = render_admitted_sample_range(
            scene,
            &lighting,
            camera,
            exposure,
            cx,
            settings,
            shutter,
            from,
            to,
            sobol.as_ref(),
            path_sobol.as_ref(),
        )?;
        merge_block(&mut film, &mut strategies, block)?;
        from = to;
    }
    film.time_mode = time_mode;
    Ok(BidirectionalRenderOutput { film, strategies })
}

struct BidirectionalSampleKernel<'run, 'assets> {
    scene: &'assets Scene,
    lighting: &'run AdmittedLighting<'assets>,
    camera: &'assets AnimatedCamera,
    exposure: CameraExposure<'assets>,
    settings: &'run Settings,
    shutter: ShutterInterval,
    sobol: Option<&'run Sobol>,
    path_sobol: Option<&'run Sobol>,
    results: &'run Mutex<Vec<Option<BidirectionalRenderOutput>>>,
    failures: &'run Mutex<Option<(u64, RenderTileFailure)>>,
    sample_blocks: u64,
    compute_ns: &'run AtomicU64,
}

impl BidirectionalSampleKernel<'_, '_> {
    fn fail(&self, block: u64, failure: RenderTileFailure) -> ControlFlow<Cancelled, ()> {
        let mut recorded = self
            .failures
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if recorded
            .as_ref()
            .is_none_or(|(recorded_block, _)| block < *recorded_block)
        {
            *recorded = Some((block, failure));
        }
        ControlFlow::Break(Cancelled)
    }
}

impl TileKernel for BidirectionalSampleKernel<'_, '_> {
    type Out = ();

    fn tiles(&self) -> TilePlan {
        TilePlan::new(BIDIRECTIONAL_SAMPLE_KERNEL, self.sample_blocks)
    }

    fn run(&self, block: u64, cx: &Cx<'_>) -> ControlFlow<Cancelled, Self::Out> {
        if cx.checkpoint().is_err() {
            return ControlFlow::Break(Cancelled);
        }
        let Some(from) = u32::try_from(block)
            .ok()
            .and_then(|block| block.checked_mul(BIDIRECTIONAL_SAMPLE_BLOCK))
        else {
            return self.fail(
                block,
                RenderTileFailure::Internal("BDPT sample-block identity overflow"),
            );
        };
        let to = from
            .saturating_add(BIDIRECTIONAL_SAMPLE_BLOCK)
            .min(self.settings.spp);
        let started = Instant::now();
        let result = render_admitted_sample_range(
            self.scene,
            self.lighting,
            self.camera,
            self.exposure,
            cx,
            self.settings,
            self.shutter,
            from,
            to,
            self.sobol,
            self.path_sobol,
        );
        atomic_saturating_add(self.compute_ns, elapsed_ns(started));
        let result = match result {
            Ok(result) => result,
            Err(TracerError::Cancelled) => return ControlFlow::Break(Cancelled),
            Err(error) => return self.fail(block, RenderTileFailure::Tracer(error)),
        };
        let mut results = self
            .results
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        let Some(slot) = usize::try_from(block)
            .ok()
            .and_then(|index| results.get_mut(index))
        else {
            return self.fail(
                block,
                RenderTileFailure::Internal("BDPT sample block outside staging plan"),
            );
        };
        if slot.is_some() {
            return self.fail(
                block,
                RenderTileFailure::Internal("BDPT sample block executed more than once"),
            );
        }
        *slot = Some(result);
        ControlFlow::Continue(())
    }
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(super) fn render_cinematic_bidirectional_execution_impl<R: RenderPoolRunner>(
    scene: &Scene,
    camera: &AnimatedCamera,
    cut_side: CutSide,
    cx: &Cx<'_>,
    settings: &Settings,
    shutter: ShutterInterval,
    execution: &RenderExecutionConfig,
    runner: &R,
) -> Result<BidirectionalRenderExecutionOutput, RenderExecutionError> {
    let setup_started = Instant::now();
    cx.checkpoint()?;
    let exposure = camera
        .admit_shutter(cx, shutter, cut_side)
        .map_err(TracerError::from)?;
    let camera_path = CameraPath::Cinematic { camera, exposure };
    let (lighting, time_mode) = preflight_render(
        scene,
        cx,
        settings,
        None,
        0,
        settings.spp,
        Some(shutter),
        camera_path,
    )?;
    let image_layout = RenderTileLayout::try_new(
        settings.width,
        settings.height,
        execution.tile_width,
        execution.tile_height,
    )
    .map_err(RenderExecutionError::Config)?;
    let sample_blocks = u64::from(settings.spp.div_ceil(BIDIRECTIONAL_SAMPLE_BLOCK));
    let active_worker_ceiling = u64::try_from(execution.workers)
        .unwrap_or(u64::MAX)
        .min(sample_blocks);
    let pixel_count = checked_pixel_len(settings.width, settings.height)?;
    let film_bytes = u64::try_from(pixel_count)
        .ok()
        .and_then(|count| count.checked_mul(size_of::<[f64; 3]>() as u64))
        .ok_or(RenderExecutionError::Config(
            RenderExecutionConfigError::InvalidImageDimensions,
        ))?;
    // Every completed block stays private until executor drain; publication
    // then owns one additional image-sized film while reducing those blocks.
    let staging_film_bytes = sample_blocks
        .checked_add(1)
        .and_then(|count| count.checked_mul(film_bytes))
        .ok_or(RenderExecutionError::Internal(
            "BDPT staging-film byte accounting overflow",
        ))?;
    let maximum_splats_per_source = u64::from(BIDIRECTIONAL_SAMPLE_BLOCK)
        .checked_mul(u64::from(settings.max_depth).saturating_add(1))
        .ok_or(RenderExecutionError::Internal(
            "BDPT splat scratch count overflow",
        ))?;
    let splat_scratch_envelope_bytes = maximum_splats_per_source
        .checked_mul(size_of::<SplatRecord>() as u64)
        .and_then(|bytes| bytes.checked_mul(active_worker_ceiling))
        .ok_or(RenderExecutionError::Internal(
            "BDPT splat scratch byte accounting overflow",
        ))?;
    let sampler_state_bytes = sampler_direction_bytes(settings.sampler);
    let lease = OperationMemoryLease::bounded(execution.memory_limit_bytes);
    let staging_charge = lease
        .reserve("bdpt-sample-block-films", staging_film_bytes)
        .map_err(RenderExecutionError::Memory)?;
    let splat_charge = lease
        .reserve("bdpt-source-splat-scratch", splat_scratch_envelope_bytes)
        .map_err(RenderExecutionError::Memory)?;
    let sampler_charge = (sampler_state_bytes != 0)
        .then(|| lease.reserve("bdpt-sobol-directions", sampler_state_bytes))
        .transpose()
        .map_err(RenderExecutionError::Memory)?;
    let mut staged = Vec::new();
    let sample_block_len = usize::try_from(sample_blocks)
        .map_err(|_| RenderExecutionError::Internal("BDPT sample-block count exceeds usize"))?;
    staged.try_reserve_exact(sample_block_len).map_err(|_| {
        RenderExecutionError::Allocation(AllocError::OutOfMemory {
            site: "bdpt-sample-block-index",
            requested_bytes: sample_block_len
                .saturating_mul(size_of::<Option<BidirectionalRenderOutput>>()),
        })
    })?;
    staged.resize_with(sample_block_len, || None);
    let sobol = pixel_sobol(settings.sampler, settings.seed);
    let path_sobol = path_sobol(settings.sampler);
    let results = Mutex::new(staged);
    let failures = Mutex::new(None);
    let compute_ns = AtomicU64::new(0);
    let kernel = BidirectionalSampleKernel {
        scene,
        lighting: &lighting,
        camera,
        exposure,
        settings,
        shutter,
        sobol: sobol.as_ref(),
        path_sobol: path_sobol.as_ref(),
        results: &results,
        failures: &failures,
        sample_blocks,
        compute_ns: &compute_ns,
    };
    let setup_ns = elapsed_ns(setup_started);
    let traversal_started = Instant::now();
    let (outcome, executor) = runner.run_render(cx, &kernel, execution.run_id, &lease);
    let traversal_ns = elapsed_ns(traversal_started);
    let failure = failures
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner)
        .clone();
    if let Some(error) = render_outcome_error(outcome, failure) {
        return Err(error);
    }
    cx.checkpoint()?;
    drop(kernel);
    let mut blocks = results.into_inner().map_err(|_| {
        RenderExecutionError::Internal("successful BDPT staging mutex was poisoned")
    })?;
    let publication_started = Instant::now();
    let mut film = Film::try_new(settings.width, settings.height)?;
    let mut strategies = BidirectionalStrategyStats::default();
    for block in &mut blocks {
        let block = block.take().ok_or(RenderExecutionError::Internal(
            "successful BDPT run omitted a sample block",
        ))?;
        merge_block(&mut film, &mut strategies, block)?;
    }
    film.time_mode = time_mode;
    let publication_ns = elapsed_ns(publication_started);
    let block_compute_ns = compute_ns.load(Ordering::Relaxed);
    let active_workers = executor.tiles_by_worker.len();
    let idle_worker_ns = traversal_ns
        .saturating_mul(active_workers as u64)
        .saturating_sub(block_compute_ns);
    drop(blocks);
    drop(sobol);
    drop(path_sobol);
    drop(sampler_charge);
    drop(splat_charge);
    drop(staging_charge);
    let memory = lease.receipt();
    Ok(BidirectionalRenderExecutionOutput {
        film,
        strategies,
        report: BidirectionalExecutionReport {
            image_layout,
            sample_blocks,
            samples_per_block: BIDIRECTIONAL_SAMPLE_BLOCK,
            requested_workers: execution.workers,
            workers: active_workers,
            staging_film_bytes,
            splat_scratch_envelope_bytes,
            sampler_state_bytes,
            setup_ns,
            traversal_ns,
            block_compute_ns,
            publication_ns,
            idle_worker_ns,
            executor,
            memory,
        },
    })
}

/// Render a complete cinematic BDPT frame with fixed sample-block parallelism.
#[allow(clippy::too_many_arguments)]
pub fn render_cinematic_bidirectional_with_execution(
    scene: &Scene,
    camera: &AnimatedCamera,
    cut_side: CutSide,
    cx: &Cx<'_>,
    settings: &Settings,
    shutter: ShutterInterval,
    execution: &RenderExecutionConfig,
) -> Result<BidirectionalRenderExecutionOutput, RenderExecutionError> {
    let pool = build_render_pool(execution, cx.mode(), settings.seed);
    render_cinematic_bidirectional_execution_impl(
        scene, camera, cut_side, cx, settings, shutter, execution, &pool,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::camera::Aperture;
    use crate::motion::{ShotTimeBounds, ShutterConvention, ShutterDistribution};
    use crate::spectral::lift_rgb;
    use fs_exec::{CancelGate, StreamKey};

    fn with_test_cx<R>(operation: impl FnOnce(&Cx<'_>) -> R) -> R {
        let gate = CancelGate::new();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 0x4244_5054,
                    kernel_id: 1,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            operation(&cx)
        })
    }

    fn emitter_quad(z: f64) -> TriMesh {
        TriMesh::new(
            vec![
                [-1.0, -1.0, z],
                [1.0, -1.0, z],
                [1.0, 1.0, z],
                [-1.0, 1.0, z],
            ],
            vec![[0, 1, 2], [0, 2, 3]],
        )
    }

    fn direct_emitter_scene(camera_above: bool) -> (Scene, AnimatedCamera, ShutterInterval) {
        let white = lift_rgb([1.0, 1.0, 1.0]);
        let (eye_z, forward_z, up_y) = if camera_above {
            (2.0, -1.0, 1.0)
        } else {
            (-2.0, 1.0, 1.0)
        };
        let camera_spec = Camera {
            eye: Point3::new(0.0, 0.0, eye_z),
            forward: Vec3::new(0.0, 0.0, forward_z),
            up: Vec3::new(0.0, up_y, 0.0),
            half_tan: 0.2,
        };
        let physical = PhysicalCamera::try_legacy_compatible(
            camera_spec.eye,
            camera_spec.forward,
            camera_spec.up,
            camera_spec.half_tan,
            2.0,
            Aperture::try_circular(0.0).unwrap(),
        )
        .unwrap();
        let scene = Scene {
            primitives: vec![Primitive {
                shape: Shape::Mesh(emitter_quad(0.0)),
                material: Material::Lambertian { reflectance: white },
                emission: Some((white, 2.0)),
            }],
            lights: vec![RectLight {
                corner: Point3::new(-1.0, -1.0, 0.0),
                edge_u: Vec3::new(2.0, 0.0, 0.0),
                edge_v: Vec3::new(0.0, 2.0, 0.0),
                prim: 0,
                emission: (white, 2.0),
            }],
            environment: None,
            camera: camera_spec,
        };
        let camera = AnimatedCamera::try_static(1, 0.0, 1.0, physical).unwrap();
        let shutter = ShutterInterval::resolve(
            0.5,
            0.0,
            ShutterConvention::Centered,
            ShutterDistribution::StratifiedCounterV1 { strata: 4 },
            ShotTimeBounds::try_new(0.0, 1.0).unwrap(),
        )
        .unwrap();
        (scene, camera, shutter)
    }

    #[test]
    fn g0_density_conversion_obeys_area_jacobian() {
        let source = Vertex::camera(Point3::new(0.0, 0.0, 0.0));
        let mut target = Vertex::camera(Point3::new(0.0, 0.0, 2.0));
        target.geometric_normal = Some(Vec3::new(0.0, 0.0, -1.0));
        assert_eq!(
            convert_density(&source, 4.0, &target).to_bits(),
            1.0_f64.to_bits()
        );
    }

    #[test]
    fn g5_full_path_qmc_is_replayable_and_domain_separated() {
        let settings = Settings {
            width: 1,
            height: 1,
            spp: 16,
            max_depth: 1,
            sampler: Sampler::OwenSobolFullPath,
            strategy: DirectStrategy::Mis,
            seed: 0x514d_4342,
        };
        let sobol = path_sobol(settings.sampler).unwrap();
        let draw = |domain| SubpathRng::new(&settings, 7, 5, domain, Some(&sobol)).next4();
        assert_eq!(draw(CAMERA_WALK_DOMAIN), draw(CAMERA_WALK_DOMAIN));
        assert_ne!(draw(CAMERA_WALK_DOMAIN), draw(LIGHT_WALK_DOMAIN));
        assert_ne!(draw(LIGHT_WALK_DOMAIN), draw(CONNECTION_DOMAIN));
    }

    #[test]
    fn g0_medium_identity_is_not_only_an_ior_comparison() {
        let glass = DielectricGlass::representative_borosilicate();
        let left = MediumEntry {
            boundary_primitive: 3,
            glass,
        };
        let right = MediumEntry {
            boundary_primitive: 4,
            glass,
        };
        assert!(!same_medium(Some(left), Some(right)));
        assert!(same_medium(Some(left), Some(left)));
    }

    #[test]
    fn g0_two_sided_emitter_launches_share_the_emissive_hit_radiance() {
        let emission = (LiftedSpectrum { c: [0.7, 0.6, 0.5] }, 3.0);
        let light = RectLight {
            corner: Point3::new(-0.5, -0.5, 1.0),
            edge_u: Vec3::new(1.0, 0.0, 0.0),
            edge_v: Vec3::new(0.0, 1.0, 0.0),
            prim: 0,
            emission,
        };
        let lights = [light];
        let admitted = AdmittedLighting::try_new(&lights, None).unwrap();
        let front = admitted
            .sample_rectangle_emission(0.3, 0.4, 0.125, 0.7)
            .unwrap();
        let back = admitted
            .sample_rectangle_emission(0.3, 0.4, 0.625, 0.7)
            .unwrap();
        assert!(front.normal.dot(front.direction) > 0.0);
        assert!(back.normal.dot(back.direction) < 0.0);

        let scene = Scene {
            primitives: Vec::new(),
            lights: Vec::new(),
            environment: None,
            camera: Camera {
                eye: Point3::new(0.0, 0.0, 0.0),
                forward: Vec3::new(0.0, 0.0, 1.0),
                up: Vec3::new(0.0, 1.0, 0.0),
                half_tan: 1.0,
            },
        };
        let wavelength_nm = 550.0;
        let front_radiance =
            Vertex::light(front, wavelength_nm).emitted_radiance(&scene, wavelength_nm);
        let back_radiance =
            Vertex::light(back, wavelength_nm).emitted_radiance(&scene, wavelength_nm);
        let hit_radiance = emission.0.eval(wavelength_nm) * emission.1;
        assert_eq!(front_radiance.to_bits(), hit_radiance.to_bits());
        assert_eq!(back_radiance.to_bits(), hit_radiance.to_bits());
    }

    #[test]
    fn g5_splat_publication_is_independent_of_worker_arrival_order() {
        let records = [
            SplatRecord {
                target_pixel: 1,
                sample: 2,
                strategy_s: 3,
                xyz: [1.0e16, 2.0, 3.0],
            },
            SplatRecord {
                target_pixel: 1,
                sample: 0,
                strategy_s: 2,
                xyz: [-1.0e16, 5.0, 7.0],
            },
            SplatRecord {
                target_pixel: 1,
                sample: 1,
                strategy_s: 4,
                xyz: [1.0, 11.0, 13.0],
            },
        ];
        let mut forward = Film::new(2, 1);
        let mut reverse = Film::new(2, 1);
        let mut forward_records = records.to_vec();
        let mut reverse_records = records.into_iter().rev().collect();
        publish_source_splats(&mut forward, &mut forward_records);
        publish_source_splats(&mut reverse, &mut reverse_records);
        assert_eq!(forward.xyz, reverse.xyz);
        assert_eq!(forward.xyz[1][0].to_bits(), 0.0_f64.to_bits());
    }

    #[test]
    fn g0_t1_global_sample_count_closes_two_strategy_mis() {
        let width = 4;
        let height = 2;
        let light_tracing_samples = strategy_sample_multiplicity(1, width, height).unwrap();
        assert_eq!(light_tracing_samples, 8.0);
        assert_eq!(strategy_sample_multiplicity(2, width, height).unwrap(), 1.0);

        // For equal path densities, balance-heuristic weights are proportional
        // to the number of samples drawn by each strategy.  This also proves
        // the two directional ratios close to one rather than merely checking
        // a sum that could hide a strategy permutation.
        let t2_over_t1 = ratio_with_sample_multiplicity(1.0, 2, 1, width, height).unwrap();
        let t1_over_t2 = ratio_with_sample_multiplicity(1.0, 1, 2, width, height).unwrap();
        assert_eq!(t2_over_t1, 1.0 / 8.0);
        assert_eq!(t1_over_t2, 8.0);
        let t1_weight = 1.0 / (1.0 + t2_over_t1);
        let t2_weight = 1.0 / (1.0 + t1_over_t2);
        assert!((t1_weight - 8.0 / 9.0).abs() <= f64::EPSILON);
        assert!((t2_weight - 1.0 / 9.0).abs() <= f64::EPSILON);
        assert!((t1_weight + t2_weight - 1.0).abs() <= f64::EPSILON);
        assert!(strategy_sample_multiplicity(1, 0, height).is_err());
    }

    #[test]
    fn g2_front_and_back_primary_emitter_hits_have_identical_energy() {
        let settings = Settings {
            width: 1,
            height: 1,
            spp: 4,
            max_depth: 0,
            sampler: Sampler::Iid,
            strategy: DirectStrategy::Mis,
            seed: 17,
        };
        with_test_cx(|cx| {
            let (front_scene, front_camera, front_shutter) = direct_emitter_scene(true);
            let front = render_cinematic_bidirectional(
                &front_scene,
                &front_camera,
                CutSide::After,
                cx,
                &settings,
                front_shutter,
            )
            .unwrap();
            let (back_scene, back_camera, back_shutter) = direct_emitter_scene(false);
            let back = render_cinematic_bidirectional(
                &back_scene,
                &back_camera,
                CutSide::After,
                cx,
                &settings,
                back_shutter,
            )
            .unwrap();
            assert_eq!(front.film.xyz, back.film.xyz);
            assert!(front.film.xyz[0].iter().all(|value| *value > 0.0));
            assert!(front.strategies.nonzero > 0);
        });
    }

    #[test]
    fn g5_parallel_sample_blocks_match_serial_and_worker_counts_bitwise() {
        let settings = Settings {
            width: 1,
            height: 1,
            spp: BIDIRECTIONAL_SAMPLE_BLOCK + 1,
            max_depth: 0,
            sampler: Sampler::Iid,
            strategy: DirectStrategy::Mis,
            seed: 29,
        };
        with_test_cx(|cx| {
            let (scene, camera, shutter) = direct_emitter_scene(true);
            let serial = render_cinematic_bidirectional(
                &scene,
                &camera,
                CutSide::After,
                cx,
                &settings,
                shutter,
            )
            .unwrap();
            let one_worker =
                RenderExecutionConfig::try_new(1, 1, 1, 16 * 1024 * 1024, RunId(0x4244_5054))
                    .unwrap();
            let two_workers =
                RenderExecutionConfig::try_new(1, 1, 2, 16 * 1024 * 1024, RunId(0x4244_5055))
                    .unwrap();
            let one = render_cinematic_bidirectional_with_execution(
                &scene,
                &camera,
                CutSide::After,
                cx,
                &settings,
                shutter,
                &one_worker,
            )
            .unwrap();
            let two = render_cinematic_bidirectional_with_execution(
                &scene,
                &camera,
                CutSide::After,
                cx,
                &settings,
                shutter,
                &two_workers,
            )
            .unwrap();
            assert_eq!(serial.film, one.film);
            assert_eq!(one.film, two.film);
            assert_eq!(serial.strategies, one.strategies);
            assert_eq!(one.strategies, two.strategies);
            assert_eq!(one.report.sample_blocks, 2);
            assert_eq!(two.report.sample_blocks, 2);
            assert_eq!(one.report.workers, 1);
            assert_eq!(two.report.workers, 2);
        });
    }
}
