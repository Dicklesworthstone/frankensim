//! Animated Euler-disc composition for the opt-in spectral tracer.
//!
//! The mechanics chart remains the source of truth. Because the current
//! axisymmetric chart intentionally carries no certified sphere-tracing claim,
//! this bridge derives one bounded, deterministic surface-of-revolution mesh
//! from its exact retained line/arc meridian. The mesh is a visualization
//! approximation, built once in center-of-mass coordinates; contact, support,
//! mass, and inertia continue to use the original chart.

use core::fmt;
use std::collections::BTreeSet;

use fs_blake3::{ContentHash, DomainHasher, hash_domain};
use fs_exec::Cx;
use fs_geom::{Aabb, Chart, Point3, Vec3};
use fs_mbd::{MassProperties, Pose as MbdPose, Vec3 as MbdVec3};
use fs_render::animated_instances::{
    AnimatedGeometryInstance, AnimatedInstanceError, RigidTransformTrajectory, TransformKeyframe,
};
use fs_render::camera::{AnimatedCamera, CameraError, CutSide, KeyframeFocus};
use fs_render::charts::TriMesh;
use fs_render::conductor::{ConductorOptics, ConductorSurface};
use fs_render::dielectric::{
    BeerLambertParameters, DielectricGlass, DielectricSurface, GlassProvenance,
};
use fs_render::instances::{GeometryInstance, InstanceError, RigidTransform, SharedGeometry};
use fs_render::motion::{ShotTimeBounds, ShutterConvention, ShutterDistribution, ShutterInterval};
use fs_render::motion_bounds::{
    FiniteLocalAabb, MotionBoundsError, conservative_trajectory_swept_aabb,
};
use fs_render::spectral::lift_rgb;
use fs_render::tracer::{
    AdaptiveRenderOutput, AdaptiveSamplingConfig, Camera, DirectStrategy, Film, Material,
    ParkedRenderScope, PendingAdaptiveRender, PendingRender, Primitive, RectLight,
    RenderExecutionConfig, RenderExecutionError, RenderExecutionOutput, Sampler, Scene, Settings,
    Shape, TracerError, film_to_exr, render, render_cinematic,
    render_cinematic_adaptive_with_execution, render_cinematic_with_execution,
};
use fs_rep_frep::{MeridianPoint, MeridianSegment};

use crate::contact_dynamics::profile_contact_geometry;
use crate::render_motion_bridge::{
    EulerRenderMotionBridge, EulerShutterPartition, RenderMotionBridgeError,
};
use crate::render_trajectory::{
    RenderContactBranch, RenderTrajectoryAuthority, RenderUnitSystem, RenderWorldFrame,
};
use crate::render_trajectory_codec::EulerRenderTrajectoryArtifact;
use crate::specimen::{ResolvedDiscProfile, ResolvedDiscProfileIdentities};
use crate::timeline_resampling::{EventEvaluationSide, ExposureEventPolicy};

/// Version of the scene-composition, asset-binding, and preview-mesh policy.
pub const EULER_RENDER_SCENE_BRIDGE_VERSION: u16 = 1;
/// Domain for the COM-centered preview mesh's complete canonical input.
pub const EULER_PREVIEW_MESH_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-euler-disc-e2e.euler-preview-mesh.v1";
/// Domain for the complete trajectory/configuration scene identity.
pub const EULER_RENDER_SCENE_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-euler-disc-e2e.euler-render-scene.v1";
/// Domain for the complete admitted scene-builder configuration.
pub const EULER_RENDER_CONFIGURATION_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-euler-disc-e2e.euler-render-configuration.v1";

/// Hard topology ceilings for the one-time render-only lathe conversion.
pub const MAX_EULER_PREVIEW_VERTICES: usize = 1_000_000;
/// Hard topology ceiling paired with [`MAX_EULER_PREVIEW_VERTICES`].
pub const MAX_EULER_PREVIEW_TRIANGLES: usize = 2_000_000;
/// Maximum fixed azimuthal samples accepted by v1.
pub const MAX_EULER_AZIMUTHAL_SEGMENTS: u32 = 4_096;
/// Maximum chord subdivisions accepted for one circular meridian arc.
pub const MAX_EULER_ARC_SUBDIVISIONS: u32 = 1_024;

const CONTACT_BASE_ALIGNMENT_TOLERANCE_M: f64 = 1.0e-10;
const CONTACT_NORMAL_ALIGNMENT_TOLERANCE: f64 = 1.0e-10;

/// Stable default object identity for the animated disc.
pub const EULER_DISC_OBJECT_ID: u64 = 0x4555_4c45_525f_0001;
/// Stable default object identity for the moving plate.
pub const EULER_BASE_PLATE_OBJECT_ID: u64 = 0x4555_4c45_525f_0002;
/// Stable default object identity for the static housing.
pub const EULER_HOUSING_OBJECT_ID: u64 = 0x4555_4c45_525f_0003;
/// Stable default object identity for an optional diagnostic marker.
pub const EULER_DEBUG_MARKER_OBJECT_ID: u64 = 0x4555_4c45_525f_00ff;

/// Binding no-claim for the current scene bridge.
pub const EULER_RENDER_SCENE_NO_CLAIM: &str = "scene composition and chordal preview geometry do not validate Euler mechanics, measured material parameters, calibrated lighting, or a real apparatus";

/// Declared length convention of scene configuration values.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum EulerSceneLengthUnit {
    /// SI metres, the only v1 admission.
    Metres,
    /// Deliberately unsupported configuration used to catch unit cross-wiring.
    Millimetres,
}

/// Fixed deterministic tessellation controls for the derived preview mesh.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct EulerTessellationConfig {
    /// Samples around a complete revolution; must be at least eight.
    pub azimuthal_segments: u32,
    /// Equal-angle chords used for each true circular meridian arc.
    pub arc_subdivisions_per_arc: u32,
}

impl EulerTessellationConfig {
    /// Bounded interactive-quality default; 4K qualification may select a
    /// finer admitted configuration and receives its own identity.
    pub const DEFAULT: Self = Self {
        azimuthal_segments: 256,
        arc_subdivisions_per_arc: 32,
    };
}

impl Default for EulerTessellationConfig {
    fn default() -> Self {
        Self::DEFAULT
    }
}

/// Explicit material appearance for one Euler-scene object.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum EulerMaterialStyle {
    /// Ideal diffuse reflection with bounded linear-RGB reflectance.
    Lambertian {
        /// Linear-RGB reflectance, each component in `[0, 1]`.
        linear_rgb: [f64; 3],
    },
    /// Opaque isotropic GGX reflection.
    Ggx {
        /// Linear-RGB reflectance, each component in `[0, 1]`.
        linear_rgb: [f64; 3],
        /// Positive isotropic microfacet roughness.
        alpha: f64,
    },
    /// Opaque spectral metal with exact complex-IOR Fresnel and an isotropic
    /// single-scattering GGX surface.
    Conductor {
        /// Validated complex-index table and its source provenance.
        optics: ConductorOptics,
        /// Validated isotropic GGX roughness.
        surface: ConductorSurface,
    },
    /// Homogeneous spectral dielectric with smooth or rough transmission.
    Dielectric {
        /// Validated interior optical medium and source provenance.
        glass: DielectricGlass,
        /// Smooth-delta or isotropic-GGX boundary treatment.
        surface: DielectricSurface,
    },
}

/// Representative base and housing dimensions in metres.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EulerBaseVisualSpec {
    /// Plate width along local x.
    pub plate_width_m: f64,
    /// Plate depth along local y.
    pub plate_depth_m: f64,
    /// Plate thickness below its contact plane.
    pub plate_thickness_m: f64,
    /// Housing width along local x.
    pub housing_width_m: f64,
    /// Housing depth along local y.
    pub housing_depth_m: f64,
    /// Housing height.
    pub housing_height_m: f64,
    /// Nominal vertical clearance from plate bottom to housing top.
    pub housing_gap_m: f64,
}

/// Stable scene instance identities.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EulerSceneObjectIds {
    /// Animated disc object ID.
    pub disc: u64,
    /// Animated base-plate object ID.
    pub base_plate: u64,
    /// Static housing object ID.
    pub housing: u64,
    /// Optional debug marker object ID.
    pub debug_marker: u64,
}

impl EulerSceneObjectIds {
    /// Reference stable IDs.
    pub const DEFAULT: Self = Self {
        disc: EULER_DISC_OBJECT_ID,
        base_plate: EULER_BASE_PLATE_OBJECT_ID,
        housing: EULER_HOUSING_OBJECT_ID,
        debug_marker: EULER_DEBUG_MARKER_OBJECT_ID,
    };
}

impl Default for EulerSceneObjectIds {
    fn default() -> Self {
        Self::DEFAULT
    }
}

/// Single static rectangular softbox admitted by tracer v1.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EulerRectLightSpec {
    /// World-space corner (m).
    pub corner_world_m: Point3,
    /// First world-space edge (m).
    pub edge_u_world_m: Vec3,
    /// Second world-space edge (m).
    pub edge_v_world_m: Vec3,
    /// Linear RGB emitter spectrum input.
    pub linear_rgb: [f64; 3],
    /// Positive radiance scale.
    pub radiance_scale: f64,
}

/// Optional diagnostic geometry kept outside the beauty-object set.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum EulerDebugOverlay {
    /// Beauty scene only.
    None,
    /// A static octahedral marker at one retained closed-contact sample.
    ContactMarker {
        /// Source sample index.
        sample_index: usize,
        /// Marker radius (m).
        radius_m: f64,
    },
}

/// Complete configuration consumed by the v1 scene builder.
#[derive(Clone, Debug, PartialEq)]
pub struct EulerSceneConfig {
    /// Explicit configuration length unit.
    pub length_unit: EulerSceneLengthUnit,
    /// Preview conversion resolution.
    pub tessellation: EulerTessellationConfig,
    /// Representative product-base dimensions.
    pub base: EulerBaseVisualSpec,
    /// Stable object IDs.
    pub object_ids: EulerSceneObjectIds,
    /// Disc appearance.
    pub disc_material: EulerMaterialStyle,
    /// Spectral base-plate appearance. The reference configuration uses a
    /// representative, explicitly non-measured glass preset.
    pub plate_material: EulerMaterialStyle,
    /// Housing appearance.
    pub housing_material: EulerMaterialStyle,
    /// One static softbox.
    pub light: EulerRectLightSpec,
    /// Validated animated physical camera.
    pub camera: AnimatedCamera,
    /// Positive declared near-depth requirement (m). The ray tracer itself has
    /// no clipping plane; this is an admission check over subject bounds.
    pub camera_near_m: f64,
    /// Declared far-depth requirement (m), larger than `camera_near_m`.
    pub camera_far_m: f64,
    /// Largest admitted endpoint angular-speed times sample duration (rad).
    /// This prevents an obviously aliased shortest-arc interpolation.
    pub maximum_angular_step_rad: f64,
    /// Optional isolated diagnostic layer.
    pub debug_overlay: EulerDebugOverlay,
}

impl EulerSceneConfig {
    /// Reference geometry and spectral-material setup around a supplied
    /// admitted camera. Final conductor, engraving, and studio look development
    /// remain separate downstream capabilities.
    #[must_use]
    pub fn reference(camera: AnimatedCamera) -> Self {
        Self {
            length_unit: EulerSceneLengthUnit::Metres,
            tessellation: EulerTessellationConfig::DEFAULT,
            base: EulerBaseVisualSpec {
                plate_width_m: 0.18,
                plate_depth_m: 0.18,
                plate_thickness_m: 0.010,
                housing_width_m: 0.20,
                housing_depth_m: 0.20,
                housing_height_m: 0.025,
                housing_gap_m: 0.002,
            },
            object_ids: EulerSceneObjectIds::DEFAULT,
            disc_material: EulerMaterialStyle::Ggx {
                linear_rgb: [0.72, 0.74, 0.78],
                alpha: 0.12,
            },
            plate_material: EulerMaterialStyle::Dielectric {
                glass: DielectricGlass::representative_crown(),
                surface: DielectricSurface::POLISHED,
            },
            housing_material: EulerMaterialStyle::Ggx {
                linear_rgb: [0.055, 0.065, 0.08],
                alpha: 0.28,
            },
            light: EulerRectLightSpec {
                corner_world_m: Point3::new(-0.08, 0.05, 0.16),
                edge_u_world_m: Vec3::new(0.16, 0.0, 0.0),
                edge_v_world_m: Vec3::new(0.0, -0.10, 0.0),
                linear_rgb: [1.0, 0.96, 0.90],
                radiance_scale: 24.0,
            },
            camera,
            camera_near_m: 0.01,
            camera_far_m: 2.0,
            maximum_angular_step_rad: core::f64::consts::FRAC_PI_2,
            debug_overlay: EulerDebugOverlay::None,
        }
    }
}

/// Audit receipt for the render-only axisymmetric conversion.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EulerPreviewMeshReceipt {
    /// Exact source chart identity.
    pub source_chart_identity: ContentHash,
    /// Canonical COM-centered triangle geometry identity.
    pub mesh_identity: ContentHash,
    /// Active azimuthal resolution.
    pub azimuthal_segments: u32,
    /// Active per-arc meridian resolution.
    pub arc_subdivisions_per_arc: u32,
    /// Published vertex count.
    pub vertex_count: usize,
    /// Published triangle count.
    pub triangle_count: usize,
    /// Maximum circular-meridian chord sagitta (m).
    pub maximum_meridian_chord_error_m: f64,
    /// Maximum azimuthal chord sagitta at the source outer radius (m).
    pub maximum_azimuthal_chord_error_m: f64,
    /// COM-centered local mesh bounds.
    pub local_bounds_m: Aabb,
    /// Deterministic BVH-layout diagnostic, not a content address.
    pub bvh_fingerprint: u64,
}

/// Stable locations of semantic primitives in the emitted scene.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EulerScenePrimitiveIndices {
    /// Animated disc primitive.
    pub disc: usize,
    /// Animated plate primitive.
    pub base_plate: usize,
    /// Static housing primitive.
    pub housing: usize,
    /// Emissive softbox primitive referenced by the sole v1 `Scene::lights` entry.
    pub light: usize,
}

/// Identity and source binding of the optional non-beauty diagnostic layer.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EulerDebugLayerReceipt {
    /// Identity bound to the beauty scene, marker geometry, and placement.
    pub layer_identity: ContentHash,
    /// Stable marker object identity.
    pub object_id: u64,
    /// Retained closed-contact sample selected by the configuration.
    pub source_sample_index: usize,
    /// Declared octahedral marker radius (m).
    pub radius_m: f64,
}

/// Exact scene transforms reconstructed for one admitted timeline time.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EulerScenePose {
    /// COM-centered disc-mesh placement.
    pub disc: RigidTransform,
    /// Plate-local placement including reduced-base displacement.
    pub base_plate: RigidTransform,
}

/// One explicit event-delimited shutter segment.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EulerPreparedFrameSegment {
    scene_identity: ContentHash,
    shutter: ShutterInterval,
    duration_weight: f64,
}

impl EulerPreparedFrameSegment {
    /// Segment-local shutter passed to the tracer.
    #[must_use]
    pub const fn shutter(&self) -> ShutterInterval {
        self.shutter
    }

    /// Segment duration divided by the complete positive exposure.
    #[must_use]
    pub const fn duration_weight(&self) -> f64 {
        self.duration_weight
    }
}

/// Event-aware prepared frame. Multi-segment films remain separate so the
/// bridge never fabricates tracer sample-count provenance while combining them.
#[derive(Clone, Debug, PartialEq)]
pub struct EulerPreparedFrame {
    scene_identity: ContentHash,
    cut_side: CutSide,
    segments: Vec<EulerPreparedFrameSegment>,
}

impl EulerPreparedFrame {
    /// Beauty scene that resolved this event-aware frame.
    #[must_use]
    pub const fn scene_identity(&self) -> ContentHash {
        self.scene_identity
    }

    /// Camera-shot ownership used for an exact cut-boundary exposure.
    #[must_use]
    pub const fn cut_side(&self) -> CutSide {
        self.cut_side
    }

    /// Ordered explicit shutter segments.
    #[must_use]
    pub fn segments(&self) -> &[EulerPreparedFrameSegment] {
        &self.segments
    }
}

/// Request for one frame exposure inside the accepted trajectory horizon.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EulerFrameRequest {
    /// Reference frame time (s).
    pub frame_time_s: f64,
    /// Nonnegative exposure duration (s).
    pub exposure_duration_s: f64,
    /// Placement of the shutter around `frame_time_s`.
    pub convention: ShutterConvention,
    /// Deterministic time-sampling distribution.
    pub distribution: ShutterDistribution,
    /// Contact/event crossing policy.
    pub event_policy: ExposureEventPolicy,
    /// Camera-shot ownership at an exact cut.
    pub cut_side: CutSide,
}

/// Complete owned render scene bound to one decoded trajectory artifact.
pub struct EulerCinematicScene<'artifact> {
    artifact: &'artifact EulerRenderTrajectoryArtifact,
    scene: Scene,
    camera: AnimatedCamera,
    source_configuration_identity: ContentHash,
    scene_identity: ContentHash,
    preview_mesh: EulerPreviewMeshReceipt,
    primitive_indices: EulerScenePrimitiveIndices,
    subject_bounds_m: Aabb,
    debug_layer: Option<EulerDebugLayer>,
}

struct EulerDebugLayer {
    receipt: EulerDebugLayerReceipt,
    instance: GeometryInstance,
    material: Material,
}

impl<'artifact> EulerCinematicScene<'artifact> {
    /// Build the complete reference scene from one admitted artifact and the
    /// exact resolved specimen asset named by its metadata.
    pub fn try_build(
        artifact: &'artifact EulerRenderTrajectoryArtifact,
        specimen: &ResolvedDiscProfile,
        config: EulerSceneConfig,
        cx: &Cx<'_>,
    ) -> Result<Self, EulerSceneError> {
        checkpoint(cx)?;
        validate_config(&config)?;
        if !artifact.declared_discontinuities().is_empty() {
            return Err(EulerSceneError::DeclaredDiscontinuityUnsupported);
        }
        let identities = validate_asset_binding(artifact, specimen)?;
        validate_angular_sampling(artifact, config.maximum_angular_step_rad)?;
        validate_contact_base_alignment(artifact, specimen, cx)?;
        validate_camera_coverage(artifact, &config.camera, cx)?;

        let (disc_mesh, preview_mesh) =
            tessellate_disc(specimen, identities, config.tessellation, cx)?;
        let disc_geometry = SharedGeometry::mesh(disc_mesh);
        let disc_trajectory = disc_transform_trajectory(artifact, cx)?;
        let disc = AnimatedGeometryInstance::try_new(
            config.object_ids.disc,
            preview_mesh.mesh_identity,
            disc_geometry,
            disc_trajectory,
        )?;

        let plate_mesh = box_mesh(
            0.5 * config.base.plate_width_m,
            0.5 * config.base.plate_depth_m,
            -config.base.plate_thickness_m,
            0.0,
        );
        let plate_geometry_identity = box_identity(
            b"animated-base-plate",
            [
                config.base.plate_width_m,
                config.base.plate_depth_m,
                config.base.plate_thickness_m,
                0.0,
            ],
        );
        let plate = AnimatedGeometryInstance::try_new(
            config.object_ids.base_plate,
            plate_geometry_identity,
            SharedGeometry::mesh(plate_mesh),
            base_transform_trajectory(artifact, cx)?,
        )?;

        let housing_top = -config.base.plate_thickness_m - config.base.housing_gap_m;
        let housing_bottom = housing_top - config.base.housing_height_m;
        let housing_geometry_identity = box_identity(
            b"static-housing",
            [
                config.base.housing_width_m,
                config.base.housing_depth_m,
                housing_bottom,
                housing_top,
            ],
        );
        let housing = GeometryInstance::try_new(
            config.object_ids.housing,
            housing_geometry_identity,
            SharedGeometry::mesh(box_mesh(
                0.5 * config.base.housing_width_m,
                0.5 * config.base.housing_depth_m,
                housing_bottom,
                housing_top,
            )),
            nominal_base_transform(artifact)?,
        )?;

        let subject_bounds_m =
            subject_bounds(artifact, preview_mesh, &config.base, &disc, &plate, cx)?;
        validate_camera_depths(
            artifact,
            &config.camera,
            subject_bounds_m,
            config.camera_near_m,
            config.camera_far_m,
            cx,
        )?;
        let source_configuration_identity = configuration_identity(&config);
        let scene_identity = scene_identity(
            artifact,
            identities,
            preview_mesh,
            plate_geometry_identity,
            housing_geometry_identity,
            &config,
        );
        let debug_layer = build_debug_layer(
            artifact,
            config.debug_overlay,
            config.object_ids.debug_marker,
            scene_identity,
        )?;

        let disc_material = material(config.disc_material);
        let plate_material = material(config.plate_material);
        let housing_material = material(config.housing_material);
        let light_spectrum = lift_rgb(config.light.linear_rgb);
        let light_emission = (light_spectrum, config.light.radiance_scale);

        let mut primitives = vec![
            Primitive {
                shape: Shape::AnimatedInstance(disc),
                material: disc_material,
                emission: None,
            },
            Primitive {
                shape: Shape::AnimatedInstance(plate),
                material: plate_material,
                emission: None,
            },
            Primitive {
                shape: Shape::Instance(housing),
                material: housing_material,
                emission: None,
            },
        ];
        let light_index = primitives.len();
        primitives.push(Primitive {
            shape: Shape::Mesh(rect_mesh(config.light)),
            material: Material::Lambertian {
                reflectance: light_spectrum,
            },
            emission: Some(light_emission),
        });

        checkpoint(cx)?;

        let first_time = artifact.trajectory().samples()[0].input().time_s;
        let first_camera = config.camera.evaluate(cx, first_time, CutSide::After)?;
        let legacy_camera = legacy_camera(&first_camera);
        let primitive_indices = EulerScenePrimitiveIndices {
            disc: 0,
            base_plate: 1,
            housing: 2,
            light: light_index,
        };
        checkpoint(cx)?;

        Ok(Self {
            artifact,
            scene: Scene {
                primitives,
                lights: vec![RectLight {
                    corner: config.light.corner_world_m,
                    edge_u: config.light.edge_u_world_m,
                    edge_v: config.light.edge_v_world_m,
                    prim: light_index,
                    emission: light_emission,
                }],
                environment: None,
                camera: legacy_camera,
            },
            camera: config.camera,
            source_configuration_identity,
            scene_identity,
            preview_mesh,
            primitive_indices,
            subject_bounds_m,
            debug_layer,
        })
    }

    /// Complete deterministic scene identity.
    #[must_use]
    pub const fn scene_identity(&self) -> ContentHash {
        self.scene_identity
    }

    /// Complete identity of the admitted scene-builder configuration.
    ///
    /// Unlike [`Self::scene_identity`], this does not include the source
    /// trajectory or resolved specimen. It is the configuration component of
    /// durable renderer-checkpoint provenance.
    #[must_use]
    pub const fn source_configuration_identity(&self) -> ContentHash {
        self.source_configuration_identity
    }

    /// Source trajectory artifact identity.
    #[must_use]
    pub fn source_trajectory_identity(&self) -> ContentHash {
        self.artifact.receipt().artifact_identity()
    }

    /// Unchanged scientific authority ceiling.
    #[must_use]
    pub fn source_authority(&self) -> RenderTrajectoryAuthority {
        self.artifact.trajectory().metadata().authority
    }

    /// Derived preview-mesh audit receipt.
    #[must_use]
    pub const fn preview_mesh_receipt(&self) -> EulerPreviewMeshReceipt {
        self.preview_mesh
    }

    /// Stable semantic primitive indices.
    #[must_use]
    pub const fn primitive_indices(&self) -> EulerScenePrimitiveIndices {
        self.primitive_indices
    }

    /// Optional diagnostic-layer receipt. The marker is never present in the
    /// beauty scene returned by [`Self::scene`] or any frame-render API.
    #[must_use]
    pub fn debug_layer_receipt(&self) -> Option<EulerDebugLayerReceipt> {
        self.debug_layer.as_ref().map(|layer| layer.receipt)
    }

    /// Conservative subject bounds over accepted source samples.
    #[must_use]
    pub const fn subject_bounds_m(&self) -> Aabb {
        self.subject_bounds_m
    }

    /// Underlying beauty-only tracer scene for orchestration and inspection.
    #[must_use]
    pub const fn scene(&self) -> &Scene {
        &self.scene
    }

    /// Validated physical/keyframed camera.
    #[must_use]
    pub const fn camera(&self) -> &AnimatedCamera {
        &self.camera
    }

    /// Reconstruct disc and plate transforms through the trajectory's public
    /// resampler rather than treating renderer interpolation as solver state.
    pub fn pose_at(
        &self,
        time_s: f64,
        side: EventEvaluationSide,
    ) -> Result<EulerScenePose, EulerSceneError> {
        let motion = EulerRenderMotionBridge::new(self.artifact.trajectory());
        let sample = motion.sample_at_time(time_s, side)?;
        let base_plate = base_transform(
            self.artifact,
            sample.timeline_sample().base_mode.displacement_m,
        )?;
        Ok(EulerScenePose {
            disc: sample.transform(),
            base_plate,
        })
    }

    /// Resolve event-delimited shutter segments before any ray is evaluated.
    pub fn prepare_frame(
        &self,
        request: EulerFrameRequest,
    ) -> Result<EulerPreparedFrame, EulerSceneError> {
        let motion = EulerRenderMotionBridge::new(self.artifact.trajectory());
        let prepared = motion.resolve_shutter(
            request.frame_time_s,
            request.exposure_duration_s,
            request.convention,
            request.distribution,
            request.event_policy,
        )?;
        let segment_count = match prepared.partition() {
            EulerShutterPartition::Static { .. } => 1,
            EulerShutterPartition::EventDelimited(partition) => partition.segments.len(),
        };
        let mut segments = Vec::new();
        segments
            .try_reserve_exact(segment_count)
            .map_err(|_| EulerSceneError::Capacity("prepared shutter segments"))?;
        for index in 0..segment_count {
            let segment = motion.shutter_segment(&prepared, index)?;
            segments.push(EulerPreparedFrameSegment {
                scene_identity: self.scene_identity,
                shutter: segment.shutter(),
                duration_weight: segment.duration_weight(),
            });
        }
        Ok(EulerPreparedFrame {
            scene_identity: self.scene_identity,
            cut_side: request.cut_side,
            segments,
        })
    }

    /// Render one explicit prepared segment. Multi-segment weighted image
    /// composition remains visible to orchestration instead of falsifying the
    /// tracer film's sample-count/time-mode provenance.
    pub fn render_segment(
        &self,
        prepared: &EulerPreparedFrame,
        segment_index: usize,
        settings: &Settings,
        cx: &Cx<'_>,
    ) -> Result<Film, EulerSceneError> {
        let segment = self.prepared_segment(prepared, segment_index)?;
        Ok(render_cinematic(
            &self.scene,
            &self.camera,
            prepared.cut_side,
            cx,
            settings,
            segment.shutter,
        )?)
    }

    /// Tile-parallel render of one explicit prepared segment under a caller-
    /// supplied, replayable execution policy. The returned output retains the
    /// executor, memory-admission, and tile-layout report.
    pub fn render_segment_with_execution(
        &self,
        prepared: &EulerPreparedFrame,
        segment_index: usize,
        settings: &Settings,
        execution: &RenderExecutionConfig,
        cx: &Cx<'_>,
    ) -> Result<RenderExecutionOutput, EulerSceneError> {
        let segment = self.prepared_segment(prepared, segment_index)?;
        Ok(render_cinematic_with_execution(
            &self.scene,
            &self.camera,
            prepared.cut_side,
            cx,
            settings,
            segment.shutter,
            execution,
        )?)
    }

    /// Render one explicit prepared segment using an already parked worker
    /// crew. Reusing the scope across segments or frames avoids creating and
    /// joining a fresh crew for every animation job while retaining each
    /// job's explicit execution policy and run identity.
    pub fn render_segment_with_parked_scope(
        &self,
        parked: &ParkedRenderScope<'_>,
        prepared: &EulerPreparedFrame,
        segment_index: usize,
        settings: &Settings,
        execution: &RenderExecutionConfig,
        cx: &Cx<'_>,
    ) -> Result<RenderExecutionOutput, EulerSceneError> {
        let segment = self.prepared_segment(prepared, segment_index)?;
        Ok(parked.render_cinematic(
            &self.scene,
            &self.camera,
            prepared.cut_side,
            cx,
            settings,
            segment.shutter,
            execution,
        )?)
    }

    /// Deterministic adaptive render of one explicit prepared segment. The
    /// returned film retains raw sums, estimator moments, sample counts, and
    /// terminal decisions; it does not convert the error proxy into a physical
    /// or perceptual quality claim.
    pub fn render_segment_adaptive_with_execution(
        &self,
        prepared: &EulerPreparedFrame,
        segment_index: usize,
        settings: &Settings,
        adaptive: AdaptiveSamplingConfig,
        execution: &RenderExecutionConfig,
        cx: &Cx<'_>,
    ) -> Result<AdaptiveRenderOutput, EulerSceneError> {
        let segment = self.prepared_segment(prepared, segment_index)?;
        Ok(render_cinematic_adaptive_with_execution(
            &self.scene,
            &self.camera,
            prepared.cut_side,
            cx,
            settings,
            adaptive,
            segment.shutter,
            execution,
        )?)
    }

    /// Adaptive render of one explicit prepared segment on an already parked
    /// animation crew.
    pub fn render_segment_adaptive_with_parked_scope(
        &self,
        parked: &ParkedRenderScope<'_>,
        prepared: &EulerPreparedFrame,
        segment_index: usize,
        settings: &Settings,
        adaptive: AdaptiveSamplingConfig,
        execution: &RenderExecutionConfig,
        cx: &Cx<'_>,
    ) -> Result<AdaptiveRenderOutput, EulerSceneError> {
        let segment = self.prepared_segment(prepared, segment_index)?;
        Ok(parked.render_cinematic_adaptive(
            &self.scene,
            &self.camera,
            prepared.cut_side,
            cx,
            settings,
            adaptive,
            segment.shutter,
            execution,
        )?)
    }

    /// Begin an opaque single-film segment job whose committed row prefixes
    /// survive cancellation without exposing a partial film. Resume it with
    /// `PendingRender::resume_on_parked` inside an animation crew scope, or
    /// with `PendingRender::resume` on a one-shot lane.
    pub fn begin_segment_render(
        &self,
        prepared: &EulerPreparedFrame,
        segment_index: usize,
        settings: Settings,
        execution: RenderExecutionConfig,
        cx: &Cx<'_>,
    ) -> Result<PendingRender<'_>, EulerSceneError> {
        let segment = self.prepared_segment(prepared, segment_index)?;
        Ok(PendingRender::begin_cinematic(
            &self.scene,
            &self.camera,
            prepared.cut_side,
            cx,
            settings,
            segment.shutter,
            execution,
        )?)
    }

    /// Begin an opaque adaptive segment render whose complete-row prefixes and
    /// statistical AOVs survive in-process cancellation.
    pub fn begin_segment_adaptive_render(
        &self,
        prepared: &EulerPreparedFrame,
        segment_index: usize,
        settings: Settings,
        adaptive: AdaptiveSamplingConfig,
        execution: RenderExecutionConfig,
        cx: &Cx<'_>,
    ) -> Result<PendingAdaptiveRender<'_>, EulerSceneError> {
        let segment = self.prepared_segment(prepared, segment_index)?;
        Ok(PendingAdaptiveRender::begin_cinematic(
            &self.scene,
            &self.camera,
            prepared.cut_side,
            cx,
            settings,
            adaptive,
            segment.shutter,
            execution,
        )?)
    }

    /// Convenience render for a frame requiring no event subdivision.
    pub fn render_frame(
        &self,
        request: EulerFrameRequest,
        settings: &Settings,
        cx: &Cx<'_>,
    ) -> Result<Film, EulerSceneError> {
        let prepared = self.prepare_frame(request)?;
        if prepared.segments.len() != 1 {
            return Err(EulerSceneError::ExposureNeedsComposition {
                segment_count: prepared.segments.len(),
            });
        }
        self.render_segment(&prepared, 0, settings, cx)
    }

    /// Tile-parallel convenience render for a frame requiring no event
    /// subdivision. Event-delimited multi-film composition remains explicit.
    pub fn render_frame_with_execution(
        &self,
        request: EulerFrameRequest,
        settings: &Settings,
        execution: &RenderExecutionConfig,
        cx: &Cx<'_>,
    ) -> Result<RenderExecutionOutput, EulerSceneError> {
        let prepared = self.prepare_frame(request)?;
        if prepared.segments.len() != 1 {
            return Err(EulerSceneError::ExposureNeedsComposition {
                segment_count: prepared.segments.len(),
            });
        }
        self.render_segment_with_execution(&prepared, 0, settings, execution, cx)
    }

    /// Parked-crew convenience render for a frame requiring no event
    /// subdivision. Call this repeatedly inside
    /// [`fs_render::tracer::RenderWorkerPool::with_parked_crew_local`] to reuse
    /// one worker crew across an animation batch. Event-delimited multi-film
    /// composition remains explicit.
    pub fn render_frame_with_parked_scope(
        &self,
        parked: &ParkedRenderScope<'_>,
        request: EulerFrameRequest,
        settings: &Settings,
        execution: &RenderExecutionConfig,
        cx: &Cx<'_>,
    ) -> Result<RenderExecutionOutput, EulerSceneError> {
        let prepared = self.prepare_frame(request)?;
        if prepared.segments.len() != 1 {
            return Err(EulerSceneError::ExposureNeedsComposition {
                segment_count: prepared.segments.len(),
            });
        }
        self.render_segment_with_parked_scope(parked, &prepared, 0, settings, execution, cx)
    }

    /// Deterministic adaptive convenience render for a frame requiring no
    /// event subdivision.
    pub fn render_frame_adaptive_with_execution(
        &self,
        request: EulerFrameRequest,
        settings: &Settings,
        adaptive: AdaptiveSamplingConfig,
        execution: &RenderExecutionConfig,
        cx: &Cx<'_>,
    ) -> Result<AdaptiveRenderOutput, EulerSceneError> {
        let prepared = self.prepare_frame(request)?;
        if prepared.segments.len() != 1 {
            return Err(EulerSceneError::ExposureNeedsComposition {
                segment_count: prepared.segments.len(),
            });
        }
        self.render_segment_adaptive_with_execution(&prepared, 0, settings, adaptive, execution, cx)
    }

    /// Parked-crew adaptive convenience render for a frame requiring no event
    /// subdivision.
    pub fn render_frame_adaptive_with_parked_scope(
        &self,
        parked: &ParkedRenderScope<'_>,
        request: EulerFrameRequest,
        settings: &Settings,
        adaptive: AdaptiveSamplingConfig,
        execution: &RenderExecutionConfig,
        cx: &Cx<'_>,
    ) -> Result<AdaptiveRenderOutput, EulerSceneError> {
        let prepared = self.prepare_frame(request)?;
        if prepared.segments.len() != 1 {
            return Err(EulerSceneError::ExposureNeedsComposition {
                segment_count: prepared.segments.len(),
            });
        }
        self.render_segment_adaptive_with_parked_scope(
            parked, &prepared, 0, settings, adaptive, execution, cx,
        )
    }

    /// Begin an opaque resumable job for a frame requiring no event
    /// subdivision. Event-delimited multi-film composition remains explicit.
    pub fn begin_frame_render(
        &self,
        request: EulerFrameRequest,
        settings: Settings,
        execution: RenderExecutionConfig,
        cx: &Cx<'_>,
    ) -> Result<PendingRender<'_>, EulerSceneError> {
        let prepared = self.prepare_frame(request)?;
        if prepared.segments.len() != 1 {
            return Err(EulerSceneError::ExposureNeedsComposition {
                segment_count: prepared.segments.len(),
            });
        }
        self.begin_segment_render(&prepared, 0, settings, execution, cx)
    }

    /// Begin an opaque adaptive job for a frame requiring no event
    /// subdivision.
    pub fn begin_frame_adaptive_render(
        &self,
        request: EulerFrameRequest,
        settings: Settings,
        adaptive: AdaptiveSamplingConfig,
        execution: RenderExecutionConfig,
        cx: &Cx<'_>,
    ) -> Result<PendingAdaptiveRender<'_>, EulerSceneError> {
        let prepared = self.prepare_frame(request)?;
        if prepared.segments.len() != 1 {
            return Err(EulerSceneError::ExposureNeedsComposition {
                segment_count: prepared.segments.len(),
            });
        }
        self.begin_segment_adaptive_render(&prepared, 0, settings, adaptive, execution, cx)
    }

    fn prepared_segment<'prepared>(
        &self,
        prepared: &'prepared EulerPreparedFrame,
        segment_index: usize,
    ) -> Result<&'prepared EulerPreparedFrameSegment, EulerSceneError> {
        if prepared.scene_identity != self.scene_identity {
            return Err(EulerSceneError::PreparedFrameMismatch);
        }
        let segment = prepared.segments.get(segment_index).ok_or(
            EulerSceneError::InvalidPreparedSegment {
                index: segment_index,
                segment_count: prepared.segments.len(),
            },
        )?;
        if segment.scene_identity != self.scene_identity {
            return Err(EulerSceneError::PreparedFrameMismatch);
        }
        Ok(segment)
    }

    /// Revalidate a prepared segment at the L6 sharding boundary without
    /// exposing the bridge's scene-binding fields. A decoded render plan is
    /// never sufficient authority to manufacture a shutter: the coordinator
    /// must present the original scene-bound prepared frame again.
    pub(crate) fn prepared_segment_shard_binding(
        &self,
        prepared: &EulerPreparedFrame,
        segment_index: usize,
    ) -> Result<(ShutterInterval, CutSide), EulerSceneError> {
        let segment = self.prepared_segment(prepared, segment_index)?;
        Ok((segment.shutter, prepared.cut_side))
    }

    /// Render and encode one non-subdivided frame as linear floating-point EXR.
    pub fn render_frame_exr(
        &self,
        request: EulerFrameRequest,
        settings: &Settings,
        cx: &Cx<'_>,
    ) -> Result<Vec<u8>, EulerSceneError> {
        let film = self.render_frame(request, settings, cx)?;
        film_to_exr(&film).map_err(|_| EulerSceneError::ImageEncoding)
    }

    /// Materialize the same composed scene at one exact renderer time. This is
    /// a focused zero-width/reference seam, not a replacement for motion blur.
    pub fn static_scene_at(
        &self,
        time_s: f64,
        cut_side: CutSide,
        cx: &Cx<'_>,
    ) -> Result<Scene, EulerSceneError> {
        checkpoint(cx)?;
        let mut primitives = Vec::new();
        primitives
            .try_reserve_exact(self.scene.primitives.len())
            .map_err(|_| EulerSceneError::Capacity("static scene primitives"))?;
        for primitive in &self.scene.primitives {
            checkpoint(cx)?;
            let shape = match &primitive.shape {
                Shape::AnimatedInstance(instance) => {
                    Shape::Instance(instance.instance_at(cx, time_s)?)
                }
                Shape::Instance(instance) => Shape::Instance(instance.clone()),
                Shape::Mesh(mesh) => Shape::Mesh(mesh.clone()),
                Shape::Chart(_) => return Err(EulerSceneError::UnexpectedSceneShape),
            };
            primitives.push(Primitive {
                shape,
                material: primitive.material,
                emission: primitive.emission,
            });
        }
        let physical = self.camera.evaluate(cx, time_s, cut_side)?;
        Ok(Scene {
            primitives,
            lights: vec![RectLight {
                corner: self.scene.lights[0].corner,
                edge_u: self.scene.lights[0].edge_u,
                edge_v: self.scene.lights[0].edge_v,
                prim: self.scene.lights[0].prim,
                emission: self.scene.lights[0].emission,
            }],
            environment: None,
            camera: legacy_camera(&physical),
        })
    }

    /// Materialize a static diagnostic scene by explicitly adding the optional
    /// marker after every beauty primitive. This is the only API that includes
    /// configured debug geometry.
    pub fn static_scene_with_debug_at(
        &self,
        time_s: f64,
        cut_side: CutSide,
        cx: &Cx<'_>,
    ) -> Result<Scene, EulerSceneError> {
        let mut scene = self.static_scene_at(time_s, cut_side, cx)?;
        if let Some(layer) = &self.debug_layer {
            scene
                .primitives
                .try_reserve_exact(1)
                .map_err(|_| EulerSceneError::Capacity("debug scene primitive"))?;
            scene.primitives.push(layer.primitive());
        }
        Ok(scene)
    }

    /// Render the exact static scene seam without assigning shutter provenance.
    pub fn render_static_at(
        &self,
        time_s: f64,
        cut_side: CutSide,
        settings: &Settings,
        cx: &Cx<'_>,
    ) -> Result<Film, EulerSceneError> {
        let scene = self.static_scene_at(time_s, cut_side, cx)?;
        Ok(render(&scene, cx, settings)?)
    }

    /// Render the explicit static diagnostic layer, when configured, together
    /// with the beauty scene. Cinematic beauty rendering never calls this API.
    pub fn render_static_with_debug_at(
        &self,
        time_s: f64,
        cut_side: CutSide,
        settings: &Settings,
        cx: &Cx<'_>,
    ) -> Result<Film, EulerSceneError> {
        let scene = self.static_scene_with_debug_at(time_s, cut_side, cx)?;
        Ok(render(&scene, cx, settings)?)
    }
}

impl EulerDebugLayer {
    fn primitive(&self) -> Primitive {
        Primitive {
            shape: Shape::Instance(self.instance.clone()),
            material: self.material,
            emission: None,
        }
    }
}

/// Structured fail-closed scene composition diagnostics.
#[derive(Debug)]
pub enum EulerSceneError {
    /// Execution scope requested cancellation.
    Cancelled,
    /// Scene configuration did not use SI metres.
    UnsupportedLengthUnit,
    /// A named configuration field was outside its finite admitted domain.
    InvalidConfig(&'static str),
    /// Artifact metadata named a different resolved asset.
    AssetIdentityMismatch(&'static str),
    /// Artifact mechanics properties differed from the resolved profile.
    MassPropertiesMismatch,
    /// V1 cannot safely interpolate a producer-declared discontinuity.
    DeclaredDiscontinuityUnsupported,
    /// The requested tessellation exceeded a hard topology ceiling.
    TessellationBudgetExceeded {
        /// Projected vertices.
        vertices: usize,
        /// Projected triangles.
        triangles: usize,
    },
    /// Deterministic mesh construction produced no surface.
    EmptyPreviewMesh,
    /// Bounded allocation failed.
    Capacity(&'static str),
    /// Closed contact did not agree with the displaced base plane.
    ContactBaseMismatch(usize),
    /// Retained support gap/contact did not agree with the bound exact chart.
    ContactSpecimenMismatch(usize),
    /// A configured diagnostic marker selected no closed contact geometry.
    MissingDebugContact(usize),
    /// Camera shots do not continuously cover the accepted trajectory horizon.
    CameraCoverage,
    /// Subject bounds violated the declared near/far depth requirement.
    CameraDepthRange,
    /// Endpoint angular speed and sample duration make shortest-arc motion
    /// obviously vulnerable to temporal aliasing.
    AngularSamplingAmbiguous {
        /// Left source interval index.
        interval: usize,
        /// Conservative endpoint diagnostic (rad).
        angular_step_rad: f64,
        /// Configured maximum (rad).
        maximum_rad: f64,
    },
    /// A frame with multiple event-delimited shutters needs explicit weighted
    /// film composition.
    ExposureNeedsComposition {
        /// Number of event-delimited films requiring weighted composition.
        segment_count: usize,
    },
    /// Prepared frame came from a different scene identity.
    PreparedFrameMismatch,
    /// Segment selection was outside a prepared frame.
    InvalidPreparedSegment {
        /// Refused segment index.
        index: usize,
        /// Number of prepared segments available.
        segment_count: usize,
    },
    /// A scene built by this module unexpectedly contained an unsupported shape.
    UnexpectedSceneShape,
    /// EXR encoding refused.
    ImageEncoding,
    /// Renderer instance admission/evaluation refused.
    Instance(InstanceError),
    /// Animated-instance admission/evaluation refused.
    Animated(AnimatedInstanceError),
    /// Camera admission/evaluation refused.
    Camera(CameraError),
    /// Euler shutter or timeline motion admission refused.
    Motion(RenderMotionBridgeError),
    /// Spectral tracer refused.
    Tracer(TracerError),
    /// Explicit tile-render execution, memory admission, or tracer execution
    /// refused while retaining its structured renderer diagnostic.
    RenderExecution(RenderExecutionError),
    /// Conservative animated-instance bounds refused.
    MotionBounds(MotionBoundsError),
}

impl fmt::Display for EulerSceneError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for EulerSceneError {}

impl From<InstanceError> for EulerSceneError {
    fn from(error: InstanceError) -> Self {
        Self::Instance(error)
    }
}

impl From<AnimatedInstanceError> for EulerSceneError {
    fn from(error: AnimatedInstanceError) -> Self {
        Self::Animated(error)
    }
}

impl From<CameraError> for EulerSceneError {
    fn from(error: CameraError) -> Self {
        Self::Camera(error)
    }
}

impl From<RenderMotionBridgeError> for EulerSceneError {
    fn from(error: RenderMotionBridgeError) -> Self {
        Self::Motion(error)
    }
}

impl From<TracerError> for EulerSceneError {
    fn from(error: TracerError) -> Self {
        Self::Tracer(error)
    }
}

impl From<RenderExecutionError> for EulerSceneError {
    fn from(error: RenderExecutionError) -> Self {
        Self::RenderExecution(error)
    }
}

impl From<MotionBoundsError> for EulerSceneError {
    fn from(error: MotionBoundsError) -> Self {
        Self::MotionBounds(error)
    }
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), EulerSceneError> {
    cx.checkpoint().map_err(|_| EulerSceneError::Cancelled)
}

fn validate_config(config: &EulerSceneConfig) -> Result<(), EulerSceneError> {
    if config.length_unit != EulerSceneLengthUnit::Metres {
        return Err(EulerSceneError::UnsupportedLengthUnit);
    }
    if !(8..=MAX_EULER_AZIMUTHAL_SEGMENTS).contains(&config.tessellation.azimuthal_segments) {
        return Err(EulerSceneError::InvalidConfig("azimuthal_segments"));
    }
    if !(1..=MAX_EULER_ARC_SUBDIVISIONS).contains(&config.tessellation.arc_subdivisions_per_arc) {
        return Err(EulerSceneError::InvalidConfig("arc_subdivisions_per_arc"));
    }
    for (field, value) in [
        ("plate_width_m", config.base.plate_width_m),
        ("plate_depth_m", config.base.plate_depth_m),
        ("plate_thickness_m", config.base.plate_thickness_m),
        ("housing_width_m", config.base.housing_width_m),
        ("housing_depth_m", config.base.housing_depth_m),
        ("housing_height_m", config.base.housing_height_m),
    ] {
        if !value.is_finite() || value <= 0.0 {
            return Err(EulerSceneError::InvalidConfig(field));
        }
    }
    if !config.base.housing_gap_m.is_finite() || config.base.housing_gap_m < 0.0 {
        return Err(EulerSceneError::InvalidConfig("housing_gap_m"));
    }
    if !config.camera_near_m.is_finite()
        || !config.camera_far_m.is_finite()
        || config.camera_near_m <= 0.0
        || config.camera_far_m <= config.camera_near_m
    {
        return Err(EulerSceneError::InvalidConfig("camera depth range"));
    }
    if !config.maximum_angular_step_rad.is_finite()
        || config.maximum_angular_step_rad <= 0.0
        || config.maximum_angular_step_rad >= core::f64::consts::PI
    {
        return Err(EulerSceneError::InvalidConfig("maximum_angular_step_rad"));
    }
    validate_style(config.disc_material, "disc_material")?;
    validate_style(config.plate_material, "plate_material")?;
    validate_style(config.housing_material, "housing_material")?;
    validate_light(config.light)?;
    let beauty_ids = [
        config.object_ids.disc,
        config.object_ids.base_plate,
        config.object_ids.housing,
    ];
    let mut unique = BTreeSet::new();
    if beauty_ids.iter().any(|id| *id == 0 || !unique.insert(*id)) {
        return Err(EulerSceneError::InvalidConfig("object_ids"));
    }
    if let EulerDebugOverlay::ContactMarker { radius_m, .. } = config.debug_overlay {
        if !radius_m.is_finite() || radius_m <= 0.0 {
            return Err(EulerSceneError::InvalidConfig("debug marker radius_m"));
        }
        if config.object_ids.debug_marker == 0 || !unique.insert(config.object_ids.debug_marker) {
            return Err(EulerSceneError::InvalidConfig("debug marker object_id"));
        }
    }
    Ok(())
}

fn validate_style(style: EulerMaterialStyle, field: &'static str) -> Result<(), EulerSceneError> {
    match style {
        EulerMaterialStyle::Lambertian { linear_rgb } => {
            validate_linear_rgb(linear_rgb, field)?;
        }
        EulerMaterialStyle::Ggx { linear_rgb, alpha } => {
            validate_linear_rgb(linear_rgb, field)?;
            if !alpha.is_finite() || alpha <= 0.0 || alpha > 1.0 {
                return Err(EulerSceneError::InvalidConfig(field));
            }
        }
        EulerMaterialStyle::Conductor { .. } => {
            // Conductor fields are private validated value types, so no
            // malformed optical state can be constructed here.
        }
        EulerMaterialStyle::Dielectric { .. } => {
            // Dielectric fields are private validated value types, so no
            // malformed optical state can be constructed here.
        }
    }
    Ok(())
}

fn validate_linear_rgb(linear_rgb: [f64; 3], field: &'static str) -> Result<(), EulerSceneError> {
    if linear_rgb
        .iter()
        .any(|value| !value.is_finite() || !(0.0..=1.0).contains(value))
    {
        Err(EulerSceneError::InvalidConfig(field))
    } else {
        Ok(())
    }
}

fn validate_light(light: EulerRectLightSpec) -> Result<(), EulerSceneError> {
    let values = [
        light.corner_world_m.x,
        light.corner_world_m.y,
        light.corner_world_m.z,
        light.edge_u_world_m.x,
        light.edge_u_world_m.y,
        light.edge_u_world_m.z,
        light.edge_v_world_m.x,
        light.edge_v_world_m.y,
        light.edge_v_world_m.z,
        light.radiance_scale,
    ];
    if values.iter().any(|value| !value.is_finite())
        || light.radiance_scale <= 0.0
        || light
            .linear_rgb
            .iter()
            .any(|value| !value.is_finite() || !(0.0..=1.0).contains(value))
        || geom_cross(light.edge_u_world_m, light.edge_v_world_m).norm() <= 0.0
    {
        return Err(EulerSceneError::InvalidConfig("light"));
    }
    Ok(())
}

fn validate_asset_binding(
    artifact: &EulerRenderTrajectoryArtifact,
    specimen: &ResolvedDiscProfile,
) -> Result<ResolvedDiscProfileIdentities, EulerSceneError> {
    let identities = specimen.content_identities();
    let metadata = artifact.trajectory().metadata();
    if metadata.world_frame != RenderWorldFrame::RightHandedZUp {
        return Err(EulerSceneError::AssetIdentityMismatch("world_frame"));
    }
    if metadata.units != RenderUnitSystem::SiRadians {
        return Err(EulerSceneError::AssetIdentityMismatch("units"));
    }
    if metadata.specimen_chart_identity != identities.chart {
        return Err(EulerSceneError::AssetIdentityMismatch(
            "specimen_chart_identity",
        ));
    }
    if metadata.specimen_profile_identity != identities.profile {
        return Err(EulerSceneError::AssetIdentityMismatch(
            "specimen_profile_identity",
        ));
    }
    if metadata.mass_properties.identity != identities.mass_properties {
        return Err(EulerSceneError::AssetIdentityMismatch(
            "mass_properties.identity",
        ));
    }
    let expected = MassProperties::new(
        specimen.mass_properties.mass,
        MbdVec3::ZERO,
        MbdVec3::new(
            specimen.mass_properties.principal_inertia.transverse,
            specimen.mass_properties.principal_inertia.transverse,
            specimen.mass_properties.principal_inertia.axial,
        ),
    )
    .map_err(|_| EulerSceneError::MassPropertiesMismatch)?;
    if metadata.mass_properties.properties != expected
        || specimen.mass_properties.center_of_mass.x.to_bits() != 0.0_f64.to_bits()
        || specimen.mass_properties.center_of_mass.y.to_bits() != 0.0_f64.to_bits()
    {
        return Err(EulerSceneError::MassPropertiesMismatch);
    }
    Ok(identities)
}

fn validate_angular_sampling(
    artifact: &EulerRenderTrajectoryArtifact,
    maximum_rad: f64,
) -> Result<(), EulerSceneError> {
    let samples = artifact.trajectory().samples();
    let mass = artifact.trajectory().metadata().mass_properties.properties;
    for (interval, pair) in samples.windows(2).enumerate() {
        let duration_s = pair[1].input().time_s - pair[0].input().time_s;
        let left = stable_mbd_norm(
            mass.angular_velocity_body_checked(pair[0].state().angular_momentum_body())
                .map_err(|_| EulerSceneError::MassPropertiesMismatch)?,
        )
        .ok_or(EulerSceneError::MassPropertiesMismatch)?;
        let right = stable_mbd_norm(
            mass.angular_velocity_body_checked(pair[1].state().angular_momentum_body())
                .map_err(|_| EulerSceneError::MassPropertiesMismatch)?,
        )
        .ok_or(EulerSceneError::MassPropertiesMismatch)?;
        let left_orientation = pair[0].state().pose().orientation().components();
        let right_orientation = pair[1].state().pose().orientation().components();
        let quaternion_dot = left_orientation
            .iter()
            .zip(right_orientation)
            .fold(0.0, |sum, (left, right)| left.mul_add(right, sum))
            .abs()
            .clamp(0.0, 1.0);
        let shortest_arc_rad = 2.0 * quaternion_dot.acos();
        let angular_step_rad = (duration_s * left.max(right)).max(shortest_arc_rad);
        if !angular_step_rad.is_finite() || angular_step_rad > maximum_rad {
            return Err(EulerSceneError::AngularSamplingAmbiguous {
                interval,
                angular_step_rad,
                maximum_rad,
            });
        }
    }
    Ok(())
}

fn validate_contact_base_alignment(
    artifact: &EulerRenderTrajectoryArtifact,
    specimen: &ResolvedDiscProfile,
    cx: &Cx<'_>,
) -> Result<(), EulerSceneError> {
    let metadata = artifact.trajectory().metadata();
    let normal = metadata
        .base_frame
        .orientation_base_to_world
        .rotate_body_to_world(MbdVec3::new(0.0, 0.0, 1.0));
    for (index, sample) in artifact.trajectory().samples().iter().enumerate() {
        checkpoint(cx)?;
        let displacement = sample
            .input()
            .base_mode
            .ok_or(EulerSceneError::ContactBaseMismatch(index))?
            .displacement_m;
        let plane_point = metadata
            .base_frame
            .origin_world_m
            .add(normal.scale(displacement));
        let state_pose = sample.state().pose();
        let relative_pose = MbdPose::new(
            state_pose.position_world().sub(plane_point),
            state_pose.orientation(),
        )
        .map_err(|_| EulerSceneError::ContactSpecimenMismatch(index))?;
        let derived =
            profile_contact_geometry(&specimen.chart, specimen.mass_properties, relative_pose, cx)
                .map_err(|_| EulerSceneError::ContactSpecimenMismatch(index))?;
        let derived_gap_m = derived.contact.gap_m;
        let derived_scale_m = derived_gap_m.abs().max(sample.input().signed_gap_m.abs());
        let derived_tolerance_m =
            CONTACT_BASE_ALIGNMENT_TOLERANCE_M.max(256.0 * f64::EPSILON * derived_scale_m);
        if (derived_gap_m - sample.input().signed_gap_m).abs() > derived_tolerance_m {
            return Err(EulerSceneError::ContactSpecimenMismatch(index));
        }
        if sample.input().contact_branch != RenderContactBranch::Closed {
            continue;
        }
        let Some(contact) = sample.input().contact_geometry else {
            return Err(EulerSceneError::ContactBaseMismatch(index));
        };
        let observed_gap = contact.point_world_m.sub(plane_point).dot(normal);
        let scale_m = observed_gap.abs().max(sample.input().signed_gap_m.abs());
        let tolerance_m = CONTACT_BASE_ALIGNMENT_TOLERANCE_M.max(256.0 * f64::EPSILON * scale_m);
        if (observed_gap - sample.input().signed_gap_m).abs() > tolerance_m
            || contact.normal_world.dot(normal) < 1.0 - CONTACT_NORMAL_ALIGNMENT_TOLERANCE
        {
            return Err(EulerSceneError::ContactBaseMismatch(index));
        }
        let derived_point_world = derived.contact.point_world_m.add(plane_point);
        let retained_point_delta = contact.point_world_m.sub(derived_point_world);
        if stable_mbd_norm(retained_point_delta)
            .is_none_or(|distance| distance > CONTACT_BASE_ALIGNMENT_TOLERANCE_M)
            || contact.support_feature
                != crate::render_trajectory::RenderSupportFeature::ProfileFeature(
                    derived.support_source_feature,
                )
        {
            return Err(EulerSceneError::ContactSpecimenMismatch(index));
        }
    }
    Ok(())
}

fn validate_camera_coverage(
    artifact: &EulerRenderTrajectoryArtifact,
    camera: &AnimatedCamera,
    cx: &Cx<'_>,
) -> Result<(), EulerSceneError> {
    let samples = artifact.trajectory().samples();
    let first = samples[0].input().time_s;
    let last = samples[samples.len() - 1].input().time_s;
    let mut covered_until = first;
    for shot in camera.shots() {
        if shot.end_s() < first {
            continue;
        }
        if shot.start_s() > last {
            break;
        }
        if shot.start_s() > covered_until {
            return Err(EulerSceneError::CameraCoverage);
        }
        covered_until = covered_until.max(shot.end_s());
        if covered_until >= last {
            break;
        }
    }
    if covered_until < last {
        return Err(EulerSceneError::CameraCoverage);
    }
    for sample in samples {
        checkpoint(cx)?;
        camera
            .evaluate(cx, sample.input().time_s, CutSide::After)
            .map_err(|_| EulerSceneError::CameraCoverage)?;
    }
    Ok(())
}

fn tessellate_disc(
    specimen: &ResolvedDiscProfile,
    identities: ResolvedDiscProfileIdentities,
    config: EulerTessellationConfig,
    cx: &Cx<'_>,
) -> Result<(TriMesh, EulerPreviewMeshReceipt), EulerSceneError> {
    checkpoint(cx)?;
    let mut meridian = Vec::new();
    let mut maximum_meridian_chord_error_m = 0.0_f64;
    for (segment_index, segment) in specimen.chart.segments().iter().copied().enumerate() {
        checkpoint(cx)?;
        let start = segment_start(segment);
        if segment_index == 0 {
            meridian.push(start);
        } else if meridian.last().copied() != Some(start) {
            return Err(EulerSceneError::InvalidConfig("disconnected meridian"));
        }
        match segment {
            MeridianSegment::Line { end, .. } => meridian.push(end),
            MeridianSegment::Arc {
                start,
                end,
                center,
                clockwise,
            } => {
                let subdivisions = config.arc_subdivisions_per_arc;
                let start_angle = (start.axial - center.axial).atan2(start.radius - center.radius);
                let end_angle = (end.axial - center.axial).atan2(end.radius - center.radius);
                let sweep = directed_sweep(start_angle, end_angle, clockwise);
                let radius = (start.radius - center.radius).hypot(start.axial - center.axial);
                let half_chord_angle = 0.5 * sweep.abs() / f64::from(subdivisions);
                maximum_meridian_chord_error_m =
                    maximum_meridian_chord_error_m.max(radius * (1.0 - half_chord_angle.cos()));
                for step in 1..=subdivisions {
                    if step % 64 == 0 {
                        checkpoint(cx)?;
                    }
                    if step == subdivisions {
                        meridian.push(end);
                    } else {
                        let angle = start_angle + sweep * f64::from(step) / f64::from(subdivisions);
                        meridian.push(MeridianPoint::new(
                            center.radius + radius * angle.cos(),
                            center.axial + radius * angle.sin(),
                        ));
                    }
                }
            }
        }
    }
    if meridian.len() < 2 {
        return Err(EulerSceneError::EmptyPreviewMesh);
    }
    let azimuthal = usize::try_from(config.azimuthal_segments)
        .map_err(|_| EulerSceneError::InvalidConfig("azimuthal_segments"))?;
    let vertex_count = meridian.len().checked_mul(azimuthal).ok_or(
        EulerSceneError::TessellationBudgetExceeded {
            vertices: usize::MAX,
            triangles: usize::MAX,
        },
    )?;
    let mut triangle_count = 0usize;
    for pair in meridian.windows(2) {
        triangle_count = triangle_count
            .checked_add(if pair[0].radius == 0.0 && pair[1].radius == 0.0 {
                0
            } else if pair[0].radius == 0.0 || pair[1].radius == 0.0 {
                azimuthal
            } else {
                2 * azimuthal
            })
            .ok_or(EulerSceneError::TessellationBudgetExceeded {
                vertices: vertex_count,
                triangles: usize::MAX,
            })?;
    }
    if vertex_count > MAX_EULER_PREVIEW_VERTICES
        || triangle_count > MAX_EULER_PREVIEW_TRIANGLES
        || vertex_count > usize::try_from(u32::MAX).unwrap_or(usize::MAX)
    {
        return Err(EulerSceneError::TessellationBudgetExceeded {
            vertices: vertex_count,
            triangles: triangle_count,
        });
    }
    let com = specimen.mass_properties.center_of_mass;
    let mut vertices = Vec::new();
    vertices
        .try_reserve_exact(vertex_count)
        .map_err(|_| EulerSceneError::Capacity("preview vertices"))?;
    for (ring, point) in meridian.iter().enumerate() {
        if ring % 64 == 0 {
            checkpoint(cx)?;
        }
        for azimuth in 0..azimuthal {
            let azimuth_u32 = u32::try_from(azimuth)
                .map_err(|_| EulerSceneError::InvalidConfig("azimuthal_segments"))?;
            let angle = core::f64::consts::TAU * f64::from(azimuth_u32)
                / f64::from(config.azimuthal_segments);
            vertices.push([
                point.radius * angle.cos() - com.x,
                point.radius * angle.sin() - com.y,
                point.axial - com.z,
            ]);
        }
    }
    let mut triangles = Vec::new();
    triangles
        .try_reserve_exact(triangle_count)
        .map_err(|_| EulerSceneError::Capacity("preview triangles"))?;
    for (strip, pair) in meridian.windows(2).enumerate() {
        if strip % 64 == 0 {
            checkpoint(cx)?;
        }
        let first = strip * azimuthal;
        let second = (strip + 1) * azimuthal;
        for azimuth in 0..azimuthal {
            let next = (azimuth + 1) % azimuthal;
            let a = u32::try_from(first + azimuth).map_err(|_| {
                EulerSceneError::TessellationBudgetExceeded {
                    vertices: vertex_count,
                    triangles: triangle_count,
                }
            })?;
            let a_next = u32::try_from(first + next).map_err(|_| {
                EulerSceneError::TessellationBudgetExceeded {
                    vertices: vertex_count,
                    triangles: triangle_count,
                }
            })?;
            let b = u32::try_from(second + azimuth).map_err(|_| {
                EulerSceneError::TessellationBudgetExceeded {
                    vertices: vertex_count,
                    triangles: triangle_count,
                }
            })?;
            let b_next = u32::try_from(second + next).map_err(|_| {
                EulerSceneError::TessellationBudgetExceeded {
                    vertices: vertex_count,
                    triangles: triangle_count,
                }
            })?;
            match (pair[0].radius == 0.0, pair[1].radius == 0.0) {
                (true, true) => {}
                (true, false) => triangles.push([a, b_next, b]),
                (false, true) => triangles.push([a, a_next, b]),
                (false, false) => {
                    triangles.push([a, a_next, b_next]);
                    triangles.push([a, b_next, b]);
                }
            }
        }
    }
    if vertices.is_empty() || triangles.is_empty() {
        return Err(EulerSceneError::EmptyPreviewMesh);
    }
    let local_bounds_m = bounds_of_vertices(&vertices)?;
    let mut mesh_identity = DomainHasher::new(EULER_PREVIEW_MESH_IDENTITY_DOMAIN);
    mesh_identity.update(identities.chart.as_bytes());
    mesh_identity.update(&EULER_RENDER_SCENE_BRIDGE_VERSION.to_le_bytes());
    mesh_identity.update(&config.azimuthal_segments.to_le_bytes());
    mesh_identity.update(&config.arc_subdivisions_per_arc.to_le_bytes());
    let vertex_count_u64 = u64::try_from(vertices.len())
        .map_err(|_| EulerSceneError::Capacity("preview vertex identity count"))?;
    let triangle_count_u64 = u64::try_from(triangles.len())
        .map_err(|_| EulerSceneError::Capacity("preview triangle identity count"))?;
    mesh_identity.update(&vertex_count_u64.to_le_bytes());
    mesh_identity.update(&triangle_count_u64.to_le_bytes());
    for vertex in &vertices {
        for component in vertex {
            mesh_identity.update(&component.to_bits().to_le_bytes());
        }
    }
    for triangle in &triangles {
        for index in triangle {
            mesh_identity.update(&index.to_le_bytes());
        }
    }
    let mesh_identity = mesh_identity.finalize();
    let mesh = TriMesh::new(vertices, triangles);
    let source_support = specimen.chart.support();
    let source_outer_radius_m = source_support.max.x.abs().max(source_support.min.x.abs());
    let receipt = EulerPreviewMeshReceipt {
        source_chart_identity: identities.chart,
        mesh_identity,
        azimuthal_segments: config.azimuthal_segments,
        arc_subdivisions_per_arc: config.arc_subdivisions_per_arc,
        vertex_count,
        triangle_count,
        maximum_meridian_chord_error_m,
        maximum_azimuthal_chord_error_m: source_outer_radius_m
            * (1.0 - (core::f64::consts::PI / f64::from(config.azimuthal_segments)).cos()),
        local_bounds_m,
        bvh_fingerprint: mesh.bvh_fingerprint(),
    };
    checkpoint(cx)?;
    Ok((mesh, receipt))
}

fn segment_start(segment: MeridianSegment) -> MeridianPoint {
    match segment {
        MeridianSegment::Line { start, .. } | MeridianSegment::Arc { start, .. } => start,
    }
}

fn directed_sweep(start: f64, end: f64, clockwise: bool) -> f64 {
    if clockwise {
        -((start - end).rem_euclid(core::f64::consts::TAU))
    } else {
        (end - start).rem_euclid(core::f64::consts::TAU)
    }
}

fn bounds_of_vertices(vertices: &[[f64; 3]]) -> Result<Aabb, EulerSceneError> {
    let first = vertices.first().ok_or(EulerSceneError::EmptyPreviewMesh)?;
    let mut min = *first;
    let mut max = *first;
    for vertex in vertices {
        if vertex.iter().any(|value| !value.is_finite()) {
            return Err(EulerSceneError::InvalidConfig("preview vertex"));
        }
        for axis in 0..3 {
            min[axis] = min[axis].min(vertex[axis]);
            max[axis] = max[axis].max(vertex[axis]);
        }
    }
    Ok(Aabb::new(
        Point3::new(min[0], min[1], min[2]),
        Point3::new(max[0], max[1], max[2]),
    ))
}

fn disc_transform_trajectory(
    artifact: &EulerRenderTrajectoryArtifact,
    cx: &Cx<'_>,
) -> Result<RigidTransformTrajectory, EulerSceneError> {
    let trajectory = artifact.trajectory();
    let motion = EulerRenderMotionBridge::new(trajectory);
    let mass = trajectory.metadata().mass_properties.properties;
    let mut keyframes = Vec::new();
    keyframes
        .try_reserve_exact(trajectory.samples().len())
        .map_err(|_| EulerSceneError::Capacity("disc keyframes"))?;
    for sample in trajectory.samples() {
        checkpoint(cx)?;
        let time_s = sample.input().time_s;
        let mapped = motion.sample_at_time(time_s, EventEvaluationSide::RightLimit)?;
        let velocity = sample
            .state()
            .center_of_mass_velocity_world(mass)
            .map_err(|_| EulerSceneError::MassPropertiesMismatch)?;
        keyframes.push(TransformKeyframe::try_new(
            time_s,
            mapped.transform(),
            [velocity.x, velocity.y, velocity.z],
        )?);
    }
    Ok(RigidTransformTrajectory::try_new(keyframes)?)
}

fn base_transform_trajectory(
    artifact: &EulerRenderTrajectoryArtifact,
    cx: &Cx<'_>,
) -> Result<RigidTransformTrajectory, EulerSceneError> {
    let mut keyframes = Vec::new();
    keyframes
        .try_reserve_exact(artifact.trajectory().samples().len())
        .map_err(|_| EulerSceneError::Capacity("base keyframes"))?;
    let orientation = artifact
        .trajectory()
        .metadata()
        .base_frame
        .orientation_base_to_world;
    for sample in artifact.trajectory().samples() {
        checkpoint(cx)?;
        let base = sample
            .input()
            .base_mode
            .ok_or(EulerSceneError::InvalidConfig("missing base mode"))?;
        let velocity =
            orientation.rotate_body_to_world(MbdVec3::new(0.0, 0.0, base.velocity_m_per_s));
        keyframes.push(TransformKeyframe::try_new(
            sample.input().time_s,
            base_transform(artifact, base.displacement_m)?,
            [velocity.x, velocity.y, velocity.z],
        )?);
    }
    Ok(RigidTransformTrajectory::try_new(keyframes)?)
}

fn base_transform(
    artifact: &EulerRenderTrajectoryArtifact,
    displacement_m: f64,
) -> Result<RigidTransform, EulerSceneError> {
    let frame = artifact.trajectory().metadata().base_frame;
    let [w, x, y, z] = frame.orientation_base_to_world.components();
    let displacement_world = frame
        .orientation_base_to_world
        .rotate_body_to_world(MbdVec3::new(0.0, 0.0, displacement_m));
    let origin = frame.origin_world_m.add(displacement_world);
    Ok(RigidTransform::try_new(
        [x, y, z, w],
        [origin.x, origin.y, origin.z],
    )?)
}

fn nominal_base_transform(
    artifact: &EulerRenderTrajectoryArtifact,
) -> Result<RigidTransform, EulerSceneError> {
    base_transform(artifact, 0.0)
}

fn material(style: EulerMaterialStyle) -> Material {
    match style {
        EulerMaterialStyle::Lambertian { linear_rgb } => Material::Lambertian {
            reflectance: lift_rgb(linear_rgb),
        },
        EulerMaterialStyle::Ggx { linear_rgb, alpha } => Material::Ggx {
            reflectance: lift_rgb(linear_rgb),
            alpha,
        },
        EulerMaterialStyle::Conductor { optics, surface } => {
            Material::Conductor { optics, surface }
        }
        EulerMaterialStyle::Dielectric { glass, surface } => {
            Material::Dielectric { glass, surface }
        }
    }
}

fn box_mesh(half_x: f64, half_y: f64, z_min: f64, z_max: f64) -> TriMesh {
    TriMesh::new(
        vec![
            [-half_x, -half_y, z_min],
            [half_x, -half_y, z_min],
            [half_x, half_y, z_min],
            [-half_x, half_y, z_min],
            [-half_x, -half_y, z_max],
            [half_x, -half_y, z_max],
            [half_x, half_y, z_max],
            [-half_x, half_y, z_max],
        ],
        vec![
            [0, 2, 1],
            [0, 3, 2],
            [4, 5, 6],
            [4, 6, 7],
            [0, 1, 5],
            [0, 5, 4],
            [1, 2, 6],
            [1, 6, 5],
            [2, 3, 7],
            [2, 7, 6],
            [3, 0, 4],
            [3, 4, 7],
        ],
    )
}

fn box_identity(label: &[u8], values: [f64; 4]) -> ContentHash {
    let mut hasher = DomainHasher::new("org.frankensim.fs-euler-disc-e2e.scene-box.v1");
    hasher.update(label);
    for value in values {
        hasher.update(&value.to_bits().to_le_bytes());
    }
    hasher.finalize()
}

fn rect_mesh(light: EulerRectLightSpec) -> TriMesh {
    let corner = light.corner_world_m;
    let u = light.edge_u_world_m;
    let v = light.edge_v_world_m;
    let uv = Vec3::new(u.x + v.x, u.y + v.y, u.z + v.z);
    TriMesh::new(
        vec![
            [corner.x, corner.y, corner.z],
            [corner.x + u.x, corner.y + u.y, corner.z + u.z],
            [corner.x + uv.x, corner.y + uv.y, corner.z + uv.z],
            [corner.x + v.x, corner.y + v.y, corner.z + v.z],
        ],
        vec![[0, 1, 2], [0, 2, 3]],
    )
}

fn geom_cross(left: Vec3, right: Vec3) -> Vec3 {
    Vec3::new(
        left.y.mul_add(right.z, -(left.z * right.y)),
        left.z.mul_add(right.x, -(left.x * right.z)),
        left.x.mul_add(right.y, -(left.y * right.x)),
    )
}

fn stable_mbd_norm(vector: MbdVec3) -> Option<f64> {
    if !vector.is_finite() {
        return None;
    }
    let scale = vector.x.abs().max(vector.y.abs()).max(vector.z.abs());
    if scale == 0.0 {
        return Some(0.0);
    }
    let scaled = MbdVec3::new(vector.x / scale, vector.y / scale, vector.z / scale);
    let magnitude = scale * scaled.dot(scaled).sqrt();
    magnitude.is_finite().then_some(magnitude)
}

fn build_debug_layer(
    artifact: &EulerRenderTrajectoryArtifact,
    overlay: EulerDebugOverlay,
    object_id: u64,
    beauty_scene_identity: ContentHash,
) -> Result<Option<EulerDebugLayer>, EulerSceneError> {
    let EulerDebugOverlay::ContactMarker {
        sample_index,
        radius_m,
    } = overlay
    else {
        return Ok(None);
    };
    let contact = artifact
        .trajectory()
        .samples()
        .get(sample_index)
        .and_then(|sample| sample.input().contact_geometry)
        .ok_or(EulerSceneError::MissingDebugContact(sample_index))?;
    let mesh = octahedron_mesh(radius_m);
    let mut identity_bytes = Vec::with_capacity(16);
    identity_bytes.extend_from_slice(&sample_index.to_le_bytes());
    identity_bytes.extend_from_slice(&radius_m.to_bits().to_le_bytes());
    let geometry_identity = hash_domain(
        "org.frankensim.fs-euler-disc-e2e.contact-marker.v1",
        &identity_bytes,
    );
    let transform = RigidTransform::try_new(
        [0.0, 0.0, 0.0, 1.0],
        [
            contact.point_world_m.x,
            contact.point_world_m.y,
            contact.point_world_m.z,
        ],
    )?;
    let instance = GeometryInstance::try_new(
        object_id,
        geometry_identity,
        SharedGeometry::mesh(mesh),
        transform,
    )?;
    let mut layer_identity =
        DomainHasher::new("org.frankensim.fs-euler-disc-e2e.euler-debug-contact-layer.v1");
    layer_identity.update(beauty_scene_identity.as_bytes());
    layer_identity.update(&object_id.to_le_bytes());
    let sample_index_u64 = u64::try_from(sample_index)
        .map_err(|_| EulerSceneError::Capacity("debug source sample index"))?;
    layer_identity.update(&sample_index_u64.to_le_bytes());
    layer_identity.update(&radius_m.to_bits().to_le_bytes());
    layer_identity.update(instance.frame_identity().as_bytes());
    Ok(Some(EulerDebugLayer {
        receipt: EulerDebugLayerReceipt {
            layer_identity: layer_identity.finalize(),
            object_id,
            source_sample_index: sample_index,
            radius_m,
        },
        instance,
        material: Material::Lambertian {
            reflectance: lift_rgb([0.95, 0.08, 0.03]),
        },
    }))
}

fn octahedron_mesh(radius: f64) -> TriMesh {
    TriMesh::new(
        vec![
            [radius, 0.0, 0.0],
            [-radius, 0.0, 0.0],
            [0.0, radius, 0.0],
            [0.0, -radius, 0.0],
            [0.0, 0.0, radius],
            [0.0, 0.0, -radius],
        ],
        vec![
            [4, 0, 2],
            [4, 2, 1],
            [4, 1, 3],
            [4, 3, 0],
            [5, 2, 0],
            [5, 1, 2],
            [5, 3, 1],
            [5, 0, 3],
        ],
    )
}

fn legacy_camera(camera: &fs_render::camera::PhysicalCamera) -> Camera {
    Camera {
        eye: camera.eye(),
        forward: camera.forward(),
        up: camera.up(),
        half_tan: camera.projection().vertical_half_tan(),
    }
}

fn subject_bounds(
    artifact: &EulerRenderTrajectoryArtifact,
    preview: EulerPreviewMeshReceipt,
    base: &EulerBaseVisualSpec,
    disc: &AnimatedGeometryInstance,
    plate: &AnimatedGeometryInstance,
    cx: &Cx<'_>,
) -> Result<Aabb, EulerSceneError> {
    let samples = artifact.trajectory().samples();
    let first_time_s = samples[0].input().time_s;
    let last_time_s = samples[samples.len() - 1].input().time_s;
    let horizon = ShotTimeBounds::try_new(first_time_s, last_time_s)
        .map_err(|_| EulerSceneError::InvalidConfig("trajectory horizon"))?;
    let shutter = ShutterInterval::resolve(
        first_time_s,
        last_time_s - first_time_s,
        ShutterConvention::FrontLoaded,
        ShutterDistribution::UniformCounterV1,
        horizon,
    )
    .map_err(|_| EulerSceneError::InvalidConfig("trajectory horizon"))?;
    checkpoint(cx)?;
    let mut bounds = conservative_trajectory_swept_aabb(
        FiniteLocalAabb::try_new(preview.local_bounds_m)?,
        disc.trajectory(),
        shutter,
    )?;
    let plate_local = Aabb::new(
        Point3::new(
            -0.5 * base.plate_width_m,
            -0.5 * base.plate_depth_m,
            -base.plate_thickness_m,
        ),
        Point3::new(0.5 * base.plate_width_m, 0.5 * base.plate_depth_m, 0.0),
    );
    bounds = bounds.union(&conservative_trajectory_swept_aabb(
        FiniteLocalAabb::try_new(plate_local)?,
        plate.trajectory(),
        shutter,
    )?);
    checkpoint(cx)?;
    let housing_top = -base.plate_thickness_m - base.housing_gap_m;
    let housing_bottom = housing_top - base.housing_height_m;
    Ok(bounds.union(&transformed_box_bounds(
        nominal_base_transform(artifact)?,
        0.5 * base.housing_width_m,
        0.5 * base.housing_depth_m,
        housing_bottom,
        housing_top,
    )))
}

fn transformed_box_bounds(
    transform: RigidTransform,
    half_x: f64,
    half_y: f64,
    z_min: f64,
    z_max: f64,
) -> Aabb {
    let mut points = Vec::with_capacity(8);
    for x in [-half_x, half_x] {
        for y in [-half_y, half_y] {
            for z in [z_min, z_max] {
                points.push(transform.transform_point(Point3::new(x, y, z)));
            }
        }
    }
    let mut min = points[0];
    let mut max = points[0];
    for point in points {
        min.x = min.x.min(point.x);
        min.y = min.y.min(point.y);
        min.z = min.z.min(point.z);
        max.x = max.x.max(point.x);
        max.y = max.y.max(point.y);
        max.z = max.z.max(point.z);
    }
    Aabb::new(min, max)
}

fn validate_camera_depths(
    artifact: &EulerRenderTrajectoryArtifact,
    camera: &AnimatedCamera,
    bounds: Aabb,
    near_m: f64,
    far_m: f64,
    cx: &Cx<'_>,
) -> Result<(), EulerSceneError> {
    for sample in artifact.trajectory().samples() {
        checkpoint(cx)?;
        let physical = camera.evaluate(cx, sample.input().time_s, CutSide::After)?;
        validate_camera_depth(&physical, bounds, near_m, far_m)?;
    }
    let first = artifact.trajectory().samples()[0].input().time_s;
    let last = artifact.trajectory().samples()[artifact.trajectory().samples().len() - 1]
        .input()
        .time_s;
    for shot in camera.shots() {
        for keyframe in shot.keyframes() {
            if keyframe.absolute_time_s() >= first && keyframe.absolute_time_s() <= last {
                checkpoint(cx)?;
                validate_camera_depth(keyframe.camera(), bounds, near_m, far_m)?;
            }
        }
    }
    Ok(())
}

fn validate_camera_depth(
    physical: &fs_render::camera::PhysicalCamera,
    bounds: Aabb,
    near_m: f64,
    far_m: f64,
) -> Result<(), EulerSceneError> {
    for x in [bounds.min.x, bounds.max.x] {
        for y in [bounds.min.y, bounds.max.y] {
            for z in [bounds.min.z, bounds.max.z] {
                let depth = physical
                    .forward()
                    .dot(Point3::new(x, y, z).delta_from(physical.eye()));
                if !depth.is_finite() || depth < near_m || depth > far_m {
                    return Err(EulerSceneError::CameraDepthRange);
                }
            }
        }
    }
    Ok(())
}

fn scene_identity(
    artifact: &EulerRenderTrajectoryArtifact,
    identities: ResolvedDiscProfileIdentities,
    preview: EulerPreviewMeshReceipt,
    plate_geometry_identity: ContentHash,
    housing_geometry_identity: ContentHash,
    config: &EulerSceneConfig,
) -> ContentHash {
    let mut hasher = DomainHasher::new(EULER_RENDER_SCENE_IDENTITY_DOMAIN);
    hasher.update(&EULER_RENDER_SCENE_BRIDGE_VERSION.to_le_bytes());
    hasher.update(artifact.receipt().artifact_identity().as_bytes());
    hasher.update(identities.chart.as_bytes());
    hasher.update(identities.profile.as_bytes());
    hasher.update(identities.mass_properties.as_bytes());
    hasher.update(preview.mesh_identity.as_bytes());
    hasher.update(plate_geometry_identity.as_bytes());
    hasher.update(housing_geometry_identity.as_bytes());
    for id in [
        config.object_ids.disc,
        config.object_ids.base_plate,
        config.object_ids.housing,
    ] {
        hasher.update(&id.to_le_bytes());
    }
    hash_base_config(&mut hasher, config);
    hash_material(&mut hasher, config.disc_material);
    hash_material(&mut hasher, config.plate_material);
    hash_material(&mut hasher, config.housing_material);
    hash_light(&mut hasher, config.light);
    hash_camera(&mut hasher, &config.camera);
    hasher.update(&config.camera_near_m.to_bits().to_le_bytes());
    hasher.update(&config.camera_far_m.to_bits().to_le_bytes());
    hasher.update(&config.maximum_angular_step_rad.to_bits().to_le_bytes());
    hasher.finalize()
}

fn configuration_identity(config: &EulerSceneConfig) -> ContentHash {
    let mut hasher = DomainHasher::new(EULER_RENDER_CONFIGURATION_IDENTITY_DOMAIN);
    hasher.update(&EULER_RENDER_SCENE_BRIDGE_VERSION.to_le_bytes());
    hasher.update(&[match config.length_unit {
        EulerSceneLengthUnit::Metres => 0,
        EulerSceneLengthUnit::Millimetres => 1,
    }]);
    hash_base_config(&mut hasher, config);
    for id in [
        config.object_ids.disc,
        config.object_ids.base_plate,
        config.object_ids.housing,
        config.object_ids.debug_marker,
    ] {
        hasher.update(&id.to_le_bytes());
    }
    hash_material(&mut hasher, config.disc_material);
    hash_material(&mut hasher, config.plate_material);
    hash_material(&mut hasher, config.housing_material);
    hash_light(&mut hasher, config.light);
    hash_camera(&mut hasher, &config.camera);
    hasher.update(&config.camera_near_m.to_bits().to_le_bytes());
    hasher.update(&config.camera_far_m.to_bits().to_le_bytes());
    hasher.update(&config.maximum_angular_step_rad.to_bits().to_le_bytes());
    match config.debug_overlay {
        EulerDebugOverlay::None => hasher.update(&[0]),
        EulerDebugOverlay::ContactMarker {
            sample_index,
            radius_m,
        } => {
            hasher.update(&[1]);
            hasher.update(
                &u64::try_from(sample_index)
                    .unwrap_or(u64::MAX)
                    .to_le_bytes(),
            );
            hasher.update(&radius_m.to_bits().to_le_bytes());
        }
    }
    hasher.finalize()
}

fn hash_base_config(hasher: &mut DomainHasher, config: &EulerSceneConfig) {
    hasher.update(&config.tessellation.azimuthal_segments.to_le_bytes());
    hasher.update(&config.tessellation.arc_subdivisions_per_arc.to_le_bytes());
    for value in [
        config.base.plate_width_m,
        config.base.plate_depth_m,
        config.base.plate_thickness_m,
        config.base.housing_width_m,
        config.base.housing_depth_m,
        config.base.housing_height_m,
        config.base.housing_gap_m,
    ] {
        hasher.update(&value.to_bits().to_le_bytes());
    }
}

fn hash_material(hasher: &mut DomainHasher, style: EulerMaterialStyle) {
    match style {
        EulerMaterialStyle::Lambertian { linear_rgb } => {
            // Preserve the v1 opaque-material byte order so existing scene
            // identities do not change merely because dielectric tag 2 exists.
            hasher.update(&[0]);
            hash_linear_rgb(hasher, linear_rgb);
        }
        EulerMaterialStyle::Ggx { linear_rgb, alpha } => {
            hasher.update(&[1]);
            hash_linear_rgb(hasher, linear_rgb);
            hasher.update(&alpha.to_bits().to_le_bytes());
        }
        EulerMaterialStyle::Dielectric { glass, surface } => {
            hasher.update(&[2]);
            for coefficient in glass.ior().coefficients() {
                hasher.update(&coefficient.to_bits().to_le_bytes());
            }
            match glass.absorption().parameters() {
                BeerLambertParameters::Clear => hasher.update(&[0]),
                BeerLambertParameters::Constant { extinction_per_m } => {
                    hasher.update(&[1]);
                    hasher.update(&extinction_per_m.to_bits().to_le_bytes());
                }
                BeerLambertParameters::ReferenceRgb {
                    linear_rgb,
                    distance_m,
                } => {
                    hasher.update(&[2]);
                    hash_linear_rgb(hasher, linear_rgb);
                    hasher.update(&distance_m.to_bits().to_le_bytes());
                }
            }
            hasher.update(&[match glass.provenance() {
                GlassProvenance::Custom => 0,
                GlassProvenance::RepresentativeBorosilicateV1 => 1,
                GlassProvenance::RepresentativeCrownV1 => 2,
            }]);
            match surface.roughness_alpha() {
                None => hasher.update(&[0]),
                Some(alpha) => {
                    hasher.update(&[1]);
                    hasher.update(&alpha.to_bits().to_le_bytes());
                }
            }
        }
        EulerMaterialStyle::Conductor { optics, surface } => {
            // Tags 0--2 remain frozen for existing Euler configurations.
            // Bind the generic L5 material identity so changes to optical,
            // Fresnel, or roughness semantics invalidate this L6 scene.
            hasher.update(&[3]);
            hasher.update(
                Material::Conductor { optics, surface }
                    .content_identity()
                    .as_bytes(),
            );
        }
    }
}

fn hash_linear_rgb(hasher: &mut DomainHasher, linear_rgb: [f64; 3]) {
    for value in linear_rgb {
        hasher.update(&value.to_bits().to_le_bytes());
    }
}

fn hash_light(hasher: &mut DomainHasher, light: EulerRectLightSpec) {
    for value in [
        light.corner_world_m.x,
        light.corner_world_m.y,
        light.corner_world_m.z,
        light.edge_u_world_m.x,
        light.edge_u_world_m.y,
        light.edge_u_world_m.z,
        light.edge_v_world_m.x,
        light.edge_v_world_m.y,
        light.edge_v_world_m.z,
        light.linear_rgb[0],
        light.linear_rgb[1],
        light.linear_rgb[2],
        light.radiance_scale,
    ] {
        hasher.update(&value.to_bits().to_le_bytes());
    }
}

fn hash_camera(hasher: &mut DomainHasher, camera: &AnimatedCamera) {
    hasher.update(
        &u64::try_from(camera.shots().len())
            .unwrap_or(u64::MAX)
            .to_le_bytes(),
    );
    for shot in camera.shots() {
        hasher.update(&shot.shot_id().to_le_bytes());
        hasher.update(&shot.start_s().to_bits().to_le_bytes());
        hasher.update(&shot.end_s().to_bits().to_le_bytes());
        hasher.update(
            &u64::try_from(shot.keyframes().len())
                .unwrap_or(u64::MAX)
                .to_le_bytes(),
        );
        for keyframe in shot.keyframes() {
            hasher.update(&keyframe.absolute_time_s().to_bits().to_le_bytes());
            let camera = keyframe.camera();
            let projection = camera.projection();
            for value in [
                camera.eye().x,
                camera.eye().y,
                camera.eye().z,
                camera.forward().x,
                camera.forward().y,
                camera.forward().z,
                camera.up().x,
                camera.up().y,
                camera.up().z,
                projection.vertical_half_tan(),
                camera.focus_distance_m(),
                camera.aperture().radius_m(),
                camera.exposure_metadata().sensitivity_iso(),
                camera.exposure_metadata().compensation_ev(),
            ] {
                hasher.update(&value.to_bits().to_le_bytes());
            }
            if let (Some(focal_length_m), Some(sensor_height_m)) =
                (projection.focal_length_m(), projection.sensor_height_m())
            {
                hasher.update(&[0]);
                hasher.update(&focal_length_m.to_bits().to_le_bytes());
                hasher.update(&sensor_height_m.to_bits().to_le_bytes());
            } else if let Some(vertical_fov_rad) = projection.vertical_fov_rad() {
                hasher.update(&[1]);
                hasher.update(&vertical_fov_rad.to_bits().to_le_bytes());
            } else {
                hasher.update(&[2]);
            }
            match camera.aperture().blades() {
                Some(blades) => {
                    hasher.update(&[1, blades]);
                    hasher.update(
                        &camera
                            .aperture()
                            .rotation_rad()
                            .unwrap_or(0.0)
                            .to_bits()
                            .to_le_bytes(),
                    );
                }
                None => hasher.update(&[0]),
            }
            match keyframe.focus() {
                KeyframeFocus::AxialDistance => hasher.update(&[0]),
                KeyframeFocus::WorldPoint(point) => {
                    hasher.update(&[1]);
                    for value in [point.x, point.y, point.z] {
                        hasher.update(&value.to_bits().to_le_bytes());
                    }
                }
            }
        }
    }
}

/// Small deterministic render settings suitable for the scene-bridge E2E.
#[must_use]
pub const fn euler_scene_smoke_settings(width: u32, height: u32) -> Settings {
    Settings {
        width,
        height,
        spp: 4,
        max_depth: 6,
        sampler: Sampler::OwenSobol,
        strategy: DirectStrategy::Mis,
        seed: 0x4555_4c45_525f_5343,
    }
}
