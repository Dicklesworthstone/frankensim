//! Deterministic scientific overlay and per-frame sidecar generation.
//!
//! The output is a transparent vector-command layer. It is deliberately
//! separate from beauty rendering: callers may composite the commands for a
//! diagnostic cut, while the source trajectory, beauty scene, film, and EXR
//! identities remain untouched.

use core::fmt;
use std::fmt::Write as _;

use fs_blake3::{ContentHash, hash_domain};
use fs_exec::Cx;
use fs_geom::Point3;
use fs_mbd::Vec3 as MbdVec3;
use fs_render::camera::{
    AnimatedCamera, CameraError, CutSide, OpticalCenterProjection, PhysicalCamera,
};

use crate::coupled_runner::{ChannelWrench, ContactTransitionKind};
use crate::render_trajectory::{
    DerivedEulerQois, RenderContactBranch, RenderNormalForceSampling, RenderSupportFeature,
    RenderTrajectoryAuthority, RenderTrajectorySampleInput,
};
use crate::render_trajectory_codec::EulerRenderTrajectoryArtifact;
use crate::timeline_resampling::{
    DeclaredDiscontinuityKind, EventEvaluationSide, TimelineEvent, TimelineResampler,
    TimelineResamplingError, TimelineSampleSource,
};

/// Schema and command semantics version.
pub const EULER_SCIENTIFIC_OVERLAY_SCHEMA_VERSION: u16 = 1;
/// Domain for the canonical sidecar bytes.
pub const EULER_SCIENTIFIC_OVERLAY_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-euler-disc-e2e.scientific-overlay-sidecar.v1";
/// Hard ceiling for retained contact-orbit points in one frame.
pub const MAX_SCIENTIFIC_OVERLAY_ORBIT_POINTS: usize = 65_536;
/// Explicit authority boundary carried by every sidecar.
pub const SCIENTIFIC_OVERLAY_NO_CLAIM: &str = "visualization-only derived from simulation evidence; not measurement, calibration, force reconstruction, uncertainty certification, or beauty-image authority";

/// A fixed color-blind-safe overlay color.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OverlayColor {
    /// Linear display-independent RGB label value in `[0, 1]`.
    pub rgb: [f32; 3],
    /// Overlay opacity in `[0, 1]`.
    pub alpha: f32,
}

impl OverlayColor {
    /// Okabe-Ito sky blue.
    pub const AXIS: Self = Self::rgb8(86, 180, 233);
    /// Okabe-Ito orange.
    pub const CONTACT: Self = Self::rgb8(230, 159, 0);
    /// Okabe-Ito blue.
    pub const FORCE: Self = Self::rgb8(0, 114, 178);
    /// Okabe-Ito reddish purple.
    pub const TORQUE: Self = Self::rgb8(204, 121, 167);
    /// Okabe-Ito bluish green.
    pub const EVENT: Self = Self::rgb8(0, 158, 115);
    /// High-contrast neutral label color.
    pub const LABEL: Self = Self::rgb8(255, 255, 255);

    const fn rgb8(red: u8, green: u8, blue: u8) -> Self {
        Self {
            rgb: [
                red as f32 / 255.0,
                green as f32 / 255.0,
                blue as f32 / 255.0,
            ],
            alpha: 1.0,
        }
    }
}

/// Raster-safe placement rectangle in continuous pixel-edge coordinates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OverlaySafeArea {
    /// Inclusive left edge.
    pub left_px: f64,
    /// Inclusive top edge.
    pub top_px: f64,
    /// Exclusive right edge.
    pub right_px: f64,
    /// Exclusive bottom edge.
    pub bottom_px: f64,
}

impl OverlaySafeArea {
    fn clamp(self, point: [f64; 2]) -> [f64; 2] {
        [
            point[0].clamp(self.left_px, previous_float(self.right_px)),
            point[1].clamp(self.top_px, previous_float(self.bottom_px)),
        ]
    }
}

/// Explicit frame-local scaling and layout controls.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ScientificOverlayConfig {
    /// Beauty raster width.
    pub width: u32,
    /// Beauty raster height.
    pub height: u32,
    /// Equal safe-title inset from every raster edge.
    pub safe_inset_px: u32,
    /// Rendered world length of the unit disc-axis indicator.
    pub axis_length_m: f64,
    /// World metres used per displayed newton.
    pub force_scale_m_per_n: f64,
    /// World metres used per displayed newton-metre.
    pub torque_scale_m_per_nm: f64,
    /// Contact marker radius in raster pixels.
    pub contact_marker_radius_px: f64,
    /// Vector/orbit stroke width in raster pixels.
    pub line_width_px: f64,
    /// Maximum exact closed-contact samples retained in the orbit.
    pub maximum_orbit_points: usize,
}

impl ScientificOverlayConfig {
    /// Conservative defaults for a supplied beauty raster.
    #[must_use]
    pub const fn reference(width: u32, height: u32) -> Self {
        Self {
            width,
            height,
            safe_inset_px: 48,
            axis_length_m: 0.04,
            force_scale_m_per_n: 0.002,
            torque_scale_m_per_nm: 0.02,
            contact_marker_radius_px: 8.0,
            line_width_px: 2.0,
            maximum_orbit_points: 16_384,
        }
    }

    fn validate(self) -> Result<OverlaySafeArea, ScientificOverlayError> {
        let doubled_inset = self
            .safe_inset_px
            .checked_mul(2)
            .ok_or(ScientificOverlayError::InvalidConfig("safe_inset_px"))?;
        if self.width == 0
            || self.height == 0
            || doubled_inset >= self.width
            || doubled_inset >= self.height
        {
            return Err(ScientificOverlayError::InvalidConfig("raster/safe area"));
        }
        for (name, value) in [
            ("axis_length_m", self.axis_length_m),
            ("force_scale_m_per_n", self.force_scale_m_per_n),
            ("torque_scale_m_per_nm", self.torque_scale_m_per_nm),
            ("contact_marker_radius_px", self.contact_marker_radius_px),
            ("line_width_px", self.line_width_px),
        ] {
            if !value.is_finite() || value <= 0.0 {
                return Err(ScientificOverlayError::InvalidConfig(name));
            }
        }
        if self.maximum_orbit_points == 0
            || self.maximum_orbit_points > MAX_SCIENTIFIC_OVERLAY_ORBIT_POINTS
        {
            return Err(ScientificOverlayError::InvalidConfig(
                "maximum_orbit_points",
            ));
        }
        Ok(OverlaySafeArea {
            left_px: f64::from(self.safe_inset_px),
            top_px: f64::from(self.safe_inset_px),
            right_px: f64::from(self.width - self.safe_inset_px),
            bottom_px: f64::from(self.height - self.safe_inset_px),
        })
    }
}

/// Request for one exact diagnostic frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ScientificOverlayRequest {
    /// Timeline time in seconds.
    pub frame_time_s: f64,
    /// Camera ownership when the time lies exactly on a hard cut.
    pub cut_side: CutSide,
    /// Timeline event-side semantics at an exact mechanics event.
    pub event_side: EventEvaluationSide,
}

/// Projection retained for world-to-screen audit and off-screen explanation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum OverlayProjection {
    /// Finite point in front of the optical centre.
    InFront {
        /// NDC coordinate, `+x` right and `+y` up.
        ndc_xy: [f64; 2],
        /// Pixel-edge coordinate, `+x` right and `+y` down.
        pixel_xy: [f64; 2],
        /// Positive axial camera depth in metres.
        depth_m: f64,
        /// Whether the raw point belongs to the half-open beauty raster.
        in_frame: bool,
    },
    /// Point on or behind the optical-centre plane.
    BehindCamera {
        /// Nonpositive axial camera depth in metres.
        signed_depth_m: f64,
    },
}

impl OverlayProjection {
    /// Raw continuous pixel coordinate when the point is in front.
    #[must_use]
    pub const fn pixel_xy(self) -> Option<[f64; 2]> {
        match self {
            Self::InFront { pixel_xy, .. } => Some(pixel_xy),
            Self::BehindCamera { .. } => None,
        }
    }
}

/// One retained point in SI world space and beauty-raster coordinates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ScientificOverlayPoint {
    /// World coordinate in metres.
    pub world_m: [f64; 3],
    /// Optical-centre projection.
    pub projection: OverlayProjection,
}

/// Stable vector vocabulary shared by commands and the sidecar.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ScientificVectorKind {
    /// Unit body symmetry axis, scaled by `axis_length_m`.
    DiscAxis,
    /// Independently retained base-normal scalar.
    NormalContactForce,
    /// Gravity-channel force.
    GravityForce,
    /// Aggregate-contact-channel force.
    ContactForce,
    /// Rolling-resistance-channel force.
    RollingForce,
    /// Reduced-base-channel force.
    BaseForce,
    /// Exterior-gas-channel force.
    GasForce,
    /// Gravity-channel torque.
    GravityTorque,
    /// Aggregate-contact-channel torque.
    ContactTorque,
    /// Rolling-resistance-channel torque.
    RollingTorque,
    /// Reduced-base-channel torque.
    BaseTorque,
    /// Exterior-gas-channel torque.
    GasTorque,
}

impl ScientificVectorKind {
    fn name(self) -> &'static str {
        match self {
            Self::DiscAxis => "disc-axis",
            Self::NormalContactForce => "normal-contact-force",
            Self::GravityForce => "gravity-force",
            Self::ContactForce => "contact-force",
            Self::RollingForce => "rolling-force",
            Self::BaseForce => "base-force",
            Self::GasForce => "gas-force",
            Self::GravityTorque => "gravity-torque",
            Self::ContactTorque => "contact-torque",
            Self::RollingTorque => "rolling-torque",
            Self::BaseTorque => "base-torque",
            Self::GasTorque => "gas-torque",
        }
    }

    fn unit(self) -> &'static str {
        match self {
            Self::DiscAxis => "1",
            Self::NormalContactForce
            | Self::GravityForce
            | Self::ContactForce
            | Self::RollingForce
            | Self::BaseForce
            | Self::GasForce => "N",
            Self::GravityTorque
            | Self::ContactTorque
            | Self::RollingTorque
            | Self::BaseTorque
            | Self::GasTorque => "N m",
        }
    }

    fn color(self) -> OverlayColor {
        match self {
            Self::DiscAxis => OverlayColor::AXIS,
            Self::GravityTorque
            | Self::ContactTorque
            | Self::RollingTorque
            | Self::BaseTorque
            | Self::GasTorque => OverlayColor::TORQUE,
            _ => OverlayColor::FORCE,
        }
    }
}

/// One scientific vector and the exact display scale used for it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ScientificVectorDiagnostic {
    /// Semantic channel and quantity.
    pub kind: ScientificVectorKind,
    /// Vector origin in world metres.
    pub origin_world_m: [f64; 3],
    /// Physical vector components in the declared `unit`.
    pub value_si: [f64; 3],
    /// Display metres per physical unit.
    pub display_scale_m_per_unit: f64,
    /// Projected origin.
    pub start: OverlayProjection,
    /// Projected scaled endpoint.
    pub end: OverlayProjection,
}

/// Current unilateral-contact information.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ScientificContactDiagnostic {
    /// One-sided contact branch at the frame time.
    pub branch: RenderContactBranch,
    /// Exact source geometry, unavailable for an interpolated or open sample.
    pub point: Option<ScientificOverlayPoint>,
    /// Exact support feature paired with `point`.
    pub support_feature: Option<RenderSupportFeature>,
}

/// Exact-source energy data. Interpolated frames do not fabricate it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ScientificEnergyDiagnostic {
    /// Total declared mechanical energy [J].
    pub mechanical_energy_j: f64,
    /// Declared energy-closure defect [J].
    pub energy_defect_j: f64,
    /// Gravity/contact/rolling/base/gas channel work [J].
    pub channel_work_j: [f64; 5],
}

/// Event kind used by the time strip and sidecar.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ScientificEventKind {
    /// Contact opens.
    ContactOpening,
    /// Contact reimpacts.
    ContactReimpact,
    /// Terminal inclination threshold.
    TerminalInclination,
    /// Declared continuation seam.
    ContinuationSeam,
    /// Other producer-declared discontinuity.
    ProducerDeclared,
}

impl ScientificEventKind {
    fn name(self) -> &'static str {
        match self {
            Self::ContactOpening => "contact-opening",
            Self::ContactReimpact => "contact-reimpact",
            Self::TerminalInclination => "terminal-inclination",
            Self::ContinuationSeam => "continuation-seam",
            Self::ProducerDeclared => "producer-declared",
        }
    }
}

/// Retained event time and numerical bracket.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ScientificEventDiagnostic {
    /// Stable event kind.
    pub kind: ScientificEventKind,
    /// Retained event time [s].
    pub time_s: f64,
    /// Inclusive numerical bracket start [s].
    pub bracket_start_s: f64,
    /// Inclusive numerical bracket end [s].
    pub bracket_end_s: f64,
}

/// Machine-readable values used to build one overlay frame.
#[derive(Clone, Debug, PartialEq)]
pub struct ScientificOverlaySidecar {
    /// Source trajectory artifact identity.
    pub source_trajectory_identity: ContentHash,
    /// Unchanged beauty-scene identity supplied by composition.
    pub beauty_scene_identity: ContentHash,
    /// Unchanged source authority ceiling.
    pub source_authority: RenderTrajectoryAuthority,
    /// Frame time [s].
    pub frame_time_s: f64,
    /// Stable camera shot identity.
    pub camera_shot_id: u64,
    /// Raster width.
    pub width: u32,
    /// Raster height.
    pub height: u32,
    /// Deterministic safe placement rectangle.
    pub safe_area: OverlaySafeArea,
    /// Exact or interpolated timeline provenance.
    pub timeline_source: TimelineSampleSource,
    /// Continuous one-mode base displacement [m].
    pub base_displacement_m: f64,
    /// Continuous one-mode base velocity [m/s].
    pub base_velocity_m_per_s: f64,
    /// Current contact state.
    pub contact: ScientificContactDiagnostic,
    /// Exact closed-contact history through this frame.
    pub contact_orbit: Vec<ScientificOverlayPoint>,
    /// Disc axis, always available from the continuous pose.
    pub disc_axis: ScientificVectorDiagnostic,
    /// Exact-source force and torque channels; empty for interpolated frames.
    pub vectors: Vec<ScientificVectorDiagnostic>,
    /// Exact-source Euler QoIs; unavailable for interpolated frames.
    pub qois: Option<DerivedEulerQois>,
    /// Exact-source energy channels; unavailable for interpolated frames.
    pub energy: Option<ScientificEnergyDiagnostic>,
    /// Events belonging to the source interval of this frame.
    pub events: Vec<ScientificEventDiagnostic>,
}

impl ScientificOverlaySidecar {
    /// Canonical compact JSON. Field order and numeric formatting are fixed by
    /// this schema version; every numeric input was admitted as finite upstream.
    #[must_use]
    pub fn to_canonical_json(&self) -> String {
        let mut out = String::new();
        write!(
            out,
            "{{\"schema_version\":{},\"visualization_only\":true,\"transparent_background\":true,\"source_authority\":\"simulation-evidence\",\"source_trajectory_identity\":\"{}\",\"beauty_scene_identity\":\"{}\",\"frame_time_s\":{},\"camera_shot_id\":{},\"raster\":{{\"width\":{},\"height\":{}}},\"safe_area_px\":[{},{},{},{}],\"timeline_source\":",
            EULER_SCIENTIFIC_OVERLAY_SCHEMA_VERSION,
            self.source_trajectory_identity,
            self.beauty_scene_identity,
            self.frame_time_s,
            self.camera_shot_id,
            self.width,
            self.height,
            self.safe_area.left_px,
            self.safe_area.top_px,
            self.safe_area.right_px,
            self.safe_area.bottom_px,
        )
        .expect("String writes cannot fail");
        write_timeline_source(&mut out, self.timeline_source);
        write!(
            out,
            ",\"base_mode\":{{\"displacement_m\":{},\"velocity_m_per_s\":{}}},\"contact\":{{\"branch\":\"{}\",\"point\":",
            self.base_displacement_m,
            self.base_velocity_m_per_s,
            contact_branch_name(self.contact.branch),
        )
        .expect("String writes cannot fail");
        write_optional_point(&mut out, self.contact.point);
        out.push_str(",\"support_feature\":");
        write_support_feature(&mut out, self.contact.support_feature);
        out.push_str("},\"contact_orbit\":[");
        write_joined(&mut out, &self.contact_orbit, write_point);
        out.push_str("],\"disc_axis\":");
        write_vector(&mut out, &self.disc_axis);
        out.push_str(",\"vectors\":[");
        write_joined(&mut out, &self.vectors, write_vector);
        out.push_str("],\"qois\":");
        write_qois(&mut out, self.qois);
        out.push_str(",\"energy\":");
        write_energy(&mut out, self.energy);
        out.push_str(",\"events\":[");
        write_joined(&mut out, &self.events, write_event);
        write!(out, "],\"no_claim\":\"{}\"}}", SCIENTIFIC_OVERLAY_NO_CLAIM)
            .expect("String writes cannot fail");
        out
    }
}

/// Renderable transparent vector commands. Labels carry explicit rectangles so
/// a downstream compositor need not guess placement or overlap policy.
#[derive(Clone, Debug, PartialEq)]
pub enum ScientificOverlayPrimitive {
    /// Clipped line segment.
    Line {
        /// Semantic role.
        role: &'static str,
        /// Safe-area-clipped endpoints.
        from_px: [f64; 2],
        /// Safe-area-clipped endpoints.
        to_px: [f64; 2],
        /// Stroke color.
        color: OverlayColor,
        /// Stroke width.
        width_px: f64,
    },
    /// Current-contact point marker or off-screen edge indicator.
    Marker {
        /// Semantic role.
        role: &'static str,
        /// Safe-area-clamped centre.
        center_px: [f64; 2],
        /// Marker radius.
        radius_px: f64,
        /// Whether the raw projection was outside the beauty raster.
        clipped: bool,
        /// Stroke color.
        color: OverlayColor,
    },
    /// Event time and its retained numerical bracket on the frame time strip.
    EventBracket {
        /// Event kind.
        kind: ScientificEventKind,
        /// Bracket start x.
        start_x_px: f64,
        /// Event x.
        event_x_px: f64,
        /// Bracket end x.
        end_x_px: f64,
        /// Stacked strip y.
        y_px: f64,
        /// Stroke color.
        color: OverlayColor,
    },
    /// Deterministically placed annotation.
    Label {
        /// Non-overlapping safe-area rectangle `[left, top, right, bottom]`.
        bounds_px: [f64; 4],
        /// Text including explicit SI units.
        text: String,
        /// Text color.
        color: OverlayColor,
    },
}

/// Complete transparent command layer and its byte-identical sidecar payload.
#[derive(Clone, Debug, PartialEq)]
pub struct ScientificOverlayFrame {
    sidecar: ScientificOverlaySidecar,
    sidecar_json: String,
    sidecar_identity: ContentHash,
    primitives: Vec<ScientificOverlayPrimitive>,
}

impl ScientificOverlayFrame {
    /// Typed sidecar used to derive every command.
    #[must_use]
    pub const fn sidecar(&self) -> &ScientificOverlaySidecar {
        &self.sidecar
    }

    /// Canonical machine-readable bytes.
    #[must_use]
    pub fn sidecar_json(&self) -> &str {
        &self.sidecar_json
    }

    /// Domain-separated identity of `sidecar_json`.
    #[must_use]
    pub const fn sidecar_identity(&self) -> ContentHash {
        self.sidecar_identity
    }

    /// Transparent vector commands in deterministic painter order.
    #[must_use]
    pub fn primitives(&self) -> &[ScientificOverlayPrimitive] {
        &self.primitives
    }
}

/// Structured fail-closed overlay diagnostics.
#[derive(Clone, Debug, PartialEq)]
pub enum ScientificOverlayError {
    /// Caller-owned scope requested cancellation.
    Cancelled,
    /// Invalid finite resource or display configuration.
    InvalidConfig(&'static str),
    /// Closed-contact orbit exceeded its explicit point ceiling.
    OrbitPointBudgetExceeded {
        /// Required retained points.
        required: usize,
        /// Configured ceiling.
        limit: usize,
    },
    /// Timeline reconstruction refused.
    Timeline(TimelineResamplingError),
    /// Camera evaluation or projection refused.
    Camera(CameraError),
}

impl fmt::Display for ScientificOverlayError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for ScientificOverlayError {}

impl From<TimelineResamplingError> for ScientificOverlayError {
    fn from(error: TimelineResamplingError) -> Self {
        Self::Timeline(error)
    }
}

impl From<CameraError> for ScientificOverlayError {
    fn from(error: CameraError) -> Self {
        Self::Camera(error)
    }
}

/// Build one diagnostic frame without mutating or re-encoding any beauty/raw
/// artifact. Exact interval forces, energies, QoIs, and contact geometry are
/// emitted only when `frame_time_s` is an exact admitted source sample.
pub fn build_scientific_overlay(
    cx: &Cx<'_>,
    artifact: &EulerRenderTrajectoryArtifact,
    camera: &AnimatedCamera,
    beauty_scene_identity: ContentHash,
    request: ScientificOverlayRequest,
    config: ScientificOverlayConfig,
) -> Result<ScientificOverlayFrame, ScientificOverlayError> {
    checkpoint(cx)?;
    let safe_area = config.validate()?;
    let evaluated = camera.evaluate_with_shot(cx, request.frame_time_s, request.cut_side)?;
    let physical_camera = evaluated.camera();
    let trajectory = artifact.trajectory();
    let timeline =
        TimelineResampler::new(trajectory).sample(request.frame_time_s, request.event_side)?;
    let exact_input = match timeline.source {
        TimelineSampleSource::ExactSample { index } => Some(trajectory.samples()[index].input()),
        TimelineSampleSource::Interpolated { .. } => None,
    };

    let current_contact = exact_input.and_then(|input| input.contact_geometry);
    let contact = ScientificContactDiagnostic {
        branch: timeline.contact_branch,
        point: current_contact
            .map(|geometry| project_point(physical_camera, geometry.point_world_m, config)),
        support_feature: current_contact.map(|geometry| geometry.support_feature),
    };

    let mut contact_orbit = Vec::new();
    for (index, sample) in trajectory.samples().iter().enumerate() {
        if sample.input().time_s > request.frame_time_s {
            break;
        }
        if let Some(geometry) = sample.input().contact_geometry {
            if contact_orbit.len() == config.maximum_orbit_points {
                return Err(ScientificOverlayError::OrbitPointBudgetExceeded {
                    required: contact_orbit.len() + 1,
                    limit: config.maximum_orbit_points,
                });
            }
            if index % 128 == 0 {
                checkpoint(cx)?;
            }
            contact_orbit.push(project_point(
                physical_camera,
                geometry.point_world_m,
                config,
            ));
        }
    }

    let pose = timeline.state.pose();
    let centre = pose.position_world();
    let axis = pose
        .orientation()
        .rotate_body_to_world(MbdVec3::new(0.0, 0.0, 1.0));
    let disc_axis = project_vector(
        physical_camera,
        config,
        ScientificVectorKind::DiscAxis,
        centre,
        axis,
        config.axis_length_m,
    );
    let mut vectors = Vec::new();
    if let Some(input) = exact_input {
        append_exact_vectors(
            &mut vectors,
            physical_camera,
            config,
            trajectory.metadata().channel_availability,
            input,
            centre,
        );
    }
    let qois = exact_input.map(|input| input.qois);
    let energy = exact_input.map(|input| ScientificEnergyDiagnostic {
        mechanical_energy_j: input.mechanical_energy_j,
        energy_defect_j: input.energy_defect_j,
        channel_work_j: [
            input.channels.gravity.work_j,
            input.channels.contact.work_j,
            input.channels.rolling.work_j,
            input.channels.base.work_j,
            input.channels.gas.work_j,
        ],
    });
    let events = timeline
        .interval_events
        .iter()
        .copied()
        .map(event_diagnostic)
        .collect::<Vec<_>>();
    let sidecar = ScientificOverlaySidecar {
        source_trajectory_identity: artifact.receipt().artifact_identity(),
        beauty_scene_identity,
        source_authority: trajectory.metadata().authority,
        frame_time_s: request.frame_time_s,
        camera_shot_id: evaluated.shot_id(),
        width: config.width,
        height: config.height,
        safe_area,
        timeline_source: timeline.source,
        base_displacement_m: timeline.base_mode.displacement_m,
        base_velocity_m_per_s: timeline.base_mode.velocity_m_per_s,
        contact,
        contact_orbit,
        disc_axis,
        vectors,
        qois,
        energy,
        events,
    };
    let sidecar_json = sidecar.to_canonical_json();
    let sidecar_identity = hash_domain(
        EULER_SCIENTIFIC_OVERLAY_IDENTITY_DOMAIN,
        sidecar_json.as_bytes(),
    );
    let primitives = build_primitives(cx, &sidecar, config)?;
    checkpoint(cx)?;
    Ok(ScientificOverlayFrame {
        sidecar,
        sidecar_json,
        sidecar_identity,
        primitives,
    })
}

fn append_exact_vectors(
    vectors: &mut Vec<ScientificVectorDiagnostic>,
    camera: &PhysicalCamera,
    config: ScientificOverlayConfig,
    availability: crate::render_trajectory::RenderChannelAvailability,
    input: &RenderTrajectorySampleInput,
    centre: MbdVec3,
) {
    if availability.normal_force_sampling != RenderNormalForceSampling::Unavailable {
        if let Some(contact) = input.contact_geometry {
            vectors.push(project_vector(
                camera,
                config,
                ScientificVectorKind::NormalContactForce,
                contact.point_world_m,
                contact.normal_world.scale(input.interval_normal_force_n),
                config.force_scale_m_per_n,
            ));
        }
    }
    for (available, force_kind, torque_kind, wrench) in [
        (
            availability.gravity,
            ScientificVectorKind::GravityForce,
            ScientificVectorKind::GravityTorque,
            input.channels.gravity,
        ),
        (
            availability.contact,
            ScientificVectorKind::ContactForce,
            ScientificVectorKind::ContactTorque,
            input.channels.contact,
        ),
        (
            availability.rolling,
            ScientificVectorKind::RollingForce,
            ScientificVectorKind::RollingTorque,
            input.channels.rolling,
        ),
        (
            availability.base,
            ScientificVectorKind::BaseForce,
            ScientificVectorKind::BaseTorque,
            input.channels.base,
        ),
        (
            availability.gas,
            ScientificVectorKind::GasForce,
            ScientificVectorKind::GasTorque,
            input.channels.gas,
        ),
    ] {
        if available {
            append_wrench(
                vectors,
                camera,
                config,
                force_kind,
                torque_kind,
                centre,
                wrench,
            );
        }
    }
}

fn append_wrench(
    vectors: &mut Vec<ScientificVectorDiagnostic>,
    camera: &PhysicalCamera,
    config: ScientificOverlayConfig,
    force_kind: ScientificVectorKind,
    torque_kind: ScientificVectorKind,
    origin: MbdVec3,
    wrench: ChannelWrench,
) {
    vectors.push(project_vector(
        camera,
        config,
        force_kind,
        origin,
        wrench.force_world_n,
        config.force_scale_m_per_n,
    ));
    vectors.push(project_vector(
        camera,
        config,
        torque_kind,
        origin,
        wrench.torque_world_nm,
        config.torque_scale_m_per_nm,
    ));
}

fn project_vector(
    camera: &PhysicalCamera,
    config: ScientificOverlayConfig,
    kind: ScientificVectorKind,
    origin: MbdVec3,
    value: MbdVec3,
    scale: f64,
) -> ScientificVectorDiagnostic {
    let endpoint = origin.add(value.scale(scale));
    ScientificVectorDiagnostic {
        kind,
        origin_world_m: vec_array(origin),
        value_si: vec_array(value),
        display_scale_m_per_unit: scale,
        start: project(camera, origin, config),
        end: project(camera, endpoint, config),
    }
}

fn project_point(
    camera: &PhysicalCamera,
    point: MbdVec3,
    config: ScientificOverlayConfig,
) -> ScientificOverlayPoint {
    ScientificOverlayPoint {
        world_m: vec_array(point),
        projection: project(camera, point, config),
    }
}

fn project(
    camera: &PhysicalCamera,
    point: MbdVec3,
    config: ScientificOverlayConfig,
) -> OverlayProjection {
    match camera
        .project_from_optical_center(
            Point3::new(point.x, point.y, point.z),
            f64::from(config.width) / f64::from(config.height),
        )
        .expect("admitted trajectory and camera must retain finite projection")
    {
        OpticalCenterProjection::BehindCamera { signed_depth_m } => {
            OverlayProjection::BehindCamera { signed_depth_m }
        }
        OpticalCenterProjection::InFront { ndc_xy, depth_m } => {
            let pixel_xy = [
                (ndc_xy[0] + 1.0) * 0.5 * f64::from(config.width),
                (1.0 - ndc_xy[1]) * 0.5 * f64::from(config.height),
            ];
            OverlayProjection::InFront {
                ndc_xy,
                pixel_xy,
                depth_m,
                in_frame: pixel_xy[0] >= 0.0
                    && pixel_xy[0] < f64::from(config.width)
                    && pixel_xy[1] >= 0.0
                    && pixel_xy[1] < f64::from(config.height),
            }
        }
    }
}

fn event_diagnostic(event: TimelineEvent) -> ScientificEventDiagnostic {
    match event {
        TimelineEvent::Contact(event) => ScientificEventDiagnostic {
            kind: match event.kind {
                ContactTransitionKind::Opening => ScientificEventKind::ContactOpening,
                ContactTransitionKind::Reimpact => ScientificEventKind::ContactReimpact,
            },
            time_s: event.time_s,
            bracket_start_s: event.bracket_start_s,
            bracket_end_s: event.bracket_end_s,
        },
        TimelineEvent::TerminalInclination(event) => ScientificEventDiagnostic {
            kind: ScientificEventKind::TerminalInclination,
            time_s: event.time_s,
            bracket_start_s: event.bracket_start_s,
            bracket_end_s: event.bracket_end_s,
        },
        TimelineEvent::Declared(event) => ScientificEventDiagnostic {
            kind: match event.kind {
                DeclaredDiscontinuityKind::ContinuationSeam => {
                    ScientificEventKind::ContinuationSeam
                }
                DeclaredDiscontinuityKind::ProducerDeclared => {
                    ScientificEventKind::ProducerDeclared
                }
            },
            time_s: event.time_s,
            bracket_start_s: event.time_s,
            bracket_end_s: event.time_s,
        },
    }
}

fn build_primitives(
    cx: &Cx<'_>,
    sidecar: &ScientificOverlaySidecar,
    config: ScientificOverlayConfig,
) -> Result<Vec<ScientificOverlayPrimitive>, ScientificOverlayError> {
    let safe = sidecar.safe_area;
    let mut primitives = Vec::new();
    for (index, pair) in sidecar.contact_orbit.windows(2).enumerate() {
        if index % 128 == 0 {
            checkpoint(cx)?;
        }
        append_projected_line(
            &mut primitives,
            "contact-orbit",
            pair[0].projection,
            pair[1].projection,
            OverlayColor::CONTACT,
            config.line_width_px,
            safe,
        );
    }
    if let Some(point) = sidecar.contact.point {
        if let OverlayProjection::InFront {
            pixel_xy, in_frame, ..
        } = point.projection
        {
            primitives.push(ScientificOverlayPrimitive::Marker {
                role: "current-contact",
                center_px: safe.clamp(pixel_xy),
                radius_px: config.contact_marker_radius_px,
                clipped: !in_frame || safe.clamp(pixel_xy) != pixel_xy,
                color: OverlayColor::CONTACT,
            });
        }
    }
    append_vector_primitive(&mut primitives, &sidecar.disc_axis, config, safe);
    for vector in &sidecar.vectors {
        append_vector_primitive(&mut primitives, vector, config, safe);
    }

    let first_time = sidecar.frame_time_s
        - sidecar
            .events
            .iter()
            .map(|event| sidecar.frame_time_s - event.bracket_start_s)
            .fold(0.0_f64, f64::max);
    let strip_span = (sidecar.frame_time_s - first_time).max(f64::EPSILON);
    for (index, event) in sidecar.events.iter().enumerate() {
        let map_x = |time_s: f64| {
            safe.left_px
                + ((time_s - first_time) / strip_span).clamp(0.0, 1.0)
                    * (safe.right_px - safe.left_px)
        };
        let row = u32::try_from(index.min(15)).expect("bounded event row");
        primitives.push(ScientificOverlayPrimitive::EventBracket {
            kind: event.kind,
            start_x_px: map_x(event.bracket_start_s),
            event_x_px: map_x(event.time_s),
            end_x_px: map_x(event.bracket_end_s),
            y_px: safe.bottom_px - 20.0 - 8.0 * f64::from(row),
            color: OverlayColor::EVENT,
        });
    }
    append_labels(&mut primitives, sidecar)?;
    Ok(primitives)
}

fn append_vector_primitive(
    primitives: &mut Vec<ScientificOverlayPrimitive>,
    vector: &ScientificVectorDiagnostic,
    config: ScientificOverlayConfig,
    safe: OverlaySafeArea,
) {
    if vector.value_si == [0.0; 3] {
        return;
    }
    append_projected_line(
        primitives,
        vector.kind.name(),
        vector.start,
        vector.end,
        vector.kind.color(),
        config.line_width_px,
        safe,
    );
}

fn append_projected_line(
    primitives: &mut Vec<ScientificOverlayPrimitive>,
    role: &'static str,
    start: OverlayProjection,
    end: OverlayProjection,
    color: OverlayColor,
    width_px: f64,
    safe: OverlaySafeArea,
) {
    let (Some(start), Some(end)) = (start.pixel_xy(), end.pixel_xy()) else {
        return;
    };
    if let Some((from_px, to_px)) = clip_line(start, end, safe) {
        primitives.push(ScientificOverlayPrimitive::Line {
            role,
            from_px,
            to_px,
            color,
            width_px,
        });
    }
}

fn append_labels(
    primitives: &mut Vec<ScientificOverlayPrimitive>,
    sidecar: &ScientificOverlaySidecar,
) -> Result<(), ScientificOverlayError> {
    let contact = match sidecar.contact.point {
        Some(_) => match sidecar.contact.support_feature {
            Some(RenderSupportFeature::CylinderRim) => "contact: closed, cylinder rim".to_owned(),
            Some(RenderSupportFeature::ProfileFeature(index)) => {
                format!("contact: closed, profile feature {index}")
            }
            None => "contact: closed, feature unavailable".to_owned(),
        },
        None if sidecar.contact.branch == RenderContactBranch::Open => "contact: open".to_owned(),
        None => "contact: closed, interpolated geometry unavailable".to_owned(),
    };
    let mut labels = vec![
        format!("t = {:.6} s", sidecar.frame_time_s),
        contact,
        format!(
            "base q = {:.6e} m; qdot = {:.6e} m/s",
            sidecar.base_displacement_m, sidecar.base_velocity_m_per_s
        ),
        format!(
            "axis scale = {:.6e} m; force scale = declared per vector",
            sidecar.disc_axis.display_scale_m_per_unit
        ),
    ];
    if let Some(qois) = sidecar.qois {
        labels.push(format!(
            "inclination = {:.6e} rad; precession = {:.6e} rad/s",
            qois.inclination_rad, qois.precession_rad_per_s
        ));
        labels.push(format!("spin = {:.6e} rad/s", qois.spin_rad_per_s));
    } else {
        labels.push("QoIs: unavailable at interpolated frame".to_owned());
    }
    if let Some(energy) = sidecar.energy {
        labels.push(format!(
            "energy = {:.6e} J; defect = {:.6e} J",
            energy.mechanical_energy_j, energy.energy_defect_j
        ));
    } else {
        labels.push("energy: unavailable at interpolated frame".to_owned());
    }
    labels.push("SIMULATION EVIDENCE / VISUALIZATION ONLY".to_owned());

    let safe = sidecar.safe_area;
    let line_height_px = 28.0;
    let panel_width_px = 0.46 * (safe.right_px - safe.left_px);
    let left = safe.right_px - panel_width_px;
    let needed_height = line_height_px
        * f64::from(
            u32::try_from(labels.len())
                .map_err(|_| ScientificOverlayError::InvalidConfig("label count"))?,
        );
    if needed_height > safe.bottom_px - safe.top_px {
        return Err(ScientificOverlayError::InvalidConfig("label safe area"));
    }
    for (index, text) in labels.into_iter().enumerate() {
        let row = u32::try_from(index)
            .map_err(|_| ScientificOverlayError::InvalidConfig("label count"))?;
        let top = safe.top_px + line_height_px * f64::from(row);
        primitives.push(ScientificOverlayPrimitive::Label {
            bounds_px: [left, top, safe.right_px, top + line_height_px],
            text,
            color: OverlayColor::LABEL,
        });
    }
    Ok(())
}

fn clip_line(
    start: [f64; 2],
    end: [f64; 2],
    rect: OverlaySafeArea,
) -> Option<([f64; 2], [f64; 2])> {
    let delta = [end[0] - start[0], end[1] - start[1]];
    let mut lower = 0.0_f64;
    let mut upper = 1.0_f64;
    for (p, q) in [
        (-delta[0], start[0] - rect.left_px),
        (delta[0], previous_float(rect.right_px) - start[0]),
        (-delta[1], start[1] - rect.top_px),
        (delta[1], previous_float(rect.bottom_px) - start[1]),
    ] {
        if p == 0.0 {
            if q < 0.0 {
                return None;
            }
            continue;
        }
        let ratio = q / p;
        if p < 0.0 {
            lower = lower.max(ratio);
        } else {
            upper = upper.min(ratio);
        }
        if lower > upper {
            return None;
        }
    }
    Some((
        [start[0] + lower * delta[0], start[1] + lower * delta[1]],
        [start[0] + upper * delta[0], start[1] + upper * delta[1]],
    ))
}

fn write_timeline_source(out: &mut String, source: TimelineSampleSource) {
    match source {
        TimelineSampleSource::ExactSample { index } => {
            write!(out, "{{\"kind\":\"exact\",\"index\":{index}}}")
                .expect("String writes cannot fail");
        }
        TimelineSampleSource::Interpolated {
            left_index,
            right_index,
            alpha,
            ..
        } => {
            write!(out, "{{\"kind\":\"interpolated-cubic-hermite-slerp-v1\",\"left_index\":{left_index},\"right_index\":{right_index},\"alpha\":{alpha}}}")
                .expect("String writes cannot fail");
        }
    }
}

fn write_optional_point(out: &mut String, point: Option<ScientificOverlayPoint>) {
    match point {
        Some(point) => write_point(out, &point),
        None => out.push_str("null"),
    }
}

fn write_point(out: &mut String, point: &ScientificOverlayPoint) {
    out.push_str("{\"world_m\":");
    write_f64_array(out, point.world_m);
    out.push_str(",\"projection\":");
    write_projection(out, point.projection);
    out.push('}');
}

fn write_projection(out: &mut String, projection: OverlayProjection) {
    match projection {
        OverlayProjection::InFront {
            ndc_xy,
            pixel_xy,
            depth_m,
            in_frame,
        } => {
            write!(out, "{{\"kind\":\"in-front\",\"ndc_xy\":[{},{}],\"pixel_xy\":[{},{}],\"depth_m\":{},\"in_frame\":{}}}", ndc_xy[0], ndc_xy[1], pixel_xy[0], pixel_xy[1], depth_m, in_frame)
                .expect("String writes cannot fail");
        }
        OverlayProjection::BehindCamera { signed_depth_m } => {
            write!(
                out,
                "{{\"kind\":\"behind-camera\",\"signed_depth_m\":{signed_depth_m}}}"
            )
            .expect("String writes cannot fail");
        }
    }
}

fn write_support_feature(out: &mut String, feature: Option<RenderSupportFeature>) {
    match feature {
        None => out.push_str("null"),
        Some(RenderSupportFeature::CylinderRim) => out.push_str("{\"kind\":\"cylinder-rim\"}"),
        Some(RenderSupportFeature::ProfileFeature(index)) => {
            write!(out, "{{\"kind\":\"profile-feature\",\"index\":{index}}}")
                .expect("String writes cannot fail");
        }
    }
}

fn write_vector(out: &mut String, vector: &ScientificVectorDiagnostic) {
    write!(
        out,
        "{{\"kind\":\"{}\",\"unit\":\"{}\",\"origin_world_m\":",
        vector.kind.name(),
        vector.kind.unit()
    )
    .expect("String writes cannot fail");
    write_f64_array(out, vector.origin_world_m);
    out.push_str(",\"value_si\":");
    write_f64_array(out, vector.value_si);
    write!(
        out,
        ",\"display_scale_m_per_unit\":{},\"start\":",
        vector.display_scale_m_per_unit
    )
    .expect("String writes cannot fail");
    write_projection(out, vector.start);
    out.push_str(",\"end\":");
    write_projection(out, vector.end);
    out.push('}');
}

fn write_qois(out: &mut String, qois: Option<DerivedEulerQois>) {
    match qois {
        None => out.push_str("null"),
        Some(qois) => {
            write!(out, "{{\"inclination_rad\":{},\"precession_rad_per_s\":{},\"spin_rad_per_s\":{},\"precession_acceleration_rad_per_s2\":{}}}", qois.inclination_rad, qois.precession_rad_per_s, qois.spin_rad_per_s, qois.precession_acceleration_rad_per_s2)
                .expect("String writes cannot fail");
        }
    }
}

fn write_energy(out: &mut String, energy: Option<ScientificEnergyDiagnostic>) {
    match energy {
        None => out.push_str("null"),
        Some(energy) => {
            write!(
                out,
                "{{\"mechanical_energy_j\":{},\"energy_defect_j\":{},\"channel_work_j\":",
                energy.mechanical_energy_j, energy.energy_defect_j
            )
            .expect("String writes cannot fail");
            write_f64_array(out, energy.channel_work_j);
            out.push('}');
        }
    }
}

fn write_event(out: &mut String, event: &ScientificEventDiagnostic) {
    write!(
        out,
        "{{\"kind\":\"{}\",\"time_s\":{},\"bracket_start_s\":{},\"bracket_end_s\":{}}}",
        event.kind.name(),
        event.time_s,
        event.bracket_start_s,
        event.bracket_end_s
    )
    .expect("String writes cannot fail");
}

fn write_joined<T>(out: &mut String, values: &[T], write_value: fn(&mut String, &T)) {
    for (index, value) in values.iter().enumerate() {
        if index > 0 {
            out.push(',');
        }
        write_value(out, value);
    }
}

fn write_f64_array<const N: usize>(out: &mut String, values: [f64; N]) {
    out.push('[');
    for (index, value) in values.into_iter().enumerate() {
        if index > 0 {
            out.push(',');
        }
        write!(out, "{value}").expect("String writes cannot fail");
    }
    out.push(']');
}

fn contact_branch_name(branch: RenderContactBranch) -> &'static str {
    match branch {
        RenderContactBranch::Open => "open",
        RenderContactBranch::Closed => "closed",
    }
}

fn vec_array(vector: MbdVec3) -> [f64; 3] {
    [vector.x, vector.y, vector.z]
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), ScientificOverlayError> {
    cx.checkpoint()
        .map_err(|_| ScientificOverlayError::Cancelled)
}

fn previous_float(value: f64) -> f64 {
    f64::from_bits(value.to_bits() - 1)
}
