//! Focused G0/G3/E2E coverage for the animated Euler scene bridge.

#![cfg(feature = "cinematic-render")]

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::{ContentHash, hash_domain};
use fs_euler_disc_e2e::coupled_runner::{ChannelOwnership, ContactTransitionKind};
use fs_euler_disc_e2e::profile_contact_geometry;
use fs_euler_disc_e2e::render_scene_bridge::{
    EulerCinematicScene, EulerDebugOverlay, EulerFrameRequest, EulerMaterialStyle, EulerSceneError,
    EulerSceneLengthUnit, EulerTessellationConfig, euler_scene_smoke_settings,
};
use fs_euler_disc_e2e::specimen::{DiscProfileSpec, ResolvedDiscProfile};
use fs_euler_disc_e2e::{
    DeclaredDiscontinuityKind, DeclaredTimelineDiscontinuity, DerivedEulerQois,
    EULER_RENDER_TRAJECTORY_SCHEMA_VERSION, EulerRenderTrajectoryArtifact, EventEvaluationSide,
    ExposureEventPolicy, RenderBaseFrame, RenderBaseModeState, RenderChannelAvailability,
    RenderContactBranch, RenderContactGeometry, RenderContactTransition, RenderMassProperties,
    RenderSampleDisposition, RenderSupportFeature, RenderTrajectory, RenderTrajectoryAuthority,
    RenderTrajectoryCodecBudget, RenderTrajectoryMetadata, RenderTrajectorySampleInput,
    RenderUnitSystem, RenderWorldFrame,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3 as GeomVec3};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_render::camera::{
    AnimatedCamera, Aperture, CameraKeyframe, CameraProjection, CameraShot, CutSide, PhysicalCamera,
};
use fs_render::dielectric::{DielectricGlass, DielectricSurface, GlassProvenance};
use fs_render::instances::SharedGeometry;
use fs_render::motion::{ShutterConvention, ShutterDistribution};
use fs_render::tracer::{Material, Shape};
use fs_rep_frep::SquatDiscEdgeTreatment;

const END_TIME_S: f64 = 0.02;
const STEP_S: f64 = 0.01;

fn assert_close(actual: f64, expected: f64, tolerance: f64, context: &str) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{context}: actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:.3e}"
    );
}

fn with_cx<R>(cancelled: bool, operation: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    if cancelled {
        gate.request();
    }
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x4555_4c45_525f_5343,
                kernel_id: 0x4252_4944_4745,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        operation(&cx)
    })
}

fn identity(label: &str) -> ContentHash {
    hash_domain(
        "org.frankensim.test.euler-render-scene.v1",
        label.as_bytes(),
    )
}

fn specimen(cx: &Cx<'_>) -> ResolvedDiscProfile {
    DiscProfileSpec::SolidCylinder {
        outer_radius_m: 0.038,
        thickness_m: 0.006,
        edge_treatment: SquatDiscEdgeTreatment::CircularFillet { radius: 0.001 },
    }
    .resolve(7_800.0, cx)
    .expect("real circular-filleted steel specimen")
}

fn render_mass(specimen: &ResolvedDiscProfile) -> MassProperties {
    MassProperties::new(
        specimen.mass_properties.mass,
        Vec3::ZERO,
        Vec3::new(
            specimen.mass_properties.principal_inertia.transverse,
            specimen.mass_properties.principal_inertia.transverse,
            specimen.mass_properties.principal_inertia.axial,
        ),
    )
    .expect("resolved specimen mass properties")
}

fn state_at(time_s: f64, mass: MassProperties) -> RigidBodyState {
    let orientation = UnitQuaternion::from_axis_angle(Vec3::new(1.0, 0.0, 0.0), 1.0 + 0.2 * time_s)
        .expect("finite tilted orientation");
    let inertia = mass.principal_inertia_body();
    RigidBodyState::new(
        Pose::new(Vec3::new(0.005 * time_s, 0.0, 0.045), orientation).expect("finite pose"),
        Vec3::new(0.005 * mass.mass(), 0.0, 0.0),
        Vec3::new(0.2 * inertia.x, 0.0, 4.0 * inertia.z),
    )
    .expect("finite rigid-body state")
}

fn sample_at(
    interval_start_time_s: f64,
    time_s: f64,
    specimen: &ResolvedDiscProfile,
    mass: MassProperties,
    transition: bool,
    final_sample: bool,
    cx: &Cx<'_>,
) -> RenderTrajectorySampleInput {
    let mut state = state_at(time_s, mass);
    let base_mode = RenderBaseModeState {
        displacement_m: 0.0002 * time_s,
        velocity_m_per_s: 0.0002,
    };
    if transition {
        let provisional =
            profile_contact_geometry(&specimen.chart, specimen.mass_properties, state.pose(), cx)
                .expect("provisional exact chart support");
        let old_pose = state.pose();
        let old_position = old_pose.position_world();
        let corrected_pose = Pose::new(
            Vec3::new(
                old_position.x,
                old_position.y,
                old_position.z + base_mode.displacement_m - provisional.contact.gap_m,
            ),
            old_pose.orientation(),
        )
        .expect("contact-corrected pose");
        state = RigidBodyState::new(
            corrected_pose,
            state.linear_momentum_world(),
            state.angular_momentum_body(),
        )
        .expect("contact-corrected state");
    }
    let orientation = state.pose().orientation();
    let exact_contact =
        profile_contact_geometry(&specimen.chart, specimen.mass_properties, state.pose(), cx)
            .expect("exact chart support");
    let exact_gap_m = exact_contact.contact.gap_m - base_mode.displacement_m;
    let (contact_branch, contact_geometry, signed_gap_m, interval_contact_active, transitions) =
        if transition {
            (
                RenderContactBranch::Closed,
                Some(RenderContactGeometry {
                    point_world_m: exact_contact.contact.point_world_m,
                    normal_world: Vec3::new(0.0, 0.0, 1.0),
                    support_feature: RenderSupportFeature::ProfileFeature(
                        exact_contact.support_source_feature,
                    ),
                }),
                exact_gap_m,
                true,
                vec![RenderContactTransition {
                    kind: ContactTransitionKind::Reimpact,
                    time_s: STEP_S,
                    bracket_start_s: 0.009,
                    bracket_end_s: 0.011,
                }],
            )
        } else {
            (
                RenderContactBranch::Open,
                None,
                exact_gap_m,
                false,
                Vec::new(),
            )
        };
    RenderTrajectorySampleInput {
        interval_start_time_s,
        time_s,
        world_frame: RenderWorldFrame::RightHandedZUp,
        units: RenderUnitSystem::SiRadians,
        center_of_mass_world_m: state.pose().position_world(),
        orientation_body_to_world: orientation.components(),
        linear_momentum_world_kg_m_per_s: state.linear_momentum_world(),
        angular_momentum_body_kg_m2_per_s: state.angular_momentum_body(),
        symmetry_axis_world: orientation.rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0)),
        contact_branch,
        contact_geometry,
        signed_gap_m,
        interval_contact_active,
        interval_normal_force_n: if transition { 1.0 } else { 0.0 },
        contact_transitions: transitions,
        base_mode: Some(base_mode),
        channels: ChannelOwnership::default(),
        mechanical_energy_j: 1.0,
        energy_defect_j: 0.0,
        qois: DerivedEulerQois::from_state(state, mass, 0.0).expect("finite Euler QoIs"),
        disposition: if final_sample {
            RenderSampleDisposition::HorizonCensored
        } else {
            RenderSampleDisposition::Continue
        },
        terminal_event: None,
    }
}

fn trajectory(
    specimen: &ResolvedDiscProfile,
    chart_identity_override: Option<ContentHash>,
    with_transition: bool,
    linear_velocity_x_mps: Option<&[f64]>,
    cx: &Cx<'_>,
) -> RenderTrajectory {
    let mass = render_mass(specimen);
    let inputs = if with_transition {
        vec![
            sample_at(0.0, 0.0, specimen, mass, false, false, cx),
            sample_at(0.0, END_TIME_S, specimen, mass, true, true, cx),
        ]
    } else {
        vec![
            sample_at(0.0, 0.0, specimen, mass, false, false, cx),
            sample_at(0.0, STEP_S, specimen, mass, false, false, cx),
            sample_at(STEP_S, END_TIME_S, specimen, mass, false, true, cx),
        ]
    };
    let mut inputs = inputs;
    if let Some(velocities) = linear_velocity_x_mps {
        assert_eq!(velocities.len(), inputs.len());
        for (input, velocity_mps) in inputs.iter_mut().zip(velocities) {
            input.linear_momentum_world_kg_m_per_s =
                Vec3::new(velocity_mps * mass.mass(), 0.0, 0.0);
        }
    }
    let identities = specimen.content_identities();
    let first = &inputs[0];
    let initial_orientation = UnitQuaternion::new(
        first.orientation_body_to_world[0],
        first.orientation_body_to_world[1],
        first.orientation_body_to_world[2],
        first.orientation_body_to_world[3],
    )
    .expect("initial orientation");
    let initial_state = RigidBodyState::new(
        Pose::new(first.center_of_mass_world_m, initial_orientation).expect("initial pose"),
        first.linear_momentum_world_kg_m_per_s,
        first.angular_momentum_body_kg_m2_per_s,
    )
    .expect("initial state");
    RenderTrajectory::try_new(
        RenderTrajectoryMetadata {
            schema_version: EULER_RENDER_TRAJECTORY_SCHEMA_VERSION,
            world_frame: RenderWorldFrame::RightHandedZUp,
            units: RenderUnitSystem::SiRadians,
            specimen_profile_identity: identities.profile,
            specimen_chart_identity: chart_identity_override.unwrap_or(identities.chart),
            mass_properties: RenderMassProperties {
                identity: identities.mass_properties,
                properties: mass,
            },
            initial_state,
            initial_base_mode: first.base_mode.expect("fixture base state"),
            base_model_identity: identity("base"),
            base_frame: RenderBaseFrame {
                origin_world_m: Vec3::ZERO,
                orientation_base_to_world: UnitQuaternion::IDENTITY,
            },
            model_identity: identity("model"),
            channel_availability: RenderChannelAvailability::NONE_AVAILABLE,
            configuration_identity: identity("configuration"),
            configuration_fingerprint: 0x4555_4c45_525f_4532,
            timestep_s: if with_transition { END_TIME_S } else { STEP_S },
            producer_version: "scene-bridge-test-v1".into(),
            applicability: "deterministic visualization bridge fixture only".into(),
            no_claims: vec!["rendering does not validate the reduced mechanics".into()],
            authority: RenderTrajectoryAuthority::SimulationEvidence,
        },
        inputs,
    )
    .expect("valid render trajectory fixture")
}

fn artifact(
    specimen: &ResolvedDiscProfile,
    chart_identity_override: Option<ContentHash>,
    with_transition: bool,
    discontinuities: Vec<DeclaredTimelineDiscontinuity>,
    cx: &Cx<'_>,
) -> EulerRenderTrajectoryArtifact {
    EulerRenderTrajectoryArtifact::try_from_trajectory(
        identity(if with_transition {
            "transition-campaign"
        } else {
            "smooth-campaign"
        }),
        trajectory(specimen, chart_identity_override, with_transition, None, cx),
        discontinuities,
        RenderTrajectoryCodecBudget::DEFAULT,
        cx,
    )
    .expect("canonical trajectory artifact")
}

fn overshoot_artifact(
    specimen: &ResolvedDiscProfile,
    cx: &Cx<'_>,
) -> EulerRenderTrajectoryArtifact {
    EulerRenderTrajectoryArtifact::try_from_trajectory(
        identity("hermite-overshoot-campaign"),
        trajectory(specimen, None, false, Some(&[100.0, -100.0, 100.0]), cx),
        Vec::new(),
        RenderTrajectoryCodecBudget::DEFAULT,
        cx,
    )
    .expect("overshoot trajectory artifact")
}

fn physical_camera() -> PhysicalCamera {
    let eye = Point3::new(0.24, -0.30, 0.18);
    let target = Point3::new(0.0, 0.0, 0.025);
    let focus_distance_m = target.delta_from(eye).norm();
    PhysicalCamera::try_look_at(
        eye,
        target,
        GeomVec3::new(0.0, 0.0, 1.0),
        CameraProjection::try_half_tangent(0.48).expect("projection"),
        focus_distance_m,
        Aperture::try_circular(0.0).expect("pinhole"),
    )
    .expect("scene camera")
}

fn camera() -> AnimatedCamera {
    AnimatedCamera::try_static(0x4555_4c45_525f_4341, 0.0, END_TIME_S, physical_camera())
        .expect("static shot")
}

fn config() -> fs_euler_disc_e2e::render_scene_bridge::EulerSceneConfig {
    let mut config = fs_euler_disc_e2e::render_scene_bridge::EulerSceneConfig::reference(camera());
    config.tessellation = EulerTessellationConfig {
        azimuthal_segments: 32,
        arc_subdivisions_per_arc: 8,
    };
    config
}

fn frame_request(event_policy: ExposureEventPolicy) -> EulerFrameRequest {
    EulerFrameRequest {
        frame_time_s: STEP_S,
        exposure_duration_s: 0.0,
        convention: ShutterConvention::Centered,
        distribution: ShutterDistribution::UniformCounterV1,
        event_policy,
        cut_side: CutSide::After,
    }
}

#[test]
fn g0_real_filleted_asset_builds_one_deterministic_com_centered_scene() {
    with_cx(false, |cx| {
        let specimen = specimen(cx);
        let artifact = artifact(&specimen, None, false, Vec::new(), cx);
        let scene_config = config();
        let first = EulerCinematicScene::try_build(&artifact, &specimen, scene_config.clone(), cx)
            .expect("first scene build");
        let replay = EulerCinematicScene::try_build(&artifact, &specimen, scene_config.clone(), cx)
            .expect("deterministic scene replay");

        assert_eq!(first.scene_identity(), replay.scene_identity());
        assert_eq!(first.preview_mesh_receipt(), replay.preview_mesh_receipt());
        assert_eq!(
            first.source_authority(),
            RenderTrajectoryAuthority::SimulationEvidence
        );
        assert_eq!(
            first.source_trajectory_identity(),
            artifact.receipt().artifact_identity()
        );

        let receipt = first.preview_mesh_receipt();
        assert_eq!(
            receipt.source_chart_identity,
            specimen.content_identities().chart
        );
        assert!(receipt.vertex_count > 0 && receipt.triangle_count > 0);
        assert!(receipt.maximum_meridian_chord_error_m > 0.0);
        assert!(receipt.maximum_meridian_chord_error_m < 0.001);
        assert!(receipt.maximum_azimuthal_chord_error_m > 0.0);
        assert_eq!(
            receipt.local_bounds_m.min.x.to_bits(),
            (-0.038_f64).to_bits()
        );
        assert_eq!(receipt.local_bounds_m.max.x.to_bits(), 0.038_f64.to_bits());
        assert_close(
            receipt.local_bounds_m.min.z,
            -0.003,
            2.0e-18,
            "lower cap bound",
        );
        assert_close(
            receipt.local_bounds_m.max.z,
            0.003,
            2.0e-18,
            "upper cap bound",
        );

        let indices = first.primitive_indices();
        assert_eq!(indices.disc, 0);
        assert_eq!(indices.base_plate, 1);
        assert_eq!(indices.housing, 2);
        assert_eq!(indices.light, 3);
        assert_eq!(first.debug_layer_receipt(), None);
        assert_eq!(first.scene().primitives.len(), 4);
        assert!(matches!(
            first.scene().primitives[indices.disc].material,
            Material::Ggx { alpha, .. } if alpha.to_bits() == 0.12_f64.to_bits()
        ));
        assert!(matches!(
            first.scene().primitives[indices.base_plate].material,
            Material::Dielectric { glass, surface }
                if glass.provenance() == GlassProvenance::RepresentativeCrownV1
                    && surface.roughness_alpha().map(f64::to_bits)
                        == Some(0.06_f64.to_bits())
        ));
        assert!(matches!(
            first.scene().primitives[indices.housing].material,
            Material::Ggx { alpha, .. } if alpha.to_bits() == 0.28_f64.to_bits()
        ));
        assert!(first.scene().primitives[indices.light].emission.is_some());
        let Shape::AnimatedInstance(disc) = &first.scene().primitives[indices.disc].shape else {
            assert!(
                matches!(
                    &first.scene().primitives[indices.disc].shape,
                    Shape::AnimatedInstance(_)
                ),
                "disc must remain an animated instance"
            );
            return;
        };
        assert_eq!(disc.object_id(), scene_config.object_ids.disc);
        assert_eq!(disc.geometry_identity(), receipt.mesh_identity);
        let SharedGeometry::Mesh(mesh) = disc.geometry() else {
            assert!(
                matches!(disc.geometry(), SharedGeometry::Mesh(_)),
                "uncertified axisymmetric chart must be rendered through the derived mesh"
            );
            return;
        };
        assert_eq!(mesh.vertices.len(), receipt.vertex_count);
        assert_eq!(mesh.triangles.len(), receipt.triangle_count);
        assert!(
            mesh.vertices.iter().any(|vertex| {
                let radius = vertex[0].hypot(vertex[1]);
                radius > 0.037 && radius < 0.038 && vertex[2] > 0.002 && vertex[2] < 0.003
            }),
            "the preview must retain chord samples on the true circular rim fillet"
        );

        let first_pose = first
            .pose_at(0.0, EventEvaluationSide::RightLimit)
            .expect("first source pose");
        assert_eq!(first_pose.disc.translation_m(), [0.0, 0.0, 0.045]);
        assert_eq!(first_pose.base_plate.translation_m(), [0.0, 0.0, 0.0]);
        let last_pose = first
            .pose_at(END_TIME_S, EventEvaluationSide::RightLimit)
            .expect("last source pose");
        assert_eq!(
            last_pose.disc.translation_m(),
            [0.005 * END_TIME_S, 0.0, 0.045]
        );
        assert_eq!(
            last_pose.base_plate.translation_m(),
            [0.0, 0.0, 0.0002 * END_TIME_S]
        );

        let mut changed_style = scene_config.clone();
        changed_style.disc_material = EulerMaterialStyle::Ggx {
            linear_rgb: [0.61, 0.74, 0.78],
            alpha: 0.12,
        };
        let restyled = EulerCinematicScene::try_build(&artifact, &specimen, changed_style, cx)
            .expect("restyled scene");
        assert_ne!(first.scene_identity(), restyled.scene_identity());
        assert_eq!(
            first.preview_mesh_receipt().mesh_identity,
            restyled.preview_mesh_receipt().mesh_identity
        );

        let mut changed_glass = scene_config;
        changed_glass.plate_material = EulerMaterialStyle::Dielectric {
            glass: DielectricGlass::representative_borosilicate(),
            surface: DielectricSurface::SMOOTH,
        };
        let reglassed = EulerCinematicScene::try_build(&artifact, &specimen, changed_glass, cx)
            .expect("alternate glass scene");
        assert_ne!(first.scene_identity(), reglassed.scene_identity());
        assert_eq!(
            first.preview_mesh_receipt().mesh_identity,
            reglassed.preview_mesh_receipt().mesh_identity
        );

        let irrelevant_camera = PhysicalCamera::try_legacy_compatible(
            Point3::new(10.0, 0.0, 0.0),
            GeomVec3::new(1.0, 0.0, 0.0),
            GeomVec3::new(0.0, 0.0, 1.0),
            0.5,
            1.0,
            Aperture::try_circular(0.0).expect("pinhole"),
        )
        .expect("deliberately irrelevant camera");
        let prior_shot = CameraShot::try_new(
            0x4555_4c45_525f_5052,
            -2.0,
            -1.0,
            vec![CameraKeyframe::try_new(-2.0, irrelevant_camera).expect("prior keyframe")],
        )
        .expect("prior shot");
        let current_shot = camera().shots()[0].clone();
        let mut extra_history = config();
        extra_history.camera =
            AnimatedCamera::try_new(vec![prior_shot, current_shot]).expect("shot history");
        EulerCinematicScene::try_build(&artifact, &specimen, extra_history, cx)
            .expect("camera gaps and poses outside the trajectory horizon are irrelevant");
    });
}

#[test]
fn g3_asset_units_aliasing_discontinuities_and_cancellation_fail_closed() {
    let (specimen, correct_artifact) = with_cx(false, |cx| {
        let specimen = specimen(cx);
        let artifact = artifact(&specimen, None, false, Vec::new(), cx);
        (specimen, artifact)
    });
    with_cx(false, |cx| {
        let wrong_asset = artifact(
            &specimen,
            Some(identity("wrong-chart")),
            false,
            Vec::new(),
            cx,
        );
        assert!(matches!(
            EulerCinematicScene::try_build(&wrong_asset, &specimen, config(), cx),
            Err(EulerSceneError::AssetIdentityMismatch(
                "specimen_chart_identity"
            ))
        ));

        let mut wrong_units = config();
        wrong_units.length_unit = EulerSceneLengthUnit::Millimetres;
        assert!(matches!(
            EulerCinematicScene::try_build(&correct_artifact, &specimen, wrong_units, cx),
            Err(EulerSceneError::UnsupportedLengthUnit)
        ));

        let mut clipped_depth = config();
        clipped_depth.camera_far_m = 0.05;
        assert!(matches!(
            EulerCinematicScene::try_build(&correct_artifact, &specimen, clipped_depth, cx),
            Err(EulerSceneError::CameraDepthRange)
        ));

        let mut aliased = config();
        aliased.maximum_angular_step_rad = 0.01;
        assert!(matches!(
            EulerCinematicScene::try_build(&correct_artifact, &specimen, aliased, cx),
            Err(EulerSceneError::AngularSamplingAmbiguous { interval: 0, .. })
        ));

        let valid_contact = trajectory(&specimen, None, true, None, cx);
        let contact_metadata = valid_contact.metadata().clone();
        let mut floating_inputs = valid_contact
            .samples()
            .iter()
            .map(|sample| sample.input().clone())
            .collect::<Vec<_>>();
        floating_inputs[1].center_of_mass_world_m.z += 0.01;
        let floating_contact = EulerRenderTrajectoryArtifact::try_from_trajectory(
            identity("floating-contact-campaign"),
            RenderTrajectory::try_new(contact_metadata, floating_inputs)
                .expect("trajectory metadata alone admits the inconsistent chart contact"),
            Vec::new(),
            RenderTrajectoryCodecBudget::DEFAULT,
            cx,
        )
        .expect("floating-contact artifact");
        assert!(matches!(
            EulerCinematicScene::try_build(&floating_contact, &specimen, config(), cx),
            Err(EulerSceneError::ContactSpecimenMismatch(1))
        ));

        let discontinuous = artifact(
            &specimen,
            None,
            false,
            vec![DeclaredTimelineDiscontinuity {
                time_s: STEP_S,
                kind: DeclaredDiscontinuityKind::ContinuationSeam,
            }],
            cx,
        );
        assert!(matches!(
            EulerCinematicScene::try_build(&discontinuous, &specimen, config(), cx),
            Err(EulerSceneError::DeclaredDiscontinuityUnsupported)
        ));
    });
    with_cx(true, |cx| {
        assert!(matches!(
            EulerCinematicScene::try_build(&correct_artifact, &specimen, config(), cx),
            Err(EulerSceneError::Cancelled)
        ));
    });
}

#[test]
fn g3_subject_bounds_include_interior_cubic_hermite_translation_overshoot() {
    with_cx(false, |cx| {
        let specimen = specimen(cx);
        let artifact = overshoot_artifact(&specimen, cx);
        let eye = Point3::new(0.0, -2.0, 0.6);
        let target = Point3::new(0.0, 0.0, 0.025);
        let physical = PhysicalCamera::try_look_at(
            eye,
            target,
            GeomVec3::new(0.0, 0.0, 1.0),
            CameraProjection::try_half_tangent(0.48).expect("projection"),
            target.delta_from(eye).norm(),
            Aperture::try_circular(0.0).expect("pinhole"),
        )
        .expect("overshoot camera");
        let mut scene_config = config();
        scene_config.camera =
            AnimatedCamera::try_static(0x4555_4c45_525f_4f56, 0.0, END_TIME_S, physical)
                .expect("overshoot shot");
        scene_config.camera_far_m = 5.0;
        let scene = EulerCinematicScene::try_build(&artifact, &specimen, scene_config, cx)
            .expect("overshoot-bounded scene");
        let midpoint = scene
            .pose_at(0.5 * STEP_S, EventEvaluationSide::RightLimit)
            .expect("Hermite midpoint")
            .disc
            .translation_m();
        assert!(
            midpoint[0] > 0.2,
            "fixture must actually overshoot endpoints"
        );
        assert!(scene.subject_bounds_m().contains(Point3::new(
            midpoint[0],
            midpoint[1],
            midpoint[2]
        )));
        assert!(
            scene.subject_bounds_m().max.x > 0.3,
            "swept bounds must include Bezier controls, not only source positions"
        );
    });
}

#[test]
fn g3_event_crossing_stays_as_two_explicit_weighted_films() {
    with_cx(false, |cx| {
        let specimen = specimen(cx);
        let artifact = artifact(&specimen, None, true, Vec::new(), cx);
        let beauty_config = config();
        let scene = EulerCinematicScene::try_build(&artifact, &specimen, beauty_config.clone(), cx)
            .expect("event-bearing scene");
        let mut debug_config = beauty_config;
        debug_config.debug_overlay = EulerDebugOverlay::ContactMarker {
            sample_index: 1,
            radius_m: 0.002,
        };
        let debug_scene = EulerCinematicScene::try_build(&artifact, &specimen, debug_config, cx)
            .expect("separate debug layer");
        assert_eq!(scene.scene_identity(), debug_scene.scene_identity());
        assert_eq!(scene.scene().primitives.len(), 4);
        assert_eq!(debug_scene.scene().primitives.len(), 4);
        let debug_receipt = debug_scene
            .debug_layer_receipt()
            .expect("contact marker receipt");
        assert_eq!(debug_receipt.source_sample_index, 1);
        assert_eq!(debug_receipt.radius_m, 0.002);
        assert_eq!(
            debug_scene
                .static_scene_with_debug_at(END_TIME_S, CutSide::After, cx)
                .expect("explicit diagnostic scene")
                .primitives
                .len(),
            5
        );
        let request = EulerFrameRequest {
            frame_time_s: STEP_S,
            exposure_duration_s: END_TIME_S,
            event_policy: ExposureEventPolicy::Subdivide,
            ..frame_request(ExposureEventPolicy::Subdivide)
        };
        let prepared = scene.prepare_frame(request).expect("event partition");
        assert_eq!(prepared.segments().len(), 2);
        assert_eq!(prepared.segments()[0].shutter().open_s(), 0.0);
        assert_eq!(prepared.segments()[0].shutter().close_s(), STEP_S);
        assert_eq!(prepared.segments()[1].shutter().open_s(), STEP_S);
        assert_eq!(prepared.segments()[1].shutter().close_s(), END_TIME_S);
        assert_eq!(prepared.segments()[0].duration_weight(), 0.5);
        assert_eq!(prepared.segments()[1].duration_weight(), 0.5);

        let settings = euler_scene_smoke_settings(2, 2);
        let before = scene
            .render_segment(&prepared, 0, &settings, cx)
            .expect("pre-event film");
        let after = scene
            .render_segment(&prepared, 1, &settings, cx)
            .expect("post-event film");
        assert_eq!(before.spp_done, settings.spp);
        assert_eq!(after.spp_done, settings.spp);
        assert!(matches!(
            scene.render_frame(request, &settings, cx),
            Err(EulerSceneError::ExposureNeedsComposition { segment_count: 2 })
        ));

        let refused = EulerFrameRequest {
            event_policy: ExposureEventPolicy::Refuse,
            ..request
        };
        assert!(scene.prepare_frame(refused).is_err());

        let beauty_request = EulerFrameRequest {
            frame_time_s: END_TIME_S,
            exposure_duration_s: 0.0,
            event_policy: ExposureEventPolicy::Refuse,
            ..request
        };
        assert_eq!(
            scene
                .render_frame(beauty_request, &settings, cx)
                .expect("beauty without configured marker")
                .xyz,
            debug_scene
                .render_frame(beauty_request, &settings, cx)
                .expect("beauty excludes configured marker")
                .xyz
        );
    });
}

#[test]
fn e2e_zero_width_frame_traces_real_pixels_and_round_trips_exr_deterministically() {
    with_cx(false, |cx| {
        let specimen = specimen(cx);
        let artifact = artifact(&specimen, None, false, Vec::new(), cx);
        let scene = EulerCinematicScene::try_build(&artifact, &specimen, config(), cx)
            .expect("renderable scene");
        let settings = euler_scene_smoke_settings(12, 8);
        let request = frame_request(ExposureEventPolicy::Refuse);

        let cinematic = scene
            .render_frame(request, &settings, cx)
            .expect("zero-width cinematic render");
        let static_film = scene
            .render_static_at(STEP_S, CutSide::After, &settings, cx)
            .expect("materialized static render");
        for (pixel_index, (cinematic_pixel, static_pixel)) in
            cinematic.xyz.iter().zip(&static_film.xyz).enumerate()
        {
            for channel in 0..3 {
                let scale = cinematic_pixel[channel]
                    .abs()
                    .max(static_pixel[channel].abs())
                    .max(1.0);
                assert_close(
                    cinematic_pixel[channel],
                    static_pixel[channel],
                    2.0e-15 * scale,
                    &format!("static/cinematic pixel {pixel_index} channel {channel}"),
                );
            }
        }
        assert_eq!(cinematic.spp_done, settings.spp);

        let first_bytes = scene
            .render_frame_exr(request, &settings, cx)
            .expect("first EXR");
        let replay_bytes = scene
            .render_frame_exr(request, &settings, cx)
            .expect("replayed EXR");
        assert_eq!(first_bytes, replay_bytes);

        let decoded = fs_img::read_exr(&first_bytes).expect("decode generated EXR");
        assert_eq!((decoded.width, decoded.height), (12, 8));
        assert_eq!(
            decoded
                .channels
                .iter()
                .map(|channel| channel.name.as_str())
                .collect::<Vec<_>>(),
            vec!["B", "G", "R"]
        );
        assert!(decoded.channels.iter().all(|channel| {
            channel.data.len() == 12 * 8 && channel.data.iter().all(|value| value.is_finite())
        }));
        assert!(
            decoded.channels.iter().any(|channel| {
                channel.data.iter().any(|value| value.abs() > 1.0e-8)
                    && channel
                        .data
                        .windows(2)
                        .any(|pair| pair[0].to_bits() != pair[1].to_bits())
            }),
            "the E2E must contain illuminated, spatially varying traced pixels"
        );
        assert_eq!(
            fs_img::write_exr(decoded.width, decoded.height, &decoded.channels)
                .expect("re-encode decoded EXR"),
            first_bytes
        );
    });
}
