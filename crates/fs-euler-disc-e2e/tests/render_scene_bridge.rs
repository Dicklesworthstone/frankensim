//! Focused G0/G3/E2E coverage for the animated Euler scene bridge.

#![cfg(feature = "cinematic-render")]

use std::collections::BTreeMap;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::{ContentHash, hash_domain};
use fs_euler_disc_e2e::coupled_runner::{ChannelOwnership, ContactTransitionKind};
use fs_euler_disc_e2e::profile_contact_geometry;
use fs_euler_disc_e2e::render_scene_bridge::{
    EULER_STUDIO_ENVIRONMENT_HEIGHT, EULER_STUDIO_ENVIRONMENT_WIDTH, EulerCinematicScene,
    EulerDebugOverlay, EulerEnvironmentStyle, EulerFrameRequest, EulerMaterialStyle,
    EulerSceneError, EulerSceneLengthUnit, EulerStudioEnvironmentSpec, EulerSupportSurfaceSpec,
    EulerTessellationConfig, MAX_EULER_STUDIO_ENVIRONMENT_RADIANCE_SCALE,
    euler_scene_smoke_settings,
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
use fs_exec::{Budget, CancelGate, Cx, ExecMode, RunId, StreamKey};
use fs_geom::{Point3, Vec3 as GeomVec3};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_render::camera::{
    AnimatedCamera, Aperture, CameraKeyframe, CameraProjection, CameraShot, CutSide, PhysicalCamera,
};
use fs_render::charts::TriMesh;
use fs_render::conductor::{ConductorOptics, ConductorProvenance, ConductorSurface};
use fs_render::dielectric::{CauchyIor, DielectricGlass, DielectricSurface, GlassProvenance};
use fs_render::instances::SharedGeometry;
use fs_render::motion::{ShutterConvention, ShutterDistribution};
use fs_render::tracer::{
    AdaptiveFilm, AdaptiveSamplingConfig, Film, Material, RenderExecutionConfig,
    RenderExecutionError, RenderWorkerPool, Shape,
};
use fs_rep_frep::SquatDiscEdgeTreatment;

const END_TIME_S: f64 = 0.02;
const STEP_S: f64 = 0.01;

fn assert_close(actual: f64, expected: f64, tolerance: f64, context: &str) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{context}: actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:.3e}"
    );
}

fn assert_film_bits_eq(actual: &Film, expected: &Film, context: &str) {
    assert_eq!(
        (
            actual.width,
            actual.height,
            actual.spp_done,
            actual.time_mode
        ),
        (
            expected.width,
            expected.height,
            expected.spp_done,
            expected.time_mode,
        ),
        "{context}: film metadata differed"
    );
    assert_eq!(actual.xyz.len(), expected.xyz.len(), "{context}");
    for (pixel, (actual, expected)) in actual.xyz.iter().zip(&expected.xyz).enumerate() {
        for channel in 0..3 {
            assert_eq!(
                actual[channel].to_bits(),
                expected[channel].to_bits(),
                "{context}: pixel={pixel} channel={channel} actual={:#018x} expected={:#018x}",
                actual[channel].to_bits(),
                expected[channel].to_bits(),
            );
        }
    }
}

fn assert_adaptive_film_bits_eq(actual: &AdaptiveFilm, expected: &AdaptiveFilm, context: &str) {
    assert_eq!(
        (
            actual.width(),
            actual.height(),
            actual.maximum_samples(),
            actual.policy(),
            actual.sampler(),
            actual.stream_seed(),
            actual.semantics_version(),
            actual.time_mode(),
        ),
        (
            expected.width(),
            expected.height(),
            expected.maximum_samples(),
            expected.policy(),
            expected.sampler(),
            expected.stream_seed(),
            expected.semantics_version(),
            expected.time_mode(),
        ),
        "{context}: adaptive film identity differed"
    );
    assert_eq!(
        actual.sample_counts(),
        expected.sample_counts(),
        "{context}"
    );
    assert_eq!(actual.decisions(), expected.decisions(), "{context}");
    for pixel in 0..actual.xyz_sums().len() {
        for (label, actual, expected) in [
            ("sum", actual.xyz_sums()[pixel], expected.xyz_sums()[pixel]),
            (
                "mean",
                actual.running_means_xyz()[pixel],
                expected.running_means_xyz()[pixel],
            ),
            ("m2", actual.m2_xyz()[pixel], expected.m2_xyz()[pixel]),
        ] {
            for channel in 0..3 {
                assert_eq!(
                    actual[channel].to_bits(),
                    expected[channel].to_bits(),
                    "{context}: {label} pixel={pixel} channel={channel} actual={:#018x} expected={:#018x}",
                    actual[channel].to_bits(),
                    expected[channel].to_bits(),
                );
            }
        }
    }
}

fn assert_closed_outward_mesh(mesh: &TriMesh, expected_min_m: [f64; 3], expected_max_m: [f64; 3]) {
    assert!(!mesh.vertices.is_empty());
    assert!(!mesh.triangles.is_empty());
    let mut actual_min_m = [f64::INFINITY; 3];
    let mut actual_max_m = [f64::NEG_INFINITY; 3];
    for vertex in &mesh.vertices {
        assert!(vertex.iter().all(|coordinate| coordinate.is_finite()));
        for axis in 0..3 {
            actual_min_m[axis] = actual_min_m[axis].min(vertex[axis]);
            actual_max_m[axis] = actual_max_m[axis].max(vertex[axis]);
        }
    }
    for axis in 0..3 {
        assert_eq!(actual_min_m[axis].to_bits(), expected_min_m[axis].to_bits());
        assert_eq!(actual_max_m[axis].to_bits(), expected_max_m[axis].to_bits());
    }

    let center_m = [
        0.5 * (expected_min_m[0] + expected_max_m[0]),
        0.5 * (expected_min_m[1] + expected_max_m[1]),
        0.5 * (expected_min_m[2] + expected_max_m[2]),
    ];
    // Production bevel meshes deliberately duplicate vertices between face
    // patches. Canonicalize edges by exact coordinates rather than raw vertex
    // index so this checks the represented surface for cracks, not one storage
    // layout. All generators here derive shared coordinates through identical
    // arithmetic, so exact bit keys are intentional.
    let vertex_key = |vertex: [f64; 3]| {
        vertex.map(|coordinate| {
            if coordinate == 0.0 {
                0.0_f64.to_bits()
            } else {
                coordinate.to_bits()
            }
        })
    };
    let mut undirected_edge_incidence = BTreeMap::<([u64; 3], [u64; 3]), usize>::new();
    for (triangle_index, triangle) in mesh.triangles.iter().enumerate() {
        for index in triangle {
            assert!(
                (*index as usize) < mesh.vertices.len(),
                "triangle {triangle_index} has out-of-range vertex {index}"
            );
        }
        for [from_index, to_index] in [
            [triangle[0], triangle[1]],
            [triangle[1], triangle[2]],
            [triangle[2], triangle[0]],
        ] {
            let from = vertex_key(mesh.vertices[from_index as usize]);
            let to = vertex_key(mesh.vertices[to_index as usize]);
            let edge = if from < to { (from, to) } else { (to, from) };
            *undirected_edge_incidence.entry(edge).or_default() += 1;
        }

        let a = mesh.vertices[triangle[0] as usize];
        let b = mesh.vertices[triangle[1] as usize];
        let c = mesh.vertices[triangle[2] as usize];
        let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let ac = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
        let normal = [
            ab[1] * ac[2] - ab[2] * ac[1],
            ab[2] * ac[0] - ab[0] * ac[2],
            ab[0] * ac[1] - ab[1] * ac[0],
        ];
        let centroid_from_center = [
            (a[0] + b[0] + c[0]) / 3.0 - center_m[0],
            (a[1] + b[1] + c[1]) / 3.0 - center_m[1],
            (a[2] + b[2] + c[2]) / 3.0 - center_m[2],
        ];
        let outward_dot = normal[0] * centroid_from_center[0]
            + normal[1] * centroid_from_center[1]
            + normal[2] * centroid_from_center[2];
        assert!(
            outward_dot > 0.0,
            "triangle {triangle_index} is degenerate or not outward-wound: dot={outward_dot:.17e}"
        );
    }
    assert!(
        undirected_edge_incidence.values().all(|count| *count == 2),
        "every represented surface edge must have exactly two incident triangles: {undirected_edge_incidence:?}"
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
        // This rendering-only fixture declares every mechanical channel
        // unavailable. Contact branch/event geometry is independent of force
        // magnitude, so retain the required exact zero instead of inventing an
        // interval-mean normal-force authority for these scene tests.
        interval_normal_force_n: 0.0,
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
        assert_eq!(scene_config.environment, EulerEnvironmentStyle::None);
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
        assert_eq!(indices.support_surface, None);
        assert_eq!(indices.light, 3);
        assert_eq!(indices.spin_fiducial, None);
        assert_eq!(first.debug_layer_receipt(), None);
        assert_eq!(first.scene().primitives.len(), 4);
        assert_eq!(
            first.scene().lights.len(),
            1,
            "the v1 Euler cinematic rig must retain exactly one rectangular emitter"
        );
        assert_eq!(first.scene().lights[0].prim, indices.light);
        assert!(
            first.scene().environment.is_none(),
            "the v1 Euler cinematic rig must not silently acquire an environment emitter"
        );
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

        let Shape::AnimatedInstance(plate) = &first.scene().primitives[indices.base_plate].shape
        else {
            assert!(
                matches!(
                    &first.scene().primitives[indices.base_plate].shape,
                    Shape::AnimatedInstance(_)
                ),
                "the dielectric base plate must remain an animated instance"
            );
            return;
        };
        assert_eq!(plate.object_id(), scene_config.object_ids.base_plate);
        let SharedGeometry::Mesh(plate_mesh) = plate.geometry() else {
            assert!(
                matches!(plate.geometry(), SharedGeometry::Mesh(_)),
                "the dielectric base plate must remain an explicit closed mesh"
            );
            return;
        };
        assert_closed_outward_mesh(
            plate_mesh,
            [
                -0.5 * scene_config.base.plate_width_m,
                -0.5 * scene_config.base.plate_depth_m,
                -scene_config.base.plate_thickness_m,
            ],
            [
                0.5 * scene_config.base.plate_width_m,
                0.5 * scene_config.base.plate_depth_m,
                0.0,
            ],
        );

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
fn g0_opt_in_studio_environment_is_identity_bound_and_propagates_to_static_scenes() {
    with_cx(false, |cx| {
        let specimen = specimen(cx);
        let artifact = artifact(&specimen, None, false, Vec::new(), cx);
        let reference_config = config();
        let reference =
            EulerCinematicScene::try_build(&artifact, &specimen, reference_config.clone(), cx)
                .expect("reference black-world scene");
        assert!(reference.scene().environment.is_none());

        let mut studio_config = reference_config;
        studio_config.environment =
            EulerEnvironmentStyle::StudioGradient(EulerStudioEnvironmentSpec::SOFT_NEUTRAL);
        let studio =
            EulerCinematicScene::try_build(&artifact, &specimen, studio_config.clone(), cx)
                .expect("native studio environment");
        let replay =
            EulerCinematicScene::try_build(&artifact, &specimen, studio_config.clone(), cx)
                .expect("deterministic studio replay");
        assert_ne!(reference.scene_identity(), studio.scene_identity());
        assert_ne!(
            reference.source_configuration_identity(),
            studio.source_configuration_identity()
        );
        assert_eq!(studio.scene_identity(), replay.scene_identity());
        assert_eq!(
            studio.source_configuration_identity(),
            replay.source_configuration_identity()
        );

        let environment = studio
            .scene()
            .environment
            .as_ref()
            .expect("configured environment must be present");
        assert_eq!(environment.width(), EULER_STUDIO_ENVIRONMENT_WIDTH);
        assert_eq!(environment.height(), EULER_STUDIO_ENVIRONMENT_HEIGHT);
        assert!(!environment.is_black());
        assert_eq!(
            environment.semantic_hash(),
            replay
                .scene()
                .environment
                .as_ref()
                .expect("replayed environment")
                .semantic_hash()
        );
        let static_scene = studio
            .static_scene_at(STEP_S, CutSide::After, cx)
            .expect("static scene with studio environment");
        assert_eq!(
            static_scene
                .environment
                .as_ref()
                .expect("static environment")
                .semantic_hash(),
            environment.semantic_hash(),
            "static materialization must retain the exact animated-scene environment"
        );

        for (label, changed_spec) in [
            (
                "overhead",
                EulerStudioEnvironmentSpec {
                    overhead_linear_rgb: [0.21, 0.26, 0.36],
                    ..EulerStudioEnvironmentSpec::SOFT_NEUTRAL
                },
            ),
            (
                "horizon",
                EulerStudioEnvironmentSpec {
                    horizon_linear_rgb: [0.60, 0.49, 0.34],
                    ..EulerStudioEnvironmentSpec::SOFT_NEUTRAL
                },
            ),
            (
                "floor",
                EulerStudioEnvironmentSpec {
                    floor_linear_rgb: [0.025, 0.018, 0.015],
                    ..EulerStudioEnvironmentSpec::SOFT_NEUTRAL
                },
            ),
            (
                "scale",
                EulerStudioEnvironmentSpec {
                    radiance_scale: 0.80,
                    ..EulerStudioEnvironmentSpec::SOFT_NEUTRAL
                },
            ),
        ] {
            let mut changed_config = studio_config.clone();
            changed_config.environment = EulerEnvironmentStyle::StudioGradient(changed_spec);
            let changed = EulerCinematicScene::try_build(&artifact, &specimen, changed_config, cx)
                .unwrap_or_else(|error| panic!("changed {label} environment: {error}"));
            assert_ne!(
                studio.source_configuration_identity(),
                changed.source_configuration_identity(),
                "{label} must affect configuration identity"
            );
            assert_ne!(
                studio.scene_identity(),
                changed.scene_identity(),
                "{label} must affect scene identity"
            );
        }
    });
}

#[test]
fn g0_studio_environment_rejects_nonfinite_negative_and_unbounded_inputs() {
    with_cx(false, |cx| {
        let specimen = specimen(cx);
        let artifact = artifact(&specimen, None, false, Vec::new(), cx);
        for (label, invalid_spec) in [
            (
                "nonfinite overhead",
                EulerStudioEnvironmentSpec {
                    overhead_linear_rgb: [f64::NAN, 0.26, 0.36],
                    ..EulerStudioEnvironmentSpec::SOFT_NEUTRAL
                },
            ),
            (
                "negative horizon",
                EulerStudioEnvironmentSpec {
                    horizon_linear_rgb: [-0.01, 0.48, 0.34],
                    ..EulerStudioEnvironmentSpec::SOFT_NEUTRAL
                },
            ),
            (
                "unbounded floor",
                EulerStudioEnvironmentSpec {
                    floor_linear_rgb: [0.025, 0.018, 1.01],
                    ..EulerStudioEnvironmentSpec::SOFT_NEUTRAL
                },
            ),
            (
                "zero radiance",
                EulerStudioEnvironmentSpec {
                    radiance_scale: 0.0,
                    ..EulerStudioEnvironmentSpec::SOFT_NEUTRAL
                },
            ),
            (
                "unbounded radiance",
                EulerStudioEnvironmentSpec {
                    radiance_scale: MAX_EULER_STUDIO_ENVIRONMENT_RADIANCE_SCALE + 0.01,
                    ..EulerStudioEnvironmentSpec::SOFT_NEUTRAL
                },
            ),
        ] {
            let mut invalid_config = config();
            invalid_config.environment = EulerEnvironmentStyle::StudioGradient(invalid_spec);
            let result = EulerCinematicScene::try_build(&artifact, &specimen, invalid_config, cx);
            assert!(
                matches!(result, Err(EulerSceneError::InvalidConfig("environment"))),
                "{label} must refuse before allocating or rendering"
            );
        }
    });
}

#[test]
fn g0_optional_support_surface_is_physical_identity_bound_geometry() {
    with_cx(false, |cx| {
        let specimen = specimen(cx);
        let artifact = artifact(&specimen, None, false, Vec::new(), cx);
        let reference_config = config();
        let reference =
            EulerCinematicScene::try_build(&artifact, &specimen, reference_config.clone(), cx)
                .expect("reference scene");

        let mut supported_config = reference_config;
        supported_config.support_surface = Some(EulerSupportSurfaceSpec {
            width_m: 1.2,
            depth_m: 1.0,
            thickness_m: 0.02,
            gap_below_housing_m: 0.001,
            material: EulerMaterialStyle::Lambertian {
                linear_rgb: [0.07, 0.055, 0.045],
            },
        });
        let supported =
            EulerCinematicScene::try_build(&artifact, &specimen, supported_config.clone(), cx)
                .expect("scene on finite support surface");
        let replay =
            EulerCinematicScene::try_build(&artifact, &specimen, supported_config.clone(), cx)
                .expect("support-surface replay");

        let indices = supported.primitive_indices();
        let support_index = indices.support_surface.expect("support primitive");
        assert_eq!(support_index, 3);
        assert_eq!(indices.light, 4);
        assert_eq!(supported.scene().primitives.len(), 5);
        assert_ne!(reference.scene_identity(), supported.scene_identity());
        assert_ne!(
            reference.source_configuration_identity(),
            supported.source_configuration_identity()
        );
        assert_eq!(supported.scene_identity(), replay.scene_identity());
        assert_eq!(
            supported.source_configuration_identity(),
            replay.source_configuration_identity()
        );

        let Shape::Instance(support) = &supported.scene().primitives[support_index].shape else {
            panic!("support surface must be a static geometry instance");
        };
        assert_eq!(
            support.object_id(),
            supported_config.object_ids.support_surface
        );
        let SharedGeometry::Mesh(mesh) = support.geometry() else {
            panic!("support surface must use explicit closed geometry");
        };
        let housing_bottom = -supported_config.base.plate_thickness_m
            - supported_config.base.housing_gap_m
            - supported_config.base.housing_height_m;
        let support_top = housing_bottom - 0.001;
        assert_closed_outward_mesh(
            mesh,
            [-0.6, -0.5, support_top - 0.02],
            [0.6, 0.5, support_top],
        );
        assert!(matches!(
            supported.scene().primitives[support_index].material,
            Material::Lambertian { .. }
        ));

        let mut duplicate_id = supported_config.clone();
        duplicate_id.object_ids.support_surface = duplicate_id.object_ids.housing;
        assert!(matches!(
            EulerCinematicScene::try_build(&artifact, &specimen, duplicate_id, cx),
            Err(EulerSceneError::InvalidConfig("support surface object_id"))
        ));

        let admitted_support = supported_config.support_surface.unwrap();
        let invalid_supports = [
            (
                EulerSupportSurfaceSpec {
                    width_m: 0.0,
                    ..admitted_support
                },
                "support width_m",
            ),
            (
                EulerSupportSurfaceSpec {
                    depth_m: f64::NAN,
                    ..admitted_support
                },
                "support depth_m",
            ),
            (
                EulerSupportSurfaceSpec {
                    thickness_m: -0.02,
                    ..admitted_support
                },
                "support thickness_m",
            ),
            (
                EulerSupportSurfaceSpec {
                    gap_below_housing_m: -0.001,
                    ..admitted_support
                },
                "support gap_below_housing_m",
            ),
            (
                EulerSupportSurfaceSpec {
                    material: EulerMaterialStyle::Lambertian {
                        linear_rgb: [0.07, f64::INFINITY, 0.045],
                    },
                    ..admitted_support
                },
                "support material",
            ),
        ];
        for (support_surface, expected_field) in invalid_supports {
            let mut invalid_config = supported_config.clone();
            invalid_config.support_surface = Some(support_surface);
            assert!(matches!(
                EulerCinematicScene::try_build(&artifact, &specimen, invalid_config, cx),
                Err(EulerSceneError::InvalidConfig(field)) if field == expected_field
            ));
        }

        for support_surface in [
            EulerSupportSurfaceSpec {
                width_m: f64::from_bits(1),
                ..admitted_support
            },
            EulerSupportSurfaceSpec {
                gap_below_housing_m: f64::MAX,
                ..admitted_support
            },
            EulerSupportSurfaceSpec {
                thickness_m: f64::from_bits(1),
                ..admitted_support
            },
        ] {
            let mut invalid_config = supported_config.clone();
            invalid_config.support_surface = Some(support_surface);
            assert!(matches!(
                EulerCinematicScene::try_build(&artifact, &specimen, invalid_config, cx),
                Err(EulerSceneError::InvalidConfig("support derived bounds"))
            ));
        }
    });
}

#[test]
fn g0_spin_fiducial_is_disc_local_identity_bound_and_default_off() {
    with_cx(false, |cx| {
        let specimen = specimen(cx);
        let artifact = artifact(&specimen, None, false, Vec::new(), cx);
        let reference_config = config();
        let reference =
            EulerCinematicScene::try_build(&artifact, &specimen, reference_config.clone(), cx)
                .expect("reference scene");
        assert_eq!(reference.primitive_indices().spin_fiducial, None);
        assert_eq!(reference.primitive_indices().support_surface, None);
        assert_eq!(reference.scene().primitives.len(), 4);
        assert_eq!(reference.primitive_indices().light, 3);

        let mut marked_config = reference_config;
        marked_config.show_spin_fiducial = true;
        let marked =
            EulerCinematicScene::try_build(&artifact, &specimen, marked_config.clone(), cx)
                .expect("scene with spin fiducial");
        let replay = EulerCinematicScene::try_build(&artifact, &specimen, marked_config, cx)
            .expect("replayed scene with spin fiducial");
        let indices = marked.primitive_indices();
        let marker_index = indices.spin_fiducial.expect("enabled marker primitive");
        assert_eq!(marker_index, 3);
        assert_eq!(indices.light, 4);
        assert_eq!(marked.scene().primitives.len(), 5);
        assert_ne!(reference.scene_identity(), marked.scene_identity());
        assert_ne!(
            reference.source_configuration_identity(),
            marked.source_configuration_identity()
        );
        assert_eq!(marked.scene_identity(), replay.scene_identity());
        assert_eq!(
            marked.source_configuration_identity(),
            replay.source_configuration_identity()
        );

        let Shape::AnimatedInstance(disc) = &marked.scene().primitives[indices.disc].shape else {
            panic!("disc must retain animated placement");
        };
        let Shape::AnimatedInstance(marker) = &marked.scene().primitives[marker_index].shape else {
            panic!("spin fiducial must share an animated placement");
        };
        assert_eq!(marker.object_id(), marked_config_object_id());
        assert_eq!(
            marker.trajectory().keyframes(),
            disc.trajectory().keyframes(),
            "the fiducial must use the exact disc trajectory rather than a sampled approximation"
        );
        for time_s in [0.0, 0.005, END_TIME_S] {
            assert_eq!(
                marker
                    .trajectory()
                    .evaluate(time_s)
                    .expect("marker trajectory"),
                disc.trajectory().evaluate(time_s).expect("disc trajectory"),
                "the fiducial must remain rigidly co-moving with the disc"
            );
        }
        let SharedGeometry::Mesh(marker_mesh) = marker.geometry() else {
            panic!("spin fiducial must be an explicit mesh");
        };
        assert!(
            marker_mesh
                .vertices
                .iter()
                .all(|vertex| vertex[2] > marked.preview_mesh_receipt().local_bounds_m.max.z),
            "positive local lift must prevent fiducial z-fighting with the derived disc mesh"
        );
        assert!(matches!(
            marked.scene().primitives[marker_index].material,
            Material::Lambertian { .. }
        ));
    });
}

fn marked_config_object_id() -> u64 {
    fs_euler_disc_e2e::render_scene_bridge::EULER_SPIN_FIDUCIAL_OBJECT_ID
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
        let mut rough_contact_inputs = valid_contact
            .samples()
            .iter()
            .map(|sample| sample.input().clone())
            .collect::<Vec<_>>();
        rough_contact_inputs[1].signed_gap_m = -2.0e-6;
        let rough_contact = EulerRenderTrajectoryArtifact::try_from_trajectory(
            identity("rough-contact-campaign"),
            RenderTrajectory::try_new(contact_metadata.clone(), rough_contact_inputs)
                .expect("resolved rough gap is valid independently of the smooth chart gap"),
            Vec::new(),
            RenderTrajectoryCodecBudget::DEFAULT,
            cx,
        )
        .expect("rough-contact artifact");
        EulerCinematicScene::try_build(&rough_contact, &specimen, config(), cx)
            .expect("rough-contact normal-law gap must not be equated to smooth chart height");

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
fn e2e_representative_conductors_bind_identity_and_render_distinctly() {
    with_cx(false, |cx| {
        let specimen = specimen(cx);
        let artifact = artifact(&specimen, None, false, Vec::new(), cx);
        let surface = ConductorSurface::try_rough(0.12).expect("valid conductor roughness");
        let tungsten_optics = ConductorOptics::representative_tungsten();
        let stainless_optics = ConductorOptics::representative_stainless_steel();
        assert_ne!(
            tungsten_optics.provenance(),
            stainless_optics.provenance(),
            "the representative presets must retain distinct provenance"
        );

        let mut tungsten_config = config();
        tungsten_config.light.linear_rgb = [1.0, 1.0, 1.0];
        tungsten_config.disc_material = EulerMaterialStyle::Conductor {
            optics: tungsten_optics,
            surface,
        };
        let mut stainless_config = tungsten_config.clone();
        stainless_config.disc_material = EulerMaterialStyle::Conductor {
            optics: stainless_optics,
            surface,
        };

        let tungsten =
            EulerCinematicScene::try_build(&artifact, &specimen, tungsten_config.clone(), cx)
                .expect("representative tungsten scene");
        let tungsten_replay =
            EulerCinematicScene::try_build(&artifact, &specimen, tungsten_config, cx)
                .expect("replayed representative tungsten scene");
        let stainless =
            EulerCinematicScene::try_build(&artifact, &specimen, stainless_config.clone(), cx)
                .expect("representative stainless-steel scene");
        let stainless_replay =
            EulerCinematicScene::try_build(&artifact, &specimen, stainless_config, cx)
                .expect("replayed representative stainless-steel scene");

        assert_eq!(tungsten.scene_identity(), tungsten_replay.scene_identity());
        assert_eq!(
            stainless.scene_identity(),
            stainless_replay.scene_identity()
        );
        assert_eq!(
            tungsten.source_configuration_identity(),
            tungsten_replay.source_configuration_identity()
        );
        assert_eq!(
            stainless.source_configuration_identity(),
            stainless_replay.source_configuration_identity()
        );
        assert_ne!(
            tungsten.source_configuration_identity(),
            stainless.source_configuration_identity(),
            "changing only the conductor optical preset must invalidate configuration identity"
        );
        assert_ne!(
            tungsten.scene_identity(),
            stainless.scene_identity(),
            "changing only the conductor optical preset must invalidate scene identity"
        );
        assert_eq!(
            tungsten.source_trajectory_identity(),
            stainless.source_trajectory_identity(),
            "material look development must not alter the mechanics trajectory"
        );
        assert_eq!(
            tungsten.preview_mesh_receipt(),
            stainless.preview_mesh_receipt(),
            "material look development must not alter the derived disc geometry"
        );

        let tungsten_indices = tungsten.primitive_indices();
        let stainless_indices = stainless.primitive_indices();
        let tungsten_material = tungsten.scene().primitives[tungsten_indices.disc].material;
        let stainless_material = stainless.scene().primitives[stainless_indices.disc].material;
        let tungsten_provenance: ConductorProvenance = match tungsten_material {
            Material::Conductor {
                optics,
                surface: found_surface,
            } => {
                assert_eq!(optics, tungsten_optics);
                assert_eq!(found_surface, surface);
                optics.provenance()
            }
            other => panic!("Euler disc did not retain tungsten conductor material: {other:?}"),
        };
        let stainless_provenance: ConductorProvenance = match stainless_material {
            Material::Conductor {
                optics,
                surface: found_surface,
            } => {
                assert_eq!(optics, stainless_optics);
                assert_eq!(found_surface, surface);
                optics.provenance()
            }
            other => {
                panic!("Euler disc did not retain stainless-steel conductor material: {other:?}")
            }
        };
        assert_eq!(tungsten_provenance, tungsten_optics.provenance());
        assert_eq!(stainless_provenance, stainless_optics.provenance());
        assert_ne!(
            tungsten_material.content_identity(),
            stainless_material.content_identity(),
            "different optical tables must produce different tracer-material identities"
        );
        for (label, tungsten_index, stainless_index) in [
            (
                "base plate",
                tungsten_indices.base_plate,
                stainless_indices.base_plate,
            ),
            (
                "housing",
                tungsten_indices.housing,
                stainless_indices.housing,
            ),
            ("light", tungsten_indices.light, stainless_indices.light),
        ] {
            assert_eq!(
                tungsten.scene().primitives[tungsten_index]
                    .material
                    .content_identity(),
                stainless.scene().primitives[stainless_index]
                    .material
                    .content_identity(),
                "the {label} material changed in a conductor-only comparison"
            );
        }

        let settings = euler_scene_smoke_settings(12, 8);
        let request = frame_request(ExposureEventPolicy::Refuse);
        let tungsten_film = tungsten
            .render_frame(request, &settings, cx)
            .expect("production tungsten render");
        let tungsten_film_replay = tungsten_replay
            .render_frame(request, &settings, cx)
            .expect("replayed production tungsten render");
        let stainless_film = stainless
            .render_frame(request, &settings, cx)
            .expect("production stainless-steel render");
        let stainless_film_replay = stainless_replay
            .render_frame(request, &settings, cx)
            .expect("replayed production stainless-steel render");
        assert_film_bits_eq(
            &tungsten_film,
            &tungsten_film_replay,
            "representative tungsten render was not exactly replayable",
        );
        assert_film_bits_eq(
            &stainless_film,
            &stainless_film_replay,
            "representative stainless-steel render was not exactly replayable",
        );

        for (label, film) in [
            ("tungsten", &tungsten_film),
            ("stainless steel", &stainless_film),
        ] {
            let mut absolute_energy = 0.0_f64;
            for (pixel_index, pixel) in film.xyz.iter().enumerate() {
                for (channel, value) in pixel.iter().copied().enumerate() {
                    assert!(
                        value.is_finite(),
                        "{label} render produced a non-finite value at pixel={pixel_index} channel={channel}: {value}"
                    );
                    absolute_energy += value.abs();
                }
            }
            assert!(
                absolute_energy > 1.0e-8,
                "{label} conductor render was black; absolute_energy={absolute_energy:.17e}"
            );
        }

        let mut differing_pixels = 0_usize;
        let mut total_absolute_delta = 0.0_f64;
        let mut comparison_energy = 0.0_f64;
        for (tungsten_pixel, stainless_pixel) in tungsten_film.xyz.iter().zip(&stainless_film.xyz) {
            let pixel_delta = tungsten_pixel
                .iter()
                .zip(stainless_pixel)
                .map(|(tungsten, stainless)| (tungsten - stainless).abs())
                .sum::<f64>();
            let pixel_energy = tungsten_pixel
                .iter()
                .zip(stainless_pixel)
                .map(|(tungsten, stainless)| tungsten.abs().max(stainless.abs()))
                .sum::<f64>();
            assert!(pixel_delta.is_finite() && pixel_energy.is_finite());
            if tungsten_pixel
                .iter()
                .zip(stainless_pixel)
                .any(|(tungsten, stainless)| tungsten.to_bits() != stainless.to_bits())
            {
                differing_pixels += 1;
            }
            total_absolute_delta += pixel_delta;
            comparison_energy += pixel_energy;
        }
        let normalized_l1_delta = total_absolute_delta / comparison_energy.max(f64::MIN_POSITIVE);
        assert!(
            differing_pixels > 0 && normalized_l1_delta > 1.0e-8,
            "representative tungsten and stainless steel must render distinguishably under the same neutral light; differing_pixels={differing_pixels}, normalized_l1_delta={normalized_l1_delta:.17e}"
        );
    });
}

#[test]
fn e2e_zero_width_frame_traces_real_pixels_and_round_trips_exr_deterministically() {
    with_cx(false, |cx| {
        let specimen = specimen(cx);
        let artifact = artifact(&specimen, None, false, Vec::new(), cx);
        let reference_config = config();
        let scene =
            EulerCinematicScene::try_build(&artifact, &specimen, reference_config.clone(), cx)
                .expect("renderable scene");
        let settings = euler_scene_smoke_settings(12, 8);
        let request = frame_request(ExposureEventPolicy::Refuse);

        let cinematic = scene
            .render_frame(request, &settings, cx)
            .expect("zero-width cinematic render");
        let execution =
            RenderExecutionConfig::try_new(5, 3, 4, 64 << 20, RunId(0x4555_4c45_525f_5041))
                .expect("valid Euler cinematic tile policy");
        let parallel = scene
            .render_frame_with_execution(request, &settings, &execution, cx)
            .expect("tile-parallel zero-width cinematic render");
        assert_film_bits_eq(
            &parallel.film,
            &cinematic,
            "Euler cinematic parallel render differed from serial oracle",
        );
        assert_eq!(parallel.report.requested_workers, 4);
        assert_eq!(parallel.report.workers, 4);
        assert_eq!(parallel.report.layout.tile_count(), 9);
        assert_eq!(parallel.report.executor.declared_run, execution.run_id());
        assert_eq!(
            (
                parallel.report.executor.completed,
                parallel.report.executor.total,
            ),
            (9, 9)
        );
        assert_eq!(parallel.report.memory.used_bytes, 0);

        let next_request = EulerFrameRequest {
            frame_time_s: END_TIME_S,
            ..request
        };
        let next_serial = scene
            .render_frame(next_request, &settings, cx)
            .expect("next zero-width cinematic render");
        let next_execution =
            RenderExecutionConfig::try_new(4, 5, 4, 64 << 20, RunId(0x4555_4c45_525f_5042))
                .expect("valid second Euler cinematic tile policy");
        let adaptive =
            AdaptiveSamplingConfig::try_new(2, 2, 0.0, 0.0, 0.0).expect("adaptive policy");
        let adaptive_reference = scene
            .render_frame_adaptive_with_execution(request, &settings, adaptive, &execution, cx)
            .expect("one-shot adaptive Euler cinematic render");
        let pending = scene
            .begin_frame_render(request, settings, execution.clone(), cx)
            .expect("admitted opaque Euler frame job");
        let adaptive_pending = scene
            .begin_frame_adaptive_render(request, settings, adaptive, execution.clone(), cx)
            .expect("admitted opaque adaptive Euler frame job");
        let worker_pool = RenderWorkerPool::new(&execution, cx.mode(), 0x4555_4c45_525f_4352);
        let (parked_first, parked_next, owned, parked_adaptive, owned_adaptive) = worker_pool
            .with_parked_crew_local(|parked| {
                let first = scene
                    .render_frame_with_parked_scope(parked, request, &settings, &execution, cx)
                    .expect("first render on parked Euler cinematic crew");
                let next = scene
                    .render_frame_with_parked_scope(
                        parked,
                        next_request,
                        &settings,
                        &next_execution,
                        cx,
                    )
                    .expect("second render on same parked Euler cinematic crew");
                let owned = pending
                    .resume_on_parked(parked, cx)
                    .expect("opaque Euler job on parked crew");
                let parked_adaptive = scene
                    .render_frame_adaptive_with_parked_scope(
                        parked, request, &settings, adaptive, &execution, cx,
                    )
                    .expect("adaptive Euler frame on parked crew");
                let owned_adaptive = adaptive_pending
                    .resume_on_parked(parked, cx)
                    .expect("opaque adaptive Euler job on parked crew");
                (first, next, owned, parked_adaptive, owned_adaptive)
            });
        assert_film_bits_eq(
            &parked_first.film,
            &cinematic,
            "first parked Euler frame differed from serial oracle",
        );
        assert_film_bits_eq(
            &parked_next.film,
            &next_serial,
            "second parked Euler frame differed from serial oracle",
        );
        assert_film_bits_eq(
            &owned.film,
            &cinematic,
            "owned single-film Euler frame differed from serial oracle",
        );
        assert!(owned.report.progress_state_bytes > 0);
        assert_eq!(owned.report.memory.used_bytes, 0);
        assert_adaptive_film_bits_eq(
            &parked_adaptive.film,
            &adaptive_reference.film,
            "parked adaptive Euler frame differed from one-shot reference",
        );
        assert_adaptive_film_bits_eq(
            &owned_adaptive.film,
            &adaptive_reference.film,
            "opaque adaptive Euler frame differed from one-shot reference",
        );
        assert!(owned_adaptive.report.progress_state_bytes > 0);
        assert_eq!(owned_adaptive.report.memory.used_bytes, 0);
        assert_eq!(
            adaptive_reference.report.executor.declared_run,
            execution.run_id()
        );
        assert_eq!(
            parked_adaptive.report.executor.declared_run,
            execution.run_id()
        );
        assert_eq!(
            owned_adaptive.report.executor.declared_run,
            execution.run_id()
        );
        assert_eq!(
            parked_first.report.executor.declared_run,
            execution.run_id()
        );
        assert_eq!(
            parked_next.report.executor.declared_run,
            next_execution.run_id()
        );
        assert_ne!(
            parked_first.report.executor.declared_run,
            parked_next.report.executor.declared_run
        );
        assert_eq!(
            (
                parked_first.report.executor.completed,
                parked_first.report.executor.total,
            ),
            (9, 9)
        );
        assert_eq!(
            (
                parked_next.report.executor.completed,
                parked_next.report.executor.total,
            ),
            (6, 6)
        );

        let one_byte_execution =
            RenderExecutionConfig::try_new(5, 3, 4, 1, RunId(0x4555_4c45_525f_4f4d))
                .expect("one-byte limit is structurally valid");
        assert!(matches!(
            scene.render_frame_with_execution(request, &settings, &one_byte_execution, cx),
            Err(EulerSceneError::RenderExecution(
                RenderExecutionError::Memory(_)
            ))
        ));
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

        let crown = DielectricGlass::representative_crown();
        let [_, crown_b_um2, crown_c_um4] = crown.ior().coefficients();
        let low_index_ior = CauchyIor::try_new(1.0, crown_b_um2, crown_c_um4)
            .expect("same-dispersion-class low-index control");
        assert!(crown.ior().is_dispersive() && low_index_ior.is_dispersive());
        let mut low_index_config = reference_config;
        low_index_config.plate_material = EulerMaterialStyle::Dielectric {
            glass: DielectricGlass::new(low_index_ior, crown.absorption(), GlassProvenance::Custom),
            surface: DielectricSurface::POLISHED,
        };
        let low_index_scene =
            EulerCinematicScene::try_build(&artifact, &specimen, low_index_config, cx)
                .expect("matched-absorption low-index scene");
        let low_index_film = low_index_scene
            .render_frame(request, &settings, cx)
            .expect("low-index comparison render");

        let mut materially_changed_pixels = 0_usize;
        let mut total_absolute_delta = 0.0_f64;
        let mut comparison_energy = 0.0_f64;
        for (glass_pixel, low_index_pixel) in cinematic.xyz.iter().zip(&low_index_film.xyz) {
            let pixel_delta = glass_pixel
                .iter()
                .zip(low_index_pixel)
                .map(|(glass, low_index)| (glass - low_index).abs())
                .sum::<f64>();
            let pixel_energy = glass_pixel
                .iter()
                .zip(low_index_pixel)
                .map(|(glass, low_index)| glass.abs().max(low_index.abs()))
                .sum::<f64>();
            assert!(pixel_delta.is_finite() && pixel_energy.is_finite());
            if pixel_delta > 1.0e-9 * pixel_energy.max(1.0) {
                materially_changed_pixels += 1;
            }
            total_absolute_delta += pixel_delta;
            comparison_energy += pixel_energy;
        }
        let normalized_l1_delta = total_absolute_delta / comparison_energy.max(f64::MIN_POSITIVE);
        assert!(
            materially_changed_pixels >= 4 && normalized_l1_delta > 1.0e-5,
            "changing only the plate's phase-index magnitude within the same dispersive estimator class must materially change multiple traced pixels; changed_pixels={materially_changed_pixels}, normalized_l1_delta={normalized_l1_delta:.17e}"
        );

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
