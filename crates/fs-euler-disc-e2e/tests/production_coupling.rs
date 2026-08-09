//! Synthetic G0 composition fixture; it establishes plumbing, not calibration.

use fs_contact::normal_patch::{
    ApplicabilityInput, ApplicabilityLimits, InputUncertainty, NormalPatchEmbedState,
    NormalPatchReceipt,
};
use fs_couple::StableId;
use fs_euler_disc_e2e::air::{
    AirFilmDiscretization, AirFilmTransactionState, AirFilmWorkOwnership, AirVec3,
    ContactExclusion, PrescribedPlaneBase, TILTED_DISC_GAS_FILM_ADAPTER_ID, TiltedDiscAirFilmInput,
    TiltedDiscKinematics,
};
use fs_euler_disc_e2e::base_response::{ReducedBasePort, ReducedBasePortIdentity};
use fs_euler_disc_e2e::external_air::{
    EulerDiscBodyFrame, EulerDiscExteriorGeometry, EulerDiscExteriorState, EulerExternalAirInput,
    EulerExternalAirWorkState, ExteriorAirPressure, ExternalAirDomain, ExternalAirIdentity,
};
use fs_euler_disc_e2e::normal_contact::{
    EulerNormalContactInput, EulerNormalContactOutcome, EulerNormalGeometry,
    NORMAL_CONTACT_ADAPTER_ID, NormalContactIdentity, NormalContactIntegrationRegime,
    NormalMaterialInterface, NormalRateResponse, evaluate_normal_contact,
};
use fs_euler_disc_e2e::patch_kinematics::{
    CurvatureMetadata, MovingOneModeBaseState, MovingOneModePatchBridgeInput,
    MovingOneModePatchKinematicsInput, OrderedSurfacePair, PatchGeometryMetadata,
    PatchKinematicThresholds, ProfileSupportKinematics, SurfaceOrder, TangentGaugeInput,
    compute_moving_one_mode_patch_kinematics,
};
use fs_euler_disc_e2e::production_coupling::{
    GasChannelReceipt, GasChannelState, GasChannelStepInput, ProductionCouplingError,
    ProductionCouplingIdentity, ProductionCouplingModel, ProductionCouplingStepInput,
    ProductionEventTrajectoryTermination, ProductionOpenFlightStepInput,
    ProductionSurfaceExcitationError, ProductionSurfaceExcitationStepInput,
    ProductionSurfaceTraceStepInput, ProductionTrajectoryBranch, ProductionTrajectoryStepReceipt,
    SmoothContactTrajectoryTermination,
};
use fs_euler_disc_e2e::rolling_contact::{
    ROLLING_CONTACT_ADAPTER_ID, RollingContactIdentity, RollingContactInput, RollingContactState,
};
use fs_euler_disc_e2e::specimen::DiscProfileSpec;
use fs_euler_disc_e2e::tangential_contact::{
    EulerTangentialContactAdapter, TangentialContactLane, TangentialContactRequest,
};
use fs_euler_disc_e2e::{
    BaseGeometryScope, BaseResponseInput, ContactLoadScope, LevelSupportInput, MovingContactLoad,
    RenderNormalForceSampling, RenderSampleDisposition, RenderTrajectory, profile_contact_geometry,
    state_at_profile_ground_contact,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_flux::Vec3 as FluxVec3;
use fs_flux::{
    ApplicabilityEnvelope, ClosedRange, ContributionFamily, CorrelationIdentity,
    CorrelationUncertainty, FormDrag, GasFilmApplicability, GasFilmBoundaryTopology, GasFilmBudget,
    GasFilmInputAuthority, GasFilmUncertainty, GasPropertyCard, IsothermalIdealGas,
    ReducedAeroComponents, ReducedAeroModel, RoughnessPolicy, SlipPolicy, SurfaceRoughness,
};
use fs_mbd::{Gravity, MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_rep_frep::{AxisymmetricCurvatureError, AxisymmetricMassProperties, SquatDiscEdgeTreatment};
use fs_solid::{
    AssemblyBudget, DampingModel, ShellIdentity, ShellMaterial, ShellNode, ShellPlate, ShellSupport,
};
use fs_tribo::partial_slip::{
    GeneralizedWorkOwnership, NormalPatchAuthority, NormalPatchView, PARTIAL_SLIP_MODEL_ID,
    PartialSlipInterface, PartialSlipLaw, PartialSlipParameters,
};
use fs_tribo::rolling_loss::{
    CoulombContourCard, LEINE_STYLE_CONTOUR_LAW_ID, RollingLossApplicability, RollingLossChannel,
    RollingLossLaw, RollingWorkOwnership,
};
use fs_tribo::surface_excitation::{
    PeriodicHarmonicSurface, PeriodicSurfaceHarmonic, SurfaceTraceMotion,
};
use fs_tribo::{ApplicabilityRange, InputAuthority, InterfaceMedium, InterfaceSystemRef};

fn id(value: &str) -> StableId {
    StableId::new(value).expect("synthetic fixture identity")
}

fn with_cx<R>(operation: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        operation(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x5052_4f46_494c_45,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn profile_mbd_mass(profile: AxisymmetricMassProperties) -> MassProperties {
    MassProperties::new(
        profile.mass,
        Vec3::ZERO,
        Vec3::new(
            profile.principal_inertia.transverse,
            profile.principal_inertia.transverse,
            profile.principal_inertia.axial,
        ),
    )
    .expect("profile mass is a valid centroidal rigid-body mass")
}

fn mass() -> MassProperties {
    MassProperties::new(0.2, Vec3::ZERO, Vec3::new(1.0e-4, 1.0e-4, 1.0e-4))
        .expect("synthetic fixture mass")
}

fn disc_state() -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(Vec3::new(0.0, 0.0, -1.0e-6), UnitQuaternion::IDENTITY).expect("pose"),
        Vec3::new(0.2, 0.0, 0.0),
        Vec3::ZERO,
    )
    .expect("state")
}

fn thin_gap_disc_state() -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(Vec3::new(0.0, 0.0, 1.0e-3), UnitQuaternion::IDENTITY).expect("thin pose"),
        Vec3::new(0.2, 0.0, 0.0),
        Vec3::ZERO,
    )
    .expect("thin state")
}

fn thresholds() -> PatchKinematicThresholds {
    PatchKinematicThresholds {
        threshold_identity: id("synthetic/production-thresholds"),
        tie_break_identity: id("synthetic/production-ties"),
        separation_gap_m: 1.0e-3,
        touching_gap_m: 1.0e-4,
        support_point_coincidence_tolerance_m: 1.0e-12,
        tangent_counterpart_coincidence_tolerance_m: 1.0e-12,
        approach_speed_m_per_s: 1.0e-3,
        impact_candidate_speed_m_per_s: 1.0,
        stationary_tangent_speed_m_per_s: 1.0e-12,
        minimum_reference_rolling_speed_m_per_s: 1.0e-6,
        gauge_degeneracy_norm: 1.0e-12,
    }
}

fn patch_input() -> MovingOneModePatchKinematicsInput {
    MovingOneModePatchKinematicsInput {
        bridge: MovingOneModePatchBridgeInput {
            profile_support: ProfileSupportKinematics {
                disc_arm_world_m: Vec3::ZERO,
                disc_point_world_m: Vec3::new(0.0, 0.0, -1.0e-6),
                gap_m: -1.0e-6,
                source_feature: 0,
                support_authority: fs_rep_frep::AxisymmetricSupportAuthority::Estimate,
            },
            disc_state: disc_state(),
            disc_mass_properties: mass(),
            base_mode: MovingOneModeBaseState {
                undeformed_contact_point_world_m: Vec3::ZERO,
                vertical_displacement_m: 0.0,
                vertical_velocity_m_per_s: 0.0,
            },
            normal_world: Vec3::new(0.0, 0.0, 1.0),
            tangent_gauge: TangentGaugeInput {
                reference_world: Vec3::new(1.0, 0.0, 0.0),
                rotation_rad: 0.0,
            },
            thresholds: thresholds(),
        },
        surfaces: OrderedSurfacePair::try_new(id("disc"), id("base"), SurfaceOrder::DiscThenBase)
            .expect("ordered pair"),
        patch: PatchGeometryMetadata {
            patch_identity: id("synthetic/rim-patch"),
            source_feature: 0,
            gap_uncertainty_m: 0.0,
            curvature: CurvatureMetadata::Known {
                curvature_identity: id("synthetic/relative-gap-curvature"),
                authority: InputAuthority::SyntheticFixture,
                first_principal_m_inverse: 25.0,
                second_principal_m_inverse: 25.0,
                uncertainty_m_inverse: 1.0e-6,
            },
        },
        tangent_effort_probe_world_n: None,
    }
}

fn normal_material() -> NormalMaterialInterface {
    NormalMaterialInterface {
        material_card_id: "synthetic/steel-card".into(),
        model_id: "synthetic/hertz-model".into(),
        source_id: "synthetic/material-source".into(),
        interface: InterfaceSystemRef::new(
            "disc->base",
            "synthetic/normal-history",
            "synthetic/interface",
            InputAuthority::SyntheticFixture,
            InterfaceMedium::Dry,
        )
        .expect("interface"),
        reduced_modulus_pa: 2.0e9,
        rate_response: NormalRateResponse::ElasticHertz,
        applicability: ApplicabilityInput {
            half_space_depth_m: 1.0,
            layer_thickness_m: 1.0,
            yield_strength_pa: 1.0e12,
            characteristic_rate_m_per_s: 1.0,
            temperature_k: 293.15,
            adhesion_energy_j_per_m2: 0.0,
        },
        limits: ApplicabilityLimits {
            max_patch_to_radius: 0.2,
            max_strain: 0.1,
            max_patch_to_depth: 0.2,
            max_patch_to_layer: 0.2,
            max_pressure_to_yield: 0.2,
            max_rate_ratio: 1.0,
            min_temperature_k: 200.0,
            max_temperature_k: 400.0,
        },
        uncertainty: InputUncertainty {
            radius_relative: 0.001,
            modulus_relative: 0.02,
            load_relative: 0.03,
        },
    }
}

fn normal_input(
    kinematics: fs_euler_disc_e2e::patch_kinematics::PatchKinematics,
) -> EulerNormalContactInput {
    EulerNormalContactInput {
        identity: NormalContactIdentity {
            case_id: "synthetic/production-case".into(),
            adapter_id: NORMAL_CONTACT_ADAPTER_ID.into(),
            solver_id: "synthetic/smooth-solver".into(),
            contact_id: "synthetic/rim-contact".into(),
            sample_id: "synthetic/sample-1".into(),
        },
        kinematics,
        material: normal_material(),
        geometry: EulerNormalGeometry::SpherePlane,
        integration_regime: NormalContactIntegrationRegime::SmoothQuasistatic,
        state: NormalPatchEmbedState::new(0.0, 1.0).expect("normal state"),
        time_s: 0.0,
        iteration: 1,
        step_s: 1.0e-6,
        converged: true,
    }
}

fn base_port() -> ReducedBasePort {
    let support = ShellSupport {
        node_indices: [0, 1, 2],
        normal: [0.0, 0.0, 1.0],
    };
    ReducedBasePort::build(
        ReducedBasePortIdentity {
            model_id: "synthetic/base-model".into(),
            configuration_id: "synthetic/base-config".into(),
        },
        BaseResponseInput {
            plate: ShellPlate {
                nodes: vec![
                    ShellNode {
                        position_m: [-0.1, -0.08, 0.0],
                    },
                    ShellNode {
                        position_m: [0.1, -0.08, 0.0],
                    },
                    ShellNode {
                        position_m: [0.0, 0.12, 0.0],
                    },
                    ShellNode {
                        position_m: [0.0, 0.0, 0.0],
                    },
                ],
                triangles: vec![[0, 1, 3], [1, 2, 3], [2, 0, 3]],
                thickness_m: 0.004,
                material: ShellMaterial {
                    youngs_modulus_pa: 70.0e9,
                    poisson_ratio: 0.33,
                    density_kg_m3: 2700.0,
                },
                identity: ShellIdentity {
                    model_id: "synthetic/base-plate".into(),
                    source_id: "synthetic/base-source".into(),
                    state_id: "initial".into(),
                },
                support: Some(support),
                damping: DampingModel::None,
                budget: AssemblyBudget::default(),
            },
            level_support: LevelSupportInput {
                support,
                level_normal: [0.0, 0.0, 1.0],
                maximum_tilt_rad: 1.0e-6,
            },
            geometry_scope: BaseGeometryScope::FlatSinglePatch,
            contact_scope: ContactLoadScope::NodalNormalLoad,
            load: MovingContactLoad {
                start_node: 3,
                end_node: 3,
                normal_force_n: 1.0,
            },
            initial_modal_displacement_m: 0.0,
            initial_modal_velocity_m_per_s: 0.0,
            timestep_s: 1.0e-6,
            steps: 1,
        },
        4,
    )
    .expect("synthetic base port")
}

fn tangential_law() -> PartialSlipLaw {
    PartialSlipLaw::new(
        PARTIAL_SLIP_MODEL_ID,
        "synthetic/friction-card",
        PartialSlipParameters {
            static_mu: 0.8,
            kinetic_mu: 0.4,
            tangential_stiffness_n_per_m: 10_000.0,
            torsional_stiffness_nm_per_rad: 100.0,
            torsional_capacity_factor: 0.5,
            partial_slip_onset_fraction: 0.5,
            partial_slip_hardening_fraction: 0.4,
        },
    )
    .expect("synthetic partial-slip law")
}

fn tangential_adapter() -> EulerTangentialContactAdapter {
    EulerTangentialContactAdapter::new(
        "synthetic/tangential-adapter",
        "synthetic/tangential-source",
        TangentialContactLane::PartialSlip {
            law: tangential_law(),
        },
    )
    .expect("adapter")
}

fn aero() -> EulerExternalAirInput {
    let range = |a, b| ClosedRange::try_new(a, b).expect("range");
    let model = ReducedAeroModel::try_new(
        CorrelationIdentity::try_new("synthetic/exterior-a", "v1", "synthetic/aero-source")
            .expect("identity"),
        ApplicabilityEnvelope {
            translational_reynolds: range(0.0, 1.0e9),
            rotational_reynolds: range(0.0, 1.0e9),
            relative_roughness: range(0.0, 1.0),
            maximum_tip_mach: 1.0,
        },
        CorrelationUncertainty {
            source_id: "synthetic/aero-uncertainty".into(),
            coefficient_relative_half_width: 0.1,
        },
        ReducedAeroComponents {
            form_drag: Some(FormDrag { coefficient: 1.0 }),
            rotational_skin_friction: None,
            edge_flow: None,
            orientation_rate_damping: None,
        },
        &[ContributionFamily::TranslationalFormDrag],
    )
    .expect("aero model");
    EulerExternalAirInput {
        domain: ExternalAirDomain::ExteriorFreeGas,
        identity: ExternalAirIdentity {
            case_id: "synthetic/production-case".into(),
            world_frame_id: "synthetic/world".into(),
            body_frame_id: "synthetic/body".into(),
            geometry_source_id: "synthetic/geometry".into(),
            state_source_id: "synthetic/state".into(),
            domain_source_id: "synthetic/exterior".into(),
        },
        geometry: EulerDiscExteriorGeometry {
            radius_m: 0.04,
            exterior_thickness_m: 0.003,
        },
        state: EulerDiscExteriorState {
            center_world_m: FluxVec3::new(0.0, 0.0, -1.0e-6),
            // Deliberately stale: production coupling must replace all of this
            // pose/rate record from the accepted fs-mbd checkpoint.
            center_velocity_world_m_per_s: FluxVec3::ZERO,
            angular_velocity_world_rad_per_s: FluxVec3::ZERO,
            body_frame: EulerDiscBodyFrame {
                x_world: FluxVec3::new(1.0, 0.0, 0.0),
                z_world: FluxVec3::new(0.0, 0.0, 1.0),
            },
        },
        gas: GasPropertyCard {
            source_id: "synthetic/gas".into(),
            density_kg_per_m3: Some(1.2),
            dynamic_viscosity_pa_s: Some(1.8e-5),
            speed_of_sound_m_per_s: Some(340.0),
            velocity_world_m_per_s: fs_flux::Vec3::ZERO,
        },
        pressure: ExteriorAirPressure {
            absolute_pressure_pa: 101_325.0,
            source_id: "synthetic/pressure".into(),
        },
        exterior_roughness: SurfaceRoughness {
            source_id: "synthetic/roughness".into(),
            height_m: 1.0e-5,
        },
        alternatives: vec![model],
    }
}

/// Synthetic, applicable thin-gap card. Its stale disc record is intentional:
/// production coupling must replace it from the fs-mbd checkpoint.
fn thin_gap_air() -> TiltedDiscAirFilmInput {
    let pressure = 101_325.0;
    TiltedDiscAirFilmInput {
        identity: fs_euler_disc_e2e::air::AirFilmIdentity {
            case_id: "synthetic/production-case".into(),
            adapter_model_id: TILTED_DISC_GAS_FILM_ADAPTER_ID.into(),
            frame_id: "synthetic/world".into(),
            base_motion_id: "synthetic/prescribed-plane".into(),
            gas_species_id: "synthetic-air".into(),
            eos_id: "synthetic-isothermal-ideal-gas".into(),
            viscosity_source_id: "synthetic-viscosity".into(),
            thermal_model_id: "synthetic-isothermal".into(),
            configuration_id: "synthetic/production-config".into(),
            deterministic_seed: 17,
            authority: GasFilmInputAuthority::SyntheticFixture,
        },
        disc_radius_m: 1.0e-2,
        disc_half_thickness_m: 1.0e-4,
        disc: TiltedDiscKinematics {
            center_world_m: AirVec3::ZERO,
            normal_away_from_base_world: AirVec3::new(0.0, 0.0, 1.0),
            center_velocity_world_m_per_s: AirVec3::ZERO,
            angular_velocity_world_rad_per_s: AirVec3::ZERO,
        },
        base: PrescribedPlaneBase {
            height_m: 0.0,
            vertical_velocity_m_per_s: 0.0,
        },
        discretization: AirFilmDiscretization {
            azimuthal_sectors: 4,
            radial_cells: 4,
        },
        contact_exclusion: ContactExclusion {
            handoff_gap_m: 1.0e-6,
        },
        gas: IsothermalIdealGas {
            specific_gas_constant_j_kg_k: 287.05,
            temperature_k: 300.0,
            dynamic_viscosity_pa_s: 1.8e-5,
            declared_density_kg_m3: pressure / (287.05 * 300.0),
            declared_specific_enthalpy_j_kg: 300_000.0,
        },
        boundary: GasFilmBoundaryTopology::Sealed,
        slip_policy: SlipPolicy::NoSlipContinuum {
            source_id: "synthetic/no-slip".into(),
        },
        roughness_policy: RoughnessPolicy::ResolvedSmooth {
            source_id: "synthetic/smooth".into(),
            maximum_roughness_m: 1.0e-8,
        },
        applicability: GasFilmApplicability {
            mean_free_path_m: 65.0e-9,
            maximum_knudsen_number: 0.01,
            maximum_gap_slope: 0.2,
            speed_of_sound_m_per_s: 347.0,
            maximum_mach_number: 0.3,
        },
        uncertainty: GasFilmUncertainty {
            viscosity_relative_bound: 0.0,
            gap_relative_bound: 0.0,
            pressure_relative_bound: 0.0,
        },
        initial_absolute_pressure_pa: pressure,
        gauge_reference_absolute_pressure_pa: pressure,
        timestep_s: 1.0e-6,
        budget: GasFilmBudget {
            maximum_iterations: 2_000,
            mass_residual_tolerance_kg_m2_s: 1.0e-9,
            relaxation: 0.8,
        },
    }
}

fn rolling_law() -> RollingLossLaw {
    RollingLossLaw::CoulombContour(
        CoulombContourCard::new(
            LEINE_STYLE_CONTOUR_LAW_ID,
            "synthetic/rolling-card",
            InputAuthority::SyntheticFixture,
            0.1,
            RollingLossApplicability::new(
                ApplicabilityRange::new(250.0, 350.0).unwrap(),
                ApplicabilityRange::new(0.0, 100.0).unwrap(),
            )
            .unwrap(),
        )
        .expect("rolling law"),
    )
}

/// Rebuilds synthetic cards from the accepted state exactly as a real profile
/// resolver would, while keeping this G0 fixture deliberately uncalibrated.
fn request_for_checkpoint(
    template: &ProductionCouplingStepInput,
    checkpoint: &fs_euler_disc_e2e::production_coupling::ProductionCouplingCheckpoint,
) -> ProductionCouplingStepInput {
    let mut input = template.clone();
    let version = checkpoint.committed_version;
    let disc_center = checkpoint.disc_state.pose().position_world();
    let disc_point = disc_center.add(input.patch.bridge.profile_support.disc_arm_world_m);
    input.expected_checkpoint_version = version;
    input.time_s = version as f64 * input.duration_s;
    input.normal.iteration = version + 1;
    input.patch.bridge.profile_support.disc_point_world_m = disc_point;
    input
        .patch
        .bridge
        .base_mode
        .undeformed_contact_point_world_m = Vec3::new(disc_point.x, disc_point.y, 0.0);
    input.tangential.request_id = format!("synthetic/tangent-step-{}", version + 1);
    input.tangential.work_ownership = GeneralizedWorkOwnership::new(
        "synthetic/rim-patch",
        format!("synthetic/tangent-interval-{}", version + 1),
        "long",
        "lat",
        "spin",
    )
    .expect("synthetic tangential work ownership");
    input.rolling.ownership = RollingWorkOwnership::new(
        "synthetic/rim-patch",
        format!("synthetic/rolling-interval-{}", version + 1),
        "contour",
        RollingLossChannel::ContourDeformation,
    )
    .expect("synthetic rolling work ownership");
    match &mut input.gas_channel {
        GasChannelStepInput::ExteriorFreeGas { exchange_key, .. }
        | GasChannelStepInput::ThinGap { exchange_key, .. } => *exchange_key = version + 1,
    }
    input.base_step_id = format!("synthetic/base-step-{}", version + 1);
    input
}

#[allow(clippy::too_many_lines)] // One fixture keeps every cross-channel identity mutually consistent.
fn production_request_template() -> (
    ProductionCouplingStepInput,
    NormalPatchView,
    PartialSlipInterface,
) {
    let kinematics = compute_moving_one_mode_patch_kinematics(patch_input())
        .expect("synthetic moving-base patch");
    let normal = normal_input(kinematics.clone());
    let normal_outcome = evaluate_normal_contact(&normal).expect("synthetic normal");
    let EulerNormalContactOutcome::Active(active) = normal_outcome else {
        panic!("penetrated synthetic patch must be active");
    };
    let NormalPatchReceipt::Point(point) = active.generic.receipt else {
        panic!("synthetic fixture must retain point-normal units");
    };
    let normal_view = NormalPatchView::new(
        "synthetic/rim-patch",
        "synthetic/steel-card",
        "synthetic/material-source",
        NormalPatchAuthority::SyntheticFixture,
        point.normal_force_n,
        point.patch_radius_m,
        point.patch_radius_m,
        point.pressure.second_moment_m2,
    )
    .expect("normal view");
    let interface = PartialSlipInterface::new(
        "disc->base",
        "synthetic/tangent-history",
        "synthetic/tangent-interface",
        NormalPatchAuthority::SyntheticFixture,
    )
    .expect("tangent interface");
    let request = ProductionCouplingStepInput {
        expected_checkpoint_version: 0,
        duration_s: 1.0e-6,
        time_s: 0.0,
        patch: patch_input(),
        normal,
        surface_excitation: None,
        tangential: TangentialContactRequest {
            request_id: "synthetic/tangent-step".into(),
            expected_state_version: 0,
            patch_kinematics: kinematics,
            normal_patch: normal_view.clone(),
            interface: interface.clone(),
            work_ownership: GeneralizedWorkOwnership::new(
                "synthetic/rim-patch",
                "synthetic/interval-1",
                "long",
                "lat",
                "spin",
            )
            .expect("ownership"),
            dt_s: 1.0e-6,
        },
        rolling: RollingContactInput {
            identity: RollingContactIdentity {
                case_id: "synthetic/production-case".into(),
                adapter_id: ROLLING_CONTACT_ADAPTER_ID.into(),
                world_frame_id: "synthetic/world".into(),
                port_id: "synthetic/rolling-port".into(),
                domain_id: "synthetic/rim-domain".into(),
            },
            patch: fs_tribo::rolling_loss::RollingPatchReceipt::new(
                "synthetic/rim-patch",
                "synthetic/normal",
                "synthetic/source",
                InputAuthority::SyntheticFixture,
                1.0,
                1.0e-4,
                fs_tribo::rolling_loss::PatchCurvature::Principal {
                    first_per_m: 25.0,
                    second_per_m: 25.0,
                },
            )
            .expect("placeholder overwritten"),
            interface: InterfaceSystemRef::new(
                "disc->base",
                "synthetic/rolling-history",
                "synthetic/interface",
                InputAuthority::SyntheticFixture,
                InterfaceMedium::Dry,
            )
            .expect("rolling interface"),
            law: rolling_law(),
            state: fs_tribo::rolling_loss::RollingLossState::zero(),
            checkpoint: None,
            ownership: RollingWorkOwnership::new(
                "synthetic/rim-patch",
                "synthetic/rolling-interval-1",
                "contour",
                RollingLossChannel::ContourDeformation,
            )
            .expect("rolling ownership"),
            partial_slip_ownership: None,
            contact_arm_world_m: Vec3::ZERO,
            contour_tangent_axis_world: Vec3::new(1.0, 0.0, 0.0),
            rolling_axis_world: Vec3::new(0.0, 1.0, 0.0),
            contour_speed_mps: 0.5,
            rolling_rate_rad_s: 0.0,
            spin_rate_rad_s: 0.0,
            temperature_kelvin: 293.15,
            excitation_frequency_hz: 1.0,
            interval_s: 1.0e-6,
        },
        gas_channel: GasChannelStepInput::ExteriorFreeGas {
            input: aero(),
            selected_correlation_id: "synthetic/exterior-a".into(),
            exchange_key: 1,
        },
        base_step_id: "synthetic/base-step-1".into(),
        base_load_progress_start: 0.0,
        base_load_progress_end: 0.0,
    };
    (request, normal_view, interface)
}

#[test]
fn synthetic_g0_multistep_composes_real_adapters_and_refuses_without_mutation() {
    let (request, normal_view, interface) = production_request_template();
    let adapter = tangential_adapter();
    let model = ProductionCouplingModel {
        identity: ProductionCouplingIdentity {
            case_id: "synthetic/production-case".into(),
            configuration_id: "synthetic/production-config".into(),
            world_frame_id: "synthetic/world".into(),
        },
        disc_mass_properties: mass(),
        gravity: Gravity::ZERO,
        base_port: base_port(),
        tangential_adapter: adapter.clone(),
    };
    let checkpoint = model
        .initial_checkpoint(
            disc_state(),
            NormalPatchEmbedState::new(0.0, 1.0).expect("normal checkpoint"),
            adapter
                .initial_state(&normal_view, &interface, 4)
                .expect("tangent checkpoint"),
            RollingContactState::zero(),
            GasChannelState::ExteriorFreeGas(
                EulerExternalAirWorkState::new("synthetic/exterior-work", 4)
                    .expect("air checkpoint"),
            ),
        )
        .expect("production checkpoint");
    let open_input = ProductionOpenFlightStepInput {
        expected_checkpoint_version: 0,
        duration_s: request.duration_s,
        gas_channel: request.gas_channel.clone(),
        base_step_id: "synthetic/open-base-step-1".into(),
        base_load_progress_start: 0.0,
        base_load_progress_end: 0.0,
    };
    let (open_next, open_receipt) = model
        .step_open_flight(&checkpoint, &open_input)
        .expect("open flight composes gravity, gas, and unforced support dynamics");
    assert_eq!(open_next.committed_version, 1);
    assert_eq!(
        open_receipt
            .base
            .receipt()
            .compressive_normal_force_on_base_n,
        0.0,
        "an open branch cannot fabricate a zero-force contact receipt"
    );
    assert_eq!(
        open_receipt.total_force_world_n,
        match &open_receipt.gas_channel {
            GasChannelReceipt::ExteriorFreeGas { candidate, .. } => {
                Vec3::new(
                    candidate.world_wrench.force_world_n.x,
                    candidate.world_wrench.force_world_n.y,
                    candidate.world_wrench.force_world_n.z,
                )
            }
            GasChannelReceipt::ThinGap { .. } => panic!("fixture selected exterior gas"),
        },
        "open-flight fs-mbd force must be the exact selected gas force"
    );
    assert_ne!(open_next.disc_state, checkpoint.disc_state);
    let contact_after_open = request_for_checkpoint(&request, &open_next);
    model
        .step(&open_next, &contact_after_open)
        .expect("open flight retains contact checkpoints for a later compliant branch");

    let mut invalid_open = open_input.clone();
    invalid_open.expected_checkpoint_version = 1;
    assert!(matches!(
        model.step_open_flight(&checkpoint, &invalid_open),
        Err(ProductionCouplingError::CheckpointVersionMismatch {
            expected: 1,
            observed: 0
        })
    ));
    model
        .validate_checkpoint(&checkpoint)
        .expect("a refused open proposal cannot mutate the shared checkpoint");

    let (next, receipt) = model
        .step(&checkpoint, &request)
        .expect("synthetic composition");
    assert_eq!(next.committed_version, 1);
    assert!(receipt.estimate_only);
    assert!(receipt.total_force_world_n.is_finite());
    assert!(
        matches!(
            &receipt.gas_channel,
            GasChannelReceipt::ExteriorFreeGas { candidate, .. }
                if candidate.world_wrench.force_world_n.x < 0.0
        ),
        "exterior force must use the checkpoint's +x velocity, not the stale card state"
    );

    let texture_a = PeriodicHarmonicSurface::new(
        "synthetic/disc-track",
        "synthetic/disc-profilometer",
        InputAuthority::SyntheticFixture,
        32.0e-5,
        32,
        vec![PeriodicSurfaceHarmonic {
            cycles_per_track: 1,
            cosine_amplitude_m: 1.0e-9,
            sine_amplitude_m: 0.0,
        }],
    )
    .and_then(|surface| surface.realize())
    .expect("synthetic disc spectrum resolves to a physical trace");
    let texture_b = PeriodicHarmonicSurface::new(
        "synthetic/base-track",
        "synthetic/base-profilometer",
        InputAuthority::SyntheticFixture,
        32.0e-5,
        32,
        Vec::new(),
    )
    .and_then(|surface| surface.realize())
    .expect("an empty spectrum is an exact smooth trace");
    let normal_interface = InterfaceSystemRef::new(
        "disc->base",
        "synthetic/normal-history",
        "synthetic/interface",
        InputAuthority::SyntheticFixture,
        InterfaceMedium::Dry,
    )
    .expect("normal interface identity");
    let excitation = receipt
        .evaluate_surface_excitation(
            &normal_interface,
            SurfaceTraceMotion {
                trace: &texture_a,
                path_coordinate_m: 0.0,
                path_speed_m_per_s: 0.5,
            },
            SurfaceTraceMotion {
                trace: &texture_b,
                path_coordinate_m: 2.0e-4,
                path_speed_m_per_s: 0.0,
            },
            0.0,
            0.01,
        )
        .expect("accepted contact drives reusable surface excitation");
    let NormalPatchReceipt::Point(normal_receipt) = &receipt.normal.generic.receipt else {
        panic!("production step admitted only point normal contact")
    };
    let expected_major_axis_m = normal_receipt
        .elliptic_patch_axes
        .map_or(normal_receipt.patch_radius_m, |axes| axes.semi_major_axis_m);
    assert_eq!(
        excitation.projected_half_width_m.to_bits(),
        expected_major_axis_m.to_bits(),
        "surface forcing must consume the accepted physical patch, not a sound preset"
    );
    assert!(
        excitation.normal_force_perturbation_n > 0.0,
        "a positive cosine crest must raise the admitted point-normal force"
    );

    let mut textured_request = request.clone();
    textured_request.surface_excitation = Some(ProductionSurfaceExcitationStepInput {
        interface: normal_interface.clone(),
        surface_a: ProductionSurfaceTraceStepInput {
            trace: texture_a.clone(),
            path_coordinate_m: 0.0,
            path_speed_m_per_s: 0.5,
        },
        surface_b: ProductionSurfaceTraceStepInput {
            trace: texture_b.clone(),
            path_coordinate_m: 2.0e-4,
            path_speed_m_per_s: 0.0,
        },
        travel_angle_from_patch_major_rad: 0.0,
        maximum_linearized_height_fraction: 0.01,
    });
    let (textured_next, textured_receipt) = model
        .step(&checkpoint, &textured_request)
        .expect("surface excitation participates in the atomic mechanics step");
    let coupled_excitation = textured_receipt
        .surface_excitation
        .as_ref()
        .expect("accepted textured step publishes its physical perturbation");
    assert_eq!(coupled_excitation, &excitation);
    assert!(
        textured_receipt
            .base
            .receipt()
            .compressive_normal_force_on_base_n
            > receipt.base.receipt().compressive_normal_force_on_base_n,
        "the support must receive the same increased action/reaction load"
    );
    assert!(
        textured_receipt.total_force_world_n.z > receipt.total_force_world_n.z,
        "the topography force must enter fs-mbd rather than remain a post-hoc audio signal"
    );
    assert_ne!(
        textured_next.disc_state, next.disc_state,
        "the shared evolving mechanics state must respond to admitted topography"
    );
    let wrong_interface = InterfaceSystemRef::new(
        "base->disc",
        "synthetic/normal-history",
        "synthetic/interface",
        InputAuthority::SyntheticFixture,
        InterfaceMedium::Dry,
    )
    .expect("reversed interface fixture");
    assert!(matches!(
        receipt.evaluate_surface_excitation(
            &wrong_interface,
            SurfaceTraceMotion {
                trace: &texture_a,
                path_coordinate_m: 0.0,
                path_speed_m_per_s: 0.0,
            },
            SurfaceTraceMotion {
                trace: &texture_b,
                path_coordinate_m: 0.0,
                path_speed_m_per_s: 0.0,
            },
            0.0,
            0.01,
        ),
        Err(ProductionSurfaceExcitationError::InterfaceIdentityMismatch { .. })
    ));
    let mut wrong_surface_request = textured_request.clone();
    wrong_surface_request
        .surface_excitation
        .as_mut()
        .expect("textured request retains a channel")
        .interface = wrong_interface;
    assert!(matches!(
        model.step(&checkpoint, &wrong_surface_request),
        Err(ProductionCouplingError::SurfaceExcitation(
            ProductionSurfaceExcitationError::InterfaceIdentityMismatch { .. }
        ))
    ));
    assert_eq!(checkpoint.committed_version, 0);
    model
        .validate_checkpoint(&checkpoint)
        .expect("a rejected topography proposal cannot mutate accepted mechanics state");

    let thin_air = thin_gap_air();
    let thin_channel_state = AirFilmTransactionState::new(
        "synthetic/thin-gap-transaction",
        AirFilmWorkOwnership {
            owner_id: "synthetic/thin-gap-wall-work".into(),
        },
        &thin_air,
        4,
    )
    .expect("synthetic thin-gap state");
    let mut thin_template = request.clone();
    thin_template.patch.bridge.profile_support.disc_arm_world_m = Vec3::new(0.0, 0.0, -1.001e-3);
    thin_template.gas_channel = GasChannelStepInput::ThinGap {
        input: thin_air,
        exchange_key: 1,
    };
    let mut thin_patch = thin_template.patch.clone();
    thin_patch.bridge.disc_state = thin_gap_disc_state();
    thin_patch.bridge.profile_support.disc_point_world_m = Vec3::new(0.0, 0.0, -1.0e-6);
    let thin_kinematics =
        compute_moving_one_mode_patch_kinematics(thin_patch).expect("synthetic thin patch");
    let thin_outcome =
        evaluate_normal_contact(&normal_input(thin_kinematics)).expect("synthetic thin normal");
    let EulerNormalContactOutcome::Active(thin_active) = thin_outcome else {
        return;
    };
    let NormalPatchReceipt::Point(thin_point) = thin_active.generic.receipt else {
        return;
    };
    let thin_normal_view = NormalPatchView::new(
        "synthetic/rim-patch",
        "synthetic/steel-card",
        "synthetic/material-source",
        NormalPatchAuthority::SyntheticFixture,
        thin_point.normal_force_n,
        thin_point.patch_radius_m,
        thin_point.patch_radius_m,
        thin_point.pressure.second_moment_m2,
    )
    .expect("thin normal view");
    let thin_checkpoint = model
        .initial_checkpoint(
            thin_gap_disc_state(),
            NormalPatchEmbedState::new(0.0, 1.0).expect("thin normal checkpoint"),
            adapter
                .initial_state(&thin_normal_view, &interface, 4)
                .expect("thin tangent checkpoint"),
            RollingContactState::zero(),
            GasChannelState::ThinGap(thin_channel_state),
        )
        .expect("thin production checkpoint");
    let thin_request = request_for_checkpoint(&thin_template, &thin_checkpoint);
    let (thin_next, thin_receipt) = model
        .step(&thin_checkpoint, &thin_request)
        .expect("applicable synthetic thin-gap composition");
    assert!(matches!(
        &thin_receipt.gas_channel,
        GasChannelReceipt::ThinGap { .. }
    ));
    let GasChannelReceipt::ThinGap { proposal } = &thin_receipt.gas_channel else {
        return;
    };
    assert!(proposal.step.receipt.wrench.force_world_n.x.is_finite());
    assert!(
        proposal
            .step
            .samples
            .iter()
            .any(|sample| sample.radial_relative_velocity_m_per_s.abs() > 0.0),
        "thin-gap card's stale zero velocity must be replaced by the checkpoint's +x velocity"
    );
    assert_eq!(thin_next.committed_version, 1);
    let mut inapplicable_thin = thin_request.clone();
    let GasChannelStepInput::ThinGap { input, .. } = &mut inapplicable_thin.gas_channel else {
        return;
    };
    input.applicability.mean_free_path_m = 1.0e-5;
    let thin_before_refusal = thin_checkpoint.clone();
    assert!(matches!(
        model.step(&thin_checkpoint, &inapplicable_thin),
        Err(ProductionCouplingError::AirFilm(_))
    ));
    assert_eq!(thin_checkpoint, thin_before_refusal);

    let mut stale_version = request.clone();
    stale_version.expected_checkpoint_version = 1;
    assert!(matches!(
        model.step(&checkpoint, &stale_version),
        Err(ProductionCouplingError::CheckpointVersionMismatch {
            expected: 1,
            observed: 0
        })
    ));
    let mut forged = checkpoint.clone();
    forged.disc_state = thin_gap_disc_state();
    assert!(matches!(
        model.step(&forged, &request),
        Err(ProductionCouplingError::CheckpointIntegrityMismatch)
    ));
    assert_eq!(
        checkpoint.committed_version, 0,
        "forgery cannot mutate the source checkpoint"
    );
    let mut overlapping = request.clone();
    overlapping.rolling.ownership = RollingWorkOwnership::new(
        "synthetic/rim-patch",
        "synthetic/interval-1",
        "contour",
        RollingLossChannel::ContourDeformation,
    )
    .expect("deliberately overlapping ownership");
    assert!(matches!(
        model.step(&checkpoint, &overlapping),
        Err(ProductionCouplingError::Rolling(_))
    ));
    let whole = model.run_smooth_contact_trajectory(checkpoint.clone(), 3, |state| {
        Ok(request_for_checkpoint(&request, state))
    });
    assert_eq!(whole.accepted_steps.len(), 3, "trajectory={whole:#?}");
    assert_eq!(
        whole.termination,
        SmoothContactTrajectoryTermination::StepLimitReached {
            maximum_accepted_steps: 3
        }
    );
    assert_ne!(
        whole.last_accepted_checkpoint.disc_state,
        checkpoint.disc_state
    );

    let replay = model.run_smooth_contact_trajectory(checkpoint.clone(), 3, |state| {
        Ok(request_for_checkpoint(&request, state))
    });
    assert_eq!(whole, replay, "deterministic inputs must replay exactly");

    let prefix = model.run_smooth_contact_trajectory(checkpoint.clone(), 1, |state| {
        Ok(request_for_checkpoint(&request, state))
    });
    let resumed =
        model.run_smooth_contact_trajectory(prefix.last_accepted_checkpoint.clone(), 2, |state| {
            Ok(request_for_checkpoint(&request, state))
        });
    assert_eq!(
        resumed.last_accepted_checkpoint,
        whole.last_accepted_checkpoint
    );
    let mut resumed_receipts = prefix.accepted_steps.clone();
    resumed_receipts.extend(resumed.accepted_steps.clone());
    assert_eq!(resumed_receipts, whole.accepted_steps);

    let eventful = model.run_eventful_compliant_trajectory(checkpoint.clone(), 3, |state| {
        let mut rebuilt = request_for_checkpoint(&request, state);
        rebuilt.patch.bridge.profile_support.gap_m = match state.committed_version {
            1 => 2.0e-3,
            _ => -1.0e-6,
        };
        Ok(rebuilt)
    });
    assert_eq!(
        eventful.termination,
        ProductionEventTrajectoryTermination::StepLimitReached {
            maximum_accepted_steps: 3
        }
    );
    assert_eq!(eventful.accepted_steps.len(), 3);
    assert!(matches!(
        eventful.accepted_steps[0].receipt,
        ProductionTrajectoryStepReceipt::CompliantContact(_)
    ));
    assert!(matches!(
        eventful.accepted_steps[1].receipt,
        ProductionTrajectoryStepReceipt::OpenFlight(_)
    ));
    assert!(matches!(
        eventful.accepted_steps[2].receipt,
        ProductionTrajectoryStepReceipt::CompliantContact(_)
    ));
    assert_eq!(
        eventful
            .accepted_steps
            .iter()
            .map(|step| step.branch)
            .collect::<Vec<_>>(),
        vec![
            ProductionTrajectoryBranch::CompliantContact,
            ProductionTrajectoryBranch::OpenFlight,
            ProductionTrajectoryBranch::CompliantContact,
        ]
    );
    assert_eq!(eventful.transitions.len(), 2);
    assert_eq!(
        (eventful.transitions[0].from, eventful.transitions[0].to),
        (
            ProductionTrajectoryBranch::CompliantContact,
            ProductionTrajectoryBranch::OpenFlight,
        )
    );
    assert_eq!(
        (eventful.transitions[1].from, eventful.transitions[1].to),
        (
            ProductionTrajectoryBranch::OpenFlight,
            ProductionTrajectoryBranch::CompliantContact,
        )
    );
    for transition in &eventful.transitions {
        assert_eq!(
            transition.bracket_end_s - transition.bracket_start_s,
            request.duration_s,
            "event timing authority is exactly one declared fixed-grid bracket"
        );
    }

    let terminal_crossing =
        model.run_eventful_compliant_trajectory(checkpoint.clone(), 1, |state| {
            let mut rebuilt = request_for_checkpoint(&request, state);
            if state.committed_version == 1 {
                rebuilt.patch.bridge.profile_support.gap_m = 2.0e-3;
            }
            Ok(rebuilt)
        });
    assert_eq!(terminal_crossing.accepted_steps.len(), 1);
    assert_eq!(terminal_crossing.transitions.len(), 1);
    assert_eq!(
        (
            terminal_crossing.transitions[0].from,
            terminal_crossing.transitions[0].to,
        ),
        (
            ProductionTrajectoryBranch::CompliantContact,
            ProductionTrajectoryBranch::OpenFlight,
        ),
        "a crossing during the final step must not vanish at the step budget"
    );
    assert_eq!(
        terminal_crossing.transitions[0].bracket_end_s,
        terminal_crossing.last_accepted_checkpoint.elapsed_time_s()
    );

    let stopped = model.run_smooth_contact_trajectory(checkpoint.clone(), 3, |state| {
        let mut rebuilt = request_for_checkpoint(&request, state);
        if state.committed_version == 1 {
            rebuilt.patch.bridge.profile_support.gap_m = 2.0e-3;
        }
        Ok(rebuilt)
    });
    assert_eq!(stopped.accepted_steps.len(), 1);
    assert_eq!(
        stopped.last_accepted_checkpoint, prefix.last_accepted_checkpoint,
        "the deliberate out-of-envelope proposal must not commit"
    );
    assert!(matches!(
        stopped.termination,
        SmoothContactTrajectoryTermination::Refused {
            attempted_checkpoint_version: 1,
            error: ProductionCouplingError::UnsupportedMechanism { .. }
        }
    ));

    let mut thin_gap = request;
    let GasChannelStepInput::ExteriorFreeGas { input, .. } = &mut thin_gap.gas_channel else {
        return;
    };
    input.domain = ExternalAirDomain::ThinGap;
    assert!(matches!(
        model.step(&checkpoint, &thin_gap),
        Err(ProductionCouplingError::ExteriorAir(_))
    ));
    assert_eq!(
        checkpoint.committed_version, 0,
        "refusals cannot mutate checkpoint"
    );
}

#[test]
#[allow(clippy::too_many_lines)] // One vertical slice proves the same resolved patch reaches both laws.
fn profile_native_fillet_curvature_reaches_production_normal_and_rolling_receipts() {
    with_cx(|cx| {
        let profile = DiscProfileSpec::SolidCylinder {
            outer_radius_m: 0.038,
            thickness_m: 0.006,
            edge_treatment: SquatDiscEdgeTreatment::CircularFillet { radius: 0.001 },
        }
        .resolve(7_800.0, cx)
        .expect("resolved physical 1 mm fillet");
        let disc_mass_properties = profile_mbd_mass(profile.mass_properties);
        let orientation = UnitQuaternion::from_axis_angle(Vec3::new(0.0, 1.0, 0.0), 0.55)
            .expect("finite tilted orientation");
        let grounded = state_at_profile_ground_contact(
            &profile.chart,
            7_800.0,
            orientation,
            // The partial-slip lane requires a nonzero reference rolling
            // direction. Give this G0 composition fixture a small declared
            // tangential velocity instead of asking the law to choose one at
            // a completely stationary contact.
            Vec3::new(profile.mass_properties.mass * 0.01, 0.0, 0.0),
            Vec3::ZERO,
            cx,
        )
        .expect("profile-native ground state");
        let penetrated_pose = Pose::new(
            grounded
                .pose()
                .position_world()
                // Keep this synthetic preload inside the normal-law card's
                // declared small-patch applicability envelope for a 1 mm
                // fillet. A 1 micrometre overlap produces a/R > 0.2 here and
                // is correctly refused by the model rather than extrapolated.
                .sub(Vec3::new(0.0, 0.0, 1.0e-7)),
            grounded.pose().orientation(),
        )
        .expect("finite compliant-contact pose");
        let penetrated_state = RigidBodyState::new(
            penetrated_pose,
            grounded.linear_momentum_world(),
            grounded.angular_momentum_body(),
        )
        .expect("finite compliant-contact state");
        let adapter = tangential_adapter();
        let model = ProductionCouplingModel {
            identity: ProductionCouplingIdentity {
                case_id: "synthetic/production-case".into(),
                configuration_id: "synthetic/production-config".into(),
                world_frame_id: "synthetic/world".into(),
            },
            disc_mass_properties,
            gravity: Gravity::ZERO,
            base_port: base_port(),
            tangential_adapter: adapter.clone(),
        };
        let (mut request, _, interface) = production_request_template();
        request.normal.geometry = EulerNormalGeometry::EllipticParaboloid;
        let texture_a = PeriodicHarmonicSurface::new(
            "synthetic/profile-disc-track",
            "synthetic/profile-disc-profilometer",
            InputAuthority::SyntheticFixture,
            32.0e-5,
            32,
            vec![PeriodicSurfaceHarmonic {
                cycles_per_track: 1,
                cosine_amplitude_m: 1.0e-12,
                sine_amplitude_m: 0.0,
            }],
        )
        .and_then(|surface| surface.realize())
        .expect("tiny resolved profile texture");
        let texture_b = PeriodicHarmonicSurface::new(
            "synthetic/profile-base-track",
            "synthetic/profile-base-profilometer",
            InputAuthority::SyntheticFixture,
            32.0e-5,
            32,
            Vec::new(),
        )
        .and_then(|surface| surface.realize())
        .expect("exact smooth counterpart trace");
        let surface_excitation = ProductionSurfaceExcitationStepInput {
            interface: InterfaceSystemRef::new(
                "disc->base",
                "synthetic/profile-normal-history",
                "synthetic/profile-interface",
                InputAuthority::SyntheticFixture,
                InterfaceMedium::Dry,
            )
            .expect("profile normal interface identity"),
            surface_a: ProductionSurfaceTraceStepInput {
                trace: texture_a,
                path_coordinate_m: 0.0,
                path_speed_m_per_s: 0.5,
            },
            surface_b: ProductionSurfaceTraceStepInput {
                trace: texture_b,
                path_coordinate_m: 0.0,
                path_speed_m_per_s: 0.0,
            },
            travel_angle_from_patch_major_rad: 0.0,
            maximum_linearized_height_fraction: 0.01,
        };
        request.surface_excitation = Some(surface_excitation.clone());
        let resolved = model
            .bind_initial_horizontal_plane_axisymmetric_profile_contact(
                &mut request,
                &profile,
                penetrated_state,
                cx,
            )
            .expect("smooth profile binds production patch input");
        assert_eq!(
            resolved.support.support_source_feature,
            resolved.curvature.source_feature
        );
        assert!(
            (resolved.curvature.meridional_m_inverse - 1_000.0).abs() <= 1.0e-9,
            "1 mm circular fillet must supply its actual meridional curvature"
        );
        assert_eq!(
            request.rolling.contact_arm_world_m,
            resolved.support.contact.radius_world_m
        );
        let CurvatureMetadata::Known {
            curvature_identity,
            authority,
            first_principal_m_inverse,
            second_principal_m_inverse,
            uncertainty_m_inverse,
        } = &request.patch.patch.curvature
        else {
            panic!("smooth profile must publish curvature metadata");
        };
        assert!(
            curvature_identity
                .as_str()
                .contains("principal-curvature-v1")
        );
        assert_eq!(*authority, InputAuthority::Estimated);
        assert_eq!(
            first_principal_m_inverse.to_bits(),
            resolved.curvature.meridional_m_inverse.to_bits()
        );
        assert_eq!(
            second_principal_m_inverse.to_bits(),
            resolved.curvature.azimuthal_m_inverse.to_bits()
        );
        assert_eq!(
            uncertainty_m_inverse.to_bits(),
            resolved.curvature.uncertainty_m_inverse.to_bits()
        );

        let kinematics = compute_moving_one_mode_patch_kinematics(request.patch.clone())
            .expect("profile-native moving-base patch");
        let mut initial_normal = request.normal.clone();
        initial_normal.kinematics = kinematics;
        let EulerNormalContactOutcome::Active(active) =
            evaluate_normal_contact(&initial_normal).expect("elliptic profile normal response")
        else {
            panic!("penetrated profile patch must have active normal response");
        };
        let NormalPatchReceipt::Point(point) = &active.generic.receipt else {
            panic!("elliptic profile contact must retain point-resultant units");
        };
        let (longitudinal, lateral) = point
            .elliptic_patch_axes
            .map_or((point.patch_radius_m, point.patch_radius_m), |axes| {
                (axes.semi_major_axis_m, axes.semi_minor_axis_m)
            });
        let normal_view = NormalPatchView::new(
            request.patch.patch.patch_identity.as_str(),
            request.normal.material.material_card_id.clone(),
            request.normal.material.source_id.clone(),
            NormalPatchAuthority::SyntheticFixture,
            point.normal_force_n,
            longitudinal,
            lateral,
            point.pressure.second_moment_m2,
        )
        .expect("profile-derived normal view");
        let checkpoint = model
            .initial_checkpoint(
                penetrated_state,
                NormalPatchEmbedState::new(0.0, 1.0).expect("normal checkpoint"),
                adapter
                    .initial_state(&normal_view, &interface, 4)
                    .expect("profile tangential checkpoint"),
                RollingContactState::zero(),
                GasChannelState::ExteriorFreeGas(
                    EulerExternalAirWorkState::new("synthetic/exterior-work", 4)
                        .expect("air checkpoint"),
                ),
            )
            .expect("profile production checkpoint");
        let first_request = request.clone();
        let mut repeated_request = production_request_template().0;
        repeated_request.normal.geometry = EulerNormalGeometry::EllipticParaboloid;
        repeated_request.surface_excitation = Some(surface_excitation);
        let repeated_resolved = model
            .bind_horizontal_plane_axisymmetric_profile_contact(
                &mut repeated_request,
                &profile,
                &checkpoint,
                cx,
            )
            .expect("deterministic repeated profile binding");
        assert_eq!(resolved, repeated_resolved);
        assert_eq!(first_request, repeated_request);

        let production_prefix =
            model.run_smooth_contact_trajectory(checkpoint.clone(), 1, |accepted| {
                Ok(request_for_checkpoint(&first_request, accepted))
            });
        assert!(
            !production_prefix.accepted_steps.is_empty(),
            "production prefix refused before its first commit: {:?}",
            production_prefix.termination
        );
        let render_prefix = RenderTrajectory::from_production_coupling_prefix(
            &model,
            &production_prefix,
            &profile,
            first_request.duration_s,
            cx,
        )
        .expect("accepted production prefix enters the common render contract");
        assert_eq!(render_prefix.samples().len(), 1);
        assert_eq!(
            render_prefix
                .metadata()
                .channel_availability
                .normal_force_sampling,
            RenderNormalForceSampling::AppliedSubstepZeroOrderHold
        );
        assert!(!render_prefix.metadata().channel_availability.contact);
        let rendered = render_prefix.samples()[0].input();
        assert_eq!(
            rendered.disposition,
            RenderSampleDisposition::HorizonCensored
        );
        assert_eq!(
            rendered
                .contact_geometry
                .expect("closed patch")
                .support_feature,
            fs_euler_disc_e2e::RenderSupportFeature::ProfileFeature(
                first_request.patch.patch.source_feature,
            )
        );
        assert_eq!(
            rendered.mechanical_energy_j.to_bits(),
            production_prefix.accepted_steps[0]
                .rigid_step
                .diagnostics_after
                .mechanical_energy
                .to_bits()
        );
        let accepted = &production_prefix.accepted_steps[0];
        let nominal_force_n = match &accepted.normal.generic.receipt {
            NormalPatchReceipt::Point(point) => point.normal_force_n,
            NormalPatchReceipt::Line(_) => panic!("profile fixture requires point normal units"),
        };
        assert!(
            rendered.interval_normal_force_n > nominal_force_n,
            "render/audio control must retain the topography force that changed mechanics"
        );
        assert_eq!(
            rendered.interval_normal_force_n.to_bits(),
            accepted
                .base
                .receipt()
                .compressive_normal_force_on_base_n
                .to_bits(),
            "trajectory forcing must be the exact accepted action/reaction load"
        );
        let endpoint_contact = profile_contact_geometry(
            &profile.chart,
            profile.mass_properties,
            accepted.next_disc_state.pose(),
            cx,
        )
        .expect("accepted endpoint support");
        let rendered_contact = rendered.contact_geometry.expect("closed endpoint contact");
        assert_eq!(
            rendered_contact.point_world_m, endpoint_contact.contact.point_world_m,
            "render geometry must be evaluated at the accepted endpoint pose"
        );
        assert_eq!(
            rendered.signed_gap_m.to_bits(),
            (endpoint_contact.contact.gap_m - accepted.base.receipt().modal_displacement_end_m)
                .to_bits()
        );

        let event_source =
            model.run_eventful_compliant_trajectory(checkpoint.clone(), 1, |accepted| {
                let mut input = request_for_checkpoint(&first_request, accepted);
                input.normal.integration_regime =
                    NormalContactIntegrationRegime::CompliantTransient;
                Ok(input)
            });
        assert!(
            matches!(
                event_source.termination,
                ProductionEventTrajectoryTermination::StepLimitReached {
                    maximum_accepted_steps: 1
                }
            ),
            "event-aware production path refused: {:?}",
            event_source.termination
        );
        let event_render = RenderTrajectory::from_production_event_trajectory(
            &model,
            &event_source,
            &profile,
            first_request.duration_s,
            cx,
        )
        .expect("event-aware contact path enters the common render/audio contract");
        let event_sample = event_render.samples()[0].input();
        let ProductionTrajectoryStepReceipt::CompliantContact(event_receipt) =
            &event_source.accepted_steps[0].receipt
        else {
            panic!("penetrated profile must select compliant contact")
        };
        assert!(event_sample.interval_contact_active);
        assert_eq!(
            event_sample.interval_normal_force_n.to_bits(),
            event_receipt
                .base
                .receipt()
                .compressive_normal_force_on_base_n
                .to_bits(),
            "event render/audio forcing must be the exact accepted support reaction"
        );
        assert_eq!(
            event_sample.contact_geometry, rendered.contact_geometry,
            "the smooth and event-aware bridges must query the same endpoint profile geometry"
        );
        let mut forged_prefix = production_prefix.clone();
        assert_ne!(
            forged_prefix.last_accepted_checkpoint.disc_state,
            checkpoint.disc_state
        );
        forged_prefix.last_accepted_checkpoint.disc_state = checkpoint.disc_state;
        assert!(matches!(
            RenderTrajectory::from_production_coupling_prefix(
                &model,
                &forged_prefix,
                &profile,
                first_request.duration_s,
                cx,
            ),
            Err(fs_euler_disc_e2e::RenderTrajectoryError::ProductionPrefixModelMismatch)
        ));

        let first = model
            .step(&checkpoint, &first_request)
            .expect("profile-native production step");
        let replay = model
            .step(&checkpoint, &repeated_request)
            .expect("deterministic profile-native replay");
        assert_eq!(first, replay);
        let (_, receipt) = first;
        assert_eq!(
            receipt.patch_kinematics.patch.curvature,
            first_request.patch.patch.curvature
        );
        assert_eq!(
            receipt.normal.curvature.authority,
            InputAuthority::Estimated
        );
        assert_eq!(
            receipt.rolling.step.generic.patch_authority,
            InputAuthority::Estimated
        );
        assert!(matches!(
            receipt.normal.generic.receipt,
            NormalPatchReceipt::Point(_)
        ));

        let mut stale_mass_profile = profile.clone();
        stale_mass_profile.mass_properties.mass *= 2.0;
        stale_mass_profile
            .mass_properties
            .principal_inertia
            .transverse *= 2.0;
        stale_mass_profile.mass_properties.principal_inertia.axial *= 2.0;
        let stale_mass_model = ProductionCouplingModel {
            identity: model.identity.clone(),
            disc_mass_properties: profile_mbd_mass(stale_mass_profile.mass_properties),
            gravity: Gravity::ZERO,
            base_port: base_port(),
            tangential_adapter: adapter.clone(),
        };
        let mut stale_mass_request = production_request_template().0;
        let stale_mass_request_before = stale_mass_request.clone();
        assert_eq!(
            stale_mass_model.bind_initial_horizontal_plane_axisymmetric_profile_contact(
                &mut stale_mass_request,
                &stale_mass_profile,
                penetrated_state,
                cx,
            ),
            Err(ProductionCouplingError::ResolvedProfileMassMismatch),
            "even a model copied from forged cached mass must not make it authoritative"
        );
        assert_eq!(stale_mass_request, stale_mass_request_before);

        let mut stale_density_profile = profile.clone();
        stale_density_profile.density_kg_per_m3 *= 0.5;
        let mut stale_density_request = production_request_template().0;
        let stale_density_request_before = stale_density_request.clone();
        assert_eq!(
            model.bind_initial_horizontal_plane_axisymmetric_profile_contact(
                &mut stale_density_request,
                &stale_density_profile,
                penetrated_state,
                cx,
            ),
            Err(ProductionCouplingError::ResolvedProfileMassMismatch),
            "retained density and cached mass must describe the same resolved solid"
        );
        assert_eq!(stale_density_request, stale_density_request_before);

        let sharp = DiscProfileSpec::SolidCylinder {
            outer_radius_m: 0.038,
            thickness_m: 0.006,
            edge_treatment: SquatDiscEdgeTreatment::Sharp,
        }
        .resolve(7_800.0, cx)
        .expect("resolved physical sharp profile");
        let sharp_model = ProductionCouplingModel {
            identity: model.identity.clone(),
            disc_mass_properties: profile_mbd_mass(sharp.mass_properties),
            gravity: Gravity::ZERO,
            base_port: base_port(),
            tangential_adapter: adapter,
        };
        let mut sharp_request = production_request_template().0;
        let sharp_request_before = sharp_request.clone();
        assert!(matches!(
            sharp_model.bind_initial_horizontal_plane_axisymmetric_profile_contact(
                &mut sharp_request,
                &sharp,
                penetrated_state,
                cx,
            ),
            Err(ProductionCouplingError::ProfileContact(
                fs_euler_disc_e2e::ContactDynamicsError::ProfileCurvatureRefusal {
                    detail: AxisymmetricCurvatureError::NonsmoothFeatureBoundary { .. }
                }
            ))
        ));
        assert_eq!(
            sharp_request, sharp_request_before,
            "nonsmooth profile refusal must not partially rewrite the step input"
        );
    });
}
