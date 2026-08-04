//! Synthetic G0 composition fixture; it establishes plumbing, not calibration.

use fs_contact::normal_patch::{
    ApplicabilityInput, ApplicabilityLimits, InputUncertainty, NormalPatchEmbedState,
    NormalPatchReceipt,
};
use fs_couple::StableId;
use fs_euler_disc_e2e::base_response::{ReducedBasePort, ReducedBasePortIdentity};
use fs_euler_disc_e2e::external_air::{
    EulerDiscBodyFrame, EulerDiscExteriorGeometry, EulerDiscExteriorState, EulerExternalAirInput,
    EulerExternalAirWorkState, ExteriorAirPressure, ExternalAirDomain, ExternalAirIdentity,
};
use fs_euler_disc_e2e::normal_contact::{
    EulerNormalContactInput, EulerNormalContactOutcome, EulerNormalGeometry,
    NORMAL_CONTACT_ADAPTER_ID, NormalContactIdentity, NormalMaterialInterface,
    evaluate_normal_contact,
};
use fs_euler_disc_e2e::patch_kinematics::{
    CurvatureMetadata, MovingOneModeBaseState, MovingOneModePatchBridgeInput,
    MovingOneModePatchKinematicsInput, OrderedSurfacePair, PatchGeometryMetadata,
    PatchKinematicThresholds, ProfileSupportKinematics, SurfaceOrder, TangentGaugeInput,
    compute_moving_one_mode_patch_kinematics,
};
use fs_euler_disc_e2e::production_coupling::{
    ProductionCouplingError, ProductionCouplingIdentity, ProductionCouplingModel,
    ProductionCouplingStepInput, SmoothContactTrajectoryTermination,
};
use fs_euler_disc_e2e::rolling_contact::{
    ROLLING_CONTACT_ADAPTER_ID, RollingContactIdentity, RollingContactInput, RollingContactState,
};
use fs_euler_disc_e2e::tangential_contact::{
    EulerTangentialContactAdapter, TangentialContactLane, TangentialContactRequest,
};
use fs_euler_disc_e2e::{
    BaseGeometryScope, BaseResponseInput, ContactLoadScope, LevelSupportInput, MovingContactLoad,
};
use fs_flux::Vec3 as FluxVec3;
use fs_flux::{
    ApplicabilityEnvelope, ClosedRange, ContributionFamily, CorrelationIdentity,
    CorrelationUncertainty, FormDrag, GasPropertyCard, ReducedAeroComponents, ReducedAeroModel,
    SurfaceRoughness,
};
use fs_mbd::{Gravity, MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
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
use fs_tribo::{ApplicabilityRange, InputAuthority, InterfaceMedium, InterfaceSystemRef};

fn id(value: &str) -> StableId {
    StableId::new(value).expect("synthetic fixture identity")
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
        hunt_crossley_dissipation_s_per_m: None,
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
    let disc_position = checkpoint.disc_state.pose().position_world();
    input.expected_checkpoint_version = version;
    input.time_s = version as f64 * input.duration_s;
    input.normal.iteration = version + 1;
    input.patch.bridge.profile_support.disc_point_world_m = disc_position;
    input
        .patch
        .bridge
        .base_mode
        .undeformed_contact_point_world_m = Vec3::new(disc_position.x, disc_position.y, 0.0);
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
    input.exterior_exchange_key = version + 1;
    input.base_step_id = format!("synthetic/base-step-{}", version + 1);
    input
}

#[test]
fn synthetic_g0_one_substep_composes_real_adapters_and_refuses_without_mutation() {
    let kinematics = compute_moving_one_mode_patch_kinematics(patch_input())
        .expect("synthetic moving-base patch");
    let normal = normal_input(kinematics.clone());
    let normal_outcome = evaluate_normal_contact(&normal).expect("synthetic normal");
    assert!(matches!(
        normal_outcome,
        EulerNormalContactOutcome::Active(_)
    ));
    let EulerNormalContactOutcome::Active(active) = normal_outcome else {
        return;
    };
    assert!(matches!(
        active.generic.receipt,
        NormalPatchReceipt::Point(_)
    ));
    let NormalPatchReceipt::Point(point) = active.generic.receipt else {
        return;
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
            EulerExternalAirWorkState::new("synthetic/exterior-work", 4).expect("air checkpoint"),
        )
        .expect("production checkpoint");
    let request = ProductionCouplingStepInput {
        expected_checkpoint_version: 0,
        duration_s: 1.0e-6,
        time_s: 0.0,
        patch: patch_input(),
        normal,
        tangential: TangentialContactRequest {
            request_id: "synthetic/tangent-step".into(),
            expected_state_version: 0,
            patch_kinematics: kinematics,
            normal_patch: normal_view,
            interface,
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
        exterior_air: aero(),
        selected_exterior_correlation_id: "synthetic/exterior-a".into(),
        exterior_exchange_key: 1,
        base_step_id: "synthetic/base-step-1".into(),
        base_load_progress_start: 0.0,
        base_load_progress_end: 0.0,
    };
    let (next, receipt) = model
        .step(&checkpoint, &request)
        .expect("synthetic composition");
    assert_eq!(next.committed_version, 1);
    assert!(receipt.estimate_only);
    assert!(receipt.total_force_world_n.is_finite());
    assert!(
        receipt.exterior_air.world_wrench.force_world_n.x < 0.0,
        "exterior force must use the checkpoint's +x velocity, not the stale card state"
    );
    let mut stale_version = request.clone();
    stale_version.expected_checkpoint_version = 1;
    assert!(matches!(
        model.step(&checkpoint, &stale_version),
        Err(ProductionCouplingError::CheckpointVersionMismatch {
            expected: 1,
            observed: 0
        })
    ));
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
    assert_eq!(whole.accepted_steps.len(), 3);
    assert!(matches!(
        whole.termination,
        SmoothContactTrajectoryTermination::StepLimitReached {
            maximum_accepted_steps: 3
        }
    ));
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
    thin_gap.exterior_air.domain = ExternalAirDomain::ThinGap;
    assert!(matches!(
        model.step(&checkpoint, &thin_gap),
        Err(ProductionCouplingError::ExteriorAir(_))
    ));
    assert_eq!(
        checkpoint.committed_version, 0,
        "refusals cannot mutate checkpoint"
    );
}
