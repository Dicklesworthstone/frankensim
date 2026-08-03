use fs_euler_disc_e2e::{
    CONTACT_NO_CLAIM_BOUNDARY, ContactDiscGeometry as DiscGeometry, ContactDynamicsError,
    ContactDynamicsInput, ContactTermination, ProfileContactDynamicsInput, contact_geometry,
    profile_contact_geometry, refine_profile_timestep_by_two, refine_timestep_by_two,
    run_contact_dynamics, run_profile_contact_dynamics, small_angle_rolling_profile_initializer,
    state_at_ground_contact, state_at_profile_ground_contact,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_mbd::{Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_rep_frep::{AxisymmetricChart, AxisymmetricMassError, SquatDiscEdgeTreatment};
use fs_tribo::{InputAuthority, InterfaceMedium, InterfaceSystemRef};

fn assert_close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "expected {expected:?}, got {actual:?}, tolerance {tolerance:?}"
    );
}

fn with_cx<R>(cancelled: bool, operation: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    if cancelled {
        gate.request();
    }
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x434f_4e54_4143_54,
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

fn tilted_input(
    static_friction_coefficient: f64,
    linear_velocity_world_m_per_s: Vec3,
    spin_rad_s: f64,
    timestep_s: f64,
    maximum_steps: u32,
) -> ContactDynamicsInput {
    let geometry = DiscGeometry {
        radius_m: 0.038,
        thickness_m: 0.006,
        mass_kg: 0.12,
    };
    let orientation = UnitQuaternion::from_axis_angle(Vec3::new(0.0, 1.0, 0.0), 0.55)
        .expect("finite tilted orientation");
    let inertia = geometry
        .mass_properties()
        .expect("physical cylinder inertia");
    let initial_state = state_at_ground_contact(
        geometry,
        orientation,
        linear_velocity_world_m_per_s.scale(geometry.mass_kg),
        Vec3::new(0.0, 0.0, inertia.principal_inertia_body().z * spin_rad_s),
    )
    .expect("tilted rim geometry admits a ground state");
    ContactDynamicsInput {
        geometry,
        initial_state,
        gravity_m_per_s2: 9.806_65,
        static_friction_coefficient,
        interface: InterfaceSystemRef::new(
            "euler-disc-e2e/disc->plane",
            "euler-disc-e2e/contact-dynamics-v1",
            "synthetic-fixture/contact-dynamics-e2e",
            InputAuthority::SyntheticFixture,
            InterfaceMedium::Dry,
        )
        .expect("declared dry interface"),
        timestep_s,
        maximum_steps,
        contact_tolerance_m: 1.0e-9,
        maximum_initial_penetration_m: 1.0e-10,
        release_speed_tolerance_m_per_s: 1.0e-8,
    }
}

#[test]
fn homogeneous_inertia_and_tilted_rim_gap_are_geometry_derived() {
    let input = tilted_input(8.0, Vec3::ZERO, 0.0, 1.0e-4, 1);
    let inertia = input.geometry.mass_properties().expect("physical cylinder");
    let radius_squared = input.geometry.radius_m * input.geometry.radius_m;
    let thickness_squared = input.geometry.thickness_m * input.geometry.thickness_m;
    assert_close(
        inertia.principal_inertia_body().x,
        input.geometry.mass_kg * (3.0 * radius_squared + thickness_squared) / 12.0,
        1.0e-18,
    );
    assert_close(
        inertia.principal_inertia_body().z,
        input.geometry.mass_kg * radius_squared / 2.0,
        1.0e-18,
    );
    let contact = contact_geometry(input.geometry, input.initial_state.pose())
        .expect("tilted cylinder has a unique lowest rim point");
    assert_close(contact.gap_m, 0.0, 1.0e-14);
    assert!(contact.radius_world_m.z < 0.0);
}

#[test]
fn gravity_reaction_and_static_cone_evolve_full_rigid_body_state() {
    let input = tilted_input(100.0, Vec3::ZERO, 0.05, 1.0e-4, 8);
    let run = run_contact_dynamics(&input).expect("admitted high-friction dynamic contact run");
    assert_eq!(run.termination, ContactTermination::HorizonReached);
    assert_eq!(run.steps.len(), 8);
    let first = run.steps.first().expect("completed first dynamic step");
    assert!(first.normal_impulse_ns > 0.0);
    assert!(first.stick.feasible);
    assert!(first.stick.normal_impulse_ns > 0.0);
    assert!(first.stick.static_capacity_impulse_ns >= first.stick.required_tangential_impulse_ns);
    assert_eq!(
        first.stick.input_authority,
        InputAuthority::SyntheticFixture
    );
    assert!(first.contact_after.gap_m.abs() <= input.contact_tolerance_m);
    assert!(first.energy.mechanical_energy_before_j.is_finite());
    assert!(first.energy.mechanical_energy_after_j.is_finite());
    assert!(first.energy.mechanical_balance_residual_j.is_finite());
    assert!(
        first
            .post_impulse_contact_velocity_world_m_per_s
            .is_finite()
    );
    assert_close(
        first.post_impulse_contact_velocity_residual_world_m_per_s.x,
        0.0,
        1.0e-10,
    );
    assert_close(
        first.post_impulse_contact_velocity_residual_world_m_per_s.y,
        0.0,
        1.0e-10,
    );
    assert_close(
        first.post_impulse_contact_velocity_residual_world_m_per_s.z,
        0.0,
        1.0e-10,
    );
    assert_close(
        first.post_impulse_contact_velocity_world_m_per_s.x,
        0.0,
        1.0e-10,
    );
    assert_close(
        first.post_impulse_contact_velocity_world_m_per_s.y,
        0.0,
        1.0e-10,
    );
    assert_close(
        first.post_impulse_contact_velocity_world_m_per_s.z,
        0.0,
        1.0e-10,
    );
    assert_ne!(
        first.state_before.pose().orientation().components(),
        first.state_after.pose().orientation().components(),
        "orientation must be integrated from state, not prescribed by theta"
    );
}

#[test]
fn hostile_initial_penetration_is_refused_instead_of_projected_away() {
    let mut input = tilted_input(4.0, Vec3::ZERO, 0.0, 1.0e-4, 1);
    let initial = input.initial_state;
    let penetrated_pose = Pose::new(
        initial
            .pose()
            .position_world()
            .sub(Vec3::new(0.0, 0.0, 1.0e-6)),
        initial.pose().orientation(),
    )
    .expect("finite hostile pose");
    input.initial_state = RigidBodyState::new(
        penetrated_pose,
        initial.linear_momentum_world(),
        initial.angular_momentum_body(),
    )
    .expect("finite hostile state");
    match run_contact_dynamics(&input) {
        Err(ContactDynamicsError::InitialPenetrationExceeded { gap_m, tolerance_m }) => {
            assert!(gap_m < -tolerance_m);
            assert_eq!(tolerance_m, input.maximum_initial_penetration_m);
        }
        other => panic!("expected penetration refusal, got {other:?}"),
    }
}

#[test]
fn geometry_relative_gap_tolerances_refuse_unbounded_projection_policy() {
    let mut input = tilted_input(4.0, Vec3::ZERO, 0.0, 1.0e-4, 1);
    input.contact_tolerance_m = 1.0e-3;
    match run_contact_dynamics(&input) {
        Err(ContactDynamicsError::InvalidInput {
            field: "geometry_relative_contact_tolerance",
        }) => {}
        other => panic!("expected geometry-relative tolerance refusal, got {other:?}"),
    }
}

#[test]
fn finite_gap_excess_reports_constraint_residual_instead_of_a_projection_claim() {
    let mut input = tilted_input(100.0, Vec3::ZERO, 0.0, 1.0e-4, 1);
    input.contact_tolerance_m = 1.0e-15;
    match run_contact_dynamics(&input) {
        Err(ContactDynamicsError::ConstraintResidualExceeded {
            field: "projected_contact_gap_m",
            residual,
            tolerance,
        }) => {
            assert!(residual.is_finite());
            assert!(residual > tolerance);
        }
        other => panic!("expected finite projected-gap residual, got {other:?}"),
    }
}

#[test]
fn horizontal_cylinder_line_contact_is_refused() {
    let input = tilted_input(4.0, Vec3::ZERO, 0.0, 1.0e-4, 1);
    let horizontal =
        UnitQuaternion::from_axis_angle(Vec3::new(0.0, 1.0, 0.0), core::f64::consts::FRAC_PI_2)
            .expect("finite horizontal orientation");
    match contact_geometry(
        input.geometry,
        Pose::new(Vec3::ZERO, horizontal).expect("finite pose"),
    ) {
        Err(ContactDynamicsError::UnsupportedLineContact) => {}
        other => panic!("expected line-contact refusal, got {other:?}"),
    }
}

#[test]
fn impulse_work_uses_free_contact_velocity_and_projection_is_separate() {
    let input = tilted_input(4.0, Vec3::ZERO, 0.0, 1.0e-4, 1);
    let run = run_contact_dynamics(&input).expect("one resting-contact step");
    let receipt = run.steps.first().expect("completed contact step");
    assert!(
        receipt.energy.contact_impulse_work_estimate_j < 0.0,
        "gravity makes the free contact velocity negative before the restoring impulse"
    );
    let recomposed = receipt.energy.gravity_work_j
        + receipt.energy.contact_impulse_work_estimate_j
        + receipt.energy.geometric_projection_work_j
        + receipt.energy.mechanical_balance_residual_j;
    assert_close(
        receipt.energy.mechanical_energy_delta_j,
        recomposed,
        1.0e-15,
    );
    assert_close(
        receipt.energy.geometric_projection_work_j,
        receipt.energy.projection_potential_shift_j,
        1.0e-15,
    );
}

#[test]
fn separating_contact_terminates_without_a_fake_reaction() {
    let input = tilted_input(2.0, Vec3::new(0.0, 0.0, 0.2), 0.0, 1.0e-4, 4);
    let run = run_contact_dynamics(&input).expect("separation is a terminal result, not refusal");
    assert!(run.steps.is_empty());
    match run.termination {
        ContactTermination::ContactLost {
            step_index,
            gap_m,
            normal_velocity_m_per_s,
        } => {
            assert_eq!(step_index, 0);
            assert_close(gap_m, 0.0, 1.0e-14);
            assert!(normal_velocity_m_per_s > input.release_speed_tolerance_m_per_s);
        }
        other => panic!("expected contact loss, got {other:?}"),
    }
}

#[test]
fn initially_separated_contact_reports_its_present_normal_velocity() {
    let mut input = tilted_input(2.0, Vec3::new(0.0, 0.0, 0.2), 0.0, 1.0e-4, 4);
    let state = input.initial_state;
    let separated_pose = Pose::new(
        state
            .pose()
            .position_world()
            .add(Vec3::new(0.0, 0.0, 1.0e-3)),
        state.pose().orientation(),
    )
    .expect("finite separated pose");
    input.initial_state = RigidBodyState::new(
        separated_pose,
        state.linear_momentum_world(),
        state.angular_momentum_body(),
    )
    .expect("finite separated state");
    let run = run_contact_dynamics(&input).expect("separation is a terminal result");
    match run.termination {
        ContactTermination::ContactLost {
            step_index,
            gap_m,
            normal_velocity_m_per_s,
        } => {
            assert_eq!(step_index, 0);
            assert!(gap_m > input.contact_tolerance_m);
            assert_close(normal_velocity_m_per_s, 0.2, 1.0e-12);
        }
        other => panic!("expected already-separated contact loss, got {other:?}"),
    }
}

#[test]
fn infeasible_static_friction_terminates_instead_of_prescribing_no_slip() {
    let input = tilted_input(0.001, Vec3::new(0.8, 0.0, 0.0), 0.0, 1.0e-4, 4);
    let run = run_contact_dynamics(&input).expect("cone failure is an honest terminal result");
    assert!(run.steps.is_empty());
    match run.termination {
        ContactTermination::StickInfeasible { step_index, stick } => {
            assert_eq!(step_index, 0);
            assert!(!stick.feasible);
            assert!(stick.required_tangential_impulse_ns > stick.static_capacity_impulse_ns);
            assert!(stick.friction_cone_margin_ns < 0.0);
        }
        other => panic!("expected stick infeasibility, got {other:?}"),
    }
}

#[test]
fn fixed_horizon_refinement_is_deterministic_and_reports_endpoint_difference() {
    let input = tilted_input(4.0, Vec3::ZERO, 0.0, 2.0e-4, 6);
    let first = run_contact_dynamics(&input).expect("coarse run");
    let repeated = run_contact_dynamics(&input).expect("repeated coarse run");
    assert_eq!(first, repeated);
    let refinement = refine_timestep_by_two(&input).expect("same-terminal-class refinement");
    assert_eq!(refinement.coarse, first);
    assert_eq!(
        refinement.coarse.termination,
        ContactTermination::HorizonReached
    );
    assert_eq!(
        refinement.fine.termination,
        ContactTermination::HorizonReached
    );
    assert_eq!(
        refinement.reference.termination,
        ContactTermination::HorizonReached
    );
    assert!(refinement.final_position_difference_m.is_finite());
    assert!(
        refinement
            .final_linear_momentum_difference_kg_m_per_s
            .is_finite()
    );
    assert!(
        refinement
            .final_angular_momentum_difference_kg_m2_per_s
            .is_finite()
    );
    assert!(refinement.final_mechanical_energy_difference_j.is_finite());
    assert!(
        refinement
            .coarse_reference_position_difference_m
            .is_finite()
    );
    assert!(refinement.fine_reference_position_difference_m.is_finite());
    assert!(
        refinement
            .coarse_reference_linear_momentum_difference_kg_m_per_s
            .is_finite()
    );
    assert!(
        refinement
            .fine_reference_linear_momentum_difference_kg_m_per_s
            .is_finite()
    );
    assert!(
        refinement
            .coarse_reference_angular_momentum_difference_kg_m2_per_s
            .is_finite()
    );
    assert!(
        refinement
            .fine_reference_angular_momentum_difference_kg_m2_per_s
            .is_finite()
    );
    assert!(
        refinement
            .coarse_reference_orientation_angle_rad
            .is_finite()
    );
    assert!(refinement.fine_reference_orientation_angle_rad.is_finite());
    assert!(CONTACT_NO_CLAIM_BOUNDARY.contains("No sliding"));
    assert!(CONTACT_NO_CLAIM_BOUNDARY.contains("convergence-order claim"));
}

#[test]
fn sharp_and_one_millimetre_fillet_have_distinct_oblique_profile_contacts_and_runs() {
    with_cx(false, |cx| {
        let sharp = AxisymmetricChart::squat_disc(0.038, 0.006, SquatDiscEdgeTreatment::Sharp)
            .expect("physical sharp disc");
        let filleted = AxisymmetricChart::squat_disc(
            0.038,
            0.006,
            SquatDiscEdgeTreatment::CircularFillet { radius: 0.001 },
        )
        .expect("physical 1 mm fillet");
        let density_kg_per_m3 = 7_800.0;
        let sharp_mass = sharp
            .mass_properties(density_kg_per_m3, cx)
            .expect("sharp profile mass");
        let filleted_mass = filleted
            .mass_properties(density_kg_per_m3, cx)
            .expect("filleted profile mass");
        let orientation = UnitQuaternion::from_axis_angle(Vec3::new(0.0, 1.0, 0.0), 0.55)
            .expect("finite tilted orientation");
        let zero_pose = Pose::new(Vec3::ZERO, orientation).expect("finite support pose");
        let sharp_contact = profile_contact_geometry(&sharp, sharp_mass, zero_pose, cx)
            .expect("oblique ground direction has a unique sharp rim vertex");
        let filleted_contact = profile_contact_geometry(&filleted, filleted_mass, zero_pose, cx)
            .expect("unique filleted support");
        let support_delta = sharp_contact
            .contact
            .radius_world_m
            .sub(filleted_contact.contact.radius_world_m);
        assert!(
            support_delta.dot(support_delta).sqrt() > 1.0e-6,
            "sharp and 1 mm fillet must supply distinct oblique contact inputs"
        );
        assert!(
            (sharp_mass.mass - filleted_mass.mass).abs() > 1.0e-9,
            "the profile mass must not silently retain a sharp-cylinder value"
        );
        assert!(
            (sharp_mass.principal_inertia.axial - filleted_mass.principal_inertia.axial).abs()
                > 1.0e-12,
            "the profile inertia must enter the rigid-body model"
        );
        assert_eq!(
            filleted_contact.mass_properties, filleted_mass,
            "the admitted contact records the mass properties it actually uses"
        );

        let inclination_rad = 0.05;
        let sharp_initial = small_angle_rolling_profile_initializer(
            &sharp,
            density_kg_per_m3,
            inclination_rad,
            9.806_65,
            cx,
        )
        .expect("sharp oblique rolling initialization");
        let filleted_initial = small_angle_rolling_profile_initializer(
            &filleted,
            density_kg_per_m3,
            inclination_rad,
            9.806_65,
            cx,
        )
        .expect("filleted oblique rolling initialization");
        let mut sharp_controls = tilted_input(100.0, Vec3::ZERO, 0.0, 1.0e-6, 2);
        sharp_controls.geometry.mass_kg = sharp_mass.mass;
        sharp_controls.initial_state = sharp_initial.state;
        let mut filleted_controls = tilted_input(100.0, Vec3::ZERO, 0.0, 1.0e-6, 2);
        filleted_controls.geometry.mass_kg = filleted_mass.mass;
        filleted_controls.initial_state = filleted_initial.state;
        let sharp_run = run_profile_contact_dynamics(
            &ProfileContactDynamicsInput {
                chart: sharp,
                density_kg_per_m3,
                controls: sharp_controls,
            },
            cx,
        )
        .expect("sharp oblique profile evolution");
        let filleted_run = run_profile_contact_dynamics(
            &ProfileContactDynamicsInput {
                chart: filleted,
                density_kg_per_m3,
                controls: filleted_controls,
            },
            cx,
        )
        .expect("filleted oblique profile evolution");
        assert!(!sharp_run.steps.is_empty());
        assert!(!filleted_run.steps.is_empty());
    });
}

#[test]
fn declared_small_angle_profile_roll_is_zero_slip_and_has_refinement_diagnostics() {
    with_cx(false, |cx| {
        let chart = AxisymmetricChart::squat_disc(
            0.038,
            0.006,
            SquatDiscEdgeTreatment::CircularFillet { radius: 0.001 },
        )
        .expect("physical 1 mm fillet");
        let density_kg_per_m3 = 7_800.0;
        let mass = chart
            .mass_properties(density_kg_per_m3, cx)
            .expect("profile mass");
        let inclination_rad = 0.05;
        let rolling = small_angle_rolling_profile_initializer(
            &chart,
            density_kg_per_m3,
            inclination_rad,
            9.806_65,
            cx,
        )
        .expect("caller-declared rolling-compatible profile state");
        let body_axis_world = rolling
            .state
            .pose()
            .orientation()
            .rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0));
        assert_close(body_axis_world.x, inclination_rad.sin(), 1.0e-14);
        assert_close(body_axis_world.y, 0.0, 1.0e-14);
        assert_close(body_axis_world.z, inclination_rad.cos(), 1.0e-14);
        assert_close(
            rolling.declared_precession_rate_rad_per_s.powi(2) * inclination_rad.sin(),
            4.0 * 9.806_65 / 0.038,
            1.0e-10,
        );
        assert_close(
            rolling.angular_velocity_body_rad_per_s.x,
            -rolling.declared_precession_rate_rad_per_s * inclination_rad.sin(),
            1.0e-13,
        );
        assert_close(rolling.angular_velocity_body_rad_per_s.y, 0.0, 1.0e-15);
        assert_close(rolling.angular_velocity_body_rad_per_s.z, 0.0, 1.0e-15);
        assert_close(rolling.contact.contact.gap_m, 0.0, 1.0e-14);
        assert_close(
            rolling.initial_contact_velocity_world_m_per_s.x,
            0.0,
            1.0e-12,
        );
        assert_close(
            rolling.initial_contact_velocity_world_m_per_s.y,
            0.0,
            1.0e-12,
        );
        assert_close(
            rolling.initial_contact_velocity_world_m_per_s.z,
            0.0,
            1.0e-12,
        );
        let mut controls = tilted_input(100.0, Vec3::ZERO, 0.0, 1.0e-6, 4);
        controls.geometry.mass_kg = mass.mass;
        controls.initial_state = rolling.state;
        let input = ProfileContactDynamicsInput {
            chart,
            density_kg_per_m3,
            controls,
        };
        let first = run_profile_contact_dynamics(&input, cx).expect("profile contact evolution");
        let repeated = run_profile_contact_dynamics(&input, cx).expect("repeat profile evolution");
        assert_eq!(first, repeated);
        assert!(
            !first.steps.is_empty(),
            "the full profile contact solver must evolve the declared initial state"
        );
        let refinement = refine_profile_timestep_by_two(&input, cx).expect("profile refinement");
        assert_eq!(refinement.coarse, first);
        assert!(
            refinement
                .coarse_reference_position_difference_m
                .is_finite()
        );
        assert!(refinement.fine_reference_position_difference_m.is_finite());
        assert!(
            refinement
                .coarse_reference_linear_momentum_difference_kg_m_per_s
                .is_finite()
        );
        assert!(
            refinement
                .fine_reference_linear_momentum_difference_kg_m_per_s
                .is_finite()
        );
        assert!(
            refinement
                .coarse_reference_angular_momentum_difference_kg_m2_per_s
                .is_finite()
        );
        assert!(
            refinement
                .fine_reference_angular_momentum_difference_kg_m2_per_s
                .is_finite()
        );
        assert!(
            refinement
                .coarse_reference_orientation_angle_rad
                .is_finite()
        );
        assert!(refinement.fine_reference_orientation_angle_rad.is_finite());
    });
}

#[test]
fn profile_initializer_has_zero_fillet_gap_and_cylinder_initializer_does_not() {
    with_cx(false, |cx| {
        let chart = AxisymmetricChart::squat_disc(
            0.038,
            0.006,
            SquatDiscEdgeTreatment::CircularFillet { radius: 0.001 },
        )
        .expect("physical 1 mm fillet");
        let density_kg_per_m3 = 7_800.0;
        let mass = chart
            .mass_properties(density_kg_per_m3, cx)
            .expect("profile mass");
        let orientation = UnitQuaternion::from_axis_angle(Vec3::new(0.0, 1.0, 0.0), 0.55)
            .expect("finite tilted orientation");
        let profile_state = state_at_profile_ground_contact(
            &chart,
            density_kg_per_m3,
            orientation,
            Vec3::ZERO,
            Vec3::ZERO,
            cx,
        )
        .expect("profile ground initialization");
        let profile_gap = profile_contact_geometry(&chart, mass, profile_state.pose(), cx)
            .expect("unique profile contact")
            .contact
            .gap_m;
        assert_close(profile_gap, 0.0, 1.0e-14);

        let cylinder_geometry = DiscGeometry {
            radius_m: 0.038,
            thickness_m: 0.006,
            mass_kg: mass.mass,
        };
        let cylinder_state =
            state_at_ground_contact(cylinder_geometry, orientation, Vec3::ZERO, Vec3::ZERO)
                .expect("legacy cylinder ground initialization");
        let cylinder_gap = profile_contact_geometry(&chart, mass, cylinder_state.pose(), cx)
            .expect("unique profile contact")
            .contact
            .gap_m;
        assert!(
            cylinder_gap.abs() > 1.0e-6,
            "a cylinder initializer must not be accepted as a fillet initializer"
        );
    });
}

#[test]
fn profile_run_rejects_inconsistent_legacy_geometry_declaration() {
    with_cx(false, |cx| {
        let chart = AxisymmetricChart::squat_disc(
            0.038,
            0.006,
            SquatDiscEdgeTreatment::CircularFillet { radius: 0.001 },
        )
        .expect("physical 1 mm fillet");
        let density_kg_per_m3 = 7_800.0;
        let mass = chart
            .mass_properties(density_kg_per_m3, cx)
            .expect("profile mass");
        let orientation = UnitQuaternion::from_axis_angle(Vec3::new(0.0, 1.0, 0.0), 0.55)
            .expect("finite tilted orientation");
        let state = state_at_profile_ground_contact(
            &chart,
            density_kg_per_m3,
            orientation,
            Vec3::ZERO,
            Vec3::ZERO,
            cx,
        )
        .expect("profile ground initialization");
        let mut controls = tilted_input(100.0, Vec3::ZERO, 0.0, 1.0e-4, 1);
        controls.initial_state = state;
        assert_ne!(controls.geometry.mass_kg, mass.mass);
        let input = ProfileContactDynamicsInput {
            chart,
            density_kg_per_m3,
            controls,
        };
        match run_profile_contact_dynamics(&input, cx) {
            Err(ContactDynamicsError::ProfileControlMismatch {
                field: "controls.geometry.mass_kg",
                declared,
                derived,
            }) => assert_ne!(declared, derived),
            other => panic!("expected profile-control declaration refusal, got {other:?}"),
        }
        input.controls.geometry.mass_kg = mass.mass;
        input.controls.geometry.radius_m = 0.037;
        match run_profile_contact_dynamics(&input, cx) {
            Err(ContactDynamicsError::ProfileControlMismatch {
                field: "controls.geometry.radius_m",
                declared,
                derived,
            }) => assert_ne!(declared, derived),
            other => panic!("expected profile-dimension declaration refusal, got {other:?}"),
        }
    });
}

#[test]
fn cancelled_profile_mass_query_is_a_typed_refusal() {
    let chart = AxisymmetricChart::squat_disc(
        0.038,
        0.006,
        SquatDiscEdgeTreatment::CircularFillet { radius: 0.001 },
    )
    .expect("physical profile");
    let input = ProfileContactDynamicsInput {
        chart,
        density_kg_per_m3: 7_800.0,
        controls: tilted_input(4.0, Vec3::ZERO, 0.0, 1.0e-4, 1),
    };
    with_cx(true, |cx| match run_profile_contact_dynamics(&input, cx) {
        Err(ContactDynamicsError::ProfileMassRefusal {
            detail: AxisymmetricMassError::Cancelled,
        }) => {}
        other => panic!("expected typed cancelled profile-mass refusal, got {other:?}"),
    });
}
