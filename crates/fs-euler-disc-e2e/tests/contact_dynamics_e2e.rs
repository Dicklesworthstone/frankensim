#[path = "../src/contact_dynamics.rs"]
mod contact_dynamics;

use contact_dynamics::{
    ContactDynamicsError, ContactDynamicsInput, ContactTermination, DiscGeometry,
    NO_CLAIM_BOUNDARY, contact_geometry, refine_timestep_by_two, run_contact_dynamics,
    state_at_ground_contact,
};
use fs_mbd::{Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_tribo::{InputAuthority, InterfaceMedium, InterfaceSystemRef};

fn assert_close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "expected {expected:?}, got {actual:?}, tolerance {tolerance:?}"
    );
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
    assert!(first.stick.normal_reaction_n > 0.0);
    assert!(first.stick.static_capacity_n >= first.stick.required_tangential_reaction_n);
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
fn impulse_work_uses_free_contact_velocity_and_projection_is_separate() {
    let input = tilted_input(4.0, Vec3::ZERO, 0.0, 1.0e-4, 1);
    let run = run_contact_dynamics(&input).expect("one resting-contact step");
    let receipt = run.steps.first().expect("completed contact step");
    assert!(
        receipt.energy.contact_impulse_work_estimate_j < 0.0,
        "gravity makes the free contact velocity negative before the restoring impulse"
    );
    let recomposed = receipt.energy.contact_impulse_work_estimate_j
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
fn infeasible_static_friction_terminates_instead_of_prescribing_no_slip() {
    let input = tilted_input(0.001, Vec3::new(0.8, 0.0, 0.0), 0.0, 1.0e-4, 4);
    let run = run_contact_dynamics(&input).expect("cone failure is an honest terminal result");
    assert!(run.steps.is_empty());
    match run.termination {
        ContactTermination::StickInfeasible { step_index, stick } => {
            assert_eq!(step_index, 0);
            assert!(!stick.feasible);
            assert!(stick.required_tangential_reaction_n > stick.static_capacity_n);
            assert!(stick.friction_cone_margin_n < 0.0);
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
    assert!(NO_CLAIM_BOUNDARY.contains("No sliding"));
    assert!(NO_CLAIM_BOUNDARY.contains("convergence-order claim"));
}
