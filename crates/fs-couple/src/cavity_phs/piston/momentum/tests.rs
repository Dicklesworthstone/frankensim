//! Native source gas, rigid inertia and continued wall motion after impact.

use super::*;
use crate::cavity_phs::piston::double_acting::DoubleActingGasPiston;
use crate::cavity_phs::piston::tests::{
    controls as step_controls, drive, gas, real_gas::real_parcel,
};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_ports::PortedModule;

fn near(actual: f64, expected: f64) {
    assert!((actual - expected).abs() <= 1e-12, "{actual} != {expected}");
}

fn contact(arm: Vec3, restitution: f64) -> RigidPistonImpactContact {
    RigidPistonImpactContact {
        rigid_arm_body_m: arm,
        piston_axis_origin_world_m: Vec3::ZERO,
        expansion_axis_world: Vec3::new(0.0, 0.0, 1.0),
        normal_sign: 1.0,
        activation_tolerance_m: 1e-9,
        restitution,
    }
}

fn rigid<G: PistonGas>(wall: &GasPiston<G>, arm: Vec3) -> (MassProperties, RigidBodyState) {
    let properties = MassProperties::new(1.0, Vec3::ZERO, Vec3::new(1.0, 2.0, 3.0)).unwrap();
    let pose = Pose::new(
        Vec3::new(0.0, 0.0, wall.position_m()).sub(arm),
        UnitQuaternion::IDENTITY,
    )
    .unwrap();
    let body = RigidBodyState::new(pose, Vec3::new(0.0, 0.0, -3.0), Vec3::ZERO).unwrap();
    (properties, body)
}

fn same_ports(left: &impl PortedModule, right: &impl PortedModule) {
    assert_eq!(left.port_schemas(), right.port_schemas());
    assert_eq!(
        left.port_coordinate_charts(),
        right.port_coordinate_charts()
    );
    assert_eq!(left.junction_schemas(), right.junction_schemas());
    assert_eq!(
        left.last_window_duration_s(),
        right.last_window_duration_s()
    );
    assert_eq!(left.last_window_transfers(), right.last_window_transfers());
    assert_eq!(
        left.last_window_transfer_error(),
        right.last_window_transfer_error()
    );
    assert_eq!(left.audit_last_window(), right.audit_last_window());
}

#[test]
fn reciprocal_two_mass_restitution_preserves_gas_and_previous_port_window() {
    for restitution in [0.0, 0.5, 1.0] {
        let mut wall = GasPiston::new(gas(0.75, 0.001), 2.0, 0.01, 0.0).unwrap();
        wall.step(drive(0.0), step_controls(), || false).unwrap();
        let before = wall.clone();
        let (properties, mut body) = rigid(&wall, Vec3::ZERO);
        let pose_before = body.pose();
        let receipt = wall
            .apply_rigid_normal_impact(
                &mut body,
                properties,
                contact(Vec3::ZERO, restitution),
                RigidPistonImpactControls::default(),
                || false,
            )
            .unwrap();
        let impulse = 2.0 * (1.0 + restitution);
        near(receipt.contact.normal_impulse, impulse);
        near(receipt.inverse_effective_mass_kg_inv, 1.5);
        near(wall.velocity_m_s(), -0.5 * impulse);
        near(body.linear_momentum_world().z, -3.0 + impulse);
        near(receipt.relative_normal_speed_after_m_s, 3.0 * restitution);
        near(
            receipt.contact.dissipated_energy,
            3.0 * (1.0 - restitution * restitution),
        );
        near(
            receipt.kinetic_energy_change_j,
            -receipt.contact.dissipated_energy,
        );
        near(
            body.linear_momentum_world().z + 2.0 * wall.velocity_m_s(),
            -3.0,
        );
        near(receipt.energy_residual_j, 0.0);
        assert_eq!(body.pose(), pose_before);
        assert_eq!(wall.gas(), before.gas());
        same_ports(&wall, &before);
    }
}

#[test]
fn off_center_native_rigid_impact_retains_rotational_recoil_and_energy() {
    let mut wall = GasPiston::new(gas(0.75, 0.001), 2.0, 0.01, 0.0).unwrap();
    let before = wall.clone();
    let arm = Vec3::new(1.0, 0.0, 0.0);
    let (properties, mut body) = rigid(&wall, arm);
    let receipt = wall
        .apply_rigid_normal_impact(
            &mut body,
            properties,
            contact(arm, 1.0),
            RigidPistonImpactControls::default(),
            || false,
        )
        .unwrap();
    near(receipt.inverse_effective_mass_kg_inv, 2.0);
    near(receipt.contact.normal_impulse, 3.0);
    near(body.linear_momentum_world().z, 0.0);
    near(body.angular_momentum_body().y, -3.0);
    near(
        receipt.rigid.point_kinematics_after.angular_velocity_body.y,
        -1.5,
    );
    near(
        receipt.rigid.point_kinematics_after.point_velocity_world.z,
        1.5,
    );
    near(wall.velocity_m_s(), -1.5);
    near(receipt.relative_normal_speed_after_m_s, 3.0);
    near(receipt.kinetic_energy_change_j, 0.0);
    near(receipt.energy_residual_j, 0.0);
    assert_eq!(wall.gas(), before.gas());
}

#[test]
fn fixed_wall_rebounds_rigid_body_and_returns_zero_work_support_reaction() {
    let mut wall = GasPiston::new(gas(0.75, 0.001), 2.0, 0.01, 0.0)
        .unwrap()
        .with_fixed_support()
        .unwrap();
    let before = wall.clone();
    let (properties, mut body) = rigid(&wall, Vec3::ZERO);
    let receipt = wall
        .apply_rigid_normal_impact(
            &mut body,
            properties,
            contact(Vec3::ZERO, 1.0),
            RigidPistonImpactControls::default(),
            || false,
        )
        .unwrap();
    near(receipt.contact.normal_impulse, 6.0);
    near(body.linear_momentum_world().z, 3.0);
    near(receipt.piston.impulse_n_s, -6.0);
    near(receipt.piston.support_impulse_n_s, 6.0);
    assert_eq!(receipt.piston.support_work_j, 0.0);
    assert_eq!(receipt.piston.external_work_j, 0.0);
    assert_eq!(wall, before);
    same_ports(&wall, &before);
}

#[test]
fn late_cancellation_and_endpoint_refusal_leave_both_owners_and_history_intact() {
    let mut wall = GasPiston::new(gas(0.75, 0.001), 2.0, 0.01, 0.0).unwrap();
    wall.step(drive(0.0), step_controls(), || false).unwrap();
    let before = wall.clone();
    let (properties, mut body) = rigid(&wall, Vec3::ZERO);
    let body_before = body;
    let mut polls = 0;
    let result = wall.apply_rigid_normal_impact(
        &mut body,
        properties,
        contact(Vec3::ZERO, 1.0),
        RigidPistonImpactControls::default(),
        || {
            polls += 1;
            polls == 2
        },
    );
    assert!(matches!(
        result,
        Err(RigidPistonImpactError::Piston(
            PistonImpulseError::Cancelled
        ))
    ));
    assert_eq!(polls, 2);
    assert_eq!(wall, before);
    assert_eq!(body, body_before);
    same_ports(&wall, &before);

    let controls = RigidPistonImpactControls {
        piston: PistonImpulseControls {
            maximum_abs_velocity_m_s: 0.1,
            ..PistonImpulseControls::default()
        },
        ..RigidPistonImpactControls::default()
    };
    let result = wall.apply_rigid_normal_impact(
        &mut body,
        properties,
        contact(Vec3::ZERO, 1.0),
        controls,
        || false,
    );
    assert!(matches!(
        result,
        Err(RigidPistonImpactError::Piston(
            PistonImpulseError::Budget { .. }
        ))
    ));
    assert_eq!(wall, before);
    assert_eq!(body, body_before);
    same_ports(&wall, &before);
}

#[test]
fn impact_drives_subsequent_native_gas_compression_without_instant_heat_credit() {
    let mut wall = GasPiston::new(gas(0.75, 0.001), 2.0, 0.01, 0.0).unwrap();
    let gas_before = wall.gas().clone();
    let (properties, mut body) = rigid(&wall, Vec3::ZERO);
    let receipt = wall
        .apply_rigid_normal_impact(
            &mut body,
            properties,
            contact(Vec3::ZERO, 0.0),
            RigidPistonImpactControls::default(),
            || false,
        )
        .unwrap();
    near(receipt.contact.dissipated_energy, 3.0);
    near(wall.velocity_m_s(), -1.0);
    assert_eq!(wall.gas(), &gas_before);
    let step = wall.step(drive(0.0), step_controls(), || false).unwrap();
    assert!(step.displacement_m < 0.0);
    assert!(wall.gas().volume_m3() < gas_before.volume_m3());
    assert!(wall.gas().state().state().temperature > gas_before.state().state().temperature);
    assert!(step.gas.work_by_gas_j < 0.0);
    assert!(step.total_energy_residual_j.abs() <= step_controls().energy_tolerance_j);
}

#[test]
fn prescribed_braking_preserves_original_ideal_and_real_gas_owners() {
    fn stop<G: PistonGas>(gas: G) {
        let mut wall = GasPiston::new(gas, 2.0, 0.01, 2.0).unwrap();
        let before = wall.clone();
        assert!(before.clone().with_fixed_support().is_err());
        let receipt = wall
            .apply_impulse(-4.0, PistonImpulseControls::default(), || false)
            .unwrap();
        near(receipt.represented_momentum_change_n_s, -4.0);
        near(receipt.kinetic_energy_change_j, -4.0);
        near(receipt.external_work_j, -4.0);
        assert_eq!(wall.gas(), before.gas());
        same_ports(&wall, &before);
        let held = wall.with_fixed_support().unwrap();
        assert_eq!(held.velocity_m_s(), 0.0);
    }
    stop(gas(0.75, 0.001));
    stop(real_parcel(300.0, 5e6));
}

#[test]
fn double_acting_impulse_changes_one_shared_mass_and_neither_source_gas() {
    let primary = GasPiston::new(gas(0.75, 0.001), 2.0, 0.01, 0.0).unwrap();
    let mut wall = DoubleActingGasPiston::new(primary, gas(0.25, 0.0), 0.008).unwrap();
    let before = wall.clone();
    let receipt = wall
        .apply_impulse(6.0, PistonImpulseControls::default(), || false)
        .unwrap();
    near(wall.primary().velocity_m_s(), 3.0);
    near(receipt.kinetic_energy_change_j, 9.0);
    assert_eq!(wall.primary().gas(), before.primary().gas());
    assert_eq!(wall.opposing(), before.opposing());
    same_ports(&wall, &before);
}

#[test]
fn common_large_world_translation_does_not_activate_a_separated_contact() {
    fn separated(offset: f64) -> RigidPistonImpactReceipt {
        let mut wall = GasPiston::new(gas(0.75, 0.001), 2.0, 0.01, 0.0).unwrap();
        let before = wall.clone();
        let arm = Vec3::new(0.0, 0.0, 0.15);
        let origin = Vec3::new(offset, offset, offset);
        let properties = MassProperties::new(1.0, Vec3::ZERO, Vec3::new(1.0, 2.0, 3.0)).unwrap();
        let mut body = RigidBodyState::new(
            Pose::new(origin, UnitQuaternion::IDENTITY).unwrap(),
            Vec3::new(0.0, 0.0, -3.0),
            Vec3::ZERO,
        )
        .unwrap();
        let body_before = body;
        let receipt = wall
            .apply_rigid_normal_impact(
                &mut body,
                properties,
                RigidPistonImpactContact {
                    piston_axis_origin_world_m: origin,
                    activation_tolerance_m: 0.01,
                    ..contact(arm, 1.0)
                },
                RigidPistonImpactControls::default(),
                || false,
            )
            .unwrap();
        assert!(!receipt.contact.active);
        assert!(!receipt.contact.impacted);
        assert_eq!(wall, before);
        assert_eq!(body, body_before);
        receipt
    }
    let local = separated(0.0);
    let translated = separated(18_014_398_509_481_984.0); // 2^54: a 0.05 m gap vanishes in absolute points.
    near(local.gap_m, 0.05);
    assert_eq!(translated.gap_m, local.gap_m);
    assert_eq!(translated.contact, local.contact);
    assert_eq!(translated.piston, local.piston);
}
