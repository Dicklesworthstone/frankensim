//! G0 physical moment reduction, G1 native trajectory limits, G3 origin/boost
//! covariance, and G4 bounded refusal/cancellation. No interval-energy claim.

use super::*;
use crate::UnitQuaternion;
use crate::rocket::rigid::tensor::TensorRocketInertia;
use crate::rocket::{RocketForcing, RocketState};

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.16e} != {b:.16e}");
}

fn vector_close(a: Vec3, b: Vec3, tolerance: f64) {
    close(a.x, b.x, tolerance);
    close(a.y, b.y, tolerance);
    close(a.z, b.z, tolerance);
}

fn controls(step: f64) -> RigidRocketSolve {
    RigidRocketSolve {
        maximum_substep_s: step,
        maximum_substeps: 10_000,
    }
}

fn fixture(
    propellant: f64,
    propellant_inertia: SpecificInertia,
    omega: Vec3,
    orientation: UnitQuaternion,
) -> TensorRocketState {
    let translation = RocketState {
        position_m: Vec3::new(3.0, -5.0, 2.0),
        velocity_m_s: Vec3::new(0.2, -0.3, 0.4),
        dry_mass_kg: 10.0,
        propellant_mass_kg: propellant,
    };
    let inertia = TensorRocketInertia::new(
        Mat3 {
            m: [2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0],
        },
        propellant_inertia,
    )
    .unwrap();
    let spin = inertia
        .mass_properties(translation)
        .unwrap()
        .angular_momentum_body_checked(omega)
        .unwrap();
    TensorRocketState::new(translation, orientation, spin, inertia).unwrap()
}

fn pair(arm: Vec3, positive: Vec3, negative: Vec3) -> BalancedNozzleGroup {
    BalancedNozzleGroup::OpposedPair {
        mass_flow_per_nozzle_kg_s: 0.5,
        arm_body_m: arm,
        positive_exhaust_velocity_body_m_s: positive,
        negative_exhaust_velocity_body_m_s: negative,
    }
}

#[test]
fn g0_same_physical_ports_determine_force_couple_full_inertia_and_power() {
    let nozzles = PreparedNozzleArray::new(
        &[
            BalancedNozzleGroup::OpposedPair {
                mass_flow_per_nozzle_kg_s: 2.0,
                arm_body_m: Vec3::new(1.0, 2.0, 3.0),
                positive_exhaust_velocity_body_m_s: Vec3::new(0.0, 5.0, 0.0),
                negative_exhaust_velocity_body_m_s: Vec3::new(0.0, -5.0, 0.0),
            },
            BalancedNozzleGroup::Central {
                mass_flow_kg_s: 1.0,
                exhaust_velocity_body_m_s: Vec3::new(-10.0, 0.0, 0.0),
            },
        ],
        || false,
    )
    .unwrap();
    close(nozzles.mass_flow_kg_s(), 5.0, 0.0);
    vector_close(
        nozzles.reaction_force_body_n(),
        Vec3::new(10.0, 0.0, 0.0),
        1e-13,
    );
    vector_close(
        nozzles.reaction_couple_body_n_m(),
        Vec3::new(60.0, 0.0, -20.0),
        1e-13,
    );
    close(nozzles.relative_kinetic_power_w(), 100.0, 1e-13);
    let k = nozzles.co_rotating_specific_inertia_m2().tensor_body_m2();
    for (actual, expected) in
        k.m.into_iter()
            .zip([10.4, -1.6, -2.4, -1.6, 8.0, -4.8, -2.4, -4.8, 4.0])
    {
        close(actual, expected, 5e-14);
    }
    let forcing = nozzles.forcing(Vec3::ZERO);
    assert_eq!(forcing.thrust_application_arm_body_m, Vec3::ZERO);
    assert_eq!(forcing.external_torque_body_n_m, Vec3::ZERO);
    vector_close(
        forcing.exhaust_swirl_specific_angular_momentum_body_m2_s,
        Vec3::new(-12.0, 0.0, 4.0),
        1e-13,
    );
}

#[test]
fn g0_g3_absolute_transport_matches_individual_rotating_particles_and_origin_shift() {
    let r = Vec3::new(1.0, 2.0, -0.5);
    let plus = Vec3::new(-4.0, 3.0, 1.0);
    let minus = Vec3::new(2.0, -5.0, 3.0);
    let central = Vec3::new(-2.0, 1.0, -4.0);
    let nozzles = PreparedNozzleArray::new(
        &[
            pair(r, plus, minus),
            BalancedNozzleGroup::Central {
                mass_flow_kg_s: 0.25,
                exhaust_velocity_body_m_s: central,
            },
        ],
        || false,
    )
    .unwrap();
    let omega = Vec3::new(0.4, -0.7, 0.2);
    let orientation = UnitQuaternion::from_axis_angle(Vec3::new(1.0, -2.0, 3.0), 0.7).unwrap();
    let body = fixture(2.0, SpecificInertia::ZERO, omega, orientation);
    let origin = Vec3::new(-6.0, 2.0, 5.0);
    let flux = nozzles.flux_at(&body, origin).unwrap();
    let (mut momentum, mut angular, mut spin) = (Vec3::ZERO, Vec3::ZERO, Vec3::ZERO);
    // Independent physical particle definition, with no prepared forcing/tensor.
    for (rate, arm, velocity) in [
        (0.5, r, plus),
        (0.5, r.scale(-1.0), minus),
        (0.25, Vec3::ZERO, central),
    ] {
        let world_arm = orientation.rotate_body_to_world(arm);
        let world_velocity = body
            .translation()
            .velocity_m_s
            .add(orientation.rotate_body_to_world(omega.cross(arm).add(velocity)));
        let port_momentum = world_velocity.scale(rate);
        momentum = momentum.add(port_momentum);
        spin = spin.add(world_arm.cross(port_momentum));
        angular = angular.add(
            body.translation()
                .position_m
                .add(world_arm)
                .sub(origin)
                .cross(port_momentum),
        );
    }
    vector_close(flux.linear_momentum_rate_world_n, momentum, 2e-13);
    vector_close(flux.spin_angular_momentum_rate_world_n_m, spin, 2e-13);
    vector_close(
        flux.angular_momentum_rate_about_origin_world_n_m,
        angular,
        3e-13,
    );
    let shift = Vec3::new(2.0, -3.0, 1.0);
    let moved = nozzles.flux_at(&body, origin.add(shift)).unwrap();
    vector_close(
        moved.angular_momentum_rate_about_origin_world_n_m,
        angular.sub(shift.cross(momentum)),
        3e-13,
    );

    let boost = Vec3::new(4.0, 1.0, -2.0);
    let boosted = TensorRocketState::new(
        RocketState {
            velocity_m_s: body.translation().velocity_m_s.add(boost),
            ..*body.translation()
        },
        orientation,
        body.angular_momentum_body_kg_m2_s(),
        body.inertia(),
    )
    .unwrap();
    let changed = nozzles.flux_at(&boosted, origin).unwrap();
    vector_close(
        changed.linear_momentum_rate_world_n,
        momentum.add(boost.scale(1.25)),
        2e-13,
    );
    vector_close(changed.spin_angular_momentum_rate_world_n_m, spin, 2e-13);
    close(
        changed.relative_kinetic_power_w,
        flux.relative_kinetic_power_w,
        0.0,
    );
}

#[test]
fn g1_pure_couple_refines_analytic_depleting_rotation_and_closes_at_burnout() {
    let nozzles = PreparedNozzleArray::new(
        &[pair(
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 2.0, 0.0),
            Vec3::new(0.0, -2.0, 0.0),
        )],
        || false,
    )
    .unwrap();
    let body = fixture(
        2.0,
        nozzles.co_rotating_specific_inertia_m2(),
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
    );
    // q=1, j_z=2, I_z=4-t: I*omega'=-2 until t=2.
    // After burnout I=2 and the remaining second is exact free coast.
    let expected_omega = 2.0 * 0.5_f64.ln();
    let powered_angle = 2.0 * (-2.0 * 0.5_f64.ln() - 2.0);
    let expected_angle = powered_angle + expected_omega;
    let mut errors = [0.0; 3];
    for (slot, count) in [32, 64, 128].into_iter().enumerate() {
        let report = body
            .advance_nozzles(
                &nozzles,
                Vec3::ZERO,
                3.0,
                controls(2.0 / f64::from(count)),
                || false,
            )
            .unwrap();
        let final_state = report.state_after;
        let omega = final_state
            .mass_properties()
            .unwrap()
            .angular_velocity_body_checked(final_state.angular_momentum_body_kg_m2_s())
            .unwrap();
        let axis = final_state
            .orientation()
            .rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0));
        errors[slot] = (omega.z - expected_omega).abs()
            + (axis.x - expected_angle.cos()).abs()
            + (axis.y - expected_angle.sin()).abs();
        close(report.burned_propellant_mass_kg, 2.0, 1e-13);
        close(report.powered_duration_s, 2.0, 1e-13);
        close(report.coast_duration_s, 1.0, 1e-13);
        vector_close(
            final_state.translation().velocity_m_s,
            body.translation().velocity_m_s,
            0.0,
        );
        vector_close(report.thrust_impulse_world_n_s, Vec3::ZERO, 0.0);
        vector_close(
            report.applied_angular_impulse_world_kg_m2_s,
            Vec3::ZERO,
            0.0,
        );
        vector_close(
            report.angular_momentum_residual_world_kg_m2_s,
            Vec3::ZERO,
            3e-11,
        );
        assert_eq!(
            nozzles.flux_at(&final_state, Vec3::ZERO).unwrap(),
            NozzleFlux::ZERO
        );
    }
    assert!(errors[0] / errors[1] > 3.5, "{errors:?}");
    assert!(errors[1] / errors[2] > 3.5, "{errors:?}");
    assert!(errors[2] < 8e-5, "{errors:?}");
}

#[test]
fn g1_corotating_off_axis_vent_preserves_spin_speed_when_exit_matches_reservoir() {
    let nozzles = PreparedNozzleArray::new(
        &[pair(Vec3::new(1.0, 2.0, 3.0), Vec3::ZERO, Vec3::ZERO)],
        || false,
    )
    .unwrap();
    let omega = Vec3::new(2.0, -1.0, 0.0); // Perpendicular to the port arm.
    let body = fixture(
        2.0,
        nozzles.co_rotating_specific_inertia_m2(),
        omega,
        UnitQuaternion::IDENTITY,
    );
    let report = body
        .advance_nozzles(&nozzles, Vec3::ZERO, 0.5, controls(0.025), || false)
        .unwrap();
    let result = report.state_after;
    let after = result
        .mass_properties()
        .unwrap()
        .angular_velocity_body_checked(result.angular_momentum_body_kg_m2_s())
        .unwrap();
    vector_close(after, omega, 2e-10);
    vector_close(
        result.translation().velocity_m_s,
        body.translation().velocity_m_s,
        0.0,
    );
    vector_close(
        report.angular_momentum_residual_world_kg_m2_s,
        Vec3::ZERO,
        2e-10,
    );
    assert!(
        report
            .exported_angular_momentum_world_kg_m2_s
            .norm_squared()
            > 1.0
    );
    close(nozzles.relative_kinetic_power_w(), 0.0, 0.0);
}

#[test]
fn g1_central_physical_exhaust_reuses_native_analytic_translation_and_mass_floor() {
    let nozzles = PreparedNozzleArray::new(
        &[BalancedNozzleGroup::Central {
            mass_flow_kg_s: 2.0,
            exhaust_velocity_body_m_s: Vec3::new(-100.0, 0.0, 0.0),
        }],
        || false,
    )
    .unwrap();
    let body = fixture(
        2.0,
        SpecificInertia::ZERO,
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
    );
    let gravity = Vec3::new(0.0, 0.0, -3.0);
    let actual = body
        .advance_nozzles(&nozzles, gravity, 2.0, controls(2.0), || false)
        .unwrap();
    let expected = body
        .translation()
        .advance(
            RocketForcing {
                mass_flow_kg_s: 2.0,
                effective_exhaust_velocity_m_s: Vec3::new(100.0, 0.0, 0.0),
                gravity_m_s2: gravity,
                external_force_n: Vec3::ZERO,
            },
            2.0,
            || false,
        )
        .unwrap();
    vector_close(
        actual.state_after.translation().position_m,
        expected.state_after.position_m,
        1e-12,
    );
    vector_close(
        actual.state_after.translation().velocity_m_s,
        expected.state_after.velocity_m_s,
        1e-12,
    );
    assert_eq!(actual.state_after.translation().propellant_mass_kg, 0.0);
    assert_eq!(actual.state_after.translation().dry_mass_kg, 10.0);
    assert_eq!(actual.substeps, 2);
}

#[test]
fn g4_preparation_is_bounded_cancellable_and_rejects_invalid_inactive_ports() {
    assert_eq!(
        PreparedNozzleArray::new(&[], || false).unwrap(),
        PreparedNozzleArray::ZERO
    );
    let central = BalancedNozzleGroup::Central {
        mass_flow_kg_s: 0.0,
        exhaust_velocity_body_m_s: Vec3::new(f64::NAN, 0.0, 0.0),
    };
    assert!(PreparedNozzleArray::new(&[central], || false).is_err());
    assert!(
        PreparedNozzleArray::new(
            &vec![pair(Vec3::ZERO, Vec3::ZERO, Vec3::ZERO); MAX_NOZZLE_GROUPS + 1],
            || false,
        )
        .is_err()
    );
    let groups = [pair(Vec3::new(1.0, 2.0, 3.0), Vec3::ZERO, Vec3::ZERO); 65];
    let mut polls = 0;
    assert!(matches!(
        PreparedNozzleArray::new(&groups, || {
            polls += 1;
            polls == 5 // During the second pass, after the total rate is known.
        }),
        Err(RigidRocketError::Cancelled)
    ));
    assert_eq!(polls, 5);
    assert!(
        PreparedNozzleArray::new(
            &[BalancedNozzleGroup::Central {
                mass_flow_kg_s: -1.0,
                exhaust_velocity_body_m_s: Vec3::ZERO,
            },],
            || false
        )
        .is_err()
    );

    // A large speed with a small rate has finite physical power; forming its
    // unweighted speed squared first would spuriously overflow.
    let high_speed = PreparedNozzleArray::new(
        &[BalancedNozzleGroup::Central {
            mass_flow_kg_s: 1e-100,
            exhaust_velocity_body_m_s: Vec3::new(1e200, 0.0, 0.0),
        }],
        || false,
    )
    .unwrap();
    assert!((high_speed.relative_kinetic_power_w() / 5e299 - 1.0).abs() < 1e-14);
    let body = fixture(
        2.0,
        SpecificInertia::ZERO,
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
    );
    assert!(matches!(
        body.advance_nozzles(&high_speed, Vec3::ZERO, 0.1, controls(0.01), || true,),
        Err(RigidRocketError::Cancelled)
    ));
}
