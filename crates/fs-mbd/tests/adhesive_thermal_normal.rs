//! G1/G3/G4: actual native spin follows the thermally changed normal disk.
#[path = "../../fs-contact/tests/adhesive_temperature/common.rs"]
mod common;

use core::f64::consts::PI;
use fs_contact::interface_binding::adhesive_shear::SourceAdhesiveShear;
use fs_contact::normal_patch::adhesive::temperature::AdhesiveTemperatures;
use fs_exec::CancelGate;
use fs_material::state_point::MaterialPropertySelection;
use fs_mbd::tangential_friction::{
    RigidContactBinding,
    adhesive::{thermal_normal::*, *},
};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_tribo::adhesive::spinning::contact::AdhesiveWrenchConfig;
use fs_tribo::partial_slip::TangentFrame;

const DT: f64 = 0.05;
const INERTIA: [f64; 2] = [0.08, 0.12];

fn temperatures() -> AdhesiveTemperatures {
    AdhesiveTemperatures {
        surface_a_k: 310.0,
        surface_b_k: 290.0,
        interface_k: 315.0,
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}
fn controls() -> RigidAdhesiveControls {
    RigidAdhesiveControls {
        maximum_bodies: 2,
        geometry_tolerance_m: 1e-12,
        velocity_tolerance_m_s: 1e-10,
        angular_velocity_tolerance_rad_s: 1e-9,
        linear_momentum_tolerance_ns: 1e-12,
        angular_momentum_tolerance_nms: 1e-12,
        energy_tolerance_j: 1e-10,
        contact: AdhesiveWrenchConfig {
            maximum_static_iterations: 96,
            maximum_iterations: 256,
            maximum_evaluations: 512,
            force_tolerance_n: 1e-9,
            moment_tolerance_nm: 1e-11,
            displacement_tolerance_m: 1e-12,
            energy_tolerance_j: 1e-12,
            relative_matrix_tolerance: 1e-13,
            maximum_force_n: 100.0,
            maximum_spin_moment_nm: 10.0,
        },
    }
}
fn bodies(spin: f64) -> ([RigidBodyState; 2], [MassProperties; 2]) {
    let properties = std::array::from_fn(|i| {
        MassProperties::new(
            [2.0, 3.0][i],
            Vec3::ZERO,
            Vec3::new(INERTIA[i], INERTIA[i], INERTIA[i]),
        )
        .unwrap()
    });
    let bodies = std::array::from_fn(|i| {
        RigidBodyState::new(
            Pose::new(
                Vec3::new(0.0, 0.0, [1.0, -1.0][i]),
                UnitQuaternion::IDENTITY,
            )
            .unwrap(),
            Vec3::ZERO,
            Vec3::new(0.0, 0.0, if i == 0 { INERTIA[i] * spin } else { 0.0 }),
        )
        .unwrap()
    });
    (bodies, properties)
}
fn binding() -> RigidContactBinding {
    RigidContactBinding {
        body_a: 0,
        body_b: 1,
        arm_a_body: Vec3::new(0.0, 0.0, -1.0),
        arm_b_body: Vec3::new(0.0, 0.0, 1.0),
    }
}
fn frame() -> TangentFrame {
    TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap()
}

#[test]
fn unequal_bulk_and_face_temperatures_change_real_disk_native_torque_and_complete_energy() {
    let f = common::Fixture::new();
    let mut normal = f.jkr();
    let before = normal.evaluate_attached(0.0).unwrap();
    let mut state = before.state();
    let mut shear = SourceAdhesiveShear::resolve(
        normal.normal_state(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let saved_sources = (normal.clone(), state, shear.clone());
    let refreshed = normal
        .refresh_at_held_indentation(state, f.cards(), temperatures())
        .unwrap();
    let parameter_work = refreshed.parameter_work_j();
    let (mut bodies, properties) = bodies(4.0);
    let saved_bodies = bodies;
    let old_patch = shear
        .sliding_patch_jkr(normal.normal_state(), &before)
        .unwrap();
    let candidate = prepare_rigid_adhesive_with_normal_refresh(
        &bodies,
        &properties,
        RigidAdhesiveContact {
            binding: binding(),
            frame: frame(),
            source: old_patch,
        },
        &normal,
        state,
        refreshed.clone(),
        DT,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(bodies, saved_bodies);
    assert_eq!((normal.clone(), state, shear.clone()), saved_sources);

    let boosted_bodies = core::array::from_fn::<_, 2, _>(|i| {
        RigidBodyState::new(
            bodies[i].pose(),
            bodies[i].linear_momentum_world().add(Vec3::new(
                0.0,
                0.0,
                properties[i].mass() * 2.0_f64.powi(30),
            )),
            bodies[i].angular_momentum_body(),
        )
        .unwrap()
    });
    let boosted = prepare_rigid_adhesive_with_normal_refresh(
        &boosted_bodies,
        &properties,
        RigidAdhesiveContact {
            binding: binding(),
            frame: frame(),
            source: old_patch,
        },
        &normal,
        state,
        refreshed.clone(),
        DT,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(
        boosted.report.mechanical.kinetic_before_j,
        boosted.report.mechanical.kinetic_after_j
    );
    assert_eq!(
        boosted.report.represented_storage_change_j,
        candidate.report.represented_storage_change_j
    );
    assert_eq!(
        boosted.report.energy_residual_j,
        candidate.report.energy_residual_j
    );

    // Independent zero-approach JKR relation: a^3 = 2*pi*W*R^2/E*.
    // Each E uses its own bulk temperature; W and tau use the actual face T.
    let reduced = 1.0 / (0.9375 / 0.9e6 + 0.9375 / 2.2e6);
    let radius = (2.0 * PI * 0.55 / reduced).cbrt();
    let moment = (2.0 / 3.0) * 1700.0 * PI * radius.powi(3);
    let inverse_inertia = 1.0 / INERTIA[0] + 1.0 / INERTIA[1];
    let heat = DT * moment * (4.0 - 0.5 * DT * moment * inverse_inertia);
    let report = &candidate.report;
    close(report.normal_after.contact_radius_m, radius, 1e-14);
    assert!(radius < before.contact_radius_m);
    assert!(report.normal_after.normal_force_n < 0.0);
    close(
        report.mechanical.spin_angular_impulse_world_nms.z,
        -DT * moment,
        1e-12,
    );
    close(
        report.mechanical.scaled_mobility_per_kg.unwrap()[2][2],
        radius * radius * inverse_inertia,
        1e-14,
    );
    close(report.mechanical.dissipated_energy_j, heat, 1e-12);
    close(report.parameter_work_j, parameter_work, 1e-14);
    close(report.energy_residual_j, 0.0, 1e-12);
    close(
        report.represented_storage_change_j,
        parameter_work - heat,
        1e-12,
    );
    close(
        candidate.bodies_after[0].angular_momentum_body().z,
        INERTIA[0] * 4.0 - DT * moment,
        1e-12,
    );
    close(
        candidate.bodies_after[1].angular_momentum_body().z,
        DT * moment,
        1e-12,
    );
    assert_eq!(
        candidate.shear_source.resolved().card_identity(),
        shear.resolved().card_identity()
    );
    assert_eq!(
        candidate.shear_source.resolved().state_coordinate("T"),
        Some(315.0)
    );
    assert_eq!(
        candidate.shear_source.resolved().state_coordinate("p"),
        Some(2e5)
    );

    let published = apply_rigid_adhesive_with_normal_refresh(
        &mut bodies,
        &properties,
        binding(),
        frame(),
        &mut normal,
        &mut state,
        &mut shear,
        refreshed,
        DT,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(bodies, candidate.bodies_after);
    assert_eq!(normal, candidate.normal_source);
    assert_eq!(state, candidate.normal_state);
    assert_eq!(shear, candidate.shear_source);
    assert_eq!(published, candidate.report);
    // A second refresh uses the accepted hot owners and still cannot switch cards.
    let again = normal
        .refresh_at_held_indentation(state, f.cards(), temperatures())
        .unwrap();
    assert_eq!(again.parameter_work_j(), 0.0);
    assert_eq!(
        shear.rebind_jkr_temperature(&normal, &again).unwrap(),
        shear
    );
}

#[test]
fn stationary_temperature_change_supplies_signed_normal_work_without_friction_heat() {
    let f = common::Fixture::new();
    let mut normal = f.jkr();
    let mut state = normal.evaluate_attached(0.0).unwrap().state();
    let mut shear = SourceAdhesiveShear::resolve(
        normal.normal_state(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let refreshed = normal
        .refresh_at_held_indentation(state, f.cards(), temperatures())
        .unwrap();
    let (mut bodies, properties) = bodies(0.0);
    let saved = bodies;
    let report = apply_rigid_adhesive_with_normal_refresh(
        &mut bodies,
        &properties,
        binding(),
        frame(),
        &mut normal,
        &mut state,
        &mut shear,
        refreshed,
        DT,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(bodies, saved);
    assert_eq!(report.mechanical.dissipated_energy_j, 0.0);
    assert_eq!(report.mechanical.midpoint_impulse_work_j, 0.0);
    assert_ne!(report.parameter_work_j, 0.0);
    assert_eq!(report.normal_potential_change_j, report.parameter_work_j);
    close(report.energy_residual_j, 0.0, 1e-16);
}

#[test]
fn wrong_history_force_servo_and_cancellation_preserve_every_body_and_source_owner() {
    let f = common::Fixture::new();
    for case in 0..3 {
        let mut normal = f.jkr();
        let mut state = normal.evaluate_attached(0.0).unwrap().state();
        let mut shear = SourceAdhesiveShear::resolve(
            normal.normal_state(),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let refreshed = if case == 1 {
            let refresh = normal
                .refresh_at_held_load(0.0, f.cards(), temperatures())
                .unwrap();
            state = refresh.before().state();
            refresh
        } else {
            normal
                .refresh_at_held_indentation(state, f.cards(), temperatures())
                .unwrap()
        };
        if case == 0 {
            state = normal.evaluate_attached(1e-4).unwrap().state();
        }
        let gate = CancelGate::new();
        if case == 2 {
            gate.request();
        }
        let (mut bodies, properties) = bodies(4.0);
        let saved = (bodies, normal.clone(), state, shear.clone());
        assert!(
            apply_rigid_adhesive_with_normal_refresh(
                &mut bodies,
                &properties,
                binding(),
                frame(),
                &mut normal,
                &mut state,
                &mut shear,
                refreshed,
                DT,
                controls(),
                &gate,
            )
            .is_err()
        );
        assert_eq!((bodies, normal, state, shear), saved);
    }
}
