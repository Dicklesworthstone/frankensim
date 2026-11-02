//! Prepared original-body response with genuine temperature-dependent junction claims.
#[path = "../../fs-contact/tests/adhesive_temperature/common.rs"]
mod source_fixture;

use core::f64::consts::PI;
use fs_contact::interface_binding::adhesive_shear::SourceAdhesiveShear;
use fs_exec::CancelGate;
use fs_material::state_point::MaterialPropertySelection;
use fs_mbd::tangential_friction::RigidContactBinding;
use fs_mbd::tangential_friction::adhesive::transient::{
    AdhesiveTransientControls, PreparedRigidAdhesiveTransient, RigidAdhesiveTransientError,
    apply_rigid_adhesive_transient, prepare_rigid_adhesive_transient,
};
use fs_mbd::tangential_friction::adhesive::{RigidAdhesiveContact, RigidAdhesiveControls};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_tribo::adhesive::spinning::contact::AdhesiveWrenchConfig;
use fs_tribo::partial_slip::TangentFrame;

fn controls() -> RigidAdhesiveControls {
    RigidAdhesiveControls {
        maximum_bodies: 4,
        geometry_tolerance_m: 1e-12,
        velocity_tolerance_m_s: 1e-10,
        angular_velocity_tolerance_rad_s: 1e-9,
        linear_momentum_tolerance_ns: 1e-11,
        angular_momentum_tolerance_nms: 1e-11,
        energy_tolerance_j: 1e-10,
        contact: AdhesiveWrenchConfig {
            maximum_static_iterations: 64,
            maximum_iterations: 128,
            maximum_evaluations: 1024,
            force_tolerance_n: 1e-10,
            moment_tolerance_nm: 1e-11,
            displacement_tolerance_m: 1e-12,
            energy_tolerance_j: 1e-11,
            relative_matrix_tolerance: 1e-13,
            maximum_force_n: 100.0,
            maximum_spin_moment_nm: 10.0,
        },
    }
}

fn transient() -> AdhesiveTransientControls {
    AdhesiveTransientControls {
        maximum_step_s: 0.02,
        minimum_step_s: 1e-12,
        maximum_attempts: 50_000,
        maximum_evaluations: 200_000,
        scaled_rate_tolerance_m_s: 1e-11,
        scaled_impulse_tolerance_ns: 1e-11,
        scaled_travel_tolerance_m: 1e-11,
        heat_tolerance_j: 1e-11,
        stopping_time_tolerance_s: 1e-9,
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}

struct Pair {
    bodies: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    binding: RigidContactBinding,
    frame: TangentFrame,
}

impl Pair {
    fn new() -> Self {
        let properties = [
            MassProperties::new(2.0, Vec3::ZERO, Vec3::new(0.06, 0.07, 0.08)).unwrap(),
            MassProperties::new(3.0, Vec3::ZERO, Vec3::new(0.09, 0.11, 0.12)).unwrap(),
        ];
        let arms = [Vec3::new(0.0, 0.0, -0.0625), Vec3::new(0.0, 0.0, 0.0625)];
        let bodies = std::array::from_fn(|i| {
            RigidBodyState::new(
                Pose::new(arms[i].scale(-1.0), UnitQuaternion::IDENTITY).unwrap(),
                Vec3::ZERO,
                Vec3::new(
                    0.0,
                    0.0,
                    [4.0, 1.0][i] * properties[i].principal_inertia_body().z,
                ),
            )
            .unwrap()
        });
        Self {
            bodies,
            properties,
            binding: RigidContactBinding {
                body_a: 0,
                body_b: 1,
                arm_a_body: arms[0],
                arm_b_body: arms[1],
            },
            frame: TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap(),
        }
    }
}

#[test]
fn prepared_thermal_trials_change_arrest_and_heat_with_finite_recoil_from_one_snapshot() {
    let fixture = source_fixture::Fixture::new();
    let normal = fixture.jkr();
    let shear = SourceAdhesiveShear::resolve(
        normal.normal_state(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let response = normal.evaluate_attached(0.0).unwrap();
    assert!(response.normal_force_n < 0.0);
    let source = shear
        .sliding_patch_jkr(normal.normal_state(), &response)
        .unwrap();
    let pair = Pair::new();
    let original = pair.bodies;
    let contact = RigidAdhesiveContact {
        binding: pair.binding,
        frame: pair.frame,
        source,
    };
    let a = response.contact_radius_m;
    let area = PI * a * a;
    let mobility = 1.0 / 0.08 + 1.0 / 0.12;
    let initial_spin = 3.0;
    let original_moment = (2.0 / 3.0) * 2000.0 * area * a;
    let dt = 1.05 * initial_spin / (mobility * original_moment);
    let prepared = PreparedRigidAdhesiveTransient::new(
        &pair.bodies,
        &pair.properties,
        contact,
        dt,
        controls(),
        transient(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(prepared.contact_points_world(), [Vec3::ZERO; 2]);
    let mut hot_result = None;
    for (temperature, strength) in [(315.0, 1700.0), (280.0, 2400.0), (300.0, 2000.0)] {
        let refreshed = source.at_temperature("T", temperature).unwrap();
        let candidate = prepared
            .advance_thermal(&refreshed, &CancelGate::new())
            .unwrap();
        let moment = (2.0 / 3.0) * strength * area * a;
        let stop = initial_spin / (mobility * moment);
        let moving = dt.min(stop);
        let impulse = -moment * moving;
        let spin_a = 4.0 + impulse / 0.08;
        let spin_b = 1.0 - impulse / 0.12;
        let heat = moment * (initial_spin * moving - 0.5 * mobility * moment * moving * moving);
        close(
            candidate.bodies_after[0].angular_momentum_body().z,
            0.08 * spin_a,
            1e-11,
        );
        close(
            candidate.bodies_after[1].angular_momentum_body().z,
            0.12 * spin_b,
            1e-11,
        );
        close(candidate.report.spin_rad_s[1], spin_a - spin_b, 1e-9);
        close(
            candidate.report.spin_angular_impulse_on_a_world_nms.z,
            impulse,
            1e-11,
        );
        close(candidate.report.dissipated_energy_j, heat, 1e-10);
        close(candidate.report.kinetic_energy_change_j, -heat, 1e-10);
        let kinetic_before = 0.5 * (0.08 * 4.0_f64.powi(2) + 0.12);
        let kinetic_after = 0.5 * (0.08 * spin_a * spin_a + 0.12 * spin_b * spin_b);
        close(kinetic_after - kinetic_before + heat, 0.0, 1e-12);
        assert_eq!(candidate.report.force_impulse_on_a_world_ns, Vec3::ZERO);
        assert_eq!(
            candidate.report.audit.unwrap().normal_velocity_bound_m_s,
            0.0
        );
        let path = candidate.report.transient.as_ref().unwrap();
        assert!(path.exact_ray);
        if stop < dt {
            let bracket = path.stopping_time_bracket_s.unwrap();
            close(bracket[0], stop, 1e-10);
            close(bracket[1], stop, 1e-10);
            close(heat, 0.5 * initial_spin * initial_spin / mobility, 1e-12);
        } else {
            assert!(path.stopping_time_bracket_s.is_none());
            assert!(candidate.report.spin_rad_s[1] > 0.0);
            hot_result = Some(candidate.clone());
        }
        for i in 0..2 {
            assert_eq!(candidate.bodies_after[i].pose(), pair.bodies[i].pose());
        }
        assert_eq!(refreshed.original(), source);
        assert_eq!(
            refreshed.friction_source().state_coordinate("T"),
            Some(temperature)
        );
        assert_eq!(
            normal
                .normal_state()
                .interface_state()
                .state_coordinate("T"),
            Some(300.0)
        );
    }
    assert_eq!(pair.bodies, original);
    let hot = source.at_temperature("T", 315.0).unwrap();
    assert_eq!(
        prepared.advance_thermal(&hot, &CancelGate::new()).unwrap(),
        hot_result.unwrap()
    );
    let unchanged = prepared.advance(&CancelGate::new()).unwrap();
    assert_eq!(
        unchanged,
        prepared
            .advance_thermal(
                &source.at_temperature("T", 300.0).unwrap(),
                &CancelGate::new()
            )
            .unwrap()
    );
    assert_eq!(
        unchanged,
        prepare_rigid_adhesive_transient(
            &pair.bodies,
            &pair.properties,
            contact,
            dt,
            controls(),
            transient(),
            &CancelGate::new()
        )
        .unwrap()
    );
    let mut accepted = pair.bodies;
    let report = apply_rigid_adhesive_transient(
        &mut accepted,
        &pair.properties,
        contact,
        dt,
        controls(),
        transient(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(accepted, unchanged.bodies_after);
    assert_eq!(report, unchanged.report);
}

#[test]
fn complete_original_patch_identity_cancellation_and_detached_cases_preserve_the_prepared_state() {
    let fixture = source_fixture::Fixture::new();
    let normal = fixture.jkr();
    let shear = SourceAdhesiveShear::resolve(
        normal.normal_state(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let response = normal.evaluate_attached(0.0).unwrap();
    let source = shear
        .sliding_patch_jkr(normal.normal_state(), &response)
        .unwrap();
    let pair = Pair::new();
    let contact = RigidAdhesiveContact {
        binding: pair.binding,
        frame: pair.frame,
        source,
    };
    let prepared = PreparedRigidAdhesiveTransient::new(
        &pair.bodies,
        &pair.properties,
        contact,
        0.1,
        controls(),
        transient(),
        &CancelGate::new(),
    )
    .unwrap();
    let hot = source.at_temperature("T", 315.0).unwrap();
    let baseline = prepared.advance_thermal(&hot, &CancelGate::new()).unwrap();
    let other_response = normal.evaluate_attached(1e-4).unwrap();
    let other_patch = shear
        .sliding_patch_jkr(normal.normal_state(), &other_response)
        .unwrap();
    assert!(matches!(
        prepared.advance_thermal(
            &other_patch.at_temperature("T", 315.0).unwrap(),
            &CancelGate::new()
        ),
        Err(RigidAdhesiveTransientError::Invalid(_))
    ));

    // Identical physical constants and radius do not identify the original card.
    let foreign_fixture = source_fixture::Fixture::with_options(1.0, 1.0, true, 320.0);
    let foreign = foreign_fixture.jkr();
    let foreign_shear = SourceAdhesiveShear::resolve(
        foreign.normal_state(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let foreign_response = foreign.evaluate_attached(0.0).unwrap();
    assert_eq!(foreign_response.contact_radius_m, response.contact_radius_m);
    let foreign_patch = foreign_shear
        .sliding_patch_jkr(foreign.normal_state(), &foreign_response)
        .unwrap();
    assert!(matches!(
        prepared.advance_thermal(
            &foreign_patch.at_temperature("T", 315.0).unwrap(),
            &CancelGate::new()
        ),
        Err(RigidAdhesiveTransientError::Invalid(_))
    ));
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(prepared.advance_thermal(&hot, &cancelled).is_err());
    assert!(prepared.advance(&cancelled).is_err());
    assert_eq!(
        prepared.advance_thermal(&hot, &CancelGate::new()).unwrap(),
        baseline
    );

    let detached = normal
        .advance_monotone(normal.detached_state(-1e-3).unwrap(), -1e-3)
        .unwrap()
        .response;
    let empty = shear
        .sliding_patch_jkr(normal.normal_state(), &detached)
        .unwrap();
    let no_contact = PreparedRigidAdhesiveTransient::new(
        &pair.bodies,
        &pair.properties,
        RigidAdhesiveContact {
            source: empty,
            ..contact
        },
        0.1,
        controls(),
        transient(),
        &CancelGate::new(),
    )
    .unwrap();
    let noop = no_contact.advance(&CancelGate::new()).unwrap();
    assert_eq!(noop.bodies_after, pair.bodies);
    assert!(noop.report.transient.is_none());
    assert_eq!(noop.report.dissipated_energy_j, 0.0);
    assert!(
        no_contact
            .advance_thermal(&hot, &CancelGate::new())
            .is_err()
    );
    assert!(empty.at_temperature("T", 315.0).is_err());
}

#[test]
fn cached_normal_response_bound_uses_the_refreshed_strength_and_keeps_work_budgets() {
    let fixture = source_fixture::Fixture::new();
    let normal = fixture.jkr();
    let shear = SourceAdhesiveShear::resolve(
        normal.normal_state(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let response = normal.evaluate_attached(0.0).unwrap();
    let source = shear
        .sliding_patch_jkr(normal.normal_state(), &response)
        .unwrap();
    let mut pair = Pair::new();
    let arm = Vec3::new(0.2, -0.3, 0.15).scale(response.contact_radius_m);
    pair.binding.arm_a_body = arm;
    let spin = Vec3::new(0.0, 0.0, 4.0);
    pair.bodies[0] = RigidBodyState::new(
        Pose::new(arm.scale(-1.0), UnitQuaternion::IDENTITY).unwrap(),
        spin.cross(arm).scale(-pair.properties[0].mass()),
        pair.bodies[0].angular_momentum_body(),
    )
    .unwrap();
    let dt = 1e-3;
    // For a diagonal inertia and r=(x,y,z), the original normal response to
    // [Fx,Fy,M/a] is [-xz/Iy,-yz/Ix,0]. The axial second body contributes zero.
    let inertia = pair.properties[0].principal_inertia_body();
    let row_norm = (arm.x * arm.z / inertia.y).hypot(arm.y * arm.z / inertia.x);
    let reference_bound = row_norm * 2000.0 * response.contact_area_m2 * 13.0_f64.sqrt() / 3.0 * dt;
    let gates = RigidAdhesiveControls {
        velocity_tolerance_m_s: reference_bound,
        ..controls()
    };
    let contact = RigidAdhesiveContact {
        binding: pair.binding,
        frame: pair.frame,
        source,
    };
    let prepared = PreparedRigidAdhesiveTransient::new(
        &pair.bodies,
        &pair.properties,
        contact,
        dt,
        gates,
        transient(),
        &CancelGate::new(),
    )
    .unwrap();
    let hot = source.at_temperature("T", 315.0).unwrap();
    let cold = source.at_temperature("T", 280.0).unwrap();
    let warm = prepared.advance_thermal(&hot, &CancelGate::new()).unwrap();
    assert!(
        warm.report
            .transient
            .as_ref()
            .unwrap()
            .stopping_time_bracket_s
            .is_none()
    );
    close(
        warm.report.audit.unwrap().normal_velocity_bound_m_s,
        0.85 * reference_bound,
        1e-12 * reference_bound,
    );
    assert!(matches!(
        prepared.advance_thermal(&cold, &CancelGate::new()),
        Err(RigidAdhesiveTransientError::Invalid(
            "curved traction cannot preserve the held normal-motion bound"
        ))
    ));
    assert_eq!(
        prepared.advance_thermal(&hot, &CancelGate::new()).unwrap(),
        warm
    );
    let bounded = PreparedRigidAdhesiveTransient::new(
        &pair.bodies,
        &pair.properties,
        contact,
        dt,
        gates,
        AdhesiveTransientControls {
            maximum_evaluations: 2,
            ..transient()
        },
        &CancelGate::new(),
    )
    .unwrap();
    assert!(bounded.advance_thermal(&hot, &CancelGate::new()).is_err());
}
