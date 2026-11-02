//! Native source-backed normal/shear recoil, work and physical-time refinement.
#[path = "../../fs-contact/tests/source_point/common.rs"]
mod common;
#[path = "native_mindlin/driven.rs"]
mod driven;

use fs_contact::interface_binding::mindlin::loading::{
    SourceHertzMindlinContact, SourceHertzMindlinState,
};
use fs_contact::interface_binding::{NormalContactModelSelection, mindlin::SourceMindlinContact};
use fs_contact::normal_patch::source::{PointContactGeometry, state::PointContactMaterialCards};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    InterfaceSystemCard, InterpolationPolicy, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, UncertaintyModel,
};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
};
use fs_mbd::tangential_friction::mindlin::*;
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_qty::Dims;
use fs_tribo::{
    InputAuthority, TangentialLaw, TangentialSolveConfig, TangentialStepInput,
    partial_slip::TangentFrame,
};

fn source(young: f64, mu: f64) -> SourceHertzMindlinContact {
    let a = common::solid("native-Mindlin-A", 1000.0, young);
    let b = common::solid("native-Mindlin-B", 2000.0, 2e7);
    let normal = common::bind(&a, &b, None, 0.0);
    let original = normal.normal_state().elastic().interface().card();
    let mut claims = original.claims().clone();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("local-friction", Dims::NONE),
            value: PropertyValue::Scalar {
                value: mu,
                dims: Dims::NONE,
            },
            validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic native friction fixture".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let interface = InterfaceSystemCard::assemble(
        original.surface_a().clone(),
        original.surface_b().clone(),
        original.context().clone(),
        claims,
        original.models().to_vec(),
    )
    .unwrap();
    let selection = MaterialPropertySelection::SingleClaimOnly;
    let normal = PointContactMaterialCards {
        surface_a: &a,
        surface_b: &b,
        interface: &interface,
        bulk_selection: [selection.clone(), selection.clone()],
        interface_selection: selection.clone(),
        model_selection: NormalContactModelSelection::SingleSupported,
        authority: InputAuthority::SyntheticFixture,
    }
    .resolve(
        &common::query(),
        PointContactGeometry::SpherePlane {
            effective_radius_m: 0.1,
        },
        common::envelope(),
    )
    .unwrap();
    SourceHertzMindlinContact::from_fixed_source(
        SourceMindlinContact::resolve(
            normal,
            2e-5,
            ScalarPropertyRequirement::try_new(
                "local-friction",
                Dims::NONE,
                ScalarAdmissibility::NonNegative,
            )
            .unwrap(),
            selection,
            64,
        )
        .unwrap(),
    )
    .unwrap()
}

fn controls() -> RigidMindlinControls {
    RigidMindlinControls {
        solve: TangentialSolveConfig {
            max_iterations: 96,
            maximum_force_n: 100.0,
            force_absolute_tolerance_n: 1e-13,
            force_relative_tolerance: 1e-13,
        },
        maximum_normal_iterations: 64,
        maximum_step_s: 1e-3,
        maximum_slip_speed_m_s: 1.0,
        geometry_tolerance_m: 1e-12,
        displacement_tolerance_m: 1e-16,
        velocity_tolerance_m_s: 1e-10,
        angular_velocity_tolerance_rad_s: 1e-10,
        linear_momentum_tolerance_ns: 1e-12,
        angular_momentum_tolerance_nms: 1e-12,
        energy_tolerance_j: 1e-13,
    }
}

fn pair(
    d: f64,
    tangent_speed: f64,
    normal_speed: f64,
) -> ([RigidBodyState; 2], [MassProperties; 2], [Vec3; 2]) {
    let arms = [
        Vec3::new(0.0, 0.2 - 2.0 * d / 3.0, 0.0),
        Vec3::new(0.0, -(0.2 - d / 3.0), 0.0),
    ];
    let properties = [1.0, 2.0].map(|mass| {
        MassProperties::new(
            mass,
            Vec3::ZERO,
            Vec3::new(0.016 * mass, 0.016 * mass, 0.016 * mass),
        )
        .unwrap()
    });
    let bodies = [
        RigidBodyState::new(
            Pose::new(arms[0].scale(-1.0), UnitQuaternion::IDENTITY).unwrap(),
            Vec3::new(tangent_speed, normal_speed, 0.0),
            Vec3::ZERO,
        )
        .unwrap(),
        RigidBodyState::new(
            Pose::new(arms[1].scale(-1.0), UnitQuaternion::IDENTITY).unwrap(),
            Vec3::ZERO,
            Vec3::ZERO,
        )
        .unwrap(),
    ];
    (bodies, properties, arms)
}

fn contact(law: &SourceHertzMindlinContact, arms_body: [Vec3; 2]) -> RigidMindlinContact<'_> {
    RigidMindlinContact {
        law,
        arms_body,
        frame: TangentFrame::new([0.0, 1.0, 0.0], [1.0, 0.0, 0.0]).unwrap(),
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}
fn advance(
    law: &SourceHertzMindlinContact,
    history: &SourceHertzMindlinState,
    bodies: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    arms: [Vec3; 2],
    dt: f64,
) -> RigidMindlinCandidate {
    advance_rigid_mindlin(
        bodies,
        properties,
        contact(law, arms),
        history,
        dt,
        controls(),
        &mut || false,
    )
    .unwrap()
}

#[test]
fn native_normal_and_rotational_shear_recoil_use_the_same_source_work() {
    let law = source(1e7, 0.6);
    let state = law.initial_state();
    let (bodies, properties, arms) = pair(state.indentation_m(), 0.01, 0.005);
    let dt = 1e-3;
    let result = advance(&law, &state, bodies, properties, arms, dt);
    let report = &result.report;
    close(report.mobility_kg_inverse[0], 1.5, 1e-14);
    close(
        report.mobility_kg_inverse[1],
        1.5 + arms[0].y.powi(2) / 0.016 + arms[1].y.powi(2) / 0.032,
        1e-13,
    );
    let q = result.contact.tangential.traction_n[0];
    let n = report.normal_average_force_n;
    assert!(q > 0.0 && n > state.normal_force_n());
    assert!(result.contact.tangential.state.indentation_m() > state.indentation_m());
    assert!(result.contact.tangential.dissipated_j > 0.0);
    for (i, sign) in [(0, -1.0), (1, 1.0)] {
        close(
            result.bodies_after[i].linear_momentum_world().x,
            bodies[i].linear_momentum_world().x + sign * dt * q,
            1e-15,
        );
        close(
            result.bodies_after[i].linear_momentum_world().y,
            bodies[i].linear_momentum_world().y + sign * dt * n,
            1e-15,
        );
        close(
            result.bodies_after[i].angular_momentum_body().z,
            -sign * arms[i].y * dt * q,
            1e-15,
        );
        assert_eq!(result.bodies_after[i].pose(), bodies[i].pose());
    }
    close(report.total_energy_residual_j, 0.0, 1e-13);
    close(
        report.angular_momentum_change_nms.norm_squared(),
        0.0,
        1e-25,
    );
    let mut lighter_rotation = properties;
    lighter_rotation[0] =
        MassProperties::new(1.0, Vec3::ZERO, Vec3::new(0.012, 0.012, 0.008)).unwrap();
    let changed = advance(&law, &state, bodies, lighter_rotation, arms, dt);
    assert!(changed.report.mobility_kg_inverse[1] > report.mobility_kg_inverse[1]);
    assert!(
        changed.contact.tangential.slip_increment_m[0]
            < result.contact.tangential.slip_increment_m[0]
    );
    let stiffer = source(3e7, 0.6);
    let changed = advance(
        &stiffer,
        &stiffer.initial_state(),
        bodies,
        properties,
        arms,
        dt,
    );
    assert!(changed.report.normal_average_force_n > n);
    assert!(
        changed.contact.tangential.state.indentation_m()
            < result.contact.tangential.state.indentation_m()
    );
}

#[test]
fn unloading_preloaded_shear_drives_bodies_without_a_fictitious_heat_source() {
    let law = source(1e7, 0.6);
    let initial = law.initial_state();
    let critical = law.fixed_source().law().gross_slip_displacement_m();
    let loaded = law
        .advance(
            &initial,
            &TangentialStepInput {
                normal_n: initial.normal_force_n(),
                slip_increment_m: [0.9 * critical, 0.0],
                dt_s: 1e-3,
                contact_temperature_k: Some(law.temperature_kelvin()),
                nominal_pressure_pa: None,
                compliance_m_per_n: 0.0,
                solve: controls().solve,
            },
        )
        .unwrap()
        .state;
    let (bodies, properties, arms) = pair(loaded.indentation_m(), 0.0, -0.003);
    let result = advance(&law, &loaded, bodies, properties, arms, 1e-3);
    assert!(result.contact.tangential.traction_n[0] > 0.0);
    assert!(result.contact.tangential.slip_increment_m[0] < 0.0);
    assert!(result.contact.tangential.stored_change_j < 0.0);
    assert!(result.contact.normal_work_j < 0.0);
    assert!(result.contact.tangential.dissipated_j >= 0.0);
    assert!(result.report.kinetic_after_j > result.report.kinetic_before_j);
    close(
        result.report.kinetic_after_j - result.report.kinetic_before_j
            + result.contact.normal_stored_change_j
            + result.contact.tangential.stored_change_j
            + result.contact.tangential.dissipated_j,
        0.0,
        1e-13,
    );
    let relaxed = advance(&law, &initial, bodies, properties, arms, 1e-3);
    assert_eq!(relaxed.contact.tangential.traction_n[0], 0.0);
    assert_eq!(relaxed.bodies_after[0].linear_momentum_world().x, 0.0);
}

#[test]
fn native_normal_motion_has_second_order_time_refinement_against_independent_ode() {
    let law = source(1e7, 0.0);
    let start = law.initial_state();
    let (initial, properties, _) = pair(start.indentation_m(), 0.0, 0.003);
    let stiffness = law.fixed_source().normal_source().stiffness_n_per_m_3_2();
    let total_time = 0.002;
    // Independent classical RK4 solution of m_eff*d'' = -K*d^(3/2).
    // This uses instantaneous Hertz force, not the implementation's path mean.
    let mut reference = [start.indentation_m(), 0.003];
    let step = total_time / 2048.0;
    let rhs = |s: [f64; 2]| [s[1], -1.5 * stiffness * s[0] * s[0].sqrt()];
    let shifted = |s: [f64; 2], k: [f64; 2], h: f64| [s[0] + h * k[0], s[1] + h * k[1]];
    for _ in 0..2048 {
        let a = rhs(reference);
        let b = rhs(shifted(reference, a, 0.5 * step));
        let c = rhs(shifted(reference, b, 0.5 * step));
        let d = rhs(shifted(reference, c, step));
        for i in 0..2 {
            reference[i] += step * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    let mut errors = Vec::new();
    for count in [8, 16, 32] {
        let mut bodies = initial;
        let mut history = start.clone();
        let dt = total_time / count as f64;
        for _ in 0..count {
            let point = bodies[0]
                .pose()
                .position_world()
                .scale(0.5)
                .add(bodies[1].pose().position_world().scale(0.5));
            let arms = bodies.map(|body| point.sub(body.pose().position_world()));
            let result = advance(&law, &history, bodies, properties, arms, dt);
            assert_eq!(result.contact.tangential.dissipated_j, 0.0);
            close(result.report.total_energy_residual_j, 0.0, 1e-13);
            for i in 0..2 {
                let old = result.bodies_after[i];
                bodies[i] = RigidBodyState::new(
                    Pose::new(
                        old.pose()
                            .position_world()
                            .add(result.report.body_displacement_world_m[i]),
                        old.pose().orientation(),
                    )
                    .unwrap(),
                    old.linear_momentum_world(),
                    old.angular_momentum_body(),
                )
                .unwrap();
            }
            history = result.contact.tangential.state;
            let geometric = start.indentation_m() + initial[1].pose().position_world().y
                - initial[0].pose().position_world().y
                - (bodies[1].pose().position_world().y - bodies[0].pose().position_world().y);
            close(geometric, history.indentation_m(), 2e-15);
        }
        let velocity = bodies[0]
            .center_of_mass_velocity_world(properties[0])
            .unwrap()
            .y
            - bodies[1]
                .center_of_mass_velocity_world(properties[1])
                .unwrap()
                .y;
        errors.push(
            (history.indentation_m() - reference[0]).abs()
                + total_time * (velocity - reference[1]).abs(),
        );
    }
    assert!(
        errors[1] < 0.27 * errors[0] && errors[2] < 0.27 * errors[1],
        "{errors:?}"
    );
    assert!(errors[2] < 1e-10, "{errors:?}");
}

#[test]
fn native_endpoint_rates_cannot_hide_outside_source_domains_behind_small_means() {
    let law = source(1e7, 0.6);
    let history = law.initial_state();
    let (bodies, properties, arms) = pair(history.indentation_m(), 0.01, 0.005);
    let ordinary = advance(&law, &history, bodies, properties, arms, 1e-3);
    let mean = ordinary.report.relative_displacement_m[1] / 1e-3;
    let mut bound = controls();
    bound.maximum_slip_speed_m_s = 0.5 * (mean + 0.01);
    assert!(mean < bound.maximum_slip_speed_m_s);
    assert!(
        advance_rigid_mindlin(
            bodies,
            properties,
            contact(&law, arms),
            &history,
            1e-3,
            bound,
            &mut || false
        )
        .is_err()
    );

    // Choose actual native masses giving d1=d0/4 from rest. The normal mean
    // speed is -0.15 m/s, inside the original 0.2 m/s source domain, while
    // the native endpoint reaches -0.30 m/s and must refuse.
    let dt = 1e-4;
    let target = 0.25 * history.indentation_m();
    let normal = law.fixed_source().normal_source();
    let average = normal
        .elastic_path_average(history.indentation_m(), target, dt)
        .unwrap();
    let inverse_mass = 2.0 * (history.indentation_m() - target) / (dt * dt * average.mean_force_n);
    let base_mass = 1.5 / inverse_mass;
    let light = [base_mass, 2.0 * base_mass].map(|mass| {
        MassProperties::new(
            mass,
            Vec3::ZERO,
            Vec3::new(0.016 * mass, 0.016 * mass, 0.016 * mass),
        )
        .unwrap()
    });
    let resting =
        bodies.map(|body| RigidBodyState::new(body.pose(), Vec3::ZERO, Vec3::ZERO).unwrap());
    assert!(matches!(
        advance_rigid_mindlin(
            resting,
            light,
            contact(&law, arms),
            &history,
            dt,
            controls(),
            &mut || false
        ),
        Err(RigidMindlinError::Normal(_))
    ));
}

#[test]
fn unsupported_motion_separation_and_refusals_preserve_native_history() {
    let law = source(1e7, 0.6);
    let history = law.initial_state();
    let (bodies, properties, arms) = pair(history.indentation_m(), 0.01, 0.005);
    let saved = bodies;
    let saved_history = history.clone();
    let mut noncentral = arms;
    for arm in &mut noncentral {
        arm.x = 1e-3;
    }
    assert!(
        advance_rigid_mindlin(
            bodies,
            properties,
            contact(&law, noncentral),
            &history,
            1e-3,
            controls(),
            &mut || false
        )
        .is_err()
    );
    for (linear, angular) in [
        (Vec3::new(0.01, 0.005, 0.01), Vec3::ZERO),
        (bodies[0].linear_momentum_world(), Vec3::new(0.0, 0.01, 0.0)),
    ] {
        let mut moving = bodies;
        moving[0] = RigidBodyState::new(bodies[0].pose(), linear, angular).unwrap();
        assert!(
            advance_rigid_mindlin(
                moving,
                properties,
                contact(&law, arms),
                &history,
                1e-3,
                controls(),
                &mut || false
            )
            .is_err()
        );
    }
    let (opening, _, _) = pair(history.indentation_m(), 0.0, -0.03);
    assert!(matches!(
        advance_rigid_mindlin(
            opening,
            properties,
            contact(&law, arms),
            &history,
            1e-3,
            controls(),
            &mut || false
        ),
        Err(RigidMindlinError::Separation)
    ));
    let mut short = controls();
    short.maximum_normal_iterations = 1;
    assert!(matches!(
        advance_rigid_mindlin(
            bodies,
            properties,
            contact(&law, arms),
            &history,
            1e-3,
            short,
            &mut || false
        ),
        Err(RigidMindlinError::NormalSolve { .. })
    ));
    let mut polls = 0;
    assert!(
        advance_rigid_mindlin(
            bodies,
            properties,
            contact(&law, arms),
            &history,
            1e-3,
            controls(),
            &mut || {
                polls += 1;
                polls == 5
            }
        )
        .is_err()
    );
    let foreign = source(3e7, 0.6);
    assert!(
        advance_rigid_mindlin(
            bodies,
            properties,
            contact(&foreign, arms),
            &history,
            1e-3,
            controls(),
            &mut || false
        )
        .is_err()
    );
    assert_eq!(bodies, saved);
    assert_eq!(history, saved_history);
    assert_eq!(
        advance(&law, &history, bodies, properties, arms, 1e-3),
        advance(&law, &history, bodies, properties, arms, 1e-3)
    );
}
