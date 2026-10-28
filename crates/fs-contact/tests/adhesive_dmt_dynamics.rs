//! G1/G3: actual inertial snap, original source work and physical time refinement.
#[path = "source_point/common.rs"]
mod common;

use fs_contact::normal_patch::dmt::{dynamics::*, fixture::DmtSpringFixture, *};
use fs_tribo::InputAuthority;

fn contact(scale: f64) -> SourceDmtSphereContact {
    let a = common::solid("A", 1000.0, 1e11);
    let b = common::solid("B", 2000.0, 2e11);
    SourceDmtSphereContact::new(
        common::bind(&a, &b, None, 0.01 * scale)
            .normal_state()
            .clone(),
        1e-6 * scale,
        1e-9 * scale,
        DmtEnvelope {
            maximum_indentation_m: 1e-9 * scale,
            half_space_depth_m: 1e-3 * scale,
            layer_thickness_m: 1e-3 * scale,
            maximum_interaction_radius_to_radius: 0.1,
            maximum_interaction_radius_to_depth: 0.1,
            maximum_interaction_radius_to_layer: 0.1,
            maximum_nominal_strain: 0.1,
            maximum_repulsive_pressure_to_yield: 0.3,
            maximum_adhesive_stress_to_yield: 0.1,
            maximum_interaction_range_to_radius: 0.01,
            maximum_tabor_parameter: 0.1,
        },
    )
    .unwrap()
}
fn scales(c: &SourceDmtSphereContact, mass: f64) -> (f64, f64, f64) {
    let z = c.interaction_range_m();
    let p = -c.force_controlled_pull_off_n();
    (z, p, (mass * z / p).sqrt())
}
fn model(scale: f64, spring_ratio: f64, damping_ratio: f64) -> DmtInertialContact {
    let c = contact(scale);
    let mass = 1e-12 * scale.powi(3);
    let (z, p, time) = scales(&c, mass);
    let stiffness = p / z;
    DmtInertialContact::new(
        c,
        mass,
        spring_ratio * stiffness,
        damping_ratio * (mass * stiffness).sqrt(),
        DmtDynamicsBudget {
            maximum_step_s: 0.05 * time,
            maximum_iterations: 80,
            minimum_indentation_m: -20.0 * z,
            maximum_indentation_step_m: 0.1 * z,
            maximum_abs_stage_displacement_m: 20.0 * z,
            maximum_abs_velocity_m_s: 20.0 * z / time,
            maximum_abs_force_n: 100.0 * p,
            force_absolute_tolerance_n: 1e-13 * p,
            force_relative_tolerance: 1e-13,
            energy_absolute_tolerance_j: 1e-12 * p * z,
            energy_relative_tolerance: 1e-12,
        },
    )
    .unwrap()
}
fn close(a: f64, b: f64, scale: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= scale * tolerance,
        "{a:e} != {b:e}, scaled error {:e}",
        (a - b).abs() / scale
    );
}
fn independent_energy(m: &DmtInertialContact, s: DmtDynamicsState) -> f64 {
    0.5 * m.mass_kg() * s.velocity_m_s().powi(2)
        + m.contact()
            .evaluate(s.indentation_m())
            .unwrap()
            .reversible_energy_j
        + 0.5 * m.stiffness_n_per_m() * (s.stage_displacement_m() - s.indentation_m()).powi(2)
}
fn check_step(
    m: &DmtInertialContact,
    before: DmtDynamicsState,
    step: &DmtDynamicsStep,
    dt: f64,
    load: f64,
) {
    let (z, p, time) = scales(m.contact(), m.mass_kg());
    let after = step.state();
    let dx = after.indentation_m() - before.indentation_m();
    let velocity = 0.5 * (before.velocity_m_s() + after.velocity_m_s());
    close(dx, dt * velocity, z, 2e-15);
    close(
        m.mass_kg() * (after.velocity_m_s() - before.velocity_m_s()),
        dt * (load + step.spring_force_n - step.discrete_contact_force_n + step.damping_force_n),
        p * dt,
        2e-11,
    );
    close(
        step.contact_work_j,
        m.contact()
            .evaluate(after.indentation_m())
            .unwrap()
            .reversible_energy_j
            - m.contact()
                .evaluate(before.indentation_m())
                .unwrap()
                .reversible_energy_j,
        p * z,
        1e-13,
    );
    close(
        step.damping_dissipation_j,
        m.damping_n_s_per_m() * velocity * velocity * dt,
        p * z,
        1e-13,
    );
    close(
        independent_energy(m, after) - independent_energy(m, before) + step.damping_dissipation_j,
        step.prescribed_stage_work_j + step.applied_load_work_j,
        p * z,
        2e-12,
    );
    close(
        step.total_energy_j,
        independent_energy(m, after),
        p * z,
        1e-13,
    );
    close(after.time_s() - before.time_s(), dt, time, 2e-13);
    assert_eq!(after.steps(), before.steps() + 1);
    assert_eq!(
        step.contact,
        m.contact().evaluate(after.indentation_m()).unwrap()
    );
    assert_eq!(step.contact.authority, InputAuthority::SyntheticFixture);
    assert_eq!(
        step.contact.normal_state_identity(),
        m.contact().normal_state().identity()
    );
}

#[test]
fn detached_mass_obeys_exact_impulse_and_applied_work() {
    let m = model(1.0, 0.0, 0.0);
    let (z, p, time) = scales(m.contact(), m.mass_kg());
    for load in [-0.25 * p, 0.0, 0.25 * p] {
        let mut state = m.initialize(-4.0 * z, -0.2 * z / time, 0.0).unwrap();
        let start = state;
        for _ in 0..100 {
            let dt = 0.001 * time;
            let step = m.advance(state, dt, 0.0, load).unwrap();
            check_step(&m, state, &step, dt, load);
            assert_eq!(step.contact.regime, DmtRegime::Separated);
            assert_eq!(step.discrete_contact_force_n, 0.0);
            assert_eq!(step.damping_dissipation_j, 0.0);
            state = step.state();
        }
        close(
            state.velocity_m_s(),
            start.velocity_m_s() + load / m.mass_kg() * state.time_s(),
            z / time,
            1e-11,
        );
        close(
            state.indentation_m(),
            start.indentation_m()
                + start.velocity_m_s() * state.time_s()
                + 0.5 * load / m.mass_kg() * state.time_s().powi(2),
            z,
            1e-12,
        );
    }
}

#[test]
fn attachment_and_pull_off_resolve_into_kinetic_energy_without_invented_heat() {
    let m = model(1.0, 0.25, 0.0);
    let (z, p, time) = scales(m.contact(), m.mass_kg());
    for (initial, stage, expected_regime) in [
        (-0.999 * z, -z, DmtRegime::Contact),
        (0.0, -4.1 * z, DmtRegime::Separated),
    ] {
        let mut state = m.initialize(initial, 0.0, stage).unwrap();
        let initial_energy = independent_energy(&m, state);
        let mut maximum_kinetic = 0.0_f64;
        let mut reached = false;
        let mut real_contact = false;
        for _ in 0..2400 {
            let dt = 0.005 * time;
            let step = m.advance(state, dt, stage, 0.0).unwrap();
            check_step(&m, state, &step, dt, 0.0);
            assert_eq!(step.damping_dissipation_j, 0.0);
            assert_eq!(step.prescribed_stage_work_j, 0.0);
            assert_eq!(step.applied_load_work_j, 0.0);
            maximum_kinetic = maximum_kinetic.max(step.kinetic_energy_j);
            reached |= step.contact.regime == expected_regime;
            real_contact |= step.contact.contact_area_m2 > 0.0;
            let shear = step.contact.friction_contact_state().unwrap();
            if step.contact.regime != DmtRegime::Contact {
                assert_eq!(shear.area_m2(), 0.0);
                assert_eq!(shear.repulsive_normal_force_n(), 0.0);
            }
            state = step.state();
            close(independent_energy(&m, state), initial_energy, p * z, 2e-10);
        }
        assert!(
            reached,
            "the actual trajectory must complete the instability"
        );
        assert!(
            maximum_kinetic > 0.1 * p * z,
            "snap must release finite kinetic energy"
        );
        if expected_regime == DmtRegime::Contact {
            assert!(real_contact);
        }
    }
}

#[test]
fn physical_damping_settles_snap_to_original_finite_spring_equilibrium() {
    let m = model(1.0, 0.25, 0.8);
    let (z, p, time) = scales(m.contact(), m.mass_kg());
    let stage = -0.9 * z;
    let equilibrium = DmtSpringFixture::new(m.contact().clone(), m.stiffness_n_per_m())
        .unwrap()
        .evaluate_contact(stage)
        .unwrap();
    let mut state = m.initialize(-z, 0.0, stage).unwrap();
    let initial = independent_energy(&m, state);
    let mut loss = 0.0;
    for _ in 0..10000 {
        let dt = 0.005 * time;
        let step = m.advance(state, dt, stage, 0.0).unwrap();
        check_step(&m, state, &step, dt, 0.0);
        assert!(step.damping_dissipation_j >= 0.0);
        loss += step.damping_dissipation_j;
        state = step.state();
    }
    assert!(loss > 0.1 * p * z);
    close(independent_energy(&m, state) + loss, initial, p * z, 1e-9);
    close(
        state.indentation_m(),
        equilibrium.contact.indentation_m,
        z,
        1e-8,
    );
    close(state.velocity_m_s(), 0.0, z / time, 1e-7);
}

#[test]
fn moving_stage_and_held_load_have_distinct_conjugate_work_ports() {
    let m = model(1.0, 2.0, 0.2);
    let (z, p, time) = scales(m.contact(), m.mass_kg());
    let mut state = m.initialize(-2.0 * z, 0.0, -2.0 * z).unwrap();
    let initial = independent_energy(&m, state);
    let (mut stage_work, mut load_work, mut loss) = (0.0, 0.0, 0.0);
    let mut attached = false;
    for i in 1..=1600 {
        let dt = 0.0025 * time;
        let stage = (-2.0 + 0.6 * i as f64 * dt / time) * z;
        let load = 0.07 * p * (i as f64 * dt / time).sin();
        let step = m.advance(state, dt, stage, load).unwrap();
        check_step(&m, state, &step, dt, load);
        let extension0 = state.stage_displacement_m() - state.indentation_m();
        let extension1 = stage - step.state().indentation_m();
        close(
            step.prescribed_stage_work_j,
            0.5 * m.stiffness_n_per_m()
                * (extension0 + extension1)
                * (stage - state.stage_displacement_m()),
            p * z,
            1e-13,
        );
        close(
            step.applied_load_work_j,
            load * (step.state().indentation_m() - state.indentation_m()),
            p * z,
            1e-13,
        );
        attached |= step.contact.contact_area_m2 > 0.0;
        stage_work += step.prescribed_stage_work_j;
        load_work += step.applied_load_work_j;
        loss += step.damping_dissipation_j;
        state = step.state();
    }
    assert!(attached);
    assert!(stage_work.abs() > 0.01 * p * z && load_work.abs() > 0.01 * p * z && loss > 0.0);
    close(
        independent_energy(&m, state) - initial + loss,
        stage_work + load_work,
        p * z,
        1e-9,
    );
}

#[test]
fn unstable_gap_motion_converges_quadratically_to_the_analytic_trajectory() {
    let m = model(1.0, 0.25, 0.0);
    let (z, p, time) = scales(m.contact(), m.mass_kg());
    let stage = -2.0 * z;
    let x0 = -0.65 * z;
    let v0 = 0.01 * z / time;
    let growth = ((p / z - m.stiffness_n_per_m()) / m.mass_kg()).sqrt();
    let unstable = -(p + m.stiffness_n_per_m() * stage) / (p / z - m.stiffness_n_per_m());
    let end = 0.5 * time;
    let expected_x =
        unstable + (x0 - unstable) * (growth * end).cosh() + v0 / growth * (growth * end).sinh();
    let expected_v = growth * (x0 - unstable) * (growth * end).sinh() + v0 * (growth * end).cosh();
    let mut errors = Vec::new();
    for count in [20, 40, 80] {
        let mut state = m.initialize(x0, v0, stage).unwrap();
        for _ in 0..count {
            let step = m.advance(state, end / count as f64, stage, 0.0).unwrap();
            assert_eq!(step.contact.regime, DmtRegime::AttractiveGap);
            state = step.state();
        }
        errors.push(
            (state.indentation_m() - expected_x).abs() / z
                + (state.velocity_m_s() - expected_v).abs() / (z / time),
        );
    }
    assert!(
        errors[0] / errors[1] > 3.9 && errors[1] / errors[2] > 3.9,
        "{errors:?}"
    );
    assert!(errors[2] < 1e-7, "{errors:?}");
}

#[test]
fn source_preserving_dynamics_scales_with_length_mass_force_and_time() {
    let reference = model(1.0, 0.25, 0.2);
    let (z, p, time) = scales(reference.contact(), reference.mass_kg());
    for length in [0.01, 10.0, 1000.0] {
        let scaled = model(length, 0.25, 0.2);
        let mut a = reference.initialize(-0.9 * z, 0.03 * z / time, -z).unwrap();
        let mut b = scaled
            .initialize(-0.9 * z * length, 0.03 * z / time, -z * length)
            .unwrap();
        for i in 1..=80 {
            let stage = (-1.0 + 0.002 * i as f64) * z;
            let x = reference.advance(a, 0.005 * time, stage, 0.02 * p).unwrap();
            let y = scaled
                .advance(
                    b,
                    0.005 * time * length,
                    stage * length,
                    0.02 * p * length.powi(2),
                )
                .unwrap();
            close(
                y.state().indentation_m() / length,
                x.state().indentation_m(),
                z,
                2e-10,
            );
            close(
                y.state().velocity_m_s(),
                x.state().velocity_m_s(),
                z / time,
                2e-10,
            );
            close(
                y.contact.contact_area_m2 / length.powi(2),
                x.contact.contact_area_m2,
                reference.contact().effective_radius_m() * z,
                2e-10,
            );
            close(
                y.total_energy_j / length.powi(3),
                x.total_energy_j,
                p * z,
                2e-10,
            );
            close(
                y.damping_dissipation_j / length.powi(3),
                x.damping_dissipation_j,
                p * z,
                2e-10,
            );
            a = x.state();
            b = y.state();
        }
    }
}

#[test]
fn refusals_and_immutable_successors_preserve_history_without_static_branch_rules() {
    let m = model(1.0, 0.25, 0.0);
    let (z, p, time) = scales(m.contact(), m.mass_kg());
    let initial = m.initialize(-0.9 * z, 0.0, -z).unwrap();
    let step = m.advance(initial, 0.005 * time, -z, 0.0).unwrap();
    assert_eq!(step, m.advance(initial, 0.005 * time, -z, 0.0).unwrap());
    let mut altered = step;
    altered.contact.indentation_m = 1e9;
    assert_eq!(
        altered.state(),
        step.state(),
        "diagnostics cannot mutate accepted physics"
    );
    assert!(m.advance(initial, 0.1 * time, -z, 0.0).is_err());
    assert!(m.advance(initial, 0.005 * time, -z, 101.0 * p).is_err());
    assert!(
        model(1.0, 0.5, 0.0)
            .advance(initial, 0.005 * time, -z, 0.0)
            .is_err()
    );
    let mut budget = m.budget();
    budget.maximum_iterations = 1;
    let one = DmtInertialContact::new(
        m.contact().clone(),
        m.mass_kg(),
        m.stiffness_n_per_m(),
        0.0,
        budget,
    )
    .unwrap();
    assert!(matches!(
        one.advance(initial, 0.005 * time, -z, 0.0),
        Err(DmtDynamicsError::NoConvergence { .. })
    ));
    budget.maximum_iterations = 80;
    budget.maximum_step_s = 4.0 * time;
    let large =
        DmtInertialContact::new(m.contact().clone(), m.mass_kg(), 0.0, 0.0, budget).unwrap();
    let s = large.initialize(-0.5 * z, 0.0, 0.0).unwrap();
    assert!(matches!(
        large.advance(s, 3.0 * time, 0.0, 0.0),
        Err(DmtDynamicsError::NonuniqueStep)
    ));
    budget.maximum_indentation_step_m = 1e-8 * z;
    let narrow = DmtInertialContact::new(
        m.contact().clone(),
        m.mass_kg(),
        m.stiffness_n_per_m(),
        0.0,
        budget,
    )
    .unwrap();
    assert_eq!(
        narrow.advance(initial, 0.005 * time, -z, 0.0),
        Err(DmtDynamicsError::OutsideStepRange)
    );
    // K=S refuses as a massless equilibrium fixture, but finite-mass evolution
    // remains well posed. No massless jump or neutral-branch rule is imported.
    let neutral = model(1.0, 1.0, 0.0);
    let s = neutral.initialize(-0.5 * z, 0.01 * z / time, -z).unwrap();
    neutral.advance(s, 0.005 * time, -z, 0.0).unwrap();
    assert_eq!(initial.steps(), 0);
    assert_eq!(step, m.advance(initial, 0.005 * time, -z, 0.0).unwrap());
}
