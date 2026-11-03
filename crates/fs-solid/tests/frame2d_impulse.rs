//! G0/G1/G3: native generalized impulses, finite-force endpoints and work.

use fs_solid::frame2d::FrameMassFormulation;
use fs_solid::{Frame2d, Geometry, Rayleigh, SectionModel};

const UNDAMPED: Rayleigh = Rayleigh { a0: 0.0, a1: 0.0 };
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:.3e}"
    );
}
fn vector_close(actual: &[f64], expected: &[f64], tolerance: f64) {
    assert_eq!(actual.len(), expected.len());
    for (&a, &b) in actual.iter().zip(expected) {
        close(a, b, tolerance);
    }
}
fn section() -> SectionModel {
    SectionModel::Elastic {
        ea: 4000.0,
        ei: 20.0,
    }
}
fn transverse(formulation: FrameMassFormulation, spin: f64) -> Frame2d {
    let mut frame = Frame2d::with_mass_formulation(Geometry::Linear, formulation);
    frame.add_node(0.0, 0.0);
    frame.add_node(2.0, 0.0);
    frame.fix(0, [true; 3]);
    frame.fix(1, [true, false, false]);
    frame
        .add_element_with_rotary_inertia(0, 1, 5, 3.0, spin, &section)
        .unwrap();
    frame
}
fn axial(formulation: FrameMassFormulation) -> Frame2d {
    let mut frame = Frame2d::with_mass_formulation(Geometry::Linear, formulation);
    frame.add_node(0.0, 0.0);
    frame.add_node(1.0, 0.0);
    frame.fix(0, [true; 3]);
    frame.fix(1, [false, true, true]);
    frame.add_element(0, 1, 3, 6.0, &section).unwrap();
    frame.add_mass(1, 2.0);
    frame
}
fn point_mass(mass: f64) -> Frame2d {
    let mut frame = Frame2d::new(Geometry::Linear);
    frame.add_node(0.0, 0.0);
    frame.add_mass(0, mass);
    frame.fix(0, [false, true, true]);
    frame
}

#[test]
fn g1_consistent_impulse_retains_cross_terms_support_reaction_and_signed_work() {
    let mut frame = transverse(FrameMassFormulation::ConsistentReference, 0.0);
    let before_u = frame.displacements().to_vec();
    let before_sections = frame.section_deformations(0).to_vec();
    let before_force = frame.element_forces(0);
    let impulse = [0.0, 0.0, 0.0, 0.0, 1.0, 0.0];
    let event = frame
        .apply_impulse(&impulse, &[0.0; 6], UNDAMPED, 0.01)
        .unwrap();
    // For L=2 and line mass=3, exact Hermite integration gives
    // M_tip=[[78,-22],[-22,8]]/35 and its inverse [[2,5.5],[5.5,19.5]].
    // The independently integrated support rows then give Jy=-1/2, Jtheta=-1/5.
    vector_close(frame.velocities(), &[0.0, 0.0, 0.0, 0.0, 2.0, 5.5], 3e-13);
    vector_close(
        &event.generalized_momentum_change,
        &[0.0, -0.5, -0.2, 0.0, 1.0, 0.0],
        3e-13,
    );
    vector_close(
        &event.support_impulse,
        &[0.0, -0.5, -0.2, 0.0, 0.0, 0.0],
        3e-13,
    );
    for i in 0..6 {
        close(
            event.generalized_momentum_change[i],
            impulse[i] + event.support_impulse[i],
            3e-13,
        );
    }
    close(event.kinetic_energy_before_j, 0.0, 0.0);
    close(event.kinetic_energy_after_j, 1.0, 3e-13);
    close(event.kinetic_energy_change_j, 1.0, 3e-13);
    close(event.impulse_work_j, 1.0, 3e-13);
    assert_eq!(frame.displacements(), before_u);
    assert_eq!(frame.section_deformations(0), before_sections);
    assert_eq!(frame.element_forces(0), before_force);
    assert_eq!(frame.accelerations(), &[0.0; 6]);

    let extraction = [0.0, 0.0, 0.0, 0.0, -0.5, 0.0];
    let removed = frame
        .apply_impulse(&extraction, &[0.0; 6], UNDAMPED, 0.01)
        .unwrap();
    vector_close(frame.velocities(), &[0.0, 0.0, 0.0, 0.0, 1.0, 2.75], 3e-13);
    close(removed.kinetic_energy_change_j, -0.75, 3e-13);
    close(removed.impulse_work_j, -0.75, 3e-13);
    close(removed.kinetic_energy_after_j, 0.25, 3e-13);
}

#[test]
fn g1_finite_post_event_force_and_damping_set_the_right_limit_and_newmark_continuation() {
    let mut frame = axial(FrameMassFormulation::ConsistentReference);
    let preload = [0.0, 0.0, 0.0, 4.0, 0.0, 0.0];
    frame.static_load(&preload, 1).unwrap();
    let before_u = frame.displacements().to_vec();
    let before_sections = frame.section_deformations(0).to_vec();
    let before_force = frame.element_forces(0);
    let damping = Rayleigh { a0: 0.2, a1: 0.001 };
    let impulse = [0.0, 0.0, 0.0, 2.0, 0.0, 0.0];
    let load = [0.0, 0.0, 0.0, 8.0, 0.0, 0.0];
    let h = 0.001;
    let event = frame.apply_impulse(&impulse, &load, damping, h).unwrap();
    let (m, k, c) = (4.0, 4000.0, 4.8);
    let q0 = before_u[3];
    let v0 = 0.5;
    let a0 = (8.0 - k * q0 - c * v0) / m;
    close(frame.velocities()[3], v0, 1e-14);
    close(frame.accelerations()[3], a0, 1e-12);
    close(event.impulse_work_j, 0.5, 1e-14);
    assert_eq!(frame.displacements(), before_u);
    assert_eq!(frame.section_deformations(0), before_sections);
    assert_eq!(frame.element_forces(0), before_force);

    // Independent one-step average-acceleration equation with that actual
    // post-event rate. The finite load never masquerades as an impulse.
    let a1 = (8.0 - k * (q0 + h * v0 + 0.25 * h * h * a0) - c * (v0 + 0.5 * h * a0))
        / (m + 0.5 * h * c + 0.25 * h * h * k);
    let q1 = q0 + h * v0 + 0.25 * h * h * (a0 + a1);
    let v1 = v0 + 0.5 * h * (a0 + a1);
    let history = frame
        .newmark_with_load(&[0.0, 0.0], h, (0.0, 0.0), damping, &load)
        .unwrap();
    close(
        history.kinetic_energy[0],
        event.kinetic_energy_after_j,
        1e-14,
    );
    close(frame.displacements()[3], q1, 2e-15);
    close(frame.velocities()[3], v1, 2e-14);
    close(frame.accelerations()[3], a1, 2e-11);
}

#[test]
fn g1_lumped_section_spin_and_impulse_are_independent_of_continuation_clock() {
    let original = transverse(FrameMassFormulation::Lumped, 1.5);
    let impulse = [1.0, -2.0, 3.0, 0.0, -3.0, 4.5];
    let mut first = original.clone();
    let mut second = original;
    let a = first
        .apply_impulse(&impulse, &[0.0; 6], UNDAMPED, 1e-6)
        .unwrap();
    let b = second
        .apply_impulse(&impulse, &[0.0; 6], UNDAMPED, 0.25)
        .unwrap();
    assert_eq!(a, b);
    assert_eq!(first.velocities(), second.velocities());
    assert_eq!(first.accelerations(), second.accelerations());
    vector_close(first.velocities(), &[0.0, 0.0, 0.0, 0.0, -1.0, 3.0], 1e-14);
    vector_close(&a.support_impulse, &[-1.0, 2.0, -3.0, 0.0, 0.0, 0.0], 1e-14);
    close(a.impulse_work_j, 8.25, 1e-14);
    close(a.kinetic_energy_change_j, 8.25, 1e-14);
}

#[test]
fn g3_consistent_event_rotates_with_reference_geometry_and_generalized_loads() {
    let build = |c: f64, s: f64| {
        let mut frame = Frame2d::with_mass_formulation(
            Geometry::Linear,
            FrameMassFormulation::ConsistentReference,
        );
        frame.add_node(0.0, 0.0);
        frame.add_node(2.0 * c, 2.0 * s);
        frame
            .add_element_with_rotary_inertia(0, 1, 5, 3.0, 0.5, &section)
            .unwrap();
        frame
    };
    let transform = |values: &[f64]| {
        let mut out = values.to_vec();
        for base in [0, 3] {
            out[base] = 0.6 * values[base] - 0.8 * values[base + 1];
            out[base + 1] = 0.8 * values[base] + 0.6 * values[base + 1];
        }
        out
    };
    let j = [0.25, -0.5, 0.125, 1.0, 0.75, -0.25];
    let force = [1.0, -2.0, 0.5, -0.5, 1.5, 0.25];
    let damping = Rayleigh { a0: 0.125, a1: 0.0 };
    let mut original = build(1.0, 0.0);
    let mut rotated = build(0.6, 0.8);
    let a = original.apply_impulse(&j, &force, damping, 0.001).unwrap();
    let b = rotated
        .apply_impulse(&transform(&j), &transform(&force), damping, 0.001)
        .unwrap();
    vector_close(
        rotated.velocities(),
        &transform(original.velocities()),
        2e-12,
    );
    vector_close(
        rotated.accelerations(),
        &transform(original.accelerations()),
        2e-12,
    );
    vector_close(
        &b.generalized_momentum_change,
        &transform(&a.generalized_momentum_change),
        2e-12,
    );
    assert_eq!(a.support_impulse, vec![0.0; 6]);
    assert_eq!(b.support_impulse, vec![0.0; 6]);
    close(a.kinetic_energy_change_j, b.kinetic_energy_change_j, 2e-12);
    close(a.impulse_work_j, b.impulse_work_j, 2e-12);
}

#[test]
fn g0_support_only_events_preserve_fixed_displacement_and_material_state_for_both_masses() {
    for formulation in [
        FrameMassFormulation::Lumped,
        FrameMassFormulation::ConsistentReference,
    ] {
        let mut frame = axial(formulation);
        frame
            .static_load(&[0.0, 0.0, 0.0, 4.0, 0.0, 0.0], 1)
            .unwrap();
        frame.fix(1, [true; 3]);
        let u = frame.displacements().to_vec();
        let section = frame.section_deformations(0).to_vec();
        let force = frame.element_forces(0);
        assert!(u[3] > 0.0);
        let j = [1.0, -2.0, 3.0, -4.0, 5.0, -6.0];
        let event = frame
            .apply_impulse(&j, &[0.0; 6], Rayleigh { a0: 0.5, a1: 0.01 }, 0.01)
            .unwrap();
        assert_eq!(event.support_impulse, j.map(|v| -v));
        assert_eq!(event.generalized_momentum_change, vec![0.0; 6]);
        assert_eq!(event.kinetic_energy_change_j, 0.0);
        assert_eq!(event.impulse_work_j, 0.0);
        assert_eq!(frame.velocities(), &[0.0; 6]);
        assert_eq!(frame.accelerations(), &[0.0; 6]);
        assert_eq!(frame.displacements(), u);
        assert_eq!(frame.section_deformations(0), section);
        assert_eq!(frame.element_forces(0), force);
    }
}

#[test]
fn g0_free_massless_coordinates_refuse_even_without_direct_rotational_impulse() {
    let mut frame = transverse(FrameMassFormulation::Lumped, 0.0);
    let before = format!("{frame:?}");
    let impulse = [0.0, 0.0, 0.0, 0.0, 1.0, 0.0];
    let error = frame
        .apply_impulse(&impulse, &[0.0; 6], UNDAMPED, 0.01)
        .unwrap_err();
    assert!(error.to_string().contains("massless"));
    assert_eq!(format!("{frame:?}"), before);
    frame.fix(1, [false, false, true]);
    let event = frame
        .apply_impulse(&impulse, &[0.0; 6], UNDAMPED, 0.01)
        .unwrap();
    close(event.velocity_change[4], 1.0 / 3.0, 1e-15);
    assert_eq!(frame.velocities()[5], 0.0);
}

#[test]
fn g0_unrepresentable_velocity_jump_and_late_acceleration_failure_are_atomic() {
    let mut frame = point_mass(1.0);
    frame
        .apply_impulse(&[1e20, 0.0, 0.0], &[0.0; 3], UNDAMPED, 0.01)
        .unwrap();
    let before = format!("{frame:?}");
    assert!(
        frame
            .apply_impulse(&[1.0, 0.0, 0.0], &[0.0; 3], UNDAMPED, 0.01)
            .is_err()
    );
    assert_eq!(
        format!("{frame:?}"),
        before,
        "an unresolved jump cannot silently discard momentum"
    );

    let mut tiny = point_mass(1e-300);
    let before = format!("{tiny:?}");
    // The proposed velocity and kinetic work are finite, but post-event
    // acceleration overflows. Even the lazy Rayleigh reference must roll back.
    assert!(
        tiny.apply_impulse(&[0.5e-300, 0.0, 0.0], &[1e100, 0.0, 0.0], UNDAMPED, 0.01)
            .is_err()
    );
    assert_eq!(format!("{tiny:?}"), before);
    for dt in [0.0, f64::NAN] {
        assert!(
            tiny.apply_impulse(&[0.5e-300, 0.0, 0.0], &[0.0; 3], UNDAMPED, dt)
                .is_err()
        );
        assert_eq!(format!("{tiny:?}"), before);
    }
    assert!(
        tiny.apply_impulse(&[f64::INFINITY, 0.0, 0.0], &[0.0; 3], UNDAMPED, 0.01)
            .is_err()
    );
    assert!(
        tiny.apply_impulse(&[0.0; 2], &[0.0; 3], UNDAMPED, 0.01)
            .is_err()
    );
    assert_eq!(format!("{tiny:?}"), before);
}

#[test]
fn g1_small_resolved_event_reports_work_without_subtracting_large_kinetic_energies() {
    let mut frame = point_mass(1.0);
    frame
        .apply_impulse(&[1e8, 0.0, 0.0], &[0.0; 3], UNDAMPED, 0.01)
        .unwrap();
    let jump = 2.0_f64.powi(-26); // One represented velocity ULP here.
    let event = frame
        .apply_impulse(&[jump, 0.0, 0.0], &[0.0; 3], UNDAMPED, 0.01)
        .unwrap();
    let expected = 1e8 * jump + 0.5 * jump * jump;
    close(event.velocity_change[0], jump, 0.0);
    close(event.kinetic_energy_change_j, expected, 4e-16);
    close(event.impulse_work_j, expected, 4e-16);
    assert!(
        (event.kinetic_energy_after_j - event.kinetic_energy_before_j - expected).abs() > 0.1,
        "the fixture must expose absolute-energy cancellation"
    );
}
