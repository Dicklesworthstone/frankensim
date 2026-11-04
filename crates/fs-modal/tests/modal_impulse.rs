//! G1/G3 momentum jumps, reciprocal port work and transactional admission.

use fs_math::c64::C64;
use fs_modal::time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeError,
    ModalAcousticTimeModel, ModalAcousticWorkspace,
};

fn state(q: f64, v: f64) -> ModalAcousticState {
    ModalAcousticState {
        displacement_m_sqrt_kg: q,
        velocity_m_sqrt_kg_per_s: v,
    }
}

fn mode(omega: f64, damping: f64, pressure: f64) -> ModalAcousticMode {
    ModalAcousticMode {
        angular_frequency_rad_s: omega,
        damping_ratio: damping,
        pressure_per_modal_velocity: C64::new(pressure, 0.0),
    }
}

fn close(actual: f64, expected: f64) {
    let scale = actual.abs().max(expected.abs()).max(1.0);
    assert!(
        actual.is_finite() && (actual - expected).abs() <= 64.0 * f64::EPSILON * scale,
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

#[test]
fn g1_impulse_is_independent_of_clock_stiffness_and_damping() {
    for sample_rate in [8_000, 48_000, 96_000] {
        for damping in [0.0, 0.1, 1.0, 5.0] {
            let mut model = ModalAcousticTimeModel::try_new(
                sample_rate,
                vec![mode(1_000.0, damping, 2.0)],
                ModalAcousticTimeBudget::audible_reference(),
            )
            .unwrap();
            model.restore_states(&[state(0.003, -0.4)]).unwrap();
            let frame = model.apply_impulse(&[1.2]).unwrap();
            assert_eq!(model.states()[0].displacement_m_sqrt_kg, 0.003);
            close(model.states()[0].velocity_m_sqrt_kg_per_s, 0.8);
            close(frame.input_work_j, 0.24);
            close(frame.total_modal_energy_j, 0.5 * (9.0 + 0.64));
            close(frame.observer_pressure_pa, 1.6);
            assert_eq!(frame.viscous_dissipation_j, 0.0);
            assert_eq!(frame.modal_viscous_losses[0].energy_j, 0.0);
        }
    }
}

#[test]
fn g1_reciprocal_port_impulse_has_the_physical_effective_mass_and_work() {
    let shapes = [0.5, 1.0, -0.2];
    let initial = [state(0.0, -2.0), state(0.001, 0.4), state(-0.002, 0.3)];
    let mut model = ModalAcousticTimeModel::try_new_with_free_coordinates(
        48_000,
        vec![
            mode(0.0, 0.0, 0.0),
            mode(200.0, 0.1, 0.0),
            mode(900.0, 1.0, 0.0),
        ],
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    model.restore_states(&initial).unwrap();
    let mut generalized = model.clone();
    let before_velocity: f64 = shapes
        .iter()
        .zip(initial)
        .map(|(shape, s)| shape * s.velocity_m_sqrt_kg_per_s)
        .sum();
    let impulse = 0.6;
    let frame = model.apply_port_impulse(&shapes, impulse).unwrap();
    let after_velocity: f64 = shapes
        .iter()
        .zip(model.states())
        .map(|(shape, s)| shape * s.velocity_m_sqrt_kg_per_s)
        .sum();
    close(after_velocity - before_velocity, impulse * 1.29);
    close(
        frame.input_work_j,
        impulse * 0.5 * (before_velocity + after_velocity),
    );
    assert!(
        frame.input_work_j < 0.0,
        "this impulse extracts kinetic energy"
    );
    let projected = shapes.map(|shape| shape * impulse);
    let projected_frame = generalized.apply_impulse(&projected).unwrap();
    assert_eq!(model.states(), generalized.states());
    assert_eq!(frame, projected_frame);
}

#[test]
fn g1_free_mass_obeys_impulse_momentum_and_can_be_stopped() {
    let mut model = ModalAcousticTimeModel::try_free_mass(
        48_000,
        4.0,
        0.2,
        3.0,
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    // phi=1/sqrt(m); J=-m*v stops the physical four-kilogram mass.
    let frame = model.apply_port_impulse(&[0.5], -12.0).unwrap();
    assert_eq!(model.states(), &[state(0.4, 0.0)]);
    assert_eq!(frame.input_work_j, -18.0);
    assert_eq!(frame.total_modal_energy_j, 0.0);
    assert_eq!(frame.viscous_dissipation_j, 0.0);
}

#[test]
fn g1_impulse_excites_the_analytic_damped_ringdown() {
    let omega = 800.0;
    let damping = 0.04;
    let impulse = 0.7;
    let duration = 0.001;
    let mut model = ModalAcousticTimeModel::try_new(
        48_000,
        vec![ModalAcousticMode {
            pressure_per_modal_velocity: C64::new(2.0, -0.3),
            ..mode(omega, damping, 0.0)
        }],
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    let mut workspace = ModalAcousticWorkspace::new(&model);
    let state_buffer = model.states().as_ptr();
    let kick = model
        .apply_impulse_into(&[impulse], &mut workspace)
        .unwrap();
    let energy_buffer = kick.modal_energy_j.as_ptr();
    let loss_buffer = kick.modal_viscous_losses.as_ptr();
    let injected_energy = kick.input_work_j;
    assert_eq!(injected_energy, 0.5 * impulse * impulse);

    let frame = model
        .step_duration_into(&[0.0], duration, &mut workspace)
        .unwrap();
    let decay = damping * omega;
    let damped_omega = omega * (1.0 - damping * damping).sqrt();
    let (sine, cosine) = (damped_omega * duration).sin_cos();
    let envelope = impulse * (-decay * duration).exp();
    let displacement = envelope * sine / damped_omega;
    let velocity = envelope * (cosine - decay * sine / damped_omega);
    close(model.states()[0].displacement_m_sqrt_kg, displacement);
    close(model.states()[0].velocity_m_sqrt_kg_per_s, velocity);
    close(
        frame.observer_pressure_pa,
        2.0 * velocity - 0.3 * omega * displacement,
    );
    close(
        frame.total_modal_energy_j + frame.viscous_dissipation_j,
        injected_energy,
    );
    assert!(frame.viscous_dissipation_j > 0.0);
    assert_eq!(model.states().as_ptr(), state_buffer);
    assert_eq!(frame.modal_energy_j.as_ptr(), energy_buffer);
    assert_eq!(frame.modal_viscous_losses.as_ptr(), loss_buffer);
}

#[test]
fn g3_impulse_subdivision_and_inverse_preserve_displacement_and_work() {
    let mut whole = ModalAcousticTimeModel::try_new(
        48_000,
        vec![mode(600.0, 0.3, 0.1)],
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    whole.restore_states(&[state(0.01, -0.5)]).unwrap();
    let mut split = whole.clone();
    let expected = whole.apply_impulse(&[1.0]).unwrap();
    let mut workspace = ModalAcousticWorkspace::new(&split);
    let mut work = 0.0;
    for _ in 0..4 {
        work += split
            .apply_impulse_into(&[0.25], &mut workspace)
            .unwrap()
            .input_work_j;
    }
    assert_eq!(whole.states(), split.states());
    close(work, expected.input_work_j);
    let reverse = split.apply_impulse(&[-1.0]).unwrap();
    assert_eq!(split.states(), &[state(0.01, -0.5)]);
    close(work + reverse.input_work_j, 0.0);
}

#[test]
fn g0_impulse_refusals_leave_the_entire_mode_bank_unchanged() {
    let mut model = ModalAcousticTimeModel::try_new(
        48_000,
        vec![mode(100.0, 0.1, 0.0), mode(200.0, 0.1, 1.0)],
        ModalAcousticTimeBudget {
            maximum_abs_velocity_m_sqrt_kg_per_s: 2.0,
            maximum_total_energy_j: 2.0,
            maximum_abs_pressure_pa: 0.5,
            ..ModalAcousticTimeBudget::audible_reference()
        },
    )
    .unwrap();
    let initial = model.states().to_vec();
    for impulses in [
        [0.1, f64::NAN],
        [0.1, f64::INFINITY],
        [0.1, 3.0],
        [0.1, 0.6],
        [2.0, 0.1],
    ] {
        assert!(model.apply_impulse(&impulses).is_err());
        assert_eq!(model.states(), initial);
    }
    assert!(matches!(
        model.apply_impulse(&[1.0]),
        Err(ModalAcousticTimeError::ImpulseCountMismatch {
            expected: 2,
            found: 1
        })
    ));
    for (shapes, impulse) in [
        ([1.0, f64::NAN], 0.0),
        ([1.0, f64::MAX], 2.0),
        ([1.0, 1.0], f64::NAN),
    ] {
        assert!(model.apply_port_impulse(&shapes, impulse).is_err());
        assert_eq!(model.states(), initial);
    }
    assert!(model.apply_port_impulse(&[1.0], 1.0).is_err());
    let other = ModalAcousticTimeModel::try_free_mass(
        48_000,
        1.0,
        0.0,
        0.0,
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    let mut wrong_workspace = ModalAcousticWorkspace::new(&other);
    assert!(
        model
            .apply_impulse_into(&[0.1, 0.1], &mut wrong_workspace)
            .is_err()
    );
    assert_eq!(model.states(), initial);
}

#[test]
fn g1_small_signed_work_is_not_erased_by_large_stored_elastic_energy() {
    let mut model = ModalAcousticTimeModel::try_new(
        48_000,
        vec![mode(1_000.0, 0.5, 0.0)],
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    model.restore_states(&[state(1.0, 1.0)]).unwrap();
    let frame = model.apply_impulse(&[1.0e-12]).unwrap();
    let after_velocity = model.states()[0].velocity_m_sqrt_kg_per_s;
    assert_eq!(
        frame.input_work_j,
        1.0e-12 * f64::midpoint(1.0, after_velocity)
    );
    assert!(frame.input_work_j > 0.0);
    assert_eq!(frame.total_modal_energy_j - 500_000.5, 0.0);
    assert_eq!(frame.viscous_dissipation_j, 0.0);
}

fn wide_free_model(velocities: &[f64]) -> ModalAcousticTimeModel {
    let mut model = ModalAcousticTimeModel::try_new_with_free_coordinates(
        48_000,
        vec![mode(0.0, 0.0, 0.0); velocities.len()],
        ModalAcousticTimeBudget {
            maximum_abs_displacement_m_sqrt_kg: 1.0e20,
            maximum_abs_velocity_m_sqrt_kg_per_s: 1.0e20,
            maximum_total_energy_j: 1.0e40,
            ..ModalAcousticTimeBudget::audible_reference()
        },
    )
    .unwrap();
    model
        .restore_states(
            &velocities
                .iter()
                .map(|&v| state(0.0, v))
                .collect::<Vec<_>>(),
        )
        .unwrap();
    model
}

#[test]
fn g1_signed_work_retains_endpoint_mean_and_product_tails_before_cancellation() {
    // The first mean is halfway between f64 values in the first case; the
    // means are exact but each product loses two joules in the second case.
    for (large, kick, expected_work) in [
        ((1_u64 << 55) as f64, 8.0, 64.5),
        ((1_u64 << 52) as f64, 6.0, 36.5),
    ] {
        for physical_port in [false, true] {
            let mut model = wide_free_model(&[large, 0.0, large]);
            let initial = model.states().to_vec();
            let impulses = [kick, 1.0, -kick];
            let mut workspace = ModalAcousticWorkspace::new(&model);
            let frame = if physical_port {
                model.apply_port_impulse_into(&impulses, 1.0, &mut workspace)
            } else {
                model.apply_impulse_into(&impulses, &mut workspace)
            }
            .unwrap();
            // v dot J cancels, leaving |J|^2 / 2 with unit modal masses.
            assert_eq!(frame.input_work_j, expected_work);
            assert_eq!(frame.viscous_dissipation_j, 0.0);
            let reverse = model
                .apply_impulse_into(&impulses.map(|j| -j), &mut workspace)
                .unwrap();
            assert_eq!(reverse.input_work_j, -expected_work);
            assert_eq!(model.states(), initial);
        }
    }
}

#[test]
fn g1_small_signed_work_survives_opposite_large_modal_contributions() {
    let large = (1_u64 << 28) as f64;
    for (small_before, small_impulse, expected_work) in [(0.0, 1.0, 0.5), (1.0, -1.0, -0.5)] {
        let mut model = wide_free_model(&[0.0, small_before, large]);
        let frame = model
            .apply_impulse(&[large, small_impulse, -large])
            .unwrap();
        // The exact works in basis order are +2^55, +/-0.5, -2^55 joules.
        assert_eq!(frame.input_work_j, expected_work);
        assert_eq!(frame.viscous_dissipation_j, 0.0);
        assert_eq!(
            model.states(),
            &[
                state(0.0, large),
                state(0.0, small_before + small_impulse),
                state(0.0, 0.0)
            ]
        );
    }
}

#[test]
fn g0_swallowed_momentum_refuses_before_any_mode_is_published() {
    let mut model = wide_free_model(&[0.0, 1.0e16]);
    let initial = model.states().to_vec();
    let mut workspace = ModalAcousticWorkspace::new(&model);
    for kick in [1.0, -1.0] {
        for physical_port in [false, true] {
            let result = if physical_port {
                model.apply_port_impulse_into(&[0.25, kick], 1.0, &mut workspace)
            } else {
                model.apply_impulse_into(&[0.25, kick], &mut workspace)
            };
            assert!(matches!(
                result,
                Err(ModalAcousticTimeError::InvalidInput {
                    what: "modal impulse is below representable velocity resolution",
                })
            ));
            assert_eq!(model.states(), initial);
        }
    }
}

#[test]
fn g0_nonzero_port_projection_cannot_silently_underflow() {
    let mut model = wide_free_model(&[0.0, 0.0]);
    let initial = model.states().to_vec();
    assert!(matches!(
        model.apply_port_impulse(&[1.0, 1.0e-200], 1.0e-200),
        Err(ModalAcousticTimeError::InvalidInput {
            what: "modal impulse projection loses nonzero momentum",
        })
    ));
    assert_eq!(model.states(), initial);
    // A physically zero mode shape remains a legitimate uncoupled coordinate.
    model.apply_port_impulse(&[1.0, 0.0], 1.0e-200).unwrap();
    assert_eq!(model.states(), &[state(0.0, 1.0e-200), state(0.0, 0.0)]);
}

#[test]
fn g0_pressure_limit_applies_before_damping_can_hide_the_jump() {
    let mut model = ModalAcousticTimeModel::try_new(
        2_000,
        vec![mode(700.0, 1.0, 1.0)],
        ModalAcousticTimeBudget {
            maximum_abs_pressure_pa: 0.5,
            ..ModalAcousticTimeBudget::audible_reference()
        },
    )
    .unwrap();
    let mut endpoint_only = model.clone();
    endpoint_only.restore_states(&[state(0.0, 1.0)]).unwrap();
    assert!(
        endpoint_only
            .step(&[0.0])
            .unwrap()
            .observer_pressure_pa
            .abs()
            < 0.5
    );
    assert!(matches!(
        model.apply_impulse(&[1.0]),
        Err(ModalAcousticTimeError::BudgetExceeded {
            what: "absolute observer pressure",
            ..
        })
    ));
    assert_eq!(model.states(), &[state(0.0, 0.0)]);
}
