//! Every mode's response from rest to a unit held force, as a step's gains.
//!
//! `ModalAcousticTimeModel::unit_force_response_into` is formed by the code that forms the
//! mobility `held_force_for_port_velocity` inverts. These tests pin the two to the bit, at the
//! full and half sample and at an off-clock duration, for a bank of damping ratios and one of
//! explicit damping rates; they check that the response is the change a held force makes to a
//! step from a state that is not rest; and they check the refusals.
//! Run with `cargo test -p fs-modal --test unit_force_response`.

use fs_math::c64::C64;
use fs_modal::time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeError,
    ModalAcousticTimeModel, ModalAcousticWorkspace,
};

const RATE_HZ: u32 = 48_000;

fn mode(omega: f64, ratio: f64) -> ModalAcousticMode {
    ModalAcousticMode {
        angular_frequency_rad_s: omega,
        damping_ratio: ratio,
        pressure_per_modal_velocity: C64::ZERO,
    }
}

fn state(q: f64, v: f64) -> ModalAcousticState {
    ModalAcousticState {
        displacement_m_sqrt_kg: q,
        velocity_m_sqrt_kg_per_s: v,
    }
}

fn budget() -> ModalAcousticTimeBudget {
    ModalAcousticTimeBudget {
        maximum_abs_displacement_m_sqrt_kg: 100.0,
        maximum_total_energy_j: 1.0e20,
        ..ModalAcousticTimeBudget::audible_reference()
    }
}

/// Damping ratios: lightly damped, undamped, strongly damped.
fn ratios() -> ModalAcousticTimeModel {
    let tau = core::f64::consts::TAU;
    ModalAcousticTimeModel::try_new(
        RATE_HZ,
        vec![
            mode(tau * 110.0, 0.01),
            mode(tau * 440.0, 0.0),
            mode(tau * 1500.0, 0.3),
        ],
        budget(),
    )
    .expect("ratio bank")
}

/// Explicit damping rates: a damped free coordinate, a damped mode, an undamped one.
fn rates() -> ModalAcousticTimeModel {
    let tau = core::f64::consts::TAU;
    ModalAcousticTimeModel::try_new_with_damping_rates(
        RATE_HZ,
        vec![
            mode(0.0, 0.0),
            mode(tau * 300.0, 0.0),
            mode(tau * 2000.0, 0.0),
        ],
        &[2.5, 0.7, 0.0],
        budget(),
    )
    .expect("rate bank")
}

fn durations(model: &ModalAcousticTimeModel) -> [f64; 3] {
    let full = model.sample_period_s();
    [full, 0.5 * full, 0.37 * full]
}

const SHAPES: [f64; 3] = [0.4, -1.3, 0.7];
const TARGET_M_S: f64 = 0.21;

#[test]
fn the_response_is_the_mobility_the_held_force_solve_inverts_to_the_bit() {
    for model in [ratios(), rates()] {
        for duration_s in durations(&model) {
            let mut responses = vec![ModalAcousticState::default(); 3];
            model
                .unit_force_response_into(duration_s, &mut responses)
                .expect("response");
            // The solve's own sum, in its order; at rest the free motion is zero.
            let mut mobility = 0.0_f64;
            for (phi, response) in SHAPES.iter().zip(&responses) {
                mobility += phi * phi * response.velocity_m_sqrt_kg_per_s;
            }
            let solved = model
                .held_force_for_port_velocity(&SHAPES, TARGET_M_S, duration_s)
                .expect("solve");
            assert_eq!(
                solved.to_bits(),
                (TARGET_M_S / mobility).to_bits(),
                "duration {duration_s:e}: solve {solved:e}, response {:e}",
                TARGET_M_S / mobility
            );
        }
    }
}

#[test]
fn from_a_moving_state_the_solve_inverts_the_same_response_at_the_full_sample() {
    for mut model in [ratios(), rates()] {
        let moving = [
            state(2.0e-4, -0.3),
            state(-1.0e-5, 0.05),
            state(3.0e-6, 0.12),
        ];
        model.restore_states(&moving).expect("restore");
        let full = model.sample_period_s();
        let mut responses = vec![ModalAcousticState::default(); 3];
        model
            .unit_force_response_into(full, &mut responses)
            .expect("response");
        let (mut free, mut mobility) = (0.0_f64, 0.0_f64);
        for (index, (phi, response)) in SHAPES.iter().zip(&responses).enumerate() {
            let unforced = model.advance_exact_zoh_step(index, moving[index], 0.0);
            free += phi * unforced.velocity_m_sqrt_kg_per_s;
            mobility += phi * phi * response.velocity_m_sqrt_kg_per_s;
        }
        let solved = model
            .held_force_for_port_velocity(&SHAPES, TARGET_M_S, full)
            .expect("solve");
        assert_eq!(solved.to_bits(), ((TARGET_M_S - free) / mobility).to_bits());
        assert_eq!(model.states(), &moving, "neither call moves the state");
    }
}

#[test]
fn the_response_is_what_a_held_force_adds_to_a_step_from_a_moving_state() {
    let forces = [0.8, -2.5, 1.7];
    for model in [ratios(), rates()] {
        let moving = [
            state(2.0e-4, -0.3),
            state(-1.0e-5, 0.05),
            state(3.0e-6, 0.12),
        ];
        for duration_s in durations(&model) {
            let mut responses = vec![ModalAcousticState::default(); 3];
            model
                .unit_force_response_into(duration_s, &mut responses)
                .expect("response");
            let step = |held: &[f64]| {
                let mut stepped = model.clone();
                stepped.restore_states(&moving).expect("restore");
                let mut workspace = ModalAcousticWorkspace::new(&stepped);
                stepped
                    .step_duration_into(held, duration_s, &mut workspace)
                    .expect("step");
                stepped.states().to_vec()
            };
            let (forced, free) = (step(&forces), step(&[0.0; 3]));
            for k in 0..3 {
                for (end, unforced, gain) in [
                    (
                        forced[k].displacement_m_sqrt_kg,
                        free[k].displacement_m_sqrt_kg,
                        responses[k].displacement_m_sqrt_kg,
                    ),
                    (
                        forced[k].velocity_m_sqrt_kg_per_s,
                        free[k].velocity_m_sqrt_kg_per_s,
                        responses[k].velocity_m_sqrt_kg_per_s,
                    ),
                ] {
                    // The change rounds at the size of the end points.
                    let allowed = 8.0 * f64::EPSILON * (end.abs() + unforced.abs());
                    let off = ((end - unforced) - forces[k] * gain).abs();
                    assert!(
                        off <= allowed,
                        "duration {duration_s:e}, mode {k}: change {:e} against {:e}",
                        end - unforced,
                        forces[k] * gain
                    );
                }
            }
        }
    }
}

#[test]
fn a_duration_or_a_slice_that_does_not_fit_is_refused_and_nothing_moves() {
    let model = ratios();
    let mut responses = vec![state(9.0, 9.0); 3];
    for bad in [0.0, -1.0e-3, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            model.unit_force_response_into(bad, &mut responses),
            Err(ModalAcousticTimeError::InvalidInput { .. })
        ));
    }
    assert_eq!(
        responses,
        vec![state(9.0, 9.0); 3],
        "a refusal writes nothing"
    );
    let mut short = vec![ModalAcousticState::default(); 2];
    assert_eq!(
        model.unit_force_response_into(model.sample_period_s(), &mut short),
        Err(ModalAcousticTimeError::StateCountMismatch {
            expected: 3,
            found: 2
        })
    );
    assert!(
        model
            .states()
            .iter()
            .all(|s| *s == ModalAcousticState::default())
    );
}
