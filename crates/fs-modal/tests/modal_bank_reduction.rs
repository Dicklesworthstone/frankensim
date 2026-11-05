//! G1/G3 modal-bank reductions: resolved ringdown, coherent pressure and signed work.

use fs_math::c64::C64;
use fs_modal::time::{
    MAX_TIME_DOMAIN_ACOUSTIC_MODES, ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget,
    ModalAcousticTimeModel, ModalAcousticWorkspace,
};

fn state(q: f64, v: f64) -> ModalAcousticState {
    ModalAcousticState {
        displacement_m_sqrt_kg: q,
        velocity_m_sqrt_kg_per_s: v,
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        actual.is_finite() && (actual - expected).abs() <= 8.0 * f64::EPSILON * expected.abs(),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

#[test]
fn g1_full_modal_bank_preserves_high_q_ringdown_and_independent_heat() {
    // A resolved 4.377 kHz, Q ~= 7232 oscillator. Replicating independent
    // coordinates multiplies energy and heat; it cannot invalidate their motion.
    let mode = ModalAcousticMode {
        angular_frequency_rad_s: 27_501.440_936_882_4,
        damping_ratio: 6.913_320_946_726_21e-5,
        pressure_per_modal_velocity: C64::ZERO,
    };
    let initial = state(1.672_451_219_462_013_7e-6, 0.079_970_662_945_226_45);
    let count = MAX_TIME_DOMAIN_ACOUSTIC_MODES;
    let limits = ModalAcousticTimeBudget::audible_reference();
    let mut single = ModalAcousticTimeModel::try_new(48_000, vec![mode], limits).unwrap();
    let mut bank = ModalAcousticTimeModel::try_new(48_000, vec![mode; count], limits).unwrap();
    single.restore_states(&[initial]).unwrap();
    bank.restore_states(&vec![initial; count]).unwrap();
    let zero_forces = vec![0.0; count];
    let mut workspace = ModalAcousticWorkspace::new(&bank);
    for _ in 0..16 {
        let reference = single.step(&[0.0]).unwrap();
        let frame = bank.step_into(&zero_forces, &mut workspace).unwrap();
        close(
            frame.total_modal_energy_j,
            count as f64 * reference.total_modal_energy_j,
        );
        close(
            frame.viscous_dissipation_j,
            count as f64 * reference.viscous_dissipation_j,
        );
        assert_eq!(frame.input_work_j, 0.0);
        assert!(frame.viscous_dissipation_j > 0.0);
        assert!(bank.states().iter().all(|s| *s == single.states()[0]));
    }
}

#[test]
fn g1_coherent_pressure_retains_quiet_modes_and_quadratures() {
    // Conditioning stress: the quiet pressure may sit either between modes
    // or in the displacement quadrature of a large velocity contribution.
    for transfers in [
        vec![
            C64::new(1.0e16, 0.0),
            C64::new(1.0, 0.0),
            C64::new(-1.0e16, 0.0),
        ],
        vec![C64::new(1.0e16, 1.0), C64::new(-1.0e16, 0.0)],
    ] {
        let modes: Vec<_> = transfers
            .iter()
            .map(|&transfer| ModalAcousticMode {
                angular_frequency_rad_s: 1_024.0,
                damping_ratio: 0.0,
                pressure_per_modal_velocity: transfer,
            })
            .collect();
        let initial = vec![state(1.0 / 1_024.0, 1.0); modes.len()];
        let zero_forces = vec![0.0; modes.len()];
        let limits = ModalAcousticTimeBudget::audible_reference();
        let mut bank = ModalAcousticTimeModel::try_new(48_000, modes.clone(), limits).unwrap();
        bank.restore_states(&initial).unwrap();
        assert_eq!(
            bank.observer_pressure_with_transfers(&transfers).unwrap(),
            1.0
        );
        assert_eq!(
            bank.observer_pressure_with_transfers_about_static_equilibrium(
                &transfers,
                &zero_forces,
            )
            .unwrap(),
            1.0,
        );
        assert_eq!(
            bank.apply_impulse(&zero_forces)
                .unwrap()
                .observer_pressure_pa,
            1.0
        );
        let frame = bank.step(&zero_forces).unwrap();
        let expected = if modes.len() == 3 {
            bank.states()[0].velocity_m_sqrt_kg_per_s
        } else {
            1_024.0 * bank.states()[0].displacement_m_sqrt_kg
        };
        close(frame.observer_pressure_pa, expected);

        // Cancellation cannot erase a real pressure-budget crossing.
        let mut limited = ModalAcousticTimeModel::try_new(
            48_000,
            modes,
            ModalAcousticTimeBudget {
                maximum_abs_pressure_pa: 0.5,
                ..limits
            },
        )
        .unwrap();
        limited.restore_states(&initial).unwrap();
        assert!(limited.apply_impulse(&zero_forces).is_err());
        assert_eq!(limited.states(), initial);
    }
}

#[test]
fn g1_signed_force_work_retains_small_transfer_between_opposite_sources() {
    let mode = ModalAcousticMode {
        angular_frequency_rad_s: 0.0,
        damping_ratio: 0.0,
        pressure_per_modal_velocity: C64::ZERO,
    };
    let limits = ModalAcousticTimeBudget {
        maximum_abs_displacement_m_sqrt_kg: 2.0,
        ..ModalAcousticTimeBudget::audible_reference()
    };
    for linear in [false, true] {
        let mut bank =
            ModalAcousticTimeModel::try_new_with_free_coordinates(1, vec![mode; 3], limits)
                .unwrap();
        bank.restore_states(&[state(0.0, 1.0), state(0.0, 0.0), state(0.0, -1.25)])
            .unwrap();
        // Exact opposing works +/- 0.28125 J surround a nonzero 5e-19 J transfer.
        let forces = [0.25, 1.0e-9, 0.25];
        let frame = if linear {
            bank.step_linear(&forces, &forces).unwrap()
        } else {
            bank.step(&forces).unwrap()
        };
        close(frame.input_work_j, 0.5 * forces[1] * forces[1]);
        assert_eq!(frame.viscous_dissipation_j, 0.0);
    }
}

#[test]
fn g0_energy_admission_counts_small_modes_beside_large_stored_energy() {
    let mode = ModalAcousticMode {
        angular_frequency_rad_s: 1.0,
        damping_ratio: 0.0,
        pressure_per_modal_velocity: C64::ZERO,
    };
    let limits = ModalAcousticTimeBudget {
        maximum_abs_displacement_m_sqrt_kg: 3.0,
        maximum_total_energy_j: 2.0 + 8.0 * f64::EPSILON,
        ..ModalAcousticTimeBudget::audible_reference()
    };
    let mut bank = ModalAcousticTimeModel::try_new(48_000, vec![mode; 65], limits).unwrap();
    let mut forces = vec![2.0_f64.powi(-26); 65];
    forces[0] = 2.0;
    let equilibrium: Vec<_> = forces.iter().map(|&q| state(q, 0.0)).collect();
    // Each small energy rounds away if independently added to 2 J, but the
    // bank contains 2 + 32*epsilon J and exceeds the supplied physical budget.
    assert!(bank.restore_states(&equilibrium).is_err());
    assert!(bank.initialize_static_equilibrium(&forces).is_err());
    assert!(
        bank.states()
            .iter()
            .all(|s| *s == ModalAcousticState::default())
    );
}
