use super::super::tests::{all_states, fixture};
use super::*;
use crate::modal_acoustic_time::ModalAcousticState;

const DECAY: [f64; 2] = [450.0, 850.0];
const GAIN: [f64; 2] = [3.0, -2.0];
const DIRECT: [f64; 2] = [0.4, -0.2];

fn budget() -> CoupledToneBudget {
    CoupledToneBudget {
        max_dimension: 256,
        max_dense_bytes: 128 * 1_048_576,
        max_multiply_terms: 100_000_000_000,
    }
}
fn controls() -> BroadbandModalControls {
    BroadbandModalControls {
        fit_order: 1,
        fit_iterations: 6,
        delay_kernel_half_len: 2,
        max_held_out_rel_error: 1e-8,
    }
}
fn analytic_transfer(mode: usize, omega: f64) -> C64 {
    C64::from_re(DIRECT[mode])
        + C64::from_re(GAIN[mode] * DECAY[mode]) / C64::new(DECAY[mode], -omega)
}
fn samples(delay_s: f64) -> ModalTransferSamples {
    let omega = vec![50.0, 180.0, 400.0, 650.0, 900.0, 1_400.0, 2_100.0];
    let transfer = |mode, w: f64| {
        analytic_transfer(mode, w) * C64::new((w * delay_s).cos(), (w * delay_s).sin())
    };
    ModalTransferSamples {
        source_id: "analytic coupled receiver; signed first-order velocity transfers".into(),
        per_mode: (0..2)
            .map(|mode| omega.iter().map(|&w| transfer(mode, w)).collect())
            .collect(),
        omega_rad_s: omega,
        held_out: [110.0, 300.0, 700.0, 1_100.0, 1_800.0]
            .into_iter()
            .map(|w| (w, (0..2).map(|mode| transfer(mode, w)).collect()))
            .collect(),
        propagation_delay_s: delay_s,
        producer_relative_error_estimate: 2e-4,
    }
}
fn source(rate: u32, frequencies: &[f64]) -> CoupledModalTones {
    CoupledModalTones::new(
        fixture(rate, true, 1e-5),
        frequencies,
        budget(),
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}
fn receiver(rate: u32, frequencies: &[f64], delay_s: f64) -> BroadbandCoupledModalTones {
    BroadbandCoupledModalTones::build(
        source(rate, frequencies),
        &samples(delay_s),
        controls(),
        rate,
        budget(),
        &CancelGate::new_clock_free(),
    )
    .unwrap()
    .0
}
fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 2e-8 * a.abs().max(b.abs()) + 2e-13,
        "{a:e} != {b:e}"
    );
}
fn rotate(receiver: &BroadbandCoupledModalTones, cosine: &mut [f64], sine: &mut [f64]) {
    let tones = receiver.source().tone_count();
    for (index, (c, s)) in cosine.iter_mut().zip(sine).enumerate() {
        (*c, *s) = receiver.advance_force_phase(index % tones, *c, *s).unwrap();
    }
}

#[test]
fn startup_filter_uses_true_frequency_coupled_motion_and_preserves_source_accounts() {
    let gate = CancelGate::new_clock_free();
    let omega = 700.0;
    let force = [C64::new(0.4, -0.2), C64::new(-0.15, 0.3)];
    let mut network = fixture(4_000, true, 0.0);
    let response = network
        .harmonic_response(omega, &force, 1_000_000, &gate)
        .unwrap();
    network.models[0]
        .restore_states(
            &(0..2)
                .map(|mode| ModalAcousticState {
                    displacement_m_sqrt_kg: response.displacement_m_sqrt_kg[mode].re,
                    velocity_m_sqrt_kg_per_s: response.velocity_m_sqrt_kg_per_s[mode].re,
                })
                .collect::<Vec<_>>(),
        )
        .unwrap();
    for (x, z) in network
        .material_memory
        .iter_mut()
        .zip(&response.material_memory)
    {
        *x = z.re;
    }
    for (x, z) in network
        .radiation_state
        .iter_mut()
        .zip(&response.radiation_state)
    {
        *x = z.re;
    }
    let original = all_states(&network);
    let source = CoupledModalTones::new(network, &[omega], budget(), &gate).unwrap();
    let (mut receiver, report) = BroadbandCoupledModalTones::build(
        source,
        &samples(0.0),
        controls(),
        4_000,
        budget(),
        &gate,
    )
    .unwrap();
    assert_eq!(report.producer_relative_error_estimate, 2e-4);
    assert!(report.held_out_max_rel_error < 1e-8);
    assert_eq!(original, all_states(receiver.network()));
    assert!(
        receiver
            .network()
            .components()
            .iter()
            .flat_map(|model| model.modes())
            .all(|mode| mode.pressure_per_modal_velocity == C64::ZERO)
    );
    let mut cosine: Vec<_> = force.iter().map(|z| z.re).collect();
    let mut sine: Vec<_> = force.iter().map(|z| z.im).collect();
    for tick in 1..=80 {
        let actual = receiver.step(&cosine, &sine).unwrap();
        let time = f64::from(tick) / 4_000.0;
        let phase = C64::new((omega * time).cos(), -(omega * time).sin());
        let mut expected = 0.0;
        for mode in 0..2 {
            let velocity = response.velocity_m_sqrt_kg_per_s[mode];
            let steady_filter = velocity / C64::new(DECAY[mode], -omega);
            let state = (steady_filter * phase).re - steady_filter.re * (-DECAY[mode] * time).exp();
            expected += GAIN[mode] * DECAY[mode] * state + DIRECT[mode] * (velocity * phase).re;
            let fitted = receiver.transfer_at(mode, omega).unwrap();
            close(fitted.re, analytic_transfer(mode, omega).re);
            close(fitted.im, analytic_transfer(mode, omega).im);
        }
        close(actual, expected);
        let frame = receiver.network().last_frame().unwrap();
        assert_eq!(frame.observer_pressure_pa, 0.0);
        assert!(
            frame.component_dissipation_j > 0.0
                && frame.connection_dissipation_j > 0.0
                && frame.material_dissipation_j > 0.0
                && frame.radiation_dissipation_j > 0.0
        );
        if tick == 1 {
            let held_endpoint = (0..2)
                .map(|mode| {
                    let velocity = (response.velocity_m_sqrt_kg_per_s[mode] * phase).re;
                    (GAIN[mode] * (1.0 - (-DECAY[mode] * time).exp()) + DIRECT[mode]) * velocity
                })
                .sum::<f64>();
            assert!(
                (actual - held_endpoint).abs() > 1e-5 * actual.abs(),
                "fixture must distinguish the actual trajectory from held endpoint velocity"
            );
        }
        rotate(&receiver, &mut cosine, &mut sine);
    }
}

#[test]
fn interacting_tones_and_filter_startup_obey_continuous_time_partition() {
    let frequencies = [700.0, 701.3];
    let mut full = receiver(4_000, &frequencies, 0.0);
    let mut half = receiver(8_000, &frequencies, 0.0);
    let mut unobserved = source(4_000, &frequencies);
    let mut cosine = [0.4, -0.39, -0.2, 0.18];
    let mut sine = [-0.2, 0.21, 0.3, -0.29];
    for _ in 0..12 {
        let pressure = full.step(&cosine, &sine).unwrap();
        let frame = full.network().last_frame().unwrap().clone();
        let original = unobserved.step(&cosine, &sine).unwrap().clone();
        assert_eq!(all_states(full.network()), all_states(unobserved.network()));
        assert_eq!(
            frame.external_work_j.to_bits(),
            original.external_work_j.to_bits()
        );
        assert_eq!(
            frame.total_dissipation_j().to_bits(),
            original.total_dissipation_j().to_bits()
        );
        assert_eq!(
            frame.total_energy_j().to_bits(),
            original.total_energy_j().to_bits()
        );
        half.step(&cosine, &sine).unwrap();
        let first = half.network().last_frame().unwrap().clone();
        rotate(&half, &mut cosine, &mut sine);
        let subdivided = half.step(&cosine, &sine).unwrap();
        let second = half.network().last_frame().unwrap();
        close(pressure, subdivided);
        for (a, b) in all_states(full.network())
            .into_iter()
            .zip(all_states(half.network()))
        {
            close(a, b);
        }
        for (&a, &b) in full.state.iter().zip(&half.state) {
            close(a, b);
        }
        close(
            frame.external_work_j,
            first.external_work_j + second.external_work_j,
        );
        close(
            frame.total_dissipation_j(),
            first.total_dissipation_j() + second.total_dissipation_j(),
        );
        rotate(&half, &mut cosine, &mut sine);
    }
}

#[test]
fn integer_delay_preserves_exact_history_and_fractional_delay_keeps_its_gate() {
    let delay = 5.0 / 4_000.0;
    let mut delayed = receiver(4_000, &[700.0], delay);
    let mut raw = receiver(4_000, &[700.0], delay);
    // Keep the identical fitted operator and source. Only remove this reference
    // receiver's sampled propagation delay, using the same existing owner.
    raw.delay = Delay::new(0.0, 2).unwrap();
    let mut history = Vec::new();
    let mut cosine = [0.4, -0.2];
    let mut sine = [-0.2, 0.3];
    for tick in 0..24 {
        let actual = delayed.step(&cosine, &sine).unwrap();
        history.push(raw.step(&cosine, &sine).unwrap());
        let expected = if tick < 5 { 0.0 } else { history[tick - 5] };
        assert_eq!(actual.to_bits(), expected.to_bits());
        assert_eq!(all_states(delayed.network()), all_states(raw.network()));
        rotate(&delayed, &mut cosine, &mut sine);
    }
    let fractional = samples(4.375 / 4_000.0);
    let admitted = BroadbandModalControls {
        max_held_out_rel_error: 0.05,
        ..controls()
    };
    let (_, report) = BroadbandCoupledModalTones::build(
        source(4_000, &[700.0]),
        &fractional,
        admitted,
        4_000,
        budget(),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    assert!(report.held_out_max_rel_error > 1e-8 && report.held_out_max_rel_error < 0.05);
    assert!(matches!(
        BroadbandCoupledModalTones::build(
            source(4_000, &[700.0]),
            &fractional,
            controls(),
            4_000,
            budget(),
            &CancelGate::new_clock_free()
        ),
        Err(BroadbandCoupledModalError::Observation(
            BroadbandModalError::HeldOut { .. }
        ))
    ));
}

#[test]
fn late_delayed_pressure_refusal_preserves_source_filter_and_delay_together() {
    let make = || {
        let mut source = source(4_000, &[700.0]);
        source.network.config.maximum_total_energy_j = 1e12;
        source.network.config.maximum_abs_pressure_pa = 0.05;
        BroadbandCoupledModalTones::build(
            source,
            &samples(4.0 / 4_000.0),
            controls(),
            4_000,
            budget(),
            &CancelGate::new_clock_free(),
        )
        .unwrap()
        .0
    };
    let mut trial = make();
    let mut reference = make();
    for tick in 0..4 {
        let forces = if tick == 0 { [1e5, 0.0] } else { [0.0, 0.0] };
        assert_eq!(
            trial.step(&forces, &[0.0; 2]).unwrap().to_bits(),
            reference.step(&forces, &[0.0; 2]).unwrap().to_bits()
        );
    }
    let states = all_states(trial.network());
    let filter = trial.state.clone();
    let frame = trial.network().last_frame().unwrap().clone();
    let pressure = trial.last_pressure_pa();
    assert!(matches!(
        trial.step(&[0.0; 2], &[0.0; 2]),
        Err(BroadbandCoupledModalError::Modal(
            ModalCouplingError::Budget {
                what: "broadband receiver pressure",
                ..
            }
        ))
    ));
    assert_eq!(states, all_states(trial.network()));
    assert_eq!(filter, trial.state);
    assert_eq!(frame, *trial.network().last_frame().unwrap());
    assert_eq!(pressure, trial.last_pressure_pa());
    // After lifting only the observation ceiling, the rejected attempt must
    // leave the same future delayed signal as a receiver that never tried it.
    trial.source.network.config.maximum_abs_pressure_pa = 1e6;
    reference.source.network.config.maximum_abs_pressure_pa = 1e6;
    for _ in 0..8 {
        assert_eq!(
            trial.step(&[0.0; 2], &[0.0; 2]).unwrap().to_bits(),
            reference.step(&[0.0; 2], &[0.0; 2]).unwrap().to_bits()
        );
        assert_eq!(all_states(trial.network()), all_states(reference.network()));
        assert_eq!(trial.state, reference.state);
    }
}

#[test]
fn receiver_admits_complete_basis_clock_storage_and_held_out_operator() {
    let gate = CancelGate::new_clock_free();
    assert!(
        BroadbandCoupledModalTones::build(
            source(4_000, &[700.0]),
            &samples(0.0),
            controls(),
            8_000,
            budget(),
            &gate
        )
        .is_err()
    );
    let mut bad = samples(0.0);
    bad.held_out[2].1[0] = bad.held_out[2].1[0].scale(2.0);
    assert!(matches!(
        BroadbandCoupledModalTones::build(
            source(4_000, &[700.0]),
            &bad,
            controls(),
            4_000,
            budget(),
            &gate
        ),
        Err(BroadbandCoupledModalError::Observation(
            BroadbandModalError::HeldOut { .. }
        ))
    ));
    assert!(matches!(
        BroadbandCoupledModalTones::build(
            source(4_000, &[700.0]),
            &samples(10.0),
            controls(),
            4_000,
            CoupledToneBudget {
                max_dense_bytes: 64_000,
                ..budget()
            },
            &gate
        ),
        Err(BroadbandCoupledModalError::Observation(
            BroadbandModalError::Cap {
                what: "additional coupled receiver bytes",
                ..
            }
        ))
    ));
}
