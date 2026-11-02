use super::super::{
    ConstantLossConnection, ModalColumnConnection, ModalCouplingConfig, PreparedModalBudget,
};
use super::*;
use crate::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use crate::render::plate::impact::radiation::{Model, Pole};
use fs_math::c64::C64;
use fs_ports::PortedModule;

fn config() -> ModalCouplingConfig {
    ModalCouplingConfig {
        max_modes: 2,
        max_connections: 2,
        max_setup_terms: 1_000_000,
        nyquist_guard_fraction: 0.9,
        maximum_total_energy_j: 100.0,
        maximum_abs_pressure_pa: 1e6,
        maximum_abs_connection_force_n: 1e8,
        solve_relative_tolerance: 1e-11,
        energy_absolute_tolerance_j: 1e-12,
        energy_relative_tolerance: 1e-9,
    }
}
fn budget() -> CoupledToneBudget {
    CoupledToneBudget {
        max_dimension: 256,
        max_dense_bytes: 128 * 1_048_576,
        max_multiply_terms: 100_000_000_000,
    }
}
fn fixture(rate: u32, anchored: bool, rest: f64) -> CoupledModalSystem {
    let gate = CancelGate::new_clock_free();
    let mut model = ModalAcousticTimeModel::try_new(
        rate,
        [(500.0, 0.01, 2.0), (900.0, 0.04, -3.0)]
            .into_iter()
            .map(|(w, zeta, gain)| ModalAcousticMode {
                angular_frequency_rad_s: w,
                damping_ratio: zeta,
                pressure_per_modal_velocity: C64::new(gain, 0.0),
            })
            .collect(),
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    model
        .restore_states(&[
            ModalAcousticState {
                displacement_m_sqrt_kg: 1e-4,
                velocity_m_sqrt_kg_per_s: 0.02,
            },
            ModalAcousticState {
                displacement_m_sqrt_kg: -2e-4,
                velocity_m_sqrt_kg_per_s: -0.01,
            },
        ])
        .unwrap();
    let links = vec![
        ModalColumnConnection {
            column: vec![1.0, -0.3],
            stiffness_n_m: 4e4,
            damping_n_s_m: 20.0,
            rest_extension_m: rest,
        },
        ModalColumnConnection {
            column: vec![0.2, 0.8],
            stiffness_n_m: 1e4,
            damping_n_s_m: 10.0,
            rest_extension_m: -0.7 * rest,
        },
    ];
    let network = if anchored {
        CoupledModalSystem::new_with_columns_prepared(
            vec![model],
            links,
            config(),
            PreparedModalBudget {
                max_modes: 2,
                max_dense_bytes: 128 * 1_048_576,
                max_multiply_terms: 100_000_000_000,
            },
            &gate,
        )
        .unwrap()
        .with_constant_loss_connections(
            &[ConstantLossConnection {
                connection: 0,
                reference_omega_rad_s: 500.0,
            }],
            &gate,
        )
        .unwrap()
    } else {
        CoupledModalSystem::new_with_columns(vec![model], links, config(), &gate).unwrap()
    };
    network
        .with_material_relaxation(
            vec![fs_phs::RelaxationBranch {
                projection: vec![1.0, 0.0, -0.4, 0.0],
                stiffness: 3e4,
                relaxation_time_s: 0.002,
            }],
            Some(vec![0.007]),
            &gate,
        )
        .unwrap()
        .with_radiation_load(
            &Model {
                ports: 2,
                poles: vec![Pole {
                    omega: 700.0,
                    zeta: 0.1,
                    coupling: vec![300.0, -200.0],
                }],
            },
            &[0, 1],
            1,
            &gate,
        )
        .unwrap()
}
fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 3e-9 * a.abs().max(b.abs()) + 2e-14,
        "{a:e} != {b:e}"
    );
}
fn all_states(network: &CoupledModalSystem) -> Vec<f64> {
    let mut states: Vec<_> = network
        .components()
        .iter()
        .flat_map(|model| model.states())
        .flat_map(|x| [x.displacement_m_sqrt_kg, x.velocity_m_sqrt_kg_per_s])
        .collect();
    states.extend_from_slice(network.material_memory());
    states.extend_from_slice(network.radiation_state());
    states
}
fn rotate(runtime: &CoupledModalTones, cosine: &mut [f64], sine: &mut [f64]) {
    let tones = runtime.tone_count();
    for (index, (c, s)) in cosine.iter_mut().zip(sine).enumerate() {
        (*c, *s) = runtime.advance_force_phase(index % tones, *c, *s).unwrap();
    }
}

#[test]
fn true_frequency_motion_and_each_heat_match_complete_coupled_harmonic_solution() {
    // The steady response independently solves the original dynamic stiffness.
    // Its off-natural drive and mutual air reactions challenge component-wise
    // propagation; an anchored dashpot must retain its distinct loss direction.
    let gate = CancelGate::new_clock_free();
    let omega = core::f64::consts::TAU * 125.0;
    let force = [C64::new(0.4, -0.2), C64::new(-0.15, 0.3)];
    let mut network = fixture(4_000, true, 0.0);
    let response = network
        .harmonic_response(omega, &force, 1_000_000, &gate)
        .unwrap();
    network.models[0]
        .restore_states(
            &(0..2)
                .map(|i| ModalAcousticState {
                    displacement_m_sqrt_kg: response.displacement_m_sqrt_kg[i].re,
                    velocity_m_sqrt_kg_per_s: response.velocity_m_sqrt_kg_per_s[i].re,
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
    let initial_energy = network.total_energy_j().unwrap();
    let mut runtime = CoupledModalTones::new(network, &[omega], budget(), &gate).unwrap();
    let mut cosine: Vec<_> = force.iter().map(|z| z.re).collect();
    let mut sine: Vec<_> = force.iter().map(|z| z.im).collect();
    let mut totals = [0.0; 5];
    for tick in 1..=32 {
        let frame = runtime.step(&cosine, &sine).unwrap().clone();
        let phase = omega * f64::from(tick) / 4_000.0;
        for (index, state) in runtime.network().components()[0]
            .states()
            .iter()
            .enumerate()
        {
            let q = response.displacement_m_sqrt_kg[index];
            let v = response.velocity_m_sqrt_kg_per_s[index];
            close(
                state.displacement_m_sqrt_kg,
                q.re * phase.cos() + q.im * phase.sin(),
            );
            close(
                state.velocity_m_sqrt_kg_per_s,
                v.re * phase.cos() + v.im * phase.sin(),
            );
        }
        for (value, z) in runtime
            .network()
            .material_memory()
            .iter()
            .zip(&response.material_memory)
        {
            close(*value, z.re * phase.cos() + z.im * phase.sin());
        }
        for (value, z) in runtime
            .network()
            .radiation_state()
            .iter()
            .zip(&response.radiation_state)
        {
            close(*value, z.re * phase.cos() + z.im * phase.sin());
        }
        let states = runtime.network().components()[0].states();
        close(
            frame.observer_pressure_pa,
            2.0 * states[0].velocity_m_sqrt_kg_per_s - 3.0 * states[1].velocity_m_sqrt_kg_per_s,
        );
        for (sum, value) in totals.iter_mut().zip([
            frame.external_work_j,
            frame.component_dissipation_j,
            frame.connection_dissipation_j,
            frame.material_dissipation_j,
            frame.radiation_dissipation_j,
        ]) {
            *sum += value;
        }
        rotate(&runtime, &mut cosine, &mut sine);
    }
    let period = 1.0 / 125.0;
    for (integral, power) in totals.into_iter().zip([
        response.external_power_w,
        response.component_dissipation_w,
        response.connection_dissipation_w,
        response.material_dissipation_w,
        response.radiation_dissipation_w,
    ]) {
        assert!(integral > 0.0);
        close(integral, power * period);
    }
    close(runtime.network().total_energy_j().unwrap(), initial_energy);
}

#[test]
fn inharmonic_interference_rest_load_and_original_port_work_obey_time_partition() {
    let gate = CancelGate::new_clock_free();
    let frequencies = [700.0, 700.3, 0.0];
    let mut whole =
        CoupledModalTones::new(fixture(4_000, false, 1e-5), &frequencies, budget(), &gate).unwrap();
    let mut halves =
        CoupledModalTones::new(fixture(8_000, false, 1e-5), &frequencies, budget(), &gate).unwrap();
    let mut cosine = [0.4, 0.1, -0.2, -0.05, 0.3, 0.15];
    let mut sine = [-0.25, 0.2, 3.0, 0.1, -0.1, -99.0];
    let full = whole.step(&cosine, &sine).unwrap().clone();
    let first = halves.step(&cosine, &sine).unwrap().clone();
    let first_ports: Vec<_> = halves
        .network()
        .last_window_transfers()
        .unwrap()
        .iter()
        .map(|v| v.energy_into().value())
        .collect();
    rotate(&halves, &mut cosine, &mut sine);
    let second = halves.step(&cosine, &sine).unwrap().clone();
    for (a, b) in all_states(whole.network())
        .into_iter()
        .zip(all_states(halves.network()))
    {
        close(a, b);
    }
    for (a, b, c) in [
        (
            full.external_work_j,
            first.external_work_j,
            second.external_work_j,
        ),
        (
            full.component_dissipation_j,
            first.component_dissipation_j,
            second.component_dissipation_j,
        ),
        (
            full.connection_dissipation_j,
            first.connection_dissipation_j,
            second.connection_dissipation_j,
        ),
        (
            full.material_dissipation_j,
            first.material_dissipation_j,
            second.material_dissipation_j,
        ),
        (
            full.radiation_dissipation_j,
            first.radiation_dissipation_j,
            second.radiation_dissipation_j,
        ),
    ] {
        close(a, b + c);
    }
    for ((a, b), c) in whole
        .network()
        .last_window_transfers()
        .unwrap()
        .iter()
        .zip(halves.network().last_window_transfers().unwrap())
        .zip(first_ports)
    {
        close(a.energy_into().value(), b.energy_into().value() + c);
    }
    let fs_ports::ModuleWindowAudit::Available(audit) =
        whole.network().audit_last_window().unwrap()
    else {
        panic!("original coupled account was lost");
    };
    assert!(audit.violations().is_empty());
}

#[test]
fn exactly_coherent_large_sources_cancel_before_the_interval_force_bound() {
    let gate = CancelGate::new_clock_free();
    let mut grouped = CoupledModalTones::new(
        fixture(4_000, false, 1e-5),
        &[700.0, 700.0, 700.0],
        budget(),
        &gate,
    )
    .unwrap();
    let mut single =
        CoupledModalTones::new(fixture(4_000, false, 1e-5), &[700.0], budget(), &gate).unwrap();
    let a = grouped
        .step(
            &[1e16, 1.0, -1e16, 1e16, -0.5, -1e16],
            &[1e16, -0.25, -1e16, 1e16, 0.75, -1e16],
        )
        .unwrap()
        .clone();
    let b = single.step(&[1.0, -0.5], &[-0.25, 0.75]).unwrap().clone();
    assert_eq!(all_states(grouped.network()), all_states(single.network()));
    assert_eq!(a, b);
}

#[test]
fn late_pressure_refusal_retains_original_states_clock_frame_and_ports() {
    let gate = CancelGate::new_clock_free();
    let mut network = fixture(4_000, false, 0.0);
    network.models[0]
        .restore_states(
            &[ModalAcousticState {
                displacement_m_sqrt_kg: 0.0,
                velocity_m_sqrt_kg_per_s: 0.0,
            }; 2],
        )
        .unwrap();
    network.material_memory.fill(0.0);
    network.radiation_state.fill(0.0);
    network.config.maximum_abs_pressure_pa = 0.01;
    network.config.maximum_total_energy_j = 1e12;
    let mut runtime = CoupledModalTones::new(network, &[700.0], budget(), &gate).unwrap();
    runtime.step(&[1e-6, 0.0], &[0.0, 0.0]).unwrap();
    let states = all_states(runtime.network());
    let frame = runtime.network().last_frame().unwrap().clone();
    let ports = runtime.network().last_window_transfers().unwrap().to_vec();
    assert!(matches!(
        runtime.step(&[1e5, 0.0], &[0.0, 0.0]),
        Err(ModalCouplingError::Budget {
            what: "observer pressure",
            ..
        })
    ));
    assert_eq!(states, all_states(runtime.network()));
    assert_eq!(frame, *runtime.network().last_frame().unwrap());
    assert_eq!(ports, runtime.network().last_window_transfers().unwrap());
    assert!(runtime.step(&[f64::NAN, 0.0], &[0.0, 0.0]).is_err());
    assert_eq!(states, all_states(runtime.network()));
    assert_eq!(frame, *runtime.network().last_frame().unwrap());
}

#[test]
fn unresolved_tones_and_natural_mode_quadrature_refuse() {
    let gate = CancelGate::new_clock_free();
    assert!(
        CoupledModalTones::new(fixture(4_000, false, 0.0), &[20_000.0], budget(), &gate,).is_err()
    );
    assert!(
        CoupledModalTones::new(
            fixture(4_000, false, 0.0),
            &[700.0],
            CoupledToneBudget {
                max_dimension: 8,
                ..budget()
            },
            &gate,
        )
        .is_err()
    );
    let model = ModalAcousticTimeModel::try_new(
        4_000,
        vec![ModalAcousticMode {
            angular_frequency_rad_s: 500.0,
            damping_ratio: 0.01,
            pressure_per_modal_velocity: C64::new(0.0, 1.0),
        }],
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    let network = CoupledModalSystem::new(vec![model], vec![], config(), &gate).unwrap();
    assert!(CoupledModalTones::new(network, &[700.0], budget(), &gate).is_err());
}
