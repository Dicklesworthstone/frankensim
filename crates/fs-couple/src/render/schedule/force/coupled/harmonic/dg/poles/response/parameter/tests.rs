use super::super::super::super::frequency::tests::{loaded, oscillator};
use super::super::tests::{options, physical_pole, pole_options, preparation, vector_error};
use super::*;

const COLUMN: [f64; 1] = [0.003];
const STEP: f64 = 0.0002;

// Reconstruct the original physical owners at every material/connection sample.
// No tangent matrix or finite-difference mode is fed to the production adapter.
fn varied_oscillator(k: f64, c: f64, strength: f64, tau: f64) -> CoupledModalSystem {
    use crate::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeBudget};
    use crate::render::plate::impact::radiation::{Model, Pole};
    let gate = CancelGate::new_clock_free();
    let model = ModalAcousticTimeModel::try_new(
        48_000,
        vec![ModalAcousticMode {
            angular_frequency_rad_s: 1800.0,
            damping_ratio: 0.003,
            pressure_per_modal_velocity: C64::ZERO,
        }],
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    let config = oscillator().config;
    CoupledModalSystem::new_with_columns(
        vec![model],
        vec![ModalColumnConnection {
            column: vec![0.6],
            stiffness_n_m: k,
            damping_n_s_m: c,
            rest_extension_m: 0.0,
        }],
        config,
        &gate,
    )
    .unwrap()
    .with_material_relaxation(
        vec![fs_phs::RelaxationBranch {
            projection: vec![1.0, 0.0],
            stiffness: strength,
            relaxation_time_s: tau,
        }],
        Some(vec![0.004]),
        &gate,
    )
    .unwrap()
    .with_radiation_load(
        &Model {
            ports: 1,
            poles: vec![Pole {
                omega: 2200.0,
                zeta: 0.02,
                coupling: vec![100.0],
            }],
        },
        &[0],
        1,
        &gate,
    )
    .unwrap()
}
fn changed(direction: DgPoleParameterDirection<'_>, theta: f64) -> (CoupledModalSystem, Vec<f64>) {
    let (mut k, mut c, mut strength, mut tau) = (20_000.0, 4.0, 30_000.0, 0.002);
    let mut column = COLUMN.to_vec();
    match direction {
        DgPoleParameterDirection::Connection {
            stiffness_n_m,
            damping_n_s_m,
            ..
        } => {
            k += theta * stiffness_n_m;
            c += theta * damping_n_s_m;
        }
        DgPoleParameterDirection::Maxwell {
            stiffness,
            relaxation_time_s,
            ..
        } => {
            strength += theta * stiffness;
            tau += theta * relaxation_time_s;
        }
        DgPoleParameterDirection::InletFlowColumn(values) => {
            for (c, v) in column.iter_mut().zip(values) {
                *c += theta * v;
            }
        }
    }
    (varied_oscillator(k, c, strength, tau), column)
}
fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(actual.abs()) + 1e-25,
        "{actual:?} != {expected:?}"
    );
}
fn difference(a: &[C64], b: &[C64]) -> Vec<C64> {
    a.iter()
        .zip(b)
        .map(|(&x, &y)| (x - y).scale(0.5 / STEP))
        .collect()
}

#[test]
fn complete_loaded_parameter_shapes_match_independently_rebuilt_physical_sources() {
    let gate = CancelGate::new_clock_free();
    let directions = [
        DgPoleParameterDirection::Connection {
            connection: 0,
            stiffness_n_m: 20_000.0,
            damping_n_s_m: 0.0,
        },
        DgPoleParameterDirection::Connection {
            connection: 0,
            stiffness_n_m: 0.0,
            damping_n_s_m: 4.0,
        },
        DgPoleParameterDirection::Maxwell {
            arm: 0,
            stiffness: 30_000.0,
            relaxation_time_s: 0.0,
        },
        DgPoleParameterDirection::Maxwell {
            arm: 0,
            stiffness: 0.0,
            relaxation_time_s: 0.002,
        },
        DgPoleParameterDirection::InletFlowColumn(&COLUMN),
    ];
    let system = oscillator();
    for circular in [false, true] {
        let load = loaded(circular, circular);
        let saved = load.state().clone();
        let source = system
            .prepare_resonance_poles_with_load(&load, &COLUMN, preparation(1800.0), &gate)
            .unwrap();
        let pole = physical_pole(&source, &gate);
        let prepared = source
            .prepare_pole_response(&pole, options(), &gate)
            .unwrap();
        let layout = source.coordinates();
        let mut observation = vec![C64::ZERO; source.dimension()];
        observation[layout.material_memory.start] = C64::new(0.3, 0.7);
        let derivatives = prepared
            .parameter_derivatives(&directions, &[&observation], &gate)
            .unwrap();
        for (&direction, derivative) in directions.iter().zip(derivatives) {
            let mut samples = Vec::new();
            for theta in [STEP, -STEP] {
                let (changed, column) = changed(direction, theta);
                let source = changed
                    .prepare_resonance_poles_with_load(&load, &column, preparation(1800.0), &gate)
                    .unwrap();
                let mut pole = source
                    .resonance_pole(pole.laplace_pole_per_s, &pole.state, &pole_options(), &gate)
                    .unwrap();
                let pivot = pole.state[prepared.normalization_pivot()];
                for value in &mut pole.state {
                    *value = *value / pivot;
                }
                for value in &mut pole.modal_velocity_m_sqrt_kg_per_s {
                    *value = *value / pivot;
                }
                for value in &mut pole.radiation_velocity {
                    *value = *value / pivot;
                }
                pole.inlet_flow_m3_s = pole.inlet_flow_m3_s / pivot;
                pole.inlet_pressure_pa = pole.inlet_pressure_pa / pivot;
                samples.push(pole);
            }
            let (plus, minus) = (&samples[0], &samples[1]);
            close(
                derivative.laplace_pole_per_s,
                (plus.laplace_pole_per_s - minus.laplace_pole_per_s).scale(0.5 / STEP),
                2e-5,
            );
            let expected = difference(&plus.state, &minus.state);
            let mut groups = vec![
                layout.modal_displacement.clone(),
                layout.radiation_displacement.clone(),
                layout.duct_energy.clone(),
                layout.terminal_state.clone(),
                layout.material_memory.clone(),
            ];
            groups.extend(layout.branch_energy.iter().map(|(_, span)| span.clone()));
            groups.extend(
                layout
                    .channels
                    .iter()
                    .filter(|c| !c.velocity_cells.is_empty())
                    .map(|c| c.velocity_cells.clone()),
            );
            for group in groups {
                // The fixed pivot has exactly zero derivative; finite difference
                // of that constrained coordinate is checked separately below.
                if group.clone().all(|i| i == prepared.normalization_pivot()) {
                    continue;
                }
                assert!(
                    vector_error(&derivative.state[group.clone()], &expected[group.clone()]) < 4e-5,
                    "{direction:?} retained state group {group:?}"
                );
            }
            assert_eq!(derivative.state[prepared.normalization_pivot()], C64::ZERO);
            assert!(expected[prepared.normalization_pivot()].abs() < 1e-10);
            for &(port, index) in &layout.terminal_flow {
                if load
                    .terminal_channels()
                    .iter()
                    .any(|c| c.input().port == port && c.channel().is_none())
                {
                    assert!(derivative.state[index].abs() < 1e-15);
                } else {
                    close(derivative.state[index], expected[index], 4e-5);
                }
            }
            for &(_, index) in &layout.junction_flow {
                close(derivative.state[index], expected[index], 4e-5);
            }
            for index in layout.channels.iter().filter_map(|c| c.pressure_drop) {
                close(derivative.state[index], expected[index], 4e-5);
            }
            assert!(
                vector_error(
                    &derivative.modal_velocity_m_sqrt_kg_per_s,
                    &difference(
                        &plus.modal_velocity_m_sqrt_kg_per_s,
                        &minus.modal_velocity_m_sqrt_kg_per_s
                    )
                ) < 4e-5
            );
            assert!(
                vector_error(
                    &derivative.radiation_velocity,
                    &difference(&plus.radiation_velocity, &minus.radiation_velocity)
                ) < 4e-5
            );
            close(
                derivative.inlet_flow_m3_s,
                (plus.inlet_flow_m3_s - minus.inlet_flow_m3_s).scale(0.5 / STEP),
                4e-5,
            );
            close(
                derivative.inlet_pressure_pa,
                (plus.inlet_pressure_pa - minus.inlet_pressure_pa).scale(0.5 / STEP),
                4e-5,
            );
            close(
                derivative.observations[0],
                C64::new(0.3, 0.7) * derivative.state[layout.material_memory.start],
                1e-13,
            );
            assert!(
                derivative.linear_residual < 1e-11 && derivative.original_relative_residual < 1e-11
            );
        }
        assert_eq!(*load.state(), saved);
    }
}

#[test]
fn native_parameter_units_fixed_row_scaling_and_replay_preserve_the_complete_tangent() {
    let gate = CancelGate::new_clock_free();
    let system = oscillator();
    let load = loaded(false, true);
    let mut reference: Option<(C64, Vec<C64>)> = None;
    for rate in [1800.0, 36_000.0] {
        let source = system
            .prepare_resonance_poles_with_load(&load, &COLUMN, preparation(rate), &gate)
            .unwrap();
        let pole = physical_pole(&source, &gate);
        let prepared = source
            .prepare_pole_response(&pole, options(), &gate)
            .unwrap();
        for amplitude in [1.0, 1e-100, 1e100] {
            let direction = DgPoleParameterDirection::Maxwell {
                arm: 0,
                stiffness: 30_000.0 * amplitude,
                relaxation_time_s: 0.0004 * amplitude,
            };
            let derivative = prepared
                .parameter_derivatives(&[direction], &[], &gate)
                .unwrap()
                .remove(0);
            let replay = prepared
                .parameter_derivatives(&[direction], &[], &gate)
                .unwrap()
                .remove(0);
            assert_eq!(derivative.state, replay.state);
            assert_eq!(derivative.laplace_pole_per_s, replay.laplace_pole_per_s);
            let normalized: Vec<_> = derivative
                .state
                .iter()
                .map(|&v| v.scale(1.0 / amplitude))
                .collect();
            // The original-row audit must detect a corrupted published shape,
            // independently of the numerical border's successful solve receipt.
            let (action, envelope) = prepared.parameter_action(direction, &gate).unwrap();
            let mut corrupt = derivative.state.clone();
            corrupt[0] = corrupt[0] + C64::from_re(0.001 * amplitude);
            assert!(
                source
                    .pencil
                    .parameter_residual(
                        prepared.mode.laplace_pole_per_s,
                        &prepared.mode.state,
                        derivative.laplace_pole_per_s,
                        &corrupt,
                        &action,
                        &envelope,
                        &mut || false,
                    )
                    .unwrap()
                    > 1e-6
            );
            let value = derivative.laplace_pole_per_s.scale(1.0 / amplitude);
            if let Some((prior, vector)) = &reference {
                close(value, *prior, 3e-8);
                assert!(vector_error(&normalized, vector) < 3e-8);
            } else {
                reference = Some((value, normalized));
            }
        }
        let zero = prepared
            .parameter_derivatives(
                &[DgPoleParameterDirection::InletFlowColumn(&[0.0])],
                &[],
                &gate,
            )
            .unwrap()
            .remove(0);
        assert!(zero.state.iter().all(|&v| v == C64::ZERO));
        assert_eq!(zero.laplace_pole_per_s, C64::ZERO);
        assert_eq!(zero.original_relative_residual, 0.0);
    }
}

#[test]
fn source_dependency_refusals_and_budgets_leave_the_prepared_factor_reusable() {
    let gate = CancelGate::new_clock_free();
    let mut system = oscillator();
    // This test isolates the derivative contract on an already prepared source
    // column; it does not test or replace the constant-loss matrix-function owner.
    system.connections[0].constant_loss_reference_omega_rad_s = Some(1800.0);
    system.harmonic_damping_columns[0] = Some(vec![0.7]);
    let load = loaded(false, false);
    let source = system
        .prepare_resonance_poles_with_load(&load, &COLUMN, preparation(1800.0), &gate)
        .unwrap();
    let pole = physical_pole(&source, &gate);
    let prepared = source
        .prepare_pole_response(&pole, options(), &gate)
        .unwrap();
    let damping = DgPoleParameterDirection::Connection {
        connection: 0,
        stiffness_n_m: 0.0,
        damping_n_s_m: 4.0,
    };
    let baseline = prepared
        .parameter_derivatives(&[damping], &[], &gate)
        .unwrap()
        .remove(0);
    for direction in [
        DgPoleParameterDirection::Connection {
            connection: 0,
            stiffness_n_m: 1.0,
            damping_n_s_m: 0.0,
        },
        DgPoleParameterDirection::Connection {
            connection: 1,
            stiffness_n_m: 0.0,
            damping_n_s_m: 1.0,
        },
        DgPoleParameterDirection::Maxwell {
            arm: 1,
            stiffness: 0.0,
            relaxation_time_s: 1.0,
        },
        DgPoleParameterDirection::Maxwell {
            arm: 0,
            stiffness: f64::NAN,
            relaxation_time_s: 0.0,
        },
        DgPoleParameterDirection::InletFlowColumn(&[]),
        DgPoleParameterDirection::InletFlowColumn(&[f64::INFINITY]),
    ] {
        assert!(
            prepared
                .parameter_derivatives(&[damping, direction], &[], &gate)
                .is_err()
        );
    }
    assert!(prepared.parameter_derivatives(&[], &[], &gate).is_err());
    assert!(
        prepared
            .parameter_derivatives(&[damping; 9], &[], &gate)
            .is_err()
    );
    assert!(
        prepared
            .parameter_derivatives(&[damping], &[&[]], &gate)
            .is_err()
    );
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        prepared.parameter_derivatives(&[damping], &[], &cancelled),
        Err(DgPoleResponseError::Cancelled)
    ));
    assert_eq!(
        prepared
            .parameter_derivatives(&[damping], &[], &gate)
            .unwrap()[0]
            .state,
        baseline.state
    );
    // Refuse a singular zero-strength chart rather than freezing its sqrt chain.
    drop(prepared);
    drop(source);
    system.relaxation_branches[0].stiffness = 0.0;
    let source = system
        .prepare_resonance_poles_with_load(&load, &COLUMN, preparation(1800.0), &gate)
        .unwrap();
    let pole = physical_pole(&source, &gate);
    let prepared = source
        .prepare_pole_response(&pole, options(), &gate)
        .unwrap();
    assert!(
        prepared
            .parameter_derivatives(
                &[DgPoleParameterDirection::Maxwell {
                    arm: 0,
                    stiffness: 1.0,
                    relaxation_time_s: 0.0,
                }],
                &[],
                &gate
            )
            .is_err()
    );
}
