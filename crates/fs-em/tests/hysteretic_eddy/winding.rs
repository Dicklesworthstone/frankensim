//! G1/G3/G4: actual native circuits against independent scalar play balances.
use super::*;
use fs_em::winding::{
    end::EndWindingMode,
    hysteretic::{HystereticWindingCircuit, HystereticWindingControls, HystereticWindingError},
    nonlinear::{
        NonlinearWindingConfig, NonlinearWindingInitial, NonlinearWindingPort, NonlinearWindings,
        resistors::GraphResistorMaterial,
    },
    ohmic::WindingConductor,
};
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_solver::NewtonKrylovConfig;

const END: f64 = 0.02 + 0.003 * 4.0;
const GRAPH_R: f64 = 0.3;

fn terminal_controls() -> HystereticWindingControls {
    HystereticWindingControls {
        field: controls(),
        maximum_field_trials: 1000,
        maximum_newton_iterations: 30,
        difference_step_a: 1e-5,
        terminal_voltage_tolerance_v: 1e-8,
        power_tolerance_w: 1e-8,
        energy_tolerance_j: 2e-8,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-10,
            relative_tolerance: 1e-12,
            forcing_minimum: 1e-4,
            linear_restart: 8,
            ..NewtonKrylovConfig::default()
        },
    }
}

fn circuit<'a>(
    owner: &'a PlanarHystereticMagnetostatics,
    all: &[EddyConductor],
    core: bool,
    current_drive: bool,
    polarity: f64,
    scale: f64,
) -> HystereticWindingCircuit<'a> {
    let core_cells: Vec<_> = all
        .iter()
        .filter(|s| core && matches!(s.triangle, 1 | 7))
        .cloned()
        .collect();
    let copper: Vec<_> = all
        .iter()
        .filter(|s| !core || !matches!(s.triangle, 1 | 7))
        .map(|s| WindingConductor {
            triangle: s.triangle,
            conductor: s.conductor.clone(),
        })
        .collect();
    let insulators: Vec<_> = copper.iter().map(|s| s.triangle).collect();
    let mut profile = vec![0.0; 8];
    for s in &copper {
        profile[s.triangle] =
            polarity * if core { 8.0 } else { 6.0 } * if s.triangle < 4 { 1.0 } else { -1.0 };
    }
    let initial = owner.demagnetized(0.0, &mut || false).unwrap();
    let field = PlanarHystereticEddy::new(
        &initial,
        core_cells,
        &insulators,
        &boundary(6, 0.0),
        0.0,
        controls(),
        &mut || false,
    )
    .unwrap();
    let graph = if current_drive {
        CircuitGraph {
            node_count: 2,
            branches: vec![
                (1, 0, Branch::CurrentSource { port: 0 }),
                (1, 0, Branch::Inductor { henries: 999.0 }),
            ],
            transformers: vec![],
        }
    } else {
        CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 0, Branch::VoltageSource { port: 0 }),
                (1, 2, Branch::Inductor { henries: 999.0 }),
                (2, 0, Branch::Resistor { ohms: GRAPH_R }),
            ],
            transformers: vec![],
        }
    };
    let port = if current_drive {
        NonlinearWindingPort::VoltageInput {
            initial_voltage_v: 0.0,
            residual_scale_ohm: scale,
            current_tolerance_a: 1e-8,
        }
    } else {
        NonlinearWindingPort::CurrentInput
    };
    HystereticWindingCircuit::new_with_end_modes(
        field,
        graph,
        NonlinearWindings {
            branch_indices: vec![1],
            profiles_per_m2: vec![profile],
            end_inductance_h: vec![0.02],
        },
        vec![port],
        vec![EndWindingMode {
            inductance_h: 0.003,
            turns: vec![2.0],
        }],
        copper,
        NonlinearWindingInitial {
            currents_a: vec![0.0],
            remaining_fluxes_wb: vec![],
            charges_c: vec![],
            external_inputs: vec![0.0],
        },
        NonlinearWindingConfig {
            max_windings: 4,
            max_profile_values: 32,
            current_balance_tolerance_turns: 1e-10,
        },
        &mut || false,
    )
    .unwrap()
}

fn root(residual: impl Fn(f64) -> f64) -> f64 {
    let (mut low, mut high) = (-29.0, 29.0);
    assert!(residual(low) < 0.0 && residual(high) > 0.0);
    for _ in 0..90 {
        let middle = f64::midpoint(low, high);
        if residual(middle) < 0.0 {
            low = middle;
        } else {
            high = middle;
        }
    }
    f64::midpoint(low, high)
}

#[test]
fn voltage_driven_minor_loops_retain_vector_history_and_pay_each_physical_loss_once() {
    let (depth, sigma, dt) = (0.3, 18.0, 0.2);
    let (owner, all) = problem(paired_squares(depth), MU, true, sigma);
    for polarity in [1.0, -1.0] {
        let mut state = circuit(&owner, &all, false, false, polarity, 1.0);
        let rw = 72.0 * depth / sigma;
        near(state.resistance_ohm()[0], rw);
        near(state.end_inductance().matrix_h()[0], END);
        let (mut old_i, mut old_z, mut old_b) = (0.0, 0.0, 0.0);
        let mut pinning = 0.0;
        for voltage in [20.0, 0.0, 8.0, 0.0, -20.0, 0.0] {
            // The two signed winding sections give H=polarity*I and
            // lambda=2*depth*polarity*B. Scalar play is an independent oracle.
            let i = root(|i| {
                let h = polarity * i;
                let b = MU * h + play(old_z, h);
                (rw + GRAPH_R) * i + (2.0 * depth * polarity * (b - old_b) + END * (i - old_i)) / dt
                    - voltage
            });
            let h = polarity * i;
            let z = play(old_z, h);
            let b = MU * h + z;
            let old_storage =
                2.0 * depth * storage(polarity * old_i, old_z) + 0.5 * END * old_i * old_i;
            let new_storage = 2.0 * depth * storage(h, z) + 0.5 * END * i * i;
            let q_pinning = 2.0 * depth * PINNING * (z - old_z).abs();
            let q_copper = dt * rw * i * i;
            let q_graph = dt * GRAPH_R * i * i;
            let supply = dt * voltage * i;
            let next = state
                .step_to_candidate(
                    state.time_s() + dt,
                    vec![],
                    state.conductors().to_vec(),
                    &[voltage],
                    terminal_controls(),
                    &mut || false,
                )
                .unwrap();
            let step = next.step();
            near(next.circuit().currents_a()[0], i);
            near(next.circuit().field().potential_wb_m()[6], 0.5 * b);
            near(next.circuit().field().potential_wb_m()[7], -0.5 * b);
            near(
                next.circuit().flux_linkages_wb()[0],
                2.0 * depth * polarity * b + END * i,
            );
            near(step.stored_before_j, old_storage);
            near(step.stored_after_j, new_storage);
            near(step.stored_change_j, new_storage - old_storage);
            near(step.field.hysteresis_heat_j, q_pinning);
            near(step.field.joule_heat_j, 0.0);
            near(step.cells.iter().map(|c| c.heat_j).sum(), q_copper);
            near(step.graph_resistors[0].heat_j, q_graph);
            near(step.heat_j, q_pinning + q_copper + q_graph);
            near(step.supplied_work_j, supply);
            near(
                step.numerical_dissipation_j,
                supply - (new_storage - old_storage) - step.heat_j,
            );
            near(step.terminal_voltage_v[0], voltage - GRAPH_R * i);
            near(step.energy_residual_j, 0.0);
            for (cell, memory) in next
                .circuit()
                .field()
                .magnetic_state()
                .cell_states()
                .iter()
                .enumerate()
            {
                let orientation = if cell < 4 { 1.0 } else { -1.0 };
                for axis in 0..2 {
                    near(
                        memory.polarizations_t()[0][axis],
                        orientation * z * DIRECTIONS[cell % 4][axis],
                    );
                }
            }
            pinning += step.field.hysteresis_heat_j;
            state = next.into_parts().0;
            (old_i, old_z, old_b) = (i, z, b);
        }
        assert!(pinning > 1.0 && old_z.abs() > 0.4);
    }
}

#[test]
fn current_drive_resolves_core_eddy_diffusion_pinning_and_end_induction_without_added_impedance() {
    let (depth, sigma, dt) = (0.3, 18.0, 0.2);
    let (owner, all) = problem(paired_squares(depth), MU, true, sigma);
    for scale in [0.2, 5.0] {
        let mut state = circuit(&owner, &all, true, true, 1.0, scale);
        let rw = 96.0 * depth / sigma;
        near(state.resistance_ohm()[0], rw);
        let (mut old_i, mut old_h, mut old_z, mut old_b) = (0.0, 0.0, 0.0, 0.0);
        for i in [12.0, 0.0, 6.0, 1.0, 5.5, 0.0, -12.0, 0.0] {
            // Two edge-connected cells carry core conductivity; the other six
            // are explicit stranded copper. Native P1 integration gives this
            // independent two-centre diffusion balance and its Ohmic integral.
            let h = root(|h| h + sigma / (96.0 * dt) * (MU * h + play(old_z, h) - old_b) - i);
            let z = play(old_z, h);
            let b = MU * h + z;
            let pinning = 2.0 * depth * PINNING * (z - old_z).abs();
            let core_heat = sigma * depth / (48.0 * dt) * (b - old_b) * (b - old_b);
            let copper_heat = dt * rw * i * i;
            let voltage = (2.0 * depth * (b - old_b) + END * (i - old_i)) / dt + rw * i;
            let before = 2.0 * depth * storage(old_h, old_z) + 0.5 * END * old_i * old_i;
            let after = 2.0 * depth * storage(h, z) + 0.5 * END * i * i;
            let next = state
                .step_to_candidate(
                    state.time_s() + dt,
                    state.field().conductors().to_vec(),
                    state.conductors().to_vec(),
                    &[i],
                    terminal_controls(),
                    &mut || false,
                )
                .unwrap();
            let step = next.step();
            near(next.circuit().currents_a()[0], i);
            near(next.circuit().field().potential_wb_m()[6], 0.5 * b);
            near(step.terminal_voltage_v[0], voltage);
            near(step.field.joule_heat_j, core_heat);
            near(step.field.hysteresis_heat_j, pinning);
            near(step.cells.iter().map(|cell| cell.heat_j).sum(), copper_heat);
            near(step.heat_j, pinning + core_heat + copper_heat);
            near(step.stored_change_j, after - before);
            near(step.supplied_work_j, dt * i * voltage);
            near(
                step.numerical_dissipation_j,
                dt * i * voltage - (after - before) - step.heat_j,
            );
            near(step.energy_residual_j, 0.0);
            assert!(
                step.field
                    .cells
                    .iter()
                    .all(|cell| matches!(cell.source.triangle, 1 | 7))
            );
            assert!(
                step.cells
                    .iter()
                    .all(|cell| !matches!(cell.source.triangle, 1 | 7))
            );
            state = next.into_parts().0;
            (old_i, old_h, old_z, old_b) = (i, h, z, b);
        }
    }
}

#[test]
fn branchless_field_matches_backward_euler_rl_without_a_fitted_inductance() {
    let (depth, sigma, dt) = (0.3, 18.0, 0.2);
    let (owner, all) = problem(paired_squares(depth), MU, false, sigma);
    let mut state = circuit(&owner, &all, false, false, 1.0, 1.0);
    let l = 2.0 * depth * MU + END;
    let r = 72.0 * depth / sigma + GRAPH_R;
    let mut previous = 0.0;
    for voltage in [3.0, 0.0, -2.0] {
        let expected = (l * previous + dt * voltage) / (l + dt * r);
        let candidate = state
            .step_to_candidate(
                state.time_s() + dt,
                vec![],
                state.conductors().to_vec(),
                &[voltage],
                terminal_controls(),
                &mut || false,
            )
            .unwrap();
        near(candidate.circuit().currents_a()[0], expected);
        near(candidate.step().heat_j, dt * r * expected * expected);
        near(
            candidate.step().numerical_dissipation_j,
            0.5 * l * (expected - previous) * (expected - previous),
        );
        assert_eq!(candidate.step().field.hysteresis_heat_j, 0.0);
        state = candidate.into_parts().0;
        previous = expected;
    }
}

#[test]
fn reassembly_failed_trials_and_late_cancellation_preserve_full_field_and_circuit_history() {
    let (owner, all) = problem(paired_squares(0.3), MU, true, 18.0);
    let zero = circuit(&owner, &all, false, false, 1.0, 1.0);
    let state = zero
        .step_to_candidate(
            0.2,
            vec![],
            zero.conductors().to_vec(),
            &[20.0],
            terminal_controls(),
            &mut || false,
        )
        .unwrap()
        .into_parts()
        .0;
    let before = state.field().magnetic_state().clone();
    let before_native = state.state().to_vec();
    let replacement = all[0]
        .conductor
        .clone()
        .into_uniform_conductor(9.0, 1.0)
        .unwrap();
    let rebound = state
        .with_graph_resistors(
            vec![GraphResistorMaterial {
                branch: 2,
                conductor: replacement,
            }],
            &mut || false,
        )
        .unwrap();
    assert_eq!(rebound.state(), state.state());
    assert_eq!(rebound.currents_a(), state.currents_a());
    assert_eq!(
        rebound.field().magnetic_state().cell_states(),
        before.cell_states()
    );
    near(rebound.time_s(), state.time_s());
    let mut limited = terminal_controls();
    limited.maximum_field_trials = 1;
    assert!(matches!(
        state.step_to_candidate(
            0.4,
            vec![],
            state.conductors().to_vec(),
            &[-20.0],
            limited,
            &mut || false
        ),
        Err(HystereticWindingError::TrialBudget)
    ));
    let (_, foreign) = sources(MU, true, 19.0);
    let foreign = state
        .conductors()
        .iter()
        .map(|s| WindingConductor {
            triangle: s.triangle,
            conductor: foreign.clone(),
        })
        .collect();
    assert!(
        state
            .step_to_candidate(
                0.4,
                vec![],
                foreign,
                &[-20.0],
                terminal_controls(),
                &mut || false
            )
            .is_err()
    );
    let mut checkpoints = 0;
    let expected = state
        .step_to_candidate(
            0.4,
            vec![],
            state.conductors().to_vec(),
            &[-20.0],
            terminal_controls(),
            &mut || {
                checkpoints += 1;
                false
            },
        )
        .unwrap();
    let mut visited = 0;
    assert!(matches!(
        state.step_to_candidate(
            0.4,
            vec![],
            state.conductors().to_vec(),
            &[-20.0],
            terminal_controls(),
            &mut || {
                visited += 1;
                visited >= checkpoints / 2
            }
        ),
        Err(HystereticWindingError::Cancelled)
    ));
    assert_eq!(state.state(), before_native);
    assert_eq!(state.field().potential_wb_m(), before.potential_wb_m());
    assert_eq!(
        state.field().magnetic_state().cell_states(),
        before.cell_states()
    );
    near(state.time_s(), 0.2);
    let retry = state
        .step_to_candidate(
            0.4,
            vec![],
            state.conductors().to_vec(),
            &[-20.0],
            terminal_controls(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(
        retry.circuit().currents_a(),
        expected.circuit().currents_a()
    );
    assert_eq!(
        retry.circuit().field().magnetic_state().cell_states(),
        expected.circuit().field().magnetic_state().cell_states()
    );
}
