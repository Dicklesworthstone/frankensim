//! Physical electrode work drives retained 3D vector memory and Ohmic current.
#[path = "hysteretic/circuit.rs"]
mod circuit;
#[allow(dead_code)]
#[path = "../hysteretic/common.rs"]
mod data;
#[allow(dead_code)]
#[path = "hysteretic_common.rs"]
mod reference;

use super::super::{grid, near, preparation};
use data::{DT, PIN, conductor, field_controls, resistivity, stored};
use fs_em::spatial::{
    SpatialMagneticError,
    vector_potential::induction::hysteretic::{
        SpatialHystereticInductionError, terminals::SpatialHystereticTerminalInduction,
    },
};
use fs_material::field_storage::hysteresis::vector::VectorHysteresisError;
use reference::{Endpoint, assert_memory, at_field, electrodes, joule};

#[test]
fn physical_electrodes_follow_scalar_play_work_and_preserve_memory_across_control_changes() {
    let card = data::card(500.0);
    let (positions, tets) = grid(1);
    let problem = data::problem(&card, positions, tets);
    let initial = data::initial(&problem, Some(conductor(&card, 300.0)));
    let original = SpatialHystereticTerminalInduction::new(
        initial,
        &electrodes(&problem),
        preparation(),
        &mut || false,
    )
    .unwrap();
    let mut state = original.clone();
    let mut old = Endpoint::default();
    let sigma = resistivity(300.0).recip();
    let mut returned = None;
    for (index, h) in [12.0, 0.0, 6.0, 1.0, 6.0].into_iter().enumerate() {
        let expected = at_field(old, sigma, DT, h);
        let before = (
            state.potential_wb().to_vec(),
            state.cell_states().to_vec(),
            state.time_s(),
        );
        let mut polls = 0;
        let candidate = state
            .step_voltage_candidate(
                &[expected.voltage, 0.0],
                &[[0.0; 3]; 6],
                DT,
                field_controls(),
                &mut || {
                    polls += 1;
                    false
                },
            )
            .unwrap();
        let current = state
            .clone()
            .into_current_driven(preparation(), &mut || false)
            .unwrap();
        assert_eq!(current.cell_states(), before.1);
        let dual = current
            .step_current_candidate(
                &[expected.incoming, -expected.incoming],
                &[[0.0; 3]; 6],
                DT,
                field_controls(),
                &mut || false,
            )
            .unwrap();
        assert_memory(candidate.state().cell_states(), expected);
        assert_memory(dual.state().cell_states(), expected);
        let heat = joule(old, expected, sigma, DT);
        let pinning = PIN * (expected.z - old.z).abs();
        let change = stored(expected.h, expected.z) - stored(old.h, old.z);
        let numerical = expected.h * (expected.b - old.b) - change - pinning;
        let step = candidate.step();
        near(
            step.terminals[0].incoming_current_a,
            expected.incoming,
            4e-7,
        );
        near(
            step.terminals[1].incoming_current_a,
            -expected.incoming,
            4e-7,
        );
        near(
            dual.step().terminals[1].voltage_to_reference_v,
            -expected.voltage,
            3e-8,
        );
        near(
            step.terminal_work_j,
            DT * expected.voltage * expected.incoming,
            4e-7,
        );
        near(step.joule_heat_j, heat, 4e-7);
        near(step.hysteresis_heat_j, pinning, 3e-8);
        near(step.heat_j, heat + pinning, 4e-7);
        near(dual.step().heat_j, heat + pinning, 4e-7);
        near(step.stored_change_j, change, 3e-8);
        near(step.numerical_dissipation_j, numerical, 3e-8);
        near(step.impressed_work_j, 0.0, 1e-14);
        near(
            step.energy_residual_j,
            0.0,
            field_controls().induction.induction.energy_tolerance_j,
        );
        assert_eq!(step.cells.len(), 6);
        for cell in &step.magnetic_cells {
            near(cell.irreversible_heat_j, pinning / 6.0, 3e-8);
        }
        if index == 0 {
            assert!(heat > 1.0 && pinning > 1.0 && numerical > 1.0);
        } else if index == 2 {
            returned = Some(candidate.state().cell_states().to_vec());
            let mut seen = 0;
            assert!(matches!(
                state.step_voltage_candidate(
                    &[expected.voltage, 0.0],
                    &[[0.0; 3]; 6],
                    DT,
                    field_controls(),
                    &mut || {
                        seen += 1;
                        seen >= polls - 1
                    },
                ),
                Err(
                    SpatialHystereticInductionError::Spatial(SpatialMagneticError::Cancelled)
                        | SpatialHystereticInductionError::Material(
                            VectorHysteresisError::Cancelled
                        )
                )
            ));
        } else if index == 4 {
            for (a, b) in candidate
                .state()
                .cell_states()
                .iter()
                .zip(returned.as_ref().unwrap())
            {
                for axis in 0..3 {
                    near(
                        a.polarizations_t()[0][axis],
                        b.polarizations_t()[0][axis],
                        3e-8,
                    );
                }
            }
        }
        assert_eq!(state.potential_wb(), before.0);
        assert_eq!(state.cell_states(), before.1);
        assert_eq!(state.time_s(), before.2);
        state = dual
            .into_parts()
            .0
            .into_voltage_driven(preparation(), &mut || false)
            .unwrap();
        old = expected;
    }
    assert_eq!(original.time_s(), 0.0);
    assert_memory(original.cell_states(), Endpoint::default());
}

#[test]
fn refined_hysteretic_electrodes_resolve_free_charge_and_are_common_mode_and_control_invariant() {
    let card = data::card(500.0);
    let (positions, tets) = grid(2);
    let problem = data::problem(&card, positions, tets);
    let original = SpatialHystereticTerminalInduction::new(
        data::initial(&problem, Some(conductor(&card, 300.0))),
        &electrodes(&problem),
        preparation(),
        &mut || false,
    )
    .unwrap();
    let zero = vec![[0.0; 3]; problem.mesh().tets.len()];
    let mut state = original.clone();
    for voltage in [3.0, -1.5] {
        let expected = state
            .step_voltage_candidate(&[voltage, 0.0], &zero, DT, field_controls(), &mut || false)
            .unwrap();
        let shifted = state
            .step_voltage_candidate(
                &[1073741824.0 + voltage, 1073741824.0],
                &zero,
                DT,
                field_controls(),
                &mut || false,
            )
            .unwrap();
        assert_eq!(
            expected.state().potential_wb(),
            shifted.state().potential_wb()
        );
        assert_eq!(
            expected.state().cell_states(),
            shifted.state().cell_states()
        );
        assert_eq!(expected.step().terminals, shifted.step().terminals);
        let incoming: Vec<_> = expected
            .step()
            .terminals
            .iter()
            .map(|t| t.incoming_current_a)
            .collect();
        let current = state
            .clone()
            .into_current_driven(preparation(), &mut || false)
            .unwrap();
        assert!(
            current
                .step_current_candidate(
                    &[incoming[0], incoming[1] + 1.0],
                    &zero,
                    DT,
                    field_controls(),
                    &mut || false,
                )
                .is_err()
        );
        let actual = current
            .step_current_candidate(&incoming, &zero, DT, field_controls(), &mut || false)
            .unwrap();
        for (a, b) in actual
            .state()
            .cell_states()
            .iter()
            .zip(expected.state().cell_states())
        {
            for axis in 0..3 {
                near(a.flux_density_t()[axis], b.flux_density_t()[axis], 4e-7);
                near(
                    a.polarizations_t()[0][axis],
                    b.polarizations_t()[0][axis],
                    4e-7,
                );
            }
        }
        near(
            actual.step().terminals[1].voltage_to_reference_v,
            -voltage,
            4e-7,
        );
        near(
            actual.step().terminal_work_j,
            expected.step().terminal_work_j,
            4e-7,
        );
        near(actual.step().heat_j, expected.step().heat_j, 4e-7);
        near(
            actual.step().hysteresis_heat_j,
            expected.step().hysteresis_heat_j,
            4e-7,
        );
        assert!(
            actual.step().maximum_charge_residual_a
                <= field_controls().induction.induction.current_tolerance_a
        );
        for (a, b) in actual.step().cells.iter().zip(&expected.step().cells) {
            for (x, y) in a
                .current_density_a_m2
                .iter()
                .flatten()
                .zip(b.current_density_a_m2.iter().flatten())
            {
                near(*x, *y, 4e-7);
            }
        }
        state = actual
            .into_parts()
            .0
            .into_voltage_driven(preparation(), &mut || false)
            .unwrap();
    }
    assert!(
        state
            .cell_states()
            .iter()
            .any(|s| s.polarizations_t()[0].iter().any(|v| v.abs() > 1e-6))
    );
    assert_eq!(original.time_s(), 0.0);
}
