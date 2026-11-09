//! Coupled thermal/interface refinement and bounded caloric retry oracles.

use super::*;
use super::super::tests::{
    Case, CONCENTRATION, enthalpy, entropy, interface_shares, link, local_energies, near, rows,
};
use crate::miscibility::multicomponent::adaptive::MulticomponentMiscibilityTimeControl;
use crate::miscibility::multicomponent::phase_field::tests::uniform_matrix;
use fs_thermochem::solution::caloric::BinaryCaloricErrorV1;
use fs_thermochem::solution::multicomponent::caloric::MulticomponentCaloricErrorV1;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

fn accuracy(step: f64) -> MulticomponentNonisothermalPhaseFieldTimeControl {
    MulticomponentNonisothermalPhaseFieldTimeControl {
        thermal: MulticomponentNonisothermalTimeControl {
            composition: MulticomponentMiscibilityTimeControl {
                initial_step_s: step,
                minimum_step_s: 1e-8,
                maximum_step_s: step,
                absolute_mole_fraction_tolerance: 1e-8,
                relative_mole_fraction_tolerance: 0.0,
                maximum_attempts: 4_000,
            },
            absolute_temperature_tolerance_k: 1e-5,
            relative_temperature_tolerance: 0.0,
        },
        absolute_interfacial_energy_tolerance_j: 1e-4,
        relative_interfacial_energy_tolerance: 0.0,
    }
}

fn controls() -> MulticomponentNonisothermalSolve {
    MulticomponentNonisothermalSolve {
        energy_tolerance_j: 1e-6,
        temperature_tolerance_k: 1e-11,
        ..MulticomponentNonisothermalSolve::default()
    }
}

#[test]
fn coupled_refinement_converges_in_all_compositions_temperatures_and_original_energy() {
    let case = Case {
        b: uniform_matrix(3, 300.0),
        ..Case::default()
    };
    let initial = case.field();
    let snapshot = initial.clone();
    let candidate = initial
        .advance_candidate(1.0, accuracy(1.0), controls(), || false)
        .unwrap();
    assert_eq!(initial, snapshot);
    assert_eq!(candidate.field().elapsed_s(), 1.0);
    assert!(candidate.advance().rejected_intervals > 0);
    assert!(candidate.advance().maximum_accepted_error_ratio <= 1.0);
    let mut reference = initial.clone();
    for _ in 0..1_024 {
        reference = reference
            .step_candidate(1.0 / 1_024.0, controls(), || false)
            .unwrap()
            .into_parts()
            .0;
    }
    assert_eq!(reference.elapsed_s(), 1.0);
    for (&actual, &expected) in candidate
        .field()
        .mole_fractions()
        .iter()
        .zip(reference.mole_fractions())
    {
        near(actual, expected, 2e-6);
    }
    for (actual, expected) in candidate.field().temperatures().iter().zip(reference.temperatures()) {
        near(actual.value(), expected.value(), 2e-4);
    }
    let before = rows(&initial);
    let after = rows(candidate.field());
    let old_energy = local_energies(&case, &initial);
    let new_energy = local_energies(&case, candidate.field());
    let old_interface = interface_shares(&case, &before).iter().sum::<f64>();
    let new_interface = interface_shares(&case, &after).iter().sum::<f64>();
    let reference_interface = interface_shares(&case, &rows(&reference)).iter().sum::<f64>();
    near(new_interface, reference_interface, 2e-3);
    let mut enthalpy_change = 0.0;
    let mut entropy_change = 0.0;
    for cell in 0..2 {
        let old_t = initial.temperatures()[cell].value();
        let new_t = candidate.field().temperatures()[cell].value();
        enthalpy_change +=
            enthalpy(&case, &after[cell], new_t) - enthalpy(&case, &before[cell], old_t);
        entropy_change +=
            entropy(&case, &after[cell], new_t) - entropy(&case, &before[cell], old_t);
        let signed_energy = if cell == 0 {
            -candidate.advance().transferred_energy_j[0]
        } else {
            candidate.advance().transferred_energy_j[0]
        };
        near(new_energy[cell] - old_energy[cell], signed_energy, 2e-6);
        near(after[cell].iter().sum(), initial.network().capacities_mol()[cell], 1e-12);
    }
    for species in 0..3 {
        near(after[0][species] + after[1][species], before[0][species] + before[1][species], 1e-12);
        near(
            after[1][species] - before[1][species],
            candidate.advance().transferred_mol[0][species],
            1e-12,
        );
    }
    near(candidate.advance().homogeneous_enthalpy_change_j, enthalpy_change, 1e-7);
    near(
        candidate.advance().interfacial_energy_change_j,
        new_interface - old_interface,
        1e-8,
    );
    near(
        candidate.advance().energy_residual_j,
        enthalpy_change + new_interface - old_interface,
        1e-7,
    );
    near(enthalpy_change + new_interface - old_interface, 0.0, 2e-6);
    near(candidate.advance().entropy_change_j_per_k, entropy_change, 5e-9);
    assert!(entropy_change > 0.0);
    assert!(candidate.advance().completed_caloric_iterations > 0);
}

fn thermal_ring() -> Case {
    let cells = 16;
    let mut links = (0..cells - 1).map(|cell| link(cell, cell + 1)).collect::<Vec<_>>();
    links.push(link(0, cells - 1));
    for face in &mut links {
        face.area_m2 = 0.001;
    }
    Case {
        cp: vec![40.0; 3],
        cp_slope: vec![0.0; 3],
        volumes: vec![0.001; cells],
        fractions: (0..cells)
            .map(|cell| {
                let phase = 2.0 * std::f64::consts::PI * cell as f64 / cells as f64;
                let perturbation = 0.01 * (2.0 * phase).cos() + 0.0025 * (5.0 * phase).cos();
                vec![1.0 / 3.0 + perturbation, 1.0 / 3.0 - perturbation, 1.0 / 3.0]
            })
            .collect(),
        temperatures: vec![400.0; cells],
        links,
        diffusivities: vec![uniform_matrix(3, 0.2)],
        conductivity: 0.0,
        kappa: CONCENTRATION * R * 400.0,
        ..Case::default()
    }
}

fn edge_energies(field: &MulticomponentNonisothermalPhaseField, kappa: f64) -> Vec<f64> {
    field
        .network()
        .faces()
        .iter()
        .map(|face| {
            let left = field.cell_mole_fractions(face.cells[0]).unwrap();
            let right = field.cell_mole_fractions(face.cells[1]).unwrap();
            kappa * face.area_m2 / (4.0 * face.distance_m)
                * left.iter().zip(right).map(|(a, b)| (a - b) * (a - b)).sum::<f64>()
        })
        .collect()
}

#[test]
fn interface_l1_and_kelvin_accuracy_each_force_refinement_when_composition_is_loose() {
    let case = thermal_ring();
    let initial = case.field();
    let snapshot = initial.clone();
    let coarse = initial.step_candidate(1.0, controls(), || false).unwrap();
    let first = initial.step_candidate(0.5, controls(), || false).unwrap();
    let fine = first.field().step_candidate(0.5, controls(), || false).unwrap();
    let edge_errors = edge_energies(fine.field(), case.kappa)
        .into_iter()
        .zip(edge_energies(coarse.field(), case.kappa))
        .map(|(fine, coarse)| fine - coarse)
        .collect::<Vec<_>>();
    assert!(edge_errors.iter().any(|error| *error > 0.0));
    assert!(edge_errors.iter().any(|error| *error < 0.0));
    let signed_error = edge_errors.iter().sum::<f64>().abs();
    let l1_error = edge_errors.iter().map(|error| error.abs()).sum::<f64>();
    assert!(l1_error > 1.1 * signed_error);
    let maximum_temperature_error = coarse
        .field()
        .temperatures()
        .iter()
        .zip(fine.field().temperatures())
        .map(|(coarse, fine)| (coarse.value() - fine.value()).abs())
        .fold(0.0_f64, f64::max);
    assert!(maximum_temperature_error > 1e-7);

    let mut energy_time = accuracy(1.0);
    energy_time.thermal.composition.absolute_mole_fraction_tolerance = 1.0;
    energy_time.thermal.composition.maximum_attempts = 1;
    energy_time.thermal.absolute_temperature_tolerance_k = 1e6;
    energy_time.absolute_interfacial_energy_tolerance_j = 0.5 * (signed_error + l1_error);
    near(
        interface_error_ratio(coarse.field(), fine.field(), energy_time, &mut || false).unwrap(),
        l1_error / energy_time.absolute_interfacial_energy_tolerance_j,
        1e-7,
    );
    assert!(matches!(
        initial.advance_candidate(1.0, energy_time, controls(), || false),
        Err(MulticomponentNonisothermalPhaseFieldAdvanceError::WorkLimit {
            attempted_intervals: 1,
            accepted_intervals: 0,
            elapsed_s,
        }) if elapsed_s == 0.0
    ));

    let mut temperature_time = energy_time;
    temperature_time.absolute_interfacial_energy_tolerance_j = 1e6;
    temperature_time.thermal.absolute_temperature_tolerance_k = 0.5 * maximum_temperature_error;
    assert!(matches!(
        initial.advance_candidate(1.0, temperature_time, controls(), || false),
        Err(MulticomponentNonisothermalPhaseFieldAdvanceError::WorkLimit {
            attempted_intervals: 1,
            accepted_intervals: 0,
            elapsed_s,
        }) if elapsed_s == 0.0
    ));
    assert_eq!(initial, snapshot);
    for mut time in [energy_time, temperature_time] {
        time.thermal.composition.maximum_attempts = 64;
        let candidate = initial.advance_candidate(1.0, time, controls(), || false).unwrap();
        assert!(candidate.advance().rejected_intervals > 0);
        assert!(candidate.advance().maximum_accepted_error_ratio <= 1.0);
        assert_eq!(candidate.field().elapsed_s(), 1.0);
    }
}

#[test]
fn caloric_bracket_retries_reuse_rates_and_late_cancellation_preserves_the_whole_field() {
    let case = Case {
        cp: vec![40.0; 3],
        cp_slope: vec![0.0; 3],
        volumes: vec![0.001; 2],
        fractions: vec![vec![0.2, 0.3, 0.5]; 2],
        temperatures: vec![300.0, 600.0],
        conductivity: 100.0,
        ..Case::default()
    };
    let initial = case.field();
    let snapshot = initial.clone();
    // One explicit second would demand 1050K and -150K. Original caloric
    // bracket refusal must refine the interval, never clip the temperatures.
    assert!(matches!(
        initial.step_candidate(1.0, controls(), || false),
        Err(MulticomponentNonisothermalError::Caloric(
            MulticomponentCaloricErrorV1::Caloric(BinaryCaloricErrorV1::HeatOutsideBracket { .. })
        ))
    ));
    let mut time = accuracy(1.0);
    time.thermal.composition.absolute_mole_fraction_tolerance = 1.0;
    time.thermal.absolute_temperature_tolerance_k = 1000.0;
    time.absolute_interfacial_energy_tolerance_j = 1e6;
    time.thermal.composition.maximum_attempts = 1;
    assert!(matches!(
        initial.advance_candidate(1.0, time, controls(), || false),
        Err(MulticomponentNonisothermalPhaseFieldAdvanceError::WorkLimit {
            attempted_intervals: 1,
            accepted_intervals: 0,
            elapsed_s,
        }) if elapsed_s == 0.0
    ));
    time.thermal.composition.minimum_step_s = 1.0;
    assert!(matches!(
        initial.advance_candidate(1.0, time, controls(), || false),
        Err(MulticomponentNonisothermalPhaseFieldAdvanceError::Solver(
            MulticomponentNonisothermalError::Caloric(
                MulticomponentCaloricErrorV1::Caloric(BinaryCaloricErrorV1::HeatOutsideBracket { .. })
            )
        ))
    ));
    assert_eq!(initial, snapshot);

    time.thermal.composition.minimum_step_s = 1e-8;
    time.thermal.composition.maximum_attempts = 256;
    let mut polls = 0;
    let candidate = initial
        .advance_candidate(1.0, time, controls(), || {
            polls += 1;
            false
        })
        .unwrap();
    let report = candidate.advance();
    assert_eq!(candidate.field().elapsed_s(), 1.0);
    assert!(report.rejected_intervals > 0);
    assert!(report.accepted_intervals > 1);
    assert_eq!(report.attempted_intervals, report.accepted_intervals + report.rejected_intervals);
    assert!(report.completed_steps >= 3 * report.accepted_intervals);
    assert!(report.completed_steps <= 3 * report.attempted_intervals);
    assert!(
        report.completed_maxwell_stefan_evaluations
            <= report.accepted_intervals + report.attempted_intervals
    );
    assert!(report.completed_maxwell_stefan_evaluations < 2 * report.attempted_intervals);
    assert!(report.completed_caloric_iterations > 0);
    assert_eq!(candidate.field().component_moles(), initial.component_moles());
    assert_eq!(report.transferred_mol, vec![vec![0.0; 3]]);
    assert_eq!(report.interfacial_energy_change_j, 0.0);
    let old_energy = local_energies(&case, &initial);
    let new_energy = local_energies(&case, candidate.field());
    near(new_energy.iter().sum(), old_energy.iter().sum(), 2e-6);
    near(new_energy[0] - old_energy[0], -report.transferred_energy_j[0], 2e-6);
    near(new_energy[1] - old_energy[1], report.transferred_energy_j[0], 2e-6);
    assert!(report.entropy_change_j_per_k > 0.0);
    assert!(candidate
        .field()
        .temperatures()
        .iter()
        .all(|temperature| (250.0..=700.0).contains(&temperature.value())));

    let mut cancelled_polls = 0;
    assert!(matches!(
        initial.advance_candidate(1.0, time, controls(), || {
            cancelled_polls += 1;
            cancelled_polls == polls
        }),
        Err(MulticomponentNonisothermalPhaseFieldAdvanceError::Solver(
            MulticomponentNonisothermalError::Cancelled
        ))
    ));
    assert_eq!(cancelled_polls, polls);
    assert_eq!(initial, snapshot);
    assert_eq!(
        candidate,
        initial.advance_candidate(1.0, time, controls(), || false).unwrap()
    );
}
