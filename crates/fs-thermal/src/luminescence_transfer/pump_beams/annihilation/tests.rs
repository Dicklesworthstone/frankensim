use super::super::super::pump_column::{SaturablePumpCell, SaturablePumpColumn};
use super::*;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-290),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

fn model(singlet_return: f64) -> SingletTripletDecay {
    SingletTripletDecay::try_new([4e-19, 3e-19], [singlet_return, 0.0], [0.0; 2], [0.0; 2]).unwrap()
}

fn path(indices: &[usize], depth: f64) -> SaturablePumpPath<1> {
    SaturablePumpPath {
        photon_energies_j: [4e-19],
        segments: indices
            .iter()
            .map(|&cell_index| SaturablePumpSegment {
                cell_index,
                optical_depth_per_ground: [[0.0, depth]],
                passive_optical_depth: [0.0],
            })
            .collect(),
    }
}

fn controls() -> (PumpBeamsLimits, AnnihilationLimits) {
    (
        PumpBeamsLimits {
            column: PumpColumnLimits {
                relative_tolerance: 2e-4,
                maximum_substeps: 2048,
                maximum_cell_evaluations: 300_000,
                ..PumpColumnLimits::default()
            },
            ..PumpBeamsLimits::default()
        },
        AnnihilationLimits {
            relative_tolerance: 2e-4,
            maximum_substeps: 2048,
            ..AnnihilationLimits::default()
        },
    )
}

fn steady_cell(step: &SaturableSingletTripletStep, initial: [f64; 3], pairs: f64) {
    for state in 0..3 {
        close(step.populations[state], initial[state], 2e-10);
    }
    close(step.annihilation_events, pairs, 2e-10);
    close(step.annihilation_singlets, 0.8 * pairs, 2e-10);
    close(step.annihilation_heat_j, pairs * 2.8e-19, 2e-10);
    close(step.pump_events[1], 2.0 * pairs, 2e-10);
    close(step.radiative_events[0], 0.8 * pairs, 2e-10);
    assert_eq!(step.nonradiative_events, [0.0; 2]);
    assert_eq!(step.crossing_events, [0.0; 2]);
    assert_eq!(step.reverse_crossing_heat_absorbed_j, 0.0);
}

#[test]
fn g1_pair_recovery_and_depletion_hold_independent_column_and_opposed_equilibrium() {
    let count: f64 = 1e12;
    let coefficient = 1e-12;
    let singlet_return = 1e9;
    let triplets = 0.75 * count;
    let pair_rate = coefficient * triplets * triplets;
    let singlets = 0.8 * pair_rate / singlet_return;
    let initial = [count - triplets - singlets, singlets, triplets];
    let depth = 2.0 / count;
    let transmission = (-depth * initial[0]).exp();
    // Stationary T requires absorbed photons=2aT². Stationary S requires
    // kS*S=phi*aT². Beer attenuation fixes the independent external input.
    let input = 2.0 * pair_rate / (-(-depth * initial[0]).exp_m1());
    let law = Some(TripletAnnihilation::try_new(coefficient, 1.0, 0.8).unwrap());
    let (limits, mut pair_limits) = controls();
    pair_limits.maximum_substeps = 128;
    let column = SaturablePumpColumn::try_new(
        &[SaturablePumpCell {
            model: model(singlet_return),
            optical_depth_per_ground: [[0.0, depth]],
            passive_optical_depth: [0.0],
        }],
        [4e-19],
        1,
        || false,
    )
    .unwrap();
    let column_step = column
        .advance_with_annihilation(
            &[initial],
            [input],
            1.0,
            &[law],
            limits.column,
            pair_limits,
            || false,
        )
        .unwrap();
    steady_cell(&column_step.cells[0].kinetics, initial, pair_rate);
    close(
        column_step.transmitted_photons[0] / input,
        transmission,
        2e-10,
    );

    // Two mirrored cells and opposing beams have the same local summed input
    // R*(1+exp(-tau)) and removal R*(1-exp(-2*tau)). Each original inventory
    // therefore has the SAME exact equilibrium, advanced only once.
    let opposed = SaturablePumpBeams::try_new(
        &[model(singlet_return); 2],
        &[path(&[0, 1], 1.0 / count), path(&[1, 0], 1.0 / count)],
        2,
        2,
        4,
        || false,
    )
    .unwrap();
    let opposed_step = opposed
        .advance_with_annihilation(
            &[initial; 2],
            &[[input]; 2],
            1.0,
            &[law; 2],
            limits,
            pair_limits,
            || false,
        )
        .unwrap();
    for cell in &opposed_step.cells {
        steady_cell(cell, initial, pair_rate);
    }
    for beam in &opposed_step.beams {
        close(beam.transmitted_photons[0] / input, transmission, 2e-10);
        close(beam.absorbed_photons[0], 2.0 * pair_rate, 2e-10);
        assert!(beam.photon_residual[0].abs() < 1e-11 * input);
    }
    assert!(opposed_step.coupling_iterations > 0);
    assert!(opposed_step.substeps <= 128);
    assert!(opposed_step.residual_j.abs() < 1e-17);
}

#[test]
fn g3_pair_return_changes_actual_pulse_absorption_and_is_invariant_to_beam_splitting() {
    let count = 1e12;
    let law = Some(TripletAnnihilation::try_new(5e-12, 1.0, 0.5).unwrap());
    let (limits, pair_limits) = controls();
    let column = SaturablePumpColumn::try_new(
        &[SaturablePumpCell {
            model: model(1e6),
            optical_depth_per_ground: [[0.0, 5.0 / count]],
            passive_optical_depth: [0.0],
        }],
        [4e-19],
        1,
        || false,
    )
    .unwrap();
    let initial = [[count, 0.0, 0.0]];
    let without = column
        .advance(&initial, [3e12], 0.5, limits.column, || false)
        .unwrap();
    let with = column
        .advance_with_annihilation(
            &initial,
            [3e12],
            0.5,
            &[law],
            limits.column,
            pair_limits,
            || false,
        )
        .unwrap();
    // Pair encounters replenish actual ground absorbers DURING the pulse.
    // A pair kick after the optical solve would leave these photons unchanged.
    assert!(with.absorbed_photons[0] > 1.1 * without.absorbed_photons[0]);
    assert!(with.cells[0].kinetics.annihilation_events > 0.0);
    assert!(with.cells[0].kinetics.radiative_events[0] > 0.0);
    let split = SaturablePumpBeams::try_new(
        &[model(1e6)],
        &[path(&[0], 5.0 / count), path(&[0], 5.0 / count)],
        1,
        2,
        2,
        || false,
    )
    .unwrap();
    let split_step = split
        .advance_with_annihilation(
            &initial,
            &[[0.9e12], [2.1e12]],
            0.5,
            &[law],
            limits,
            pair_limits,
            || false,
        )
        .unwrap();
    for state in 0..3 {
        close(
            split_step.cells[0].populations[state],
            with.cells[0].kinetics.populations[state],
            5e-4,
        );
    }
    close(
        split_step.cells[0].annihilation_events,
        with.cells[0].kinetics.annihilation_events,
        5e-4,
    );
    for (beam, fraction) in split_step.beams.iter().zip([0.3, 0.7]) {
        close(
            beam.absorbed_photons[0],
            fraction * with.absorbed_photons[0],
            5e-4,
        );
        close(
            beam.transmitted_photons[0],
            fraction * with.transmitted_photons[0],
            5e-4,
        );
    }
    close(split_step.cells[0].populations.iter().sum(), count, 2e-12);
    assert!(split_step.residual_j.abs() < 1e-17);
}

#[test]
fn g4_dark_shared_and_unvisited_pairs_recover_once_with_zero_law_and_failure_identity() {
    let (limits, pair_limits) = controls();
    let law = Some(TripletAnnihilation::try_new(1e-12, 1.0, 0.8).unwrap());
    let initial = [[3e11, 0.0, 7e11], [2e11, 1e11, 7e11], [1e11, 0.0, 9e11]];
    let network = SaturablePumpBeams::try_new(
        &[model(1e6); 3],
        &[path(&[0, 1], 1e-12), path(&[1, 0], 1e-12)],
        3,
        2,
        4,
        || false,
    )
    .unwrap();
    let duration = 0.1;
    let step = network
        .advance_with_annihilation(
            &initial,
            &[[0.0]; 2],
            duration,
            &[law; 3],
            limits,
            pair_limits,
            || false,
        )
        .unwrap();
    let local = SaturableSingletTripletPlan::prepare(
        model(1e6),
        [0.0; 2],
        [0.0; 2],
        duration,
        limits.column.kinetic,
        || false,
    )
    .unwrap();
    for index in 0..3 {
        let reference = local
            .advance_with_annihilation(initial[index], law, pair_limits, || false)
            .unwrap();
        for state in 0..3 {
            close(
                step.cells[index].populations[state],
                reference.kinetics.populations[state],
                5e-4,
            );
        }
        close(
            step.cells[index].annihilation_events,
            reference.kinetics.annihilation_events,
            5e-4,
        );
        assert!(step.cells[index].radiative_events[0] > 0.0);
    }
    assert!(step.substeps > 1);
    assert_eq!(step.coupling_iterations, 0);
    for beam in &step.beams {
        assert_eq!(beam.incident_photons, [0.0]);
    }
    let original = network
        .advance(&initial, &[[0.0]; 2], duration, limits, || false)
        .unwrap();
    let zero = Some(TripletAnnihilation::try_new(0.0, 1.0, 0.8).unwrap());
    assert_eq!(
        network
            .advance_with_annihilation(
                &initial,
                &[[0.0]; 2],
                duration,
                &[zero; 3],
                limits,
                pair_limits,
                || false
            )
            .unwrap(),
        original
    );
    let mut work_limit = limits;
    work_limit.column.maximum_cell_evaluations = 1;
    assert!(
        network
            .advance_with_annihilation(
                &initial,
                &[[0.0]; 2],
                duration,
                &[law; 3],
                work_limit,
                pair_limits,
                || false
            )
            .is_err()
    );
    let mut checks = 0;
    assert!(
        network
            .advance_with_annihilation(
                &initial,
                &[[0.0]; 2],
                duration,
                &[law; 3],
                limits,
                pair_limits,
                || {
                    checks += 1;
                    checks > 100
                }
            )
            .is_err()
    );
    assert_eq!(
        network
            .advance_with_annihilation(
                &initial,
                &[[0.0]; 2],
                duration,
                &[law; 3],
                limits,
                pair_limits,
                || false
            )
            .unwrap(),
        step
    );
}
